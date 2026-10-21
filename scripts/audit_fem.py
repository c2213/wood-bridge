#!/usr/bin/env python3
"""構造監査: 2D フレーム FEM（線形・座屈）で破壊荷重を見積もり、CP = P / W^1.5 を構造形式間で比較する。

標準ライブラリのみ。片側1構面（リブまたはトラス1面）を 2D でモデル化し、全体荷重 P の半分を載せる。
単位: mm, N, MPa, g。

    python3 scripts/audit_fem.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_draft as md  # noqa: E402

# ---- 材料値（文献の代表値。実測で置き換える）----
# バルサ: USDA Wood Handbook (FPL-GTR-190) 表5-4a、比重 0.16。圧縮は市販材のばらつきを見て 14.9→12 に低減。
# 接着: 低密度材の PVA 接合は木部破壊で決まる（同 Handbook 10章）。交差接合は木目直交の
#       ころがりせん断（平行せん断 2.1 MPa の約 25%）、添え板の重ね接着は繊維方向せん断を低減して使う。
MAT = dict(
    E=3400.0,       # 繊維方向ヤング率
    fc=12.0,        # 繊維方向圧縮強度
    ft=14.0,        # 繊維方向引張強度
    fm=21.6,        # 曲げ強度（MOR）
    fc90=0.6,       # 繊維直角方向の圧縮（めり込み）強度（文献推定 0.4〜1.0）
    tau=0.5,        # 交差接合（木目直交）のせん断強度（ころがりせん断 0.4〜0.6）
    tau_lap=1.0,    # 添え板の重ね接着のせん断強度
    rho=0.16e-3,    # g/mm3
)
HANGERS = {  # EA [N]、結び目込みの破断荷重 [N]（JAFS 標準直径・lb 表示から換算）
    "ナイロン3号": dict(EA=2800 * 0.0638, T=53 * 0.75),
    "PE1号": dict(EA=1900.0, T=78 * 0.40),
    "バルサ3mm角": dict(EA=9 * MAT["E"], T=60.0),
}
GLUE_PER_JOINT = 0.03     # g/接合（乾燥後の仮定）
GUSSET_L = 15.0           # 部材端の添え板（厚紙ガセット）の重ね長さ。両面接着
SPLICE_L = 30.0           # 床版縦桁（タイ）継手の添え板の重ね長さ。両面接着
DECK_OPTIONS = {"厚紙 t=0.5": 0.5 * 90 * 505 * 0.70e-3, "厚紙 t=0.3": 0.3 * 90 * 505 * 0.70e-3,
                "ケント紙 t=0.2": 0.2 * 90 * 505 * 0.70e-3, "床版なし（横桁のみ・ルール要確認）": 0.0}
CROSS_W = 5.6             # 床版横桁（現案）の重量 g


# ---------------------------------------------------------------- 線形代数（帯行列 LDLᵀ）
def ldl(A, n, bw):
    for k in range(n):
        d = A[k][k]
        if abs(d) < 1e-300:
            d = A[k][k] = 1e-300
        hi = min(n, k + bw + 1)
        for i in range(k + 1, hi):
            aik = A[i][k]
            if aik == 0.0:
                continue
            f = aik / d
            Ai = A[i]
            for j in range(k + 1, i + 1):
                ajk = A[j][k]
                if ajk != 0.0:
                    Ai[j] -= f * ajk
        for i in range(k + 1, hi):
            A[i][k] /= d
    return A


def ldl_solve(A, n, bw, b):
    y = b[:]
    for i in range(n):
        for k in range(max(0, i - bw), i):
            y[i] -= A[i][k] * y[k]
    for i in range(n):
        y[i] /= A[i][i]
    for i in range(n - 1, -1, -1):
        for k in range(i + 1, min(n, i + bw + 1)):
            y[i] -= A[k][i] * y[k]
    return y


# ---------------------------------------------------------------- 2D フレーム
class Model:
    def __init__(self):
        self.nodes, self.elems, self.fix, self.loads = [], [], {}, {}

    def node(self, x, y):
        self.nodes.append((x, y))
        return len(self.nodes) - 1

    def frame(self, i, j, EA, EI, hi=False, hj=False, **info):
        self.elems.append(dict(t="F", i=i, j=j, EA=EA, EI=EI, hi=hi, hj=hj, **info))

    def truss(self, i, j, EA, **info):
        self.elems.append(dict(t="T", i=i, j=j, EA=EA, **info))

    def beam_line(self, pts, nsub, EA, EI, **info):
        """pts を順に結ぶ連続梁を nsub 分割で張り、各 pts のノード番号を返す。"""
        ids = [pts[0]]
        for a, b in zip(pts, pts[1:]):
            (xa, ya), (xb, yb) = self.nodes[a], self.nodes[b]
            prev = a
            for k in range(1, nsub):
                t = k / nsub
                n = self.node(xa + (xb - xa) * t, ya + (yb - ya) * t)
                self.frame(prev, n, EA, EI, **info)
                prev = n
            self.frame(prev, b, EA, EI, **info)
            ids.append(b)
        return ids

    # -- DOF 番号付け（ノードを x 順に並べて帯幅を小さくする）
    def number(self):
        rot_keys = {}
        for e_id, e in enumerate(self.elems):
            if e["t"] != "F":
                continue
            for end, h in (("i", e["hi"]), ("j", e["hj"])):
                nd = e[end]
                key = ("h", e_id, end) if h else ("n", nd)
                rot_keys.setdefault(nd, [])
                if key not in rot_keys[nd]:
                    rot_keys[nd].append(key)
        order = sorted(range(len(self.nodes)), key=lambda n: (self.nodes[n][0], self.nodes[n][1]))
        dof, idx = {}, 0
        for n in order:
            dof[("u", n)], dof[("v", n)] = idx, idx + 1
            idx += 2
            for key in rot_keys.get(n, []):
                dof[key] = idx
                idx += 1
        self.dof, self.ndof = dof, idx
        for e_id, e in enumerate(self.elems):
            ds = [dof[("u", e["i"])], dof[("v", e["i"])]]
            if e["t"] == "F":
                ds.append(dof[("h", e_id, "i")] if e["hi"] else dof[("n", e["i"])])
            ds += [dof[("u", e["j"])], dof[("v", e["j"])]]
            if e["t"] == "F":
                ds.append(dof[("h", e_id, "j")] if e["hj"] else dof[("n", e["j"])])
            e["dofs"] = ds

    def geom(self, e):
        (xa, ya), (xb, yb) = self.nodes[e["i"]], self.nodes[e["j"]]
        L = math.hypot(xb - xa, yb - ya)
        return L, (xb - xa) / L, (yb - ya) / L

    def k_local(self, e, N=None):
        """N=None なら弾性剛性、N を与えると幾何剛性（引張正）。"""
        L, c, s = self.geom(e)
        if e["t"] == "T":
            if N is None:
                k = e["EA"] * (0 if e.get("off") else 1) / L
                kl = [[k, 0, -k, 0], [0, 0, 0, 0], [-k, 0, k, 0], [0, 0, 0, 0]]
            else:
                g = N / L
                kl = [[0, 0, 0, 0], [0, g, 0, -g], [0, 0, 0, 0], [0, -g, 0, g]]
            T = [[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]]
            return kl, T
        if N is None:
            EA, EI = e["EA"], e["EI"]
            a, b2, b3, b4, b5 = EA / L, 12 * EI / L ** 3, 6 * EI / L ** 2, 4 * EI / L, 2 * EI / L
            kl = [[a, 0, 0, -a, 0, 0], [0, b2, b3, 0, -b2, b3], [0, b3, b4, 0, -b3, b5],
                  [-a, 0, 0, a, 0, 0], [0, -b2, -b3, 0, b2, -b3], [0, b3, b5, 0, -b3, b4]]
        else:
            g = N / L
            kl = [[0] * 6, [0, 6 / 5 * g, L / 10 * g, 0, -6 / 5 * g, L / 10 * g],
                  [0, L / 10 * g, 2 * L * L / 15 * g, 0, -L / 10 * g, -L * L / 30 * g], [0] * 6,
                  [0, -6 / 5 * g, -L / 10 * g, 0, 6 / 5 * g, -L / 10 * g],
                  [0, L / 10 * g, -L * L / 30 * g, 0, -L / 10 * g, 2 * L * L / 15 * g]]
        T = [[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
             [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]]
        return kl, T

    def assemble(self, forces=None):
        free = [d for d in range(self.ndof) if d not in self.fixed]
        pos = {d: k for k, d in enumerate(free)}
        n = len(free)
        K = [[0.0] * n for _ in range(n)]
        bw = 0
        for e_id, e in enumerate(self.elems):
            kl, T = self.k_local(e, None if forces is None else forces[e_id]["N"])
            m = len(kl)
            # Kg = Tᵀ kl T
            kt = [[sum(kl[a][c] * T[c][b] for c in range(m)) for b in range(m)] for a in range(m)]
            kg = [[sum(T[c][a] * kt[c][b] for c in range(m)) for b in range(m)] for a in range(m)]
            ps = [pos.get(d) for d in e["dofs"]]
            for a in range(m):
                if ps[a] is None:
                    continue
                for b in range(m):
                    if ps[b] is None or ps[b] > ps[a]:
                        continue
                    K[ps[a]][ps[b]] += kg[a][b]
                    bw = max(bw, ps[a] - ps[b])
        return K, n, bw, free, pos

    def solve(self):
        self.fixed = set(self.dof[k] for k in self.fix)
        for _ in range(20):  # 圧縮になった釣り糸を外して再計算
            K, n, bw, free, pos = self.assemble()
            ldl(K, n, bw)
            b = [0.0] * n
            for key, f in self.loads.items():
                b[pos[self.dof[key]]] += f
            x = ldl_solve(K, n, bw, b)
            u = [0.0] * self.ndof
            for d, k in pos.items():
                u[d] = x[k]
            forces = self.element_forces(u)
            slack = [k for k, e in enumerate(self.elems)
                     if e.get("tonly") and not e.get("off") and forces[k]["N"] < 0]
            if not slack:
                break
            for k in slack:
                self.elems[k]["off"] = True
        self.u, self.forces = u, forces
        return u, forces

    def element_forces(self, u):
        out = []
        for e in self.elems:
            kl, T = self.k_local(e)
            m = len(kl)
            dg = [u[d] for d in e["dofs"]]
            dl = [sum(T[a][b] * dg[b] for b in range(m)) for a in range(m)]
            f = [sum(kl[a][b] * dl[b] for b in range(m)) for a in range(m)]
            if e["t"] == "T":
                out.append(dict(N=f[2], V=0.0, M=0.0))
            else:
                out.append(dict(N=f[3], V=abs(f[1]), M=max(abs(f[2]), abs(f[5]))))
        return out

    def buckling(self, lo=1e-3, hi=1e7, iters=60):
        """K + λ·Kg が初めて正定値でなくなる λ を Sturm 列（負のピボット数）の二分法で求める。"""
        def neg(lam):
            K, n, bw, _, _ = self.assemble()
            G, *_ = self.assemble(self.forces)
            for a in range(n):
                Ka, Ga = K[a], G[a]
                for b in range(max(0, a - bw), a + 1):
                    Ka[b] += lam * Ga[b]
            ldl(K, n, bw)
            return sum(1 for a in range(n) if K[a][a] < 0)
        if neg(hi) == 0:
            return float("inf")
        for _ in range(iters):
            mid = math.sqrt(lo * hi)
            if neg(mid) > 0:
                hi = mid
            else:
                lo = mid
            if hi / lo < 1.002:
                break
        return hi


# ---------------------------------------------------------------- 構造モデル
def sq(s):
    return s * s, s ** 4 / 12, s ** 3 / 6  # A, I, Z


def add_deck(m, xs, s_str=5.0, tie=True, anchor=False):
    """床版縦桁（タイ）。支点: 左ピン・右ローラー（anchor=True なら右もピン）。"""
    A, I, _ = sq(s_str)
    E = MAT["E"]
    nodes = [m.node(x, s_str / 2) for x in xs]
    ids = m.beam_line(nodes, 2, E * A * (1 if tie else 1e-4), E * I, kind="S", rows=1, s=s_str)
    m.fix[("u", ids[0])] = m.fix[("v", ids[0])] = m.fix[("v", ids[-1])] = True
    if anchor:
        m.fix[("u", ids[-1])] = True
    return ids


def davinci(hanger, s_long=5.0, s_trans=3.0, rows=(2, 1), rigid_ends=False, tie=True, anchor=False,
            nsub=3, n_seg=8, s_str=5.0):
    r, alpha, delta, chord, rise = md.solve_arch(n_seg, s_long, s_trans)
    theta = math.asin(md.SPAN / 2 / r)

    def ph(i):
        return -theta + i * alpha

    hh0 = s_trans / 2 * (abs(math.cos(ph(0))) + abs(math.sin(ph(0))))
    yc = md.STRINGER + hh0 - r * math.cos(theta)

    def P(rr, i):
        return md.SPAN / 2 + rr * math.sin(ph(i)), yc + rr * math.cos(ph(i))

    m = Model()
    E = MAT["E"]
    xs = [P(r, i)[0] for i in range(n_seg + 1)]
    deck = add_deck(m, xs, s_str=s_str, tie=tie, anchor=anchor)
    tn = [m.node(*P(r, i)) for i in range(n_seg + 1)]
    A1, I1, _ = sq(s_long)
    big = 1e3

    def arm(a, t, fam_rows):  # 部材側は剛、横材側はピン（接触＋接着）
        m.frame(a, t, E * A1 * big, E * I1 * big, hj=True, kind="J", rows=fam_rows, s=s_long)

    ends = {}  # (family, node) -> member end node（rigid_ends 用）
    mids = {}
    members = [("A", k, k + 2) for k in range(0, n_seg, 2)] + [("B", k, k + 2) for k in range(1, n_seg - 1, 2)] \
        + [("B", 0, 1), ("B", n_seg - 1, n_seg)]
    for fam, i0, i1 in members:
        nr = rows[0] if fam == "A" else rows[1]
        pts = []
        for i in (i0, i1):
            key = (fam, i)
            if rigid_ends and key in ends:
                p = ends[key]
            else:
                p = m.node(*P(r + delta, i))
                arm(p, tn[i], nr)
                ends[key] = p
            pts.append(p)
        if i1 - i0 == 2:
            pm = m.node(*P(r - delta, i0 + 1))
            arm(pm, tn[i0 + 1], nr)
            mids[i0 + 1] = pm
            pts.insert(1, pm)
        m.beam_line(pts, nsub, E * A1 * nr, E * I1 * nr, kind=fam, rows=nr, s=s_long)
    # 端部横材 T0/T8 → 縦桁
    for t, d in ((tn[0], deck[0]), (tn[-1], deck[-1])):
        m.frame(d, t, E * A1 * big, E * I1 * big, hj=True, kind="L", rows=1, s=s_long)
    # ハンガー（奇数節点は横材、偶数節点は B 系中央部から）
    h = HANGERS[hanger]
    for i in range(1, n_seg):
        top = tn[i] if i % 2 else mids[i]
        m.truss(deck[i], top, h["EA"], kind="H", tonly=hanger != "バルサ3mm角", cap=h["T"])
    loads_at(m, deck, xs)
    length = sum(rows[0 if fam == "A" else 1] * math.dist(P(r + delta, a), P(r + delta, b)) + 2 * md.OVERHANG
                 for fam, a, b in members)
    w_rib = length * s_long ** 2 * MAT["rho"]
    w_t = (8 * 22 + 5 * 134) * s_trans ** 2 * MAT["rho"] / 2  # 片側分
    joints = (sum(rows) + 1) * (n_seg + 1) + 2 * n_seg       # 交差・ハンガー接合のおおよその数（片側）
    return m, dict(w_side=w_rib + w_t, joints=joints, rise=rise, crown=P(r, n_seg // 2)[1], s_str=s_str)


def loads_at(m, deck, xs):
    for x in md.LOADS:
        i = min(range(1, len(xs) - 1), key=lambda k: abs(xs[k] - x))
        m.loads[("v", deck[i])] = m.loads.get(("v", deck[i]), 0.0) - 0.25  # 全体 P=1 → 片側 P/2 を2点
    m.load_nodes = [deck[min(range(1, len(xs) - 1), key=lambda k: abs(xs[k] - x))] for x in md.LOADS]


def polyarch(hanger, s=5.0, rows=1, rigid=True, n_seg=8, crown=125.9, tie=True, anchor=False, nsub=3, s_str=5.0):
    """通常の（編まない）タイドアーチ。同じライズの円弧上に節点を置き、節点間を1本の直材で結ぶ。"""
    half = md.SPAN / 2
    y0 = md.STRINGER + s / 2
    f = crown - y0
    r = (half * half + f * f) / (2 * f)
    theta = math.asin(half / r)
    m = Model()
    E = MAT["E"]
    A1, I1, _ = sq(s)
    pts_xy = [(half + r * math.sin(-theta + 2 * theta * i / n_seg),
               y0 - r * math.cos(theta) + r * math.cos(-theta + 2 * theta * i / n_seg)) for i in range(n_seg + 1)]
    xs = [p[0] for p in pts_xy]
    deck = add_deck(m, xs, s_str=s_str, tie=tie, anchor=anchor)
    an = [m.node(*p) for p in pts_xy]
    for a, b in zip(an, an[1:]):
        if rigid:
            m.beam_line([a, b], nsub, E * A1 * rows, E * I1 * rows, kind="A", rows=rows, s=s, gusset=True)
        else:
            m.frame(a, b, E * A1 * rows, E * I1 * rows, hi=True, hj=True, kind="A", rows=rows, s=s, gusset=True)
    for t, d in ((an[0], deck[0]), (an[-1], deck[-1])):
        m.frame(d, t, E * A1 * 1e3, E * I1 * 1e3, hj=True, kind="L", rows=1, s=s)
    h = HANGERS[hanger]
    for i in range(1, n_seg):
        m.truss(deck[i], an[i], h["EA"], kind="H", tonly=hanger != "バルサ3mm角", cap=h["T"])
    loads_at(m, deck, xs)
    length = sum(math.dist(p, q) for p, q in zip(pts_xy, pts_xy[1:])) * rows
    w_t = 5 * 134 * 9 * MAT["rho"] / 2
    return m, dict(w_side=length * s * s * MAT["rho"] + w_t, joints=2 * n_seg + 2, crown=crown, s_str=s_str)


def pratt(s_chord=5.0, s_web=3.0, h=90.0, nsub=2):
    """下路式プラットトラス。下弦材＝床版縦桁、上弦材は有効高さの上。"""
    xs = [0.0, 87.5, 175.0, 250.0, 325.0, 412.5, 500.0]
    m = Model()
    E = MAT["E"]
    deck = add_deck(m, xs, s_str=s_chord)
    yt = md.STRINGER / 2 + h
    top = {i: m.node(xs[i], yt) for i in range(1, 6)}
    Ac, Ic, _ = sq(s_chord)
    Aw, Iw, _ = sq(s_web)
    m.beam_line([top[i] for i in range(1, 6)], nsub, E * Ac, E * Ic, kind="TC", rows=1, s=s_chord)
    for a, b in ((deck[0], top[1]), (deck[6], top[5])):
        m.beam_line([a, b], nsub, E * Ac, E * Ic, kind="EP", rows=1, s=s_chord)
    web = []
    for i in range(1, 6):
        web.append((deck[i], top[i]))
    web += [(top[1], deck[2]), (top[2], deck[3]), (top[4], deck[3]), (top[5], deck[4])]
    for a, b in web:
        m.beam_line([a, b], nsub, E * Aw, E * Iw, kind="W", rows=1, s=s_web)
    loads_at(m, deck, xs)
    L = lambda a, b: math.dist(m.nodes[a], m.nodes[b])
    length_c = L(top[1], top[5]) + 2 * L(deck[0], top[1])
    length_w = sum(L(a, b) for a, b in web)
    w_t = 5 * 134 * 9 * MAT["rho"] / 2
    w = (length_c * s_chord ** 2 + length_w * s_web ** 2) * MAT["rho"] + w_t
    return m, dict(w_side=w, joints=len(web) * 2 + 6, crown=yt, s_str=s_chord)


# ---------------------------------------------------------------- 照査
def evaluate(m):
    m.number()
    u, F = m.solve()
    worst = dict(member=(0, None), joint=(0, None), hanger=(0, None), oop=(0, None))
    tie_N = 0.0
    for e, f in zip(m.elems, F):
        k = e.get("kind")
        if k in ("A", "B", "TC", "EP", "W"):
            A, I, Z = sq(e["s"])
            N, M = f["N"] / e["rows"], f["M"] / e["rows"]
            fa = MAT["ft"] if N > 0 else MAT["fc"]
            ua, ub = abs(N) / (A * fa), M / (Z * MAT["fm"])
            if ua + ub > worst["member"][0]:
                worst["member"] = (ua + ub, "%s 軸%.0f%%・曲げ%.0f%%" % (k, 100 * ua / (ua + ub), 100 * ub / (ua + ub)))
            if N < 0:  # 面外座屈（横材・節点間を座屈長、ピン）
                Lb = 75.0 if k in ("A", "B") else math.dist(m.nodes[e["i"]], m.nodes[e["j"]]) * 2
                ncr = math.pi ** 2 * MAT["E"] * I / Lb ** 2
                if -N / ncr > worst["oop"][0]:
                    worst["oop"] = (-N / ncr, "%s Lb=%.0f" % (k, Lb))
            if k in ("TC", "EP", "W") or (k == "A" and e.get("gusset")):
                cap = 2 * e["s"] * GUSSET_L * MAT["tau_lap"]   # 端部ガセット（引張のみ照査）
                if max(0.0, N) / cap > worst["joint"][0]:
                    worst["joint"] = (max(0.0, N) / cap, "%s 端部ガセット引張" % k)
        elif k == "S":
            tie_N = max(tie_N, f["N"])
            cap = 2 * e["s"] * SPLICE_L * MAT["tau_lap"]
            if max(0.0, f["N"]) / cap > worst["joint"][0]:
                worst["joint"] = (max(0.0, f["N"]) / cap, "縦桁（タイ）継手")
            A, I, Z = sq(e["s"])
            fa = MAT["ft"] if f["N"] > 0 else MAT["fc"]
            ut = abs(f["N"]) / (A * fa) + f["M"] / (Z * MAT["fm"])
            if ut > worst["member"][0]:
                worst["member"] = (ut, "床版縦桁（曲げ・軸）")
        elif k == "J":
            area = e["s"] * md.T_TRANS
            n = e["rows"]
            ub = max(0.0, -f["N"]) / n / (area * MAT["fc90"])
            us = f["V"] / n / (area * MAT["tau"])
            ut = max(0.0, f["N"]) / n / (area * MAT["tau"])
            u_ = max(ub, us, ut)
            if u_ > worst["joint"][0]:
                which = "めり込み" if u_ == ub else ("接着せん断" if u_ == us else "接着引張")
                worst["joint"] = (u_, which)
        elif k == "H" and not e.get("off"):
            if f["N"] / e["cap"] > worst["hanger"][0]:
                worst["hanger"] = (f["N"] / e["cap"], "ハンガー %.2f·P" % f["N"])
    lam = m.buckling()
    defl = -min(u[m.dof[("v", n)]] for n in m.load_nodes)
    caps = {k: (1 / v if v > 0 else float("inf"), why) for k, (v, why) in worst.items()}
    caps["buckle"] = (lam, "面内座屈")
    gov = min(caps.items(), key=lambda kv: kv[1][0])
    return dict(caps=caps, gov=gov, tie=tie_N, defl=defl)


def weight(info, deck_g):
    glue = 2 * info["joints"] * GLUE_PER_JOINT
    stringer = 2 * 514 * info["s_str"] ** 2 * MAT["rho"]
    return 2 * info["w_side"] + stringer + CROSS_W + deck_g + glue


def report(name, m, info, deck="厚紙 t=0.5"):
    r = evaluate(m)
    W = weight(info, DECK_OPTIONS[deck])
    P = r["gov"][1][0]
    c = r["caps"]
    print("\n■ %s" % name)
    print("  構造重量(両側) %.1f g ＋床版系ほか → 総重量 W ≈ %.1f g（%s）" % (2 * info["w_side"], W, deck))
    for key, label in (("member", "部材強度"), ("oop", "部材の面外座屈"), ("buckle", "面内全体座屈"),
                       ("joint", "交差接合"), ("hanger", "ハンガー")):
        v, why = c[key]
        if v != float("inf"):
            print("  %-10s P = %7.1f N  (%s)" % (label, v, why))
    print("  → 支配: %s  P ≈ %.0f N (%.1f kgf)" % (r["gov"][0], P, P / 9.81))
    print("  タイ張力 %.2f·P、載荷点たわみ %.3f mm/N" % (r["tie"], r["defl"]))
    print("  荷重/自重 = %.0f 倍   CP = P[N]/W[g]^1.5 = %.2f" % (P / 9.81 / (W / 1000), P / W ** 1.5))
    return dict(name=name, P=P, W=W, gov=r["gov"][0], caps=c, tie=r["tie"], defl=r["defl"])


def main():
    print("材料仮定:", MAT)
    res = []
    res.append(report("現案 ダビンチ（ナイロン3号ハンガー）", *davinci("ナイロン3号")))
    res.append(report("ダビンチ・PE1号ハンガー", *davinci("PE1号")))
    res.append(report("ダビンチ・バルサ3mmハンガー", *davinci("バルサ3mm角")))
    res.append(report("ダビンチ・同系端部を剛結（バルサハンガー）", *davinci("バルサ3mm角", rigid_ends=True)))
    res.append(report("ダビンチ・アンカー（タイなし）", *davinci("バルサ3mm角", tie=False, anchor=True)))
    res.append(report("通常タイドアーチ 5角1本・剛結", *polyarch("バルサ3mm角", s=5.0)))
    res.append(report("通常タイドアーチ 5角2本・剛結", *polyarch("バルサ3mm角", s=5.0, rows=2)))
    res.append(report("プラットトラス 弦5角・腹3角", *pratt()))
    res.append(report("プラットトラス 弦5角・腹4角", *pratt(s_web=4.0)))
    res.append(report("プラットトラス 弦4角・腹3角", *pratt(s_chord=4.0)))
    res.append(report("ダビンチ・アンカー＋縦桁3角（タイ省略）", *davinci("バルサ3mm角", tie=False, anchor=True, s_str=3.0)))
    res.append(report("通常タイドアーチ 5角2本・アンカー＋縦桁3角", *polyarch("バルサ3mm角", rows=2, tie=False,
                                                                 anchor=True, s_str=3.0)))
    print("\n=== 床版の重さを変えた場合")
    for deck in DECK_OPTIONS:
        for name, mi in (("ダビンチ", davinci("バルサ3mm角")), ("アーチ2本", polyarch("バルサ3mm角", rows=2)),
                         ("プラット5/4", pratt(s_web=4.0))):
            r = evaluate(mi[0])
            W = weight(mi[1], DECK_OPTIONS[deck])
            P = r["gov"][1][0]
            print("  %-22s %-8s W=%.1f g  P=%.0f N  荷重/自重=%.0f  CP=%.2f"
                  % (deck, name, W, P, P / 9.81 / W * 1000, P / W ** 1.5))
    print("\n=== 材料値の感度（支配荷重 P [N]）")
    base = dict(MAT)
    for key, vals in (("tau", (0.4, 1.0)), ("tau_lap", (0.5, 2.0)), ("fc90", (0.4, 1.0)),
                      ("E", (2100.0, 4300.0)), ("fc", (9.0, 14.9))):
        for v in vals:
            MAT[key] = v
            out = []
            for name, mi in (("ダビンチ", davinci("バルサ3mm角")), ("アーチ2本", polyarch("バルサ3mm角", rows=2)),
                             ("プラット5/4", pratt(s_web=4.0))):
                r = evaluate(mi[0])
                out.append("%s %.0f(%s)" % (name, r["gov"][1][0], r["gov"][0]))
            print("  %s=%g: %s" % (key, v, "  ".join(out)))
            MAT.update(base)


if __name__ == "__main__":
    main()
