#!/usr/bin/env python3
"""ダビンチアーチ橋（下路式・タイドアーチ、厚紙なし案）の簡易図面 SVG と確認用 STL を生成する。

標準ライブラリのみ使用。寸法の単位はすべて mm。
座標系: x = 橋軸方向（左支点 0 → 右支点 SPAN）、y = 鉛直上向き（支点頂部 0）、
        z = 幅員方向（橋軸中心 0）。STL では (x, z, y) に並べ替えて Z 上向きで出力する。

構成: 床版の両外側に立てたダビンチアーチのリブ2本から、載荷横桁を釣り糸で吊る。床版シートは使わない。
アーチの水平力は縦桁（タイ）で受け、三角支点には鉛直反力だけを伝える。
構成の根拠は docs/design-exploration.md（厚紙なし条件での探索の最良案）。

    python3 scripts/make_draft.py
"""
import math
import os
from collections import OrderedDict

# ---- 課題条件（実施要領） ----
SPAN = 500.0              # 支間（三角支点の頂部間）
MAX_MEMBER = 150.0        # 一部材の最大長
DECK_W_MIN = 80.0         # 床版有効幅の下限
CLEAR_MIN = 80.0          # 床版上有効高さの下限
LOAD_PITCH = 150.0        # 2点載荷の間隔（支間中央に対称と仮定）

# ---- 設計パラメータ ----
N_SEG = 10                # アーチ節点間の分割数（偶数: A系が両端まで通る）
T_LONG = 4.0              # アーチ縦材 断面（角材）
T_TRANS = 5.0             # アーチ横材 断面（角材）
OVERHANG = 4.0            # 縦材の端部節点からの突き出し（横材幅/2 + 1.5）
RIB_Z = 51.5              # リブ中心の幅方向位置 ±z
RIB_ROWS = (("A", 2.5), ("B", -2.5))  # リブ内の縦材配置（リブ中心から外側を正とする z）
RIB_T_HALF = 6.5          # リブ横材の半長（リブ中心から）
DECK_W = 90.0             # 有効幅（リブ内面間の純幅）
BRACE = 3.0               # リブ間つなぎ材 断面（角材）
STRINGER = 15.0           # 縦桁（タイ）の高さ（5×15 縦使い）
STRINGER_W = 5.0          # 縦桁（タイ）の幅
CROSS = 5.0               # 端横桁 断面（角材）
LOAD_BEAM = (5.0, 15.0)   # 載荷横桁（幅, 高さ: 5×15 縦使い）
HANGER = 0.5              # 釣り糸ハンガー（表示用の太さ）
RHO_BALSA = 0.16e-3       # g/mm3（バルサ密度の仮定値）
RHO_PAPER = 0.70e-3       # g/mm3（厚紙密度の仮定値）

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "work")
LOADS = (SPAN / 2 - LOAD_PITCH / 2, SPAN / 2 + LOAD_PITCH / 2)


# ---------------------------------------------------------------- 幾何
def solve_arch(n_seg, t_long, t_trans):
    """節点円の半径 R と節点角 α を求める。

    ダビンチアーチの自立条件: 縦材は両端で横材の上（半径 R+δ）、中央で隣の横材の下
    （半径 R-δ）を通る直線材なので cos α = (R-δ)/(R+δ)。これと支間条件
    n·α = 2·asin(SPAN/2/R) を同時に満たす R を二分法で解く。
    """
    delta = (t_long + t_trans) / 2
    half = SPAN / 2

    def f(r):
        return n_seg * math.acos((r - delta) / (r + delta)) - 2 * math.asin(half / r)

    lo, hi = half * 1.000001, half * 50
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    r = (lo + hi) / 2
    alpha = math.acos((r - delta) / (r + delta))
    chord = 2 * (r + delta) * math.sin(alpha)  # 縦材の端部節点間距離
    rise = r - math.sqrt(r * r - half * half)  # 節点円のライズ
    return r, alpha, delta, chord, rise


R, ALPHA, DELTA, CHORD, RISE = solve_arch(N_SEG, T_LONG, T_TRANS)
THETA = math.asin(SPAN / 2 / R)


def phi(i):
    return -THETA + i * ALPHA


def tbar_half_height(i):
    """節点 i の横材（断面を節点の半径方向に向けた角材）の鉛直方向の半高さ。"""
    a = phi(i)
    return T_TRANS / 2 * (abs(math.cos(a)) + abs(math.sin(a)))


# 端部の横材 T0 が床版縦桁の上面に載るように節点円の中心高さを決める
Y_C = STRINGER + tbar_half_height(0) - R * math.cos(THETA)


def pt(r, i):
    a = phi(i)
    return SPAN / 2 + r * math.sin(a), Y_C + r * math.cos(a)


def sub(a, b):
    return tuple(p - q for p, q in zip(a, b))


def add(a, b):
    return tuple(p + q for p, q in zip(a, b))


def mul(a, s):
    return tuple(p * s for p in a)


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return mul(a, 1 / math.sqrt(dot(a, a)))


class Box:
    """中心線 p0→p1、断面 sv×sw の直方体部材。v は断面の sv 方向。"""

    def __init__(self, mark, label, kind, material, p0, p1, v, sv, sw, rho):
        self.mark, self.label, self.kind, self.material = mark, label, kind, material
        self.p0, self.p1 = p0, p1
        u = norm(sub(p1, p0))
        v = norm(sub(v, mul(u, dot(v, u))))
        w = cross(u, v)
        self.length = math.sqrt(dot(sub(p1, p0), sub(p1, p0)))
        self.volume = self.length * sv * sw
        self.mass = self.volume * rho
        # corners[i]: bit0 = 端（p0/p1）, bit1 = ±v, bit2 = ±w
        self.corners = []
        for i in range(8):
            p = p1 if i & 1 else p0
            p = add(p, mul(v, sv / 2 if i & 2 else -sv / 2))
            p = add(p, mul(w, sw / 2 if i & 4 else -sw / 2))
            self.corners.append(p)

    def center(self):
        return mul(add(self.p0, self.p1), 0.5)

    def lo(self, k):
        return min(c[k] for c in self.corners)

    def hi(self, k):
        return max(c[k] for c in self.corners)

    def edges(self):
        for i in range(8):
            for bit in (1, 2, 4):
                if not i & bit:
                    yield self.corners[i], self.corners[i | bit]

    def triangles(self):
        c = self.corners
        ctr = self.center()
        tris = []
        for bit in (1, 2, 4):
            others = [b for b in (1, 2, 4) if b != bit]
            for side in (0, bit):
                quad = [side, side | others[0], side | others[0] | others[1], side | others[1]]
                for t in ((quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])):
                    a, b, d = c[t[0]], c[t[1]], c[t[2]]
                    n = cross(sub(b, a), sub(d, a))
                    if dot(n, sub(a, ctr)) < 0:  # 外向きにそろえる
                        b, d = d, b
                        n = mul(n, -1)
                    tris.append((norm(n), a, b, d))
        return tris


def balsa(s):
    return "バルサ {:g}mm角".format(s) if isinstance(s, float) else "バルサ {:g}×{:g}".format(*s)


def split(start, end, nodes):
    """start〜end を節点位置で継ぎ、各片を MAX_MEMBER 以下にする。"""
    cuts = [start]
    while end - cuts[-1] > MAX_MEMBER:
        cand = [x for x in nodes if cuts[-1] < x <= cuts[-1] + MAX_MEMBER]
        if not cand:
            raise SystemExit("%g〜%g を %g mm 以下に分割できません" % (start, end, MAX_MEMBER))
        cuts.append(max(cand))
    cuts.append(end)
    return cuts


def build():
    boxes = []
    xs = [pt(R, i)[0] for i in range(N_SEG + 1)]
    deck_top = STRINGER

    # ---- アーチリブ: 縦材は端部節点で横材の上、中央節点で横材の下を通る
    def longitudinal(mark, label, i0, i1, family):
        a, b = pt(R + DELTA, i0), pt(R + DELTA, i1)
        d = norm(sub(b, a))
        a = sub(a, mul(d, OVERHANG))
        b = add(b, mul(d, OVERHANG))
        for side in (-1, 1):
            for fam, dz in RIB_ROWS:
                if fam != family:
                    continue
                z = side * (RIB_Z + dz)
                boxes.append(Box(mark, label, mark, balsa(T_LONG), (a[0], a[1], z), (b[0], b[1], z),
                                 (-d[1], d[0], 0.0), T_LONG, T_LONG, RHO_BALSA))

    for n, k in enumerate(range(0, N_SEG, 2)):
        longitudinal("A", "A%d" % (n + 1), k, k + 2, "A")
    for n, k in enumerate(range(1, N_SEG - 1, 2)):
        longitudinal("B", "B%d" % (n + 1), k, k + 2, "B")
    longitudinal("BS", "BS", 0, 1, "B")          # B系の端部短材
    longitudinal("BS", "BS", N_SEG - 1, N_SEG, "B")

    # ---- 横材（各リブ内）。床版上の有効空間（高さ CLEAR_MIN）より上の節点だけ、リブ間をつなぎ材で結ぶ
    full_nodes = []
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        a = phi(i)
        v = (math.sin(a), math.cos(a), 0.0)
        for side in (-1, 1):
            zc = side * RIB_Z
            boxes.append(Box("T", "T%d" % i, "T", balsa(T_TRANS), (x, y, zc - RIB_T_HALF),
                             (x, y, zc + RIB_T_HALF), v, T_TRANS, T_TRANS, RHO_BALSA))
        if y - tbar_half_height(i) >= deck_top + CLEAR_MIN:
            full_nodes.append(i)
            zz = RIB_Z - RIB_T_HALF
            boxes.append(Box("TB", "TB%d" % i, "TB", balsa(BRACE), (x, y, -zz), (x, y, zz),
                             v, BRACE, BRACE, RHO_BALSA))

    # ---- 床版縦桁（タイ）: 下面 y=0 で三角支点に載る。節点位置で継ぐ
    rib = [b for b in boxes if b.kind in ("A", "B", "BS", "T", "TB")]
    s0, s1 = min(b.lo(0) for b in rib) - 1, max(b.hi(0) for b in rib) + 1
    cuts_s = split(s0, s1, xs[1:-1])
    for k in range(len(cuts_s) - 1):
        sym = min(k, len(cuts_s) - 2 - k) + 1
        for side in (-1, 1):
            boxes.append(Box("S%d" % sym, "S%d" % (k + 1), "S", balsa((STRINGER_W, STRINGER)),
                             (cuts_s[k], STRINGER / 2, side * RIB_Z), (cuts_s[k + 1], STRINGER / 2, side * RIB_Z),
                             (0.0, 1.0, 0.0), STRINGER, STRINGER_W, RHO_BALSA))

    # ---- 横桁: 縦桁の内側に突き付け、上面を縦桁上面にそろえる。床版シートは使わず、
    #      載荷治具は載荷点直下の載荷横桁（5×15 縦使い）に直接載る。両端は端横桁（5角）
    load_nodes = sorted(set(min(range(1, N_SEG), key=lambda i: abs(xs[i] - x)) for x in LOADS))
    zin = RIB_Z - STRINGER_W / 2
    for x, (w, h), mark, mat in ([(xs[0], (CROSS, CROSS), "CE", balsa(CROSS)),
                                  (xs[-1], (CROSS, CROSS), "CE", balsa(CROSS))]
                                 + [(x, LOAD_BEAM, "CL", balsa(LOAD_BEAM)) for x in LOADS]):
        yc = STRINGER - h / 2
        boxes.append(Box(mark, "C", "C", mat, (x, yc, -zin), (x, yc, zin), (0.0, 1.0, 0.0), h, w, RHO_BALSA))

    # ---- 釣り糸ハンガー: 荷重が入る節点の横材下面（リブ中心）→ 載荷横桁の端（縦桁上面）
    for i, x in zip(load_nodes, LOADS):
        top = pt(R, i)[1] - tbar_half_height(i)
        for side in (-1, 1):
            boxes.append(Box("H", "H%d" % i, "H", "釣り糸", (x, STRINGER, side * RIB_Z),
                             (xs[i], top, side * RIB_Z), (1.0, 0.0, 0.0), HANGER, HANGER, 0.0))

    lv = dict(deck_top=deck_top, cuts_s=cuts_s, cuts_d=[xs[0], xs[-1]], full_nodes=full_nodes,
              load_nodes=load_nodes, xs=xs)
    return boxes, lv


# ---------------------------------------------------------------- 照査
def checks(boxes, lv):
    res = []
    xs = lv["xs"]
    res.append(("支間 = %g mm（節点T0・T%d が支点直上）" % (SPAN, N_SEG), abs(xs[-1] - xs[0] - SPAN) < 1e-6,
                "%.1f mm" % (xs[-1] - xs[0])))
    solid = [b for b in boxes if b.kind != "H"]
    longest = max(solid, key=lambda b: b.length)
    res.append(("一部材の長さ ≤ %g mm（釣り糸を除く）" % MAX_MEMBER, longest.length <= MAX_MEMBER,
                "最大 %.1f mm（%s, 余裕 %.1f mm）" % (longest.length, longest.mark, MAX_MEMBER - longest.length)))
    inner = RIB_Z - RIB_T_HALF
    clear_w = min(DECK_W, 2 * inner)
    res.append(("床版有効幅 ≥ %g mm" % DECK_W_MIN, clear_w >= DECK_W_MIN,
                "リブ内側の純幅 %g mm（床版シートなし。載荷横桁が床版を兼ねる）" % (2 * inner)))
    # 床版上の有効空間（幅 DECK_W_MIN × 高さ CLEAR_MIN × 支間）に部材が入らないか
    zl, yl = DECK_W_MIN / 2, lv["deck_top"]
    hit = [b for b in boxes if b.kind != "D" and b.lo(2) < zl and b.hi(2) > -zl
           and b.lo(1) < yl + CLEAR_MIN and b.hi(1) > yl and b.lo(0) < SPAN and b.hi(0) > 0]
    over = [b.lo(1) - yl for b in boxes if b.mark == "TB"]
    res.append(("床版上有効高さ ≥ %g mm" % CLEAR_MIN, not hit,
                "つなぎ材の下端まで最小 %.1f mm（つなぎ材は T%s）"
                % (min(over), "・T".join(str(i) for i in lv["full_nodes"]))))
    near = [b for b in boxes if b.lo(0) < 15 and b.hi(0) > -15 or b.lo(0) < SPAN + 15 and b.hi(0) > SPAN - 15]
    low = min(b.lo(1) for b in near)
    res.append(("支点付近（±15 mm）で支点頂部より下に部材なし", low >= -1e-6,
                "最下面 y = %.1f mm（縦桁・端横桁の下面が支点頂部に線接触）" % low))
    beams = [b.p0[0] for b in boxes if b.mark == "CL"]
    on_beam = all(any(abs(x - bx) < 1e-6 for bx in beams) for x in LOADS)
    gaps = ["x=%g→節点T%d（x=%.1f）" % (x, i, lv["xs"][i]) for x, i in zip(LOADS, lv["load_nodes"])]
    res.append(("載荷点の直下に載荷横桁、ハンガーは直近節点へ", on_beam, "、".join(gaps)))
    return res


def bom(boxes):
    groups = OrderedDict()
    for b in boxes:
        groups.setdefault((b.mark, b.material, round(b.length, 1)), []).append(b)
    return [(mark, mat, length, len(items), sum(b.mass for b in items))
            for (mark, mat, length), items in groups.items()]


# ---------------------------------------------------------------- STL
def write_stl(boxes, path):
    lines = ["solid davinci_arch_bridge"]
    sw = lambda p: (p[0], p[2], p[1])  # (x, z, y) → Z 上向き
    for b in boxes:
        for n, *vs in b.triangles():
            vs = [sw(p) for p in vs]
            vs = [vs[0], vs[2], vs[1]]  # 軸の入替えで向きが反転するので頂点順も入れ替える
            lines.append("  facet normal %.6e %.6e %.6e" % sw(n))
            lines.append("    outer loop")
            for p in vs:
                lines.append("      vertex %.4f %.4f %.4f" % p)
            lines.append("    endloop")
            lines.append("  endfacet")
    lines.append("endsolid davinci_arch_bridge")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- SVG
def hull(pts):
    pts = sorted(set((round(x, 4), round(y, 4)) for x, y in pts))
    if len(pts) <= 2:
        return pts

    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cr(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cr(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def cut_x(box, x):
    """部材を平面 x = 一定 で切った断面の (z, y) 多角形。交差しなければ None。"""
    pts = []
    for a, b in box.edges():
        if (a[0] - x) * (b[0] - x) <= 0 and a[0] != b[0]:
            t = (x - a[0]) / (b[0] - a[0])
            p = add(a, mul(sub(b, a), t))
            pts.append((p[2], p[1]))
    return hull(pts) if len(pts) >= 3 else None


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Sheet:
    def __init__(self, w, h):
        self.w, self.h, self.items = w, h, []

    def add(self, s):
        self.items.append(s)

    def poly(self, pts, cls):
        self.add('<polygon class="%s" points="%s"/>' % (cls, " ".join("%.2f,%.2f" % p for p in pts)))

    def rect(self, x, y, w, h, cls):
        self.add('<rect class="%s" x="%.2f" y="%.2f" width="%.2f" height="%.2f"/>' % (cls, x, y, w, h))

    def line(self, x1, y1, x2, y2, cls="thin"):
        self.add('<line class="%s" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f"/>' % (cls, x1, y1, x2, y2))

    def text(self, x, y, s, size=6.0, anchor="start", cls="", rot=None):
        tr = ' transform="rotate(%g %.2f %.2f)"' % (rot, x, y) if rot else ""
        self.add('<text class="%s" x="%.2f" y="%.2f" font-size="%g" text-anchor="%s"%s>%s</text>'
                 % (cls, x, y, size, anchor, tr, esc(s)))

    def arrow(self, x1, y1, x2, y2, cls="load"):
        self.add('<line class="%s" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" marker-end="url(#ah)"/>'
                 % (cls, x1, y1, x2, y2))

    def dim_h(self, x1, x2, y, label, ext=None, size=5.5):
        if ext:
            for x, y0 in ((x1, ext[0]), (x2, ext[1])):
                self.line(x, y0, x, y + (2 if y > y0 else -2), "dim")
        self.line(x1, y, x2, y, "dim")
        for x in (x1, x2):
            self.line(x - 1.5, y + 1.5, x + 1.5, y - 1.5, "dimtick")
        self.text((x1 + x2) / 2, y - 1.5, label, size, "middle")

    def dim_v(self, x, y1, y2, label, ext=None, size=5.5):
        if ext:
            for y, x0 in ((y1, ext[0]), (y2, ext[1])):
                self.line(x0, y, x + (2 if x > x0 else -2), y, "dim")
        self.line(x, y1, x, y2, "dim")
        for y in (y1, y2):
            self.line(x - 1.5, y + 1.5, x + 1.5, y - 1.5, "dimtick")
        self.text(x - 1.5, (y1 + y2) / 2, label, size, "middle", rot=-90)

    def svg(self):
        style = """
  text { font-family: 'Hiragino Sans', 'Yu Gothic', 'Noto Sans JP', sans-serif; fill: #222; }
  .bg { fill: #fff; }
  .frame { fill: none; stroke: #222; stroke-width: 0.8; }
  .thin { stroke: #222; stroke-width: 0.3; fill: none; }
  .A  { fill: #3d74c9; fill-opacity: 0.55; stroke: #1d4a8f; stroke-width: 0.35; }
  .B  { fill: #d1493f; fill-opacity: 0.55; stroke: #8e231c; stroke-width: 0.35; }
  .BS { fill: #e8968f; fill-opacity: 0.65; stroke: #8e231c; stroke-width: 0.35; }
  .T  { fill: #f0b429; stroke: #8a6100; stroke-width: 0.3; }
  .TB { fill: #c98a00; stroke: #6b4a00; stroke-width: 0.3; }
  .H  { fill: #2e7d32; stroke: #2e7d32; stroke-width: 0.4; }
  .S  { fill: #b78a5b; fill-opacity: 0.8; stroke: #6b4a26; stroke-width: 0.3; }
  .C  { fill: #8d6e63; fill-opacity: 0.8; stroke: #4e342e; stroke-width: 0.3; }
  .D  { fill: #e9dcc3; stroke: #7a6a4a; stroke-width: 0.3; }
  .deckhidden { fill: none; stroke: #7a6a4a; stroke-width: 0.35; stroke-dasharray: 3 1.5; }
  .env { fill: #2e7d32; fill-opacity: 0.06; stroke: #2e7d32; stroke-width: 0.35; stroke-dasharray: 2 1.5; }
  .support { fill: #d9d9d9; stroke: #555; stroke-width: 0.4; }
  .cl { stroke: #555; stroke-width: 0.3; stroke-dasharray: 8 2 1.5 2; fill: none; }
  .dim { stroke: #444; stroke-width: 0.25; fill: none; }
  .dimtick { stroke: #444; stroke-width: 0.5; }
  .load { stroke: #c62828; stroke-width: 0.7; fill: none; }
  .loadtxt { fill: #c62828; } .envtxt { fill: #2e7d32; }
  .ok { fill: #2e7d32; } .ng { fill: #c62828; }
  .tblh { fill: #eee; stroke: #222; stroke-width: 0.3; }
"""
        head = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %g %g" width="%gmm" height="%gmm">\n'
                % (self.w, self.h, self.w, self.h))
        defs = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" '
                'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs>\n')
        return head + "<style>" + style + "</style>\n" + defs + "\n".join(self.items) + "\n</svg>\n"


ORDER = {"C": 1, "S": 2, "H": 3, "TB": 4, "T": 4, "B": 5, "BS": 5, "A": 6}
LEGEND = (("A", "A系縦材（偶数節点を結ぶ・リブ外側）"), ("B", "B系縦材（奇数節点を結ぶ・リブ内側）"),
          ("BS", "B系端部短材"), ("T", "横材（各リブ内）"), ("TB", "リブ間つなぎ材（有効高さより上）"),
          ("H", "釣り糸ハンガー（載荷節点のみ）"), ("S", "縦桁 5×15 縦使い＝タイ"),
          ("C", "横桁（載荷横桁5×15・端横桁5角）"))


def draw(boxes, lv, chk, rows):
    sh = Sheet(900, 920)
    sh.rect(0, 0, 900, 920, "bg")
    sh.rect(10, 10, 880, 900, "frame")
    sh.text(30, 32, "ダビンチアーチ橋（下路式・タイドアーチ、厚紙なし）一般図 — 簡易下書き", 11)
    sh.text(30, 43, "単位: mm ／ SVG座標 1単位 = 1mm（全図 1:1）／ AutoCAD 作図前の検討用", 5.5)
    draw_order = sorted(boxes, key=lambda b: ORDER[b.kind])
    deck_top = lv["deck_top"]
    crown = pt(R, N_SEG // 2)

    # ---- 側面図
    ox, oy = 80.0, 330.0
    S = lambda x, y: (ox + x, oy - y)
    sh.text(ox + SPAN / 2, 58, "側面図", 8, "middle")
    for x0 in (0.0, SPAN):
        sh.poly([S(x0, 0), S(x0 - 9, -16), S(x0 + 9, -16)], "support")
    sh.line(ox + SPAN / 2, oy + 25, ox + SPAN / 2, oy - crown[1] - 20, "cl")
    sh.text(ox + SPAN / 2 - 2, oy - crown[1] - 16, "CL", 5, "end")
    for b in draw_order:
        sh.poly([S(x, y) for x, y in hull([(c[0], c[1]) for c in b.corners])], b.kind)
    for b in boxes:   # 部材記号（手前のリブ・リブ中心側の列のみ）
        if b.kind in ("A", "B", "BS") and abs(b.p0[2] - (RIB_Z + (2.5 if b.kind == "A" else -2.5))) < 1e-6:
            mx, my = (b.p0[0] + b.p1[0]) / 2, (b.p0[1] + b.p1[1]) / 2
            p = S(mx, my + 9) if b.kind == "A" else S(mx, my - 12)
            sh.text(p[0], p[1], b.label, 4.5, "middle")
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        p = S(x, y)
        left = x < SPAN / 2
        sh.text(p[0] + (4 if left else -4), p[1] + 7, "T%d" % i, 3.8, "start" if left else "end")
    for x in LOADS:
        sh.arrow(ox + x, oy - deck_top - 26, ox + x, oy - deck_top - 0.8)
        sh.text(ox + x, oy - deck_top - 28, "P/2", 5, "middle", "loadtxt")
    sh.line(ox, oy - deck_top - CLEAR_MIN, ox + SPAN, oy - deck_top - CLEAR_MIN, "deckhidden")
    sh.text(ox + 4, oy - deck_top - CLEAR_MIN - 1.5, "床版上 有効高さ %g（幅 %g の範囲に部材なし）" % (CLEAR_MIN, DECK_W_MIN),
            4, cls="envtxt")
    yd = oy - crown[1] - 30
    sh.dim_h(ox, ox + LOADS[0], yd, "%g" % LOADS[0])
    sh.dim_h(ox + LOADS[0], ox + LOADS[1], yd, "%g（載荷間隔）" % LOAD_PITCH)
    sh.dim_h(ox + LOADS[1], ox + SPAN, yd, "%g" % (SPAN - LOADS[1]))
    sh.dim_h(ox, ox + SPAN, oy + 30, "支間 %g（三角支点の頂部間）" % SPAN, ext=(oy + 17, oy + 17))
    sh.dim_v(ox + SPAN + 40, oy, oy - crown[1], "節点ライズ %.1f" % crown[1], ext=(ox + SPAN + 10, ox + SPAN / 2 + 8))
    sh.dim_v(ox - 25, oy, oy - deck_top, "%.1f" % deck_top, ext=(ox - 12, ox - 12), size=4.2)
    sh.text(ox - 30, oy - 12, "横桁上面", 4.2, "end")
    sh.text(ox - 12, oy + 25, "支点は鉛直反力のみ（水平力は縦桁タイで処理）", 4.2, "start")
    for name, i in (("A", N_SEG // 2), ("B", 1)):
        x = ox + pt(R, i)[0] + (6 if name == "A" else 3)
        for yy in (oy - crown[1] - 12, oy + 8):
            sh.line(x, yy, x, yy + 6, "frame")
        sh.text(x + 1.5, oy - crown[1] - 13, name, 5)
        sh.text(x + 1.5, oy + 18, name, 5)

    # ---- 平面図
    px, py = ox, 470.0
    P = lambda x, z: (px + x, py - z)
    sh.text(px + SPAN / 2, py - 78, "平面図（床版シートなし・破線＝有効幅 %g）" % DECK_W, 8, "middle")
    sh.line(px - 15, py, px + SPAN + 15, py, "cl")
    for b in draw_order:
        if True:
            sh.poly([P(x, z) for x, z in hull([(c[0], c[2]) for c in b.corners])], b.kind)
    d0, d1 = lv["cuts_d"][0], lv["cuts_d"][-1]
    sh.rect(px + d0, py - DECK_W / 2, d1 - d0, DECK_W, "deckhidden")
    for x in LOADS:
        sh.line(px + x, py - DECK_W / 2, px + x, py + DECK_W / 2, "load")
        sh.text(px + x, py - 2, "載荷線", 4.2, "middle", "loadtxt")
    zt = RIB_Z + RIB_T_HALF
    zs = [zt] + [RIB_Z + dz for _, dz in sorted(RIB_ROWS, key=lambda r: -r[1])] + [RIB_Z - RIB_T_HALF,
                                                                                -(RIB_Z - RIB_T_HALF)]
    for a, b in zip(zs, zs[1:]):
        sh.dim_v(px + SPAN + 28, py - a, py - b, "%g" % (a - b), size=4.2)
    sh.dim_v(px + SPAN + 42, py - DECK_W / 2, py + DECK_W / 2, "有効幅 %g（≧%g）" % (DECK_W, DECK_W_MIN),
             size=4.2)
    sh.dim_v(px + SPAN + 56, py - zt, py + zt, "全幅 %g" % (2 * zt), size=4.2)
    for z, name in ((RIB_Z, "リブ"), (0.0, "床版")):
        sh.text(px - 12, py - z + 1.5, name, 4, "end")
    xs = lv["xs"]
    for a, b in zip(xs, xs[1:]):
        sh.dim_h(px + a, px + b, py + zt + 14, "%.1f" % (b - a), ext=(py + zt + 2, py + zt + 2), size=4.2)
    sh.text(px + SPAN / 2, py + zt + 24, "節点（横材・ハンガー・横桁）の水平間隔", 4.5, "middle")

    # ---- 断面図（1:1、平面 x = 節点位置 で切断）
    def section(cx, top_y, i, title):
        k = 1.0
        xcut = xs[i] + 0.01
        polys = [(b.kind, cut_x(b, xcut)) for b in draw_order]
        polys = [(kd, pg) for kd, pg in polys if pg]
        y_hi = max(max(p[1] for p in pg) for _, pg in polys)
        y_hi = max(y_hi, deck_top + CLEAR_MIN)
        Q = lambda z, y: (cx + k * z, top_y + k * (y_hi - y))
        a, b = Q(-DECK_W_MIN / 2, deck_top + CLEAR_MIN), Q(DECK_W_MIN / 2, deck_top)
        sh.rect(a[0], a[1], b[0] - a[0], b[1] - a[1], "env")
        sh.text(cx, a[1] + 6, "有効空間 %g×%g" % (DECK_W_MIN, CLEAR_MIN), 4, "middle", "envtxt")
        for kd, pg in polys:
            sh.poly([Q(z, y) for z, y in pg], kd)
        y_lo = min(min(p[1] for p in pg) for _, pg in polys)
        sh.line(cx, top_y - 3, cx, Q(0, y_lo)[1] + 4, "cl")
        sh.text(cx, top_y - 12, title, 7, "middle")
        zt_ = RIB_Z + RIB_T_HALF
        sh.dim_h(Q(-zt_, 0)[0], Q(zt_, 0)[0], Q(0, y_lo)[1] + 10, "全幅 %g" % (2 * zt_), size=4.2)
        return Q(0, y_lo)[1] + 10

    section(760, 85, 1, "断面B-B（節点T1）")
    section(760, 230, N_SEG // 2, "断面A-A（クラウンT%d）" % (N_SEG // 2))
    lx, ly = 660, 560
    for n, (kind, name) in enumerate(LEGEND):
        sh.rect(lx, ly + n * 8, 8, 5, kind)
        sh.text(lx + 11, ly + n * 8 + 4.3, name, 4.4)

    # ---- 照査・幾何
    tx, ty = 30, 600
    sh.text(tx, ty, "ルール照査（自動計算）", 7)
    for n, (name, ok, val) in enumerate(chk):
        y = ty + 10 + n * 15
        sh.text(tx, y, "OK" if ok else "NG", 5, cls="ok" if ok else "ng")
        sh.text(tx + 12, y, name, 4.6)
        sh.text(tx + 12, y + 6, val, 4.2)
    gy = ty + 10 + len(chk) * 15 + 6
    sh.text(tx, gy, "幾何（ダビンチアーチの自立条件 cosα=(R−δ)/(R+δ)）", 7)
    geo = ["分割数 n = %d、δ = (縦材+横材)/2 = %.1f mm" % (N_SEG, DELTA),
           "節点円半径 R = %.1f mm、節点角 α = %.2f°" % (R, math.degrees(ALPHA)),
           "縦材の節点間距離 = %.1f mm（＋突き出し %g×2）" % (CHORD, OVERHANG),
           "節点ライズ = %.1f mm、床版上面 = %.1f mm（支点頂部から）" % (crown[1], deck_top)]
    for n, s in enumerate(geo):
        sh.text(tx, gy + 10 + n * 8, s, 4.6)

    # ---- 部材表
    bx, by = 400, 600
    sh.text(bx, by, "部材表（重量は密度仮定による概算）", 7)
    cols = [("記号", 20), ("材料", 52), ("長さ", 26), ("本数", 20), ("重量g", 26)]
    y = by + 5
    x = bx
    for name, w in cols:
        sh.rect(x, y, w, 8, "tblh")
        sh.text(x + w / 2, y + 5.8, name, 4.4, "middle")
        x += w
    total = 0.0
    for mark, mat, length, count, mass in rows:
        y += 8
        x = bx
        for (name, w), val in zip(cols, (mark, mat, "%.1f" % length, "%d" % count, "%.2f" % mass)):
            sh.rect(x, y, w, 8, "thin")
            sh.text(x + w / 2, y + 5.8, val, 4.2, "middle")
            x += w
        total += mass
    y += 8
    sh.text(bx, y + 7, "合計（接着剤を除く）≈ %.1f g" % total, 5.5)
    sh.text(bx, y + 14, "密度仮定: バルサ %.2f g/cm³" % (RHO_BALSA * 1e3), 4.2)

    # ---- 表題欄
    tbx, tby = 700, 860
    for n, (k_, v_) in enumerate((("図名", "一般図（簡易下書き）"), ("構造形式", "ダビンチアーチ・下路式タイドアーチ（厚紙なし）"),
                                  ("作成日", "2026-10-21"))):
        sh.rect(tbx, tby + n * 13, 40, 13, "thin")
        sh.rect(tbx + 40, tby + n * 13, 140, 13, "thin")
        sh.text(tbx + 20, tby + n * 13 + 8.5, k_, 5, "middle")
        sh.text(tbx + 45, tby + n * 13 + 8.5, v_, 5)
    return sh.svg(), total


def compare():
    out = []
    for n_seg in (8, 10):
        for tl, tt in ((5, 3), (4, 3), (5, 4), (4, 4), (5, 5)):
            r, a, d, c, f = solve_arch(n_seg, tl, tt)
            out.append((n_seg, tl, tt, r, f, c + 2 * (tt / 2 + 1.5)))
    return out


def main():
    boxes, lv = build()
    chk = checks(boxes, lv)
    rows = bom(boxes)
    os.makedirs(OUT_DIR, exist_ok=True)
    svg, total = draw(boxes, lv, chk, rows)
    with open(os.path.join(OUT_DIR, "bridge-draft.svg"), "w") as f:
        f.write(svg)
    write_stl(boxes, os.path.join(OUT_DIR, "bridge-draft.stl"))

    print("R = %.2f mm, alpha = %.3f deg, delta = %.2f mm" % (R, math.degrees(ALPHA), DELTA))
    print("crown node y = %.2f mm, deck top = %.2f mm" % (pt(R, N_SEG // 2)[1], lv["deck_top"]))
    print("full-width transverse nodes:", lv["full_nodes"], " load beams at nodes:", lv["load_nodes"])
    print("stringer cuts:", ", ".join("%.1f" % x for x in lv["cuts_s"]))
    print("deck cuts:", ", ".join("%.1f" % x for x in lv["cuts_d"]))
    print("\n[checks]")
    for name, ok, val in chk:
        print("  %s  %s : %s" % ("OK" if ok else "NG", name, val))
    print("\n[BOM] mark | material | length | count | mass g")
    for r_ in rows:
        print("  %-3s | %s | %6.1f | %2d | %.2f" % r_)
    print("  total = %.1f g (members: %d)" % (total, len(boxes)))
    print("\n[compare] n | long | trans | R | node rise | member length")
    for n_seg, tl, tt, r, f, ln in compare():
        print("  %2d | %g | %g | %.1f | %.1f | %.1f %s" % (n_seg, tl, tt, r, f + tt / 2, ln,
                                                        "" if ln <= MAX_MEMBER else "NG(>150)"))


if __name__ == "__main__":
    main()
