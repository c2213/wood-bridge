#!/usr/bin/env python3
"""ダビンチアーチ＋ワーレン腹材（床版なし）の簡易図面 SVG と確認用 STL を生成する。

標準ライブラリのみ使用。寸法の単位はすべて mm。
座標系: x = 橋軸方向（左支点 0 → 右支点 SPAN）、y = 鉛直上向き（支点頂部 0）、
        z = 幅員方向（橋軸中心 0）。STL では (x, z, y) に並べ替えて Z 上向きで出力する。

構成: ダビンチアーチ（A・B 系の縦材を横材で編む）を上弦材とし、支点高さのタイ（3mm 角 2 本）と
ワーレン腹材（タイと同じ両側の 2 面）で結ぶ。床版はない。錘はアーチの節点 T3・T5 の横材から吊る。
水平力はタイで受け、三角支点には鉛直反力だけを伝える。
構成の根拠は docs/design-exploration.md（床版なし・2点荷重向け補強の探索）。

    python3 scripts/make_draft.py
"""
import math
import os
from collections import OrderedDict

# ---- 課題条件（実施要領） ----
SPAN = 500.0              # 支間（三角支点の頂部間）
MAX_MEMBER = 150.0        # 一部材の最大長
DECK_W_MIN = 80.0         # 床版有効幅の下限（本案は床版なし。探索スクリプトが参照）
CLEAR_MIN = 80.0          # 床版上有効高さの下限（同上）
LOAD_PITCH = 150.0        # 2点載荷の間隔（支間中央に対称と仮定）
DECK_W = 90.0             # 旧案の有効幅（探索スクリプトが参照）

# ---- 設計パラメータ ----
N_SEG = 8                 # アーチ節点間の分割数（偶数: A系が両端まで通る）
T_LONG = 3.0              # アーチ縦材 断面（角材）
T_TRANS = 5.0             # アーチ横材 断面（角材）
OVERHANG = 4.0            # 支点側の端部での縦材の突き出し（横材の半幅 + 1.5）。中間節点は突き付け
ROWS = (("A", -6.0), ("B", -2.0), ("B", 2.0), ("A", 6.0))  # 縦材の列（z）。A2 列・B2 列
TRANS_LEN = 22.0          # 横材の長さ（縦材の列幅 15 + 両側に腹材 3 + 余裕 0.5）
TIE = 3.0                 # タイ 断面（角材）
TIE_Z = 9.0               # タイと腹材の面の位置 ±z（縦材の列の外側、横材の端の下）
WEB = 3.0                 # ワーレン腹材 断面（角材）
HANGER = 0.5              # 錘の吊り糸（表示用の太さ）
RHO_BALSA = 0.16e-3       # g/mm3（バルサ密度の仮定値）
RHO_PAPER = 0.70e-3       # g/mm3（厚紙。探索スクリプトの比較用）

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


# 端部の横材 T0 がタイの上面に載るように節点円の中心高さを決める
Y_C = TIE + tbar_half_height(0) - R * math.cos(THETA)


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
    return "バルサ {:g}mm角".format(s)


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

    # ---- アーチ縦材: 端部節点で横材の上、中央節点で横材の下を通る。中間節点では同系の端を突き付ける
    def longitudinal(mark, label, i0, i1, family):
        a, b = pt(R + DELTA, i0), pt(R + DELTA, i1)
        d = norm(sub(b, a))
        if i0 == 0:
            a = sub(a, mul(d, OVERHANG))
        if i1 == N_SEG:
            b = add(b, mul(d, OVERHANG))
        for fam, z in ROWS:
            if fam == family:
                boxes.append(Box(mark, label, mark, balsa(T_LONG), (a[0], a[1], z), (b[0], b[1], z),
                                 (-d[1], d[0], 0.0), T_LONG, T_LONG, RHO_BALSA))

    for n, k in enumerate(range(0, N_SEG, 2)):
        mark = "A" if 0 < k < N_SEG - 2 else "AE"
        longitudinal(mark, "A%d" % (n + 1), k, k + 2, "A")
    for n, k in enumerate(range(1, N_SEG - 1, 2)):
        longitudinal("B", "B%d" % (n + 1), k, k + 2, "B")
    longitudinal("BS", "BS", 0, 1, "B")          # B系の端部短材
    longitudinal("BS", "BS", N_SEG - 1, N_SEG, "B")

    # ---- 横材（全節点、幅 TRANS_LEN）
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        a = phi(i)
        boxes.append(Box("T", "T%d" % i, "T", balsa(T_TRANS), (x, y, -TRANS_LEN / 2), (x, y, TRANS_LEN / 2),
                         (math.sin(a), math.cos(a), 0.0), T_TRANS, T_TRANS, RHO_BALSA))

    # ---- タイ（2 本、下面 y=0 で三角支点に載る）。節点位置で継ぐ
    arch = [b for b in boxes if b.kind in ("A", "AE", "B", "BS", "T")]
    s0, s1 = min(b.lo(0) for b in arch) - 1, max(b.hi(0) for b in arch) + 1
    cuts = split(s0, s1, xs[1:-1])
    for k in range(len(cuts) - 1):
        sym = min(k, len(cuts) - 2 - k) + 1
        for side in (-1, 1):
            boxes.append(Box("S%d" % sym, "S%d" % (k + 1), "S", balsa(TIE),
                             (cuts[k], TIE / 2, side * TIE_Z), (cuts[k + 1], TIE / 2, side * TIE_Z),
                             (0.0, 1.0, 0.0), TIE, TIE, RHO_BALSA))

    # ---- ワーレン腹材（タイ上面 ↔ 横材の下面、両側の面）
    webs = []
    for i in range(N_SEG):
        lo_i, hi_i = (i, i + 1) if i % 2 == 0 else (i + 1, i)   # (タイ側の節点, アーチ側の節点)
        top = (xs[hi_i], pt(R, hi_i)[1] - tbar_half_height(hi_i))
        bot = (xs[lo_i], TIE)
        webs.append((lo_i, hi_i))
        for side in (-1, 1):
            d = norm(sub(top, bot))
            boxes.append(Box("W%d" % (min(i, N_SEG - 1 - i) + 1), "W%d" % (i + 1), "W", balsa(WEB),
                             (bot[0], bot[1], side * TIE_Z), (top[0], top[1], side * TIE_Z),
                             (-d[1], d[0], 0.0), WEB, WEB, RHO_BALSA))

    # ---- 錘の吊り糸（表示用）: 載荷節点の横材の両端から下へ
    load_nodes = sorted(set(min(range(1, N_SEG), key=lambda i: abs(xs[i] - x)) for x in LOADS))
    for i in load_nodes:
        x, y = pt(R, i)
        for side in (-1, 1):
            boxes.append(Box("H", "H%d" % i, "H", "吊り糸（試験用）", (x, y - tbar_half_height(i), side * 10.5),
                             (x, -30.0, side * 10.5), (1.0, 0.0, 0.0), HANGER, HANGER, 0.0))

    lv = dict(cuts=cuts, load_nodes=load_nodes, xs=xs, webs=webs)
    return boxes, lv


# ---------------------------------------------------------------- 照査
def checks(boxes, lv):
    res = []
    xs = lv["xs"]
    res.append(("支間 = %g mm（節点T0・T%d が支点直上）" % (SPAN, N_SEG), abs(xs[-1] - xs[0] - SPAN) < 1e-6,
                "%.1f mm" % (xs[-1] - xs[0])))
    solid = [b for b in boxes if b.kind != "H"]
    longest = max(solid, key=lambda b: b.length)
    res.append(("一部材の長さ ≤ %g mm（吊り糸を除く）" % MAX_MEMBER, longest.length <= MAX_MEMBER,
                "最大 %.1f mm（%s, 余裕 %.1f mm）" % (longest.length, longest.mark, MAX_MEMBER - longest.length)))
    res.append(("床版の有効幅・有効高さ", None, "床版なし（ユーザー判断）。規定の扱いは担当教員に要確認"))
    near = [b for b in boxes if b.kind != "H" and (b.lo(0) < 15 and b.hi(0) > -15
                                                   or b.lo(0) < SPAN + 15 and b.hi(0) > SPAN - 15)]
    low = min(b.lo(1) for b in near)
    res.append(("支点付近（±15 mm）で支点頂部より下に部材なし", low >= -1e-6,
                "最下面 y = %.1f mm（タイの下面が支点頂部に線接触）" % low))
    gaps = ["x=%g→T%d（x=%.1f）" % (x, i, xs[i]) for x, i in zip(LOADS, lv["load_nodes"])]
    off = max(abs(x - xs[i]) for x, i in zip(LOADS, lv["load_nodes"]))
    res.append(("錘を吊る節点と載荷点のずれ ≤ 5 mm", off <= 5.0, "、".join(gaps)))
    return res


def bom(boxes):
    groups = OrderedDict()
    for b in boxes:
        if b.kind == "H":
            continue
        groups.setdefault((b.mark, b.material, round(b.length, 1)), []).append(b)
    return [(mark, mat, length, len(items), sum(b.mass for b in items))
            for (mark, mat, length), items in groups.items()]


# ---------------------------------------------------------------- STL
def write_stl(boxes, path):
    lines = ["solid davinci_arch_bridge"]
    sw = lambda p: (p[0], p[2], p[1])  # (x, z, y) → Z 上向き
    for b in boxes:
        if b.kind == "H":
            continue
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
  .A, .AE { fill: #3d74c9; fill-opacity: 0.55; stroke: #1d4a8f; stroke-width: 0.35; }
  .B  { fill: #d1493f; fill-opacity: 0.55; stroke: #8e231c; stroke-width: 0.35; }
  .BS { fill: #e8968f; fill-opacity: 0.65; stroke: #8e231c; stroke-width: 0.35; }
  .T  { fill: #f0b429; stroke: #8a6100; stroke-width: 0.3; }
  .W  { fill: #2e9d8f; fill-opacity: 0.6; stroke: #1b5e55; stroke-width: 0.3; }
  .H  { fill: #c62828; stroke: #c62828; stroke-width: 0.3; stroke-dasharray: 2 1; }
  .S  { fill: #b78a5b; fill-opacity: 0.85; stroke: #6b4a26; stroke-width: 0.3; }
  .support { fill: #d9d9d9; stroke: #555; stroke-width: 0.4; }
  .cl { stroke: #555; stroke-width: 0.3; stroke-dasharray: 8 2 1.5 2; fill: none; }
  .dim { stroke: #444; stroke-width: 0.25; fill: none; }
  .dimtick { stroke: #444; stroke-width: 0.5; }
  .load { stroke: #c62828; stroke-width: 0.7; fill: none; }
  .loadtxt { fill: #c62828; }
  .ok { fill: #2e7d32; } .ng { fill: #c62828; } .chk { fill: #b26a00; }
  .tblh { fill: #eee; stroke: #222; stroke-width: 0.3; }
"""
        head = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %g %g" width="%gmm" height="%gmm">\n'
                % (self.w, self.h, self.w, self.h))
        defs = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" '
                'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs>\n')
        return head + "<style>" + style + "</style>\n" + defs + "\n".join(self.items) + "\n</svg>\n"


ORDER = {"S": 1, "W": 2, "H": 3, "T": 4, "B": 5, "BS": 5, "A": 6, "AE": 6}
LEGEND = (("A", "A系縦材（偶数節点を結ぶ・外側2列）"), ("B", "B系縦材（奇数節点を結ぶ・内側2列）"),
          ("BS", "B系端部短材"), ("T", "横材（全節点）"), ("W", "ワーレン腹材（両側の2面）"),
          ("S", "タイ（2本・アーチの水平力を受ける）"), ("H", "錘の吊り糸（試験用・重量外）"))


def draw(boxes, lv, chk, rows):
    W_, H_ = 900, 800
    sh = Sheet(W_, H_)
    sh.rect(0, 0, W_, H_, "bg")
    sh.rect(10, 10, W_ - 20, H_ - 20, "frame")
    sh.text(30, 32, "ダビンチアーチ＋ワーレン腹材（床版なし）一般図 — 簡易下書き", 11)
    sh.text(30, 43, "単位: mm ／ SVG座標 1単位 = 1mm（側面・平面 1:1、断面 3:1）／ AutoCAD 作図前の検討用", 5.5)
    draw_order = sorted(boxes, key=lambda b: ORDER[b.kind])
    crown = pt(R, N_SEG // 2)
    xs = lv["xs"]

    # ---- 側面図
    ox, oy = 80.0, 245.0
    S = lambda x, y: (ox + x, oy - y)
    sh.text(ox + SPAN / 2, 58, "側面図", 8, "middle")
    for x0 in (0.0, SPAN):
        sh.poly([S(x0, 0), S(x0 - 9, -16), S(x0 + 9, -16)], "support")
    sh.line(ox + SPAN / 2, oy + 25, ox + SPAN / 2, oy - crown[1] - 20, "cl")
    sh.text(ox + SPAN / 2 - 2, oy - crown[1] - 16, "CL", 5, "end")
    for b in draw_order:
        if b.kind == "H":
            continue
        sh.poly([S(x, y) for x, y in hull([(c[0], c[1]) for c in b.corners])], b.kind)
    for b in boxes:   # 部材記号（手前の列のみ）
        if b.kind in ("A", "AE", "B", "BS") and b.p0[2] in (6.0, 2.0):
            mx, my = (b.p0[0] + b.p1[0]) / 2, (b.p0[1] + b.p1[1]) / 2
            p = S(mx, my + 8) if b.kind in ("A", "AE") else S(mx, my - 10)
            sh.text(p[0], p[1], b.label, 4.5, "middle")
        if b.kind == "W" and b.p0[2] > 0:
            mx, my = (b.p0[0] + b.p1[0]) / 2, (b.p0[1] + b.p1[1]) / 2
            sh.text(*S(mx + 3, my), b.label, 3.8)
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        p = S(x, y)
        left = x < SPAN / 2
        sh.text(p[0] + (5 if left else -5), p[1] - 4, "T%d" % i, 3.8, "start" if left else "end")
    for i, x in zip(lv["load_nodes"], LOADS):
        xn, yn = pt(R, i)
        top = yn - tbar_half_height(i)
        sh.line(ox + xn, oy - top, ox + xn, oy + 8, "load")
        sh.arrow(ox + xn, oy + 8, ox + xn, oy + 20)
        sh.text(ox + xn + 2, oy + 16, "P/2（T%dから吊る）" % i, 4.5, cls="loadtxt")
    yd = oy - crown[1] - 30
    sh.dim_h(ox, ox + LOADS[0], yd, "%g" % LOADS[0])
    sh.dim_h(ox + LOADS[0], ox + LOADS[1], yd, "%g（載荷間隔）" % LOAD_PITCH)
    sh.dim_h(ox + LOADS[1], ox + SPAN, yd, "%g" % (SPAN - LOADS[1]))
    sh.dim_h(ox, ox + SPAN, oy + 36, "支間 %g（三角支点の頂部間）" % SPAN, ext=(oy + 17, oy + 17))
    sh.dim_v(ox + SPAN + 40, oy, oy - crown[1], "節点ライズ %.1f" % crown[1], ext=(ox + SPAN + 10, ox + SPAN / 2 + 8))
    sh.text(ox - 12, oy + 48, "支点は鉛直反力のみ（水平力はタイで処理）", 4.2, "start")
    for name, i in (("A", N_SEG // 2), ("B", 1)):
        x = ox + pt(R, i)[0] + (6 if name == "A" else 3)
        for yy in (oy - crown[1] - 12, oy + 22):
            sh.line(x, yy, x, yy + 6, "frame")
        sh.text(x + 1.5, oy - crown[1] - 13, name, 5)
        sh.text(x + 1.5, oy + 32, name, 5)

    # ---- 平面図
    px, py = ox, 345.0
    P = lambda x, z: (px + x, py - z * 2)
    sh.text(px + SPAN / 2, py - 32, "平面図（幅方向のみ 2 倍に拡大）", 8, "middle")
    sh.line(px - 15, py, px + SPAN + 15, py, "cl")
    for b in draw_order:
        if b.kind != "H":
            sh.poly([P(x, z) for x, z in hull([(c[0], c[2]) for c in b.corners])], b.kind)
    zs = [TRANS_LEN / 2, 6 + T_LONG / 2, -(6 + T_LONG / 2), -TRANS_LEN / 2]
    for a, b in zip(zs, zs[1:]):
        sh.dim_v(px + SPAN + 26, py - 2 * a, py - 2 * b, "%g" % (a - b), size=4.2)
    for a, b in zip(xs, xs[1:]):
        sh.dim_h(px + a, px + b, py + 38, "%.1f" % (b - a), ext=(py + 24, py + 24), size=4.2)
    sh.text(px + SPAN / 2, py + 48, "節点（横材）の水平間隔", 4.5, "middle")

    # ---- 断面図（3 倍、平面 x = 節点位置 で切断）
    def section(cx, top_y, i, title):
        k = 3.0
        xcut = xs[i] + 0.01
        polys = [(b.kind, cut_x(b, xcut)) for b in draw_order if b.kind != "H"]
        polys = [(kd, pg) for kd, pg in polys if pg]
        ys = [p[1] for _, pg in polys for p in pg]
        # タイ付近とアーチ付近の 2 か所だけを描く（間は省略）
        bands = [(-1, TIE + 6), (pt(R, i)[1] - 10, pt(R, i)[1] + 8)]
        y_cur = top_y
        for lo_, hi_ in reversed(bands):
            Q = lambda z, y: (cx + k * z, y_cur + k * (hi_ - y))
            for kd, pg in polys:
                if all(lo_ <= p[1] <= hi_ for p in pg):
                    sh.poly([Q(z, y) for z, y in pg], kd)
                elif any(lo_ <= p[1] <= hi_ for p in pg):   # 帯をまたぐ部材（腹材）は帯で切って描く
                    clip = [(z, min(max(y, lo_), hi_)) for z, y in pg]
                    sh.poly([Q(z, y) for z, y in clip], kd)
            sh.text(cx + k * 12, y_cur + k * (hi_ - lo_) / 2, "y=%.0f〜%.0f" % (lo_, hi_), 4)
            y_cur += k * (hi_ - lo_) + 8
            if lo_ > 0:
                sh.text(cx, y_cur - 3, "（中間省略）", 4, "middle")
        sh.text(cx, top_y - 8, title, 7, "middle")
        sh.dim_h(cx - k * TRANS_LEN / 2, cx + k * TRANS_LEN / 2, y_cur + 6, "横材 %g" % TRANS_LEN, size=4.2)
        return y_cur

    section(715, 80, N_SEG // 2, "断面A-A（クラウンT%d）" % (N_SEG // 2))
    section(715, 260, 1, "断面B-B（節点T1）")
    lx, ly = 640, 440
    for n, (kind, name) in enumerate(LEGEND):
        sh.rect(lx, ly + n * 8, 8, 5, kind)
        sh.text(lx + 11, ly + n * 8 + 4.3, name, 4.4)

    # ---- 照査・幾何
    tx, ty = 30, 470
    sh.text(tx, ty, "ルール照査（自動計算）", 7)
    for n, (name, ok, val) in enumerate(chk):
        y = ty + 10 + n * 15
        lab, cls = ("要確認", "chk") if ok is None else (("OK", "ok") if ok else ("NG", "ng"))
        sh.text(tx, y, lab, 5, cls=cls)
        sh.text(tx + 20, y, name, 4.6)
        sh.text(tx + 20, y + 6, val, 4.2)
    gy = ty + 10 + len(chk) * 15 + 6
    sh.text(tx, gy, "幾何（ダビンチアーチの自立条件 cosα=(R−δ)/(R+δ)）", 7)
    geo = ["分割数 n = %d、δ = (縦材+横材)/2 = %.1f mm" % (N_SEG, DELTA),
           "節点円半径 R = %.1f mm、節点角 α = %.2f°" % (R, math.degrees(ALPHA)),
           "縦材の節点間距離 = %.1f mm（中間節点は突き付け、支点側のみ突き出し %g）" % (CHORD, OVERHANG),
           "クラウン節点高さ = %.1f mm（支点頂部から）" % crown[1]]
    for n, s in enumerate(geo):
        sh.text(tx, gy + 10 + n * 8, s, 4.6)

    # ---- 部材表
    bx, by = 345, 470
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
    tbx, tby = 700, H_ - 60
    for n, (k_, v_) in enumerate((("図名", "一般図（簡易下書き）"), ("構造形式", "ダビンチアーチ＋ワーレン腹材（床版なし）"),
                                  ("作成日", "2026-10-21"))):
        sh.rect(tbx, tby + n * 13, 40, 13, "thin")
        sh.rect(tbx + 40, tby + n * 13, 140, 13, "thin")
        sh.text(tbx + 20, tby + n * 13 + 8.5, k_, 5, "middle")
        sh.text(tbx + 45, tby + n * 13 + 8.5, v_, 5)
    return sh.svg(), total


def compare():
    out = []
    for n_seg in (8, 10):
        for tl, tt in ((3, 5), (3, 4), (4, 4), (4, 5), (5, 5)):
            r, a, d, c, f = solve_arch(n_seg, tl, tt)
            out.append((n_seg, tl, tt, r, f, c + (tt / 2 + 1.5)))
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

    print("R = %.2f mm, alpha = %.3f deg, delta = %.2f mm, chord = %.2f mm" % (R, math.degrees(ALPHA), DELTA, CHORD))
    print("crown node y = %.2f mm" % pt(R, N_SEG // 2)[1])
    print("load nodes:", lv["load_nodes"], " tie cuts:", ", ".join("%.1f" % x for x in lv["cuts"]))
    print("\n[checks]")
    for name, ok, val in chk:
        print("  %s  %s : %s" % ("要確認" if ok is None else "OK" if ok else "NG", name, val))
    print("\n[BOM] mark | material | length | count | mass g")
    for r_ in rows:
        print("  %-3s | %s | %6.1f | %2d | %.2f" % r_)
    print("  total = %.1f g (members: %d)" % (total, len([b for b in boxes if b.kind != "H"])))
    print("\n[compare] n | long | trans | R | node rise | member length (butt joints)")
    for n_seg, tl, tt, r, f, ln in compare():
        print("  %2d | %g | %g | %.1f | %.1f | %.1f %s" % (n_seg, tl, tt, r, f + tt / 2, ln,
                                                        "" if ln <= MAX_MEMBER else "NG(>150)"))


if __name__ == "__main__":
    main()
