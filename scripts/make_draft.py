#!/usr/bin/env python3
"""ダビンチアーチ橋（上路式）の簡易図面 SVG と確認用 STL を生成する。

標準ライブラリのみ使用。寸法の単位はすべて mm。
座標系: x = 橋軸方向（左支点 0 → 右支点 SPAN）、y = 鉛直上向き（支点上面 0）、
        z = 幅員方向（橋軸中心 0）。STL では (x, z, y) に並べ替えて Z 上向きで出力する。

    python3 scripts/make_draft.py
"""
import math
import os
from collections import OrderedDict

# ---- 課題条件（実施要領） ----
SPAN = 500.0              # 支間
MAX_MEMBER = 150.0        # 一部材の最大長
DECK_W_MIN = 80.0         # 床版有効幅の下限
CLEAR_MIN = 80.0          # 床版上有効高さの下限
LOAD_PITCH = 150.0        # 2点載荷の間隔（支間中央に対称と仮定）

# ---- 設計パラメータ ----
N_SEG = 8                 # アーチ節点間の分割数（偶数: A系が両端まで通る）
T_LONG = 5.0              # アーチ縦材 断面（角材）
T_TRANS = 3.0             # アーチ横材 断面（角材）
OVERHANG = 3.0            # 縦材の端部節点からの突き出し
ROWS_A = (-30.0, 30.0)    # A系縦材の幅方向位置 z
ROWS_B = (-10.0, 10.0)    # B系縦材の幅方向位置 z
TRANS_LEN = 90.0          # 横材の長さ
DECK_W = 90.0             # 床版（厚紙）幅
DECK_T = 0.5              # 床版（厚紙）厚
STRINGER = 5.0            # 床版縦桁 断面（角材）
STRINGER_Z = 40.0         # 床版縦桁の幅方向位置 ±z
POST = 4.0                # 鉛直材 断面（角材）
RHO_BALSA = 0.16e-3       # g/mm3（バルサ密度の仮定値）
RHO_PAPER = 0.70e-3       # g/mm3（厚紙密度の仮定値）

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "work")


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
Y_C = -R * math.cos(THETA) + T_TRANS / 2  # 節点円の中心 y（横材T0の下面を支点上面 y=0 に合わせる）


def phi(i):
    return -THETA + i * ALPHA


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

    def ymax(self):
        return max(c[1] for c in self.corners)

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


def build():
    boxes = []
    balsa = lambda s: "バルサ {:g}mm角".format(s)

    # アーチ縦材: 端部節点では横材の上（R+δ）、中央節点では横材の下を通る
    def longitudinal(mark, label, i0, i1, rows):
        a, b = pt(R + DELTA, i0), pt(R + DELTA, i1)
        d = norm(sub(b, a))
        a = sub(a, mul(d, OVERHANG))
        b = add(b, mul(d, OVERHANG))
        for z in rows:
            boxes.append(Box(mark, label, mark, balsa(T_LONG), (a[0], a[1], z), (b[0], b[1], z),
                             (-d[1], d[0], 0.0), T_LONG, T_LONG, RHO_BALSA))

    for n, k in enumerate(range(0, N_SEG, 2)):
        longitudinal("A", "A%d" % (n + 1), k, k + 2, ROWS_A)
    for n, k in enumerate(range(1, N_SEG - 1, 2)):
        longitudinal("B", "B%d" % (n + 1), k, k + 2, ROWS_B)
    longitudinal("BS", "BS", 0, 1, ROWS_B)          # B系の端部短材
    longitudinal("BS", "BS", N_SEG - 1, N_SEG, ROWS_B)

    # アーチ横材
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        a = phi(i)
        boxes.append(Box("T", "T%d" % i, "T", balsa(T_TRANS), (x, y, -TRANS_LEN / 2), (x, y, TRANS_LEN / 2),
                         (math.sin(a), math.cos(a), 0.0), T_TRANS, T_TRANS, RHO_BALSA))

    arch_top = max(b.ymax() for b in boxes if b.kind in ("A", "B", "BS"))
    trans_top = max(b.ymax() for b in boxes if b.kind == "T")
    plate_bot = max(arch_top, trans_top + STRINGER)
    str_bot = plate_bot - STRINGER

    # 鉛直材（横材の端部上 → 床版縦桁下面）
    for i in range(N_SEG + 1):
        tb = [b for b in boxes if b.label == "T%d" % i][0]
        bottom = tb.ymax()
        if str_bot - bottom < 2.0:
            continue  # 縦桁が横材に直接（または厚紙スペーサを介して）載る
        x = pt(R, i)[0]
        for z in (-STRINGER_Z, STRINGER_Z):
            boxes.append(Box("P%d" % min(i, N_SEG - i), "P%d" % i, "P", balsa(POST),
                             (x, bottom, z), (x, str_bot, z), (1.0, 0.0, 0.0), POST, POST, RHO_BALSA))

    # 床版縦桁・床版（厚紙）: 節点位置で継ぎ、各片を MAX_MEMBER 以下にする
    xs = [pt(R, i)[0] for i in range(N_SEG + 1)]
    start, end = xs[0] - POST / 2, xs[-1] + POST / 2
    cuts = [start]
    while end - cuts[-1] > MAX_MEMBER:
        cand = [x for x in xs[1:-1] if cuts[-1] < x <= cuts[-1] + MAX_MEMBER]
        if not cand:
            raise SystemExit("床版縦桁を %g mm 以下に分割できません" % MAX_MEMBER)
        cuts.append(max(cand))
    cuts.append(end)
    n_piece = len(cuts) - 1
    for k in range(n_piece):
        sym = min(k, n_piece - 1 - k) + 1
        x0, x1 = cuts[k], cuts[k + 1]
        for z in (-STRINGER_Z, STRINGER_Z):
            boxes.append(Box("S%d" % sym, "S%d" % (k + 1), "S", balsa(STRINGER),
                             (x0, str_bot + STRINGER / 2, z), (x1, str_bot + STRINGER / 2, z),
                             (0.0, 1.0, 0.0), STRINGER, STRINGER, RHO_BALSA))
        boxes.append(Box("D%d" % sym, "D%d" % (k + 1), "D", "厚紙 t={:g}".format(DECK_T),
                         (x0, plate_bot + DECK_T / 2, 0.0), (x1, plate_bot + DECK_T / 2, 0.0),
                         (0.0, 1.0, 0.0), DECK_T, DECK_W, RHO_PAPER))

    levels = dict(arch_top=arch_top, trans_top=trans_top, plate_bot=plate_bot,
                  str_bot=str_bot, deck_top=plate_bot + DECK_T, cuts=cuts)
    return boxes, levels


# ---------------------------------------------------------------- 照査
def checks(boxes, lv):
    res = []
    xs = [pt(R, i)[0] for i in range(N_SEG + 1)]
    res.append(("支間 = %g mm" % SPAN, abs(xs[-1] - xs[0] - SPAN) < 1e-6, "%.1f mm" % (xs[-1] - xs[0])))
    longest = max(boxes, key=lambda b: b.length)
    res.append(("一部材の長さ ≤ %g mm" % MAX_MEMBER, longest.length <= MAX_MEMBER,
                "最大 %.1f mm（%s, 余裕 %.1f mm）" % (longest.length, longest.mark, MAX_MEMBER - longest.length)))
    res.append(("床版有効幅 ≥ %g mm" % DECK_W_MIN, DECK_W >= DECK_W_MIN, "%g mm" % DECK_W))
    above = [b for b in boxes if b.kind != "D" and min(c[1] for c in b.corners) > lv["deck_top"] - 1e-6]
    res.append(("床版上有効高さ ≥ %g mm" % CLEAR_MIN, not above, "上路式のため床版上に部材なし"))
    loads = (SPAN / 2 - LOAD_PITCH / 2, SPAN / 2 + LOAD_PITCH / 2)
    on_deck = all(lv["cuts"][0] <= x <= lv["cuts"][-1] for x in loads)
    res.append(("載荷点が床版上", on_deck, "x = %g, %g mm" % loads))
    # 中央（z=0）から錘を吊る場合に吊り糸が横材と干渉しないか
    slot = min(abs(z) for z in ROWS_A + ROWS_B) * 2 - T_LONG
    gaps = []
    for x in loads:
        g = []
        for b in boxes:
            if b.kind != "T":
                continue
            lo, hi = min(c[0] for c in b.corners), max(c[0] for c in b.corners)
            g.append(0.0 if lo <= x <= hi else min(abs(x - lo), abs(x - hi)))
        gaps.append(min(g))
    res.append(("中央吊り（z=0）時の吊り糸と横材の離れ", min(gaps) >= 2.0,
                "最小 %.1f mm（縦材間スリット幅 %.0f mm）" % (min(gaps), slot)))
    side = TRANS_LEN / 2
    res.append(("側方吊り（|z| > %g）時の干渉" % side, max(DECK_W, TRANS_LEN) / 2 <= side + 1e-6,
                "床版・横材とも |z| ≤ %g mm" % side))
    return res, loads


def bom(boxes):
    groups = OrderedDict()
    for b in boxes:
        groups.setdefault((b.mark, b.material, round(b.length, 1)), []).append(b)
    rows = []
    for (mark, mat, length), items in groups.items():
        rows.append((mark, mat, length, len(items), sum(b.mass for b in items)))
    return rows


# ---------------------------------------------------------------- STL
def write_stl(boxes, path):
    lines = ["solid davinci_arch_bridge"]
    for b in boxes:
        for n, *vs in b.triangles():
            sw = lambda p: (p[0], p[2], p[1])  # (x, z, y) → Z 上向き
            n = sw(n)
            vs = [sw(p) for p in vs]
            # 軸の入替えで向きが反転するので頂点順も入れ替える
            vs = [vs[0], vs[2], vs[1]]
            lines.append("  facet normal %.6e %.6e %.6e" % n)
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
  .P  { fill: #9aa0a6; fill-opacity: 0.7; stroke: #4a4f54; stroke-width: 0.3; }
  .S  { fill: #b78a5b; fill-opacity: 0.75; stroke: #6b4a26; stroke-width: 0.3; }
  .D  { fill: #e9dcc3; stroke: #7a6a4a; stroke-width: 0.3; }
  .deckhidden { fill: none; stroke: #7a6a4a; stroke-width: 0.35; stroke-dasharray: 3 1.5; }
  .ground { fill: #d9d9d9; stroke: #555; stroke-width: 0.4; }
  .cl { stroke: #555; stroke-width: 0.3; stroke-dasharray: 8 2 1.5 2; fill: none; }
  .dim { stroke: #444; stroke-width: 0.25; fill: none; }
  .dimtick { stroke: #444; stroke-width: 0.5; }
  .load { stroke: #c62828; stroke-width: 0.7; fill: none; }
  .thrust { stroke: #2e7d32; stroke-width: 0.7; fill: none; }
  .loadtxt { fill: #c62828; }
  .ok { fill: #2e7d32; } .ng { fill: #c62828; }
  .tblh { fill: #eee; stroke: #222; stroke-width: 0.3; }
"""
        head = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %g %g" width="%gmm" height="%gmm">\n'
                % (self.w, self.h, self.w, self.h))
        defs = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" '
                'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs>\n')
        return head + "<style>" + style + "</style>\n" + defs + "\n".join(self.items) + "\n</svg>\n"


ORDER = {"D": 0, "S": 1, "P": 2, "T": 3, "B": 4, "BS": 4, "A": 5}


def draw(boxes, lv, chk, loads, rows):
    sh = Sheet(900, 800)
    sh.rect(0, 0, 900, 800, "bg")
    sh.rect(10, 10, 880, 780, "frame")
    sh.text(30, 32, "ダビンチアーチ橋（上路式）一般図 — 簡易下書き", 11)
    sh.text(30, 43, "単位: mm ／ SVG座標 1単位 = 1mm（断面図のみ2倍）／ AutoCAD 作図前の検討用", 5.5)

    draw_order = sorted(boxes, key=lambda b: ORDER[b.kind])

    # ---- 側面図
    ox, oy = 80.0, 245.0
    S = lambda x, y: (ox + x, oy - y)
    sh.text(ox + SPAN / 2, 62, "側面図", 8, "middle")
    for x0 in (-45.0, SPAN):
        sh.rect(ox + x0, oy, 45, 25, "ground")
    sh.line(ox + SPAN / 2, oy + 30, ox + SPAN / 2, oy - lv["deck_top"] - 30, "cl")
    sh.text(ox + SPAN / 2 + 2, oy - lv["deck_top"] - 26, "CL", 5)
    for b in draw_order:
        sh.poly([S(x, y) for x, y in hull([(c[0], c[1]) for c in b.corners])], b.kind)
    # 部材記号
    for b in boxes:
        if b.kind in ("A", "B", "BS") and b.p0[2] > 0:
            mx, my = (b.p0[0] + b.p1[0]) / 2, (b.p0[1] + b.p1[1]) / 2
            a = math.atan2(mx - SPAN / 2, my - Y_C)
            off = 10 if b.kind == "A" else -9
            p = S(mx + off * math.sin(a) * 0.4, my + off * math.cos(a))
            if b.kind != "A":
                p = S(mx, my - 12)
            sh.text(p[0], p[1], b.label, 4.5, "middle", "")
    for i in range(N_SEG + 1):
        x, y = pt(R, i)
        p = S(x, y)
        sh.text(p[0] + (4 if x < SPAN / 2 else -4), p[1] + 7, "T%d" % i, 3.8,
                "start" if x < SPAN / 2 else "end")
    # 載荷
    for x in loads:
        sh.arrow(ox + x, oy - lv["deck_top"] - 22, ox + x, oy - lv["deck_top"] - 0.8)
        sh.text(ox + x, oy - lv["deck_top"] - 24, "P/2", 5, "middle", "loadtxt")
    yd = oy - lv["deck_top"] - 32
    sh.dim_h(ox, ox + loads[0], yd, "%g" % loads[0])
    sh.dim_h(ox + loads[0], ox + loads[1], yd, "%g（載荷間隔）" % LOAD_PITCH)
    sh.dim_h(ox + loads[1], ox + SPAN, yd, "%g" % (SPAN - loads[1]))
    # 水平反力（スパン外アンカー）
    sh.add('<line class="thrust" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" marker-end="url(#ah)"/>'
           % (ox - 40, oy - 4, ox - 6, oy - 4))
    sh.add('<line class="thrust" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" marker-end="url(#ah)"/>'
           % (ox + SPAN + 40, oy - 4, ox + SPAN + 6, oy - 4))
    sh.text(ox - 43, oy + 31, "→ 水平反力はスパン外アンカーで受ける ←", 4.2)
    # 寸法
    sh.dim_h(ox, ox + SPAN, oy + 38, "支間 %g" % SPAN, ext=(oy + 2, oy + 2))
    sh.dim_v(ox - 52, oy, oy - lv["deck_top"], "床版上面 %.1f" % lv["deck_top"], ext=(ox - 45, ox - 2))
    crown = pt(R, N_SEG // 2)
    sh.dim_v(ox + SPAN + 55, oy, oy - crown[1], "節点ライズ %.1f" % crown[1],
             ext=(ox + SPAN + 45, ox + SPAN / 2 + 8))
    # 断面位置
    for name, i in (("A", N_SEG // 2), ("B", 1)):
        x = ox + pt(R, i)[0] + (6 if name == "A" else 0)
        for yy in (oy - lv["deck_top"] - 6, oy + 10):
            sh.line(x, yy, x, yy + 6, "frame")
        sh.text(x + 1.5, oy - lv["deck_top"] - 7, name, 5)
        sh.text(x + 1.5, oy + 22, name, 5)

    # ---- 平面図（床版を除いて表示）
    px, py = ox, 395.0
    P = lambda x, z: (px + x, py - z)
    sh.text(px + SPAN / 2, py - 62, "平面図（床版を取り除いて表示・破線＝床版外形）", 8, "middle")
    sh.line(px - 15, py, px + SPAN + 15, py, "cl")
    for b in draw_order:
        if b.kind == "D":
            continue
        sh.poly([P(x, z) for x, z in hull([(c[0], c[2]) for c in b.corners])], b.kind)
    c0, c1 = lv["cuts"][0], lv["cuts"][-1]
    sh.rect(px + c0, py - DECK_W / 2, c1 - c0, DECK_W, "deckhidden")
    for x in loads:
        sh.line(px + x, py - DECK_W / 2 - 8, px + x, py + DECK_W / 2 + 4, "load")
        sh.text(px + x, py - DECK_W / 2 - 10, "載荷線", 4.5, "middle", "loadtxt")
    # 幅方向寸法
    zs = [TRANS_LEN / 2, ROWS_A[1], ROWS_B[1], ROWS_B[0], ROWS_A[0], -TRANS_LEN / 2]
    for a, b in zip(zs, zs[1:]):
        sh.dim_v(px + SPAN + 30, py - a, py - b, "%g" % (a - b), size=4.5)
    sh.dim_v(px + SPAN + 45, py - DECK_W / 2, py + DECK_W / 2, "床版幅 %g（有効幅≧%g）" % (DECK_W, DECK_W_MIN),
             size=4.5)
    for z, name in ((ROWS_A[1], "A系"), (ROWS_B[1], "B系"), (STRINGER_Z, "縦桁・鉛直材")):
        sh.text(px - 10, py - z + 1.5, name, 4, "end")
    # 節点間隔
    xs = [pt(R, i)[0] for i in range(N_SEG + 1)]
    for a, b in zip(xs, xs[1:]):
        sh.dim_h(px + a, px + b, py + DECK_W / 2 + 16, "%.1f" % (b - a), ext=(py + 46, py + 46), size=4.2)
    sh.text(px + SPAN / 2, py + DECK_W / 2 + 26, "節点（横材）の水平間隔", 4.5, "middle")
    # 床版継手
    for x in lv["cuts"][1:-1]:
        sh.text(px + x, py + DECK_W / 2 + 4.5, "継手", 3.6, "middle")

    # ---- 断面図（2倍）
    def section(cx, base_y, i, title):
        k = 2.0
        ys = []
        rects = []  # (kind, z0, z1, y0, y1)
        rects.append(("D", -DECK_W / 2, DECK_W / 2, lv["plate_bot"], lv["deck_top"]))
        for z in (-STRINGER_Z, STRINGER_Z):
            rects.append(("S", z - STRINGER / 2, z + STRINGER / 2, lv["str_bot"], lv["plate_bot"]))
        tc = pt(R, i)[1]
        rects.append(("T", -TRANS_LEN / 2, TRANS_LEN / 2, tc - T_TRANS / 2, tc + T_TRANS / 2))
        ends, mids = (ROWS_A, ROWS_B) if i % 2 == 0 else (ROWS_B, ROWS_A)
        ye = pt(R + DELTA, i)[1]
        for z in ends:
            rects.append(("A" if ends is ROWS_A else "B", z - T_LONG / 2, z + T_LONG / 2,
                          ye - T_LONG / 2, ye + T_LONG / 2))
        if 0 < i < N_SEG:
            ym = pt(R - DELTA, i)[1]
            for z in mids:
                rects.append(("A" if mids is ROWS_A else "B", z - T_LONG / 2, z + T_LONG / 2,
                              ym - T_LONG / 2, ym + T_LONG / 2))
        posts = [b for b in boxes if b.label == "P%d" % i]
        for b in posts:
            z = b.p0[2]
            rects.append(("P", z - POST / 2, z + POST / 2, b.p0[1], b.p1[1]))
        y_lo = min(r[3] for r in rects) - 2
        Q = lambda z, y: (cx + k * z, base_y - k * (y - y_lo))
        for kind, z0, z1, y0, y1 in sorted(rects, key=lambda r: ORDER[r[0]]):
            a, b = Q(z0, y1), Q(z1, y0)
            sh.rect(a[0], a[1], b[0] - a[0], b[1] - a[1], kind)
        top = Q(0, lv["deck_top"])[1]
        sh.text(cx, top - 22, title, 7, "middle")
        sh.dim_h(cx - k * DECK_W / 2, cx + k * DECK_W / 2, top - 6, "床版幅 %g" % DECK_W, size=4.5)
        sh.line(cx, top - 2, cx, base_y + 4, "cl")
        lbl = "A系" if ends is ROWS_A else "B系"
        p = Q(TRANS_LEN / 2 + 1, ye)
        sh.text(p[0], p[1] - 3, "%s端部（横材の上）" % lbl, 3.8)
        if 0 < i < N_SEG:
            lbl2 = "B系" if ends is ROWS_A else "A系"
            p = Q(TRANS_LEN / 2 + 1, ym)
            sh.text(p[0], p[1] + 5, "%s中央部（横材の下）" % lbl2, 3.8)
        p = Q(TRANS_LEN / 2 + 1, tc)
        sh.text(p[0], p[1] + 1.5, "横材T%d" % i, 3.8)
        return base_y

    section(735, 235, 1, "断面B-B（節点T1・2倍）")
    section(735, 330, N_SEG // 2, "断面A-A（クラウンT%d・2倍）" % (N_SEG // 2))
    # 凡例
    lx, ly = 650, 360
    for n, (kind, name) in enumerate((("A", "A系縦材（偶数節点を結ぶ）"), ("B", "B系縦材（奇数節点を結ぶ）"),
                                      ("BS", "B系端部短材"), ("T", "横材"), ("P", "鉛直材"),
                                      ("S", "床版縦桁"), ("D", "床版（厚紙）"))):
        sh.rect(lx, ly + n * 8, 8, 5, kind)
        sh.text(lx + 11, ly + n * 8 + 4.3, name, 4.5)

    # ---- 照査
    tx, ty = 30, 505
    sh.text(tx, ty, "ルール照査（自動計算）", 7)
    for n, (name, ok, val) in enumerate(chk):
        y = ty + 10 + n * 8
        sh.text(tx, y, "OK" if ok else "NG", 5, cls="ok" if ok else "ng")
        sh.text(tx + 12, y, "%s：%s" % (name, val), 4.6)
    gy = ty + 10 + len(chk) * 8 + 8
    sh.text(tx, gy, "幾何（ダビンチアーチの自立条件 cosα=(R−δ)/(R+δ)）", 7)
    geo = ["分割数 n = %d、δ = (縦材+横材)/2 = %.1f mm" % (N_SEG, DELTA),
           "節点円半径 R = %.1f mm、節点角 α = %.2f°" % (R, math.degrees(ALPHA)),
           "縦材の節点間距離 = %.1f mm（＋突き出し %g×2）" % (CHORD, OVERHANG),
           "節点ライズ = %.1f mm、床版上面高さ = %.1f mm" % (RISE + T_TRANS / 2, lv["deck_top"])]
    for n, s in enumerate(geo):
        sh.text(tx, gy + 10 + n * 8, s, 4.6)

    # ---- 部材表
    bx, by = 470, 505
    sh.text(bx, by, "部材表（重量は密度仮定による概算）", 7)
    cols = [("記号", 22), ("材料", 62), ("長さ", 30), ("本数", 22), ("重量g", 30)]
    y = by + 5
    x = bx
    for name, w in cols:
        sh.rect(x, y, w, 8, "tblh")
        sh.text(x + w / 2, y + 5.8, name, 4.6, "middle")
        x += w
    total = 0.0
    for mark, mat, length, count, mass in rows:
        y += 8
        x = bx
        for (name, w), val in zip(cols, (mark, mat, "%.1f" % length, "%d" % count, "%.2f" % mass)):
            sh.rect(x, y, w, 8, "thin")
            sh.text(x + w / 2, y + 5.8, val, 4.6, "middle")
            x += w
        total += mass
    y += 8
    sh.text(bx, y + 7, "合計（接着剤を除く）≈ %.1f g" % total, 5.5)
    sh.text(bx, y + 15, "密度仮定: バルサ %.2f g/cm³、厚紙 %.2f g/cm³" % (RHO_BALSA * 1e3, RHO_PAPER * 1e3), 4.4)

    # ---- 表題欄
    tbx, tby = 700, 735
    for n, (k_, v_) in enumerate((("図名", "一般図（簡易下書き）"), ("構造形式", "ダビンチアーチ・上路式"),
                                  ("縮尺", "1:1（断面図 2:1）"), ("作成日", "2026-10-21"))):
        sh.rect(tbx, tby + n * 11, 40, 11, "thin")
        sh.rect(tbx + 40, tby + n * 11, 140, 11, "thin")
        sh.text(tbx + 20, tby + n * 11 + 7.5, k_, 5, "middle")
        sh.text(tbx + 45, tby + n * 11 + 7.5, v_, 5)
    return sh.svg(), total


def compare():
    out = []
    for n_seg in (8, 10):
        for tl, tt in ((5, 3), (4, 3), (5, 4), (4, 4), (5, 5)):
            r, a, d, c, f = solve_arch(n_seg, tl, tt)
            out.append((n_seg, tl, tt, r, f, c + 2 * OVERHANG))
    return out


def main():
    boxes, lv = build()
    chk, loads = checks(boxes, lv)
    rows = bom(boxes)
    os.makedirs(OUT_DIR, exist_ok=True)
    svg, total = draw(boxes, lv, chk, loads, rows)
    with open(os.path.join(OUT_DIR, "bridge-draft.svg"), "w") as f:
        f.write(svg)
    write_stl(boxes, os.path.join(OUT_DIR, "bridge-draft.stl"))

    print("R = %.2f mm, alpha = %.3f deg, delta = %.2f mm" % (R, math.degrees(ALPHA), DELTA))
    print("node rise = %.2f mm, deck top = %.2f mm" % (RISE + T_TRANS / 2, lv["deck_top"]))
    print("deck/stringer cuts:", ", ".join("%.1f" % x for x in lv["cuts"]))
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
