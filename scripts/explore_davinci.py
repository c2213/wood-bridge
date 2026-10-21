#!/usr/bin/env python3
"""ダビンチアーチを維持したまま、CP = P / W^1.5 が最大になる構成を総当たりで探す。

FEM・照査・材料値は audit_fem.py と共通。標準ライブラリのみ。
    python3 scripts/explore_davinci.py            # 全探索して上位を表示
"""
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit_fem as af  # noqa: E402
import make_draft as md  # noqa: E402

MAT = af.MAT
# 断面: (幅 b [z 方向 or 弧に沿う方向], 高さ h [半径方向])
SECTIONS = {"3角": (3.0, 3.0), "4角": (4.0, 4.0), "5角": (5.0, 5.0), "5×15平": (15.0, 5.0)}
ROWS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 3)]
SHEETS = {"なし": 0.0, "ケント紙0.2": 6.4, "厚紙0.5": 15.9}
GAP = 1.0           # リブ内の列間のすき間
GUSSET = dict(L=10.0, t=0.5)   # 交差部の厚紙ガセット（両側面、重ね長さ L、厚さ t）
CLEAR_EDGE = 45.0   # 下路: リブ内面の z（有効幅 80 + 余裕）


def geometry(n_seg, lng, trs, trs_load=None):
    r, alpha, delta, chord, rise = md.solve_arch(n_seg, lng[1], trs[1])
    over = max(trs[0], trs_load[0] if trs_load else 0.0) / 2 + 1.5
    theta = math.asin(md.SPAN / 2 / r)
    feasible = abs(n_seg * alpha - 2 * theta) < 1e-6   # 自立条件と支間条件を同時に満たす R が存在するか
    return dict(r=r, alpha=alpha, delta=delta, chord=chord, rise=rise, over=over,
                length=chord + 2 * over, theta=theta, feasible=feasible)


def sub2(a, b):
    return a[0] - b[0], a[1] - b[1]


def ab_joint_cap(m, ab, b, interfaces):
    """A–B 側面接着点で受け渡す力（P=1 あたり）から、その接着が決める破壊荷重を返す。"""
    worst = 0.0
    for nd, sin_t in ab:
        fx = fy = 0.0
        for e in m.elems:
            if e.get("kind") != "A" or nd not in (e["i"], e["j"]):
                continue
            kl, T = m.k_local(e)
            dg = [m.u[d] for d in e["dofs"]]
            dl = [sum(T[a][c] * dg[c] for c in range(6)) for a in range(6)]
            fl = [sum(kl[a][c] * dl[c] for c in range(6)) for a in range(6)]
            fg = [sum(T[c][a] * fl[c] for c in range(6)) for a in range(6)]
            o = 0 if e["i"] == nd else 3
            fx += fg[o]
            fy += fg[o + 1]
        cap = interfaces * (b * b / sin_t) * MAT["tau_lap"]   # 交差する2材の側面の重なり面積
        worst = max(worst, math.hypot(fx, fy) / cap)
    return 1 / worst if worst > 0 else float("inf")


def build(cfg):
    """cfg から 2D モデルと重量を作る。下路は片側リブ（荷重 P/2）、上路は全幅（荷重 P）。"""
    n, lng, trs, (nA, nB) = cfg["n"], SECTIONS[cfg["long"]], SECTIONS[cfg["trans"]], cfg["rows"]
    tl = SECTIONS[cfg["trans_load"]] if cfg.get("trans_load") else None
    g = geometry(n, lng, trs, tl)
    if cfg.get("over") is not None:   # 突き出しを指定（横材の半幅 + 0.5 mm 以上が必要）
        g["length"] += 2 * (cfg["over"] - g["over"])
        g["over"] = cfg["over"]
        if cfg["over"] < max(trs[0], tl[0] if tl else 0.0) / 2 + 0.5:
            g["feasible"] = False
    r, alpha, delta, theta = g["r"], g["alpha"], g["delta"], g["theta"]
    E, rho = MAT["E"], MAT["rho"]
    s_str = cfg.get("tie", 5.0)
    h_str = cfg.get("tie_h", s_str)

    def ph(i):
        return -theta + i * alpha

    hh = lambda i: trs[1] / 2 * abs(math.cos(ph(i))) + trs[0] / 2 * abs(math.sin(ph(i)))
    yc = h_str + hh(0) - r * math.cos(theta)

    def P(rr, i):
        return md.SPAN / 2 + rr * math.sin(ph(i)), yc + rr * math.cos(ph(i))

    m = af.Model()
    xs = [P(r, i)[0] for i in range(n + 1)]
    through = cfg["type"] == "下路"
    deck = af.add_deck(m, xs, s_str=s_str, h_str=h_str)
    if not through:   # add_deck は 1 本分の剛性。上路（全幅モデル）はタイ 2 本
        for e in m.elems:
            if e.get("kind") == "S":
                e["EA"] *= 2
                e["EI"] *= 2
                e["rows"] = 2
    tn = [m.node(*P(r, i)) for i in range(n + 1)]
    b, h = lng
    A1, I1 = b * h, b * h ** 3 / 12
    area = b * trs[0]

    jcap = area * MAT["tau"]
    if cfg["joint"] == "厚紙ガセット":   # 部材と横材の側面どうしを厚紙で両面から接ぐ（弱い側の面で決まる）
        jcap += 2 * min(h, trs[1]) * GUSSET["L"] * MAT["tau_lap"]

    load_nodes = sorted(set(min(range(1, n), key=lambda k: abs(P(r, k)[0] - x)) for x in md.LOADS))
    wide = set(load_nodes) if tl else set()   # 荷重が入る節点だけ横材を幅広（5×15 平使い）にする

    def arm(a, t, rows):
        ar, jc = (b * tl[0], jcap * tl[0] / trs[0]) if t in [tn[i] for i in wide] else (area, jcap)
        if cfg.get("soft_joints"):
            # 交差部の実際の柔らかさ: 法線方向は横材の繊維直角圧縮（E90）、接線方向は横材の
            # ころがりせん断（G_roll）が、厚さ δ の層で変形するとみなす。rows 列ぶんを並列に足す。
            L = math.dist(m.nodes[a], m.nodes[t])
            kn = MAT["E90"] * ar / delta * rows
            kt = MAT["G_roll"] * ar / delta * rows
            m.frame(a, t, kn * L, kt * L ** 3 / 3, hj=True, kind="J", rows=rows, s=b, area=ar, jcap=jc)
        else:
            m.frame(a, t, E * A1 * 1e3, E * I1 * 1e3, hj=True, kind="J", rows=rows, s=b, area=ar, jcap=jc)

    mids = {}
    members = [("A", k, k + 2) for k in range(0, n, 2)] + [("B", k, k + 2) for k in range(1, n - 1, 2)] \
        + [("B", 0, 1), ("B", n - 1, n)]
    length_total = 0.0
    ends = {}
    mems = []
    for fam, i0, i1 in members:
        nr = nA if fam == "A" else nB
        pts = []   # (部材上の位置 t, ノード)
        for t_, i in ((0.0, i0), (1.0, i1)):
            if cfg.get("rigid_ends") and (fam, i) in ends:   # 同系の端どうしを突き付け接着（連続材として扱う）
                pts.append((t_, ends[(fam, i)]))
                continue
            p = m.node(*P(r + delta, i))
            arm(p, tn[i], nr)
            ends[(fam, i)] = p
            pts.append((t_, p))
        if i1 - i0 == 2:
            pm = m.node(*P(r - delta, i0 + 1))   # 弦の中点（半径方向＝部材の法線方向）
            arm(pm, tn[i0 + 1], nr)
            mids[i0 + 1] = pm
            pts.append((0.5, pm))
        mems.append(dict(fam=fam, nr=nr, e0=P(r + delta, i0), e1=P(r + delta, i1), pts=pts))
        length_total += nr * (math.dist(P(r + delta, i0), P(r + delta, i1)) + 2 * g["over"])
    # A系と B系を密着させ、側面図で X 字に交差する点を接着する（両者を同じノードで剛に結ぶ）
    ab = []
    if cfg.get("ab_glue"):
        for ma in (x for x in mems if x["fam"] == "A"):
            for mb in (x for x in mems if x["fam"] == "B"):
                da, db = sub2(ma["e1"], ma["e0"]), sub2(mb["e1"], mb["e0"])
                den = da[0] * db[1] - da[1] * db[0]
                if abs(den) < 1e-12:
                    continue
                w0 = sub2(mb["e0"], ma["e0"])
                ta = (w0[0] * db[1] - w0[1] * db[0]) / den
                tb = (w0[0] * da[1] - w0[1] * da[0]) / den
                if 0.03 < ta < 0.97 and 0.03 < tb < 0.97:
                    nd = m.node(ma["e0"][0] + ta * da[0], ma["e0"][1] + ta * da[1])
                    ma["pts"].append((ta, nd))
                    mb["pts"].append((tb, nd))
                    cosang = abs(da[0] * db[0] + da[1] * db[1]) / math.hypot(*da) / math.hypot(*db)
                    ab.append((nd, math.sqrt(max(1e-6, 1 - cosang ** 2))))
    for x in mems:
        pts = [nd for _, nd in sorted(x["pts"])]
        x_nr = x["nr"]
        m.beam_line(pts, 2, E * A1 * x_nr, E * I1 * x_nr, kind=x["fam"], rows=x_nr, s=b, b=b, h=h,
                    lb=g["chord"] / 2)
    for t, d in ((tn[0], deck[0]), (tn[-1], deck[-1])):
        m.frame(d, t, E * A1 * 1e3, E * I1 * 1e3, hj=True, kind="L", rows=1, s=b)

    deck_top = h_str + 0.5
    w = {}
    sides = 2 if through else 1
    w["アーチ縦材"] = sides * length_total * b * h * rho
    gap = 0.0 if cfg.get("ab_glue") else GAP
    rib_w = (nA + nB) * b + (nA + nB - 1) * gap + 4.0     # リブ幅（横材は両端 2 mm ずつ突き出す）
    if through:
        z_c = CLEAR_EDGE + rib_w / 2
        full = [i for i in range(n + 1) if P(r, i)[1] - hh(i) >= deck_top + md.CLEAR_MIN]
        # 横材は各リブ内（長さ rib_w）。高い節点ではリブ間を 3mm 角のつなぎ材（長さ 2·z_c）で結ぶ
        t_len = (n + 1) * 2 * rib_w
        w["つなぎ材"] = len(full) * 2 * z_c * 9 * rho
        t_max = max(rib_w, 2 * z_c) if full else rib_w
        hset = load_nodes if cfg["hang"] == "載荷点のみ" else list(range(1, n))
        hg = af.HANGERS[cfg["hanger"]]
        hl = 0.0
        for i in hset:
            top = tn[i] if (i % 2 or i not in mids) else mids[i]
            m.truss(deck[i], top, hg["EA"], kind="H", tonly=cfg["hanger"] != "バルサ3mm角", cap=hg["T"])
            hl += math.dist(m.nodes[deck[i]], m.nodes[top])
        for i in load_nodes:
            m.loads[("v", deck[i])] = m.loads.get(("v", deck[i]), 0.0) - 0.25
        m.load_nodes = [deck[i] for i in load_nodes]
        span_c = 2 * z_c
        cross = len(load_nodes) * 75 * span_c + 2 * 25 * span_c
        cross_x = [xs[0], xs[-1]] + [xs[i] for i in load_nodes]
        if cfg["hang"] != "載荷点のみ" or cfg["sheet"] != "なし":
            cross += (n - 1 - len(load_nodes)) * 25 * span_c
            cross_x += [xs[i] for i in range(1, n) if i not in load_nodes]
        w["横桁"] = cross * rho
        w["ハンガー"] = 2 * hl * 9 * rho if cfg["hanger"] == "バルサ3mm角" else 0.0
        w["床版"] = SHEETS[cfg["sheet"]]
        joints = 2 * ((nA + nB) * (n + 1) + len(hset)) + 8
        width_ok = 2 * CLEAR_EDGE >= md.DECK_W_MIN
    else:
        t_len = (n + 1) * max(md.DECK_W, rib_w)
        t_max = max(md.DECK_W, rib_w)
        for i in load_nodes:
            m.loads[("v", tn[i])] = m.loads.get(("v", tn[i]), 0.0) - 0.5
        m.load_nodes = [tn[i] for i in load_nodes]
        span_c = md.DECK_W
        cross_x = list(xs)
        w["床版"] = 0.0
        joints = (nA + nB) * (n + 1) + 6
        width_ok = True
    w["横材"] = t_len * trs[0] * trs[1] * rho
    if tl:
        w["横材"] += len(wide) * sides * (rib_w if through else t_len / (n + 1)) * (tl[0] * tl[1] - trs[0] * trs[1]) * rho
    w["タイ（縦桁）"] = 2 * 514 * s_str * h_str * rho
    if cfg.get("load_beam"):
        w["荷重受け"] = 2 * 75 * md.DECK_W * rho
    w["接着剤"] = joints * af.GLUE_PER_JOINT
    if cfg["joint"] == "厚紙ガセット":
        w["交差部ガセット"] = (sides * (nA + nB) * (n + 1) * 2 * (GUSSET["L"] * (h + trs[1] + 4))
                          * GUSSET["t"] * md.RHO_PAPER)
    ok = g["feasible"] and g["length"] <= md.MAX_MEMBER and t_max <= md.MAX_MEMBER and span_c <= md.MAX_MEMBER and width_ok
    return m, dict(w=w, W=sum(w.values()), ok=ok, geo=g, load_nodes=load_nodes, xs=xs, t_max=t_max, ab=ab,
                   ab_if=nA + nB - 1, b=b, cross_x=sorted(cross_x), s_str=s_str, h_str=h_str,
                   crown=P(r, n // 2)[1], cross_span=span_c)


def run(cfg):
    m, info = build(cfg)
    if not info["ok"]:
        return None
    r = af.evaluate(m)
    caps = dict(r["caps"])
    if info["ab"]:
        caps["abjoint"] = (ab_joint_cap(m, info["ab"], info["b"], info["ab_if"]), "A–B 側面接着")
    if cfg["type"] == "下路":   # 載荷横桁（5×15 縦使い）の曲げ。P/2 を幅方向に等分布
        L = info["cross_span"]
        caps["cross"] = (MAT["fm"] * (5 * 15 ** 2 / 6) / (0.5 * L / 8), "載荷横桁の曲げ")
    else:                      # 上路: 載荷点の横材が P/2 を幅方向に受け、各列で支持される
        trs = SECTIONS[cfg["trans"]]
        nrow = sum(cfg["rows"])
        L = md.DECK_W / max(1, nrow)
        z = 5 * 15 ** 2 / 6 if cfg.get("load_beam") else trs[0] * trs[1] ** 2 / 6
        caps["cross"] = (MAT["fm"] * z / (0.5 * L / 8), "載荷点の横材の曲げ")
    # 縦桁の横倒れ座屈（横桁位置を横方向の支点とみなす）: M_cr = (π/Lb)·√(E·Iy·G·J)
    bs, hs = info["s_str"], info["h_str"]
    Lb = max(b_ - a_ for a_, b_ in zip(info["cross_x"], info["cross_x"][1:]))
    Iy, J = hs * bs ** 3 / 12, hs * bs ** 3 / 3 * (1 - 0.63 * min(bs, hs) / max(bs, hs))
    mcr = math.pi / Lb * math.sqrt(MAT["E"] * Iy * 0.054 * MAT["E"] * J)
    m_s = max(f["M"] for e, f in zip(m.elems, m.forces) if e.get("kind") == "S")
    caps["ltb"] = (mcr / m_s if m_s > 0 else float("inf"), "縦桁の横倒れ座屈 Lb=%.0f" % Lb)
    gov = min(caps.items(), key=lambda kv: kv[1][0])
    P = gov[1][0]
    W = info["W"]
    # アーチ分担率: ハンガーがアーチへ運ぶ荷重の割合（残りは縦桁が梁として運ぶ）
    share = (sum(f["N"] for e, f in zip(m.elems, m.forces) if e.get("kind") == "H" and not e.get("off")) / 0.5
             if cfg["type"] == "下路" else 1.0)
    return dict(cfg=cfg, P=P, W=W, CP=P / W ** 1.5, ratio=P / 9.81 / W * 1000, gov=gov, caps=caps, info=info,
                tie=r["tie"], defl=r["defl"], share=share)


def configs():
    for n, lng, trs, rows in itertools.product((8, 10, 12), SECTIONS, SECTIONS, ROWS):
        for joint in ("接着のみ", "厚紙ガセット"):
            for hanger, hang in itertools.product(("バルサ3mm角", "ナイロン3号"), ("全節点", "載荷点のみ")):
                for sheet in ("なし", "ケント紙0.2"):
                    yield dict(type="下路", n=n, long=lng, trans=trs, rows=rows, joint=joint, hanger=hanger,
                               hang=hang, sheet=sheet)
            for tie in (3.0, 5.0):
                yield dict(type="上路", n=n, long=lng, trans=trs, rows=rows, joint=joint, tie=tie, sheet="なし")


def refine_configs():
    """全探索の上位付近を、縦桁（タイ）断面と荷重受けを加えて絞り込む。"""
    ties = [(4.0, 4.0), (5.0, 5.0), (5.0, 15.0)]
    for lng, rows, hanger, hang, (tb, th) in itertools.product(("4角", "5角"), ((1, 1), (2, 1)),
                                                              ("バルサ3mm角", "ナイロン3号"), ("全節点", "載荷点のみ"), ties):
        yield dict(type="下路", n=8, long=lng, trans="3角", rows=rows, joint="厚紙ガセット", hanger=hanger,
                   hang=hang, sheet="なし", tie=tb, tie_h=th)
    for n, lng, trs, rows, (tb, th), lb in itertools.product((8, 10), ("4角", "5角"), ("3角", "4角"),
                                                            ((1, 1), (2, 1), (1, 2), (2, 2)),
                                                            [(3.0, 3.0)] + ties, (False, True)):
        yield dict(type="上路", n=n, long=lng, trans=trs, rows=rows, joint="厚紙ガセット", tie=tb, tie_h=th,
                   load_beam=lb, sheet="なし")


def nocard_configs():
    """厚紙を使わない（床版シート・交差部ガセットなし）下路案の探索。
    --soft: 交差部の柔らかさを考慮し、A–B 側面接着と突き出し 3 mm（全 5 角が収まる）も選択肢に加える。"""
    soft = "--soft" in sys.argv
    for n, lng, trs, rows, hanger, hang, th in itertools.product(
            (8, 10, 12), ("4角", "5角", "5×15平"), SECTIONS, ROWS, ("バルサ3mm角", "ナイロン3号"),
            ("全節点", "載荷点のみ"), (5.0, 15.0)):
        base = dict(type="下路", n=n, long=lng, trans=trs, rows=rows, joint="接着のみ", hanger=hanger,
                    hang=hang, sheet="なし", tie=5.0, tie_h=th)
        variants = [base]
        if trs == "5角" and not soft:
            variants.append(dict(base, trans_load="5×15平"))
        if soft:
            variants = [dict(v, soft_joints=True, ab_glue=ab, **({"over": 3.0} if ov else {}))
                        for v in variants for ab in (False, True) for ov in (False, True)]
        yield from variants


def label(c):
    s = "%s n=%d 縦材%s 横材%s%s 列A%d/B%d %s" % (c["type"], c["n"], c["long"], c["trans"],
                                                 "（載荷節点 %s）" % c["trans_load"] if c.get("trans_load") else "",
                                                 *c["rows"], c["joint"])
    if c["type"] == "下路":
        s += " %s(%s) 床版%s" % (c["hanger"], c["hang"], c["sheet"])
    else:
        s += " タイ%g×%g%s" % (c["tie"], c.get("tie_h", c["tie"]), " 荷重受け5×15" if c.get("load_beam") else "")
    if c["type"] == "下路" and "tie_h" in c:
        s += " 縦桁%g×%g" % (c["tie"], c["tie_h"])
    if c.get("ab_glue"):
        s += " A-B側面接着"
    if c.get("over") is not None:
        s += " 突出%g" % c["over"]
    return s


def main():
    results, cache = [], {}
    gen = refine_configs if "--refine" in sys.argv else nocard_configs if "--nocard" in sys.argv else configs
    for cfg in gen():
        key = tuple(sorted((k, str(v)) for k, v in cfg.items() if k != "sheet"))
        if key in cache:          # 床版シートは重量だけに効くので FEM を使い回す
            base = cache[key]
            if base is None:
                continue
            m, info = build(cfg)
            res = dict(base, cfg=cfg, W=info["W"], info=dict(base["info"], w=info["w"], W=info["W"]))
            res["CP"] = res["P"] / res["W"] ** 1.5
            res["ratio"] = res["P"] / 9.81 / res["W"] * 1000
        else:
            res = run(cfg)
            cache[key] = res
            if res is None:
                continue
        results.append(res)
    results.sort(key=lambda r: -r["CP"])
    print("評価した構成: %d（ルール適合）" % len(results))
    for typ in ("下路", "上路"):
        sub = [r for r in results if r["cfg"]["type"] == typ]
        print("\n=== %s 上位10" % typ)
        for r in sub[:10]:
            print("  CP %.2f  P %4.0f N  W %5.1f g  荷重/自重 %4.0f  アーチ分担 %3.0f%%  支配:%s  | %s"
                  % (r["CP"], r["P"], r["W"], r["ratio"], 100 * r["share"], r["gov"][1][1], label(r["cfg"])))
        arch = [r for r in sub if r["share"] >= 0.6]
        if typ == "下路" and arch:
            print("\n=== 下路のうちアーチ分担 60%% 以上（アーチが主に荷重を運ぶ案）上位10")
            for r in arch[:10]:
                print("  CP %.2f  P %4.0f N  W %5.1f g  荷重/自重 %4.0f  アーチ分担 %3.0f%%  支配:%s  | %s"
                      % (r["CP"], r["P"], r["W"], r["ratio"], 100 * r["share"], r["gov"][1][1], label(r["cfg"])))
    base = [r for r in results if r["cfg"] == dict(type="下路", n=8, long="5角", trans="3角", rows=(2, 1),
                                                   joint="接着のみ", hanger="ナイロン3号", hang="全節点",
                                                   sheet="ケント紙0.2")]
    if base:
        r = base[0]
        print("\n参考 現案相当（床版ケント紙）: CP %.2f P %.0f N W %.1f g 支配:%s" % (r["CP"], r["P"], r["W"], r["gov"][1][1]))
    return results


if __name__ == "__main__":
    main()
