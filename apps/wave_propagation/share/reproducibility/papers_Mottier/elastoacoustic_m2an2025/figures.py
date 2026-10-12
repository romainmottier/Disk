#!/usr/bin/env python3
"""
figures.py — figures et tableaux reproduits de [M2AN25] (appele par run_campaign.py / plot_mottier.py).

Conventions identifiees en comparant aux valeurs publiees (voir ../LISEZMOI.txt, section 4) :
  - normes des Figs. 5 a 7 : ||.||_HHO = (e_p^2 + e_v^2 / 2)^1/2, ||.||_dG = (e_m^2 + e_s^2 / 2)^1/2, ou e_p, e_v sont les
    erreurs L2 des inconnues de cellule contre les projections L2 de la solution exacte et e_m, e_s les erreurs L2
    (quadrature, norme L2 simple) des inconnues duales : contributions elastiques divisees par deux. La definition
    (6.3) de l'article (normes d'energie) est tracee en variante ;
  - traces aux capteurs (Figs. 10 a 12) dans les unites de Gar6more2D (source ponctuelle d'amplitude 1) : valeurs HHO
    divisees par S = c rho theta lambda^3 / (2 pi^5/2) (../gar6more.py) ;
  - Tables 1 et 2 : valeurs publiees reproduites avec les poids elastiques eta^s = 1.75 (ordre egal) et 1.04 (ordre
    mixte) au lieu des eta_* de (6.2) (1.54 et 1.38) ; les deux sont calcules.
"""
import math
import numpy as np
import matplotlib.pyplot as plt

import plot_mottier as pm
from plot_mottier import GREEN, ORANGE, DARKRED, NAVY

PAPER = "M2AN25"
INTRO = """
Reproduction des cas tests de la section 6 de [M2AN25] avec le prototype d'ordre un de wave_propagation (fichiers
d'entree Lua de ce dossier, lances par run_campaign.py). Chaque section : figure reproduite (panneaux disposes comme
dans l'article), comparaison avec les valeurs relevees sur la figure publiee (valeurs_article.json) et ecarts.
"""

K_COLORS = {1: GREEN, 2: ORANGE, 3: DARKRED, 4: NAVY}         # Figs. 5-7, 10-12
K_MARKERS = {1: "s", 2: "o", 3: "D", 4: "D"}
SPEC_COLORS = {1: ORANGE, 2: DARKRED, 3: NAVY}                # Fig. 3, Fig. 13


CONVENTIONS = {
    "plein": "(e_p^2 + e_v^2)^1/2 et (e_m^2 + e_s^2)^1/2, erreurs L2 des inconnues de cellule (Fig. 5)",
    "article": "(e_p^2 + e_v^2/2)^1/2 et (e_m^2 + e_s^2/2)^1/2, erreurs L2 (contributions elastiques divisees par deux)",
    "projection": "(e_p^2 + e_v^2)^1/2 et (|Pi m - m_T|^2 + |Pi s - s_T|^2)^1/2, erreurs dG contre les projections L2",
    "definition": "normes d'energie de la definition (6.3) : (e_p^2 + e_v^2)^1/2 et (e_m^2 + |s - s_T|_C^-1^2)^1/2"}


def norms(d, i, convention="article"):
    """(HHO, dG) de la ligne i d'un tableau de composantes (output.errors = "energy")"""
    p, v, m = d["hho_p"][i], d["hho_v"][i], d["dg_m"][i]
    if convention == "article":
        s = d["dg_s_plain"][i]
        return math.sqrt(p * p + 0.5 * v * v), math.sqrt(m * m + 0.5 * s * s)
    if convention == "plein":
        return math.hypot(p, v), math.hypot(m, d["dg_s_plain"][i])
    if convention == "projection":
        return math.hypot(p, v), math.hypot(d["dg_m_projection"][i], d["dg_s_projection_plain"][i])
    return math.hypot(p, v), math.hypot(m, d["dg_s"][i])        # definition (6.3)


def _pv(ctx, *keys):
    d = ctx.paper
    for k in keys:
        if not isinstance(d, dict) or str(k) not in d:
            return None
        d = d[str(k)]
    return d


# ------------------------------------------------------------------------------------------------------------------
def fig3(ctx):
    """Fig. 3 — rayon spectral en fonction de eta (acoustique pur, elastique pur, couplage)"""
    item = "Fig. 3"
    titles = {"acoustic": "Purely Acoustic Case", "elastic": "Purely Elasticity Case", "coupled": "Coupling Case"}
    etas = [2.0 ** w for w in range(-3, 4)]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    rows = []
    for ax, medium in zip(axs, ("acoustic", "elastic", "coupled")):
        for order, mk in (("equal", "o"), ("mixed", "d")):
            for k in (1, 2, 3):
                ys = []
                for w in range(-3, 4):
                    r = ctx.case(f"m2an_fig3_{medium}_{'eo' if order == 'equal' else 'mo'}_k{k}_w{w:+d}")
                    ys.append(r["data"]["rho_normalized"][0] if r else None)
                xs = [e for e, y in zip(etas, ys) if y is not None]
                ax.plot(xs, [y for y in ys if y is not None], mk, color=SPEC_COLORS[k], ms=7 if mk == "o" else 8)
                pv = _pv(ctx, item, medium, order, k)
                if pv:
                    ax.plot(etas, pv, mk, mfc="none", mec="k", ms=11, mew=0.8)
                    rows.append([medium, order, k, pm.fmt_ratio(pm.ratio_stats(ys, pv))])
        pm.loglog(ax, r"$\eta$", "Spectral radius")
        ax.set_xticks(etas); ax.set_xticklabels([f"$2^{{{w}}}$" for w in range(-3, 4)])
        ax.set_title(titles[medium])
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=SPEC_COLORS[k], lw=2, label=f"k = {k}") for k in (1, 2, 3)]
    h += [Line2D([], [], color="k", marker="o", ls="", label="Equal-Order"), Line2D([], [], color="k", marker="d", ls="",
          label="Mixed-Order"), Line2D([], [], color="k", marker="o", mfc="none", ls="", ms=10, label="article (releve)")]
    axs[2].legend(handles=h, loc="upper right", ncol=2)
    img = ctx.save(fig, item, "M2AN25_figure_3")
    notes = ["Maillage : 8 x 4 carres de (-1,1) x (0,1) (32 cellules, \"l1\" des fichiers de l'article), le meme pour les "
             "cas purs ; rayon spectral normalise par sqrt(#cellules) ; eta^f = eta^s = eta. Cercles vides : valeurs "
             "relevees sur la figure publiee."]
    ctx.section(item, "Rayon spectral en fonction de eta", [img], notes,
                (["milieu", "ordre", "k", "nous / article : mediane [min, max] (points)"], rows, "") if rows else None,
                status="reproduit" if rows else "")


def tab1(ctx):
    """Table 1 — rayon spectral du couplage, eta = 2^w eta_*"""
    item = "Table 1"
    ws = list(range(-3, 4))
    head = ["", "k"] + [f"{2.0 ** w:g}" if w >= 0 else f"1/{2 ** -w}" for w in ws]
    rows_t, rows_e, md = [], [], []
    stats = {"": [], "_eff": []}
    for order in ("equal", "mixed"):
        for k in (1, 2, 3):
            pv = _pv(ctx, item, order, k)
            lt, le = [order, k], [order, k]
            for w in ws:
                vals = {}
                for tag in ("", "_eff"):
                    r = ctx.case(f"m2an_tab1{tag}_{'eo' if order == 'equal' else 'mo'}_k{k}_w{w:+d}")
                    vals[tag] = r["data"]["rho_normalized"][0] if r else None
                    if vals[tag] and pv:
                        stats[tag].append(vals[tag] / pv[w + 3])
                lt.append(f"{vals['']:.1f}" if vals[""] else "-")
                le.append(f"{vals['_eff']:.1f}" if vals["_eff"] else "-")
            rows_t.append(lt); rows_e.append(le)
            if pv:
                md.append([order, k] + [f"{p:.1f}" for p in pv])
    imgs = [pm.table_image(ctx, item, "M2AN25_tableau_1", head, rows_t,
                           "Rayon spectral, poids eta_* de (6.2) : 0.88 / 1.54 (ordre egal), 0.80 / 1.38 (mixte)"),
            pm.table_image(ctx, item, "M2AN25_tableau_1_poids_effectifs", head, rows_e,
                           "Variante : poids elastiques 1.75 (ordre egal) et 1.04 (ordre mixte)")]
    if md:
        imgs.append(pm.table_image(ctx, item, "M2AN25_tableau_1_article", head, md, "Article (valeurs publiees)"))
    notes = ["Les valeurs publiees ne sont pas celles des poids eta_* annonces : en ordre egal, nos valeurs avec eta_* "
             "sont 1.14 fois les publiees pour w < 0 et 0.88 fois pour w >= 0, rapports constants qui correspondent a un "
             "poids elastique 1.54 x 1.14 = 1.75 ; en ordre mixte (0.75 a 0.84 puis 1.32), le poids elastique 1.04 "
             "reproduit le tableau. Le rayon spectral du couplage ne depend pratiquement que du poids elastique. La "
             "quadrature d'interface du code de l'article (space.interface_quadrature = 'legacy') ne change ces valeurs "
             "que de moins de 1 %."]
    tab = None
    if stats[""]:
        tab = (["poids", "nous / article : mediane [min, max] (points)"],
               [["eta_* de (6.2)", pm.fmt_ratio((float(np.median(stats[""])), min(stats[""]), max(stats[""]), len(stats[""])))],
                ["poids effectifs", pm.fmt_ratio((float(np.median(stats["_eff"])), min(stats["_eff"]), max(stats["_eff"]),
                                                  len(stats["_eff"]))) if stats["_eff"] else "-"]], "")
    ctx.section(item, "Rayon spectral du couplage en fonction de eta / eta_*", imgs, notes, tab,
                status="reproduit avec les poids effectifs (variante)" if stats["_eff"] else "")


def tab2(ctx):
    """Table 2 — rayon spectral selon la forme des cellules, k = 1..6"""
    item = "Table 2"
    geoms = ("simplexes", "quadrangles", "polygones")
    head = ["k"] + [f"{g} {o}" for g in geoms for o in ("EO", "MO")]
    out = {}
    for tag in ("", "_eff"):
        rows = []
        for k in range(1, 7):
            line = [k]
            for g in geoms:
                for o in ("eo", "mo"):
                    r = ctx.case(f"m2an_tab2{tag}_{g[:4]}_{o}_k{k}")
                    line.append(f"{r['data']['rho_normalized'][0]:.1f}" if r else "-")
            rows.append(line)
        out[tag] = rows
    imgs = [pm.table_image(ctx, item, "M2AN25_tableau_2", head, out[""], "Rayon spectral, poids eta_* de (6.2)"),
            pm.table_image(ctx, item, "M2AN25_tableau_2_poids_effectifs", head, out["_eff"],
                           "Variante : poids elastiques effectifs 1.75 (EO) / 1.04 (MO)")]
    pv = _pv(ctx, item)
    rows_cmp = []
    if pv:
        md = [[k] + [f"{pv[g][o][k - 1]:.1f}" for g in geoms for o in ("equal", "mixed")] for k in range(1, 7)]
        imgs.append(pm.table_image(ctx, item, "M2AN25_tableau_2_article", head, md, "Article (valeurs publiees)"))
        for g in geoms:
            for o, oo in (("eo", "equal"), ("mo", "mixed")):
                ours = []
                for k in range(1, 7):
                    r = ctx.case(f"m2an_tab2_eff_{g[:4]}_{o}_k{k}")
                    ours.append(r["data"]["rho_normalized"][0] if r else None)
                rows_cmp.append([g, oo, pm.fmt_ratio(pm.ratio_stats(ours, pv[g][oo]))])
    notes = ["Maillages : 8 x 4 carres (32 cellules), triangles gmsh lc = 0.35 (52 cellules) et polygones PolyMesher a "
             "32 germes par carre unite (64 cellules), symetriques par rapport a x = 0 ; les maillages de l'article ne sont "
             "pas disponibles (simplexes et polygones regeneres, voir ../meshes.py)."]
    ctx.section(item, "Rayon spectral selon la forme des cellules", imgs, notes,
                (["maillage", "ordre", "nous (poids effectifs) / article"], rows_cmp, "") if rows_cmp else None,
                status="reproduit (simplexes, quadrangles) ; polygones : maillages regeneres" if rows_cmp else "")


def fig4(ctx):
    """Fig. 4 — maillages cartesien et polygonal (l = 2)"""
    item = "Fig. 4"
    import meshes
    fig, axs = plt.subplots(1, 2, figsize=(11, 3))
    nx, ny = 16, 8
    xs, ys = np.linspace(-1, 1, nx + 1), np.linspace(0, 1, ny + 1)
    pts = np.array([[x, y] for y in ys for x in xs])
    cells = [[j * (nx + 1) + i, j * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i]
             for j in range(ny) for i in range(nx)]
    tags = [1 if pts[c].mean(axis=0)[0] < 0 else 2 for c in cells]
    pm.draw_mesh(axs[0], pts, cells, tags, {1: "#b00000", 2: "#0000ff"})
    p = ctx.paper_dir / "meshes" / "poly_128.txt"
    src = p if p.exists() else ctx.paper_dir / "meshes" / "poly_128.txt.bz2"
    if src.exists():
        pts, cells, tags = meshes.read_polygons(src)
        pm.draw_mesh(axs[1], pts, cells, tags, {1: "#b00000", 2: "#0000ff"})
    img = ctx.save(fig, item, "M2AN25_figure_4")
    ctx.section(item, "Maillages cartesien et polygonal (l = 2)", [img],
                ["Maillage polygonal regenere (PolyMesher, 128 germes par carre unite, symetrise) : le fichier de "
                 "l'article n'est pas disponible."], status="dessin")


# ------------------------------------------------------------------------------------------------------------------
def _conv_panel(ctx, ax, cases, xcol, labels, colors, markers, convention="article"):
    out = {}
    for name, key in cases:
        r = ctx.case(name)
        if not r or not r.get("data"):
            continue
        d = r["data"]
        x = d[xcol]
        hho, dg = zip(*[norms(d, i, convention) for i in range(len(x))])
        c, mk = colors[key], markers[key]
        ax.plot(x, dg, "-", color=c, marker=mk, ms=5, label=labels[key])
        ax.plot(x, hho, "--", color=c, marker=mk, ms=5)
        out[key] = (list(x), list(dg), list(hho))
    return out


def _conv_compare(ctx, item, panel, ours, pkey_fn):
    """Rapports nous / article point par point, les deux series triees par abscisse croissante (les releves ne sont pas
    tous dans le meme ordre que nos niveaux)"""
    rows = []
    for key, (x, dg, hho) in ours.items():
        pv = pkey_fn(key)
        if not pv:
            continue
        o = sorted(range(len(x)), key=lambda i: x[i])
        p = sorted(range(len(pv["x"])), key=lambda i: pv["x"][i])
        n = min(len(o), len(p))
        rows.append([panel, key, pm.fmt_ratio(pm.ratio_stats([dg[i] for i in o[:n]], [pv["dG"][i] for i in p[:n]])),
                     pm.fmt_ratio(pm.ratio_stats([hho[i] for i in o[:n]], [pv["HHO"][i] for i in p[:n]]))])
    return rows


def fig5(ctx):
    """Fig. 5 — convergence en temps (solution polynomiale en espace, k = 4, l = 1)"""
    item = "Fig. 5"
    imgs, rows = [], []
    for tag, suffix in (("", ""), ("_o1h", "_variante_erk_O1h")):
        fig, axs = plt.subplots(1, 2, figsize=(11, 5))
        for ax, scheme, stages, lab in ((axs[0], "erk", (2, 3, 4), "ERK({s})"),
                                        (axs[1], "sdirk", (1, 2, 3), "SDIRK({s},{t})")):
            cols = dict(zip(stages, (ORANGE, DARKRED, NAVY)))
            mks = dict(zip(stages, ("s", "o", "D")))
            labels = {s: lab.format(s=s, t=s + 1) for s in stages}
            case_tag = tag if scheme == "erk" else ""
            ours = _conv_panel(ctx, ax, [(f"m2an_fig5_{scheme}{s}{case_tag}", s) for s in stages], "dt", labels, cols,
                               mks, convention="plein")
            side = "gauche" if scheme == "erk" else "droite"
            for s in ours:
                pv = _pv(ctx, item, side, f"s{s}")
                if pv:
                    ax.plot(pv["x"], pv["dG"], "o", mfc="none", mec="k", ms=9, mew=0.7)
                    ax.plot(pv["x"], pv["HHO"], "s", mfc="none", mec="0.4", ms=9, mew=0.7)
            if not tag or scheme == "erk":
                rows += [[("ERK O(1/h)" if tag else ("ERK O(1)" if scheme == "erk" else "SDIRK O(1/h)"))] + r[1:]
                         for r in _conv_compare(ctx, item, scheme, ours, lambda s: _pv(ctx, item, side, f"s{s}"))]
            pm.loglog(ax, r"$\delta t$", "L$^2$-error"); pm.legend(ax, loc="lower right")
            if tag and scheme == "erk":
                ax.set_yscale("log"); ax.set_ylim(1e-15, 1e-1)
            ax.set_title(("ERK(s), O(1)" if not tag else "ERK(s), O(1/h) (legende)") if scheme == "erk"
                         else "SDIRK(s, s+1), O(1/h)")
        if tag:
            fig.suptitle("Variante : ERK(s) avec la stabilisation O(1/h) de la legende (instables aux deux plus grands pas)",
                         fontsize=9)
        imgs.append(ctx.save(fig, item, "M2AN25_figure_5" + suffix))
    notes = ["Trait plein : ||.||_dG, tirets : ||.||_HHO, normes " + CONVENTIONS["plein"] + " (cette figure : rapports "
             "1.00, sauf ||.||_HHO des schemas d'ordre 4, 1.4 a tous les pas ; la convention des Figs. 6 et 7, contributions "
             "elastiques divisees par deux, donne 0.7 a 0.8 pour ||.||_dG) ; symboles "
             "vides : valeurs relevees sur l'article. La legende de l'article inverse les panneaux (ERK a gauche dans la "
             "figure publiee). Les ERK(s) publies sont stables aux six pas de temps, ce que donne la stabilisation O(1) ; "
             "avec la O(1/h) de la legende (variante), ils explosent aux deux plus grands pas (meme erreur ailleurs)."]
    ctx.section(item, "Convergence en temps", imgs, notes,
                (["schema", "s", "dG nous / article", "HHO nous / article"], rows, "") if rows else None,
                status="reproduit (ERK avec stabilisation O(1))" if rows else "")


def _fig67(ctx, item, num, title):
    titles = {"gauche": "equal order, O(1), ERK(4)", "centre": "mixed order, O(1), ERK(4)",
              "droite": "mixed order, O(1/h), SDIRK(3,4)"}
    labels = {k: f"k = {k}" for k in (1, 2, 3, 4)}
    imgs, rows = [], []
    for conv, suffix, tag in (("article", "", ""), ("projection", "_variante_normes_projection", ""),
                              ("definition", "_variante_normes_6_3", ""), ("projection", "_variante_poids_MO_1.04", "_mo104")):
        fig, axs = plt.subplots(1, 3, figsize=(16, 5))
        any_case = False
        for ax, panel in zip(axs, ("gauche", "centre", "droite")):
            if tag and panel == "gauche":
                ax.axis("off"); continue
            cases = [(f"m2an_fig{num}_{panel}_k{k}{tag}", k) for k in (1, 2, 3, 4)]
            ours = _conv_panel(ctx, ax, cases, "h", labels, K_COLORS, K_MARKERS, convention=conv)
            any_case = any_case or bool(ours)
            for k in ours:
                pv = _pv(ctx, item, panel, f"k{k}")
                if pv:
                    ax.plot(pv["x"], pv["dG"], "o", mfc="none", mec="k", ms=9, mew=0.7)
                    ax.plot(pv["x"], pv["HHO"], "s", mfc="none", mec="0.4", ms=9, mew=0.7)
            for r in _conv_compare(ctx, item, panel, ours, lambda k: _pv(ctx, item, panel, f"k{k}")):
                rows.append([conv + (" poids MO 1.04" if tag else "")] + r)
            pm.loglog(ax, "$h$", "L$^2$-error"); pm.legend(ax, loc="lower right"); ax.set_title(titles[panel], fontsize=10)
        if not any_case:
            plt.close(fig); continue
        if conv != "article" or tag:
            fig.suptitle("Variante : " + CONVENTIONS[conv] + (" ; poids elastique 1.04 en ordre mixte" if tag else ""),
                         fontsize=9)
        imgs.append(ctx.save(fig, item, f"M2AN25_figure_{num}{suffix}"))
    notes = ["Trait plein : ||.||_dG, tirets : ||.||_HHO ; symboles vides : valeurs relevees sur l'article ; n = 8 ; poids "
             "eta_* de (6.2). Figure principale : " + CONVENTIONS["article"] + ", convention qui reproduit le panneau de "
             "gauche (ordre egal) a quelques % pres. Les panneaux en ordre mixte sont mieux reproduits par les erreurs par "
             "projection (variante) et, pour O(1), avec le poids elastique 1.04 identifie sur la Table 1 ; la definition (6.3) "
             "de l'article est tracee en variante. Rapports nous / article par convention dans le tableau."]
    ctx.section(item, title, imgs, notes,
                (["convention", "panneau", "k", "dG nous / article", "HHO nous / article"], rows, "") if rows else None,
                status="compare" if rows else "")


def fig6(ctx):
    """Fig. 6 — convergence en espace, maillages cartesiens"""
    _fig67(ctx, "Fig. 6", 6, "Convergence en espace, maillages cartesiens")


def fig7(ctx):
    """Fig. 7 — convergence en espace, maillages polygonaux"""
    _fig67(ctx, "Fig. 7", 7, "Convergence en espace, maillages polygonaux")


# ------------------------------------------------------------------------------------------------------------------
def fig8(ctx):
    """Fig. 8 — domaine du cas test de l'ondelette de Ricker"""
    item = "Fig. 8"
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 0.5, color="#f2c9c9"))
    ax.add_patch(plt.Rectangle((-0.5, 0.0), 1, 0.5, color="#c9c9f2"))
    ax.plot([-0.5, 0.5], [0, 0], "k", lw=2)
    ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fill=False, lw=2))
    ax.plot(0.0, 0.125, "o", color="#c000c0", ms=8)
    for (x, y, n) in ((-0.15, 0.1, r"$S^f$"), (-0.15, -0.1, r"$S^s$")):
        ax.plot(x, y, "o", color="teal", ms=7); ax.text(x - 0.13, y, n, color="teal", fontsize=12)
    ax.text(-0.05, 0.3, r"$\Omega^f$", color="b", fontsize=16); ax.text(-0.05, -0.35, r"$\Omega^s$", color="#b00000", fontsize=16)
    ax.set_aspect("equal"); ax.axis("off")
    img = ctx.save(fig, item, "M2AN25_figure_8")
    ctx.section(item, "Domaine de l'ondelette de Ricker", [img],
                ["(-0.5, 0.5)^2, fluide y > 0, solide y < 0, source (0, 0.125), capteurs du cas academique."], status="dessin")


def _snapshot_limits(frames, pct):
    """Echelles de couleur communes aux instantanes : pression symetrique au centile pct[0] de |p|, norme de la vitesse
    au centile pct[1] (la figure publiee ecrete aussi |v|, a 2.5)"""
    za = np.concatenate([np.abs(np.concatenate([v for (m, _), v in zip(c, vals) if m == 1] or [np.zeros(1)]))
                         for c, vals, _ in frames])
    ze = np.concatenate([np.concatenate([v for (m, _), v in zip(c, vals) if m == 0] or [np.zeros(1)])
                         for c, vals, _ in frames])
    a = float(np.percentile(za, pct[0])) or 1.0
    e = float(np.percentile(ze, pct[1])) or 1.0
    return (-a, a), (0.0, e)


def _snapshots(ctx, item, case, title, name, cols=2, notes=(), variants=(), pct=None):
    """Instantanes du cas (et des variantes [(cas, suffixe, legende)]) ; pct : centiles d'ecretage des echelles
    communes (sinon chaque instantane a son echelle)"""
    imgs, done = [], False
    for c, suffix, label in [(case, "", "")] + list(variants):
        r = ctx.case(c)
        if not r or not r.get("snapshots"):
            continue
        done = done or c == case
        d = pm.Path(r["_dir"])
        frames = [pm.read_snapshot(d / f"{c}_cells.dat", d / f) for f in r["snapshots"]]
        clim_a, clim_e = _snapshot_limits(frames, pct) if pct else (None, None)
        n = len(frames)
        rows = (n + cols - 1) // cols
        fig, axs = plt.subplots(rows, cols, figsize=(5.2 * cols, 6.0 * rows), squeeze=False)
        for ax, (cells, vals, t) in zip(axs.flat, frames):
            out = pm.draw_snapshot(ax, cells, vals, cmap_a="RdBu_r", cmap_e="RdBu_r", clim_a=clim_a, clim_e=clim_e)
            if 1 in out:
                fig.colorbar(out[1], ax=ax, location="top", shrink=0.8, label=r"$p_{\mathcal{T}^F}$")
            if 0 in out:
                fig.colorbar(out[0], ax=ax, location="bottom", shrink=0.8, label=r"$\|\mathbf{v}_{\mathcal{T}^S}\|$")
            ax.set_title(f"t = {t:.4g}", fontsize=10)
        for ax in list(axs.flat)[n:]:
            ax.axis("off")
        if label:
            fig.suptitle("Variante : " + label, fontsize=10)
        imgs.append(ctx.save(fig, item, name + suffix))
    if not imgs:
        ctx.section(item, title, [], ["Cas non disponible."], status="non calcule")
        return
    txt = ("Pression acoustique (fluide) et norme de la vitesse (solide), polynomes de cellule evalues aux sommets et au "
           "barycentre" + (f" ; echelles communes aux instantanes, ecretees aux centiles {pct[0]:g} de |p| et {pct[1]:g} "
                           "de |v| (la figure publiee ecrete |v| a 2.5)" if pct else "") + ".")
    ctx.section(item, title, imgs, [txt] + list(notes),
                status="reproduit (qualitatif)" if done else "variante seulement (cas principal non calcule)")


def fig9(ctx):
    """Fig. 9 — instantanes, cas academique (k = 1, l = 7)"""
    _snapshots(ctx, "Fig. 9", "m2an_fig9_instantanes", "Instantanes du cas academique", "M2AN25_figure_9")


def _traces(ctx, item, num, levels, refname, title, tmax, tvalid=None, extra_notes=()):
    """tvalid : fin de la fenetre de comparaison a Gar6more2D (avant les reflexions sur les bords du domaine de calcul,
    absentes des demi-espaces infinis de Gar6more2D)"""
    import run_campaign as rc
    ref = [r for r in rc.REFERENCES if r.name == refname][0]
    S = ref.scale()[1]
    quantities = (("Sf", "p", "Acoustic pressure"), ("Ss", "vx", r"$v_x$ in bedrock"), ("Ss", "vy", r"$v_y$ in bedrock"))
    fig, axs = plt.subplots(3, 2, figsize=(12, 13))
    rows = []
    for j, l in enumerate(levels):
        for i, (sensor, q, ylab) in enumerate(quantities):
            ax = axs[i, j]
            g = ctx.reference(refname, sensor)
            for k in (1, 2, 3, 4):
                r = ctx.case(f"m2an_fig{num}_l{l}_k{k}")
                if not r or not r.get("sensors") or sensor not in r["sensors"]:
                    continue
                s = r["sensors"][sensor]
                ax.plot(s["t"], np.array(s[q]) / S, color=K_COLORS[k], lw=1.6, label=f"k = {k}")
                if g:
                    ts = np.array(s["t"])
                    w = ts <= (tvalid if tvalid else ts[-1])
                    gi = np.interp(ts[w], g["t"], g[q])
                    err = np.linalg.norm(np.array(s[q])[w] - gi) / max(np.linalg.norm(gi), 1e-300)
                    rows.append([f"l={l}", sensor, q, k, f"{err:.3e}"])
            if g:
                ax.plot(g["t"], np.array(g[q]) / S, "k", lw=1.6, label="Gar6more")
            pv = _pv(ctx, item, f"l{l}", q)
            for key, cur in (pv or {}).items():
                col = "0.5" if key == "Gar6more" else K_COLORS.get(int(key[1:]), "0.5") if key[1:].isdigit() else "0.5"
                ax.plot(cur["t"], cur["v"], ls=":", color=col, lw=1.0)
            ax.set_xlim(0, tmax); ax.set_xlabel("Time"); ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(f"l = {l}")
            if i == 0 and j == 0:
                pm.legend(ax, loc="upper left")
    img = ctx.save(fig, item, f"M2AN25_figure_{num}")
    notes = [f"Unites de Gar6more2D (valeurs HHO divisees par S = {S:.6g}) ; pointilles : courbes relevees sur l'article "
             "(Gar6more en gris, k = 1 en vert). Ecarts relatifs L2 en temps a la reference Gar6more2D dans le tableau "
             "(retard du signal source de Gar6more2D : 1.5 / f0, voir run_campaign.py)."] + list(extra_notes)
    ctx.section(item, title, [img], notes,
                (["niveau", "capteur", "champ", "k", "ecart L2 relatif a Gar6more"], rows, "") if rows else None,
                status="compare" if rows else "")


def fig10(ctx):
    """Fig. 10 — ondelette de Ricker, milieux academiques, comparaison a Gar6more2D"""
    _traces(ctx, "Fig. 10", 10, (4, 5), "m2an_academique", "Ricker academique : capteurs et Gar6more2D", 0.25)


def fig11(ctx):
    """Fig. 11 — ondelette de Ricker, granite-eau, comparaison a Gar6more2D"""
    _traces(ctx, "Fig. 11", 11, (5, 6), "m2an_granite_eau", "Ricker granite-eau : capteurs et Gar6more2D", 0.5)


def fig12(ctx):
    """Fig. 12 — ondelette de Ricker, granite-air, comparaison a Gar6more2D"""
    _traces(ctx, "Fig. 12", 12, (5, 6), "m2an_granite_air", "Ricker granite-air : capteurs et Gar6more2D", 0.5,
            tvalid=0.42, extra_notes=["Ecarts calcules pour t <= 0.42 : ensuite, les ondes rapides du solide (c_p = 6.36) "
                                      "reflechies sur le bord du domaine de calcul atteignent S^s, ce que la solution de "
                                      "Gar6more2D (demi-espaces infinis) ne contient pas."])


# ------------------------------------------------------------------------------------------------------------------
def _repartition(ax, r, scale, label_scale, ls="-", lw=2.5, suffix=""):
    e = r["energy"]
    ax.plot(e["t"], np.array(e["E"]) / scale, ls, color=NAVY, lw=lw, label="Total energy" + suffix)
    ax.plot(e["t"], np.array(e["E_acoustic"]) / scale, ls, color=GREEN, lw=lw, label="Acoustic energy" + suffix)
    ax.plot(e["t"], np.array(e["E_elastic"]) / scale, ls, color=DARKRED, lw=lw, label="Elastic energy" + suffix)
    ax.set_xlabel("Time"); ax.set_ylabel(f"Energy $\\times$ {label_scale}")
    pm.legend(ax)


def fig13(ctx):
    """Fig. 13 — energie : repartition (k = 3, l = 6, t <= 10) et perte relative (k = 1..3, l = 4..6)"""
    item = "Fig. 13"
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    notes, rows = [], []
    r = ctx.case("m2an_fig13_repartition")
    if r and r.get("energy"):
        _repartition(axs[0], r, 1e-5, "$10^5$")
        e = r["energy"]
        pv = _pv(ctx, item, "repartition")
        if pv:
            for key, c in (("total", NAVY), ("acoustic", GREEN), ("elastic", DARKRED)):
                if key in pv:                           # releves en unites absolues, comme nos energies
                    axs[0].plot(pv[key]["t"], np.array(pv[key]["v"]) / 1e-5, ":", color="0.4", lw=1)
        rows.append(["repartition k=3 l=6", f"E0 = {e['E'][0]:.4e}", f"E(T) / E0 = {e['E'][-1] / e['E'][0]:.6f}"])
        if pv and "acoustic" in pv:
            t, ea = np.array(e["t"]), np.array(e["E_acoustic"])
            tp, vp = np.array(pv["acoustic"]["t"]), np.array(pv["acoustic"]["v"])
            for lo, hi in ((0.0, 2.0), (2.0, 10.0)):
                m = (tp >= lo) & (tp <= hi)
                gap = np.sqrt(np.mean((np.interp(tp[m], t, ea) - vp[m]) ** 2)) / e["E"][0]
                rows.append([f"energie acoustique, {lo:g} <= t <= {hi:g}", "ecart quadratique moyen / E0", f"{gap:.3f}"])
    ls = {4: ":", 5: "--", 6: "-"}
    for l in (4, 5, 6):
        for k in (1, 2, 3):
            r = ctx.case(f"m2an_fig13_perte_l{l}_k{k}")
            if not r or not r.get("energy"):
                continue
            t = np.array(r["energy"]["t"]); E = np.array(r["energy"]["E"])
            loss = (E[0] - E) / E[0]
            m = (t > 0) & (loss > 0)
            axs[1].plot(t[m], loss[m], ls[l], color=SPEC_COLORS[k], lw=2.5)
            pv = _pv(ctx, item, "perte", f"l{l}_k{k}")
            st = ""
            if pv:
                ours = np.interp(pv["t"], t, loss)
                st = pm.fmt_ratio(pm.ratio_stats(list(ours), pv["v"]))
            rows.append([f"perte l={l} k={k}", f"perte(1) = {loss[-1]:.3e}", st])
    pm.loglog(axs[1], "Time", "Relative energy loss")
    from matplotlib.lines import Line2D
    axs[1].legend(handles=[Line2D([], [], color=SPEC_COLORS[k], lw=2, label=f"k = {k}") for k in (1, 2, 3)]
                  + [Line2D([], [], color="k", ls=ls[l], label=f"l = {l}") for l in (4, 5, 6)], ncol=2)
    img = ctx.save(fig, item, "M2AN25_figure_13")
    ctx.section(item, "Energie : repartition et perte relative", [img],
                ["Pointilles gris : courbes relevees sur l'article. Perte relative (droite) : k = 1, l = 4 concorde (0.41 "
                 "contre 0.35 a t = 1) ; pour les discretisations plus fines nos pertes sont bien plus faibles (k = 3, "
                 "l = 6 : 2.5e-5 contre 6e-4 a t = 1) et n'augmentent pas quand l'onde atteint l'interface (t ~ 0.05), "
                 "contrairement aux courbes publiees. Avec la quadrature d'interface du code de l'article "
                 "(space.interface_quadrature = 'legacy', regle de degre k), la perte de k = 3, l = 4 a t = 1 passe de "
                 "9.4e-3 a 1.8e-2 (article 1.4e-2), mais le saut publie vers t = 0.05 n'est pas reproduit : ecart non "
                 "explique."],
                (["courbe", "valeur", "nous / article (perte)"], rows, "") if rows else None,
                status="compare" if rows else "")


REALISTIC_NOTE = ("Granite-eau realiste : 280 x 280 cellules (h_x = 17.9 m, h_y = 12.5 m) au lieu des 560 x 560 de "
                  "l'article, qui demandent environ 85 Go de memoire avec l'assembleur actuel (voir "
                  "fig14_15_granite_eau_realiste.lua ; mesh.level = 2 redonne le maillage de l'article).")


SOURCE_NOTE = ("Source en (0, 100) : les figures publiees sont celles d'une source a environ 100 m de l'interface "
               "(centre de l'anneau de pression a 90-110 m sur les deux instantanes publies, rayon exterieur c t + 130 m ; "
               "energie elastique des t ~ 0.03, maximale a t ~ 0.16) ; avec les y_c = 500 m du texte (variante), rien "
               "n'atteint l'interface avant t ~ 0.28.")


def fig14(ctx):
    """Fig. 14 — repartition de l'energie, granite-eau (l = 6) et granite-eau realiste"""
    item = "Fig. 14"
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    rows = []
    for ax, cases, scale, lab in ((axs[0], ["m2an_fig14_granite_eau"], None, None),
                                  (axs[1], ["m2an_fig15_realiste", "m2an_fig15_realiste_yc500"], 1e11, "$10^{-11}$")):
        drawn = False
        for case in cases:
            r = ctx.case(case)
            if not r or not r.get("energy"):
                continue
            E0 = r["energy"]["E"][0]
            sc = scale or 10.0 ** math.floor(math.log10(E0))
            variant = case.endswith("_yc500")
            _repartition(ax, r, sc, lab or f"$10^{{{-int(math.log10(sc))}}}$", ls=":" if variant else "-",
                         lw=1.5 if variant else 2.5, suffix=" ($y_c$ = 500 m, texte)" if variant else "")
            rows.append([case, f"E0 = {E0:.4e}", f"E(T)/E0 = {r['energy']['E'][-1] / E0:.6f}"])
            drawn = True
        if not drawn:
            ax.text(0.5, 0.5, "non calcule", ha="center")
    img = ctx.save(fig, item, "M2AN25_figure_14")
    ctx.section(item, "Repartition de l'energie, cas contrastes", [img], [REALISTIC_NOTE, SOURCE_NOTE],
                (["cas", "energie initiale", "energie finale"], rows, "") if rows else None, status="compare" if rows else "")


def fig15(ctx):
    """Fig. 15 — instantanes du cas granite-eau realiste"""
    _snapshots(ctx, "Fig. 15", "m2an_fig15_realiste", "Instantanes granite-eau realiste", "M2AN25_figure_15",
               notes=[REALISTIC_NOTE, SOURCE_NOTE], pct=(99.9, 99.0),
               variants=[("m2an_fig15_realiste_yc500", "_variante_yc500", "source en (0, 500) comme le texte")])


FIGURES = [fig3, tab1, tab2, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15]
