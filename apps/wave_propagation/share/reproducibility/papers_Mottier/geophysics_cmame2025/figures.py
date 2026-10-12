#!/usr/bin/env python3
"""
figures.py — figures et tableaux reproduits de [CMAME25] (appele par run_campaign.py / plot_mottier.py).

Conventions (voir ../LISEZMOI.txt, section 4) :
  - CFL* = c dt / h_min, h diametre des cellules, dt = T / (dernier nombre de pas stable) ;
  - erreurs de couplage (Fig. 11) : Delta_kin = |(v_F - m_T).n|, Delta_dyn = |p_F n - s_T n| au capteur S^i, en mm/s et mPa ;
  - Table 4 et Fig. 13 : la reference SEM3D de l'article n'est pas disponible ; reference semi-analytique Gar6more2D
    (valable tant qu'aucune reflexion sur les bords du domaine n'atteint les capteurs, ce qui est le cas avant 0.5 s) ;
  - Table 5 et Fig. 17 : reference = HHO k = 4 sur le maillage de triangles, Delta t / 2.
"""
import math
import numpy as np
import matplotlib.pyplot as plt

import plot_mottier as pm
from plot_mottier import GREEN, ORANGE, DARKRED, NAVY

PAPER = "CMAME25"
INTRO = """
Reproduction des cas tests de la section 5 de [CMAME25] avec le prototype d'ordre un de wave_propagation (fichiers
d'entree Lua de ce dossier, lances par run_campaign.py). Chaque section : figure ou tableau reproduit, comparaison avec
les valeurs de l'article (valeurs_article.json) et ecarts.
"""
K3 = {1: ORANGE, 2: DARKRED, 3: NAVY}


def _pv(ctx, *keys):
    d = ctx.paper
    for k in keys:
        if not isinstance(d, dict) or str(k) not in d:
            return None
        d = d[str(k)]
    return d


def _mesh_file(ctx, name):
    import meshes
    for p in (ctx.paper_dir / "meshes" / name, ctx.paper_dir / "meshes" / (name + ".bz2"), ctx.out / "meshes" / name):
        if p.exists():
            return meshes.read_polygons(p)
    return None


# ------------------------------------------------------------------------------------------------------------------
def fig5(ctx):
    """Fig. 5 — maillages cartesien, simplicial et polygonal (l = 3)"""
    item = "Fig. 5"
    fig, axs = plt.subplots(1, 3, figsize=(15, 3))
    nx, ny = 16, 8
    xs, ys = np.linspace(-1, 1, nx + 1), np.linspace(0, 1, ny + 1)
    pts = np.array([[x, y] for y in ys for x in xs])
    cells = [[j * (nx + 1) + i, j * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i + 1, (j + 1) * (nx + 1) + i]
             for j in range(ny) for i in range(nx)]
    tags = [1 if pts[c].mean(axis=0)[0] < 0 else 2 for c in cells]
    cols = {1: "#b00000", 2: "#0000ff"}
    pm.draw_mesh(axs[0], pts, cells, tags, cols)
    for ax, name in ((axs[1], "tri_lc0.15.txt"), (axs[2], "poly_64.txt")):
        m = _mesh_file(ctx, name)
        if m:
            pm.draw_mesh(ax, *m, cols)
    img = ctx.save(fig, item, "CMAME25_figure_5")
    ctx.section(item, "Maillages cartesien, simplicial et polygonal (l = 3)", [img],
                ["Triangles gmsh (lc = 0.15) et polygones PolyMesher (64 germes par carre unite) regeneres et symetrises : "
                 "les fichiers de l'article ne sont pas disponibles."], status="dessin")


def fig6(ctx):
    """Fig. 6 — CFL* des cas purs en fonction de eta"""
    item = "Fig. 6"
    etas = [2.0 ** w for w in range(-3, 4)]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.5))
    rows = []
    for ax, medium, title in ((axs[0], "acoustic", "Purely Acoustic Case"), (axs[1], "elastic", "Purely Elasticity Case")):
        for k in (1, 2, 3):
            ys = []
            for w in range(-3, 4):
                r = ctx.case(f"cmame_fig6_{medium}_k{k}_w{w:+d}")
                ys.append(r["data"]["cfl_h_min"][0] if r else None)
            ax.plot([e for e, y in zip(etas, ys) if y], [y for y in ys if y], "d", color=K3[k], ms=10, label=f"k = {k}")
            pv = _pv(ctx, item, medium, k)
            if pv:
                ax.plot(etas, pv, "d", mfc="none", mec="k", ms=13, mew=0.8)
                rows.append([medium, k, pm.fmt_ratio(pm.ratio_stats(ys, pv))])
        pm.loglog(ax, r"$\eta$", "CFL*"); ax.set_title(title)
        ax.set_xticks(etas); ax.set_xticklabels([f"$2^{{{w}}}$" for w in range(-3, 4)]); pm.legend(ax)
    img = ctx.save(fig, item, "CMAME25_figure_6")
    ctx.section(item, "CFL* des cas purs en fonction de eta", [img],
                ["ERK(4), 16 x 8 carres (niveau 2 du code de l'article : reproduit les valeurs publiees ; le niveau 4 des "
                 "tableaux donne des CFL* plus faibles), seuil de 5 % ; losanges vides : valeurs relevees sur l'article."],
                (["milieu", "k", "nous / article : mediane [min, max] (points)"], rows, "") if rows else None,
                status="reproduit" if rows else "")


def _cfl_table(ctx, item, geoms, paper_fn, title, name):
    head = ["k"] + [f"{g} s={s}" for g in geoms for s in (2, 3, 4)]
    rows, md, ratios = [], [], []
    for k in (1, 2, 3):
        line, mline = [k], [k]
        for g in geoms:
            for i, s in enumerate((2, 3, 4)):
                r = ctx.case(f"cmame_cfl_{g}_k{k}_s{s}")
                v = r["data"]["cfl_h_min"][0] if r else None
                p = paper_fn(g, k, i)
                line.append(f"{v:.3f}" if v else "-")
                mline.append(f"{v:.3f} / {p:.3f}" if v and p else "-")
                if v and p:
                    ratios.append(v / p)
        rows.append(line); md.append(mline)
    imgs = [pm.table_image(ctx, item, name, head, rows, title)]
    st = (float(np.median(ratios)), min(ratios), max(ratios), len(ratios)) if ratios else None
    return imgs, (head, md, "nous / article"), st


def tab2(ctx):
    """Table 2 — CFL* des schemas ERK(s), maillages cartesiens"""
    item = "Table 2"
    imgs, tab, st = _cfl_table(ctx, item, ("carres",), lambda g, k, i: (_pv(ctx, item, k) or [None] * 3)[i],
                               "CFL* (carres, 64 x 32)", "CMAME25_tableau_2")
    ctx.section(item, "CFL* des schemas ERK(s)", imgs,
                [f"64 x 32 carres (niveau 4 du code de l'article), poids 0.8 / 1.5 ; nous / article : {pm.fmt_ratio(st)}."],
                tab, status="reproduit" if st and abs(st[0] - 1) < 0.02 else "compare")


def tab3(ctx):
    """Table 3 — CFL* selon la forme des cellules"""
    item = "Table 3"
    key = {"triangles": "triangles", "carres": "squares", "polygones": "hexagons"}
    imgs, tab, st = _cfl_table(ctx, item, ("triangles", "carres", "polygones"),
                               lambda g, k, i: (_pv(ctx, item, key[g], k) or [None] * 3)[i], "CFL* selon la geometrie",
                               "CMAME25_tableau_3")
    ctx.section(item, "CFL* selon la forme des cellules", imgs,
                [f"Triangles gmsh lc = 0.035 et polygones PolyMesher a 1024 germes par carre unite regeneres ; nous / "
                 f"article : {pm.fmt_ratio(st)}."], tab, status="compare" if st else "")


# ------------------------------------------------------------------------------------------------------------------
# efficacite : erreur dG (convention des figures de [M2AN25] : (e_m^2 + e_s^2 / 2)^1/2) et CPU de chaque niveau,
# normalises par un point de reference comme dans l'article
MARK = {"erk2": "o", "erk3": "*", "erk4": "d", "sdirk1": "o", "sdirk2": "*", "sdirk3": "d",
        "triangles": "o", "polygones": "*", "carres": "d"}


def _series(ctx, name):
    r = ctx.case(name)
    if not r or not r.get("data"):
        return None
    d = r["data"]
    err = [math.sqrt(m * m + 0.5 * s_ * s_) for m, s_ in zip(d["dg_m"], d["dg_s_plain"])]
    return list(d["cpu"]), err


def _eff_panel(ctx, ax, series, ref, pv=None):
    """series : [(nom du cas, cle de marqueur, k, legende)] ; ref : nom du cas de normalisation (premier niveau)"""
    r0 = _series(ctx, ref)
    rows = []
    for name, *_ in ([] if r0 else series):             # serie de reference absente : premiere serie calculee
        r0 = _series(ctx, name)
        if r0:
            rows.append([f"(normalisation par {name}, {ref} absent)", "", ""])
            break
    if not r0:
        return rows
    c0, e0 = r0[0][0], r0[1][0]
    for name, mk, k, lab in series:
        sr = _series(ctx, name)
        if not sr:
            continue
        cpu, err = sr
        ax.plot(np.array(cpu) / c0, np.array(err) / e0, MARK[mk], color=K3[k], ms=9 if mk != "*" else 12, ls="none")
        rows.append([lab, " ".join(f"{c / c0:.3g}" for c in cpu), " ".join(f"{e / e0:.3g}" for e in err)])
    for key, groups in (pv or {}).items():
        kk = int(key[1:])
        for g in groups:
            ax.plot(g["cpu"], [g["error"]] * 2, "-", color=K3[kk], alpha=0.35, lw=6, solid_capstyle="round")
    pm.loglog(ax, "Relative CPU time", r"Relative $\|\cdot\|_{dG}$")
    return rows


def _eff_legend(ax, marks, labels):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color="k", marker=MARK[m], ls="none", ms=9, label=l) for m, l in zip(marks, labels)]
    h += [Line2D([], [], color=K3[k], lw=3, label=f"k = {k}") for k in (1, 2, 3)]
    h += [Line2D([], [], color="0.6", lw=6, alpha=0.4, label="article (groupes releves)")]
    ax.legend(handles=h, loc="lower left", fontsize=8, ncol=2)


EFF_NOTE = ("Niveaux l = 1..5 pour tous les k ; temps de calcul (horloge) mesure sur un fil, les cas d'efficacite etant "
            "lances un par un, seuls sur la machine ; erreur dG et CPU normalises "
            "par le premier niveau de la serie de reference, comme dans l'article. Bandes grises : groupes de marqueurs "
            "releves sur la figure publiee (marqueurs souvent superposes). Les niveaux, dt_0 et Tol_0 de l'article ne sont "
            "pas donnes : etude refaite selon la procedure decrite (voir l'entree Lua).")


def fig7(ctx):
    """Fig. 7 — efficacite ERK(s) et SDIRK(s,s+1), erreur spatiale dominante"""
    item = "Fig. 7"
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.5))
    rows = []
    rows += _eff_panel(ctx, axs[0], [(f"cmame_eff_espace_erk{s}_k{k}", f"erk{s}", k, f"ERK({s}) k={k}") for k in (1, 2, 3)
                                     for s in (2, 3, 4)], "cmame_eff_espace_erk4_k1", _pv(ctx, item, "gauche"))
    rows += _eff_panel(ctx, axs[1], [(f"cmame_eff_espace_sdirk{s}_k{k}", f"sdirk{s}", k, f"SDIRK({s},{s + 1}) k={k}")
                                     for k in (1, 2, 3) for s in (1, 2, 3)], "cmame_eff_espace_sdirk3_k2",
                       _pv(ctx, item, "droite"))
    _eff_legend(axs[0], ("erk2", "erk3", "erk4"), ("ERK(2)", "ERK(3)", "ERK(4)"))
    _eff_legend(axs[1], ("sdirk1", "sdirk2", "sdirk3"), ("SDIRK(1,2)", "SDIRK(2,3)", "SDIRK(3,4)"))
    img = ctx.save(fig, item, "CMAME25_figure_7")
    ctx.section(item, "Efficacite des schemas explicites et implicites (erreur spatiale dominante)", [img], [EFF_NOTE],
                (["serie", "CPU relatif par niveau", "erreur relative par niveau"], rows, "") if rows else None,
                status="reconstitue" if rows else "")


def fig8(ctx):
    """Fig. 8 — ERK(2) contre SDIRK(3,4) ; ERK(4) selon la forme des cellules"""
    item = "Fig. 8"
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.5))
    rows = []
    rows += _eff_panel(ctx, axs[0], [(f"cmame_eff_espace_erk2_k{k}", "erk2", k, f"ERK(2) k={k}") for k in (2, 3)]
                       + [(f"cmame_eff_espace_sdirk3_k{k}", "sdirk3", k, f"SDIRK(3,4) k={k}") for k in (2, 3)],
                       "cmame_eff_espace_erk2_k2", _pv(ctx, item, "gauche"))
    rows += _eff_panel(ctx, axs[1], [(f"cmame_eff_espace_{g}_k{k}", g, k, f"{g} k={k}") for k in (1, 2, 3)
                                     for g in ("triangles", "polygones")]
                       + [(f"cmame_eff_espace_erk4_k{k}", "carres", k, f"carres k={k}") for k in (1, 2, 3)],
                       "cmame_eff_espace_triangles_k1", _pv(ctx, item, "droite"))
    _eff_legend(axs[0], ("erk2", "sdirk3"), ("ERK(2)", "SDIRK(3,4)"))
    _eff_legend(axs[1], ("triangles", "polygones", "carres"), ("Simplices", "Polygons", "Squares"))
    img = ctx.save(fig, item, "CMAME25_figure_8")
    ctx.section(item, "ERK(2) contre SDIRK(3,4) ; ERK(4) selon la forme des cellules", [img], [EFF_NOTE],
                (["serie", "CPU relatif par niveau", "erreur relative par niveau"], rows, "") if rows else None,
                status="reconstitue" if rows else "")


def fig9(ctx):
    """Fig. 9 — efficacite ERK(s) et ERK(4) selon la forme des cellules, erreur temporelle dominante"""
    item = "Fig. 9"
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.5))
    rows = []
    rows += _eff_panel(ctx, axs[0], [(f"cmame_eff_temps_erk{s}_k{k}", f"erk{s}", k, f"ERK({s}) k={k}") for k in (1, 2, 3)
                                     for s in (2, 3, 4)], "cmame_eff_temps_erk2_k1", _pv(ctx, item, "gauche"))
    rows += _eff_panel(ctx, axs[1], [(f"cmame_eff_temps_{g}_k{k}", g, k, f"{g} k={k}") for k in (1, 2, 3)
                                     for g in ("triangles", "polygones")]
                       + [(f"cmame_eff_temps_erk4_k{k}", "carres", k, f"carres k={k}") for k in (1, 2, 3)],
                       "cmame_eff_temps_triangles_k1", _pv(ctx, item, "droite"))
    _eff_legend(axs[0], ("erk2", "erk3", "erk4"), ("ERK(2)", "ERK(3)", "ERK(4)"))
    _eff_legend(axs[1], ("triangles", "polygones", "carres"), ("Simplices", "Polygons", "Squares"))
    img = ctx.save(fig, item, "CMAME25_figure_9")
    ctx.section(item, "Efficacite, erreur temporelle dominante", [img], [EFF_NOTE],
                (["serie", "CPU relatif par niveau", "erreur relative par niveau"], rows, "") if rows else None,
                status="reconstitue" if rows else "")


# ------------------------------------------------------------------------------------------------------------------
def fig10(ctx):
    """Fig. 10 — domaine des cas tests de l'ondelette de Ricker"""
    item = "Fig. 10"
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 0.5, color="#f2c9c9"))
    ax.add_patch(plt.Rectangle((-0.5, 0.0), 1, 0.5, color="#c9c9f2"))
    ax.plot([-0.5, 0.5], [0, 0], "k", lw=2); ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fill=False, lw=2))
    ax.plot(0.0, 0.125, "o", color="#c000c0", ms=8)
    ax.plot(-0.15, 0.0, "o", color="w", mec="k", ms=8); ax.text(-0.2, -0.06, r"$S^i$", fontsize=12)
    ax.text(-0.05, 0.3, r"$\Omega^f$", color="b", fontsize=16); ax.text(-0.05, -0.35, r"$\Omega^s$", color="#b00000", fontsize=16)
    ax.set_aspect("equal"); ax.axis("off")
    img = ctx.save(fig, item, "CMAME25_figure_10")
    ctx.section(item, "Domaine de l'ondelette de Ricker", [img],
                ["Cas academique de la Fig. 11 : capteur d'interface S^i = (-0.15, 0) (croquis de l'article ; le texte "
                 "donne (-0.3, 0))."], status="dessin")


def fig11(ctx):
    """Fig. 11 — erreurs sur les conditions de couplage (k = 1, 3 ; l = 4..7)"""
    item = "Fig. 11"
    cols = {4: GREEN, 5: ORANGE, 6: DARKRED, 7: NAVY}
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    rows = []
    for i, k in enumerate((1, 3)):
        for j, (q, lab) in enumerate((("kin", r"$\Delta_{KIN}$ [mm/s]"), ("dyn", r"$\Delta_{DYN}$ [mPa]"))):
            ax = axs[i, j]
            for l in (4, 5, 6, 7):
                r = ctx.case(f"cmame_fig11_k{k}_l{l}")
                if not r or not r.get("sensors") or "Si" not in r["sensors"]:
                    continue
                s = r["sensors"]["Si"]
                t = np.array(s["t"])
                if q == "kin":
                    v = 1e3 * np.abs(np.array(s["vFy"]) - np.array(s["my"]))
                else:
                    v = 1e3 * np.hypot(np.array(s["sxy"]), np.array(s["pF"]) - np.array(s["syy"]))
                ax.plot(t, v, color=cols[l], lw=1.8, label=f"l = {l}")
                pv = _pv(ctx, item, f"{q}_k{k}", f"l{l}")
                pmax = max(pv["v"]) if pv and pv["v"] else None
                if pv:
                    ax.plot(pv["t"], pv["v"], ":", color=cols[l], lw=1.0)
                rows.append([k, l, q, f"{v.max():.3g}", f"{pmax:.3g}" if pmax else "-"])
            ax.set_xlabel("t [s]"); ax.set_ylabel(lab); ax.set_title(f"k = {k}"); pm.legend(ax, ncol=2)
    img = ctx.save(fig, item, "CMAME25_figure_11")
    ctx.section(item, "Erreurs sur les conditions de couplage", [img],
                ["S^i = (-0.15, 0), theta = 1 ; pointilles : courbes relevees sur l'article. Le maximum sur [0, 1] est compare "
                 "dans le tableau. Seule la courbe cinematique k = 1, l = 4 concorde ; nos erreurs decroissent "
                 "regulierement avec l, les publiees stagnent pour k = 3 (le pas de temps n'y est pour rien : Delta t "
                 "8 fois plus grand change nos maxima de moins de 15 %). Avec la quadrature d'interface du code de "
                 "l'article (space.interface_quadrature = 'legacy'), k = 3, l = 5 : dynamique 0.174 (article 0.214, "
                 "quadrature exacte 0.067), cinematique 0.092 (article 0.022, exacte 0.056) : ecart non explique."], (["k", "l", "erreur", "max nous", "max article"], rows, "") if rows else None,
                status="compare" if rows else "")


def _ref_error(r, sensor, cols, ref):
    s = r["sensors"].get(sensor) if r.get("sensors") else None
    if not s or ref is None:
        return None
    t = np.array(s["t"])
    num = den = 0.0
    for c in cols:
        g = np.interp(t, ref["t"], ref[c])
        num += np.sum((np.array(s[c]) - g) ** 2); den += np.sum(g ** 2)
    return math.sqrt(num / den) if den > 0 else None


def tab4(ctx):
    """Table 4 — pas de temps, CPU et erreurs, granite-eau realiste"""
    item = "Table 4"
    rf, rs = ctx.reference("cmame_granite_eau", "Sf"), ctx.reference("cmame_granite_eau", "Ss")
    head = ["schema", "solveur", "k", "dt~", "CPU [s]", "Err^f", "Err^s", "art. dt~", "art. CPU", "art. Err^f",
            "art. Err^s"]
    rows = []
    pv = _pv(ctx, item) or {}
    h = math.hypot(3800.0 / 114, 2600.0 / 78)
    for scheme, solver, key in (("SDIRK(3,4)", "direct", "sdirk_direct"), ("SDIRK(3,4)", "iteratif", "sdirk_iterative"),
                                ("ERK(4)", "-", "erk4")):
        for k in (1, 3):
            name = f"cmame_tab4_sdirk_{solver}_k{k}" if scheme.startswith("SDIRK") else f"cmame_tab4_erk4_k{k}"
            r = ctx.case(name)
            p = pv.get(key, {}).get(str(k), {})
            if r:
                dt = r["data"]["dt"][0]
                cpu = r["data"]["cpu"][0]
                ef = _ref_error(r, "Sf", ("p",), rf)
                es = _ref_error(r, "Ss", ("vx", "vy"), rs)
                rows.append([scheme, solver, k, f"{6000.0 * dt / h:.3f}", f"{cpu:.0f}", f"{ef:.2e}" if ef else "-",
                             f"{es:.2e}" if es else "-", p.get("dt", "-"), p.get("cpu", "-"), p.get("err_f", "-"),
                             p.get("err_s", "-")])
    imgs = [pm.table_image(ctx, item, "CMAME25_tableau_4", head, rows, "Granite-eau realiste", col_width=1.0)] if rows else []
    ctx.section(item, "Pas de temps, CPU et erreurs (granite-eau realiste)", imgs,
                ["dt~ = c_p^s dt / h, h = 47.1 m (diametre des carres de 33.3 m) ; CPU : temps (horloge) sur un fil de "
                 "l'assemblage et de l'integration, cas lances un par un et seuls sur la machine, comme ceux des Figs. 7 a 9 ; "
                 "erreurs relatives l2 en temps aux capteurs S^f (pression) et S^s (vitesse) par rapport a la reference "
                 "Gar6more2D (la reference SEM3D de l'article n'est pas disponible)."],
                (head, rows, "") if rows else None, status="compare" if rows else "")


def _snapshots(ctx, item, case, title, name, cols=2, size=(6.0, 4.6)):
    r = ctx.case(case)
    if not r or not r.get("snapshots"):
        ctx.section(item, title, [], ["Cas non disponible."], status="non calcule")
        return
    d = pm.Path(r["_dir"])
    files = r["snapshots"]
    n = len(files)
    nrow = (n + cols - 1) // cols
    fig, axs = plt.subplots(nrow, cols, figsize=(size[0] * cols, size[1] * nrow), squeeze=False)
    for ax, f in zip(axs.flat, files):
        cells, vals, t = pm.read_snapshot(d / f"{case}_cells.dat", d / f)
        out = pm.draw_snapshot(ax, cells, vals, cmap_a="RdBu_r", cmap_e="RdBu_r")
        if 1 in out:
            fig.colorbar(out[1], ax=ax, location="top", shrink=0.7, label="p")
        if 0 in out:
            fig.colorbar(out[0], ax=ax, location="bottom", shrink=0.7, label="|v|")
        ax.set_title(f"t = {t:.3g} s", fontsize=10)
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    img = ctx.save(fig, item, name)
    ctx.section(item, title, [img], ["Pression (fluide) et norme de la vitesse (solide), polynomes de cellule aux sommets et "
                                     "au barycentre."], status="reproduit (qualitatif)")


def fig12(ctx):
    """Fig. 12 — instantanes, granite-eau realiste (t = 0.36 et 0.5 s)"""
    _snapshots(ctx, "Fig. 12", "cmame_tab4_sdirk_direct_k3", "Instantanes granite-eau realiste (SDIRK(3,4), k = 3)",
               "CMAME25_figure_12")


def _fit(pub, ref, q, sc):
    """Facteur a minimisant || a pub - ref || aux instants releves (ref en unites de la figure) et ecart relatif restant"""
    if not pub or not ref:
        return None
    t, v = np.array(pub["t"]), np.array(pub["v"])
    g = np.interp(t, ref["t"], np.array(ref[q]) * sc)
    a = float(np.dot(v, g) / np.dot(v, v))
    return a, float(np.linalg.norm(a * v - g) / np.linalg.norm(g))


def fig13(ctx):
    """Fig. 13 — traces aux capteurs, granite-eau realiste"""
    item = "Fig. 13"
    quantities = (("Sf", "p", 1.0, "P [Pa]", (0.2, 0.5), "Sf_p"), ("Si", "vx", 1e9, r"$V_x$ [nm/s]", (0.1, 0.45), "Si_vx"),
                  ("Ss", "vx", 1e9, r"$V_x$ [nm/s]", (0.2, 0.5), "Ss_vx"))
    imgs, rows = [], []
    for variant in (False, True):
        fig, axs = plt.subplots(3, 2, figsize=(13, 13))
        for j, (col, label) in enumerate((("implicit", "SDIRK(3,4)"), ("explicit", "ERK(4)"))):
            for i, (sensor, q, sc, ylab, win, key) in enumerate(quantities):
                ax = axs[i, j]
                ref = ctx.reference("cmame_granite_eau", sensor)
                if ref:
                    ax.plot(ref["t"], np.array(ref[q]) * sc, "k", lw=2.5, label="Reference (Gar6more2D)")
                for k, style, c in ((1, ":", DARKRED), (3, "--", NAVY)):
                    name = f"cmame_tab4_sdirk_direct_k{k}" if col == "implicit" else f"cmame_tab4_erk4_k{k}"
                    r = ctx.case(name)
                    field = "vFx" if sensor == "Si" else q
                    if r and r.get("sensors") and sensor in r["sensors"]:
                        s = r["sensors"][sensor]
                        ax.plot(s["t"], np.array(s[field]) * sc, style, color=c, lw=2, label=f"k = {k}")
                pv = _pv(ctx, item, col, key)
                pub = pv.get("reference") if pv else None
                fit = _fit(pub, ref, q, sc)
                if pub:
                    a = fit[0] if (variant and fit) else 1.0
                    ax.plot(pub["t"], a * np.array(pub["v"]), ".", color="0.6", ms=2,
                            label="article (releve)" + (f" x {a:.3g}" if variant else ""))
                if fit and not variant:
                    rows.append([label, sensor, q, f"{fit[0]:.3f}", f"{fit[1]:.2e}"])
                ax.set_xlim(*win); ax.set_xlabel("t [s]"); ax.set_ylabel(ylab)
                if i == 0 or variant:
                    if i == 0:
                        ax.set_title(label)
                    pm.legend(ax, fontsize=8)
        if variant:
            if not rows:
                plt.close(fig); continue
            fig.subplots_adjust(top=0.95)
            fig.suptitle("Variante : courbes publiees multipliees par le facteur ajuste sur la reference", fontsize=10, y=0.975)
        imgs.append(ctx.save(fig, item, "CMAME25_figure_13" + ("_variante_article_mis_a_l_echelle" if variant else "")))
    notes = ["Points gris : reference SEM3D relevee sur la figure publiee ; trait noir : Gar6more2D (capteur d'interface "
             "evalue cote solide, 0.2 m sous l'interface) ; S^i : vitesse de l'inconnue de face elastique. Facteur a "
             "minimisant ||a x courbe publiee - Gar6more2D|| (tableau) : la pression en S^f concorde (a = 1) ; la courbe "
             "publiee en S^i vaut le dixieme de la vitesse physique (a = 10 : vraisemblablement un facteur 10^8 au lieu "
             "de 10^9 pour le passage en nm/s) et celle de S^s est de signe oppose (a = -1 : V_x au capteur symetrique "
             "(1000, -500) ou -V_x). Variante : courbes publiees multipliees par a, superposees a la reference."]
    ctx.section(item, "Traces aux capteurs, granite-eau realiste", imgs, notes,
                (["colonne", "capteur", "champ", "facteur a", "ecart relatif apres mise a l'echelle"], rows, "")
                if rows else None, status="compare ; S^i publie a 1/10, S^s de signe oppose" if rows else "compare")


# ------------------------------------------------------------------------------------------------------------------
BASIN_MESHES = ("triangles", "quadrangles", "hybride", "non_conforme")
BASIN_COLORS = {1: "#40c0c0", 2: "#000060", 3: "#c8c000"}


def fig14(ctx):
    """Fig. 14 — structure du bassin (maillage de triangles) et materiaux"""
    item = "Fig. 14"
    m = _mesh_file(ctx, "bassin_triangles.txt")
    fig, ax = plt.subplots(figsize=(10, 5))
    if m:
        from matplotlib.collections import PolyCollection
        pts, cells, tags = m
        polys = [pts[c] for c in cells]
        ax.add_collection(PolyCollection(polys, facecolors=[BASIN_COLORS.get(t, "w") for t in tags], edgecolors="0.3",
                                         linewidths=0.1))
        ax.set_xlim(-2500, 2500); ax.set_ylim(-1500, 2000); ax.set_aspect("equal")
    for (x, y, n) in ((-650, 1300, r"$S^f$"), (-750, 1080, r"$S^s$")):
        ax.plot(x, y, "o", color="w", mec="k", ms=7); ax.text(x + 60, y + 60, n)
    ax.plot(563, 233, "o", color="#8000c0", ms=10)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    img = ctx.save(fig, item, "CMAME25_figure_14")
    rows = [["atmosphere", "1.225", "343", "-"], ["sediments", "1300", "1600", "900"], ["socle", "2570", "5350", "3009"]]
    ctx.section(item, "Structure du bassin et materiaux", [img],
                ["Maillage de triangles de la these (gmsh), domaine x dans (-4000, 3500) m, y dans (-6000, 2500) m (l'article "
                 "decrit (-4, 4) x (-6, 1.4) km, interface vers 1.1 km)."],
                (["materiau", "rho [kg/m3]", "c_p [m/s]", "c_s [m/s]"], rows, ""), status="dessin")


def fig15(ctx):
    """Fig. 15 — zooms des maillages hybride et non conforme"""
    item = "Fig. 15"
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    for ax, name in ((axs[0], "bassin_hybride.txt"), (axs[1], "bassin_non_conforme.txt")):
        m = _mesh_file(ctx, name)
        if not m:
            continue
        from matplotlib.collections import PolyCollection
        pts, cells, tags = m
        polys = [pts[c] for c in cells]
        ax.add_collection(PolyCollection(polys, facecolors=[BASIN_COLORS.get(t, "w") for t in tags], edgecolors="k",
                                         linewidths=0.3))
        ax.set_xlim(-1400, -700); ax.set_ylim(900, 1250); ax.set_aspect("equal"); ax.set_title(name.replace(".txt", ""))
    img = ctx.save(fig, item, "CMAME25_figure_15")
    ctx.section(item, "Zooms des maillages hybride et non conforme", [img],
                ["Non conforme : maillage de quadrangles de la these raffine une fois dans l'atmosphere et le bassin "
                 "(noeuds pendants ; cellules voisines du socle traitees comme des pentagones). L'article part d'un maillage "
                 "plus grossier (aretes de 28 m) : ici le maillage non conforme a plus de cellules que celui de quadrangles."],
                status="dessin")


def tab5(ctx):
    """Table 5 — inconnues et erreurs aux capteurs, bassin"""
    item = "Table 5"
    ref = ctx.case("cmame_bassin_reference")
    pv = _pv(ctx, item) or {}
    head = ["maillage", "cellules", "inconnues", "condensees", "Err^s", "Err^f", "art. inconnues", "art. condensees",
            "art. Err^s", "art. Err^f"]
    rows = []
    for m in BASIN_MESHES:
        r = ctx.case(f"cmame_bassin_{m}")
        if not r:
            continue
        d = r["data"]
        es = ef = None
        if ref and ref.get("sensors"):
            es = _ref_error(r, "Ss", ("vx", "vy"), ref["sensors"].get("Ss"))
            ef = _ref_error(r, "Sf", ("mx", "my"), ref["sensors"].get("Sf"))
        mesh = _mesh_file(ctx, f"bassin_{m}.txt")
        n_cells = len(mesh[1]) if mesh else "-"
        rows.append([m, n_cells, f"{d['dofs'][0]:.0f}", f"{d['face_dofs'][0]:.0f}", f"{es:.2e}" if es else "-",
                     f"{ef:.2e}" if ef else "-", pv.get("dofs", {}).get(m, "-"), pv.get("dofs_condenses", {}).get(m, "-"),
                     pv.get("err_s", {}).get(m, "-"), pv.get("err_f", {}).get(m, "-")])
    imgs = [pm.table_image(ctx, item, "CMAME25_tableau_5", head, rows, "Bassin, k = 3", col_width=1.3)] if rows else []
    ctx.section(item, "Inconnues et erreurs aux capteurs (bassin)", imgs,
                ["Maillages de la these (plus grossiers que ceux de l'article : 132 439 triangles contre 185 923, 71 178 "
                 "quadrangles contre 94 019, 67 952 cellules hybrides contre 93 490 ; non conforme raffine a partir des "
                 "quadrangles : 179 448 cellules contre 84 204). Erreurs relatives l2 en temps par rapport a la solution HHO "
                 "k = 4 sur les triangles avec Delta t / 2 (reference SEM3D de l'article non disponible) : vitesse du solide "
                 "en S^s, vitesse du fluide en S^f."], (head, rows, "") if rows else None, status="compare" if rows else "")


def fig16(ctx):
    """Fig. 16 — instantanes du bassin (t = 0, 0.16, 0.6, 1 s)"""
    _snapshots(ctx, "Fig. 16", "cmame_bassin_triangles", "Instantanes du bassin (triangles, k = 3)", "CMAME25_figure_16",
               size=(7.0, 4.5))


def _rel(t, v, tr, vr, win):
    """Ecart relatif l2 de (t, v) a (tr, vr) aux instants de tr dans la fenetre win"""
    tr, vr = np.array(tr), np.array(vr)
    m = (tr >= win[0]) & (tr <= win[1])
    if not m.any():
        return None
    g = np.interp(tr[m], t, v)
    return float(np.linalg.norm(g - vr[m]) / np.linalg.norm(vr[m]))


def fig17(ctx):
    """Fig. 17 — vitesse verticale aux capteurs du bassin"""
    item = "Fig. 17"
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    cols = {"triangles": GREEN, "quadrangles": ORANGE, "hybride": DARKRED, "non_conforme": NAVY}
    rows = []
    for ax, sensor, q, win, key in ((axs[0], "Sf", "my", (0.7, 1.45), "Sf_vy"), (axs[1], "Ss", "vy", (0.15, 1.25), "Ss_vy")):
        pv = _pv(ctx, item, key)
        sem = pv.get("SEM") if pv else None
        ref = ctx.case("cmame_bassin_reference")
        rs = ref["sensors"][sensor] if ref and ref.get("sensors") and sensor in ref["sensors"] else None
        for m in BASIN_MESHES + ("reference",):
            r = ctx.case(f"cmame_bassin_{m}")
            if not (r and r.get("sensors") and sensor in r["sensors"]):
                continue
            s = r["sensors"][sensor]
            if m == "reference":
                ax.plot(s["t"], s[q], "--", color="#e0e000", lw=1.5, label="reference k = 4")
            else:
                ax.plot(s["t"], s[q], color=cols[m], lw=3, label=m)
            e_sem = _rel(s["t"], s[q], sem["t"], sem["v"], win) if sem else None
            e_ref = _rel(s["t"], s[q], rs["t"], rs[q], win) if rs and m != "reference" else None
            rows.append([m, sensor, f"{e_sem:.2e}" if e_sem is not None else "-", f"{e_ref:.2e}" if e_ref is not None else "-"])
        if sem:
            ax.plot(sem["t"], sem["v"], ".", color="0.5", ms=2, label="SEM (releve)")
        ax.set_xlim(*win); ax.set_xlabel("t [s]"); ax.set_ylabel(r"$V_y$ [m/s]"); pm.legend(ax, fontsize=8, ncol=2)
        ax.set_title("S^f (atmosphere)" if sensor == "Sf" else "S^s (bassin)")
    img = ctx.save(fig, item, "CMAME25_figure_17")
    ctx.section(item, "Vitesse verticale aux capteurs du bassin", [img],
                ["Atmosphere : vitesse du fluide m_y ; bassin : vitesse du solide v_y ; points gris : reference SEM relevee "
                 "sur la figure publiee. Ecarts relatifs l2 sur la fenetre tracee, a la courbe SEM relevee et a la "
                 "reference HHO k = 4 (triangles, Delta t / 2)."],
                (["maillage", "capteur", "ecart a SEM (releve)", "ecart a la reference k = 4"], rows, "") if rows else None,
                status="compare" if rows else "non calcule")


FIGURES = [fig5, fig6, tab2, tab3, fig7, fig8, fig9, fig10, fig11, tab4, fig12, fig13, fig14, fig15, tab5, fig16, fig17]
