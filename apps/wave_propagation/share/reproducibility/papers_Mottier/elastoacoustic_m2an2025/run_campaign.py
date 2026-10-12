#!/usr/bin/env python3
"""
run_campaign.py — reproduction des cas tests de [M2AN25] R. Mottier, A. Ern, R. Khot, L. Guillot, Hybrid high-order
methods for elasto-acoustic wave propagation in the time domain, ESAIM Math. Model. Numer. Anal. 59 (2025) 2685-2715,
section 6 : analyse spectrale (Fig. 3, Tables 1 et 2), convergence en temps et en espace (Figs. 5 a 7), ondelette de
Ricker comparee a Gar6more2D (Figs. 9 a 12), energie (Figs. 13 et 14), cas realiste granite-eau (Figs. 14 et 15) ;
dessins des maillages (Fig. 4) et du domaine (Fig. 8).

Dans le build (<build>/apps/wave_propagation/reproducibility/papers_Mottier/elastoacoustic_m2an2025) :
    python3 run_campaign.py              campagne complete (cas les plus longs d'abord), figures et rapport
    python3 run_campaign.py --quick      verification rapide (niveaux et durees reduits), resultats dans rapide/
    python3 run_campaign.py --only fig10 cas dont le nom contient fig10 ;  --list : liste des cas
    python3 run_campaign.py --no-run     figures et rapport a partir des resultats existants
Voir ../LISEZMOI.txt.
"""
import math, sys
from pathlib import Path

sys.dont_write_bytecode = True                  # pas de __pycache__ dans les sources

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from campaign_mottier import Case, main          # noqa: E402
from gar6more import Reference                   # noqa: E402

P = "M2AN25"
CASES = []
SQ3 = math.sqrt(3.0)

# poids eta_* de (6.2) et poids effectifs reproduisant les Tables 1 et 2 (voir LISEZMOI)
ETA_STAR = {"equal": (0.88, 1.54), "mixed": (0.80, 1.38)}
ETA_EFF = {"equal": (0.88, 1.75), "mixed": (0.80, 1.04)}
ORDER = {"equal": "eo", "mixed": "mo"}


def spectral(cells, k):
    """(duree en min, memoire en Go) du probleme aux valeurs propres dense : n = cellules x dim P_k, memoire en n^2 et
    duree en n^3 (mesure : 64 polygones, k = 6 : 3.9 Go, ~12 min)"""
    r = cells / 64.0 * (k + 1) * (k + 2) / 56.0
    return 0.05 + 12.0 * r ** 3, 0.3 + 4.0 * r ** 2

# ------------------------------------------------------------------------------------------------------------------
# Fig. 3 : rayon spectral en fonction de eta (cas purs et couple)
MEDIA = {"acoustic": "coupled_poly_in_time", "elastic": "coupled_poly_in_time", "coupled": "coupled_poly_in_time"}
for medium in ("acoustic", "elastic", "coupled"):
    for order in ("equal", "mixed"):
        for k in (1, 2, 3):
            for w in range(-3, 4):
                eta = 2.0 ** w
                CASES.append(Case(f"m2an_fig3_{medium}_{ORDER[order]}_k{k}_w{w:+d}", P, "Fig. 3", medium, "fig3_rayon_spectral.lua",
                                  [f"problem.medium = '{medium}'", f"space.k = {k}", f"space.order = '{order}'",
                                   f"space.weight_acoustic = {eta!r}", f"space.weight_elastic = {eta!r}"],
                                  cost=spectral(32, k)[0], mem=spectral(32, k)[1],
                                  meta={"medium": medium, "order": order, "k": k, "eta": eta}))

# Table 1 : couplage, eta = 2^w eta_* (valeurs publiees) et variante (poids effectifs)
for variant, etas in (("", ETA_STAR), ("poids effectifs", ETA_EFF)):
    tag = "" if not variant else "_eff"
    for order in ("equal", "mixed"):
        for k in (1, 2, 3):
            for w in range(-3, 4):
                ef, es = etas[order]
                CASES.append(Case(f"m2an_tab1{tag}_{ORDER[order]}_k{k}_w{w:+d}", P, "Table 1", f"{order} k={k}",
                                  "tab1_rayon_spectral_couplage.lua",
                                  [f"space.k = {k}", f"space.order = '{order}'", f"space.weight_acoustic = {ef * 2.0 ** w!r}",
                                   f"space.weight_elastic = {es * 2.0 ** w!r}"],
                                  cost=spectral(32, k)[0], mem=spectral(32, k)[1], variant=variant,
                                  meta={"order": order, "k": k, "w": w}))

# Table 2 : geometrie des cellules, k = 1..6
GEOM = {"simplexes": ["tri_lc0.35.txt"], "quadrangles": None, "polygones": ["poly_32.txt"]}
CELLS = {"simplexes": 52, "quadrangles": 32, "polygones": 64}
for variant, etas in (("", ETA_STAR), ("poids effectifs", ETA_EFF)):
    tag = "" if not variant else "_eff"
    for geom, files in GEOM.items():
        for order in ("equal", "mixed"):
            for k in range(1, 7):
                ef, es = etas[order]
                st = [f"space.k = {k}", f"space.order = '{order}'", f"space.weight_acoustic = {ef!r}",
                      f"space.weight_elastic = {es!r}"]
                if files:
                    st.append("mesh.type = 'polygons'")
                CASES.append(Case(f"m2an_tab2{tag}_{geom[:4]}_{ORDER[order]}_k{k}", P, "Table 2", f"{geom} {order} k={k}",
                                  "tab2_rayon_spectral_geometrie.lua", st, meshes=files, cost=spectral(CELLS[geom], k)[0],
                                  mem=spectral(CELLS[geom], k)[1],
                                  variant=variant, meta={"geometry": geom, "order": order, "k": k}))

# ------------------------------------------------------------------------------------------------------------------
# Fig. 5 : convergence en temps (k = 4, l = 1, ordre mixte, O(1/h)), n = 4..9 ; figure publiee : ERK a gauche,
# SDIRK a droite (la legende de l'article inverse les deux panneaux). Les ERK(s) publies sont stables aux six pas de
# temps et reproduits avec la stabilisation O(1) ; avec la O(1/h) de la legende, ils explosent aux deux plus grands pas
# (variante)
for scheme, stages in (("sdirk", (1, 2, 3)), ("erk", (2, 3, 4))):
    for s in stages:
        st = [f"time.scheme = '{scheme}'", f"time.stages = {s}"]
        CASES.append(Case(f"m2an_fig5_{scheme}{s}", P, "Fig. 5", "droite" if scheme == "sdirk" else "gauche",
                          "fig5_convergence_temps.lua", st + (["space.scaling = '1'"] if scheme == "erk" else []),
                          cost=2.0, mem=0.5, quick=["study.time_levels = 4"], meta={"scheme": scheme, "s": s, "x": "dt"}))
        if scheme == "erk":
            CASES.append(Case(f"m2an_fig5_{scheme}{s}_o1h", P, "Fig. 5", "gauche", "fig5_convergence_temps.lua", st,
                              cost=2.0, mem=0.5, quick=["study.time_levels = 4"], variant="O(1/h) de la legende",
                              meta={"scheme": scheme, "s": s, "x": "dt"}))

# Figs. 6 et 7 : convergence en espace, n = 8, l = 2..5, k = 1..4
PANELS6 = [("gauche", "equal", "1", "erk"), ("centre", "mixed", "1", "erk"), ("droite", "mixed", "1/h", "sdirk")]
POLY7 = ["poly_128.txt", "poly_512.txt", "poly_2048.txt", "poly_8192.txt"]
for fig, lua, meshes in (("Fig. 6", "fig6_convergence_espace_cartesien.lua", None),
                         ("Fig. 7", "fig7_convergence_espace_polygones.lua", POLY7)):
    num = fig.split()[-1]
    for panel, order, scaling, scheme in PANELS6:
        for k in (1, 2, 3, 4):
            ef, es = ETA_STAR[order]
            st = [f"space.k = {k}", f"space.order = '{order}'", f"space.scaling = '{scaling}'",
                  f"space.weight_acoustic = {ef!r}", f"space.weight_elastic = {es!r}", f"time.scheme = '{scheme}'",
                  f"time.stages = {4 if scheme == 'erk' else 3}"]
            th = {3: 2, 4: 4}.get(k, 1) if scheme == "sdirk" else 1    # SDIRK : solveur direct (Pardiso) multifil
            if order == "mixed":                     # variante : poids elastique 1.04 (Table 1) en ordre mixte
                CASES.append(Case(f"m2an_fig{num}_{panel}_k{k}_mo104", P, fig, panel, lua,
                                  st + ["space.weight_elastic = 1.04"], meshes=meshes,
                                  cost=(3.0 if scheme == "erk" else 6.0) * k * k * (4 if meshes else 1), threads=th,
                                  mem=1.0 + 0.5 * k * k * (2 if meshes else 1), variant="poids MO 1.04",
                                  quick=["study.last_level = math.min(study.last_level, " + ("3" if meshes is None else "2") + ")"],
                                  meta={"order": order, "scaling": scaling, "scheme": scheme, "k": k, "x": "h"}))
            CASES.append(Case(f"m2an_fig{num}_{panel}_k{k}", P, fig, panel, lua, st, meshes=meshes,
                              cost=(3.0 if scheme == "erk" else 6.0) * k * k * (4 if meshes else 1), threads=th,
                              mem=1.0 + 0.5 * k * k * (2 if meshes else 1),
                              quick=["study.last_level = math.min(study.last_level, " + ("3" if meshes is None else "2") + ")"],
                              meta={"order": order, "scaling": scaling, "scheme": scheme, "k": k, "x": "h"}))

# ------------------------------------------------------------------------------------------------------------------
# Ricker : Figs. 9 a 12
AC = {"acoustic": (1.0, 1.0), "elastic": (1.0, SQ3, 1.0)}
EAU = {"acoustic": (0.5, 0.5), "elastic": (1.3, 2.0, 1.0)}
AIR = {"acoustic": (1.0 / 2.75, 1.0 / 2.75), "elastic": (800.0, 6.36, 3.27)}


def materials(m):
    a, e = m["acoustic"], m["elastic"]
    return [f"problem.acoustic = {{ rho = {a[0]!r}, vp = {a[1]!r} }}",
            f"problem.elastic = {{ rho = {e[0]!r}, vp = {e[1]!r}, vs = {e[2]!r} }}",
            f"problem.ricker.velocity = {a[1]!r}"]


DT9 = 0.1 * 2.0 ** -9
# Fig. 9 : instantanes, k = 1, l = 7
CASES.append(Case("m2an_fig9_instantanes", P, "Fig. 9", "k=1 l=7", "fig9_10_ricker_academique.lua",
                  ["space.k = 1", "mesh.level = 7", "time.final_time = 0.32", f"time.n_steps = {round(0.32 / DT9)}",
                   "output.snapshots = { 0.0, 0.25, 0.27, 0.32 }"], cost=25.0, threads=4, mem=6.0,
                  quick=["mesh.level = 5"], meta={"times": [0.0, 0.25, 0.27, 0.32]}))
# Fig. 10 : academique, l = 4, 5 ; Figs. 11, 12 : contrastes, l = 5, 6 ; k = 1..4
for fig, lua, mats, levels, tf, name in (("Fig. 10", "fig9_10_ricker_academique.lua", AC, (4, 5), 0.25, "academique"),
                                         ("Fig. 11", "fig11_12_ricker_contraste.lua", EAU, (5, 6), 0.5, "granite_eau"),
                                         ("Fig. 12", "fig11_12_ricker_contraste.lua", AIR, (5, 6), 0.5, "granite_air")):
    num = fig.split()[-1]
    for l in levels:
        for k in (1, 2, 3, 4):
            st = materials(mats) + [f"space.k = {k}", f"mesh.level = {l}", f"time.final_time = {tf!r}",
                                    f"time.n_steps = {round(tf / DT9)}"]
            CASES.append(Case(f"m2an_fig{num}_l{l}_k{k}", P, fig, f"l={l} k={k}", lua, st,
                              cost=0.3 * 4 ** (l - 4) * k * k * tf / 0.25, threads=2 if l >= 6 else 1,
                              mem=0.5 + 0.02 * 4 ** (l - 4) * (k + 1) ** 2,     # l = 6 : 0.9 / 4.3 / 7.7 Go (k = 1, 3, 4)
                              quick=[f"mesh.level = {l - 2}"], meta={"level": l, "k": k, "reference": f"m2an_{name}"}))

# ------------------------------------------------------------------------------------------------------------------
# Fig. 13 : energie (academique) ; Fig. 14 : granite-eau
CASES.append(Case("m2an_fig13_repartition", P, "Fig. 13", "gauche", "fig13_14_energie.lua",
                  ["time.final_time = 10.0", f"time.n_steps = {round(10.0 / DT9)}"], cost=400.0, threads=6, mem=7.0,
                  quick=["mesh.level = 4", "time.final_time = 1.0", f"time.n_steps = {round(1.0 / DT9)}"],
                  meta={"panel": "repartition"}))
for l in (4, 5, 6):
    for k in (1, 2, 3):
        CASES.append(Case(f"m2an_fig13_perte_l{l}_k{k}", P, "Fig. 13", f"droite l={l} k={k}", "fig13_14_energie.lua",
                          [f"space.k = {k}", f"mesh.level = {l}"], cost=1.0 * 4 ** (l - 4) * k * k, threads=2 if l >= 6 else 1,
                          mem=0.5 + 0.02 * 4 ** (l - 4) * (k + 1) ** 2,
                          quick=[f"mesh.level = {l - 2}"], meta={"panel": "perte", "level": l, "k": k}))
CASES.append(Case("m2an_fig14_granite_eau", P, "Fig. 14", "gauche", "fig13_14_energie.lua",
                  materials(EAU) + ["time.final_time = 10.0", f"time.n_steps = {round(10.0 / DT9)}"], cost=400.0, threads=6,
                  mem=7.0,
                  quick=["mesh.level = 4", "time.final_time = 1.0", f"time.n_steps = {round(1.0 / DT9)}"],
                  meta={"panel": "granite_eau"}))
# Figs. 14 (droite) et 15 : granite-eau realiste, 280 x 280 cellules (560 x 560 dans l'article : memoire) ; source a
# 100 m de l'interface comme les figures publiees, variante avec les 500 m du texte
for tag, st, var in (("", [], ""), ("_yc500", ["problem.ricker.center = {0.0, 500.0}"], "centre (0, 500) du texte")):
    CASES.append(Case(f"m2an_fig15_realiste{tag}", P, "Fig. 15", "k=3", "fig14_15_granite_eau_realiste.lua", st,
                      cost=540.0, threads=8, mem=26.0, quick=["mesh.level = 0", "time.n_steps = 640"], variant=var,
                      meta={"times": [0.1275, 0.3825]}))

# ------------------------------------------------------------------------------------------------------------------
# references Gar6more2D
# retard du signal source de Gar6more2D : 1.5 / f0, comme les references publiees (petit defaut de l'onde conique en
# t = 0.24 pour granite-eau) ; avec 2 / f0, ce defaut tombe en t = 0.19 et fausse le premier creux de la pression (15 %)
REFERENCES = [
    Reference("m2an_academique", AC["acoustic"], AC["elastic"], (0.0, 0.125), 10.0, "acoustic_flux", 10.0,
              [("Sf", -0.15, 0.1), ("Ss", -0.15, -0.1)], 0.25, DT9, delay=1.5),
    Reference("m2an_granite_eau", EAU["acoustic"], EAU["elastic"], (0.0, 0.125), 10.0, "acoustic_flux", 10.0,
              [("Sf", -0.1, 0.05), ("Ss", -0.1, -0.05)], 0.5, DT9, delay=1.5),
    Reference("m2an_granite_air", AIR["acoustic"], AIR["elastic"], (0.0, 0.125), 10.0, "acoustic_flux", 10.0,
              [("Sf", -0.1, 0.05), ("Ss", -0.1, -0.05)], 0.5, DT9, delay=1.5),
]

if __name__ == "__main__":
    import figures
    if len(sys.argv) > 2 and sys.argv[1] == "--plot-only":
        import plot_mottier
        sys.exit(plot_mottier.make_all(HERE, Path(sys.argv[2]).resolve(), figures))
    sys.exit(main(CASES, HERE, __doc__, REFERENCES, figures))
