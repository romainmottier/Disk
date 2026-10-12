#!/usr/bin/env python3
"""
run_campaign.py — reproduction des cas tests de [CMAME25] R. Mottier, A. Ern, L. Guillot, Elasto-acoustic wave
propagation in geophysical media using hybrid high-order methods on general meshes, Comput. Methods Appl. Mech. Engrg.
(2025), section 5 : CFL des schemas ERK (Fig. 6, Tables 2 et 3), efficacite explicite / implicite (Figs. 7 a 9),
ondelette de Ricker (erreurs de couplage, Fig. 11 ; granite-eau realiste, Table 4, Figs. 12 et 13), bassin
sedimentaire (Table 5, Figs. 16 et 17) ; dessins des maillages et des domaines (Figs. 5, 10, 14, 15).

Dans le build (<build>/apps/wave_propagation/reproducibility/papers_Mottier/geophysics_cmame2025) :
    python3 run_campaign.py              campagne complete (cas les plus longs d'abord), figures et rapport
    python3 run_campaign.py --quick      verification rapide (niveaux et durees reduits), resultats dans rapide/
    python3 run_campaign.py --only tab4  cas dont le nom contient tab4 ;  --list : liste des cas
    python3 run_campaign.py --no-run     figures et rapport a partir des resultats existants
Les Figs. 7 a 9 et la Table 4 comparent des temps de calcul (horloge) : lancer ces cas un par un et seuls sur la
machine (environ 40 min), par exemple
    python3 run_campaign.py --only '^(?!cmame_(eff|tab4))'   puis   python3 run_campaign.py --only 'cmame_(eff|tab4)' --jobs 1
Voir ../LISEZMOI.txt.
"""
import math, sys
from pathlib import Path

sys.dont_write_bytecode = True                  # pas de __pycache__ dans les sources

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from campaign_mottier import Case, main          # noqa: E402
from gar6more import Reference                   # noqa: E402

P = "CMAME25"
CASES = []
SQ3 = math.sqrt(3.0)
CFL_TABLE2 = {(1, 2): 0.205, (1, 3): 0.253, (1, 4): 0.282, (2, 2): 0.099, (2, 3): 0.123, (2, 4): 0.138,
              (3, 2): 0.063, (3, 3): 0.079, (3, 4): 0.087}

# ------------------------------------------------------------------------------------------------------------------
# Fig. 6 : CFL* des cas purs en fonction de eta (ERK(4), ordre egal)
for medium in ("acoustic", "elastic"):
    for k in (1, 2, 3):
        for w in range(-3, 4):
            eta = 2.0 ** w
            w_st = [f"space.weight_acoustic = {eta!r}"] if medium == "acoustic" else [f"space.weight_elastic = {eta!r}"]
            CASES.append(Case(f"cmame_fig6_{medium}_k{k}_w{w:+d}", P, "Fig. 6", medium, "fig6_cfl_cas_purs.lua",
                              [f"problem.medium = '{medium}'", f"space.k = {k}", f"time.n_steps = {120 * k}"] + w_st,
                              cost=0.3 * k * k, meta={"medium": medium, "k": k, "eta": eta}))

# Tables 2 et 3 : CFL* du probleme couple, ERK(s), carres (niveau 4), triangles, polygones
QUICK_FILES = {"triangles": ["tri_lc0.15.txt"], "polygones": ["poly_64.txt"]}
for geom, files in (("carres", None), ("triangles", ["tri_lc0.035.txt"]), ("polygones", ["poly_1024.txt"])):
    item = "Table 2" if geom == "carres" else "Table 3"
    for k in (1, 2, 3):
        for s in (2, 3, 4):
            st = [f"space.k = {k}", f"time.stages = {s}", f"time.n_steps = {300 * k}"]
            if files:
                st += ["mesh.type = 'polygons'", "mesh.level = 1"]
            CASES.append(Case(f"cmame_cfl_{geom}_k{k}_s{s}", P, item, f"{geom} k={k} s={s}", "tab2_tab3_cfl.lua", st,
                              meshes=files, cost=1.5 * k * k, quick=["mesh.level = math.min(mesh.level, 2)"] if not files else [],
                              quick_meshes=QUICK_FILES.get(geom),
                              meta={"geometry": geom, "k": k, "s": s}))

# ------------------------------------------------------------------------------------------------------------------
# Figs. 7 a 9 : efficacite, niveaux l = 1..5 ; un fil par cas (CPU)
CFL_GEO = {"triangles": {1: 0.265, 2: 0.147, 3: 0.100}, "polygones": {1: 0.363, 2: 0.188, 3: 0.115}}
GEO_FILES = {"triangles": [f"tri_lc{lc}.txt" for lc in ("0.6", "0.3", "0.15", "0.075", "0.0375")],
             "polygones": [f"poly_{n}.txt" for n in (4, 16, 64, 256, 1024)]}
for regime, item, om in (("espace", "Fig. 7", ["problem.omega = 5", "problem.theta = math.sqrt(2)"]),
                         ("temps", "Fig. 9", ["problem.omega = 1", "problem.theta = 10"])):
    for k in (1, 2, 3):
        for s in (2, 3, 4):
            CASES.append(Case(f"cmame_eff_{regime}_erk{s}_k{k}", P, item, f"ERK({s}) k={k}", "fig7_9_efficacite.lua",
                              om + [f"space.k = {k}", f"time.stages = {s}", f"time.dt_h = {0.9 * CFL_TABLE2[(k, s)] / SQ3!r}"],
                              cost=2.0 * k * k, mem=1.0 + 0.3 * k * k, quick=["study.last_level = 2"],
                              meta={"regime": regime, "scheme": "erk", "s": s, "k": k, "geometry": "carres"}))
        for geom, files in GEO_FILES.items():
            CASES.append(Case(f"cmame_eff_{regime}_{geom}_k{k}", P, "Fig. 8" if regime == "espace" else "Fig. 9",
                              f"ERK(4) {geom} k={k}", "fig7_9_efficacite.lua",
                              om + [f"space.k = {k}", "time.stages = 4", f"time.dt_h = {0.9 * CFL_GEO[geom][k] / SQ3!r}",
                                    "mesh.type = 'polygons'", "mesh.level = 1", "study.last_level = 5"],
                              meshes=files, cost=3.0 * k * k, mem=1.0 + 0.3 * k * k, quick=["study.last_level = 3"],
                              meta={"regime": regime, "scheme": "erk", "s": 4, "k": k, "geometry": geom}))
    if regime == "espace":
        for k in (1, 2, 3):
            for s in (1, 2, 3):
                n0 = 16 if s == 1 else 8
                CASES.append(Case(f"cmame_eff_espace_sdirk{s}_k{k}", P, "Fig. 7", f"SDIRK({s},{s + 1}) k={k}",
                                  "fig7_9_efficacite.lua",
                                  om + [f"space.k = {k}", "space.order = 'mixed'", "space.scaling = '1/h'",
                                        "time.scheme = 'sdirk'", f"time.stages = {s}", f"time.n_steps = {n0}",
                                        f"time.steps_growth = {2.0 ** ((k + 1) / (s + 1))!r}", "solver.iterative = true",
                                        "solver.tolerance = 1e-2", f"solver.tolerance_growth = {2.0 ** -(k + 1)!r}"],
                                  cost=3.0 * k * k, mem=1.0 + 0.5 * k * k, quick=["study.last_level = 2"],
                                  meta={"regime": "espace", "scheme": "sdirk", "s": s, "k": k, "geometry": "carres"}))

# ------------------------------------------------------------------------------------------------------------------
# Fig. 11 : erreurs sur les conditions de couplage, k = 1, 3, l = 4..7 (h = 2^-l : mesh.level = l - 1)
for k in (1, 3):
    for l in (4, 5, 6, 7):
        CASES.append(Case(f"cmame_fig11_k{k}_l{l}", P, "Fig. 11", f"k={k} l={l}", "fig11_erreurs_couplage.lua",
                          [f"space.k = {k}", f"mesh.level = {l - 1}"], cost=0.4 * 4 ** (l - 4) * k * k,
                          threads=4 if l >= 6 else 1, mem=0.5 + 4.5 * 4 ** (l - 7) * ((k + 1) / 4) ** 2,   # l = 7, k = 3 : 4.3 Go
                          quick=[f"mesh.level = {max(l - 3, 1)}", "time.final_time = 0.25", "time.n_steps = 640"],
                          meta={"k": k, "level": l}))

# Table 4, Figs. 12 et 13 : granite-eau realiste, SDIRK(3,4) (direct, iteratif) et ERK(4), k = 1, 3, un fil (CPU)
ERK_STEPS = {1: 400, 3: 1180}
for k in (1, 3):
    for solver in ("direct", "iteratif"):
        st = [f"space.k = {k}"]
        if solver == "iteratif":
            st += ["solver.iterative = true", "solver.tolerance = 1e-2"]
        if k == 3 and solver == "direct":
            st += ["output.snapshots = { 0.36, 0.5 }"]
        CASES.append(Case(f"cmame_tab4_sdirk_{solver}_k{k}", P, "Table 4", f"SDIRK(3,4) {solver} k={k}",
                          "tab4_fig12_13_granite_eau.lua", st, cost=15.0 * k * k, mem=1.0 + 0.5 * k * k,
                          quick=["mesh.nx = 38", "mesh.ny = 26", "time.n_steps = 40"],
                          meta={"scheme": "sdirk", "solver": solver, "k": k}))
    CASES.append(Case(f"cmame_tab4_erk4_k{k}", P, "Table 4", f"ERK(4) k={k}", "tab4_fig12_13_granite_eau.lua",
                      [f"space.k = {k}", "space.order = 'equal'", "space.scaling = '1'", "time.scheme = 'erk'",
                       "time.stages = 4", f"time.n_steps = {ERK_STEPS[k]}", "solver.static_condensation = false"],
                      cost=10.0 * k * k, mem=1.0 + 0.5 * k * k, quick=["mesh.nx = 38", "mesh.ny = 26"],
                      meta={"scheme": "erk", "k": k}))

# ------------------------------------------------------------------------------------------------------------------
# Table 5, Figs. 16 et 17 : bassin sedimentaire, k = 3 sur les quatre maillages ; reference : k = 4 sur les
# triangles, Delta t / 2 (la reference SEM3D de l'article n'est pas disponible)
BASIN = {"triangles": 31.0, "quadrangles": 19.0, "hybride": 18.0, "non_conforme": 48.0}
for mesh, mem in BASIN.items():
    st = [f"space.k = 3"]
    if mesh == "triangles":
        st += ["output.snapshots = { 0.0, 0.16, 0.6, 1.0 }"]
    CASES.append(Case(f"cmame_bassin_{mesh}", P, "Table 5", mesh, "tab5_fig16_17_bassin.lua", st,
                      meshes=[f"bassin_{mesh}.txt"], cost=90.0, threads=8, mem=mem,
                      quick=["space.k = 1", "time.final_time = 1.5 * 20 / 320", "time.n_steps = 20"],
                      meta={"mesh": mesh, "k": 3}))
CASES.append(Case("cmame_bassin_reference", P, "Table 5", "reference k=4", "tab5_fig16_17_bassin.lua",
                  ["space.k = 4", "time.n_steps = 640"], meshes=["bassin_triangles.txt"], cost=360.0, threads=12, mem=55.0,
                  quick=["space.k = 1", "time.final_time = 1.5 * 20 / 320", "time.n_steps = 40"],
                  meta={"mesh": "triangles", "k": 4, "reference": True}))

# ------------------------------------------------------------------------------------------------------------------
REFERENCES = [
    Reference("cmame_granite_eau", (1025.0, 1500.0), (2690.0, 6000.0, 3000.0), (0.0, 200.0), 10.0, "acoustic_pressure", 1.0,
              [("Sf", -500.0, 200.0), ("Ss", -1000.0, -500.0), ("Si", -400.0, 0.0)], 0.5, 0.5 / 4000),
]

if __name__ == "__main__":
    import figures
    if len(sys.argv) > 2 and sys.argv[1] == "--plot-only":
        import plot_mottier
        sys.exit(plot_mottier.make_all(HERE, Path(sys.argv[2]).resolve(), figures))
    sys.exit(main(CASES, HERE, __doc__, REFERENCES, figures))
