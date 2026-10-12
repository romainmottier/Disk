#!/usr/bin/env python3
"""
run_campaign.py — reproduction des figures 1 et 2 de la section 6.1 de [JSC21] avec wave_propagation :

  [JSC21]  E. Burman, O. Duran, A. Ern, M. Steins, Convergence analysis of hybrid high-order methods for
           the wave equation, J. Sci. Comput. 87 (2021) 91, hal-02922720, section 6.1 (Fig. 1 et 2)

Les parametres de chaque figure sont dans les fichiers Lua de ce dossier, ceux qui reproduisent la figure publiee :
temps final T = 0.125 du code de l'article pour la solution (80) (Fig. 1 bas gauche), donnees de Dirichlet de t = 0
pour l'energie (Fig. 1 bas droite). La Fig. 2 publiee est faite des images de resultats acoustiques de [CAMC22]
(Fig. 3 droite, Fig. 4 droite, courbes tau = O(1) de la Fig. 5, identiques pixel pour pixel ou a 5 % pres) : les
fichiers fig2_publiee_*.lua refont ces calculs, les fichiers fig2_sdirk*.lua (cas « variante ») ceux que decrit
l'article (elastodynamique, solutions (78), (79), (80)). Les cas ci-dessous ne fixent que ce qui varie d'une courbe
a l'autre, les ordres attendus et le dessin des courbes (legende et couleurs de l'article). A lancer depuis la copie
de ce dossier dans le build, <build>/apps/wave_propagation/reproducibility/papers_Duran/elastic_jsc2021, ou les
resultats sont ecrits : figure_<n>/ (la figure n reproduite, ses panneaux, la comparaison avec les valeurs relevees
sur la figure publiee, valeurs_article.json, et un dossier par cas), rapport.md, JSC21_recapitulatif.png.
Voir ../LISEZMOI.txt.

Usage:
    python3 run_campaign.py               # campagne complete : ~ 1 h sur 24 coeurs
    python3 run_campaign.py --quick       # maillages <= niveau 4 : ~ 5 min, resultats dans rapide/
    python3 run_campaign.py --only fig2   # les cas dont le nom contient (regex) fig2
    python3 run_campaign.py --list        # liste des cas et leur duree
    python3 run_campaign.py --dry-run     # entree de chaque cas, sans executer
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True                  # pas de __pycache__ a cote de campaign_common.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from campaign_common import Case, col, QUICK_SPACE, QUICK_TIME, main

E2 = ("elastic_u", "elastic_sigma")             # second-order formulation: displacement, stress
E1 = ("elastic_v", "elastic_sigma")             # first-order formulation: velocity, stress
A1 = ("acoustic_p", "acoustic_q")               # first-order acoustic formulation ([CAMC22]): v = p_t, sigma = grad p
NEWMARK = r"Newmark-($\beta$=0.25,$\gamma$=0.5) :"


def labels(cols, k, l2, h1, sym="k", var="h"):
    """Legende d'une courbe d'erreur, comme dans l'article : L2 : k=1 O(h^3), H1 : k=1 O(h^2)"""
    return {cols[0]: rf"L$_2$ : {sym}={k} O(${var}^{{{l2}}}$)", cols[1]: rf"H$^1$ : {sym}={k} O(${var}^{{{h1}}}$)"}


def ref(panel, cols, sym, v):
    """Series de valeurs_article.json comparees aux deux colonnes d'erreur d'un cas"""
    return {"panel": panel, "series": {cols[0]: f"L2_{sym}{v}", cols[1]: f"H1_{sym}{v}"}}


def cases():
    P = "JSC21"
    sp = lambda c, l2, h1, **kw: [col(c[0], "L2", l2, **kw), col(c[1], "H1", h1, **kw)]
    DESCR = "telle que décrite : élastodynamique"
    c = []

    # ---- Fig. 1 : HHO-Newmark, ordre egal
    for k in (1, 2, 3):
        c.append(Case(f"jsc21_fig1_space_k{k}", P, "Fig. 1 (haut, gauche)",
                      r"Newmark, ordre égal, sol. (78), $\Delta t=0.1\cdot 2^{-7}$", "fig1_newmark_space.lua",
                      [f"space.k = {k}"], sp(E2, k + 2, k + 1), ("k", k), quick=QUICK_SPACE,
                      cost=(4.8, 7.2, 17.7)[k - 1], legend=NEWMARK, labels=labels(E2, k, k + 2, k + 1),
                      style={"color": k - 1, "markers": True}, ref=ref("fig_1_haut_gauche", E2, "k", k), points=6))
    c.append(Case("jsc21_fig1_time", P, "Fig. 1 (haut, droite)", r"Newmark, ordre égal, sol. (79), $k=2$, $h=2^{-6}$",
                  "fig1_newmark_time.lua", [], sp(E2, 2, 2), ("k", 2), x="dt", quick=QUICK_TIME, cost=11.2,
                  legend=r"Newmark($\beta$=0.25,$\gamma$ =0.5) :",
                  labels={E2[0]: r"L$_2$ : O($\delta t^2$)", E2[1]: r"H$^1$: O($\delta t^2$)"},
                  style={"color": 0, "markers": True, "pair": True},
                  ref={"panel": "fig_1_haut_droite", "series": {E2[0]: "L2", E2[1]: "H1"}}, points=8))
    for k in (1, 2, 3):
        # k = 1 : regime pre-asymptotique jusqu'au niveau 3 (OCV L2 1.3, 1.6 puis 3.2, comme dans l'article) :
        # ordre L2 retenu sur le dernier segment (mode rapide, niveaux 1..4)
        checks = [col(E2[0], "L2", k + 2, window=(-1, None) if k == 1 else None), col(E2[1], "H1", k + 1)]
        c.append(Case(f"jsc21_fig1_nonpoly_k{k}", P, "Fig. 1 (bas, gauche)",
                      r"Newmark, ordre égal, sol. (80), $T=0.125$, $\Delta t=0.1\cdot 2^{-10}$",
                      "fig1_newmark_space_nonpoly.lua", [f"space.k = {k}"], checks, ("k", k),
                      quick=QUICK_SPACE, cost=(2.1, 10.4, 15.6)[k - 1], legend=NEWMARK,
                      labels=labels(E2, k, k + 2, k + 1), style={"color": k - 1},
                      ref=ref("fig_1_bas_gauche", E2, "k", k), points=6,
                      note="figure publiée : temps final T = 0.125 du code de l'article pour la solution (80) "
                           "(erreurs H1 1.6 à 1.9 fois plus grandes à T = 1)"))
    for beta, gamma in ((0.25, 0.5), (0.275, 0.5), (0.275, 0.55)):
        tag = f"b{beta}_g{gamma}"
        for exact in (False, True):
            c.append(Case(f"jsc21_fig1_energy{'_exact' if exact else ''}_{tag}", P, "Fig. 1 (bas, droite)",
                          r"Newmark, ordre égal, sol. (80), $k=2$, $h=2^{-6}$, $\Delta t=0.1\cdot 2^{-7}$"
                          + ("" if exact else ", données de Dirichlet de t = 0"), "fig1_newmark_energy.lua",
                          [f"time.newmark_beta = {beta}", f"time.newmark_gamma = {gamma}"]
                          + (["problem.dirichlet = 'exact'"] if exact else []), [], ("k", 2), quick=["mesh.level = 4"],
                          energy=True, energy_check=1e-10 if (beta, gamma) == (0.25, 0.5) and not exact else None,
                          cost=6.0, legend=r"Newmark-($\beta$,$\gamma$) :",
                          labels={"energy": rf"($\beta$,$\gamma$) = ({beta},{gamma})"},
                          style={"color": [(0.25, 0.5), (0.275, 0.5), (0.275, 0.55)].index((beta, gamma)),
                                 "solid": (beta, gamma) == (0.25, 0.5), "loc": "center right"},
                          variant="données de Dirichlet exactes de (80)" if exact else "",
                          ref=None if exact else {"panel": "fig_1_bas_droite", "series": {"energy": tag}},
                          note="" if exact else
                          "figure publiée : données de Dirichlet de t = 0 (nulles), conservation à 1e-15 comme dans "
                          "l'article ; mêmes comportements que les courbes publiées, amplitudes 2 à 6 fois plus faibles "
                          "pour (β, γ) ≠ (0.25, 0.5) (écart non expliqué) ; avec les données exactes de (80), l'énergie "
                          "(69) dérive en O(Δt²) (variante)"))

    # ---- Fig. 2, telle que publiee : resultats acoustiques de [CAMC22]
    for k in range(3):
        c.append(Case(f"jsc21_fig2_publiee_space_k{k}", P, "Fig. 2 (haut, gauche)",
                      r"telle que publiée = [CAMC22] Fig. 3 (droite) : SDIRK(2,3), acoustique, ordre mixte, "
                      r"$\tau=O(h^{-1})$", "fig2_publiee_sdirk23_space.lua", [f"space.k = {k}"],
                      sp(A1, k + 2, k + 1), ("k", k), quick=QUICK_SPACE, cost=(0.4, 1.0, 2.2)[k],
                      legend="SDIRK(2,3)–HHO:", labels=labels(A1, k, k + 2, k + 1), style={"color": k, "markers": True},
                      ref=ref("fig_2_haut_gauche", A1, "k", k), points=6,
                      note="l'image publiée est celle de [CAMC22] Fig. 3 (droite), identique pixel pour pixel : "
                           "résultats acoustiques (solution (60) de [CAMC22], ordre mixte, k = 0, 1, 2)"))
        c.append(Case(f"jsc21_fig2_space_k{k + 1}", P, "Fig. 2 (haut, gauche)",
                      r"SDIRK(2,3), ordre égal, sol. (78), $\Delta t=0.1\cdot 2^{-7}$", "fig2_sdirk23_space.lua",
                      [f"space.k = {k + 1}"], sp(E1, k + 3, k + 2), ("k", k + 1), quick=QUICK_SPACE,
                      cost=(8.8, 18.8, 25.4)[k], legend="SDIRK(2,3)–HHO:", labels=labels(E1, k + 1, k + 3, k + 2),
                      style={"color": k, "markers": True}, points=6, variant=DESCR + ", sol. (78), k = 1, 2, 3"))
    for s in (1, 2, 3):
        c.append(Case(f"jsc21_fig2_publiee_time_s{s}", P, "Fig. 2 (haut, droite)",
                      r"telle que publiée = [CAMC22] Fig. 4 (droite) : SDIRK($s$,$s$+1), acoustique, $k=3$, $h=2^{-3}$",
                      "fig2_publiee_sdirk_time.lua", [f"time.stages = {s}"],
                      [col(A1[0], "L2", s + 1), col(A1[1], "H1", s + 1, tol=0.4 if s == 3 else 0.3)], ("s", s),
                      x="dt", cost=(0.1, 0.3, 0.2)[s - 1], legend="SDIRK–HHO :",
                      labels={A1[0]: rf"L$_2$ : O($\Delta t^{s + 1}$)", A1[1]: rf"H$^1$: O($\Delta t^{s + 1}$)"},
                      style={"color": s - 1, "markers": True}, ref=ref("fig_2_haut_droite", A1, "s", s), points=8,
                      note="courbes publiées identiques pixel pour pixel à celles de [CAMC22] Fig. 4 (droite) "
                           "(acoustique, solution (61) de [CAMC22], k = 3, h = 2⁻³), seule la légende diffère"))
        c.append(Case(f"jsc21_fig2_time_s{s}", P, "Fig. 2 (haut, droite)",
                      r"SDIRK($s$,$s$+1), ordre égal, sol. (79), $k=2$, $h=2^{-6}$", "fig2_sdirk_time.lua",
                      [f"time.stages = {s}"], sp(E1, s + 1, s + 1), ("s", s), x="dt", quick=QUICK_TIME,
                      cost=(14.7, 22.7, 25.8)[s - 1], legend="SDIRK–HHO :",
                      labels={E1[0]: rf"L$_2$ : O($\Delta t^{s + 1}$)", E1[1]: rf"H$^1$: O($\Delta t^{s + 1}$)"},
                      style={"color": s - 1, "markers": True}, points=8, variant=DESCR + ", sol. (79), k = 2, h = 2⁻⁶"))
    for order, side in (("equal", "gauche"), ("mixed", "droite")):
        stab = "hho" if order == "equal" else "least_squares"
        word = "égal" if order == "equal" else "mixte"
        for k in range(3):
            c.append(Case(f"jsc21_fig2_publiee_energy_{order}_k{k}", P, f"Fig. 2 (bas, {side})",
                          r"telle que publiée = [CAMC22] Fig. 5, $\tau=O(1)$ : SDIRK(3,4), acoustique, ordre " + word
                          + r", $h=2^{-6}$", "fig2_publiee_sdirk34_energy.lua",
                          [f"space.k = {k}", f"space.order = '{order}'", f"space.stabilization = '{stab}'"], [],
                          ("k", k), quick=["mesh.level = 4"], energy=True, cost=(0.7, 1.6, 3.9)[k],
                          legend="SDIRK(3,4)–HHO :", labels={"energy": f"k = {k + 1}"}, style={"color": k},
                          ref={"panel": f"fig_2_bas_{side}", "series": {"energy": f"k{k + 1}"}},
                          note="courbes publiées = courbes τ = O(1) de [CAMC22] Fig. 5 (acoustique, k = 0, 1, 2, "
                               "légendées k = 1, 2, 3)"))
            c.append(Case(f"jsc21_fig2_energy_{order}_k{k + 1}", P, f"Fig. 2 (bas, {side})",
                          r"SDIRK(3,4), ordre " + word + r", sol. (80), $h=2^{-5}$, $\Delta t=0.1\cdot 2^{-7}$",
                          "fig2_sdirk34_energy.lua",
                          [f"space.k = {k + 1}", f"space.order = '{order}'", f"space.stabilization = '{stab}'"], [],
                          ("k", k + 1), quick=["mesh.level = 4"], energy=True, cost=(1.3, 6.4, 9.4)[k],
                          legend="SDIRK(3,4)–HHO :", labels={"energy": f"k = {k + 1}"}, style={"color": k},
                          variant=DESCR + ", sol. (80), k = 1, 2, 3"))
    return c


if __name__ == "__main__":
    sys.exit(main(cases(), Path(__file__).resolve().parent, __doc__))
