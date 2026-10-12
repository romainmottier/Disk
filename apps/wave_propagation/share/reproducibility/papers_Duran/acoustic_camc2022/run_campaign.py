#!/usr/bin/env python3
"""
run_campaign.py — reproduction des figures 1 a 8 de la section 5.1 de [CAMC22] avec wave_propagation :

  [CAMC22] E. Burman, O. Duran, A. Ern, Hybrid high-order methods for the acoustic wave equation in the
           time domain, Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, section 5.1

Les parametres de chaque figure sont dans les fichiers Lua de ce dossier, ceux qui reproduisent la figure publiee ;
quand ils different de la legende de l'article (pas de temps de la Fig. 2, cos au lieu de sin dans la solution de la
Fig. 1 droite, maillage de la Fig. 5), le fichier Lua l'explique et des cas « variante » refont le calcul avec les
parametres de la legende. Les cas ci-dessous ne fixent que ce qui varie d'une courbe a l'autre (k, nombre d'etages s,
echelle de la stabilisation, (beta, gamma), maillages), les ordres attendus et le dessin des courbes (legende et
couleurs de l'article). A lancer depuis la copie de ce dossier dans le build,
<build>/apps/wave_propagation/reproducibility/papers_Duran/acoustic_camc2022, ou les resultats sont ecrits :
figure_<n>/ (la figure n reproduite, ses panneaux, la comparaison avec les valeurs relevees sur la figure publiee,
valeurs_article.json, et un dossier par cas), rapport.md, CAMC22_recapitulatif.png. Voir ../LISEZMOI.txt.

Usage:
    python3 run_campaign.py               # campagne complete : ~ 1 h sur 24 coeurs
    python3 run_campaign.py --quick       # maillages <= niveau 4 : ~ 5 min, resultats dans rapide/
    python3 run_campaign.py --only fig3   # les cas dont le nom contient (regex) fig3
    python3 run_campaign.py --list        # liste des cas et leur duree
    python3 run_campaign.py --dry-run     # entree de chaque cas, sans executer
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True                  # pas de __pycache__ a cote de campaign_common.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from campaign_common import Case, col, QUICK_SPACE, QUICK_TIME, main

C2 = ("acoustic_phi", "acoustic_grad_phi")      # second-order formulation: p, reconstructed gradient
C1 = ("acoustic_p", "acoustic_q")               # first-order formulation: v = p_t, sigma = grad p
NEWMARK = r"Newmark-($\beta$=0.25,$\gamma$=0.5) :"
POLY = [f"unit_square_polymesh_nel_{20 * 2 ** l}.txt" for l in range(9)]    # 20 x 2^l cellules, maillages de l'article


def labels(cols, k, l2, h1, sym="k", var="h"):
    """Legende d'une courbe d'erreur, comme dans l'article : L2 : k=0 O(h^2), H1 : k=0 O(h^1)"""
    return {cols[0]: rf"L$_2$ : {sym}={k} O(${var}^{{{l2}}}$)", cols[1]: rf"H$^1$ : {sym}={k} O(${var}^{{{h1}}}$)"}


def ref(panel, cols, sym, v):
    """Series de valeurs_article.json comparees aux deux colonnes d'erreur d'un cas"""
    return {"panel": panel, "series": {cols[0]: f"L2_{sym}{v}", cols[1]: f"H1_{sym}{v}"}}


def cases():
    P = "CAMC22"
    sp = lambda c, l2, h1, **kw: [col(c[0], "L2", l2, **kw), col(c[1], "H1", h1, **kw)]
    # tau = O(1) : l'article annonce un ordre L2 sous-optimal, k+1 au lieu de k+2 (perte de la supraconvergence) ;
    # pour k = 0 les pentes H1 des figures de l'article restent nettement < 1 (Fig. 6 droite, Fig. 8 bas gauche) :
    # tolerance 0.4
    sub = lambda c, k: [col(c[0], "L2", k + 1, mode="equal"), col(c[1], "H1", k + 1, tol=0.4 if k == 0 else 0.3)]
    c = []

    # ---- Fig. 1 : HHO-Newmark, ordre egal
    for k in range(3):
        c.append(Case(f"camc22_fig1_space_k{k}", P, "Fig. 1 (gauche)",
                      r"Newmark, ordre égal, sol. (60), $\Delta t=0.1\cdot 2^{-7}$", "fig1_newmark_space.lua",
                      [f"space.k = {k}"], sp(C2, k + 2, k + 1), ("k", k), quick=QUICK_SPACE, cost=(0.7, 2.1, 3.2)[k],
                      legend=NEWMARK, labels=labels(C2, k, k + 2, k + 1), style={"color": k, "markers": True},
                      ref=ref("fig_1_gauche", C2, "k", k), points=6))
    time_labels = {C2[0]: r"L$_2$ : O($\delta t^2$)", C2[1]: r"H$^1$: O($\delta t^2$)"}
    c.append(Case("camc22_fig1_time", P, "Fig. 1 (droite)",
                  r"Newmark, ordre égal, $p=\cos(\sqrt{2}\pi t)\,x(1-x)\,y(1-y)$, $k=2$, $h=2^{-6}$",
                  "fig1_newmark_time.lua", [], sp(C2, 2, 2), ("k", 2), x="dt", quick=QUICK_TIME, cost=3.6,
                  legend=r"Newmark($\beta$=0.25,$\gamma$ =0.5) :", labels=time_labels,
                  style={"color": 0, "markers": True, "pair": True},
                  ref={"panel": "fig_1_droite", "series": {C2[0]: "L2", C2[1]: "H1"}}, points=8,
                  note="figure publiée : solution en cos(√2πt) ; la solution (61) imprimée (sin) donne des erreurs "
                       "|tan(√2π)| = 3.6 fois plus petites (variante)"))
    c.append(Case("camc22_fig1_time_sin", P, "Fig. 1 (droite)",
                  r"Newmark, ordre égal, sol. (61) (sin), $k=2$, $h=2^{-6}$", "fig1_newmark_time.lua",
                  ["problem.solution = 'acoustic_quadratic_in_space'"], sp(C2, 2, 2), ("k", 2), x="dt", quick=QUICK_TIME,
                  cost=3.6, variant="(61) telle qu'imprimée, sin(√2πt)", legend=r"Newmark($\beta$=0.25,$\gamma$ =0.5) :",
                  labels=time_labels, style={"color": 0, "markers": True, "pair": True}, points=8))

    # ---- Fig. 2 : HHO-Newmark, ordre egal, sol. (62)
    for k in range(3):
        # pour k = 2 l'erreur L2 des maillages les plus fins est saturee par l'erreur en temps (comme dans l'article)
        w = (-4, -2) if k == 2 else None
        for dt7 in (False, True):
            c.append(Case(f"camc22_fig2_space{'_dt7' if dt7 else ''}_k{k}", P, "Fig. 2 (gauche)",
                          r"Newmark, ordre égal, sol. (62), $\Delta t=0.1\cdot 2^{-" + ("7" if dt7 else "6") + "}$",
                          "fig2_newmark_space.lua", [f"space.k = {k}"] + (["time.n_steps = 1280"] if dt7 else []),
                          [col(C2[0], "L2", k + 2, window=w), col(C2[1], "H1", k + 1)], ("k", k), quick=QUICK_SPACE,
                          cost=((0.5, 2.0, 3.7) if dt7 else (0.3, 0.9, 1.9))[k], legend=NEWMARK,
                          labels=labels(C2, k, k + 2, k + 1), style={"color": k}, points=6,
                          variant="Δt de la légende, 0.1·2⁻⁷" if dt7 else "",
                          ref=None if dt7 else ref("fig_2_gauche", C2, "k", k),
                          note="" if dt7 else "figure publiée : Δt = 0.1·2⁻⁶ (légende : 0.1·2⁻⁷, variante)"))
    for beta, gamma in ((0.25, 0.5), (0.275, 0.5), (0.275, 0.55)):
        tag = f"b{beta}_g{gamma}"
        for dt7 in (False, True):
            c.append(Case(f"camc22_fig2_energy{'_dt7' if dt7 else ''}_{tag}", P, "Fig. 2 (droite)",
                          r"Newmark, ordre égal, sol. (62), $k=2$, $h=2^{-6}$, $\Delta t=0.1\cdot 2^{-"
                          + ("7" if dt7 else "5") + "}$", "fig2_newmark_energy.lua",
                          [f"time.newmark_beta = {beta}", f"time.newmark_gamma = {gamma}"]
                          + (["time.n_steps = 1280"] if dt7 else []), [], ("k", 2), quick=["mesh.level = 4"],
                          energy=True, energy_check=1e-10 if (beta, gamma) == (0.25, 0.5) else None,
                          cost=2.6 if dt7 else 0.7, legend=r"Newmark-($\beta$,$\gamma$) :",
                          labels={"energy": rf"($\beta$,$\gamma$) = ({beta},{gamma})"},
                          style={"color": [(0.25, 0.5), (0.275, 0.5), (0.275, 0.55)].index((beta, gamma)),
                                 "solid": (beta, gamma) == (0.25, 0.5), "loc": "center left"},
                          variant="Δt de la légende, 0.1·2⁻⁷" if dt7 else "",
                          ref=None if dt7 else {"panel": "fig_2_droite", "series": {"energy": tag}},
                          note="" if dt7 else "figure publiée : Δt = 0.1·2⁻⁵ (légende : 0.1·2⁻⁷, variante)"))

    # ---- Fig. 3 : HHO-SDIRK(2,3), ordre mixte
    for k in range(3):
        c.append(Case(f"camc22_fig3_tau1_k{k}", P, "Fig. 3 (gauche)",
                      r"SDIRK(2,3), ordre mixte, $\tau=O(1)$, sol. (60)", "fig3_sdirk23_space.lua",
                      [f"space.k = {k}", "space.scaling = '1'"], sub(C1, k), ("k", k), quick=QUICK_SPACE,
                      cost=(0.5, 1.9, 4.1)[k], legend="SDIRK(2,3)–HHO:", labels=labels(C1, k, k + 1, k + 1),
                      style={"color": k, "markers": True}, ref=ref("fig_3_gauche", C1, "k", k), points=6))
        c.append(Case(f"camc22_fig3_tauh_k{k}", P, "Fig. 3 (droite)",
                      r"SDIRK(2,3), ordre mixte, $\tau=O(h^{-1})$, sol. (60)", "fig3_sdirk23_space.lua",
                      [f"space.k = {k}", "space.scaling = '1/h'"], sp(C1, k + 2, k + 1), ("k", k), quick=QUICK_SPACE,
                      cost=(0.7, 2.4, 4.0)[k], legend="SDIRK(2,3)–HHO:", labels=labels(C1, k, k + 2, k + 1),
                      style={"color": k, "markers": True}, ref=ref("fig_3_droite", C1, "k", k), points=6))

    # ---- Fig. 4 : HHO-SDIRK, ordre egal, tau = O(1/h)
    for k in range(3):
        c.append(Case(f"camc22_fig4_space_k{k}", P, "Fig. 4 (gauche)",
                      r"SDIRK(3,4), ordre égal, $\tau=O(h^{-1})$, sol. (62), $\Delta t=0.1\cdot 2^{-8}$",
                      "fig4_sdirk34_space.lua", [f"space.k = {k}"], sp(C1, k + 2, k + 1), ("k", k),
                      quick=QUICK_SPACE, cost=(0.9, 2.8, 4.3)[k], legend="SDIRK(3,4)–HHO:",
                      labels=labels(C1, k, k + 2, k + 1), style={"color": k}, ref=ref("fig_4_gauche", C1, "k", k),
                      points=5))
    for s in (1, 2, 3):
        # s = 3 : apres la bosse vers dt = 0.025 (visible aussi sur la figure de l'article), l'erreur sur le flux
        # rejoint l'ordre 4 lentement (OCV 3.5, 3.6, 3.7 sur les derniers pas) : tolerance 0.4, et les 8 pas de temps
        # sont necessaires, y compris en mode rapide (cas peu couteux, h = 2^-3)
        checks = [col(C1[0], "L2", s + 1), col(C1[1], "H1", s + 1, tol=0.4 if s == 3 else 0.3)]
        c.append(Case(f"camc22_fig4_time_s{s}", P, "Fig. 4 (droite)",
                      r"SDIRK($s$,$s$+1), ordre égal, $\tau=O(h^{-1})$, sol. (61), $k=3$, $h=2^{-3}$",
                      "fig4_sdirk_time.lua", [f"time.stages = {s}"], checks, ("s", s), x="dt", cost=(0.2, 0.4, 0.5)[s - 1],
                      legend="SDIRK–HHO :", labels=labels(C1, s, s + 1, s + 1, sym="s", var=r"\Delta t"),
                      style={"color": s - 1, "markers": True}, ref=ref("fig_4_droite", C1, "s", s), points=8))

    # ---- Fig. 5 : HHO-SDIRK(3,4), perte d'energie relative
    for order, side in (("equal", "gauche"), ("mixed", "droite")):
        stab = "hho" if order == "equal" else "least_squares"
        for tau in ("1", "h"):
            for k in range(3):
                c.append(Case(f"camc22_fig5_{order}_tau{tau}_k{k}", P, f"Fig. 5 ({side})",
                              r"SDIRK(3,4), ordre " + ("égal" if order == "equal" else "mixte")
                              + r", sol. (62), $h=2^{-6}$, $\Delta t=0.1\cdot 2^{-7}$", "fig5_sdirk34_energy.lua",
                              [f"space.k = {k}", f"space.order = '{order}'", f"space.stabilization = '{stab}'",
                               "space.scaling = '1'" if tau == "1" else "space.scaling = '1/h'"],
                              [], ("k", k), quick=["mesh.level = 4"], energy=True, cost=(1.6, 4.6, 6.6)[k],
                              legend=f"SDIRK(3,4)–HHO {'equal' if order == 'equal' else 'mixed'}-order:",
                              labels={"energy": f"k = {k} – O(1)" if tau == "1" else rf"k = {k} – O($h_T$$^{{-1}}$)"},
                              style={"color": k, "thick": tau == "h"},
                              ref={"panel": f"fig_5_{side}", "series": {"energy": f"k{k}_tau{tau}"}},
                              note="figure publiée : h = 2⁻⁶ (légende : 2⁻⁵) ; " + (
                                  "seule la courbe k = 1, τ = O(1) de l'ordre égal ne concorde pas (4 fois au-dessous "
                                  "de la nôtre, qui coïncide avec la courbe de l'ordre mixte)" if order == "equal"
                                  else "pertes presque identiques en ordre égal et mixte, comme dans l'article")))

    # ---- Fig. 6 : HHO-ERK(4), ordre mixte
    for k in range(3):
        c.append(Case(f"camc22_fig6_tauh_k{k}", P, "Fig. 6 (gauche)",
                      r"ERK(4), ordre mixte, $\tau=O(h^{-1})$, sol. (62), $\Delta t=h/320$", "fig6_erk4_space.lua",
                      [f"space.k = {k}", "space.scaling = '1/h'"], sp(C1, k + 2, k + 1), ("k", k), quick=QUICK_SPACE,
                      cost=(3.7, 5.5, 11.7)[k], legend="ERK(4,4)–HHO:", labels=labels(C1, k, k + 2, k + 1),
                      style={"color": k}, ref=ref("fig_6_gauche", C1, "k", k), points=5,
                      note="écart non expliqué : erreurs H1 1.2 à 1.3 fois au-dessus de l'article (L2, k = 0 : "
                           "nettement au-dessous) ; même discrétisation en espace que les Fig. 3 (droite) et 4 "
                           "(gauche), qui concordent, et ERK(4) = SDIRK(3,4) à 7 chiffres sur ces maillages"))
        c.append(Case(f"camc22_fig6_tau1_k{k}", P, "Fig. 6 (droite)",
                      r"ERK(4), ordre mixte, $\tau=O(1)$, sol. (62), $\Delta t=h/320$", "fig6_erk4_space.lua",
                      [f"space.k = {k}", "space.scaling = '1'"], sub(C1, k), ("k", k), quick=QUICK_SPACE,
                      cost=(1.7, 5.6, 10.6)[k], legend="ERK(4,4)–HHO:", labels=labels(C1, k, k + 1, k + 1),
                      style={"color": k}, ref=ref("fig_6_droite", C1, "k", k), points=5,
                      note="dans l'article aussi, pour k = 0 et tau = O(1), les pentes restent inferieures a 1 sur ces maillages"))

    # ---- Fig. 8 : maillages polygonaux PolyMesher de l'article (Fig. 7 : les quatre premiers)
    first_order = ["formulation = 'first_order'"]
    erk = first_order + ["space.order = 'mixed'", "space.stabilization = 'least_squares'", "time.scheme = 'erk'",
                         "time.stages = 4", "time.dt_h = 1/320", "study.last_level = 6"]
    panels = [
        ("newmark_equal", "Fig. 8 (haut, gauche)", r"Newmark, ordre égal, sol. (60), $\Delta t=0.1\cdot 2^{-5}$", [],
         C2, (2, 1), NEWMARK, (0.7, 1.4, 3.8)),
        ("newmark_mixed", "Fig. 8 (haut, droite)", r"Newmark, ordre mixte, sol. (60), $\Delta t=0.1\cdot 2^{-5}$",
         ["space.order = 'mixed'", "space.stabilization = 'least_squares'"], C2, (2, 1), NEWMARK, (0.7, 2.1, 3.4)),
        ("sdirk34_equal_tau1", "Fig. 8 (milieu, gauche)",
         r"SDIRK(3,4), ordre égal, $\tau=O(1)$, sol. (60), $\Delta t=0.1\cdot 2^{-5}$",
         first_order + ["time.scheme = 'sdirk'", "time.stages = 3", "space.scaling = '1'"], C1, (1, 1), "SDIRK(3,4):",
         (1.9, 5.7, 8.5)),
        ("sdirk34_equal_tauh", "Fig. 8 (milieu, droite)",
         r"SDIRK(3,4), ordre égal, $\tau=O(h^{-1})$, sol. (60), $\Delta t=0.1\cdot 2^{-5}$",
         first_order + ["time.scheme = 'sdirk'", "time.stages = 3"], C1, (2, 1), "SDIRK(3,4):", (1.8, 5.0, 8.1)),
        ("erk4_mixed_tau1", "Fig. 8 (bas, gauche)", r"ERK(4), ordre mixte, $\tau=O(1)$, sol. (60), $\Delta t=h/320$",
         erk + ["space.scaling = '1'"], C1, (1, 1), "ERK(4):", (3.3, 6.3, 14.7)),
        ("erk4_mixed_tauh", "Fig. 8 (bas, droite)", r"ERK(4), ordre mixte, $\tau=O(h^{-1})$, sol. (60), $\Delta t=h/320$",
         erk + ["space.scaling = '1/h'"], C1, (2, 1), "ERK(4):", (2.7, 6.5, 15.6)),
    ]
    side = {"haut, gauche": "haut_gauche", "haut, droite": "haut_droite", "milieu, gauche": "milieu_gauche",
            "milieu, droite": "milieu_droite", "bas, gauche": "bas_gauche", "bas, droite": "bas_droite"}
    for tag, panel, box, stmts, cc, (dl2, dh1), legend, cost in panels:
        n = 6 if tag.startswith("erk") else 9
        key = "fig_8_" + side[panel[panel.index("(") + 1:-1]]
        for k in range(3):                          # k = 0..2 comme dans l'article
            checks = sub(cc, k) if dl2 == 1 else sp(cc, k + dl2, k + dh1)
            # ERK(4), tau = O(1/h), k = 2 : limite de stabilite en h^2 (texte de l'article), depassee par dt = h/320 sur
            # les maillages 5 et 6 (valeurs non finies) ; dt = h/640, sans effet sur les erreurs (erreur en temps
            # negligeable : memes erreurs a 6 chiffres sur les maillages 1 a 4)
            dt640 = tag == "erk4_mixed_tauh" and k == 2
            # legende de l'article : ordres de tau = O(1/h) dans les six panneaux
            c.append(Case(f"camc22_fig8_{tag}_k{k}", P, panel, box + "\nmaillages polygonaux (PolyMesher)",
                          "fig8_polygonal.lua", stmts + [f"space.k = {k}"] + (["time.dt_h = 1/640"] if dt640 else []),
                          checks, ("k", k), x="h_min", quick=QUICK_SPACE, polymesher=POLY[:n],
                          cost=2 * cost[k] if dt640 else cost[k], legend=legend,
                          labels=labels(cc, k, k + 2, k + 1), style={"color": k}, ref=ref(key, cc, "k", k), points=n,
                          note=("pour k = 0 et tau = O(1), la pente H1 reste inferieure a 1 sur ces maillages, comme sur "
                                "la figure de l'article ; " if tag == "erk4_mixed_tau1" else "")
                               + ("k = 2 : dt = h/640, ERK(4) etant instable avec dt = h/320 sur les maillages 5 et 6 "
                                  "(tau = O(1/h)) ; erreurs inchangees sur les maillages 1 a 4 ; " if dt640 else "")
                               + ("k = 0 en ordre egal : notre erreur L2 decroit plus vite que celle de l'article "
                                  "(rapport nous/article de 0.9 a 0.3 du maillage le plus grossier au plus fin), ecart "
                                  "non explique ; les autres courbes concordent a 10 % ; "
                                  if k == 0 and tag in ("newmark_equal", "sdirk34_equal_tauh") else "")
                               + "maillages PolyMesher de l'article (meshes/), h = plus petit diamètre de cellule"))
    return c


if __name__ == "__main__":
    sys.exit(main(cases(), Path(__file__).resolve().parent, __doc__))
