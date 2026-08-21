#!/usr/bin/env python3
"""
Etude de la mise a l'echelle en h des valeurs propres de LHS
(LHS_spectrum_acou.hpp, fichier "lhs_scaling_summary.txt").

Pour chacune des quatre quantites :
  rho_full  : rayon spectral de LHS entier
  rho_cell  : rayon spectral restreint aux vap dont le vecteur propre est
              majoritairement dans le bloc cellule-cellule (frac_cell > 0.9)
  rho_face  : idem pour le bloc face (frac_cell < 0.1)
  rho_schur : rayon spectral du complement de Schur n_c x n_c
              S_SCHUR = Kcc - Kcf*Sff_inv*Kfc (faces eliminees par
              condensation statique) -- l'operateur cellule-cellule
              effectif reellement utilise par les schemas explicites

on fait une regression log-log :  log(rho) = p * log(h) + c
et on affiche la pente p. p ~ -1 -> rho ~ 1/h ; p ~ -2 -> rho ~ 1/h^2.

Usage:
    python3 plot_lhs_scaling.py [dir]

    dir : repertoire contenant "lhs_scaling_summary.txt"
          (par defaut : repertoire courant)

Genere "lhs_scaling.png" (log-log, avec droites de reference h^-1 et h^-2)
et affiche les pentes ajustees dans la console.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def load_summary(path):
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # colonnes : level nx ny h_max n_dof n_c n_f rho_full rho_cell rho_face rho_schur
    h         = data[:, 3]
    rho_full  = data[:, 7]
    rho_cell  = data[:, 8]
    rho_face  = data[:, 9]
    rho_schur = data[:, 10]
    order = np.argsort(h)[::-1]  # h decroissant -> maillage de plus en plus fin
    return h[order], rho_full[order], rho_cell[order], rho_face[order], rho_schur[order]


def fit_loglog(h, rho):
    """Ajuste log(rho) = p*log(h) + c par moindres carres. Ignore rho<=0."""
    mask = rho > 0
    if mask.sum() < 2:
        return None, None
    p, c = np.polyfit(np.log(h[mask]), np.log(rho[mask]), 1)
    return p, c


def plot_series(ax, h, rho, label, color):
    ax.loglog(h, rho, "o-", color=color, label=label, zorder=3)
    p, c = fit_loglog(h, rho)
    if p is not None:
        fit_line = np.exp(c) * h ** p
        ax.loglog(h, fit_line, "--", color=color, linewidth=1,
                  label=rf"{label} fit: $h^{{{p:.2f}}}$")
    return p


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    summary_path = os.path.join(directory, "lhs_scaling_summary.txt")

    if not os.path.isfile(summary_path):
        print(f"Fichier introuvable : {summary_path}")
        sys.exit(1)

    h, rho_full, rho_cell, rho_face, rho_schur = load_summary(summary_path)

    print(f"{len(h)} niveaux de maillage charges depuis '{summary_path}':")
    for hi, rf, rc, rface, rs in zip(h, rho_full, rho_cell, rho_face, rho_schur):
        print(f"  h={hi:.6g}   rho_full={rf:.6g}   rho_cell={rc:.6g}   rho_face={rface:.6g}   rho_schur={rs:.6g}")

    fig, ax = plt.subplots(figsize=(7, 6))

    series = [
        ("rho_full (spectre complet)", rho_full, "tab:blue"),
        ("rho_cell (bloc cellule-cellule)", rho_cell, "tab:red"),
        ("rho_face (bloc face)", rho_face, "tab:green"),
        ("rho_schur (complement de Schur)", rho_schur, "tab:purple"),
    ]

    print("\nPentes ajustees (log(rho) = p*log(h) + c) :")
    print("  p > 0 : rho -> 0 quand h -> 0 (rho ~ h^p)")
    print("  p < 0 : rho -> infini quand h -> 0 (rho ~ h^p = 1/h^|p|)")
    slopes = {}
    for label, rho, color in series:
        p = plot_series(ax, h, rho, label, color)
        slopes[label] = p
        if p is not None:
            candidates = {"h^0 (constante)": 0.0, "h^1": 1.0, "h^2": 2.0, "1/h": -1.0, "1/h^2": -2.0}
            best = min(candidates, key=lambda k: abs(p - candidates[k]))
            print(f"  {label:35s} p = {p:+.3f}   proche de {best}")
        else:
            print(f"  {label:35s} pas assez de points positifs pour un fit")

    # Droites de reference h^-1 et h^-2, calees sur rho_full au plus petit h
    h_ref = np.array([h.min(), h.max()])
    if rho_full[np.argmin(h)] > 0:
        anchor_h = h[np.argmin(h)]
        anchor_rho = rho_full[np.argmin(h)]
        ax.loglog(h_ref, anchor_rho * (h_ref / anchor_h) ** (-1), ":", color="grey", linewidth=1, label=r"reference $h^{-1}$")
        ax.loglog(h_ref, anchor_rho * (h_ref / anchor_h) ** (-2), ":", color="black", linewidth=1, label=r"reference $h^{-2}$")

    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"$\rho$")
    ax.set_title("Mise a l'echelle en h des rayons spectraux de LHS")
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(fontsize=8, loc="best")

    out_name = os.path.join(directory, "lhs_scaling.png")
    fig.savefig(out_name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure -> {out_name}")


if __name__ == "__main__":
    main()
