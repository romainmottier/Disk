#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Log-log convergence plot for ERK4_LTS_Lshape_MMS_conv_test.hpp: pressure L2
error against the closed-form manufactured solution (corner-singular
harmonic function times a smooth time profile), as a function of h_max, for
the actual explicit LTS-RK4 + HHO time integrator on the graded L-shaped
reentrant-corner mesh. Unlike plot_lshape_convergence.py (which compares
against a numerically-computed reference mesh), every point here is an
exact error, no reference-mesh dependency.

Usage:
    python3 plot_lshape_mms_convergence.py [file]

    file : the "lshape_mms_convergence_k_<k>.txt" log written by the C++
           prototype (default: first match of
           lshape_mms_convergence_k_*.txt in the current directory, else
           in ../results/ relative to it -- matches the
           build/apps/wave_propagation/lshape/{results,plots}/ layout).
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

SURFACE, GRID, AXIS, INK, MUTED = "#fcfcfb", "#e1e0d9", "#c3c2b7", "#0b0b0b", "#52514e"


def load(path):
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # columns: N  h_max  h_min  L2_error_pressure
    N, h_max, h_min, err = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    order = np.argsort(h_max)
    return N[order], h_max[order], h_min[order], err[order]


def find_default_input():
    matches = sorted(glob.glob("lshape_mms_convergence_k_*.txt"))
    if matches:
        return matches[0]
    matches = sorted(glob.glob("results/lshape_mms_convergence_k_*.txt"))
    if matches:
        return matches[0]
    return None


def output_path_for(path):
    abspath = os.path.abspath(path)
    if os.sep + "results" + os.sep in abspath:
        plots_path = abspath.replace(os.sep + "results" + os.sep, os.sep + "plots" + os.sep)
        os.makedirs(os.path.dirname(plots_path), exist_ok=True)
        return os.path.splitext(plots_path)[0] + ".png"
    return os.path.splitext(path)[0] + ".png"


def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = find_default_input()
        if path is None:
            print("Aucun fichier lshape_mms_convergence_k_*.txt trouve.")
            sys.exit(1)

    N, h_max, h_min, err = load(path)
    print(f"Charge depuis '{path}':")
    for ni, hi, hmi, ei in zip(N, h_max, h_min, err):
        print(f"  N={int(ni):4d}   h_max={hi:.6g}   h_min={hmi:.6g}   L2_error(pressure)={ei:.6g}")

    p_fit, c_fit = np.polyfit(np.log(h_max), np.log(err), 1)
    print(f"\nPente ajustee sur tous les points (log(err) = p*log(h_max) + c) : p = {p_fit:+.3f}")
    if len(h_max) >= 3:
        p_ends = (np.log(err[-1]) - np.log(err[0])) / (np.log(h_max[-1]) - np.log(h_max[0]))
        print(f"Pente agregee premier<->dernier point : p = {p_ends:+.3f}")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.loglog(h_max, err, "o-", color="#2a78d6", lw=2, ms=6, zorder=3,
              label="L2 error (pressure) vs solution exacte")

    h_ref = np.array([h_max.min(), h_max.max()])
    anchor_h, anchor_e = h_max[0], err[0]
    ax.loglog(h_ref, anchor_e * (h_ref / anchor_h) ** p_fit, "--", color=MUTED, lw=1,
              label=rf"fit: $h_{{max}}^{{{p_fit:.2f}}}$")

    ax.set_xlabel(r"$h_{\max}$ (taille de maille du fond)", color=INK, fontsize=11)
    ax.set_ylabel(r"$\|p_h - p_{\mathrm{exact}}\|_{L^2}$", color=INK, fontsize=11)
    ax.set_title("L-shape coin rentrant: HHO+ERK4-LTS, solution manufacturee\n(maillage gradue, erreur exacte, pas de reference numerique)",
                 color=INK, fontsize=10)
    ax.grid(True, which="both", color=GRID, lw=0.4, alpha=0.9)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend = ax.legend(fontsize=9, frameon=True, facecolor=SURFACE, edgecolor=AXIS, labelcolor=INK)
    legend.get_frame().set_linewidth(0.8)

    fig.tight_layout()
    out_path = output_path_for(path)
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nFigure -> {out_path}")


if __name__ == "__main__":
    main()
