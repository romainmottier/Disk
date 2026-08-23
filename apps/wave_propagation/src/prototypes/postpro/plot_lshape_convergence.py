#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Log-log convergence plot for ERK4_LTS_Lshape_conv_test.hpp: pressure L2
error (vs the finest/reference graded mesh) as a function of h_max, for the
HHO + ERK4-LTS solve on the graded L-shaped reentrant-corner mesh --
the first-order-formulation analogue of Fig. 10 in Grote/Michel/Sauter
(arXiv:2005.13350, Sec. 4.3).

Usage:
    python3 plot_lshape_convergence.py [file]

    file : the "lshape_convergence_k_<k>.txt" log written by the C++
           prototype (default: first match of lshape_convergence_k_*.txt
           in the current directory)
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
    # columns: N  h_min  L2_error_pressure  h_max
    N, h, err = data[:, 0], data[:, 1], data[:, 2]
    order = np.argsort(h)
    return N[order], h[order], err[order]


def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        matches = sorted(glob.glob("lshape_convergence_k_*.txt"))
        if not matches:
            print("Aucun fichier lshape_convergence_k_*.txt trouve.")
            sys.exit(1)
        path = matches[0]

    N, h, err = load(path)
    print(f"Charge depuis '{path}':")
    for ni, hi, ei in zip(N, h, err):
        print(f"  N={int(ni):4d}   h_max={hi:.6g}   L2_error(pressure)={ei:.6g}")

    p_fit, c_fit = np.polyfit(np.log(h), np.log(err), 1)
    print(f"\nPente ajustee (log(err) = p*log(h) + c) : p = {p_fit:+.3f}")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.loglog(h, err, "o-", color="#2a78d6", lw=2, ms=6, zorder=3,
              label="L2 error (pressure) vs reference")

    h_ref = np.array([h.min(), h.max()])
    anchor_h, anchor_e = h[0], err[0]
    ax.loglog(h_ref, anchor_e * (h_ref / anchor_h) ** p_fit, "--", color=MUTED, lw=1,
              label=rf"fit: $h^{{{p_fit:.2f}}}$")

    ax.set_xlabel(r"$h_{\min}$ (corner cell size)", color=INK, fontsize=11)
    ax.set_ylabel(r"$\|p_h - p_{\mathrm{ref}}\|_{L^2}$", color=INK, fontsize=11)
    ax.set_title("L-shape reentrant corner: HHO+ERK4-LTS convergence\n(graded mesh, pressure field vs finest-mesh reference)",
                 color=INK, fontsize=10)
    ax.grid(True, which="both", color=GRID, lw=0.4, alpha=0.9)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend = ax.legend(fontsize=9, frameon=True, facecolor=SURFACE, edgecolor=AXIS, labelcolor=INK)
    legend.get_frame().set_linewidth(0.8)

    fig.tight_layout()
    out_path = os.path.splitext(path)[0] + ".png"
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nFigure -> {out_path}")


if __name__ == "__main__":
    main()
