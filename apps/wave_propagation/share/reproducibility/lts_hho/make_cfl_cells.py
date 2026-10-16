#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure cfl_cells.pdf of the LTS_HHO preprint, main.tex 6.4 (what sets the CFL limit of an HHO cell): CFL limit of ERK4 (random
data, k = 3 and 1) on the 16 x 16 square meshes of mesh_generation/lts_efficiency/generate_cfl_cells.py with one
modified cell (thin rectangle, small face, cut corner) or anisotropic cells, against the smallest face, the smallest
diameter and the smallest width (altitude) of the cells. Input: the cfl.json of analyze_cfl.py (efficiency/).
Usage: python3 make_cfl_cells.py <cfl.json> <figures>
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open(sys.argv[1]))
OUT = sys.argv[2]
plt.rcParams.update({"font.family": "serif", "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.2,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "mathtext.fontset": "cm"})
GRID = dict(color="#c8c8c4", lw=0.4, ls=":")
FAM = {"face": ("small face (square cells)", "#2b6cb0", "s"), "thin": ("thin rectangle", "#b7202e", "^"),
       "corner": ("cut corner (triangle)", "#dd6b20", "D"), "aniso": ("all cells $h\\times h/a$", "#2f855a", "o"),
       "slab": ("column of thin cells", "#805ad5", "v"), "none": ("squares", "#1a1a1a", "*")}
fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), sharey=True)
for ax, key, xl in ((axs[0], "face", "smallest face"), (axs[1], "diam", "smallest diameter"),
                    (axs[2], "w", r"smallest width $w$ (altitude)")):
    for fam, (lab, c, mk) in FAM.items():
        for k, fill in ((3, True), (1, False)):
            S = sorted([r for r in R if r["mesh"].split("_")[0] == fam and r["k"] == k], key=lambda r: r[key])
            if not S:
                continue
            ax.plot([r[key] for r in S], [r["dt"] for r in S], mk + ("-" if fam != "none" else ""), color=c, ms=3.2,
                    lw=0.7, mfc=c if fill else "white", label=f"{lab}" if k == 3 else None)
    if key == "w":
        w = np.array([5e-4, 0.1])
        ax.plot(w, 0.26 * w, "--", color="#6a6a66", lw=0.6)
        ax.plot(w, 0.115 * w, ":", color="#6a6a66", lw=0.6)
        ax.text(1.2e-3, 0.26 * 1.2e-3 * 1.6, r"$0.26\,w/c$", fontsize=6, color="#4a4a46")
        ax.text(1.2e-2, 0.115 * 1.2e-2 * 0.45, r"$0.115\,w/c$", fontsize=6, color="#4a4a46")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
    ax.set_xlabel(xl, labelpad=1); ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
axs[0].set_ylabel(r"CFL limit $\Delta t_{\rm crit}$ (ERK4)", labelpad=1)
h, l = axs[0].get_legend_handles_labels()
axs[1].legend(h, l, loc="upper left", frameon=False)
axs[1].text(0.03, 0.03, "filled: $k=3$, open: $k=1$", transform=axs[1].transAxes, va="bottom", fontsize=6)
fig.tight_layout(pad=0.2, w_pad=0.5)
fig.savefig(os.path.join(OUT, "cfl_cells.pdf"), bbox_inches="tight", pad_inches=0.02)
print("ok cfl_cells.pdf")
