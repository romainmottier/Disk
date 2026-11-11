#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure cfl_width.pdf of the LTS_HHO preprint (what sets the CFL limit of an HHO cell): CFL limit of ERK4 (random data)
in the setting (1E) with k = 3 and in the setting (2E) with k = 2 (cell diameter h_T in the 1/h stabilization,
space.stabilization_size = "cell"), three tests on the 16 x 16 squares, one geometry each, each varying one size:
    (i)   small faces      a vertex inserted on an edge of a square: a face of length eps h (face_<eps>), the diameters
                           and widths unchanged
    (ii)  stretched cells  all the cells rectangles h x a h (aniso_<a>, a = 2 .. 16): width h, diameter up to 16 h
    (iii) thin cells       all the cells rectangles h x eps h (aniso_<eps>, eps = 1/2 .. 1/64): width eps h, diameter
                           about h
    one panel per test: the CFL limit as a function of the size varied by the test (smallest value over the mesh), for
    both settings ((1E) filled markers, (2E) open ones), star: the square mesh, dashed: the CFL limit proportional to
    that size through the square mesh (the sketches of the three tests above the panels are drawn in TikZ in main.tex)
Also prints, for the text, dt c / w on the meshes with an isolated thin cell (thin_<eps>, corner_<eps>, circle_<i>).
Inputs: the cfl.json of analyze_cfl.py for efficiency/run_cfl_cells.sh (setting (1E)), run_cfl_random.sh
(circle-fitted meshes, setting (1E), printed only) and run_cfl_2e.sh (setting (2E), inputs/cfl_cells2.lua).
Usage: python3 make_cfl_width.py <1E families cfl.json> <1E circles cfl.json> <2E cfl.json> <figures>
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FJ, CJ, EJ, OUT = sys.argv[1:5]


def kind(r):
    f, _, e = r["mesh"].partition("_")
    if f == "aniso":
        return "stretch" if float(e) > 1 else "thin"
    return {"none": "none", "face": "face", "thin": "isolated", "corner": "isolated", "circle": "isolated"}.get(f)


def load(files, k):
    return [dict(r, kind=kind(r)) for f in files for r in json.load(open(f)) if r["k"] == k and kind(r)]


ROWS = ((r"Setting (1E): first-order formulation, $k=3$, $\eta=0.8$", load((FJ, CJ), 3)),
        (r"Setting (2E): second-order formulation, $k=2$, $\eta=3\eta^\star$", load((EJ,), 2)))
plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False, "font.size": 8.5,
                     "axes.linewidth": 0.6, "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42})
GRID = dict(ls=":", lw=0.4, color="#b4b4b0")
h = 1 / 16
# test -> (title, color, marker, size, size varied, label)
TESTS = {
    "face":    ("(i) small faces", "#2b6cb0", "s", 4.0, "face", "smallest face"),
    "stretch": ("(ii) stretched cells", "#2f855a", "o", 4.0, "diam", r"smallest diameter $h_T$"),
    "thin":    ("(iii) thin cells", "#b7202e", "^", 4.5, "w", r"smallest width $w_T$"),
}

fig = plt.figure(figsize=(6.3, 2.15))
grid = fig.add_gridspec(1, 3, wspace=0.1)
# the CFL limit against the size varied by each test, both settings in the same panel: (1E) filled markers and star,
# (2E) open ones; dashed: proportional to the size through the square mesh of each setting
dts = [r["dt"] for _, R in ROWS for r in R if r["kind"] != "isolated"]
ylim = (min(dts) / 2.5, max(dts) * 4)
SHORT = ("(1E)", "(2E)")
for i, (f, (lab, c, mk, ms, m, xlab)) in enumerate(TESTS.items()):
    ax = fig.add_subplot(grid[0, i])
    txt = []
    for row, (title, R) in enumerate(ROWS):
        sq = next(r for r in R if r["kind"] == "none")
        S = [r for r in R if r["kind"] == f]
        xs = [r[m] for r in S + [sq]]
        x = np.geomspace(min(xs) / 2.5, max(xs) * 2.5, 50)
        ax.plot(x, sq["dt"] * x / sq[m], ls="--", color="#8a8a86", lw=0.7, zorder=1)
        ax.plot([r[m] for r in S], [r["dt"] for r in S], mk, color=c, ms=ms, mew=0.9, ls="none", zorder=3,
                mfc=c if row == 0 else "white")
        ax.plot([sq[m]], [sq["dt"]], "*", color="#1a1a1a", ms=6.5, mew=0.7, zorder=4,
                mfc="#1a1a1a" if row == 0 else "white")
        if f == "thin":
            v = np.array([r["dt"] / r["w"] for r in S + [sq]])
            txt.append(rf"{SHORT[row]}: $\Delta t_{{\rm crit}}\,c/w_T\in[{v.min():.3f},{v.max():.3f}]$")
        else:
            v = np.array([r["dt"] for r in S]) / sq["dt"]
            txt.append(rf"{SHORT[row]}: ratio to $\star$ ${v.min():.2f}$ to ${v.max():.2f}$")
    ax.text(0.95, 0.05, "\n".join(txt), transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5,
            linespacing=1.3, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#9a9a96", lw=0.5))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(*ylim); ax.set_xlim(x[0] * 1.5, x[-1] / 1.5)
    ax.set_xlabel(xlab, labelpad=1, fontsize=7.5)
    if i == 0:
        ax.set_ylabel(r"$\Delta t_{\rm crit}$", labelpad=1, fontsize=7.5)
    else:
        ax.tick_params(labelleft=False)
    ax.grid(True, which="major", **GRID); ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
# legend of the settings, of the stars and of the dashed lines, below the panels
G = "#5a5a56"
fig.legend(handles=[plt.Line2D([], [], marker="o", color=G, ms=4.5, ls="none"),
                    plt.Line2D([], [], marker="o", color=G, mfc="white", mew=0.9, ms=4.5, ls="none"),
                    plt.Line2D([], [], marker="*", color="#1a1a1a", ms=6.5, ls="none"),
                    plt.Line2D([], [], ls="--", color="#8a8a86", lw=0.7)],
           labels=[r"setting (1E): first-order, $k=3$, $\eta=0.8$",
                   r"setting (2E): second-order, $k=2$, $\eta=3\eta^\star$",
                   r"square mesh ($16\times16$ squares)",
                   r"$\Delta t_{\rm crit}\propto$ size, through the square mesh"],
           loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False, fontsize=7.5, handlelength=2.2,
           columnspacing=2.0)
fig.savefig(os.path.join(OUT, "cfl_width.pdf"), bbox_inches="tight", pad_inches=0.02)
print("ok cfl_width.pdf")
for title, R in ROWS:
    sq = next(r for r in R if r["kind"] == "none")
    print(title, "square", f"{sq['dt']:.4e}", f"dt c/w {sq['dt'] / sq['w']:.4f}")
    for f in list(TESTS) + ["isolated"]:
        T = [r for r in R if r["kind"] == f]
        v = np.array([r["dt"] / r["w"] for r in T]); u = np.array([r["dt"] for r in T]) / sq["dt"]
        print(f"  {f:9} {len(T):2} meshes: dt/dt_square in [{u.min():.3f}, {u.max():.3f}], dt c/w in [{v.min():.4f}, "
              f"{v.max():.4f}]")
    v = np.array([r["dt"] / r["w"] for r in R if r["kind"] != "isolated"])
    print(f"  three tests + square: dt c/w in [{v.min():.4f}, {v.max():.4f}], spread x{v.max() / v.min():.3f}")
