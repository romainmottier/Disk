#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure cost_examples.pdf of the LTS_HHO preprint (efficiency, cost model (eq. gains) of the section "Efficiency and
cost"): two quadtree meshes of the unit square (32 x 32 squares of side h, refined towards interfaces, balanced 2:1),
and the gains over ERK4 given by the cost model, in the settings (1E) and (2E), of MLTS with all its levels, of MLTS
with two levels (all the refined cells at the finest step) and of LI. No computation: unknowns are counted per cell.
    top     extended refinement of moderate depth: a curved layer of thickness 3h, the squares meeting it refined once,
            containing a small inclusion of radius h/5, refined three times (bands of about one square per level),
            h_max / h_min = 8: MLTS with its four levels is the fastest scheme in the setting (1E), LI in the setting
            (2E)
    bottom  localized refinement of large depth: a small inclusion of radius h/64 refined down to h/256 (bands of
            about one square per level), h_max / h_min = 256: LI is the fastest scheme in both settings
rho: cost of a fine solve in products of the setting, measured on these meshes (analyze_cost_check.py): 16 in the
setting (1E), 6 in the setting (2E), whose products also solve for the faces of the whole mesh (three to four times
the cost per cell of a product of the setting (1E), for a solve about 1.5 times as expensive).
Rule of the cost model (two levels): LI is faster than MLTS when beta p_1 > 3 (1 + rho) / 4, p_1 = P_{L-1} the depth
of the refinement, rho growing with the extent of the refined region (fill-in of the factors).
Model (per coarse step, in products): ERK4 4 P_{L-1} N; MLTS 4 N + beta sum_{l>=1} m_l P_l |T_l|, m_{L-1} = 4, m_l = 7
otherwise, T_l: band of level l and one ring (1E), levels >= l and two rings (2E), beta = 1 (1E), 2 (2E); LI
4 N + 3 (1 + rho) |S|, S: the refined cells and their rings (one in (1E), two in (2E)). Level of a cell = its depth
(squares, uniform speed); rings by vertex adjacency, as the overlap rings of the code.
Usage: python3 make_cost_examples.py <figures> [<meshes>]
"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False, "font.size": 8,
                     "axes.linewidth": 0.6, "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42})
N0 = 32
H = 1.0 / N0


# ------------------------------------------------------------------------------------------------ quadtree meshes
def build(depth, features):
    """Leaves (i, j, d), square [i, i+1] x [j, j+1] h / 2^d. features: (distance(x0, x1, y0, y1), depth, W): a square
    of side s and depth d < depth is split when it meets the feature (W = 0) or lies within W s of it."""
    leaves = set()
    stack = [(i, j, 0) for i in range(N0) for j in range(N0)]
    while stack:
        i, j, d = stack.pop()
        s = H / 2 ** d
        x0, y0 = i * s, j * s
        split = False
        for dist, D, W in features:
            if d < D:
                g = dist(x0, x0 + s, y0, y0 + s)
                if g == 0.0 or g < W * s:
                    split = True
                    break
        if split:
            stack += [(2 * i + a, 2 * j + b, d + 1) for a in (0, 1) for b in (0, 1)]
        else:
            leaves.add((i, j, d))
    return balance(leaves)


def balance(leaves):
    """2:1 balance: a leaf having a neighbour two levels deeper is split, until no change."""
    leaves = set(leaves)
    changed = True
    while changed:
        changed = False
        dm = max(d for (_, _, d) in leaves)
        n = N0 * 2 ** dm

        def leaf_at_most(px, py, depth):
            return any((px >> (dm - e), py >> (dm - e), e) in leaves for e in range(depth, -1, -1))

        for (i, j, d) in list(leaves):
            if d >= dm - 1:
                continue
            sc = 2 ** (dm - d)
            x0, y0 = i * sc, j * sc
            pts = [p for t in range(0, sc, max(1, sc // 4))
                   for p in ((x0 - 1, y0 + t), (x0 + sc, y0 + t), (x0 + t, y0 - 1), (x0 + t, y0 + sc))]
            if any(0 <= px < n and 0 <= py < n and not leaf_at_most(px, py, d + 1) for (px, py) in pts):
                leaves.remove((i, j, d))
                leaves.update((2 * i + a, 2 * j + b, d + 1) for a in (0, 1) for b in (0, 1))
                changed = True
    return leaves


def neighbours(leaves):
    """Vertex adjacency (closed squares meeting), as the overlap rings of the code."""
    dm = max(d for (_, _, d) in leaves)
    n = N0 * 2 ** dm

    def find(px, py):
        for e in range(dm + 1):
            k = (px >> (dm - e), py >> (dm - e), e)
            if k in leaves:
                return k

    nb = {}
    for (i, j, d) in leaves:
        sc = 2 ** (dm - d)
        x0, y0 = i * sc, j * sc
        pts = [(x0 - 1, y0 - 1), (x0 + sc, y0 - 1), (x0 - 1, y0 + sc), (x0 + sc, y0 + sc)]
        for t in list(range(0, sc, max(1, sc // 4))) + [sc - 1]:
            pts += [(x0 - 1, y0 + t), (x0 + sc, y0 + t), (x0 + t, y0 - 1), (x0 + t, y0 + sc)]
        nb[(i, j, d)] = {find(px, py) for (px, py) in pts if 0 <= px < n and 0 <= py < n}
    return nb


def rings(A, nb, k):
    """Cells within k adjacency steps of the set A, A excluded."""
    seen, front = set(A), set(A)
    for _ in range(k):
        front = set().union(*(nb[c] for c in front)) - seen
        seen |= front
    return seen - set(A)


def gains(leaves, rho):
    nb = neighbours(leaves)
    L = max(d for (_, _, d) in leaves) + 1
    N = len(leaves)
    ge = [{c for c in leaves if c[2] >= l} for l in range(L)]
    band = [sum(1 for c in leaves if c[2] == l) for l in range(L)]
    P = lambda l: 2 ** l
    m = lambda l: 4 if l == L - 1 else 7
    G = {}
    for name, beta, k in (("1E", 1, 1), ("2E", 2, 2)):
        T = {l: band[l] + len(rings(ge[l], nb, 1)) if name == "1E" else len(ge[l]) + len(rings(ge[l], nb, 2))
             for l in range(1, L)}
        S = len(ge[1]) + len(rings(ge[1], nb, k))
        G[name] = {"ERK4": 1.0,
                   "MLTS, 2 levels": 4 * P(L - 1) * N / (4 * N + beta * 4 * P(L - 1) * S),
                   f"MLTS, {L} levels": 4 * P(L - 1) * N / (4 * N + beta * sum(m(l) * P(l) * T[l] for l in range(1, L))),
                   "LI": P(L - 1) * N / (N + 0.75 * (1 + rho[name]) * S)}
    return G, N, L, band


# ------------------------------------------------------------------------------------------------ interfaces
def layer(yc, t):
    """Distance (vertical, sampled) from a square to the layer |y - yc(x)| <= t / 2."""
    def dist(x0, x1, y0, y1):
        xs = np.linspace(x0, x1, 9)
        return float(np.maximum(np.maximum(yc(xs) - t / 2 - y1, y0 - yc(xs) - t / 2), 0.0).min())
    return dist


def disc(cx, cy, r):
    return lambda x0, x1, y0, y1: max(0.0, math.hypot(max(x0 - cx, 0.0, cx - x1), max(y0 - cy, 0.0, cy - y1)) - r)


YC = lambda x: 0.47 + 0.07 * np.sin(2 * np.pi * (x - 0.15))
TA, RA = 3 * H, H / 5
CA = (0.5, float(YC(0.5)))
RB, CB = H / 64, (0.5, 0.5)
RHO = {"1E": 16, "2E": 6}    # cost of a fine solve in products of the setting, measured on these meshes
EXAMPLES = (
    dict(title=r"(a) layer with an inclusion, $h_{\max}/h_{\min}=8$", depth=3, zoom=4.0, centre=CA,
         features=[(layer(YC, TA), 1, 0), (disc(*CA, RA), 3, 1)]),
    dict(title=r"(b) small inclusion, $h_{\max}/h_{\min}=256$", depth=8, zoom=6.0, centre=CB,
         features=[(disc(*CB, RB), 8, 1)], inset=(256.0, (0.55, 0.55, 0.42, 0.42))),
)
COL = {"ERK4": "#1a1a1a", "MLTS, 2 levels": "#dd6b20", "MLTS": "#b7202e", "LI": "#2b6cb0"}
light = lambda c: tuple(0.55 + 0.45 * x for x in matplotlib.colors.to_rgb(c))


def draw_mesh(ax, cells, color, lo, hi, lw):
    """The cells meeting [lo, hi], coloured by level (the interfaces are not drawn: the refined cells show them)."""
    sel = [c for c in cells if c[0] + c[2] > lo[0] and c[0] < hi[0] and c[1] + c[2] > lo[1] and c[1] < hi[1]]
    polys = [[(x, y), (x + h, y), (x + h, y + h), (x, y + h)] for x, y, h, _ in sel]
    ax.add_collection(PolyCollection(polys, facecolors=[color(l) for *_, l in sel], edgecolors="#3a3a38",
                                     linewidths=lw))
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def main():
    fig = plt.figure(figsize=(6.3, 4.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.45], hspace=0.40, wspace=0.32)
    for row, E in enumerate(EXAMPLES):
        leaves = build(E["depth"], E["features"])
        G, N, L, band = gains(leaves, RHO)
        cells = [(i * H / 2 ** d, j * H / 2 ** d, H / 2 ** d, d) for (i, j, d) in leaves]
        color = lambda l: "white" if l == 0 else plt.cm.Reds(0.22 + 0.58 * (l - 1) / max(1, L - 2))
        cx, cy = E["centre"]
        # whole mesh, then zoom on the refined region (and, in (b), on the inclusion)
        ax = fig.add_subplot(gs[row, 0])
        draw_mesh(ax, cells, color, (0, 0), (1, 1), 0.15)
        ax.set_title(E["title"], fontsize=8, loc="left", pad=3)
        ax.text(0.0, -0.08, rf"$N={N}$ cells, $L={L}$", transform=ax.transAxes, fontsize=6.5, va="top")
        ax = zax = fig.add_subplot(gs[row, 1])
        w = 0.5 / E["zoom"]
        draw_mesh(ax, cells, color, (cx - w, cy - w), (cx + w, cy + w), 0.3)
        ax.set_title(rf"zoom $\times{E['zoom']:.0f}$", fontsize=7, pad=3, loc="right")
        if "inset" in E:
            z, box = E["inset"]
            ins = ax.inset_axes(box)
            wi = 0.5 / z
            draw_mesh(ins, cells, color, (cx - wi, cy - wi), (cx + wi, cy + wi), 0.3)
            ins.text(0.05, 0.95, rf"$\times{z:.0f}$", transform=ins.transAxes, fontsize=6.5, va="top",
                     bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
            mark_inset(ax, ins, loc1=2, loc2=4, fc="none", ec="#6a6a66", lw=0.5)
        # legend of the levels (the deepest ones, too small to be seen, grouped) and of the interfaces, below the zoom
        shown = list(range(min(L, 4)))
        hs = [plt.Rectangle((0, 0), 1, 1, fc=color(l), ec="#3a3a38", lw=0.4) for l in shown]
        labs = [rf"$\ell={l}$" for l in shown]
        if L > 4:
            hs.append(plt.Rectangle((0, 0), 1, 1, fc=color(L - 1), ec="#3a3a38", lw=0.4)); labs.append(rf"$\ell=4..{L - 1}$")
        order = [r * 3 + c for c in range(3) for r in range(math.ceil(len(hs) / 3)) if r * 3 + c < len(hs)]
        zax.legend([hs[i] for i in order], [labs[i] for i in order], loc="upper center", bbox_to_anchor=(0.5, -0.03),
                   ncol=3, fontsize=6, frameon=False,
                   handlelength=0.9, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.0)
        # gains of the cost model
        ax = fig.add_subplot(gs[row, 2])
        names = list(G["1E"].keys())
        x = np.arange(len(names))
        for f, (setting, off) in enumerate((("1E", -0.19), ("2E", 0.19))):
            for i, nm in enumerate(names):
                c = COL["MLTS"] if nm.startswith("MLTS,") and "2 levels" not in nm else COL[nm]
                v = G[setting][nm]
                ax.bar(i + off, v, width=0.36, color=c if f == 0 else light(c), edgecolor="#1a1a1a", linewidth=0.4)
                ax.text(i + off, v * 1.08, f"{v:.0f}" if v >= 10 else f"{v:.2f}" if v < 1 else f"{v:.1f}", ha="center",
                        va="bottom", fontsize=5.5)
        ax.set_yscale("log"); ax.axhline(1, color="#6a6a66", lw=0.6, ls="--", zorder=0)
        ax.set_xticks(x); ax.set_xticklabels([nm.replace(", ", ",\n") for nm in names], fontsize=6.5)
        ax.set_ylabel("model gain over ERK4", labelpad=1)
        top = max(max(G[s].values()) for s in G)
        ax.set_ylim(0.3, top * 4)
        ax.grid(True, axis="y", which="major", ls=":", lw=0.4, color="#b4b4b0")
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        if row == 0:
            ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#5a5a56", ec="#1a1a1a", lw=0.4, label="setting (1E)"),
                               plt.Rectangle((0, 0), 1, 1, fc=light("#5a5a56"), ec="#1a1a1a", lw=0.4, label="setting (2E)")],
                      loc="upper left", fontsize=6.5, frameon=False, handlelength=1.2)
        print(E["title"], "N", N, "L", L, "cells per level", band, "rho", RHO)
        for s in G:
            print("  ", s, {k: round(v, 2) for k, v in G[s].items()})
    fig.savefig(os.path.join(OUT, "cost_examples.pdf"), bbox_inches="tight", pad_inches=0.02)
    print("ok cost_examples.pdf")
    # optional: the two meshes in the "polygons" format of wave_propagation (hanging nodes as polygon vertices), for the
    # measurements of the gains (cost_a.txt: mesh.level 0, cost_b.txt: mesh.level 1)
    if len(sys.argv) > 2:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src",
                                        "mesh_generation", "lts_efficiency"))
        from generate_refined_quadtree import write
        os.makedirs(sys.argv[2], exist_ok=True)
        for name, E in zip(("cost_a", "cost_b"), EXAMPLES):
            path = os.path.join(sys.argv[2], name + ".txt")
            nc, per_depth, nv = write(path, build(E["depth"], E["features"]), N0, "make_cost_examples.py: " + E["title"])
            print(path, nc, "cells, per depth", dict(sorted(per_depth.items())), "at most", nv, "vertices per cell")


if __name__ == "__main__":
    main()
