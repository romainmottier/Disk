#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure cost_examples.pdf of the LTS_HHO preprint (section "Efficiency and cost"): two quadtree meshes of the unit
square (32 x 32 squares of side h, refined towards interfaces, balanced 2:1), and the gains over ERK4, in the settings
(1E) and (2E), of MLTS with all its levels, of MLTS with two levels (all the refined cells at the finest step) and of
LI. Two modes:
    --runs <cost> [<cost_r2> <cost_r3>]   gains MEASURED by the runs cost/ of article/jobs.txt (one thread, n and 2n
            coarse steps, time per coarse step = difference of the wall times: setup and factorization excluded; best
            of the repetitions given), with the peak memory relative to ERK4 below them (maximum resident set size of
            all the runs of a scheme); also writes the same split in two figures, those of the preprint since
            2026-11-13: cost_meshes.pdf (the meshes and their zooms, Fig. cost-meshes) and cost_gains.pdf (gains and
            peak memories, Fig. cost-gains)
    --formulations   with --runs, also cost_formulations.pdf (2026-11-15, Fig. cost-formulations): wall time per coarse
            step and memory of the schemes in the two settings, as absolute values (same runs, same time steps in both
            settings), for the comparison of the formulations; off by default
    --paper-notation   labels of cost_gains.pdf and cost_formulations.pdf as in the paper (2026-11-17); off by default
    without --runs   gains of the cost model below (no computation: unknowns are counted per cell)
    top     extended refinement of moderate depth: a curved layer of thickness 3h, the squares meeting it refined once,
            containing a small inclusion of radius h/5, refined three times (bands of about one square per level),
            h_max / h_min = 8: MLTS with its four levels is the fastest scheme in the setting (1E), LI in the setting
            (2E)
    bottom  localized refinement of large depth: a small inclusion of radius h/64 refined down to h/256 (bands of
            about one square per level), h_max / h_min = 256: LI is the fastest scheme in both settings
rho: cost of a fine solve in products of the setting, measured on these meshes (analyze_cost_check.py): 16 in the
setting (1E), 6 in the setting (2E), whose products also solve for the faces of the whole mesh (three to four times
the cost per cell of a product of the setting (1E), for a solve about 1.5 times as expensive). These values come from
SparseLU in (1E) and Pardiso with two refinement steps in (2E); with Pardiso in both settings and no iterative
refinement (campaign of 2026-11-15, article/pardiso), rho = 8 to 12 in (1E) and 3.5 to 5 in (2E), and a solve costs
1.2 to 2.1 times as much in (2E) as in (1E) (RHO below, used only without --runs, is kept).
Rule of the cost model (two levels): LI is faster than MLTS when beta p_1 > 3 (1 + rho) / 4, p_1 = P_{L-1} the depth
of the refinement, rho growing with the extent of the refined region (fill-in of the factors).
Model (per coarse step, in products): ERK4 4 P_{L-1} N; MLTS 4 N + beta sum_{l>=1} m_l P_l |T_l|, m_{L-1} = 4, m_l = 7
otherwise, T_l: band of level l and one ring (1E), levels >= l and two rings (2E), beta = 1 (1E), 2 (2E); LI
4 N + 3 (1 + rho) |S|, S: the refined cells and their rings (one in (1E), two in (2E)). Level of a cell = its depth
(squares, uniform speed); rings by vertex adjacency, as the overlap rings of the code.
Usage: python3 make_cost_examples.py <figures> [<meshes>] [--formulations] [--paper-notation] [--memory <algmem runs>]
                                      [--runs <cost> <cost_r2> <cost_r3>]
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
SCHEMES = (("ERK4", "erk"), ("MLTS, 2 levels", "mlts2"), ("MLTS", "mlts"), ("LI", "li"))   # name, run of jobs.txt
MEM = None   # --memory <runs>: algebraic memory (output.algebraic_memory) of the runs <runs>/cost/<run>/log.txt
PAPER = False   # --paper-notation (2026-11-17, off by default): ERK(4) and "CPU time" in cost_gains/cost_formulations.pdf


def shown(label):
    """Label as printed: ERK4 written ERK(4) with --paper-notation."""
    return label.replace("ERK4", "ERK(4)") if PAPER else label


def measured(runs, tag, L):
    """Measured gains over ERK4 and peak memories relative to ERK4 of mesh tag ("a" or "b") in the two settings, from
    the runs <rep>/o{1,2}_<tag>_<scheme>_{n,2n} (article/jobs.txt, cost/): time per coarse step = (W(2n) - W(n)) / n
    with the smallest wall times W of the repetitions, peak memory = largest maximum resident set size of the runs."""
    import re

    def runs_of(run):
        out = []
        for d in runs:
            f = os.path.join(d, run, "time.txt")
            if os.path.exists(f):
                m = re.findall(r"([0-9.]+) s ([0-9]+) KB", open(f).read())
                if m:
                    out.append((float(m[-1][0]), int(m[-1][1])))
        if not out:
            raise RuntimeError("no run " + run + " in " + " ".join(runs))
        return out

    def n_coarse(run):
        for d in runs:
            f = os.path.join(d, run, "case.lua")
            if os.path.exists(f):
                return round(float(re.findall(r"time.final_time = ([0-9.e-]+)", open(f).read())[-1]) * 400)

    G, Mem = {}, {}
    for setting, s in (("1E", "o1"), ("2E", "o2")):
        t, mem = {}, {}
        for name, sc in SCHEMES:
            a, b = runs_of(f"{s}_{tag}_{sc}_n"), runs_of(f"{s}_{tag}_{sc}_2n")
            t[name] = (min(x[0] for x in b) - min(x[0] for x in a)) / n_coarse(f"{s}_{tag}_{sc}_n")
            mem[name] = max(x[1] for x in a + b)
            if MEM:   # algebraic memory of the run (line "algebraic memory: ... total X MB" of its log.txt)
                f = os.path.join(MEM, "cost", f"{s}_{tag}_{sc}_n", "log.txt")
                mem[name] = float(re.search(r"algebraic memory: .* total ([0-9.]+) MB", open(f).read()).group(1))
        key = lambda nm: f"MLTS, {L} levels" if nm == "MLTS" else nm
        G[setting] = {key(nm): t["ERK4"] / t[nm] for nm, _ in SCHEMES}
        Mem[setting] = {key(nm): mem[nm] / mem["ERK4"] for nm, _ in SCHEMES}
    return G, Mem


def measured_absolute(runs, tag, L):
    """Absolute values behind measured(), for the comparison of the formulations (--formulations): wall time per coarse
    step in ms (same runs and same rule: (W(2n) - W(n)) / n with the smallest wall times of the repetitions) and memory in
    MiB (algebraic memory of the run with --memory, otherwise the largest maximum resident set size of the runs)."""
    import re

    def runs_of(run):
        out = []
        for d in runs:
            f = os.path.join(d, run, "time.txt")
            if os.path.exists(f):
                m = re.findall(r"([0-9.]+) s ([0-9]+) KB", open(f).read())
                if m:
                    out.append((float(m[-1][0]), int(m[-1][1])))
        if not out:
            raise RuntimeError("no run " + run + " in " + " ".join(runs))
        return out

    def n_coarse(run):
        for d in runs:
            f = os.path.join(d, run, "case.lua")
            if os.path.exists(f):
                return round(float(re.findall(r"time.final_time = ([0-9.e-]+)", open(f).read())[-1]) * 400)

    key = lambda nm: f"MLTS, {L} levels" if nm == "MLTS" else nm
    T, Mem = {}, {}
    for setting, s in (("1E", "o1"), ("2E", "o2")):
        T[setting], Mem[setting] = {}, {}
        for name, sc in SCHEMES:
            a, b = runs_of(f"{s}_{tag}_{sc}_n"), runs_of(f"{s}_{tag}_{sc}_2n")
            T[setting][key(name)] = 1e3 * (min(x[0] for x in b) - min(x[0] for x in a)) / n_coarse(f"{s}_{tag}_{sc}_n")
            mem = max(x[1] for x in a + b) / 1024.0
            if MEM:   # "algebraic memory: ... total X MB" of its log.txt (MB = 2^20 bytes)
                f = os.path.join(MEM, "cost", f"{s}_{tag}_{sc}_n", "log.txt")
                mem = float(re.search(r"algebraic memory: .* total ([0-9.]+) MB", open(f).read()).group(1))
            Mem[setting][key(name)] = mem
    return T, Mem


def draw_mesh(ax, cells, color, lo, hi, lw):
    """The cells meeting [lo, hi], coloured by level (the interfaces are not drawn: the refined cells show them)."""
    sel = [c for c in cells if c[0] + c[2] > lo[0] and c[0] < hi[0] and c[1] + c[2] > lo[1] and c[1] < hi[1]]
    polys = [[(x, y), (x + h, y), (x + h, y + h), (x, y + h)] for x, y, h, _ in sel]
    ax.add_collection(PolyCollection(polys, facecolors=[color(l) for *_, l in sel], edgecolors="#3a3a38",
                                     linewidths=lw))
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def mesh_panels(ax, zax, E, leaves, L, N, short=False):
    """Whole mesh (ax) and zoom on the refined region (zax, with an inset in (b)), colored by the levels, and the
    legend of the levels below the zoom; short: title without the ratio h_max / h_min, given below the mesh."""
    cells = [(i * H / 2 ** d, j * H / 2 ** d, H / 2 ** d, d) for (i, j, d) in leaves]
    color = lambda l: "white" if l == 0 else plt.cm.Reds(0.22 + 0.58 * (l - 1) / max(1, L - 2))
    cx, cy = E["centre"]
    draw_mesh(ax, cells, color, (0, 0), (1, 1), 0.15)
    title, ratio = E["title"].split(", ")
    ax.set_title(title if short else E["title"], fontsize=8, loc="left", pad=3)
    ax.text(0.0, -0.08, rf"$N={N}$ cells, $L={L}$" + ("\n" + ratio if short else ""), transform=ax.transAxes,
            fontsize=6.5, va="top")
    w = 0.5 / E["zoom"]
    draw_mesh(zax, cells, color, (cx - w, cy - w), (cx + w, cy + w), 0.3)
    zax.set_title(rf"zoom $\times{E['zoom']:.0f}$", fontsize=7, pad=3, loc="right")
    if "inset" in E:
        z, box = E["inset"]
        ins = zax.inset_axes(box)
        wi = 0.5 / z
        draw_mesh(ins, cells, color, (cx - wi, cy - wi), (cx + wi, cy + wi), 0.3)
        ins.text(0.05, 0.95, rf"$\times{z:.0f}$", transform=ins.transAxes, fontsize=6.5, va="top",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
        mark_inset(zax, ins, loc1=2, loc2=4, fc="none", ec="#6a6a66", lw=0.5)
    # legend of the levels (the deepest ones, too small to be seen, grouped), below the zoom
    shown = list(range(min(L, 4)))
    hs = [plt.Rectangle((0, 0), 1, 1, fc=color(l), ec="#3a3a38", lw=0.4) for l in shown]
    labs = [rf"$\ell={l}$" for l in shown]
    if L > 4:
        hs.append(plt.Rectangle((0, 0), 1, 1, fc=color(L - 1), ec="#3a3a38", lw=0.4)); labs.append(rf"$\ell=4..{L - 1}$")
    order = [r * 3 + c for c in range(3) for r in range(math.ceil(len(hs) / 3)) if r * 3 + c < len(hs)]
    zax.legend([hs[i] for i in order], [labs[i] for i in order], loc="upper center", bbox_to_anchor=(0.5, -0.03),
               ncol=3, fontsize=6, frameon=False,
               handlelength=0.9, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.0)


def main():
    global MEM, PAPER
    args, runs = sys.argv[1:], None
    formulations = "--formulations" in args   # also cost_formulations.pdf (absolute values, both settings)
    if formulations:
        args.remove("--formulations")
    if "--paper-notation" in args:   # labels in the notation of the paper (2026-11-17, off by default)
        PAPER = True
        args.remove("--paper-notation")
    if "--memory" in args:    # algebraic memory of the runs instead of the peak memory (campaign lts_work/algmem)
        k = args.index("--memory")
        MEM, args = args[k + 1], args[:k] + args[k + 2:]
    if "--runs" in args:
        k = args.index("--runs")
        args, runs = args[:k], args[k + 1:]
    out = args[0] if args else "."
    data = []
    for row, E in enumerate(EXAMPLES):
        leaves = build(E["depth"], E["features"])
        G, N, L, band = gains(leaves, RHO)
        Mem = None
        if runs:
            G, Mem = measured(runs, "ab"[row], L)
        data.append((E, leaves, G, Mem, N, L))
        print(E["title"], "N", N, "L", L, "cells per level", band, "rho", RHO, "(measured)" if runs else "(model)")
        for s in G:
            print("  ", s, "gains", {k: round(v, 2) for k, v in G[s].items()})
            if Mem:
                print("  ", s, ("algebraic" if MEM else "peak") + " memory / ERK4", {k: round(v, 2) for k, v in Mem[s].items()})
    # combined figure: one row per mesh (mesh, zoom, gains and, measured, peak memory below them)
    fig = plt.figure(figsize=(6.3, 5.0 if runs else 4.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.45], hspace=0.40, wspace=0.32)
    for row, (E, leaves, G, Mem, N, L) in enumerate(data):
        mesh_panels(fig.add_subplot(gs[row, 0]), fig.add_subplot(gs[row, 1]), E, leaves, L, N)
        if runs:
            sub = gs[row, 2].subgridspec(2, 1, height_ratios=[1.8, 1], hspace=0.08)
            ax, axm = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
        else:
            ax, axm = fig.add_subplot(gs[row, 2]), None
        bar_panels(ax, axm, G, Mem, runs, row == 0)
    fig.savefig(os.path.join(out, "cost_examples.pdf"), bbox_inches="tight", pad_inches=0.02)
    print("ok cost_examples.pdf")
    if runs:
        # the same, split in the two figures of the preprint: the meshes, then the gains and peak memories
        mesh_figure(data, os.path.join(out, "cost_meshes.pdf"))
        print("ok cost_meshes.pdf")
        gains_figure(data, os.path.join(out, "cost_gains.pdf"))
        print("ok cost_gains.pdf")
        if formulations:
            data_abs = []
            for row, (E, leaves, G, Mem, N, L) in enumerate(data):
                T, M = measured_absolute(runs, "ab"[row], L)
                data_abs.append((E, T, M))
                for s in T:
                    print("  ", E["title"].split(", ")[0], s, "ms per coarse step", {k: round(v, 3) for k, v in T[s].items()},
                          "memory (MiB)", {k: round(v, 1) for k, v in M[s].items()})
            formulations_figure(data_abs, os.path.join(out, "cost_formulations.pdf"))
            print("ok cost_formulations.pdf")
    # optional: the two meshes in the "polygons" format of wave_propagation (hanging nodes as polygon vertices), for the
    # measurements of the gains (cost_a.txt: mesh.level 0, cost_b.txt: mesh.level 1)
    if len(args) > 1:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src",
                                        "mesh_generation", "lts_efficiency"))
        from generate_refined_quadtree import write
        os.makedirs(args[1], exist_ok=True)
        for name, E in zip(("cost_a", "cost_b"), EXAMPLES):
            path = os.path.join(args[1], name + ".txt")
            nc, per_depth, nv = write(path, build(E["depth"], E["features"]), N0, "make_cost_examples.py: " + E["title"])
            print(path, nc, "cells, per depth", dict(sorted(per_depth.items())), "at most", nv, "vertices per cell")


def mesh_figure(data, path):
    """cost_meshes.pdf: for (a) and (b), the whole mesh (zoom window outlined) and the zoom, colored by the levels of
    MLTS, a short title centered over each pair and a discrete color bar of the levels below it (style of the other
    figures of the preprint: small centered titles, no text in the panels but the zoom factors)."""
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.cm import ScalarMappable
    fig = plt.figure(figsize=(6.3, 2.0))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.14, 1, 1], wspace=0.06)
    box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
    for row, (E, leaves, G, Mem, N, L) in enumerate(data):
        ax, zax = fig.add_subplot(gs[0, 3 * row]), fig.add_subplot(gs[0, 3 * row + 1])
        cells = [(i * H / 2 ** d, j * H / 2 ** d, H / 2 ** d, d) for (i, j, d) in leaves]
        color = lambda l: "white" if l == 0 else plt.cm.Reds(0.22 + 0.58 * (l - 1) / max(1, L - 2))
        cx, cy = E["centre"]
        w = 0.5 / E["zoom"]
        draw_mesh(ax, cells, color, (0, 0), (1, 1), 0.15)
        ax.add_patch(plt.Rectangle((cx - w, cy - w), 2 * w, 2 * w, fill=False, ec="#1a1a1a", lw=0.7, zorder=3))
        draw_mesh(zax, cells, color, (cx - w, cy - w), (cx + w, cy + w), 0.3)
        zax.text(0.04, 0.96, rf"$\times{E['zoom']:.0f}$", transform=zax.transAxes, fontsize=6.5, va="top", bbox=box)
        if "inset" in E:
            z, ibox = E["inset"]
            ins = zax.inset_axes(ibox)
            wi = 0.5 / z
            draw_mesh(ins, cells, color, (cx - wi, cy - wi), (cx + wi, cy + wi), 0.3)
            ins.text(0.05, 0.95, rf"$\times{z:.0f}$", transform=ins.transAxes, fontsize=6, va="top", bbox=box)
            mark_inset(zax, ins, loc1=2, loc2=4, fc="none", ec="#6a6a66", lw=0.5)
        p0, p1 = ax.get_position(), zax.get_position()
        fig.text((p0.x0 + p1.x1) / 2, p0.y1 + 0.035, E["title"].split(", ")[0], ha="center", va="bottom", fontsize=7)
        width = p1.x1 - p0.x0
        cax = fig.add_axes([p0.x0 + 0.2 * width, p0.y0 - 0.10, 0.6 * width, 0.032])
        cb = fig.colorbar(ScalarMappable(norm=BoundaryNorm(np.arange(-0.5, L), L),
                                         cmap=ListedColormap([color(l) for l in range(L)])),
                          cax=cax, orientation="horizontal", ticks=range(L))
        cb.outline.set_linewidth(0.4)
        cb.ax.tick_params(labelsize=6, length=1.5, width=0.4, pad=1)
        cb.ax.set_xlabel(r"level $\ell$", fontsize=6.5, labelpad=1)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)


def gains_figure(data, path):
    """cost_gains.pdf: gains over ERK4 (top, log scale) and peak memory relative to ERK4 (bottom), (a) left and (b)
    right, settings (1E) dark and (2E) light (style of form.pdf: small centered titles, framed legend)."""
    plt.rcParams["hatch.linewidth"] = 0.4
    hatch = {0: None, 1: "////"}      # setting (2E) hatched: the two bars stay distinct in black and white
    mmax = max(v for (E, leaves, G, Mem, N, L) in data for s in Mem for v in Mem[s].values())
    mbig = mmax > 1.4                 # memories far from that of ERK4 (algebraic memory): linear scale up to the largest
    fig, axs = plt.subplots(2, 2, figsize=(6.3, 3.9), gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.10, wspace=0.16))
    for col, (E, leaves, G, Mem, N, L) in enumerate(data):
        ax, axm = axs[0, col], axs[1, col]
        names = list(G["1E"].keys())
        x = np.arange(len(names))
        for f, (setting, off) in enumerate((("1E", -0.19), ("2E", 0.19))):
            for i, nm in enumerate(names):
                c = COL["MLTS"] if nm.startswith("MLTS,") and "2 levels" not in nm else COL[nm]
                c = c if f == 0 else light(c)
                v = G[setting][nm]
                ax.bar(i + off, v, width=0.36, color=c, edgecolor="#1a1a1a", linewidth=0.4, hatch=hatch[f])
                # --paper-notation (2026-11-19): a gain below 1 is labelled above the dashed ratio-1 line, not across it
                ax.text(i + off, (max(v, 1.0) if PAPER else v) * 1.10,
                        f"{v:.0f}" if v >= 10 else f"{v:.2f}" if v < 1 else f"{v:.1f}", ha="center", va="bottom", fontsize=6)
                v = Mem[setting][nm]
                axm.bar(i + off, v, width=0.36, color=c, edgecolor="#1a1a1a", linewidth=0.4, hatch=hatch[f])
                axm.text(i + off, v + (0.02 * mmax if mbig else 0.03), f"{v:.1f}" if mbig else f"{v:.2f}", ha="center",
                         va="bottom", fontsize=5.5)
        from matplotlib.ticker import LogLocator, FuncFormatter, NullLocator
        def plain_log(a, lo, hi):   # log axis with the ticks 1, 2, 5 per decade, written as plain numbers
            a.set_yscale("log"); a.set_ylim(lo, hi)
            a.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            a.yaxis.set_minor_locator(NullLocator())
        ax.axhline(1, color="#6a6a66", lw=0.6, ls="--", zorder=0)
        plain_log(ax, 0.5, max(max(G[s].values()) for s in G) * 3)
        ax.set_title(E["title"].split(", ")[0], fontsize=7, pad=3)
        if mbig:   # linear scale (user, 2026-11-13: no log scale on the memories)
            axm.axhline(1, color="#6a6a66", lw=0.6, ls="--", zorder=0)
            axm.set_ylim(0, mmax * 1.25)
        else:
            axm.set_ylim(0, 1.45); axm.set_yticks([0, 0.5, 1])
        for a in (ax, axm):
            a.set_xticks(x); a.set_xlim(-0.6, len(names) - 0.4)
            a.grid(True, axis="y", which="major", ls=":", lw=0.4, color="#b4b4b0")
            a.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        ax.tick_params(labelbottom=False)
        axm.set_xticklabels([shown(nm).replace(", ", ",\n") for nm in names], fontsize=7)
        if col == 0:
            ax.set_ylabel(shown("gain over ERK4"), labelpad=2)
            axm.set_ylabel(shown(("algebraic" if MEM else "peak") + " memory / ERK4"), labelpad=2)
    axs[0, 0].legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#5a5a56", ec="#1a1a1a", lw=0.4, label="setting (1E)"),
                              plt.Rectangle((0, 0), 1, 1, fc=light("#5a5a56"), ec="#1a1a1a", lw=0.4, hatch=hatch[1],
                                            label="setting (2E)")],
                     loc="upper left", fontsize=6.5, framealpha=0.92, edgecolor="#9a9a96", handlelength=1.2,
                     borderpad=0.35, labelspacing=0.25)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)


def formulations_figure(data_abs, path):
    """cost_formulations.pdf (--formulations): wall time per coarse step (top, ms, log scale) and memory (bottom, MiB,
    linear scale) of ERK4, MLTS with two and with all its levels and LI, settings (1E) dark and (2E) light hatched, meshes
    (a) left and (b) right; absolute values, for the comparison of the formulations (style of cost_gains.pdf)."""
    from matplotlib.ticker import LogLocator, FuncFormatter, NullLocator
    plt.rcParams["hatch.linewidth"] = 0.4
    hatch = {0: None, 1: "////"}
    mmax = max(v for (E, T, M) in data_abs for s in M for v in M[s].values())
    fig, axs = plt.subplots(2, 2, figsize=(6.3, 3.9), gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.10, wspace=0.20))
    # --paper-notation (2026-11-18): bar labels with three significant digits, so that the ratios quoted in the text can
    # be recovered from them; default unchanged (two significant digits). Fourth review (2026-11-18): four significant
    # digits from 100 on (238.8 / 122.5 = 1.95, the ratio 1.9 of the text, where 239 / 122 gave 1.96), units in brackets
    if PAPER:
        lab = lambda v: f"{v:.1f}" if v >= 10 else f"{v:.2f}"
    else:
        lab = lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}"
    for col, (E, T, M) in enumerate(data_abs):
        ax, axm = axs[0, col], axs[1, col]
        names = list(T["1E"].keys())
        x = np.arange(len(names))
        for f, (setting, off) in enumerate((("1E", -0.19), ("2E", 0.19))):
            for i, nm in enumerate(names):
                c = COL["MLTS"] if nm.startswith("MLTS,") and "2 levels" not in nm else COL[nm]
                c = c if f == 0 else light(c)
                v = T[setting][nm]
                ax.bar(i + off, v, width=0.36, color=c, edgecolor="#1a1a1a", linewidth=0.4, hatch=hatch[f])
                ax.text(i + off, v * 1.10, lab(v), ha="center", va="bottom", fontsize=6)
                v = M[setting][nm]
                axm.bar(i + off, v, width=0.36, color=c, edgecolor="#1a1a1a", linewidth=0.4, hatch=hatch[f])
                axm.text(i + off, v + 0.02 * mmax, lab(v), ha="center", va="bottom", fontsize=5.5)
        tmin = min(v for s in T for v in T[s].values())
        tmax = max(v for s in T for v in T[s].values())
        ax.set_yscale("log"); ax.set_ylim(tmin / 2, tmax * 4)
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_title(E["title"].split(", ")[0], fontsize=7, pad=3)
        axm.set_ylim(0, mmax * 1.25)   # linear scale (user, 2026-11-13: no log scale on the memories)
        for a in (ax, axm):
            a.set_xticks(x); a.set_xlim(-0.6, len(names) - 0.4)
            a.grid(True, axis="y", which="major", ls=":", lw=0.4, color="#b4b4b0")
            a.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        ax.tick_params(labelbottom=False)
        axm.set_xticklabels([shown(nm).replace(", ", ",\n") for nm in names], fontsize=7)
        if col == 0:
            ax.set_ylabel("CPU time per coarse step [ms]" if PAPER else "wall time per coarse step (ms)", labelpad=2)
            axm.set_ylabel(("algebraic" if MEM else "peak") + (" memory [MiB]" if PAPER else " memory (MiB)"), labelpad=2)
    axs[0, 0].legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#5a5a56", ec="#1a1a1a", lw=0.4, label="setting (1E)"),
                              plt.Rectangle((0, 0), 1, 1, fc=light("#5a5a56"), ec="#1a1a1a", lw=0.4, hatch=hatch[1],
                                            label="setting (2E)")],
                     loc="upper left", fontsize=6.5, framealpha=0.92, edgecolor="#9a9a96", handlelength=1.2,
                     borderpad=0.35, labelspacing=0.25)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)


def bar_panels(ax, axm, G, Mem, runs, legend):
    """Gains over ERK4 (ax, log scale) and, measured, peak memory relative to ERK4 below them (axm), in the settings
    (1E) (dark) and (2E) (light); legend: the legend of the settings."""
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
    if runs:
        ax.tick_params(labelbottom=False)
    ax.set_ylabel("gain over ERK4" if runs else "model gain over ERK4", labelpad=1)
    top = max(max(G[s].values()) for s in G)
    ax.set_ylim(0.3, top * 4)
    ax.grid(True, axis="y", which="major", ls=":", lw=0.4, color="#b4b4b0")
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    if legend:
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#5a5a56", ec="#1a1a1a", lw=0.4, label="setting (1E)"),
                           plt.Rectangle((0, 0), 1, 1, fc=light("#5a5a56"), ec="#1a1a1a", lw=0.4, label="setting (2E)")],
                  loc="upper left", fontsize=6.5, frameon=False, handlelength=1.2)
    if not Mem:
        return
    # peak memory relative to ERK4 (measured runs only)
    for f, (setting, off) in enumerate((("1E", -0.19), ("2E", 0.19))):
        for i, nm in enumerate(names):
            c = COL["MLTS"] if nm.startswith("MLTS,") and "2 levels" not in nm else COL[nm]
            v = Mem[setting][nm]
            axm.bar(i + off, v, width=0.36, color=c if f == 0 else light(c), edgecolor="#1a1a1a", linewidth=0.4)
            axm.text(i + off, v + 0.04, f"{v:.2f}", ha="center", va="bottom", fontsize=5)
    axm.set_xlim(ax.get_xlim())
    axm.axhline(1, color="#6a6a66", lw=0.6, ls="--", zorder=0)
    axm.set_xticks(x); axm.set_xticklabels([nm.replace(", ", ",\n") for nm in names], fontsize=6.5)
    axm.set_ylim(0, 1.7); axm.set_yticks([0, 0.5, 1, 1.5])
    axm.set_ylabel("peak memory\n/ ERK4", labelpad=1, fontsize=7)
    axm.grid(True, axis="y", which="major", ls=":", lw=0.4, color="#b4b4b0")
    axm.tick_params(which="major", length=2.5, width=0.5, pad=1.5)


if __name__ == "__main__":
    main()
