#!/usr/bin/env python3
"""make_article.py -- figures (PDF) of Sections 5 and 6 of the LTS_HHO preprint (main.tex),
with one plotting style: the one of the a priori / convergence papers on Leray-Lions problems (Computer Modern fonts at
the printed size, panels 3.0 in wide for 0.48\\textwidth, full-width figures 6.3 in, dotted major grid, observed slope of
every segment, dashed reference slopes).

Color and marker code of the time schemes, the same in every figure:
    ERK4  black, circles        LTS   orange, squares      MLTS  dark red, triangles
    LI    blue, diamonds        SDIRK grey, inverted triangles (implicit reference)

Data: <results> as written by run_article.sh (one directory per run: c.dat, log.txt, time.txt, sensors, snapshots):
    alg/, form/ 6.6    rho/, stab/, stab2/ 6.2    conv/, tconv/, conv2/, tconv2/, meshes/ 6.4    ls/, ls2/, ls13/ 6.5
    cfl/, head/, head1/ 6.7 (earlier head model)
and the reduced models (4): reduced_models/reduced_mlts4_li.npz.
Usage: python3 make_article.py <results> <output directory>"""
import glob, math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "..", "..", "src"))
LS_MESHES = os.path.join(SRC, "mesh_generation", "lshape_graded", "meshes")
HEAD_MESHES = os.path.join(SRC, "mesh_generation", "brain_mri", "meshes")
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "article_figures")

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "font.size": 8, "axes.linewidth": 0.6, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42})

SCH = {   # scheme -> (label, color, marker)
    "erk":   ("ERK4", "#1a1a1a", "o"),
    "lts":   ("LTS", "#e6750d", "s"),
    "mlts":  ("MLTS", plt.cm.Reds(0.85), "^"),
    "li":    ("LI", plt.cm.Blues(0.8), "D"),
    "sdirk": ("SDIRK", "#8c8c88", "v"),
}
GRID = dict(ls=":", lw=0.4, color="#b4b4b0")
BOX = dict(boxstyle="round,pad=0.35", fc="white", ec="#9a9a96", lw=0.5)
MISSING = []


# ----------------------------------------------------------------------------------------------------------- data
def dat(d):
    """(column names, rows) of <R>/<d>/c.dat, or None."""
    f = os.path.join(R, d, "c.dat")
    if not os.path.exists(f):
        return None
    names = open(f).readline().lstrip("#").split()
    rows = np.atleast_2d(np.loadtxt(f, comments="#"))
    return names, rows


def col(d, name):
    D = dat(d)
    if D is None or name not in D[0]:
        return None
    return D[1][:, D[0].index(name)]


def log(d):
    f = os.path.join(R, d, "log.txt")
    return re.sub(r"\x1b\[[0-9;]*m", "", open(f).read()) if os.path.exists(f) else ""


def wall(d):
    f = os.path.join(R, d, "time.txt")
    if not os.path.exists(f):
        return None, None
    m = re.match(r"([0-9.]+) s ([0-9]+) KB", open(f).read().strip())
    return (float(m.group(1)), float(m.group(2)) / 1024) if m else (None, None)


def n_cells(d):
    m = re.search(r"matrix-free operator: (\d+) cells", log(d))
    return int(m.group(1)) if m else None


def mlts_levels(d):
    m = re.search(r"MLTS: (\d+) levels, cells per level \(coarsest first\):([0-9 ]+)", log(d))
    return (int(m.group(1)), [int(x) for x in m.group(2).split()]) if m else (None, None)


def li_system(d):
    """Unknowns of the factorized fine system of LI (the fine faces when condensed) and memory of the factors (MB)."""
    m = re.search(r"LI fine system: (\d+) unknowns.*?\(~([0-9.e+]+) MB\)", log(d))
    if not m:
        return None, None
    c = re.search(r"condensed on the faces: (\d+) unknowns", log(d))
    return (int(c.group(1)) if c else int(m.group(1))), float(m.group(2))


def rates(h, e):
    return np.log(e[:-1] / e[1:]) / np.log(h[:-1] / h[1:])


def fmt(x, f="{:.2e}"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    return f.format(x)


def tex_sci(x, digits=2):
    if x is None or not np.isfinite(x):
        return "--"
    m, e = f"{x:.{digits}e}".split("e")
    return rf"${m}\,{{\times}}10^{{{int(e)}}}$"


def write(name, s):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(s)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("ok", name)


def style_log(ax, xlabel, ylabel):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="major", **GRID)
    ax.set_xlabel(xlabel, labelpad=1); ax.set_ylabel(ylabel, labelpad=1)
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    ax.tick_params(which="minor", length=1.5, width=0.4)


def mem_ticks(axis):
    """Log axis of memories (MiB): ticks 1, 2, 5 per decade, plain numbers."""
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


def legend(ax, **kw):
    kw.setdefault("loc", "lower right")
    lg = ax.legend(framealpha=0.92, edgecolor="#9a9a96", handlelength=1.8, borderpad=0.35, labelspacing=0.25,
                   fancybox=False, **kw)
    lg.get_frame().set_linewidth(0.5)
    return lg


def slope_labels(fig, ax, curves, theory, name, labelled=None):
    """Observed slope of every segment, placed greedily from the finest segment on, never on a label, a solid curve,
    the legend or a text box, avoiding the dashed reference lines when possible (make_figures.py of the a priori
    paper)."""
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    box = ax.get_window_extent(rnd)
    placed = [t.get_window_extent(rnd) for t in ax.texts]
    if ax.get_legend():
        placed.append(ax.get_legend().get_window_extent(rnd))

    def sample(points, n):
        p = ax.transData.transform(points)
        return [a + t * (b - a) for a, b in zip(p[:-1], p[1:]) for t in np.linspace(0., 1., n)]
    hard = np.array([q for h, e, _ in curves for q in sample(np.column_stack([h, e]), 30)])
    soft = np.array([q for hl, el in theory for q in sample(np.column_stack(
        [np.geomspace(hl[0], hl[1], 40), np.geomspace(el[0], el[1], 40)]), 4)]) if theory else np.zeros((0, 2))

    def hits(pts, bb):
        return len(pts) and np.any((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1))
    offsets = [(-6, 4), (6, -9), (-11, -2), (10, 2), (0, 7), (0, -11), (-15, 6), (14, -11)]
    nseg = max(len(h) - 1 for h, _, _ in curves)
    for back in range(1, nseg + 1):
        for k, (h, e, c) in enumerate(curves):
            if labelled is not None and k not in labelled:
                continue
            i = len(h) - 1 - back
            if i < 0:
                continue
            r = rates(h, e)[i]
            xm, ym = math.sqrt(h[i] * h[i + 1]), math.sqrt(e[i] * e[i + 1])
            done = False
            for strict in (True, False):
                for off in offsets:
                    t = ax.annotate(f"{r:.2f}", xy=(xm, ym), xytext=off, textcoords="offset points", fontsize=5.5,
                                    color=c, ha="center", va="center", annotation_clip=False)
                    bb = t.get_window_extent(rnd).expanded(1.08, 1.15)
                    ok = (box.x0 <= bb.x0 and bb.x1 <= box.x1 and box.y0 <= bb.y0 and bb.y1 <= box.y1
                          and not hits(hard, bb) and not (strict and hits(soft, bb))
                          and not any(bb.overlaps(o) for o in placed))
                    if ok:
                        placed.append(bb); done = True
                        break
                    t.remove()
                if done:
                    break
            if not done and back == 1:
                MISSING.append(f"{name}: finest segment of a curve colored {c}")


def conv_panel(ax, fig, series, order, xlabel, ylabel, formula, name, legend_loc="lower right", ref_index=0,
               labelled=None):
    """Convergence panel: series = [(h, e, scheme key or (label, color, marker), linestyle)], dashed line of slope
    `order` through the first point of series[ref_index], observed slopes on every segment."""
    H = np.concatenate([s[0] for s in series]); E = np.concatenate([s[1] for s in series])
    h_lo, h_hi = H.min() * 0.6, H.max() * 1.6
    ax.set_xlim(h_lo, h_hi); ax.set_ylim(E.min() * 0.2, E.max() * 8.)
    style_log(ax, xlabel, ylabel)
    curves, theory = [], []
    for h, e, key, ls in series:
        lab, c, mk = SCH[key] if isinstance(key, str) else key
        ax.plot(h, e, ls + mk, color=c, lw=1.0, ms=3.0, mfc=c if ls == "-" else "white", zorder=3, label=lab)
        curves.append((h, e, c))
    if order:
        h0, e0 = series[ref_index][0][0], series[ref_index][1][0]
        hl = np.array([h_lo, h0])
        ax.plot(hl, e0 * 0.5 * (hl / h0) ** order, ls="--", lw=0.8, color="#6a6a66", zorder=1)
        theory.append((hl, e0 * 0.5 * (hl / h0) ** order))
    if formula:
        ax.text(0.03, 0.96, formula, transform=ax.transAxes, fontsize=6.5, va="top", ha="left", bbox=BOX)
    legend(ax, loc=legend_loc)
    slope_labels(fig, ax, curves, theory, name, labelled)


def read_cells(path):
    """Cells of a c_cells.dat file: (medium, vertices)."""
    cells = []
    for l in open(path):
        if l.startswith("#"):
            continue
        v = l.split(); n = int(v[1])
        cells.append((int(v[0]), np.array(v[2:2 + 2 * n], float).reshape(n, 2)))
    return cells


def read_polygons(path):
    """Points, cells (vertex lists) and tags of a "polygons" mesh file."""
    L = [l for l in open(path) if not l.startswith("#")]
    npt, nc = map(int, L[0].split())
    P = np.array([list(map(float, l.split())) for l in L[1:1 + npt]])
    cells, tags = [], []
    for l in L[1 + npt:1 + npt + nc]:
        a = list(map(int, l.split()))
        cells.append(P[a[1:1 + a[0]]]); tags.append(a[-1])
    return P, cells, np.array(tags)


def edges(polys):
    return [[p[i], p[(i + 1) % len(p)]] for p in polys for i in range(len(p))]


# ------------------------------------------------------------------------------------ 6.6 the two formulations
def formulations():
    # top: times; bottom: peak memory (maximum resident set size, time.txt) of the same runs
    fig, A = plt.subplots(2, 2, figsize=(6.3, 4.5), gridspec_kw=dict(width_ratios=[1.15, 1]))
    axs, axm = A[0], A[1]
    ax = axs[0]
    for o, lab, c, mk in (("o1", r"setting (1E), $k=3$", "#1a1a1a", "o"), ("o2", r"setting (2E), $k=2$", "#9a9a96", "s")):
        nu, ms, mem = [], [], []
        for l in range(4):
            D = dat(f"alg/{o}_l{l}")
            if D is None:
                continue
            names, rows = D
            n_steps = round(1.0 / rows[-1, names.index("dt")])
            # integration time per step: difference of the runs with n and n + 80000 / 2^l steps (the setup cancels;
            # the second-order setup estimates omega_max and the splitting radius), best of the repetitions (_r2 to
            # _r6, those present); without the long runs, wall time per step
            extra = 80000 // 2 ** l
            per = []
            for rep in ("",) + tuple(f"_r{i}" for i in range(2, 7)):
                w, _ = wall(f"alg/{o}_l{l}{rep}")
                wl, _ = wall(f"alg/{o}_l{l}_long{rep}")
                if w is not None and wl is not None:
                    per.append((wl - w) / extra)
            if not per:
                w, _ = wall(f"alg/{o}_l{l}")
                per = [w / n_steps]
            nu.append(8 * 2 ** l); ms.append(1e3 * min(per)); mem.append(wall(f"alg/{o}_l{l}")[1] or np.nan)
        if nu:
            ax.plot(np.array(nu) ** 2, ms, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
            axm[0].plot(np.array(nu) ** 2, mem, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    style_log(ax, "cells", "wall time per ERK4 step (ms)")
    legend(ax, loc="upper left")
    style_log(axm[0], "cells", "peak memory of ERK4 (MiB)"); mem_ticks(axm[0].yaxis)
    legend(axm[0], loc="upper left")
    # right: wall time of the whole run on one thread (top) and its peak memory (bottom), the two formulations side by
    # side for each scheme (1E: full colour, 2E: light); second order: MLTS and LI with the projected condensation
    groups = (("erk", "erk", "erk"), ("mlts", "mlts", "mltsp"), ("li", "li", "lip"))     # (style, 1E run, 2E run)
    light = lambda c: tuple(0.55 + 0.45 * x for x in matplotlib.colors.to_rgb(c))
    ax = axs[1]
    for k, (st, r1, r2) in enumerate(groups):
        color = SCH[st][1]
        for f, (o, run) in enumerate((("o1", r1), ("o2", r2))):
            # best of the repetitions form/, form_r2/, form_r3/ (same jobs, run alone one after the other)
            ws = [wall(f"{rep}/{o}_{run}_q1") for rep in ("form", "form_r2", "form_r3")]
            ws = [w for w in ws if w[0] is not None]
            v, mv = min(ws) if ws else (np.nan, np.nan)      # memory of the run with the best time
            xi = k + (-0.18 if f == 0 else 0.18)
            for a, y in ((ax, v), (axm[1], mv)):
                a.bar(xi, y, color=color if f == 0 else light(color), width=0.34, edgecolor="#1a1a1a", linewidth=0.4)
            if np.isfinite(v):
                ax.text(xi, v + 1.5, f"{v:.0f}" if v >= 10 else f"{v:.1f}", ha="center", va="bottom", fontsize=6)
            if np.isfinite(mv):
                axm[1].text(xi, mv + 6, f"{mv:.0f}", ha="center", va="bottom", fontsize=6)
    xs, labels = list(range(len(groups))), [SCH[g[0]][0] for g in groups]
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=7)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, fc="#5a5a56", ec="#1a1a1a", lw=0.4, label="setting (1E), $k=3$"),
                       Rectangle((0, 0), 1, 1, fc=light("#5a5a56"), ec="#1a1a1a", lw=0.4, label="setting (2E), $k=2$")],
              loc="upper right", fontsize=6.5, frameon=False, handlelength=1.4)
    ax.set_ylabel(r"wall time on one thread (s)", labelpad=1); ax.grid(True, axis="y", which="major", **GRID)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    ax.set_ylim(0, max(b.get_height() for b in ax.patches if np.isfinite(b.get_height())) * 1.22)
    ax.set_title(r"$64\times64$ cells, central patch refined 4 times ($p=16$)", fontsize=7, pad=3)
    ax = axm[1]
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("peak memory on one thread (MiB)", labelpad=1); ax.grid(True, axis="y", which="major", **GRID)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    ax.set_ylim(0, max(b.get_height() for b in ax.patches if np.isfinite(b.get_height())) * 1.22)
    fig.tight_layout(pad=0.2, w_pad=1.2, h_pad=1.0)
    save(fig, "form.pdf")


# --------------------------------------------------------------------------------------------- 5 reduced models
# reduced models (reduced_models/reduced_mlts4_li.py): chains of four tiers, ratio p between consecutive levels, total
# ratio p^3; ERK4 with the global step, MLTS with L = 4 levels, LI (fine region: the levels >= 1)
RED_LABEL = {"mlts4": r"MLTS, $L=4$"}
RED_STYLE = {"erk": "erk", "mlts4": "mlts", "li": "li"}


def reduced():
    f = os.path.join(HERE, "reduced_models", "reduced_mlts4_li.npz")
    if not os.path.exists(f):
        print("reduced models: missing npz"); return
    D = np.load(f)
    X = D["p"].astype(float) ** 3
    models = (("order2", "model 1: second order"), ("order1_nobuf", "model 2: first order, no overlap"),
              ("order1_buf", "model 3: first order, overlap"))
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), sharey=True)
    for ax, (key, title) in zip(axs, models):
        for s in ("erk", "mlts4", "li"):
            lab, c, mk = SCH[RED_STYLE[s]]
            ax.plot(X, D[f"{key}_{s}_rho"].max(axis=0), "-" + mk, color=c, lw=1.0, ms=2.5, label=RED_LABEL.get(s, lab))
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
        ax.set_xlabel(r"$p^3$", labelpad=1); ax.set_title(title, fontsize=7, pad=3)
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        ax.axhline(1.0, color="#6a6a66", lw=0.6, ls="--")
    axs[0].set_ylabel(r"$\max_\varepsilon\,\rho(C_p(\Delta t_{\rm ref}))$", labelpad=1)
    legend(axs[0], loc="upper left")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "red_rho.pdf")
    # critical step (largest step below which every step is stable), relative to that of ERK4 at p = 1, smallest over eps
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), sharey=True)
    for ax, (key, title) in zip(axs, models):
        for s in ("erk", "mlts4", "li"):
            lab, c, mk = SCH[RED_STYLE[s]]
            ax.plot(X, D[f"{key}_{s}_dtc"].min(axis=0), "-" + mk, color=c, lw=1.0, ms=2.5, label=RED_LABEL.get(s, lab))
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
        ax.set_xlabel(r"$p^3$", labelpad=1); ax.set_title(title, fontsize=7, pad=3)
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    axs[0].set_ylabel(r"$\min_\varepsilon\,\Delta t_c(p)/\Delta t_{\rm ref}$", labelpad=1)
    legend(axs[0], loc="lower left")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "red_dtc.pdf")


# ------------------------------------------------------------------------------------------------ 6.2 stability
def stability():
    L = [0, 1, 2, 3, 4, 5]
    LAB = {"lts": "MLTS, 2 levels", "mlts": "MLTS, ratios 2"}   # no "LTS": the two-level scheme is MLTS with L = 2
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.35))
    ax = axs[0]
    for s in ("erk", "lts", "mlts", "li"):
        pp, rr = [], []
        for l in L:
            r = col(f"rho/L{l}_{s}", "rho")
            if r is not None and r[0] > 0:
                pp.append(2 ** l); rr.append(r[0])
        lab, c, mk = SCH[s]; lab = LAB.get(s, lab)
        ax.plot(pp, rr, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
    ax.axhline(1.0, color="#6a6a66", lw=0.6, ls="--")
    ax.set_xlabel(r"$p=h_{\max}/h_{\min}$", labelpad=1); ax.set_ylabel(r"$\rho(\mathcal{C})$, $\Delta t=6.4\times10^{-3}$", labelpad=1)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    legend(ax, loc="upper left")
    ax = axs[1]
    for s in ("erk", "lts", "mlts", "li"):
        pp, dd = [], []
        for l in range(6):
            d = f"stab/L{l}_{'lis' if s == 'li' and l > 0 else s}"
            v = col(d, "dt")
            if v is None or (s == "li" and l == 0):
                continue
            pp.append(2 ** l); dd.append(v[0])
        lab, c, mk = SCH[s]; lab = LAB.get(s, lab)
        ax.plot(pp, dd, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    # second-order formulation (setting (2E), stab2/): ERK4, LTS, MLTS and LI with the projected condensation (dashed,
    # hollow); LI from p = 2 (on the uniform mesh, its fine region is the whole mesh)
    second = False
    for s, key in (("erk", "erk"), ("lts", "ltsp"), ("mlts", "mltsp"), ("li", "lip")):
        pp, dd = [], []
        for l in range(6):
            v = col(f"stab2/L{l}_{key}", "dt")
            if v is None or (s == "li" and l == 0):
                continue
            pp.append(2 ** l); dd.append(v[0])
        if pp:
            lab, c, mk = SCH[s]
            ax.plot(pp, dd, "--" + mk, color=c, lw=0.9, ms=3.2, mfc="white", mew=0.8)
            second = True
    if second:
        ax.legend(handles=[Line2D([], [], color="#1a1a1a", lw=1.0, marker="o", ms=3.2, label="first order (1E)"),
                           Line2D([], [], color="#1a1a1a", lw=0.9, ls="--", marker="o", ms=3.2, mfc="white",
                                  label="second order (2E)")],
                  loc="lower left", fontsize=6, frameon=False)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
    ax.set_xlabel(r"$p=h_{\max}/h_{\min}$", labelpad=1); ax.set_ylabel(r"CFL limit $\Delta t_{\rm crit}$", labelpad=1)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    ax.set_ylim(1.5e-4, 2.5e-2)
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "stab.pdf")


# -------------------------------------------------------------------------------------- 6.4 smooth convergence
def meshes_cart():
    """One row: the uniform meshes 8 x 8 and 16 x 16 (ERK4), then the same with the central box refined three times
    (MLTS and LI)."""
    names = [("u", 0), ("u", 1), ("g", 0), ("g", 1)]
    files = [os.path.join(R, "meshes", f"cart_{k}{l}", "c_cells.dat") for k, l in names]
    if not all(os.path.exists(f) for f in files):
        print("meshes of 6.4: missing c_cells.dat"); return
    fig, axes = plt.subplots(1, 4, figsize=(6.3, 1.75))
    for ax, f, (k, l) in zip(axes, files, names):
        cells = read_cells(f)
        lw = 0.5 if len(cells) <= 600 else 0.25
        ax.add_collection(LineCollection(edges([c for _, c in cells]), colors="#1a1a1a", linewidths=lw))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        n = 8 * 2 ** l
        ax.set_title(f"${n}\\times{n}$" + (", refined" if k == "g" else "") + f", {len(cells)} cells", fontsize=7.5,
                     pad=3)
    fig.tight_layout(pad=0.2, w_pad=0.8)
    save(fig, "cart_meshes.pdf")


def convergence():
    """5.2.1, smooth solution, both settings in one figure, one panel per setting with both errors: left (1E), k = 3
    (conv/: v full markers, p hollow ones); right (2E), k = 2 (conv2/, MLTS and LI with the projected condensation: p
    full, grad p hollow). ERK4 on the uniform meshes (dashed), MLTS and LI on the meshes refined three times; dashed grey: reference
    slopes; observed slopes of MLTS on every segment (LI coincides with MLTS)."""
    P = (("conv", {"erk": "conv/erk_u", "mlts": "conv/mlts_g", "li": "conv/li_g"},
          (("acoustic_q", r"$v$", 4), ("acoustic_p", r"$p$", 4)), r"(1E), $k=3$, $\Delta t\propto h_{\max}$"),
         ("conv2", {"erk": "conv2/erk_u", "mlts": "conv2/mltsp_g", "li": "conv2/lip_g"},
          (("acoustic_phi", r"$p$", 4), ("acoustic_grad_phi", r"$\nabla p$", 3)), r"(2E), $k=2$, $\Delta t\propto h_{\max}$"))
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.7))
    for ax, (_, S, Q, tag) in zip(axs, P):
        D = {k: dat(v) for k, v in S.items()}
        if any(v is None for v in D.values()):
            print("convergence: missing runs"); return
        g = lambda k, n: D[k][1][:, D[k][0].index(n)]
        H = g("erk", "h"); E = np.concatenate([g(k, f) for k in S for f, _, _ in Q])
        h_lo, h_hi = H.min() * 0.6, H.max() * 2.4
        ax.set_xlim(h_lo, h_hi); ax.set_ylim(E.min() * 0.15, E.max() * 8.)
        style_log(ax, r"$h_{\max}$", r"error in the $L^2(\Omega)$-norm")
        curves, theory, labelled, drawn = [], [], [], set()
        upper = max(range(len(Q)), key=lambda j: g("erk", Q[j][0])[0])
        for j, (fld, qlab, order) in enumerate(Q):
            for key, mk, lab in (("erk", "o", "ERK4, uniform"), ("mlts", "^", "MLTS, refined"), ("li", "D", "LI, refined")):
                c = SCH[key][1]
                ls = "--" if key == "erk" else "-"     # first error full markers, second one hollow
                ax.plot(g(key, "h"), g(key, fld), ls + mk, color=c, lw=1.0, ms=3.0 if key != "li" else 2.6, mew=0.8,
                        mfc=c if j == 0 else "white",
                        zorder=3, label=lab if j == 0 else None)
                if key == "mlts":
                    labelled.append(len(curves))
                curves.append((g(key, "h"), g(key, fld), c))
            # dashed reference slope below the lowest bundle of each order; the name of the error at the coarse end,
            # above the upper bundle and below the lower one
            h0, e0 = g("erk", "h")[0], min(g(k, fld)[0] for k in S)
            e_top = max(g(k, fld)[0] for k in S)
            if order not in drawn or j != upper:
                hl = np.array([h_lo, h0])
                ref = min(g(k, f)[0] for k in S for f, _, o in Q if o == order)
                ax.plot(hl, ref * 0.45 * (hl / h0) ** order, ls="--", lw=0.8, color="#6a6a66", zorder=1)
                theory.append((hl, ref * 0.45 * (hl / h0) ** order)); drawn.add(order)
            # upper label ending at the coarse end (stays inside the frame), lower one starting after it
            ax.text(h0 * 1.05 if j == upper else h0 * 1.12, e_top * 1.6 if j == upper else e0 / 1.5,
                    qlab + (" (full)" if j == 0 else " (hollow)"), fontsize=7, va="bottom" if j == upper else "top",
                    ha="right" if j == upper else "left")
        ax.text(0.03, 0.96, tag, transform=ax.transAxes, fontsize=6.5, va="top", ha="left", bbox=BOX)
        legend(ax, loc="lower right")
        slope_labels(fig, ax, curves, theory, "conv_cart.pdf", labelled)
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "conv_cart.pdf")


def time_orders():
    """5.2.1, time convergence on the 16 x 16 mesh refined three times, both settings in one figure with the format of
    conv_cart.pdf, for the solution p = sin(w t)(x(1-x) + y(1-y)), w = 8 sqrt2 pi (problem.solution =
    "acoustic_quadratic_sum", problem.frequency), of degree two in space and reproduced exactly by the HHO spaces: the
    error on p against it at t_f is the time error alone (no reference solution). Runs of lts_work/tconv_poly merged
    by its merge.py into tconvp/ (6 levels and their extensions to smaller steps): left (1E), k = 3; right (2E), k = 2
    (MLTS and LI with the projected condensation, globally implicit SDIRK(3,4) for comparison)."""
    P = ((r"(1E), $k=3$, $16\times16$, $p=8$", "acoustic_p", (("mlts", "tconvp/mlts"), ("li", "tconvp/li"))),
         (r"(2E), $k=2$, $16\times16$, $p=8$", "acoustic_phi",
          (("mlts", "tconvp/2_mltsp"), ("li", "tconvp/2_lip"), ("sdirk", "tconvp/2_sdirk"))))
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.7))
    for ax, (tag, fld, runs) in zip(axs, P):
        series = []
        for key, d in runs:
            D = dat(d)
            if D is None:
                continue
            dt, e = D[1][:, D[0].index("dt")], D[1][:, D[0].index(fld)]
            ok = np.isfinite(e) & (e < 1)
            series.append((dt[ok], e[ok], key, "-"))
        if not series:
            print("time_orders: missing runs"); return
        conv_panel(ax, fig, series, 4, r"$\Delta t$", r"$\Vert p(t_f)-p_h(t_f)\Vert_{L^2(\Omega)}$", tag, "tconv.pdf",
                   legend_loc="lower right")
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "tconv.pdf")


def convergence_2e():
    """Merged into convergence() (conv_cart.pdf, right panel) on 2026-11-09."""
    convergence()


def time_orders_2e():
    """Merged into time_orders() (tconv.pdf, right panel) on 2026-11-09."""
    time_orders()


# ----------------------------------------------------------------------------------------------- 6.5 L-shape
def lshape_meshes():
    f = os.path.join(LS_MESHES, "lshape_mu1-5_n8.txt")
    if not os.path.exists(f):
        print("L-shape mesh missing"); return
    _, cells, _ = read_polygons(f)
    smin = min(np.ptp(c[:, 0]) for c in cells)
    fig, axes = plt.subplots(1, 3, figsize=(5.1, 1.7))
    fig.patch.set_alpha(0.0)
    for ax, w in zip(axes, (0.5, 0.5 / 32, 0.5 / 2 ** round(math.log2(0.5 / (6 * smin))))):
        sel = [c for c in cells if np.min(np.abs(c).max(axis=1)) <= 1.5 * w]
        ax.add_collection(LineCollection(edges(sel), colors="#1a1a1a", linewidths=0.45 if w == 0.5 else 0.4))
        ax.set_xlim(-w, w); ax.set_ylim(-w, w); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        # the removed quadrant (0, w)^2 stays empty and transparent: no frame, the L-shaped outline of the window instead
        ax.patch.set_alpha(0.0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.plot([-w, w, w, 0, 0, -w, -w], [-w, -w, 0, 0, w, w, -w], color="#1a1a1a", lw=0.7, clip_on=False,
                solid_joinstyle="miter")
        ax.set_title(f"{len(cells)} cells" if w == 0.5 else rf"zoom $\times{0.5 / w:.0f}$", fontsize=7.5, pad=3)
    axes[0].text(0.03, 0.03, rf"$h_{{\max}}/h_{{\min}}={0.0625 / smin:.0f}$", transform=axes[0].transAxes, fontsize=6, bbox=BOX)
    fig.tight_layout(pad=0.2, w_pad=3.2)
    save(fig, "lshape_meshes.pdf")


def ls_series(prefix, levels):
    h, p, q, cpu, rows = [], [], [], [], []
    for l in levels:
        d = f"ls/{prefix}_l{l}"
        D = dat(d)
        if D is None:
            continue
        c = lambda n: D[1][-1, D[0].index(n)]
        if not (np.isfinite(c("acoustic_q")) and c("acoustic_q") < 1):
            continue
        h.append(c("h")); p.append(c("acoustic_p")); q.append(c("acoustic_q")); cpu.append(c("cpu"))
        rows.append(d)
    return np.array(h), np.array(p), np.array(q), np.array(cpu), rows


def two_errors_panel(ax, fig, Q, name, erk=None, legend_loc="lower right"):
    """Errors of two quantities against the mesh size in one panel, as conv_cart.pdf: Q = [(label, order, {scheme:
    (h, e)}), ...], first quantity full markers, second hollow; dashed reference slope below each order, the name of the
    error at the coarse end; erk = [(h, e), (h, e)]: ERK4 points of the two quantities; observed slopes of MLTS."""
    H = np.concatenate([d[0] for _, _, D in Q for d in D.values()]); E = np.concatenate([d[1] for _, _, D in Q for d in D.values()])
    h_lo, h_hi = H.min() * 0.6, H.max() * 2.6
    ax.set_xlim(h_lo, h_hi); ax.set_ylim(E.min() * 0.15, E.max() * 8.)
    style_log(ax, r"$h_{\max}$", r"error in the $L^2(\Omega)$-norm")
    curves, theory, labelled, drawn = [], [], [], set()
    upper = max(range(len(Q)), key=lambda j: max(d[1][0] for d in Q[j][2].values()))
    for j, (qlab, order, D) in enumerate(Q):
        for key, (h, e) in D.items():
            lab, c, mk = SCH[key]
            ax.plot(h, e, "-" + mk, color=c, lw=1.0, ms=3.0 if key != "li" else 2.6, mew=0.8, mfc=c if j == 0 else "white",
                    zorder=3, label=lab if j == 0 else None)
            if key == "mlts":
                labelled.append(len(curves))
            curves.append((h, e, c))
        if erk is not None and len(erk[j][0]):
            ax.plot(erk[j][0], erk[j][1], "o", color=SCH["erk"][1], ms=3.5, mew=0.8, mfc=SCH["erk"][1] if j == 0 else "white",
                    zorder=4, label="ERK4" if j == 0 else None)
        h0 = max(d[0][0] for d in D.values()); e0 = min(d[1][0] for d in D.values()); e_top = max(d[1][0] for d in D.values())
        if order not in drawn or j != upper:
            ref = min(d[1][0] for _, o, DD in Q if o == order for d in DD.values())
            hl = np.array([h_lo, h0])
            ax.plot(hl, ref * 0.45 * (hl / h0) ** order, ls="--", lw=0.8, color="#6a6a66", zorder=1)
            theory.append((hl, ref * 0.45 * (hl / h0) ** order)); drawn.add(order)
        ax.text(h0 * 1.05 if j == upper else h0 * 1.12, e_top * 1.6 if j == upper else e0 / 1.5,
                qlab + (" (full)" if j == 0 else " (hollow)"), fontsize=7, va="bottom" if j == upper else "top",
                ha="right" if j == upper else "left")
    legend(ax, loc=legend_loc)
    slope_labels(fig, ax, curves, theory, name, labelled)


def lshape():
    lv = (1, 2, 3, 4, 5)
    # LI at half the step (dt = 1/(40n)): the error constant of SDIRK(3,4) is about 20 times that of RK4 (ls/li2_*)
    data = {"mlts": ls_series("mlts", lv), "li": ls_series("li2", lv)}
    erk = ls_series("erk", (1,))
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.35))
    Q = [(r"$v$", 4, {s: (data[s][0], data[s][2]) for s in ("mlts", "li") if len(data[s][0]) >= 2}),
         (r"$p$", 4, {s: (data[s][0], data[s][1]) for s in ("mlts", "li") if len(data[s][0]) >= 2})]
    two_errors_panel(axs[0], fig, Q, "lshape[errors]", erk=[(erk[0], erk[2]), (erk[0], erk[1])])
    ax = axs[1]
    for s in ("mlts", "li"):
        hh, pp, qq, cc, _ = data[s]
        lab, c, mk = SCH[s]
        ax.plot(cc, qq, "-" + mk, color=c, lw=1.0, ms=3.0, label=lab)
    if len(erk[0]):
        ax.plot(erk[3], erk[2], "o", color=SCH["erk"][1], ms=3.5, label="ERK4")
    style_log(ax, "CPU time (s)", r"$\Vert v-v_h\Vert_{L^2(\Omega)}$")
    legend(ax, loc="lower left")
    ax = axs[2]
    for s in ("mlts", "li"):
        lab, c, mk = SCH[s]
        ax.plot([wall(d)[1] for d in data[s][4]], data[s][2], "-" + mk, color=c, lw=1.0, ms=3.0, label=lab)
    if len(erk[0]):
        ax.plot([wall(erk[4][0])[1]], erk[2], "o", color=SCH["erk"][1], ms=3.5, label="ERK4")
    style_log(ax, "peak memory (MiB)", r"$\Vert v-v_h\Vert_{L^2(\Omega)}$"); mem_ticks(ax.xaxis)
    legend(ax, loc="lower left")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "lshape.pdf")


def ls2_series(prefix, levels):
    """(h, p, grad p, time) of ls2/<prefix>_l<l>: time = cpu column (construction and integration, as the setting (1E)),
    wall time of the run alone for the runs without it."""
    h, p, g, t, m = [], [], [], [], []
    for l in levels:
        d = f"ls2/{prefix}_l{l}"
        D = dat(d)
        if D is None:
            continue
        c = lambda n: D[1][-1, D[0].index(n)]
        if not (np.isfinite(c("acoustic_grad_phi")) and c("acoustic_grad_phi") < 1):
            continue
        h.append(c("h")); p.append(c("acoustic_phi")); g.append(c("acoustic_grad_phi"))
        t.append(c("cpu") if "cpu" in D[0] else (wall(d)[0] or np.nan)); m.append(wall(d)[1] or np.nan)
    return np.array(h), np.array(p), np.array(g), np.array(t), np.array(m)


def lshape_2e():
    """5.2.2, setting (2E), k = 2: MLTS at dt = 1/(25n) with the cubic forcing (ls2/mltsp_*), LI at dt = 1/(100n)
    (ls2/li4_*), both with the projected condensation (two rings)."""
    lv = (1, 2, 3, 4, 5)
    data = {"mlts": ls2_series("mltsp", lv), "li": ls2_series("li4", lv)}
    if any(len(v[0]) < 2 for v in data.values()):
        print("lshape_2e: missing runs"); return
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.35))
    Q = [(r"$p$", 4, {s: (data[s][0], data[s][1]) for s in ("mlts", "li")}),
         (r"$\nabla p$", 3, {s: (data[s][0], data[s][2]) for s in ("mlts", "li")})]
    two_errors_panel(axs[0], fig, Q, "lshape2[errors]")
    ax = axs[1]
    for s in ("mlts", "li"):
        lab, c, mk = SCH[s]
        ax.plot(data[s][3], data[s][2], "-" + mk, color=c, lw=1.0, ms=3.0, label=lab)
    style_log(ax, "CPU time (s)", r"$\Vert \nabla p-\nabla p_h\Vert_{L^2(\Omega)}$")
    legend(ax, loc="upper right")
    ax = axs[2]
    for s in ("mlts", "li"):
        lab, c, mk = SCH[s]
        ax.plot(data[s][4], data[s][2], "-" + mk, color=c, lw=1.0, ms=3.0, label=lab)
    style_log(ax, "peak memory (MiB)", r"$\Vert \nabla p-\nabla p_h\Vert_{L^2(\Omega)}$"); mem_ticks(ax.xaxis)
    legend(ax, loc="upper right")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "lshape2.pdf")


def lshape_formulations():
    """5.3, comparison of the formulations on the L-shaped domain (the runs of Figs. ls and ls2): MLTS and LI in the
    setting (1E) (k = 3, full lines) and (2E) (k = 2, dashed, hollow markers); CPU time (construction and integration, six
    threads) and peak memory as functions of the mesh size, and error on grad u (velocity in (1E), gradient of p = u in
    (2E)) as a function of both."""
    lv = (1, 2, 3, 4, 5)
    one = {"mlts": ls_series("mlts", lv), "li": ls_series("li2", lv)}
    two = {"mlts": ls2_series("mltsp", lv), "li": ls2_series("li4", lv)}
    if any(len(v[0]) < 2 for v in list(one.values()) + list(two.values())):
        print("lshape_formulations: missing runs"); return
    # (h, error on grad u, time, memory) of each series
    S = {}
    for s in ("mlts", "li"):
        h, _, q, cpu, rows = one[s]
        S[(s, 1)] = (h, q, cpu, np.array([wall(d)[1] or np.nan for d in rows]))
        h, _, g, t, m = two[s]
        S[(s, 2)] = (h, g, t, m)
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5))
    for (s, o), (h, e, t, m) in S.items():
        lab, c, mk = SCH[s]
        kw = dict(color=c, lw=1.1 if o == 1 else 1.0, ms=3.4, label=f"{lab}, setting ({o}E)")
        if o == 2:
            kw.update(mfc="white", mew=0.8)
        ls = "-" if o == 1 else "--"
        axs[0].plot(t, e, ls + mk, **kw)
        axs[1].plot(m, e, ls + mk, **kw)
    style_log(axs[0], "CPU time on six threads (s)", r"error on $\nabla u$")
    style_log(axs[1], "peak memory (MiB)", r"error on $\nabla u$"); mem_ticks(axs[1].xaxis)
    fig.legend(*axs[0].get_legend_handles_labels(), loc="upper center", ncol=4, fontsize=6.5, frameon=False,
               bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.6)
    fig.tight_layout(pad=0.2, w_pad=1.5, rect=(0, 0, 1, 0.9))
    save(fig, "lshape_form.pdf")


# --------------------------------------------------------------------------------------------- 6.7 head model
TISSUES = [(3, "water", "#cfe3f1"), (4, "scalp", "#e9b9a8"), (1, "cortical bone", "#f2eee2"),
           (2, "cancellous bone", "#c9b98f"), (5, "CSF", "#7fb2d9"), (6, "gray matter", "#9c9c9c"),
           (7, "white matter", "#efefef"), (8, "titanium", "#2b2b2b")]


def head_model():
    tf = os.path.join(HEAD_MESHES, "head_s80_plate_tissues.npy")
    mf = os.path.join(HEAD_MESHES, "head_s80_plate.txt")
    if not (os.path.exists(tf) and os.path.exists(mf)):
        print("head model files missing"); return
    T = np.load(tf)
    P, cells, tags = read_polygons(mf)
    col_of = {t: c for t, _, c in TISSUES}
    order = [t for t, _, _ in TISSUES]
    cmap = matplotlib.colors.ListedColormap([col_of[t] for t in order])
    idx = np.vectorize({t: i for i, t in enumerate(order)}.get)(T)
    ti = [c for c, t in zip(cells, tags) if t == 8]
    cx = np.mean([c.mean(axis=0) for c in ti], axis=0) if ti else np.array([132., 190.])
    zooms = [((14.0, 102.0), 40.0, "left temporal region"), ((cx[0], cx[1]), 26.0, "titanium plate and screws"),
             ((137.2, 183.4), 6.0, "screw")]
    fig, axes = plt.subplots(1, 4, figsize=(6.3, 2.05), gridspec_kw=dict(width_ratios=[196 / 238, 1, 1, 1]))
    a0 = axes[0]
    a0.pcolormesh(np.arange(T.shape[1] + 1), np.arange(T.shape[0] + 1), idx, cmap=cmap, vmin=-0.5,
                  vmax=len(order) - 0.5, shading="flat", edgecolors="none", rasterized=False)
    a0.set_xlim(0, T.shape[1]); a0.set_ylim(0, T.shape[0])
    a0.set_title("tissues (1 mm voxels)", fontsize=7, pad=2)
    fc = [col_of[t] for t in tags]
    for k, (a, (c, w, title)) in enumerate(zip(axes[1:], zooms)):
        x0, y0 = c[0] - w / 2, c[1] - w / 2
        sel = [m for m, p in enumerate(cells) if p[:, 0].max() > x0 and p[:, 0].min() < x0 + w and
               p[:, 1].max() > y0 and p[:, 1].min() < y0 + w]
        a.add_collection(PolyCollection([cells[m] for m in sel], facecolors=[fc[m] for m in sel], edgecolors="#1a1a1a",
                                        linewidths=0.25 if w > 10 else 0.3))
        a.set_xlim(x0, x0 + w); a.set_ylim(y0, y0 + w)
        a.set_title(title + f" ({w:g} mm)", fontsize=7, pad=2)
        if k < 2:
            nx, ny, nw = zooms[k + 1][0][0] - zooms[k + 1][1] / 2, zooms[k + 1][0][1] - zooms[k + 1][1] / 2, zooms[k + 1][1]
            if k == 1:
                a.add_patch(Rectangle((nx, ny), nw, nw, fill=False, ec=plt.cm.Reds(0.85), lw=0.8))
        a0.add_patch(Rectangle((x0, y0), w, w, fill=False, ec=plt.cm.Reds(0.85), lw=0.7)) if k < 2 else None
    for a in axes:
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_linewidth(0.6)
    handles = [Patch(fc=c, ec="#6a6a66", lw=0.4, label=n) for _, n, c in TISSUES]
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=6, frameon=False, handlelength=1.0,
               columnspacing=0.9, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "head_model.pdf")




def head_perf():
    """Stable steps relative to ERK4 and CPU times (1 thread over 12 us, 6 threads over 120 us)."""
    sch = (("erk", "ERK4", "erk", "erk"), ("lts", "LTS", "lts", "lts"), ("mlts", "MLTS", "mlts", "mlts16"),
           ("li", "LI", "lis", "lis9"))
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), gridspec_kw=dict(width_ratios=[1.1, 1, 1]))
    ax = axs[0]
    for k, (lv, title) in enumerate((("l1", "intact"), ("l2", "with plate"))):
        e = col(f"cfl/{lv}_erk", "dt")
        for j, (s, lab, cf, run) in enumerate(sch):
            v = col(f"cfl/{lv}_{cf}", "dt")
            if v is None or e is None:
                continue
            x = k * 5 + j
            ax.bar(x, v[0] / e[0], color=SCH[s][1], width=0.8, edgecolor="#1a1a1a", linewidth=0.4,
                   label=lab if k == 1 else None)
            ax.text(x, v[0] / e[0] * 1.12, f"{v[0] / e[0]:.0f}" if v[0] / e[0] >= 2 else f"{v[0] / e[0]:.1f}",
                    ha="center", va="bottom", fontsize=6)
    ax.set_yscale("log"); ax.set_ylim(0.5, 400)
    ax.set_xticks([1.5, 6.5]); ax.set_xticklabels(["intact head", "with plate"])
    ax.set_ylabel(r"$\Delta t_{\rm crit}/\Delta t_{\rm crit}^{\rm ERK4}$", labelpad=1)
    ax.grid(True, axis="y", which="major", **GRID); ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    legend(ax, loc="upper left")
    for ax, base, title in ((axs[1], "head1", "1 thread, $t\\in[0,12]\\,\\mu$s"), (axs[2], "head", "6 threads, $t\\in[0,120]\\,\\mu$s")):
        e = col(f"{base}/l2_erk", "cpu")
        for j, (s, lab, cf, run) in enumerate(sch):
            v = col(f"{base}/l2_{run}", "cpu")
            if v is None:
                continue
            ax.bar(j, v[-1], color=SCH[s][1], width=0.7, edgecolor="#1a1a1a", linewidth=0.4)
            ax.text(j, v[-1] * 1.03, f"{v[-1]:.0f} s" + ("" if s == "erk" else f"\n($\\times${e[-1] / v[-1]:.1f})"),
                    ha="center", va="bottom", fontsize=5.8)
        ax.set_xticks(range(len(sch))); ax.set_xticklabels([x[1] for x in sch])
        ax.set_ylim(0, max(col(f"{base}/l2_erk", "cpu")[-1] * 1.25, 1))
        ax.set_title(title, fontsize=7.5, pad=3); ax.set_ylabel("CPU time (s)", labelpad=1)
        ax.grid(True, axis="y", which="major", **GRID); ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    fig.tight_layout(pad=0.2, w_pad=0.8)
    save(fig, "head_perf.pdf")


def head_sensors():
    sens = (("water", "acoustic_pressure", "water, near the source"), ("skull_left", "elastic_velocity", "left temporal bone"),
            ("brain_center", "acoustic_pressure", "brain, centre"), ("under_plate", "acoustic_pressure", "under the plate"))
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 3.4), sharex=True)
    for ax, (name, fld, title) in zip(axes.ravel(), sens):
        for s, run, ls, lw in (("erk", "erk", "-", 1.6), ("lts", "lts", "--", 0.8), ("mlts", "mlts16", "-", 0.8), ("li", "lis9", ":", 1.0)):
            f = os.path.join(R, f"head/l2_{run}", f"c_{name}.csv")
            if not os.path.exists(f):
                continue
            Araw = [l for l in open(f) if not l.startswith("#") and l.strip()]
            Am = np.array([[float(x) for x in l.split(",")] for l in Araw[1:]])
            t = Am[:, 0]
            y = Am[:, 1] if Am.shape[1] == 2 else np.sqrt((Am[:, 1:] ** 2).sum(axis=1))
            lab, c, _ = SCH[s]
            ax.plot(t, y, ls, color=c, lw=lw, label=lab, alpha=0.45 if s == "erk" else 1.0)
        ax.grid(True, which="major", **GRID)
        ax.set_title(title, fontsize=7.5, pad=3)
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        ax.set_ylabel(r"$p$" if fld == "acoustic_pressure" else r"$|v|$", labelpad=1)
    for ax in axes[1]:
        ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1)
    legend(axes[0][0], loc="upper right")
    fig.tight_layout(pad=0.2, h_pad=0.6, w_pad=1.0)
    save(fig, "head_sensors.pdf")


def head_snapshots():
    d = os.path.join(R, "head", "l2_snap")
    files = sorted(glob.glob(os.path.join(d, "c_snapshot_*_samples.dat")), key=lambda s: int(re.search(r"_(\d+)_samples", s).group(1)))
    cf = os.path.join(d, "c_cells.dat")
    if not files or not os.path.exists(cf):
        print("head snapshots missing"); return
    cells = read_cells(cf)
    acoustic = np.array([c[0] for c in cells]) == 1
    T = np.load(os.path.join(HEAD_MESHES, "head_s80_plate_tissues.npy"))
    res = 0.125
    W, H = T.shape[1], T.shape[0]
    nx, ny = int(round(W / res)), int(round(H / res))
    fig, axes = plt.subplots(1, len(files), figsize=(6.3, 6.3 / len(files) * H / W + 0.35))
    for ax, f in zip(np.atleast_1d(axes), files):
        tt = float(re.search(r"t requested ([0-9.e+-]+)", open(f).readline()).group(1))
        img = np.full((ny, nx), np.nan)
        for k, l in enumerate(l for l in open(f) if not l.startswith("#")):
            v = np.array(l.split(), float)
            x0, x1, y0, y1 = v[:4]; vals = v[4:]
            m = int(round(math.sqrt(len(vals))))
            if not acoustic[k]:
                continue
            i0, i1 = int(round(x0 / res)), int(round(x1 / res)); j0, j1 = int(round(y0 / res)), int(round(y1 / res))
            blk = vals.reshape(m, m)
            ni, nj = i1 - i0, j1 - j0
            ri = (np.arange(nj) * m) // nj; ci = (np.arange(ni) * m) // ni
            img[j0:j1, i0:i1] = blk[np.ix_(ri, ci)]
        vmax = np.nanpercentile(np.abs(img), 99.5) or 1.0
        ax.imshow(img, origin="lower", extent=(0, W, 0, H), cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="bilinear")
        # elastic cells (bone, titanium) in grey, outline of the brain
        ax.add_collection(PolyCollection([c[1] for c in cells if c[0] == 0], facecolors="#cfc8b4", edgecolors="none"))
        brain = np.isin(T, [5, 6, 7]).astype(float)
        ax.contour(np.arange(W) + 0.5, np.arange(H) + 0.5, brain, levels=[0.5], colors="#3a3a3a", linewidths=0.3)
        ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        ax.set_title(rf"$t={tt:.0f}\,\mu$s", fontsize=7.5, pad=3)
    fig.tight_layout(pad=0.2, w_pad=0.4)
    fig.savefig(os.path.join(OUT, "head_snapshots.pdf"), bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)
    print("ok head_snapshots.pdf")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in (formulations, reduced, stability, meshes_cart, convergence, time_orders, convergence_2e, time_orders_2e,
              lshape_meshes, lshape, lshape_2e, lshape_formulations, head_model, head_perf, head_sensors, head_snapshots):
        try:
            f()
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f.__name__, "failed:", repr(e))
    for m in MISSING:
        print("WARNING: no room for the slope of", m)
    print("written to", OUT)
