#!/usr/bin/env python3
"""make_article.py -- figures (PDF) of Sections 5 and 6 of the LTS_HHO preprint (main.tex),
with one plotting style: the one of the a priori / convergence papers on Leray-Lions problems (Computer Modern fonts at
the printed size, panels 3.0 in wide for 0.48\\textwidth, full-width figures 6.3 in, dotted major grid, observed slope of
every segment, dashed reference slopes).

Color and marker code of the time schemes, the same in every figure:
    ERK4  black, circles        LTS   orange, squares      MLTS  dark red, triangles
    LI    blue, diamonds        SDIRK grey, inverted triangles (implicit reference)

Data: <results> as written by run_article.sh (one directory per run: c.dat, log.txt, time.txt, sensors, snapshots):
    alg/, form/ 6.6    rho/, stab/ 6.2    conv/, tconv/, meshes/ 6.4    ls/ 6.5    cfl/, head/, head1/ 6.7 (earlier head model)
and the reduced models (5): reduced_models/plan_extra.npz, reduced_lts_erk.npz, reduced_li_stage.npz.
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
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.35), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]
    for o, lab, c, mk in (("o1", r"first order, $k=3$", SCH["mlts"][1], "o"), ("o2", r"second order, mixed $k=2$", "#4a4a46", "s")):
        nu, ms = [], []
        for l in range(4):
            D = dat(f"alg/{o}_l{l}")
            if D is None:
                continue
            names, rows = D
            n_steps = round(1.0 / rows[-1, names.index("dt")])
            # integration time per step: difference of the runs with n and n + 80000 / 2^l steps (the setup cancels;
            # the second-order setup estimates omega_max and the splitting radius), best of the repetitions (_r2);
            # without the long runs, wall time per step
            extra = 80000 // 2 ** l
            per = []
            for rep in ("", "_r2"):
                w, _ = wall(f"alg/{o}_l{l}{rep}")
                wl, _ = wall(f"alg/{o}_l{l}_long{rep}")
                if w is not None and wl is not None:
                    per.append((wl - w) / extra)
            if not per:
                w, _ = wall(f"alg/{o}_l{l}")
                per = [w / n_steps]
            nu.append(8 * 2 ** l); ms.append(1e3 * min(per))
        if nu:
            ax.plot(np.array(nu) ** 2, ms, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    style_log(ax, "cells", "wall time per ERK4 step (ms)")
    legend(ax, loc="upper left")
    ax = axs[1]
    # wall time (the second-order runs have no CPU column), one thread (hatched) and six threads
    labels, xs = [], []
    k = 0
    # second order: MLTS with the projected condensation only (runs mltsp, faces.projected), labelled MLTS
    for o, name, schemes in (("o1", "first order", ("erk", "mlts")), ("o2", "second order", ("erk", "mltsp"))):
        for s in schemes:
            color = SCH["mlts" if s == "mltsp" else s][1]
            for j, (suf, th) in enumerate((("_q1", 1), ("_q", 6))):
                w, _ = wall(f"form/{o}_{s}{suf}")
                v = w if w is not None else np.nan
                xi = k + (j - 0.5) * 0.36
                ax.bar(xi, v, color=color, width=0.34, edgecolor="#1a1a1a", linewidth=0.4,
                       hatch="////" if th == 1 else None)
                if np.isfinite(v):
                    ax.text(xi, v * 1.12, f"{v:.0f}" if v >= 10 else f"{v:.1f}", ha="center", va="bottom", fontsize=5.5)
            labels.append(f"{name}\n{SCH['mlts' if s == 'mltsp' else s][0]}"); xs.append(k)
            k += 1
    vals = [b.get_height() for b in ax.patches]
    ax.set_yscale("log"); ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=5.8)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, fc="white", ec="#1a1a1a", lw=0.4, hatch="////", label="1 thread"),
                       Rectangle((0, 0), 1, 1, fc="white", ec="#1a1a1a", lw=0.4, label="6 threads")],
              loc="upper left", fontsize=6, frameon=False)
    ax.set_ylabel(r"wall time (s), $t_f=1$", labelpad=1); ax.grid(True, axis="y", which="major", **GRID)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    ax.set_ylim(top=max(v for v in vals if np.isfinite(v)) * 6)
    ax.set_title(r"$64\times64$ cells, central patch refined 4 times ($p=16$)", fontsize=7, pad=3)
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "form.pdf")


# --------------------------------------------------------------------------------------------- 5 reduced models
# reduced models: the two-level LTS is MLTS with L = 2 levels
RED_LABEL = {"lts": r"MLTS, $L=2$", "mlts": r"MLTS, $L=3$ ($p_1=p_2=p$)"}


def reduced():
    f1 = os.path.join(HERE, "reduced_models", "plan_extra.npz")
    f2 = os.path.join(HERE, "reduced_models", "reduced_lts_erk.npz")
    f3 = os.path.join(HERE, "reduced_models", "reduced_li_stage.npz")
    if not all(os.path.exists(f) for f in (f1, f2, f3)):
        print("reduced models: missing npz"); return
    A, B, C = (np.load(f, allow_pickle=True) for f in (f1, f2, f3))
    P = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]); P3 = np.array([1, 2, 4, 8, 16, 32, 64])
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), sharey=True)
    for ax, key, mkey, title in ((axs[0], "order2", "order2", "model 1: second order"),
                                  (axs[1], "order1_nobuf", None, "model 2: first order, no overlap"),
                                  (axs[2], "order1_buf", "order1_buf", "model 3: first order, overlap")):
        rows = [("erk", B[f"erk_{key}_rho"], P), ("lts", B[f"lts_{key}_rho"], P), ("li", C[f"li_stage_{key}_rho"], P)]
        if mkey:
            rows.insert(2, ("mlts", A[f"mlts_{mkey}_alg4_rho"], P3))
        for s, D, PP in rows:
            lab, c, mk = SCH[s]
            ax.plot(PP, D.max(axis=0), "-" + mk, color=c, lw=1.0, ms=3.0, label=RED_LABEL.get(s, lab))
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
        ax.set_xlabel(r"$p$", labelpad=1); ax.set_title(title, fontsize=7, pad=3)
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
        ax.axhline(1.0, color="#6a6a66", lw=0.6, ls="--")
    axs[0].set_ylabel(r"$\max_\varepsilon\,\rho(C_p(\Delta t_{\rm ref}))$", labelpad=1)
    legend(axs[0], loc="upper left")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "red_rho.pdf")
    # critical step (largest step below which every step is stable), relative to that of RK4 at p = 1, smallest over eps
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.2), sharey=True)
    for ax, key, title in ((axs[0], "order2", "model 1: second order"), (axs[1], "order1_buf", "model 3: first order, overlap")):
        rows = [("erk", B[f"erk_{key}_dtc"], P), ("lts", B[f"lts_{key}_dtc"], P),
                ("mlts", A[f"mlts_{key}_alg4_dtc"], P3), ("li", C[f"li_stage_{key}_dtc"], P)]
        for s, D, PP in rows:
            lab, c, mk = SCH[s]
            ax.plot(PP, D.min(axis=0), "-" + mk, color=c, lw=1.0, ms=3.0, label=RED_LABEL.get(s, lab))
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
        ax.set_xlabel(r"$p$", labelpad=1); ax.set_title(title, fontsize=7, pad=3)
        ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    axs[0].set_ylabel(r"$\min_\varepsilon\,\Delta t_c(p)/\Delta t_{\rm ref}$", labelpad=1)
    legend(axs[0], loc="lower left")
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "red_dtc.pdf")


# ------------------------------------------------------------------------------------------------ 6.2 stability
def stability():
    L = [0, 1, 2, 3, 4, 5]
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.35))
    ax = axs[0]
    for s in ("erk", "lts", "mlts", "li"):
        pp, rr = [], []
        for l in L:
            r = col(f"rho/L{l}_{s}", "rho")
            if r is not None and r[0] > 0:
                pp.append(2 ** l); rr.append(r[0])
        lab, c, mk = SCH[s]
        ax.plot(pp, rr, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.grid(True, which="major", **GRID)
    ax.axhline(1.0, color="#6a6a66", lw=0.6, ls="--")
    ax.set_xlabel(r"$p=h_{\max}/h_{\min}$", labelpad=1); ax.set_ylabel(r"$\rho(\mathcal{C})$, $\Delta t=6.4\times10^{-3}$", labelpad=1)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)
    legend(ax, loc="upper left")
    ax = axs[1]
    for s in ("erk", "lts", "mlts", "li"):
        pp, dd = [], []
        for l in range(5):
            d = f"stab/L{l}_{'lis' if s == 'li' and l > 0 else s}"
            v = col(d, "dt")
            if v is None or (s == "li" and l == 0):
                continue
            pp.append(2 ** l); dd.append(v[0])
        lab, c, mk = SCH[s]
        ax.plot(pp, dd, "-" + mk, color=c, lw=1.0, ms=3.2, label=lab)
    # second-order formulation (setting (2E), stab2/): ERK4, LTS and MLTS with the projected condensation (dashed, hollow)
    second = False
    for s, key in (("erk", "erk"), ("lts", "ltsp"), ("mlts", "mltsp")):
        pp, dd = [], []
        for l in range(5):
            v = col(f"stab2/L{l}_{key}", "dt")
            if v is None:
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
    ax.set_ylim(2.5e-4, 2.5e-2)
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "stab.pdf")


# -------------------------------------------------------------------------------------- 6.4 smooth convergence
def meshes_cart():
    names = [("u", l) for l in (0, 1, 2)] + [("g", l) for l in (0, 1, 2)]
    files = [os.path.join(R, "meshes", f"cart_{k}{l}", "c_cells.dat") for k, l in names]
    if not all(os.path.exists(f) for f in files):
        print("meshes of 6.4: missing c_cells.dat"); return
    fig, axes = plt.subplots(2, 3, figsize=(6.3, 4.3))
    for ax, f, (k, l) in zip(axes.ravel(), files, names):
        cells = read_cells(f)
        lw = 0.6 if len(cells) <= 600 else (0.3 if len(cells) <= 3000 else 0.1)
        zoom = k == "g" and l == 2
        if zoom:
            lw = 0.3
        ax.add_collection(LineCollection(edges([c for _, c in cells]), colors="#1a1a1a", linewidths=lw))
        lo, hi = (0.3, 0.7) if zoom else (0.0, 1.0)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        n = 8 * 2 ** l
        ax.set_title((f"${n}\\times{n}$" if k == "u" else f"${n}\\times{n}$, $L=3$") + f", {len(cells)} cells"
                     + (r", zoom on $[0.3,0.7]^2$" if zoom else ""), fontsize=7.5, pad=3)
    fig.tight_layout(pad=0.2, w_pad=1.0, h_pad=0.8)
    save(fig, "cart_meshes.pdf")


def convergence():
    S = {"erk": "conv/erk_u", "mlts": "conv/mlts_g", "li": "conv/li_g"}
    D = {k: dat(v) for k, v in S.items()}
    if any(v is None for v in D.values()):
        print("convergence: missing runs"); return
    g = lambda k, n: D[k][1][:, D[k][0].index(n)]
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5))
    for ax, fld, yl in ((axs[0], "acoustic_q", r"$\Vert v-v_h\Vert_{L^2(\Omega)}$"),
                        (axs[1], "acoustic_p", r"$\Vert p-p_h\Vert_{L^2(\Omega)}$")):
        series = [(g("erk", "h"), g("erk", fld), ("ERK4, uniform", SCH["erk"][1], "o"), "--"),
                  (g("mlts", "h"), g("mlts", fld), ("MLTS, refined", SCH["mlts"][1], "^"), "-"),
                  (g("li", "h"), g("li", fld), ("LI, refined", SCH["li"][1], "D"), "-")]
        # LI coincides with MLTS to three digits: slopes of ERK4 and MLTS only
        conv_panel(ax, fig, series, 4, r"$h_{\max}$", yl, r"$k=3$, $\Delta t\propto h_{\max}$", "conv_cart.pdf",
                   labelled=(0, 1))
    fig.tight_layout(pad=0.2, w_pad=1.2)
    save(fig, "conv_cart.pdf")


def time_orders():
    series = []
    for key, d in (("mlts", "tconv/mlts"), ("li", "tconv/li")):
        D = dat(d)
        if D is None:
            continue
        dt, e = D[1][:, D[0].index("dt")], D[1][:, D[0].index("mass_norm")]
        ok = np.isfinite(e) & (e < 1)
        series.append((dt[ok], e[ok], key, "-"))
    if series:
        fig, ax = plt.subplots(figsize=(3.0, 2.3))
        conv_panel(ax, fig, series, 4, r"$\Delta t$", r"$\Vert U_h(t_f)-U_h^{\rm ref}(t_f)\Vert$",
                   r"$16\times16$, $L=3$ ($p=8$)", "tconv.pdf", legend_loc="lower right")
        fig.tight_layout(pad=0.2)
        save(fig, "tconv.pdf")


# ----------------------------------------------------------------------------------------------- 6.5 L-shape
def lshape_meshes():
    f = os.path.join(LS_MESHES, "lshape_mu1-5_n8.txt")
    if not os.path.exists(f):
        print("L-shape mesh missing"); return
    _, cells, _ = read_polygons(f)
    smin = min(np.ptp(c[:, 0]) for c in cells)
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.2))
    for ax, w in zip(axes, (0.5, 0.5 / 32, 0.5 / 2 ** round(math.log2(0.5 / (6 * smin))))):
        sel = [c for c in cells if np.min(np.abs(c).max(axis=1)) <= 1.5 * w]
        ax.add_collection(LineCollection(edges(sel), colors="#1a1a1a", linewidths=0.5 if w == 0.5 else 0.4))
        ax.set_xlim(-w, w); ax.set_ylim(-w, w); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.add_patch(Rectangle((0, 0), w, w, fc="#e9e9e6", ec="none", zorder=0))
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        ax.set_title(f"{len(cells)} cells" if w == 0.5 else rf"zoom $\times{0.5 / w:.0f}$", fontsize=7.5, pad=3)
    axes[0].text(0.03, 0.03, rf"$h_{{\max}}/h_{{\min}}={0.0625 / smin:.0f}$", transform=axes[0].transAxes, fontsize=6, bbox=BOX)
    fig.tight_layout(pad=0.2, w_pad=1.0)
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


def lshape():
    lv = (1, 2, 3, 4, 5)
    # LI at half the step (dt = 1/(40n)): the error constant of SDIRK(3,4) is about 30 times that of RK4 (ls/li2_*)
    data = {"mlts": ls_series("mlts", lv), "li": ls_series("li2", lv)}
    erk = ls_series("erk", (1,))
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.35), gridspec_kw=dict(width_ratios=[1, 1, 1]))
    for ax, idx, yl in ((axs[0], 2, r"$\Vert v-v_h\Vert_{L^2(\Omega)}$"), (axs[1], 1, r"$\Vert p-p_h\Vert_{L^2(\Omega)}$")):
        series = [(data[s][0], data[s][idx], s, "-") for s in ("mlts", "li") if len(data[s][0]) >= 2]
        conv_panel(ax, fig, series, 4, r"$h_{\max}$", yl, None, f"lshape[{idx}]", legend_loc="lower right")
        if len(erk[0]):
            ax.plot(erk[0], erk[idx], "o", color=SCH["erk"][1], ms=3.5, zorder=4, label="ERK4")
            legend(ax, loc="lower right")
    ax = axs[2]
    for s in ("mlts", "li"):
        hh, pp, qq, cc, _ = data[s]
        lab, c, mk = SCH[s]
        ax.plot(cc, qq, "-" + mk, color=c, lw=1.0, ms=3.0, label=lab)
    if len(erk[0]):
        ax.plot(erk[3], erk[2], "o", color=SCH["erk"][1], ms=3.5, label="ERK4")
    style_log(ax, "CPU time (s)", r"$\Vert v-v_h\Vert_{L^2(\Omega)}$")
    legend(ax, loc="upper right")
    fig.tight_layout(pad=0.2, w_pad=0.8)
    save(fig, "lshape.pdf")


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
    for f in (formulations, reduced, stability, meshes_cart, convergence, time_orders, lshape_meshes, lshape, head_model,
              head_perf, head_sensors, head_snapshots):
        try:
            f()
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f.__name__, "failed:", repr(e))
    for m in MISSING:
        print("WARNING: no room for the slope of", m)
    print("written to", OUT)
