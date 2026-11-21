#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures of Section 6.7 of the LTS_HHO preprint (main.tex): transcranial focused ultrasound at 250 kHz in the
seven-tissue head model with and without a titanium fixation plate (inputs article/inputs/tfus_fit.lua (fitted solid
interfaces) and tfus.lua (voxel staircase), jobs fit/* and tfus/* of article/jobs.txt). Reads the runs of the directory
given as first argument (fit/plate_li, intact_li, water_li, plate_li_dt2, plate_mlts; tfus/plate_mlts, intact_mlts,
water_mlts) and writes into the second one:
    tfus_model.pdf      tissues, transducer, focus, beam axis; zoom on the plate
    tfus_model_steps.pdf  (2026-11-15) construction of the mesh: T1 slice with the tissues in transparency, quadtree
                        on the pixels (staircase) and cells split along the smoothed skull, over the T1 slice (needs the
                        MRI of ~/Downloads/mri_dataset, or the folder given by the environment variable TFUS_MRI_DATA)
    tfus_snapshots.pdf  pressure in the fluids and mean stress -tr(sigma)/2 in the solids, one colour scale
    tfus_skull.pdf      zoom on the skull and the plate: mean stress (P waves) and shear stress (S waves)
    tfus_focus.pdf      peak pressure maps (water, intact skull, plate) and profiles along and across the axis
    tfus_signals.pdf    pressure at the focus (three models; schemes and steps with the plate), discrete energy (fitted
                        and staircase interfaces), costs
    tfus_cost.txt       Table tfus-cost: best of three integration times (6 threads: fit/plate_*; 1 thread: t1/*, 8 us
                        extrapolated, the factorization of LI counted once), gains, memory
Memories: with --memory <runs>, the algebraic memory of the time integration (output.algebraic_memory: operator, state
and what the scheme allocates and keeps; <runs>/<fit|t1>/<run>/log.txt, campaign lts_work/algmem, same jobs, 3 steps);
without it, the peak memory (maximum resident set size, time.txt).
With --short (2026-11-15), the six-thread times of MLTS and LI also come from 8 us runs without outputs
(fit/plate_mlts_8us, fit/plate_li_8us), extrapolated as on one thread: the same protocol for every scheme.
With --three-snapshots (2026-11-16), tfus_snapshots.pdf has three panels: t = 4 us, before the burst reaches the skull
(run fit/plate_li_early: fit/plate_li stopped after 70 steps, job of article/pardiso/jobs_snap.txt), 20 and 44 us.
With --signals-xmin <v> (2026-11-17, off by default: 44): left panel of tfus_signals.pdf (--signals-row) from t = v us.
Usage: python3 make_tfus.py <runs> <figures> [--memory <runs>] [--short] [--three-snapshots]
"""
import math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

ARGS, MEM = sys.argv[1:], None
# --short (2026-11-15, off by default): on six threads too, every scheme timed over the first 8 us without outputs
# (fit/plate_erk_8us, fit/plate_mlts_8us, fit/plate_li_8us) and extrapolated, as on one thread (t1/)
SHORT = "--short" in ARGS
# --three-snapshots (2026-11-16, off by default): Fig. tfus-snap with t = 4 (fit/plate_li_early), 20 and 44 us
THREE_SNAPSHOTS = "--three-snapshots" in ARGS
# --skull-strip (2026-11-16, off by default): Fig. tfus-skull as strips along the plate (rotated view, plate horizontal,
# scalp up), shear stress on a common scale for the three times
SKULL_STRIP = "--skull-strip" in ARGS
# --focus-row (2026-11-16, off by default): Fig. tfus-focus in one row, the three peak-pressure maps and the profile
# along the axis in dB, from the geometric focus; --signals-row: Fig. tfus-signals without the focal signals of the
# three configurations and without the cost bars (Table tfus-cost): accuracy check and discrete energy only
FOCUS_ROW = "--focus-row" in ARGS
SIGNALS_ROW = "--signals-row" in ARGS
# --mesh-zooms (2026-11-16, off by default): tfus_model_steps.pdf with three zooms on the same window, staircase,
# quadtree with split cells and Voronoi with split cells (instead of the staircase and the split quadtree)
MESH_ZOOMS = "--mesh-zooms" in ARGS
# --fine-timing (2026-11-16, off by default): Table tfus-cost from the timed runs of the meshes twice finer
# (lts_work/head_fine/chain_T.sh): ERK(4) over 2 us (plate_erk_2us, t1/erk_2us), MLTS and LI over 278 of 2788 steps
FINE_TIMING = "--fine-timing" in ARGS
# --paper-labels (2026-11-18, off by default, fourth review of the paper): units in brackets ([us], [dB]) in
# tfus_signals, tfus_focus and tfus_polytopal; no panel titles in tfus_signals and tfus_polytopal (the captions name the
# panels); in tfus_polytopal, "pressure at the geometric focus" and slightly smaller fonts in the right panel; in
# tfus_snapshots, "mean pressure" for -tr(sigma)/2 instead of "mean stress"
PAPER_LABELS = "--paper-labels" in ARGS
ARGS = [a for a in ARGS if a not in ("--short", "--three-snapshots", "--skull-strip", "--focus-row", "--signals-row",
                                     "--mesh-zooms", "--fine-timing", "--paper-labels")]
# --paper-labels: upright micro sign (U+00B5) in the unit, as \textmu in the preprint (2026-11-19, was $\mu$, italic)
T_LABEL = "$t$ [µs]" if PAPER_LABELS else r"$t$ ($\mu$s)"


def minus(s):
    """--paper-labels (2026-11-19): typographic minus sign (U+2212) instead of the hyphen, for numbers only."""
    return s.replace("-", "−") if PAPER_LABELS else s


def t_title(t):
    """Time in the panels: "t = 4 [mu s]" with --paper-labels (2026-11-19, unit in brackets), else "t = 4 mu s"."""
    return f"$t={t:.0f}$ [µs]" if PAPER_LABELS else rf"$t={t:.0f}\,\mu$s"


DB = " [dB]" if PAPER_LABELS else " dB"   # --paper-labels (2026-11-19): unit in brackets in the titles of tfus_focus
if "--memory" in ARGS:
    _k = ARGS.index("--memory")
    MEM, ARGS = ARGS[_k + 1], ARGS[:_k] + ARGS[_k + 2:]
# --fine-meshes <dir> (2026-11-16, off by default): the meshes twice finer of <dir> (head_s80_{fit,hf,vor}2_*,
# lts_work/head_fine/meshes) instead of those of mesh_generation/brain_mri/meshes; the tissue rasters are identical
FINE = None
if "--fine-meshes" in ARGS:
    _k = ARGS.index("--fine-meshes")
    FINE, ARGS = ARGS[_k + 1], ARGS[:_k] + ARGS[_k + 2:]
# --energy-ymin <v> (2026-11-17, off by default): energy axis of Fig. tfus-signals (--signals-row) from v to 1.005,
# instead of the automatic range, which magnifies the small losses of the fitted interfaces
ENERGY_YMIN = None
if "--energy-ymin" in ARGS:
    _k = ARGS.index("--energy-ymin")
    ENERGY_YMIN, ARGS = float(ARGS[_k + 1]), ARGS[:_k] + ARGS[_k + 2:]
# --signals-xmin <v> (2026-11-17, off by default: 44): time axis of the left panel of Fig. tfus-signals (--signals-row)
# from v to 56 us instead of 44 to 56 us
SIGNALS_XMIN = 44
if "--signals-xmin" in ARGS:
    _k = ARGS.index("--signals-xmin")
    SIGNALS_XMIN, ARGS = float(ARGS[_k + 1]), ARGS[:_k] + ARGS[_k + 2:]
RUNS = ARGS[0] if len(ARGS) > 0 else "."
OUT = ARGS[1] if len(ARGS) > 1 else "."
MEMORY_LABEL = "algebraic memory" if MEM else "peak memory"
HERE = os.path.dirname(os.path.abspath(__file__))
MESHES = os.path.join(HERE, "..", "..", "..", "src", "mesh_generation", "brain_mri", "meshes")


def mesh_path(name):
    """A mesh file (or tissue raster) of the preprint, or with --fine-meshes its twice finer version."""
    if FINE:
        for a in ("fit", "hf", "vor"):
            name = name.replace(f"head_s80_{a}_", f"head_s80_{a}2_")
        return os.path.join(FINE, name)
    return os.path.join(MESHES, name)
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.5,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "mathtext.fontset": "cm", "savefig.dpi": 300})
GRID = dict(color="#c8c8c4", lw=0.4, ls=":")
SCH = {"erk": ("ERK4", "#1a1a1a"), "lts": ("LTS", "#dd6b20"), "mlts": ("MLTS", "#b7202e"), "li": ("LI", "#2b6cb0")}
# geometry of the input (tfus.lua)
FOCUS = np.array([120.5, 179.6]); AXIS = math.radians(44.3); RADIUS = 80.0; APERTURE = math.radians(25.0)
FREQ, C0 = 0.25, 1.5
TISSUES = [(3, "water", "#cfe3f1"), (4, "scalp", "#e9b9a8"), (1, "cortical bone", "#f2eee2"),
           (2, "cancellous bone", "#c9b98f"), (5, "CSF", "#7fb2d9"), (6, "gray matter", "#9c9c9c"),
           (7, "white matter", "#efefef"), (8, "titanium", "#2b2b2b")]


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("ok", name)


def run(d, kind="fit"):
    return os.path.join(RUNS, kind, d)


def first_run(kind, *names):
    """The first existing run among names, else the first one. The fine campaign (2026-11-16) has the staircase runs
    with LI (tfus/*_li, not tfus/*_mlts) and the Voronoi run with MLTS (vor/plate_mlts, not vor/plate_li)."""
    for r in names:
        if os.path.isdir(run(r, kind)):
            return run(r, kind)
    return run(names[0], kind)


def read_samples(path):
    """cells (x0, x1, y0, y1) and their m x m samples (rows bottom to top)."""
    D = np.loadtxt(path, comments="#")
    m = int(round(math.sqrt(D.shape[1] - 4)))
    return D[:, :4], D[:, 4:].reshape(-1, m, m)


def raster(boxes, vals, X0, X1, Y0, Y1, px, fill=np.nan):
    """Image of the cell samples on [X0, X1] x [Y0, Y1] with square pixels px (each pixel takes the sample of the
    sub-square of the cell containing its centre)."""
    nx, ny = int(round((X1 - X0) / px)), int(round((Y1 - Y0) / px))
    img = np.full((ny, nx), fill)
    xc = X0 + (np.arange(nx) + 0.5) * px
    yc = Y0 + (np.arange(ny) + 0.5) * px
    m = vals.shape[1]
    sel = np.nonzero((boxes[:, 1] > X0) & (boxes[:, 0] < X1) & (boxes[:, 3] > Y0) & (boxes[:, 2] < Y1))[0]
    for k in sel:
        x0, x1, y0, y1 = boxes[k]
        i0, i1 = np.searchsorted(yc, y0), np.searchsorted(yc, y1)
        j0, j1 = np.searchsorted(xc, x0), np.searchsorted(xc, x1)
        if i1 <= i0 or j1 <= j0:
            continue
        r = np.minimum(((yc[i0:i1] - y0) / (y1 - y0) * m).astype(int), m - 1)
        c = np.minimum(((xc[j0:j1] - x0) / (x1 - x0) * m).astype(int), m - 1)
        img[i0:i1, j0:j1] = vals[k][np.ix_(r, c)]
    return img


def read_cells(cells_file):
    """medium (1 elastic, 0 acoustic) and vertices of every cell, order of the samples (<table>_cells.dat)."""
    med, polys = [], []
    for line in open(cells_file):
        if line.startswith("#"):
            continue
        v = line.split()
        a, nv = int(v[0]), int(v[1])
        polys.append(np.array(v[2:2 + 2 * nv], dtype=float).reshape(nv, 2))
        med.append(1 - a)
    return np.array(med), polys


def raster_cells(polys, vals, X0, X1, Y0, Y1, px, fill=np.nan):
    """Image of the cell samples (m x m on the bounding box of each cell) on [X0, X1] x [Y0, Y1]: the axis-aligned
    rectangles are painted on their box, the cut cells only on the pixels inside their polygon (exact interfaces)."""
    from matplotlib.path import Path
    nx, ny = int(round((X1 - X0) / px)), int(round((Y1 - Y0) / px))
    img = np.full((ny, nx), fill)
    xc = X0 + (np.arange(nx) + 0.5) * px
    yc = Y0 + (np.arange(ny) + 0.5) * px
    m = vals.shape[1]
    later = []
    for k, q in enumerate(polys):
        x0, x1, y0, y1 = q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()
        if x1 <= X0 or x0 >= X1 or y1 <= Y0 or y0 >= Y1:
            continue
        rect = abs(0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1]) - (x1 - x0) * (y1 - y0)) \
            < 1e-9 * (x1 - x0) * (y1 - y0)
        if not rect:
            later.append(k)
            continue
        i0, i1 = np.searchsorted(yc, y0), np.searchsorted(yc, y1)
        j0, j1 = np.searchsorted(xc, x0), np.searchsorted(xc, x1)
        if i1 <= i0 or j1 <= j0:
            continue
        r = np.minimum(((yc[i0:i1] - y0) / (y1 - y0) * m).astype(int), m - 1)
        c = np.minimum(((xc[j0:j1] - x0) / (x1 - x0) * m).astype(int), m - 1)
        img[i0:i1, j0:j1] = vals[k][np.ix_(r, c)]
    for k in later:
        q = polys[k]
        x0, x1, y0, y1 = q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()
        i0, i1 = np.searchsorted(yc, y0), np.searchsorted(yc, y1)
        j0, j1 = np.searchsorted(xc, x0), np.searchsorted(xc, x1)
        if i1 <= i0 or j1 <= j0:
            continue
        XX, YY = np.meshgrid(xc[j0:j1], yc[i0:i1])
        inside = Path(q).contains_points(np.column_stack([XX.ravel(), YY.ravel()]), radius=1e-9).reshape(XX.shape)
        r = np.minimum(((yc[i0:i1] - y0) / (y1 - y0) * m).astype(int), m - 1)
        c = np.minimum(((xc[j0:j1] - x0) / (x1 - x0) * m).astype(int), m - 1)
        sub = img[i0:i1, j0:j1]
        sub[inside] = vals[k][np.ix_(r, c)][inside]
    return img


def solid_outline(ax, med, polys, box, color="#1a1a1a", lw=0.4):
    """Edges between elastic and acoustic cells (the fitted interfaces), from the vertices of the cells."""
    from matplotlib.collections import LineCollection
    X0, X1, Y0, Y1 = box
    owner = {}
    segs = []
    for k, q in enumerate(polys):
        if q[:, 0].max() < X0 - 2 or q[:, 0].min() > X1 + 2 or q[:, 1].max() < Y0 - 2 or q[:, 1].min() > Y1 + 2:
            continue
        n = len(q)
        for e in range(n):
            a, b = tuple(np.round(q[e], 9)), tuple(np.round(q[(e + 1) % n], 9))
            key = (min(a, b), max(a, b))
            if key in owner and owner[key] != med[k]:
                segs.append([a, b])
            owner.setdefault(key, med[k])
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw))


def tissue_underlay(ax, box, alpha=0.35):
    tis = np.load(mesh_path("head_s80_fit_plate_tissues.npy"))
    H, W = tis.shape
    order = [t for t, _, _ in TISSUES]
    cmap = ListedColormap([c for _, _, c in TISSUES])
    idx = np.vectorize({t: i for i, t in enumerate(order)}.get)(tis)
    ax.pcolormesh(np.arange(W + 1), np.arange(H + 1), idx, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5, alpha=alpha,
                  shading="flat", edgecolors="none", antialiased=False)
    ax.set_xlim(box[0], box[1]); ax.set_ylim(box[2], box[3])


def medium_raster(cells_file, X0, X1, Y0, Y1, px):
    """1 in the elastic cells, 0 in the acoustic ones (from <table>_cells.dat: medium, nv, vertices)."""
    boxes, med = [], []
    for line in open(cells_file):
        if line.startswith("#"):
            continue
        v = line.split()
        a, nv = int(v[0]), int(v[1])
        xy = np.array(v[2:2 + 2 * nv], dtype=float).reshape(nv, 2)
        boxes.append((xy[:, 0].min(), xy[:, 0].max(), xy[:, 1].min(), xy[:, 1].max()))
        med.append(1 - a)
    boxes = np.array(boxes)
    vals = np.array(med, dtype=float)[:, None, None] * np.ones((1, 1, 1))
    return raster(boxes, vals, X0, X1, Y0, Y1, px, fill=0.0), boxes, np.array(med)


def outline(ax, mask, extent, color="#1a1a1a", lw=0.35):
    ax.contour(mask, levels=[0.5], extent=extent, colors=color, linewidths=lw, origin="lower")


def transducer(ax, lw=1.6):
    th = np.linspace(AXIS - APERTURE, AXIS + APERTURE, 100)
    ax.plot(FOCUS[0] + RADIUS * np.cos(th), FOCUS[1] + RADIUS * np.sin(th), "-", color="#d69e2e", lw=lw, solid_capstyle="butt")
    ax.plot(*FOCUS, "+", color="#d69e2e", ms=6, mew=1.2)


# --------------------------------------------------------------------------------------------------------- model
def model():
    tis = np.load(mesh_path("head_s80_fit_plate_tissues.npy"))
    H, W = tis.shape
    order = [t for t, _, _ in TISSUES]
    cmap = ListedColormap([c for _, _, c in TISSUES])
    idx = np.vectorize({t: i for i, t in enumerate(order)}.get)(tis)
    fig = plt.figure(figsize=(6.3, 3.1))
    ax = fig.add_axes([0.0, 0.05, 0.42, 0.93])
    ax.pcolormesh(np.arange(W + 1), np.arange(H + 1), idx, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5, shading="flat",
                  edgecolors="none", antialiased=False)
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal")
    transducer(ax)
    e = np.array([math.cos(AXIS), math.sin(AXIS)])
    ax.plot([FOCUS[0], FOCUS[0] + RADIUS * e[0]], [FOCUS[1], FOCUS[1] + RADIUS * e[1]], ":", color="#d69e2e", lw=0.6)
    ax.add_patch(Rectangle((138, 196), 26, 26, fill=False, ec="#b7202e", lw=0.8))
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(0.03, 0.97, f"{W} x {H} mm", transform=ax.transAxes, va="top", fontsize=6)
    # zoom on the plate: cells of the mesh around it
    ax2 = fig.add_axes([0.45, 0.05, 0.53, 0.93])
    pts, cells = read_mesh(mesh_path("head_s80_fit_plate.txt"), (138, 164, 196, 222))
    from matplotlib.collections import PolyCollection
    col = {t: c for t, _, c in TISSUES}
    ax2.add_collection(PolyCollection([pts[r] for r, _ in cells], facecolors=[col[t] for _, t in cells], edgecolors="k",
                                      linewidths=0.08))
    ax2.set_xlim(138, 164); ax2.set_ylim(196, 222); ax2.set_aspect("equal"); ax2.set_xticks([]); ax2.set_yticks([])
    for s in ax2.spines.values():
        s.set_edgecolor("#b7202e"); s.set_linewidth(0.8)
    ax2.text(0.02, 0.02, "cells cut along the smoothed skull\nand the exact plate and screws", transform=ax2.transAxes,
             va="bottom", fontsize=6, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    handles = [Rectangle((0, 0), 1, 1, fc=c, ec="#555", lw=0.3, label=n) for _, n, c in TISSUES]
    ax.legend(handles=handles, loc="lower left", fontsize=5.2, frameon=True, ncol=2, handlelength=1.0, columnspacing=0.6)
    save(fig, "tfus_model.pdf")


# --------------------------------------------------------------------------------------------- model, mesh steps
# T1 MRI and FSL segmentation of the head model (OpenNeuro ds001734, subject 1, slice 80; ~/Downloads/mri_dataset)
MRI_DATA = os.environ.get("TFUS_MRI_DATA", os.path.expanduser("~/Downloads/mri_dataset/data"))
ZOOM = (138, 164, 196, 222)     # window around the plate (mm), as in model()
ZOOM_MESH = (153.5, 163.5, 199.5, 209.5)   # --mesh-zooms: 10 mm around the lower screw (158.5, 204.5)


def mri_slice():
    """T1 slice of the head model cropped as the meshes (generate_head_model.py --margin 40 --base 4, default
    thresholds), rows = y, columns = x, 1 mm pixels, with its tissue raster; the crop is checked against the tissue
    raster saved with the meshes (head_s80_fit_plate_tissues.npy)."""
    sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src", "mesh_generation", "brain_mri"))
    import generate_head_model as G
    seg = np.load(os.path.join(MRI_DATA, "segmentation", "sub-001_T1w.npy"))[80]
    t1 = G.read_nifti(os.path.join(MRI_DATA, "mri", "sub-001_T1w.nii.gz"))[80]
    tag, (r0, c0) = G.crop(G.tissues(seg, t1, 60.0, 120.0), 40, 4)
    ref = np.load(mesh_path("head_s80_fit_plate_tissues.npy"))
    if tag.shape != ref.shape or not np.array_equal(tag, ref):
        raise RuntimeError("the crop of the MRI slice differs from the tissue raster of the meshes")
    H, W = tag.shape
    img = np.zeros((H, W))
    sr0, sc0 = max(r0, 0), max(c0, 0)
    sr1, sc1 = min(r0 + H, t1.shape[0]), min(c0 + W, t1.shape[1])
    img[sr0 - r0:sr1 - r0, sc0 - c0:sc1 - c0] = t1[sr0:sr1, sc0:sc1]
    return img, tag


def model_steps():
    """tfus_model_steps.pdf (2026-11-15): construction of the mesh of the head model with the plate. Left: T1 slice
    with the tissues in transparency, transducer, its axis and focus, zoom window. Center: zoom, quadtree on the pixels
    (staircase interfaces, head_s80_hf_plate.txt) over the T1 slice. Right: same zoom, cells split along the smoothed
    skull (level 1/2 of the bone indicator smoothed by a Gaussian of 1 pixel, cubic interpolation, as in
    generate_head_model_fitted.py; red) and along the exact plate and screws (head_s80_fit_plate.txt)."""
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import to_rgba
    from scipy import ndimage as ndi
    img, tag = mri_slice()
    H, W = tag.shape
    order = [t for t, _, _ in TISSUES]
    idx = np.vectorize({t: i for i, t in enumerate(order)}.get)(tag)
    cmap = ListedColormap([c for _, _, c in TISSUES])
    lo, hi = np.percentile(img[img > 0], [1, 99.5])
    gray = dict(cmap="gray", vmin=lo, vmax=hi, origin="lower", interpolation="nearest", extent=(0, W, 0, H))
    col = {t: c for t, _, c in TISSUES}
    if MESH_ZOOMS:                    # (2026-11-16) three zooms of 10 mm around the lower screw, same window
        x0, x1, y0, y1 = ZOOM_MESH
        Hf, hgt = 1.95, 0.78
        names = (("head_s80_hf_plate.txt", "staircase"), ("head_s80_fit_plate.txt", "quadtree, split cells"),
                 ("head_s80_vor_plate.txt", "Voronoi, split cells"))
    else:
        x0, x1, y0, y1 = ZOOM
        Hf, hgt = 2.45, 0.86
        names = (("head_s80_hf_plate.txt", None), ("head_s80_fit_plate.txt", None))
    fig = plt.figure(figsize=(6.3, Hf))
    wa = hgt * Hf / 6.3 * W / H
    wz = hgt * Hf / 6.3 * (x1 - x0) / (y1 - y0)
    ax = fig.add_axes([0.0, 0.12, wa, hgt])
    ax.imshow(img, **gray)
    ax.imshow(idx, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5, alpha=0.5, origin="lower", interpolation="nearest",
              extent=(0, W, 0, H))
    transducer(ax, lw=1.2)
    e = np.array([math.cos(AXIS), math.sin(AXIS)])
    ax.plot([FOCUS[0], FOCUS[0] + RADIUS * e[0]], [FOCUS[1], FOCUS[1] + RADIUS * e[1]], ":", color="#d69e2e", lw=0.6)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="#b7202e", lw=0.8))
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    # --paper-labels (2026-11-19): multiplication sign and unit in brackets, as in the preprint
    ax.text(0.03, 0.97, f"{W} × {H} [mm]" if PAPER_LABELS else f"{W} x {H} mm", transform=ax.transAxes, va="top",
            fontsize=5.5, color="white")
    # smoothed skull: level 1/2 of the smoothed bone indicator, cubic interpolation on a fine grid of the window
    phi = ndi.gaussian_filter(np.isin(tag, [1, 2]).astype(float), 1.0) - 0.5
    coef = ndi.spline_filter(phi, order=3)
    xs, ys = np.linspace(x0, x1, 521), np.linspace(y0, y1, 521)
    XX, YY = np.meshgrid(xs, ys)
    PHI = ndi.map_coordinates(coef, [YY.ravel() - 0.5, XX.ravel() - 0.5], order=3, prefilter=False,
                              mode="nearest").reshape(XX.shape)
    for k, (name, title) in enumerate(names):
        axz = fig.add_axes([wa + 0.02 + k * (wz + 0.02), 0.12, wz, hgt])
        axz.imshow(img, **gray)
        pts, cells = read_mesh(mesh_path(name), (x0, x1, y0, y1))
        axz.add_collection(PolyCollection([pts[r] for r, _ in cells], facecolors=[to_rgba(col[t], 0.4) for _, t in cells],
                                          edgecolors=(0, 0, 0, 0.7), linewidths=0.08))
        if MESH_ZOOMS:
            # (2026-11-16) the solid-fluid interface of each mesh (edges between a solid cell, bone or titanium, and a
            # fluid cell), thin, in the three zooms: steps on the staircase, the smoothed skull on the split meshes
            from matplotlib.collections import LineCollection
            owner, segs = {}, []
            for ring, t in cells:
                s_ = t in (1, 2, 8)
                for a, b in zip(ring, ring[1:] + ring[:1]):
                    key = (min(a, b), max(a, b))
                    if key in owner and owner[key] != s_:
                        segs.append([pts[a], pts[b]])
                    owner.setdefault(key, s_)
            axz.add_collection(LineCollection(segs, colors="#b7202e", linewidths=0.45))
        elif k >= 1:
            axz.contour(XX, YY, PHI, levels=[0.0], colors="#b7202e", linewidths=0.7)
        if title:
            axz.set_title(title, fontsize=7, pad=2)
        axz.set_xlim(x0, x1); axz.set_ylim(y0, y1); axz.set_aspect("equal"); axz.set_xticks([]); axz.set_yticks([])
        for s in axz.spines.values():
            s.set_edgecolor("#b7202e"); s.set_linewidth(0.8)
    handles = [Rectangle((0, 0), 1, 1, fc=c, ec="#555", lw=0.3, label=n) for _, n, c in TISSUES]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=8, fontsize=5.5, frameon=False,
               handlelength=1.0, columnspacing=0.8)
    save(fig, "tfus_model_steps.pdf")


def read_mesh(path, box):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    npt, nc = map(int, lines[0].split())
    pts = np.array([list(map(float, l.split())) for l in lines[1:1 + npt]])
    x0, x1, y0, y1 = box
    cells = []
    for l in lines[1 + npt:1 + npt + nc]:
        v = list(map(int, l.split()))
        ring = v[1:1 + v[0]]
        P = pts[ring]
        if P[:, 0].max() < x0 or P[:, 0].min() > x1 or P[:, 1].max() < y0 or P[:, 1].min() > y1:
            continue
        cells.append((ring, v[1 + v[0]]))
    return pts, cells


# ----------------------------------------------------------------------------------------------------- snapshots
def snapshot_times(d):
    T = []
    for f in sorted(os.listdir(d)):
        m = re.match(r"c_snapshot_(\d+)_samples\.dat$", f)
        if m:
            t = float(open(os.path.join(d, f)).readline().split()[3].rstrip(":"))
            T.append((int(m.group(1)), t))
    return sorted(T)


def signed_root(v, lim):
    """sign(v) sqrt(|v| / lim), clipped to [-1, 1]: strong focus and weak transmitted waves on one scale."""
    return np.clip(np.sign(v) * np.sqrt(np.abs(v) / lim), -1, 1)


def snapshots():
    d = run("plate_li")
    if not os.path.isdir(d):
        print("snapshots: missing", d); return
    X0, X1, Y0, Y1, px = 60.0, 225.0, 120.0, 270.0, 0.2
    med, polys = read_cells(os.path.join(d, "c_cells.dat"))
    times = snapshot_times(d)
    if THREE_SNAPSHOTS:
        de = run("plate_li_early")                                    # same mesh, stopped at t = 4 us
        sel = [(de, i, t) for i, t in snapshot_times(de) if t == 4.0] + [(d, i, t) for i, t in times if t in (20.0, 44.0)]
        H = 2.3
        fig, axs = plt.subplots(1, 3, figsize=(6.3, H))
    else:
        sel = [(d, i, t) for i, t in times if t in (20.0, 28.0, 36.0, 44.0, 52.0, 68.0)]
        H = 4.1
        fig, axs = plt.subplots(2, 3, figsize=(6.3, H))
    lim = 1.0
    for ax, (ds, i, t) in zip(np.ravel(axs), sel):
        _, vals = read_samples(os.path.join(ds, f"c_snapshot_{i}_samples.dat"))
        img = signed_root(raster_cells(polys, vals, X0, X1, Y0, Y1, px), lim)
        tissue_underlay(ax, (X0, X1, Y0, Y1), alpha=0.45)
        rgba = plt.get_cmap("RdBu_r")((img + 1) / 2)
        rgba[..., 3] = np.clip(np.abs(img) * 3.0, 0, 1)                 # transparent where the field is small
        im = ax.imshow(rgba, origin="lower", extent=(X0, X1, Y0, Y1), interpolation="bilinear")
        solid_outline(ax, med, polys, (X0, X1, Y0, Y1), lw=0.3)
        transducer(ax, lw=1.0)
        ax.text(0.03, 0.96, t_title(t), transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
        ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    # colour bar: same size in inches for both layouts (4.1 in high with six panels)
    cax = fig.add_axes([0.25, 0.04 * 4.1 / H + (0.02 if THREE_SNAPSHOTS else 0.0), 0.5, 0.013 * 4.1 / H])
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-1, 1), cmap="RdBu_r"), cax=cax, orientation="horizontal")
    ticks = np.array([-1, -0.25, -0.04, 0, 0.04, 0.25, 1])
    cb.set_ticks(signed_root(ticks, lim)); cb.set_ticklabels([minus(f"{t:g}") for t in ticks])
    # with --three-snapshots, the unit amplitude of the burst is stated in the text of the preprint (Source paragraph)
    # --paper-labels (2026-11-18): "mean pressure" -tr(sigma)/2, the sign of the pressure (mean stress = +tr(sigma)/2)
    solid = "mean pressure" if PAPER_LABELS else "mean stress"
    # --paper-labels (2026-11-19): 1/2 tr s with the blackboard s := -sigma of the preprint
    trace = r" $\frac{1}{2}\mathrm{tr}\,\mathbb{s}$" if PAPER_LABELS else r" $-\frac{1}{2}\mathrm{tr}\,\sigma$"
    cb.set_label(r"pressure $p$ (fluids), " + solid + trace + r" (bone, titanium)"
                 + ("" if THREE_SNAPSHOTS else "; burst amplitude 1"), fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=(0.17 if THREE_SNAPSHOTS else 0.09), wspace=0.03,
                        hspace=0.04)
    save(fig, "tfus_snapshots.pdf")


def skull_strip(d):
    """tfus_skull.pdf with --skull-strip (2026-11-16): strips of 28 x 14 mm along the plate, in the frame of the plate
    (principal axis of the titanium cells of head_s80_fit_plate.txt, normal towards the transducer), resampled from an
    axis-aligned raster (nearest pixel of 0.04 mm). Top row: p and the mean stress (signed square root, as before); bottom
    row: shear stress tau_max in the solids, one square-root scale for the three times."""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import PowerNorm
    from scipy import ndimage as ndi
    pts, cells = read_mesh(mesh_path("head_s80_fit_plate.txt"), (130.0, 180.0, 185.0, 235.0))
    ti = [pts[r] for r, tag in cells if tag == 8]
    cent = np.array([q.mean(axis=0) for q in ti])
    area = np.array([abs(0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1])) for q in ti])
    c = (cent * area[:, None]).sum(0) / area.sum()
    t = np.linalg.eigh(((cent - c).T * area) @ (cent - c))[1][:, 1]
    n = np.array([-t[1], t[0]])
    if n @ np.array([math.cos(AXIS), math.sin(AXIS)]) < 0:                # normal towards the transducer (scalp up)
        n = -n
    if t[0] < 0:
        t = -t
    S0, S1, R0, R1, px = -14.0, 14.0, -9.0, 5.0, 0.04
    corn = np.array([c + s * t + r * n for s in (S0, S1) for r in (R0, R1)])
    X0, X1 = math.floor(corn[:, 0].min()) - 1.0, math.ceil(corn[:, 0].max()) + 1.0
    Y0, Y1 = math.floor(corn[:, 1].min()) - 1.0, math.ceil(corn[:, 1].max()) + 1.0
    sc = S0 + (np.arange(int(round((S1 - S0) / px))) + 0.5) * px
    rc = R0 + (np.arange(int(round((R1 - R0) / px))) + 0.5) * px
    SS, RR = np.meshgrid(sc, rc)
    idx = [(c[1] + SS * t[1] + RR * n[1] - Y0) / px - 0.5, (c[0] + SS * t[0] + RR * n[0] - X0) / px - 0.5]
    strip = lambda img: ndi.map_coordinates(img, idx, order=0, mode="nearest")
    med, polys = read_cells(os.path.join(d, "c_cells.dat"))
    solid = strip(raster_cells(polys, np.array(med, float)[:, None, None] * np.ones((1, 1, 1)), X0, X1, Y0, Y1, px, fill=0.0))
    owner, segs = {}, []                                                 # solid outline, in the frame of the plate
    for k, q in enumerate(polys):
        if q[:, 0].max() < X0 or q[:, 0].min() > X1 or q[:, 1].max() < Y0 or q[:, 1].min() > Y1:
            continue
        for e in range(len(q)):
            a, b = tuple(np.round(q[e], 9)), tuple(np.round(q[(e + 1) % len(q)], 9))
            key = (min(a, b), max(a, b))
            if key in owner and owner[key] != med[k]:
                segs.append([((np.array(p_) - c) @ t, (np.array(p_) - c) @ n) for p_ in (a, b)])
            owner.setdefault(key, med[k])
    times = dict(snapshot_times(d))
    sel = [i for i, tt in sorted(times.items()) if tt in (28.0, 36.0, 44.0)]
    fig, axs = plt.subplots(2, 3, figsize=(6.3, 2.3))
    ext = (S0, S1, R0, R1)
    for j, i in enumerate(sel):
        _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_samples.dat"))
        im0 = axs[0, j].imshow(signed_root(strip(raster_cells(polys, vals, X0, X1, Y0, Y1, px)), 0.6), cmap="RdBu_r",
                               vmin=-1, vmax=1, origin="lower", extent=ext, interpolation="bilinear")
        axs[0, j].set_title(t_title(times[i]), fontsize=7, pad=2)
        _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_shear.dat"))
        tau = np.where(solid > 0.5, strip(raster_cells(polys, vals, X0, X1, Y0, Y1, px)), np.nan)
        im1 = axs[1, j].imshow(tau, cmap="magma_r", norm=PowerNorm(0.5, vmin=0.0, vmax=2.0), origin="lower", extent=ext,
                               interpolation="bilinear")
        for a in (axs[0, j], axs[1, j]):
            a.add_collection(LineCollection(segs, colors="#1a1a1a", linewidths=0.5))
            a.set_xlim(S0, S1); a.set_ylim(R0, R1); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    for txt, r in (("scalp", 4.3), ("skull", -4.6)):
        axs[0, 0].text(S0 + 0.6, r, txt, fontsize=5.5, va="center", bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.6))
    fig.subplots_adjust(left=0.01, right=0.9, top=0.92, bottom=0.02, wspace=0.03, hspace=0.06)
    # --paper-labels (2026-11-19): shorter color bars, so that their end labels (-0.6 and 2) do not touch
    cb_kw = dict(shrink=0.85) if PAPER_LABELS else {}
    cb = fig.colorbar(im0, ax=axs[0, :].tolist(), fraction=0.02, pad=0.01, **cb_kw)
    ticks = np.array([-0.6, -0.15, -0.02, 0, 0.02, 0.15, 0.6])
    cb.set_ticks(signed_root(ticks, 0.6)); cb.set_ticklabels([minus(f"{x:g}") for x in ticks])
    cb.set_label(r"$p$, $\frac{1}{2}\mathrm{tr}\,\mathbb{s}$" if PAPER_LABELS else r"$p$, $-\frac{1}{2}\mathrm{tr}\,\sigma$",
                 fontsize=6.5); cb.ax.tick_params(labelsize=5.5)   # --paper-labels (2026-11-19): s := -sigma
    cb = fig.colorbar(im1, ax=axs[1, :].tolist(), fraction=0.02, pad=0.01, **cb_kw)
    cb.set_ticks([0, 0.1, 0.5, 1, 2]); cb.set_ticklabels(["0", "0.1", "0.5", "1", "2"])
    cb.set_label(r"$\tau_{\max}$", fontsize=6.5); cb.ax.tick_params(labelsize=5.5)
    save(fig, "tfus_skull.pdf")


def skull():
    d = run("plate_li")
    if not os.path.isdir(d):
        print("skull: missing", d); return
    if SKULL_STRIP:
        return skull_strip(d)
    X0, X1, Y0, Y1, px = 136.0, 170.0, 192.0, 226.0, 0.04
    med, polys = read_cells(os.path.join(d, "c_cells.dat"))
    solid = raster_cells(polys, np.array(med, float)[:, None, None] * np.ones((1, 1, 1)), X0, X1, Y0, Y1, px, fill=0.0)
    times = dict(snapshot_times(d))
    sel = [i for i, t in sorted(times.items()) if t in (28.0, 36.0, 44.0)]
    fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.25))
    for j, i in enumerate(sel):
        _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_samples.dat"))
        img = raster_cells(polys, vals, X0, X1, Y0, Y1, px)
        ax = axs[0, j]
        im0 = ax.imshow(signed_root(img, 0.6), cmap="RdBu_r", vmin=-1, vmax=1, origin="lower", extent=(X0, X1, Y0, Y1),
                        interpolation="bilinear")
        solid_outline(ax, med, polys, (X0, X1, Y0, Y1), lw=0.5)
        ax.text(0.03, 0.96, rf"$t={times[i]:.0f}\,\mu$s", transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
        _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_shear.dat"))
        img = raster_cells(polys, vals, X0, X1, Y0, Y1, px)
        img = np.where(solid > 0.5, img, np.nan)
        ax = axs[1, j]
        vmax = np.nanpercentile(img, 99)
        im1 = ax.imshow(img, cmap="magma_r", vmin=0, vmax=vmax, origin="lower", extent=(X0, X1, Y0, Y1), interpolation="bilinear")
        solid_outline(ax, med, polys, (X0, X1, Y0, Y1), lw=0.5)
        ax.text(0.03, 0.96, f"max {vmax:.2f}", transform=ax.transAxes, va="top", fontsize=6,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
        for a in (axs[0, j], axs[1, j]):
            a.set_xlim(X0, X1); a.set_ylim(Y0, Y1); a.set_xticks([]); a.set_yticks([])
    cb = fig.colorbar(im0, ax=axs[0, :].tolist(), fraction=0.025, pad=0.01)
    ticks = np.array([-0.6, -0.15, -0.02, 0, 0.02, 0.15, 0.6])
    cb.set_ticks(signed_root(ticks, 0.6)); cb.set_ticklabels([f"{t:g}" for t in ticks])
    cb.set_label(r"$p$, $-\frac{1}{2}\mathrm{tr}\,\sigma$", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    cb = fig.colorbar(im1, ax=axs[1, :].tolist(), fraction=0.025, pad=0.01)
    cb.set_label(r"shear stress $\tau_{\max}$ (scale of each panel)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    cb.set_ticks([])
    save(fig, "tfus_skull.pdf")


# ----------------------------------------------------------------------------------------------------- focusing
def focus():
    runs = [("water_li", "water"), ("intact_li", "intact skull"), ("plate_li", "skull and plate")]
    cells = {"water_li": "intact_li", "intact_li": "intact_li", "plate_li": "plate_li"}   # water: the intact mesh
    if not all(os.path.exists(os.path.join(run(r), "c_peak_samples.dat")) for r, _ in runs):
        print("focus: missing peak maps"); return
    e = np.array([math.cos(AXIS), math.sin(AXIS)]); n = np.array([-e[1], e[0]])
    X0, X1, Y0, Y1, px = 85.0, 175.0, 140.0, 230.0, 0.25
    P, geo = {}, {}
    for r, _ in runs:
        med, polys = read_cells(os.path.join(run(cells[r]), "c_cells.dat"))
        _, vals = read_samples(os.path.join(run(r), "c_peak_samples.dat"))
        P[r] = raster_cells(polys, vals, X0, X1, Y0, Y1, px)
        geo[r] = (med, polys)
    tis = np.load(mesh_path("head_s80_fit_plate_tissues.npy"))
    xc = X0 + (np.arange(P["water_li"].shape[1]) + 0.5) * px
    yc = Y0 + (np.arange(P["water_li"].shape[0]) + 0.5) * px
    XX, YY = np.meshgrid(xc, yc)
    brain = np.isin(tis[np.clip(YY.astype(int), 0, tis.shape[0] - 1), np.clip(XX.astype(int), 0, tis.shape[1] - 1)], [6, 7])
    ref = np.nanmax(np.where(brain, P["water_li"], 0))
    kw = np.unravel_index(np.nanargmax(np.where(brain, P["water_li"], 0)), XX.shape)
    fw = np.array([XX[kw], YY[kw]])
    if FOCUS_ROW:
        # one row (2026-11-16): the three maps (side s inches) and a vertical colour bar; the titles give the peak in
        # the brain and the peak at the geometric focus (sensor "focus"), both in dB re the focal peak in water
        W, s, g, x0, y0 = 6.3, 1.88, 0.06, 0.02, 0.03
        H = y0 + s + 0.34
        fig = plt.figure(figsize=(W, H))
        row_axes = [fig.add_axes([(x0 + j * (s + g)) / W, y0 / H, s / W, s / H]) for j in range(3)]
    else:
        fig = plt.figure(figsize=(6.3, 4.2))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.3, wspace=0.08)
    stats = {}
    for j, (r, lab) in enumerate(runs):
        ax = row_axes[j] if FOCUS_ROW else fig.add_subplot(gs[0, j])
        db = 20 * np.log10(np.maximum(P[r], 1e-12) / ref)
        im = ax.imshow(db, cmap="inferno", vmin=-18, vmax=0, origin="lower", extent=(X0, X1, Y0, Y1), interpolation="bilinear")
        if r != "water_li":
            solid_outline(ax, *geo[r], (X0, X1, Y0, Y1), color="white", lw=0.35)
        Pb = np.where(brain, P[r], 0)
        k = np.unravel_index(np.nanargmax(Pb), Pb.shape)
        f = np.array([XX[k], YY[k]]); dd = f - fw
        lab_, _ = ndi_label(Pb >= Pb[k] / math.sqrt(2))
        stats[r] = dict(peak=Pb[k] / ref, x=f[0], y=f[1], shift=float(np.hypot(*dd)), axial=float(dd @ e),
                        lateral=float(dd @ n), area=float(np.sum(lab_ == lab_[k]) * px * px))
        ax.plot(*FOCUS, "+", color="#4fd1c5", ms=6, mew=1.0)
        ax.plot(*f, "x", color="white", ms=4, mew=0.8)
        db = 20 * math.log10(stats[r]['peak'])
        ax.set_title(f"{lab}: {minus(f'{db:+.1f}')}{DB}", fontsize=7, pad=2)
        ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([])
    if FOCUS_ROW:
        for j, (r, lab) in enumerate(runs):
            f = os.path.join(run(r), "c_focus.csv")
            tgt = np.max(np.abs(csv(f)[:, 1])) / ref if os.path.exists(f) else float("nan")
            db, db_t = 20 * math.log10(stats[r]['peak']), 20 * math.log10(tgt)
            row_axes[j].set_title(f"{lab}\nmax {minus(f'{db:+.1f}')}{DB}, "
                                  f"target {minus(f'{db_t:+.1f}')}{DB}", fontsize=7, pad=2, linespacing=1.3)
        cax = fig.add_axes([(x0 + 3 * s + 2 * g + 0.06) / W, y0 / H, 0.09 / W, s / H])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("peak pressure [dB]" if PAPER_LABELS else "peak pressure (dB)", fontsize=6.5, labelpad=2)
        cb.ax.tick_params(labelsize=6)
        save(fig, "tfus_focus.pdf")
    else:
        focus_profiles(fig, gs, im, runs, P, ref, fw, e, n, X0, Y0, px)
    with open(os.path.join(OUT, "tfus_focus.txt"), "w") as fo:
        fo.write(f"water focus at ({fw[0]:.2f}, {fw[1]:.2f}), {np.hypot(*(fw - FOCUS)):.2f} mm from the geometric focus\n")
        for r, lab in runs:
            st = stats[r]
            fo.write(f"{lab}: peak in the brain {st['peak']:.3f} ({20 * math.log10(st['peak']):+.2f} dB) at ({st['x']:.2f}, "
                     f"{st['y']:.2f}), shift {st['shift']:.2f} mm (axial {st['axial']:+.2f}, lateral {st['lateral']:+.2f}), "
                     f"-3 dB area {st['area']:.1f} mm^2\n")
    print(open(os.path.join(OUT, "tfus_focus.txt")).read())


def focus_profiles(fig, gs, im, runs, P, ref, fw, e, n, X0, Y0, px):
    """Default layout of tfus_focus.pdf (before 2026-11-16): colour bar at the right of the maps, profiles of the peak
    pressure along the axis (from the water focus) and across it, in the second row."""
    cb = fig.colorbar(im, ax=fig.axes[:3], fraction=0.02, pad=0.01)
    cb.set_label("peak pressure (dB re focus in water)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    for j, (dirv, xl, ss) in enumerate(((e, "along the axis, from the water focus (mm)", np.linspace(-40, 30, 561)),
                                         (n, "across the axis (mm)", np.linspace(-15, 15, 241)))):
        ax = fig.add_subplot(gs[1, 0] if j == 0 else gs[1, 1:3])
        for (r, lab), c in zip(runs, ("#2b6cb0", "#2f855a", "#b7202e")):
            pts = fw[None, :] + ss[:, None] * dirv[None, :]
            ii = np.clip(((pts[:, 1] - Y0) / px).astype(int), 0, P[r].shape[0] - 1)
            jj = np.clip(((pts[:, 0] - X0) / px).astype(int), 0, P[r].shape[1] - 1)
            ax.plot(ss, P[r][ii, jj] / ref, color=c, lw=0.9, label=lab)
        ax.set_xlabel(xl, labelpad=1); ax.grid(True, **GRID)
        if j == 0:
            ax.set_ylabel("peak pressure / water focus", labelpad=1)
            ax.legend(loc="upper left", frameon=False)
    save(fig, "tfus_focus.pdf")


def ndi_label(mask):
    from scipy import ndimage
    return ndimage.label(mask)


# ------------------------------------------------------------------------------------------------- signals, cost
def csv(path):
    L = [l for l in open(path) if not l.startswith("#")]
    return np.array([[float(x) for x in l.split(",")] for l in L[1:] if l.strip()])


def wall(d):
    f = os.path.join(d, "time.txt")
    return float(open(f).read().split()[0]) if os.path.exists(f) else None


def best_integration(d):
    """Best integration wall time of d and of its repetitions d_r2, d_r3 (same job, run alone)."""
    v = [integration(x) for x in (d, d + "_r2", d + "_r3")]
    v = [x for x in v if x is not None]
    return min(v) if v else None


def integration(d):
    f = os.path.join(d, "log.txt")
    if not os.path.exists(f):
        return None
    m = re.search(r"integration wall time \(s\): ([0-9.eE+-]+)", re.sub(r"\x1b\[[0-9;]*m", "", open(f).read()))
    return float(m.group(1)) if m else None


def energy(d):
    f = os.path.join(d, "c_energy.dat")
    if not os.path.exists(f):
        return None
    E = np.loadtxt(f, comments="#")
    return E[:, 0], E[:, 1] / E[0, 1]


def signals_row():
    """tfus_signals.pdf with --signals-row (2026-11-16): accuracy check with the plate (LI at dt and dt/2, MLTS) and
    discrete energy; the costs are in Table tfus-cost, the focal peaks in Fig. tfus-focus."""
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.0))
    ax = axs[0]
    mlts_run = "plate_mlts_w" if os.path.isdir(run("plate_mlts_w")) else "plate_mlts"
    for r, lab, c, ls in (("plate_li", r"LI, $\Delta t$", SCH["li"][1], "-"), ("plate_li_dt2", r"LI, $\Delta t/2$", "#90cdf4", "--"),
                          (mlts_run, "MLTS", SCH["mlts"][1], ":")):
        f = os.path.join(run(r), "c_focus.csv")
        if os.path.exists(f):
            A = csv(f); ax.plot(A[:, 0], A[:, 1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(T_LABEL, labelpad=1); ax.set_ylabel("pressure at the geometric focus", labelpad=1)
    ax.set_xlim(SIGNALS_XMIN, 56); ax.grid(True, **GRID)
    ax.legend(loc="upper right", frameon=False, fontsize=5.5)
    if not PAPER_LABELS:
        ax.set_title("skull and plate", fontsize=7, pad=2)
    ax = axs[1]
    # --paper-labels (2026-11-19): configurations in colors outside the scheme colors (MLTS red, LI blue in the paper)
    c_water, c_plate = ("#4a4a4a", "#805ad5") if PAPER_LABELS else ("#2b6cb0", "#b7202e")
    for r, kind, lab, c, ls in (("water_li", "fit", "water", c_water, "-"), ("intact_li", "fit", "intact, fitted", "#2f855a", "-"),
                                ("plate_li", "fit", "plate, fitted", c_plate, "-"),
                                ("intact_mlts", "tfus", "intact, staircase", "#2f855a", "--"),
                                ("plate_mlts", "tfus", "plate, staircase", c_plate, "--")):
        E = energy(first_run(kind, r, r.replace("_mlts", "_li")))
        if E is not None:
            ax.plot(E[0], E[1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(T_LABEL, labelpad=1); ax.set_ylabel(r"$E(t)/E(0)$", labelpad=1); ax.grid(True, **GRID)
    if ENERGY_YMIN is not None:
        ax.set_ylim(ENERGY_YMIN, 1.005)
    ax.legend(loc="lower left", frameon=False, fontsize=5.5)
    if not PAPER_LABELS:
        ax.set_title("discrete energy", fontsize=7, pad=2)
    fig.tight_layout(pad=0.2, w_pad=1.5)
    save(fig, "tfus_signals.pdf")


def signals():
    if SIGNALS_ROW:
        return signals_row()
    # top: signals and energies; bottom: costs of Table tfus-cost (times, peak memories)
    fig = plt.figure(figsize=(6.3, 4.2))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 0.9])
    axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
           fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
    ax = axs[0]
    for (r, lab), c in zip((("water_li", "water"), ("intact_li", "intact skull"), ("plate_li", "skull and plate")),
                           ("#2b6cb0", "#2f855a", "#b7202e")):
        f = os.path.join(run(r), "c_focus.csv")
        if os.path.exists(f):
            A = csv(f); ax.plot(A[:, 0], A[:, 1], color=c, lw=0.7, label=lab)
    ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_ylabel("pressure at the focus", labelpad=1)
    ax.set_xlim(30, 80); ax.grid(True, **GRID); ax.legend(loc="lower left", frameon=False, fontsize=5.2)
    ax = axs[1]
    mlts_run = "plate_mlts_w" if os.path.isdir(run("plate_mlts_w")) else "plate_mlts"
    for r, lab, c, ls in (("plate_li", r"LI, $\Delta t$", SCH["li"][1], "-"), ("plate_li_dt2", r"LI, $\Delta t/2$", "#90cdf4", "--"),
                          (mlts_run, "MLTS", SCH["mlts"][1], ":")):
        f = os.path.join(run(r), "c_focus.csv")
        if os.path.exists(f):
            A = csv(f); ax.plot(A[:, 0], A[:, 1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_xlim(44, 56); ax.grid(True, **GRID)
    ax.legend(loc="lower left", frameon=False, fontsize=5.2); ax.set_title("skull and plate", fontsize=7, pad=2)
    ax = axs[2]
    for r, kind, lab, c, ls in (("water_li", "fit", "water", "#2b6cb0", "-"), ("intact_li", "fit", "intact, fitted", "#2f855a", "-"),
                                ("plate_li", "fit", "plate, fitted", "#b7202e", "-"),
                                ("intact_mlts", "tfus", "intact, staircase", "#2f855a", "--"),
                                ("plate_mlts", "tfus", "plate, staircase", "#b7202e", "--")):
        E = energy(run(r, kind))
        if E is not None:
            ax.plot(E[0], E[1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_ylabel(r"$E(t)/E(0)$", labelpad=1); ax.grid(True, **GRID)
    ax.legend(loc="lower left", frameon=False, fontsize=5.0); ax.set_title("discrete energy", fontsize=7, pad=2)
    T, M = cost_data()
    names = [s_ for s_ in ("erk", "mlts", "li") if T.get((s_, 6)) or T.get((s_, 1))]
    for ax, D, yl in ((axs[3], T, "time integration (min)"), (axs[4], M, f"{MEMORY_LABEL} (GiB)")):
        for i_, s_ in enumerate(names):
            for j_, th in enumerate((1, 6)):
                v = D.get((s_, th))
                if not v:
                    continue
                xi = i_ + (j_ - 0.5) * 0.38
                ax.bar(xi, v, width=0.36, color=SCH[s_][1], edgecolor="#1a1a1a", linewidth=0.4,
                       hatch="////" if th == 1 else None)
                lab = f"{v:.0f}" if D is T else f"{v:.1f}"
                if s_ != "erk" and D.get(("erk", th)):
                    r_ = D[("erk", th)] / v if D is T else v / D[("erk", th)]
                    lab += f"\n($\\times${r_:.1f})" if D is T else f"\n($\\times${r_:.2f})"
                ax.text(xi, v * (1.08 if D is T else 1.02), lab, ha="center", va="bottom", fontsize=5.2)
        ax.set_xticks(range(len(names))); ax.set_xticklabels([SCH[s_][0] for s_ in names], fontsize=6.5)
        ax.set_ylabel(yl, labelpad=1); ax.grid(True, axis="y", **GRID)
        vals = [v for v in D.values() if v]
        if D is T:
            ax.set_yscale("log"); ax.set_ylim(min(vals) / 2, max(vals) * 8)
        else:
            ax.set_ylim(0, max(vals) * 1.4)
    axs[3].legend(handles=[Rectangle((0, 0), 1, 1, fc="white", ec="#1a1a1a", lw=0.4, hatch="////", label="1 thread"),
                           Rectangle((0, 0), 1, 1, fc="white", ec="#1a1a1a", lw=0.4, label="6 threads")],
                  loc="upper right", fontsize=5.5, frameon=False)
    axs[3].set_title("time integration (best of three)", fontsize=7, pad=2)
    axs[4].set_title(MEMORY_LABEL, fontsize=7, pad=2)
    fig.tight_layout(pad=0.2, w_pad=0.6, h_pad=1.0)
    save(fig, "tfus_signals.pdf")


def polytopal():
    """Voronoi (polytopal) mesh of the head with the plate: mesh around the plate, snapshot, pressure at the focus against
    the quadtree mesh (runs vor/plate_li and fit/plate_li; fine campaign: vor/plate_mlts)."""
    d = first_run("vor", "plate_li", "plate_mlts")
    if not os.path.isdir(d):
        print("polytopal: missing", d); return
    from matplotlib.collections import PolyCollection
    fig = plt.figure(figsize=(6.3, 2.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 1.1], wspace=0.32)
    ax = fig.add_subplot(gs[0, 0])
    pts, cells = read_mesh(mesh_path("head_s80_vor_plate.txt"), (144, 160, 202, 218))
    col = {t: c for t, _, c in TISSUES}
    ax.add_collection(PolyCollection([pts[r] for r, _ in cells], facecolors=[col[t] for _, t in cells], edgecolors="k",
                                     linewidths=0.12))
    ax.set_xlim(144, 160); ax.set_ylim(202, 218); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if not PAPER_LABELS:
        ax.set_title("Voronoi cells around the plate", fontsize=6.5, pad=2)
    ax = fig.add_subplot(gs[0, 1])
    X0, X1, Y0, Y1, px = 60.0, 225.0, 120.0, 270.0, 0.25
    med, polys = read_cells(os.path.join(d, "c_cells.dat"))
    times = dict(snapshot_times(d))
    i = [k for k, t in times.items() if t == 44.0][0]
    _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_samples.dat"))
    img = signed_root(raster_cells(polys, vals, X0, X1, Y0, Y1, px), 1.0)
    tissue_underlay(ax, (X0, X1, Y0, Y1), alpha=0.45)
    rgba = plt.get_cmap("RdBu_r")((img + 1) / 2)
    rgba[..., 3] = np.clip(np.abs(img) * 3.0, 0, 1)
    ax.imshow(rgba, origin="lower", extent=(X0, X1, Y0, Y1), interpolation="bilinear")
    solid_outline(ax, med, polys, (X0, X1, Y0, Y1), lw=0.3)
    transducer(ax, lw=1.0)
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if not PAPER_LABELS:
        ax.set_title(r"$t=44\,\mu$s, Voronoi mesh", fontsize=6.5, pad=2)
    ax = fig.add_subplot(gs[0, 2])
    # --paper-labels (2026-11-19): colors of the schemes as in the other figures, MLTS (Voronoi run) red, LI (quadtree) blue
    c_quad, c_vor = ("#2b6cb0", "#b7202e") if PAPER_LABELS else ("#b7202e", "#2b6cb0")
    for r, kind, lab, c, ls in (("plate_li", "fit", "quadtree", c_quad, "-"),
                                ("plate_li", "vor", "Voronoi", c_vor, "--")):
        f = os.path.join(first_run(kind, r, "plate_mlts") if kind == "vor" else run(r, kind), "c_focus.csv")
        if os.path.exists(f):
            A = csv(f); ax.plot(A[:, 0], A[:, 1], ls, color=c, lw=0.8, label=lab)
    if PAPER_LABELS:   # smaller fonts than the defaults (labels 7, ticks 6), closer to the small panels on the left
        ax.set_xlabel(T_LABEL, labelpad=1, fontsize=7); ax.set_ylabel("pressure at the geometric focus", labelpad=1,
                                                                     fontsize=7)
        ax.tick_params(labelsize=6)
    else:
        ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_ylabel("pressure at the focus", labelpad=1)
    # --paper-labels (2026-11-19): legend at the upper right, off the dip of the signal near 45-46 us
    ax.set_xlim(40, 60); ax.grid(True, **GRID)
    ax.legend(loc="upper right" if PAPER_LABELS else "lower left", frameon=False, fontsize=5.5)
    save(fig, "tfus_polytopal.pdf")


def cost_data():
    """Costs of Table tfus-cost, best of three: integration times (min) T[(scheme, threads)] and peak memories (GiB)
    M[(scheme, threads)]. ERK4 over 8 us (10150 identical steps) x 10; on one thread, MLTS and LI over 8 us (139 of the
    1394 steps) extrapolated, the factorization of LI ("fine system") counted once."""
    def fact(d):
        m = re.search(r"fine system ([0-9.eE+-]+)", re.sub(r"\x1b\[[0-9;]*m", "", open(os.path.join(d, "log.txt")).read()))
        return float(m.group(1))
    def best(d, f):
        v = [f(x) for x in (d, d + "_r2", d + "_r3") if integration(x) is not None]
        return min(v) / 60 if v else None
    mlts_run = "plate_mlts_w" if os.path.isdir(run("plate_mlts_w")) else "plate_mlts"
    li_run, erk6, erk1 = "plate_li", "plate_erk_8us", "erk_8us"
    T = {("erk", 6): best(run("plate_erk_8us"), lambda x: 10 * integration(x)),
         ("mlts", 6): best(run(mlts_run), integration), ("li", 6): best(run("plate_li"), integration),
         ("erk", 1): best(run("erk_8us", "t1"), lambda x: 10 * integration(x)),
         ("mlts", 1): best(run("mlts_8us", "t1"), lambda x: integration(x) * 1394 / 139),
         ("li", 1): best(run("li_8us", "t1"), lambda x: (integration(x) - fact(x)) * 1394 / 139 + fact(x))}
    if SHORT:
        mlts_run, li_run = "plate_mlts_8us", "plate_li_8us"
        T[("mlts", 6)] = best(run(mlts_run), lambda x: integration(x) * 1394 / 139)
        T[("li", 6)] = best(run(li_run), lambda x: (integration(x) - fact(x)) * 1394 / 139 + fact(x))
    if FINE_TIMING:
        # --fine-timing (2026-11-16, meshes twice finer, lts_work/head_fine/chain_T.sh): ERK(4) over 2 us x 40, MLTS and
        # LI over the first 278 of the 2788 macro steps, the factorization of LI counted once, on six and one threads
        r_ = 2788 / 278
        mlts_run, li_run, erk6, erk1 = "plate_mlts_8us", "plate_li_8us", "plate_erk_2us", "erk_2us"
        T = {("erk", 6): best(run(erk6), lambda x: 40 * integration(x)),
             ("mlts", 6): best(run(mlts_run), lambda x: integration(x) * r_),
             ("li", 6): best(run(li_run), lambda x: (integration(x) - fact(x)) * r_ + fact(x)),
             ("erk", 1): best(run(erk1, "t1"), lambda x: 40 * integration(x)),
             ("mlts", 1): best(run("mlts_8us", "t1"), lambda x: integration(x) * r_),
             ("li", 1): best(run("li_8us", "t1"), lambda x: (integration(x) - fact(x)) * r_ + fact(x))}
    M = {(s_, th): os.path.join(d, "time.txt") for s_, th, d in
         (("erk", 6, run(erk6)), ("mlts", 6, run(mlts_run)), ("li", 6, run(li_run)),
          ("erk", 1, run(erk1, "t1")), ("mlts", 1, run("mlts_8us", "t1")), ("li", 1, run("li_8us", "t1")))}
    if MEM:   # algebraic memory of the same runs (line "algebraic memory: ... total X MB" of their log.txt), GiB
        M = {k: os.path.join(MEM, os.path.relpath(os.path.dirname(f), RUNS), "log.txt") for k, f in M.items()}
        M = {k: float(re.search(r"algebraic memory: .* total ([0-9.]+) MB", open(f).read()).group(1)) / 1024
             for k, f in M.items() if os.path.exists(f) and "algebraic memory" in open(f).read()}
        return T, M
    M = {k: float(open(f).read().split()[2]) / 1048576 for k, f in M.items() if os.path.exists(f)}
    return T, M


def cost_table():
    """Table tfus-cost (cost_data), memory of the runs on six threads."""
    T, M = cost_data()
    M = {s_: M[(s_, 6)] for s_ in ("erk", "mlts", "li") if (s_, 6) in M}
    with open(os.path.join(OUT, "tfus_cost.txt"), "w") as fo:
        fo.write("scheme | 6 threads [min] gain | 1 thread [min] gain | memory [GiB] ratio\n")
        for s_ in ("erk", "mlts", "li"):
            g = lambda th: T[("erk", th)] / T[(s_, th)] if T[(s_, th)] and T[("erk", th)] else float("nan")
            fo.write(f"{s_:5} | {T[(s_, 6)] or float('nan'):6.1f} {g(6):5.2f} | {T[(s_, 1)] or float('nan'):6.1f} {g(1):5.2f} | "
                     f"{M.get(s_, float('nan')):4.1f} {M.get(s_, float('nan')) / M.get('erk', float('nan')):5.2f}\n")
        if FINE_TIMING:
            # Voronoi mesh, MLTS on six threads (vor/plate_mlts_8us: first kv of the nv steps of the full run, case.lua)
            d = run("plate_mlts_8us", "vor")
            if os.path.exists(os.path.join(d, "case.lua")):
                n_ = [int(x) for x in re.findall(r"time\.n_steps = ([0-9]+)", open(os.path.join(d, "case.lua")).read())]
                v = [integration(x) for x in (d, d + "_r2", d + "_r3") if integration(x) is not None]
                if v and len(n_) >= 2:
                    fo.write(f"voronoi mlts | {min(v) * n_[0] / n_[-1] / 60:6.1f} min on 6 threads ({n_[-1]} of {n_[0]} "
                             f"steps, best of {len(v)})\n")
    print(open(os.path.join(OUT, "tfus_cost.txt")).read())


if __name__ == "__main__":
    for f in (model, model_steps, snapshots, skull, focus, signals, polytopal, cost_table):
        try:
            f()
        except Exception as ex:
            import traceback; traceback.print_exc()
            print(f.__name__, "failed:", repr(ex))
