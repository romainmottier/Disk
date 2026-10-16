#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures of Section 6.7 of the LTS_HHO preprint (main.tex): transcranial focused ultrasound at 250 kHz in the
seven-tissue head model with and without a titanium fixation plate (inputs article/inputs/tfus_fit.lua (fitted solid
interfaces) and tfus.lua (voxel staircase), jobs fit/* and tfus/* of article/jobs.txt). Reads the runs of the directory
given as first argument (fit/plate_li, intact_li, water_li, plate_li_dt2, plate_mlts; tfus/plate_mlts, intact_mlts,
water_mlts) and writes into the second one:
    tfus_model.pdf      tissues, transducer, focus, beam axis; zoom on the plate
    tfus_snapshots.pdf  pressure in the fluids and mean stress -tr(sigma)/2 in the solids, one colour scale
    tfus_skull.pdf      zoom on the skull and the plate: mean stress (P waves) and shear stress (S waves)
    tfus_focus.pdf      peak pressure maps (water, intact skull, plate) and profiles along and across the axis
    tfus_signals.pdf    pressure at the focus (three models; schemes and steps with the plate), discrete energy (fitted
                        and staircase interfaces), costs
Usage: python3 make_tfus.py <runs> <figures>
"""
import math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

RUNS = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
HERE = os.path.dirname(os.path.abspath(__file__))
MESHES = os.path.join(HERE, "..", "..", "..", "src", "mesh_generation", "brain_mri", "meshes")
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
    tis = np.load(os.path.join(MESHES, "head_s80_fit_plate_tissues.npy"))
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
    tis = np.load(os.path.join(MESHES, "head_s80_fit_plate_tissues.npy"))
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
    pts, cells = read_mesh(os.path.join(MESHES, "head_s80_fit_plate.txt"), (138, 164, 196, 222))
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
    sel = [i for i, t in times if t in (20.0, 28.0, 36.0, 44.0, 52.0, 68.0)]
    fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.1))
    lim = 1.0
    for ax, i in zip(axs.flat, sel):
        _, vals = read_samples(os.path.join(d, f"c_snapshot_{i}_samples.dat"))
        img = signed_root(raster_cells(polys, vals, X0, X1, Y0, Y1, px), lim)
        tissue_underlay(ax, (X0, X1, Y0, Y1), alpha=0.45)
        rgba = plt.get_cmap("RdBu_r")((img + 1) / 2)
        rgba[..., 3] = np.clip(np.abs(img) * 3.0, 0, 1)                 # transparent where the field is small
        im = ax.imshow(rgba, origin="lower", extent=(X0, X1, Y0, Y1), interpolation="bilinear")
        solid_outline(ax, med, polys, (X0, X1, Y0, Y1), lw=0.3)
        transducer(ax, lw=1.0)
        t = dict(times)[i]
        ax.text(0.03, 0.96, rf"$t={t:.0f}\,\mu$s", transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
        ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    cax = fig.add_axes([0.25, 0.04, 0.5, 0.013])
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-1, 1), cmap="RdBu_r"), cax=cax, orientation="horizontal")
    ticks = np.array([-1, -0.25, -0.04, 0, 0.04, 0.25, 1])
    cb.set_ticks(signed_root(ticks, lim)); cb.set_ticklabels([f"{t:g}" for t in ticks])
    cb.set_label(r"pressure $p$ (fluids), mean stress $-\frac{1}{2}\mathrm{tr}\,\sigma$ (bone, titanium); burst amplitude 1",
                 fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.09, wspace=0.03, hspace=0.04)
    save(fig, "tfus_snapshots.pdf")


def skull():
    d = run("plate_li")
    if not os.path.isdir(d):
        print("skull: missing", d); return
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
    tis = np.load(os.path.join(MESHES, "head_s80_fit_plate_tissues.npy"))
    xc = X0 + (np.arange(P["water_li"].shape[1]) + 0.5) * px
    yc = Y0 + (np.arange(P["water_li"].shape[0]) + 0.5) * px
    XX, YY = np.meshgrid(xc, yc)
    brain = np.isin(tis[np.clip(YY.astype(int), 0, tis.shape[0] - 1), np.clip(XX.astype(int), 0, tis.shape[1] - 1)], [6, 7])
    ref = np.nanmax(np.where(brain, P["water_li"], 0))
    kw = np.unravel_index(np.nanargmax(np.where(brain, P["water_li"], 0)), XX.shape)
    fw = np.array([XX[kw], YY[kw]])
    fig = plt.figure(figsize=(6.3, 4.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.3, wspace=0.08)
    stats = {}
    for j, (r, lab) in enumerate(runs):
        ax = fig.add_subplot(gs[0, j])
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
        ax.set_title(f"{lab}: {20 * math.log10(stats[r]['peak']):+.1f} dB", fontsize=7, pad=2)
        ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([])
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
    with open(os.path.join(OUT, "tfus_focus.txt"), "w") as fo:
        fo.write(f"water focus at ({fw[0]:.2f}, {fw[1]:.2f}), {np.hypot(*(fw - FOCUS)):.2f} mm from the geometric focus\n")
        for r, lab in runs:
            st = stats[r]
            fo.write(f"{lab}: peak in the brain {st['peak']:.3f} ({20 * math.log10(st['peak']):+.2f} dB) at ({st['x']:.2f}, "
                     f"{st['y']:.2f}), shift {st['shift']:.2f} mm (axial {st['axial']:+.2f}, lateral {st['lateral']:+.2f}), "
                     f"-3 dB area {st['area']:.1f} mm^2\n")
    print(open(os.path.join(OUT, "tfus_focus.txt")).read())


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


def signals():
    fig, axs = plt.subplots(1, 4, figsize=(6.3, 2.15), gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.0, 0.8]))
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
                                ("intact_mlts", "tfus", "intact, voxels", "#2f855a", "--"),
                                ("plate_mlts", "tfus", "plate, voxels", "#b7202e", "--")):
        E = energy(run(r, kind))
        if E is not None:
            ax.plot(E[0], E[1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_ylabel(r"$E(t)/E(0)$", labelpad=1); ax.grid(True, **GRID)
    ax.legend(loc="lower left", frameon=False, fontsize=5.0); ax.set_title("discrete energy", fontsize=7, pad=2)
    ax = axs[3]
    # costs on six threads: ERK4 measured over 8 us and extrapolated to 80 us (its steps are all alike), MLTS, LI; memory
    def mem(d):
        f = os.path.join(d, "time.txt")
        return float(open(f).read().split()[2]) / 1048576 if os.path.exists(f) else None
    t, m = {}, {}
    e8 = integration(run("plate_erk_8us"))
    if e8:
        t["erk"] = e8 * 80.0 / 8.0; m["erk"] = mem(run("plate_erk_8us"))
    for s_, d in (("mlts", mlts_run), ("li", "plate_li")):
        if integration(run(d)):
            t[s_] = integration(run(d)); m[s_] = mem(run(d))
    names = [s_ for s_ in ("erk", "mlts", "li") if s_ in t]
    ax.bar(range(len(names)), [t[s_] / 60 for s_ in names], color=[SCH[s_][1] for s_ in names], edgecolor="#1a1a1a", linewidth=0.4)
    for i_, s_ in enumerate(names):
        lab = f"{t[s_] / 60:.0f} min"
        if s_ != "erk" and "erk" in t:
            lab += f"\n($\\times${t['erk'] / t[s_]:.1f})"
        if m.get(s_):
            lab += f"\n{m[s_]:.1f} GiB"
        ax.text(i_, t[s_] / 60 * 1.02, lab, ha="center", va="bottom", fontsize=5.2)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([SCH[s_][0] for s_ in names], fontsize=6.5)
    ax.set_ylabel("time integration (min), 6 threads", labelpad=1)
    if names:
        ax.set_ylim(0, max(t[s_] / 60 for s_ in names) * 1.45)
    ax.grid(True, axis="y", **GRID)
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "tfus_signals.pdf")


def polytopal():
    """Voronoi (polytopal) mesh of the head with the plate: mesh around the plate, snapshot, pressure at the focus against
    the quadtree mesh (runs vor/plate_li and fit/plate_li)."""
    d = run("plate_li", "vor")
    if not os.path.isdir(d):
        print("polytopal: missing", d); return
    from matplotlib.collections import PolyCollection
    fig = plt.figure(figsize=(6.3, 2.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 1.1], wspace=0.32)
    ax = fig.add_subplot(gs[0, 0])
    pts, cells = read_mesh(os.path.join(MESHES, "head_s80_vor_plate.txt"), (144, 160, 202, 218))
    col = {t: c for t, _, c in TISSUES}
    ax.add_collection(PolyCollection([pts[r] for r, _ in cells], facecolors=[col[t] for _, t in cells], edgecolors="k",
                                     linewidths=0.12))
    ax.set_xlim(144, 160); ax.set_ylim(202, 218); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
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
    ax.set_title(r"$t=44\,\mu$s, Voronoi mesh", fontsize=6.5, pad=2)
    ax = fig.add_subplot(gs[0, 2])
    for r, kind, lab, c, ls in (("plate_li", "fit", "quadtree, cut cells", "#b7202e", "-"),
                                ("plate_li", "vor", "Voronoi, cut cells", "#2b6cb0", "--")):
        f = os.path.join(run(r, kind), "c_focus.csv")
        if os.path.exists(f):
            A = csv(f); ax.plot(A[:, 0], A[:, 1], ls, color=c, lw=0.8, label=lab)
    ax.set_xlabel(r"$t$ ($\mu$s)", labelpad=1); ax.set_ylabel("pressure at the focus", labelpad=1)
    ax.set_xlim(40, 60); ax.grid(True, **GRID); ax.legend(loc="lower left", frameon=False, fontsize=5.5)
    save(fig, "tfus_polytopal.pdf")


if __name__ == "__main__":
    for f in (model, snapshots, skull, focus, signals, polytopal):
        try:
            f()
        except Exception as ex:
            import traceback; traceback.print_exc()
            print(f.__name__, "failed:", repr(ex))
