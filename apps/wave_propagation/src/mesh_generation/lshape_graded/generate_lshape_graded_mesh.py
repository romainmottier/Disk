#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Graded quadtree-refined mesh for the L-shaped domain with a reentrant
corner, in DiSk++'s "poly_2d" text format (read by polygon_2d_mesh_reader).
===============================================================================
Domain: Omega = (0,1)^2 \\ [0.5,1) x (0.5,1], reentrant corner at
H = (0.5, 0.5).

IMPORTANT, learned the hard way (see git history / session notes): refining
ONLY a fixed-size neighbourhood of the corner while leaving the background
mesh untouched does NOT test the classical graded-mesh h->0 limit at all.
With a level-independent background, the measured L2 error converges at the
level-independent, k-INDEPENDENT rate O(h_corner^{2*gamma}) (gamma=2/3 for
this 270 deg corner, i.e. order ~1.33), which is exactly the un-graded
"quasi-uniform mesh" rate from classical FEM theory (Ciarlet-Raviart type
estimate for a solution with u in H^{1+gamma}) -- confirmed experimentally
by getting an (uncorrelated with k) order ~1.3 at BOTH k=1 and k=3 on such
a mesh. The corner refinement was doing nothing for the achievable order
because the background error dominates and never shrinks.

The correct construction couples TWO refinements per convergence level ell:
  - UNIFORM refinement of the whole background mesh, ell times (halving the
    far-field cell size each time: h_bg(ell) = h0 / 2^ell) -- this is the
    part that must shrink for the classical h->0 theory to apply at all.
  - EXTRA, corner-only quadtree refinement on top of that, so the corner
    cell size shrinks FASTER than the background by the ratio dictated by
    graded-mesh theory: h_corner ~ h_bg^{(k+1)/gamma}, i.e. (since both
    refinements are dyadic) total corner depth L(ell) = ell*(k+1)/gamma =
    1.5*(k+1)*ell, so EXTRA corner-only levels beyond the uniform part is
    (1.5*(k+1) - 1)*ell.

Every cell is still an axis-aligned quad (or a pentagon at a single
hanging-node interface) -- see cells_to_polygons() -- so there is still no
angle or aspect-ratio degeneracy anywhere; what changed is only which
cells get refined and how many levels, not the underlying quadtree/
hanging-node machinery.
===============================================================================
"""

import os
import math
import argparse
from fractions import Fraction as F

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meshes")

CORNER = F(1, 2)
GAMMA = F(2, 3)  # regularity exponent of the 270 deg reentrant corner


def _in_notch(x0, y0):
    return x0 >= CORNER and y0 >= CORNER


def _uniform_refine(cells):
    """One uniform quadtree pass on every cell (background refinement)."""
    new_cells = []
    for (x0, x1, y0, y1, lv) in cells:
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        for a0, a1, b0, b1 in [(x0, xm, y0, ym), (xm, x1, y0, ym),
                                (x0, xm, ym, y1), (xm, x1, ym, y1)]:
            if _in_notch(a0, b0):
                continue
            new_cells.append((a0, a1, b0, b1, lv + 1))
    return new_cells


def _corner_refine(cells, level, halo=1.0):
    """One corner-only quadtree pass: leaf cells at `level-1` within a
    neighbourhood of the corner get split.

    The neighbourhood is defined RELATIVE to each cell's own size (Chebyshev
    distance from the cell centre to the corner <= (0.5 + halo) * cell_size),
    not just "the single cell whose bounding box contains the corner point"
    (halo=0). This is the standard graded/geometric-mesh construction from
    the corner-singularity FEM literature (Apel, Babuska-Guo): a
    self-similar halo of cells around the singularity gets refined at every
    level, not an isolated chain of touching cells one cell wide. A wider
    halo gives a gentler size transition between the corner region and the
    background (fewer, less extreme hanging-node jumps per level), at the
    cost of more cells per level.
    """
    new_cells = []
    for (x0, x1, y0, y1, lv) in cells:
        s = x1 - x0
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        cheb = max(abs(cx - CORNER), abs(cy - CORNER))
        touches_corner = (lv == level - 1) and (cheb <= (F(1, 2) + halo) * s)
        if not touches_corner:
            new_cells.append((x0, x1, y0, y1, lv))
            continue
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        for a0, a1, b0, b1 in [(x0, xm, y0, ym), (xm, x1, y0, ym),
                                (x0, xm, ym, y1), (xm, x1, ym, y1)]:
            if _in_notch(a0, b0):
                continue
            new_cells.append((a0, a1, b0, b1, level))
    return new_cells


def corner_boost(k):
    """Extra corner-refinement ratio dictated by classical graded-mesh
    theory: total corner depth should scale like (k+1)/gamma times the
    uniform/background depth."""
    return F(k + 1, 1) / GAMMA  # = 1.5*(k+1)


def build_quadtree_cells(n0, ell, k, halo=1.0):
    """n0: base grid resolution (even). ell: convergence level (background
    refined `ell` times uniformly). k: HHO polynomial degree, used only to
    compute how many EXTRA corner-only levels are needed on top of the
    uniform part (graded-mesh theory, see module docstring). halo: width
    (in units of each cell's own size) of the neighbourhood around the
    corner refined at every corner-only pass, see _corner_refine()."""
    if n0 % 2 != 0:
        raise ValueError("n0 must be even so that 0.5 is a base grid line")
    h0 = F(1, n0)
    cells = []
    for i in range(n0):
        for j in range(n0):
            x0, x1 = i * h0, (i + 1) * h0
            y0, y1 = j * h0, (j + 1) * h0
            if _in_notch(x0, y0):
                continue
            cells.append((x0, x1, y0, y1, 0))

    for _ in range(ell):
        cells = _uniform_refine(cells)

    total_depth = int(round(float(corner_boost(k) * ell)))
    for level in range(ell + 1, total_depth + 1):
        cells = _corner_refine(cells, level, halo=halo)

    return [(x0, x1, y0, y1) for (x0, x1, y0, y1, lv) in cells]


def cells_to_polygons(cells):
    """Insert hanging-node vertices along each cell's 4 edges, using ONLY
    other cells that genuinely have an edge on the same coordinate line
    with an overlapping interval (not a global coordinate-value match --
    that earlier approach spuriously injected hanging nodes from unrelated,
    far-away refined regions into any cell that happened to share a
    coordinate value, silently splitting true single faces into several,
    a real HHO face-matching bug)."""
    from collections import defaultdict

    horiz = defaultdict(list)
    vert = defaultdict(list)
    for (x0, x1, y0, y1) in cells:
        horiz[y0].append((x0, x1))
        horiz[y1].append((x0, x1))
        vert[x0].append((y0, y1))
        vert[x1].append((y0, y1))

    def hanging_on_horizontal(y, lo, hi):
        pts = set()
        for a, b in horiz.get(y, []):
            if lo < a < hi:
                pts.add(a)
            if lo < b < hi:
                pts.add(b)
        return sorted(pts)

    def hanging_on_vertical(x, lo, hi):
        pts = set()
        for a, b in vert.get(x, []):
            if lo < a < hi:
                pts.add(a)
            if lo < b < hi:
                pts.add(b)
        return sorted(pts)

    point_index = {}
    points = []

    def get_pid(x, y):
        key = (x, y)
        if key not in point_index:
            point_index[key] = len(points)
            points.append((float(x), float(y)))
        return point_index[key]

    polygons = []
    for (x0, x1, y0, y1) in cells:
        bottom = [(x, y0) for x in [x0] + hanging_on_horizontal(y0, x0, x1) + [x1]]
        right  = [(x1, y) for y in [y0] + hanging_on_vertical(x1, y0, y1) + [y1]]
        top    = [(x, y1) for x in [x1] + hanging_on_horizontal(y1, x0, x1)[::-1] + [x0]]
        left   = [(x0, y) for y in [y1] + hanging_on_vertical(x0, y0, y1)[::-1] + [y0]]
        ring = bottom[:-1] + right[:-1] + top[:-1] + left[:-1]
        polygons.append([get_pid(x, y) for x, y in ring])

    bc_ids = set()
    for pid, (x, y) in enumerate(points):
        on_outer = x == 0.0 or x == 1.0 or y == 0.0 or y == 1.0
        on_cut_x = abs(x - 0.5) < 1e-12 and y >= 0.5 - 1e-12
        on_cut_y = abs(y - 0.5) < 1e-12 and x >= 0.5 - 1e-12
        if on_outer or on_cut_x or on_cut_y:
            bc_ids.add(pid)

    return points, polygons, sorted(bc_ids)


def generate_mesh(ell, k, n0=4, halo=1.0):
    """`ell`: convergence level (uniform background refinement count).
    `k`: HHO polynomial degree (determines the extra corner-only boost,
    see module docstring). `halo`: corner-neighbourhood width, see
    _corner_refine()."""
    cells = build_quadtree_cells(n0=n0, ell=ell, k=k, halo=halo)
    return cells_to_polygons(cells)


def write_poly2d(path, points, polys, bc_ids):
    with open(path, "w") as f:
        f.write(f"{len(points)} {len(polys)} 1\n")
        for x, y in points:
            f.write(f"{x:.15g} {y:.15g}\n")
        for poly in polys:
            f.write(f"{len(poly)} " + " ".join(str(v + 1) for v in poly) + "\n")
        f.write(" ".join(str(i + 1) for i in bc_ids) + "\n")


def h_min_max(points, polys):
    def dist(a, b):
        return math.hypot(points[a][0] - points[b][0], points[a][1] - points[b][1])

    def diam(poly):
        return max(dist(a, b) for idx, a in enumerate(poly) for b in poly[idx + 1:])

    diams = [diam(p) for p in polys]
    return min(diams), max(diams)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, required=True,
                     help="HHO polynomial degree the mesh is graded for (sets the corner boost)")
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2],
                     help="convergence level indices ell (uniform background refinement count)")
    ap.add_argument("--reference", type=int, default=3,
                     help="ell for the fine reference mesh")
    ap.add_argument("--halo", type=float, default=1.0,
                     help="width (in units of each cell's own size) of the "
                          "neighbourhood around the corner refined at every "
                          "corner-only pass; 0 = only the single cell whose "
                          "bounding box contains the corner point (old, too "
                          "narrow, behaviour), larger = wider/gentler halo")
    ap.add_argument("--plot", action="store_true", help="save a PNG plot of each mesh")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    all_levels = list(args.levels) + [args.reference]
    manifest = []
    for ell in all_levels:
        points, polys, bc_ids = generate_mesh(ell, args.k, halo=args.halo)
        fname = f"lshape_graded_k{args.k}_N{ell}.txt"
        path = os.path.join(OUTDIR, fname)
        write_poly2d(path, points, polys, bc_ids)
        hmin, hmax = h_min_max(points, polys)
        tag = "reference" if ell == args.reference else "level"
        total_depth = int(round(float(corner_boost(args.k) * ell)))
        manifest.append((ell, path, len(points), len(polys), hmin, hmax))
        print(f"[{tag:9s}] ell={ell:3d}  (bg depth={ell}, corner depth={total_depth})  "
              f"points={len(points):6d}  cells={len(polys):6d}  "
              f"h_min={hmin:.5g}  h_max={hmax:.5g}  p~{round(hmax/hmin) if hmin>0 else float('inf')}  -> {path}")

        if args.plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 6))
            for poly in polys:
                xs = [points[v][0] for v in poly] + [points[poly[0]][0]]
                ys = [points[v][1] for v in poly] + [points[poly[0]][1]]
                ax.plot(xs, ys, "-", color="#2a78d6", lw=0.5)
            ax.set_aspect("equal")
            ax.set_title(f"L-shape quadtree mesh, ell={ell}, k={args.k}")
            png_path = os.path.join(OUTDIR, f"lshape_graded_k{args.k}_N{ell}.png")
            fig.savefig(png_path, dpi=150)
            plt.close(fig)
            print(f"             plot -> {png_path}")

    manifest_path = os.path.join(OUTDIR, "manifest.txt")
    with open(manifest_path, "w") as f:
        f.write("# ell  file  n_points  n_cells  h_min  h_max\n")
        for ell, path, npt, nq, hmin, hmax in manifest:
            f.write(f"{ell} {os.path.basename(path)} {npt} {nq} {hmin:.15g} {hmax:.15g}\n")
    print(f"\nManifest -> {manifest_path}")


if __name__ == "__main__":
    main()
