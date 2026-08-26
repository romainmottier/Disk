#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Graded quadtree-refined mesh for the PLAIN UNIT SQUARE (0,1)^2, in DiSk++'s
"poly_2d" text format (read by polygon_2d_mesh_reader).

Same construction and the SAME quadtree/hanging-node machinery as
../lshape_graded/generate_lshape_graded_mesh.py, minus the notch cut-out:
grading is centred on the ORDINARY (non-reentrant) corner (0, 0) of the
plain square, matching the L-SHAPE's OWN topology (graded region touches
only 2 of the domain's boundary sides, not all 4 the way a CENTER-graded
mesh does) while still carrying no solution singularity there. Built to
test whether ../square_graded/'s CENTER-graded variant's anomalous
N=1->N=2 error INCREASE (unlike the L-shape, which converges cleanly at
every tested depth with the identical algorithm) is a topology artifact
of grading a region surrounded on all 4 sides, rather than a corner
region open on 2 sides like every other mesh validated this session.
===============================================================================
"""

import os
import math
import argparse
from fractions import Fraction as F

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meshes")

CENTER = F(1, 2)
GAMMA = F(2, 3)  # kept identical to the L-shape script purely so corner_boost(k)
                  # produces the SAME grading depth per level -- there is no
                  # regularity reason for this exponent here (the solution is
                  # smooth), it only controls how aggressively the mesh grades
                  # toward the corner, exercising the multi-level machinery to
                  # a comparable depth as the L-shape meshes.


def _uniform_refine(cells):
    """One uniform quadtree pass on every cell (background refinement)."""
    new_cells = []
    for (x0, x1, y0, y1, lv) in cells:
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        for a0, a1, b0, b1 in [(x0, xm, y0, ym), (xm, x1, y0, ym),
                                (x0, xm, ym, y1), (xm, x1, ym, y1)]:
            new_cells.append((a0, a1, b0, b1, lv + 1))
    return new_cells


def _center_refine(cells, level, halo=1.0):
    """One corner-only quadtree pass: leaf cells at `level-1` within a
    neighbourhood of the corner get split. Identical construction to the
    L-shape script's _corner_refine(), just without the notch exclusion."""
    new_cells = []
    for (x0, x1, y0, y1, lv) in cells:
        s = x1 - x0
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        cheb = max(abs(cx - CENTER), abs(cy - CENTER))
        touches_center = (lv == level - 1) and (cheb <= (F(1, 2) + halo) * s)
        if not touches_center:
            new_cells.append((x0, x1, y0, y1, lv))
            continue
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        for a0, a1, b0, b1 in [(x0, xm, y0, ym), (xm, x1, y0, ym),
                                (x0, xm, ym, y1), (xm, x1, ym, y1)]:
            new_cells.append((a0, a1, b0, b1, level))
    return new_cells


def corner_boost(k):
    """Same formula as the L-shape script, kept only to reproduce a
    comparable grading DEPTH (number of bands) per level -- see module
    docstring."""
    return F(k + 1, 1) / GAMMA  # = 1.5*(k+1)


def build_quadtree_cells(n0, ell, k, halo=1.0, extra=6):
    """n0: base grid resolution (even). ell: convergence level (background
    refined `ell` times uniformly, this is what h_max is measured against
    for the convergence-order plot). k: HHO polynomial degree (unused now,
    kept for signature compatibility). halo: width (in units of each
    cell's own size) of the neighbourhood around the corner refined at
    every corner-only pass. `extra`: FIXED number of extra corner-only
    levels beyond the uniform background depth, the SAME at every ell --
    unlike the L-shape/original variable-depth scripts (extra = ell *
    (corner_boost(k)-1), growing linearly with ell -> p_global grows like
    32^ell, N=4 already needs ~1e6 fine steps), a fixed extra keeps
    p_global = 2^extra CONSTANT across the whole N sweep, so the
    multi-level machinery is exercised at the SAME depth (same L, same
    band structure) at every N while total cost grows only with the
    background mesh size (polynomial in ell), not exponentially --
    letting a full 5+ point convergence curve be computed in practical
    time even though the errors involved are already very small."""
    if n0 % 2 != 0:
        raise ValueError("n0 must be even so that 0.5 is a base grid line")
    h0 = F(1, n0)
    cells = []
    for i in range(n0):
        for j in range(n0):
            x0, x1 = i * h0, (i + 1) * h0
            y0, y1 = j * h0, (j + 1) * h0
            cells.append((x0, x1, y0, y1, 0))

    for _ in range(ell):
        cells = _uniform_refine(cells)

    total_depth = ell + extra
    for level in range(ell + 1, total_depth + 1):
        cells = _center_refine(cells, level, halo=halo)

    return [(x0, x1, y0, y1) for (x0, x1, y0, y1, lv) in cells]


def cells_to_polygons(cells):
    """Insert hanging-node vertices along each cell's 4 edges -- identical
    logic to the L-shape script's cells_to_polygons()."""
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
        if on_outer:
            bc_ids.add(pid)

    return points, polygons, sorted(bc_ids)


def generate_mesh(ell, k, n0=4, halo=1.0, extra=6):
    cells = build_quadtree_cells(n0=n0, ell=ell, k=k, halo=halo, extra=extra)
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
                     help="HHO polynomial degree the mesh is graded for (sets the centre boost)")
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3],
                     help="convergence level indices ell (uniform background refinement count)")
    ap.add_argument("--halo", type=float, default=1.0,
                     help="width (in units of each cell's own size) of the "
                          "neighbourhood around the centre refined at every "
                          "centre-only pass")
    ap.add_argument("--extra", type=int, default=6,
                     help="FIXED number of extra corner-only refinement "
                          "levels beyond the background depth, the SAME at "
                          "every N -- keeps p_global=2^extra constant across "
                          "the whole sweep instead of growing exponentially")
    ap.add_argument("--plot", action="store_true", help="save a PNG plot of each mesh")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    manifest = []
    for ell in args.levels:
        points, polys, bc_ids = generate_mesh(ell, args.k, halo=args.halo, extra=args.extra)
        fname = f"centersquarefixed_graded_k{args.k}_N{ell}.txt"
        path = os.path.join(OUTDIR, fname)
        write_poly2d(path, points, polys, bc_ids)
        hmin, hmax = h_min_max(points, polys)
        total_depth = ell + args.extra
        manifest.append((ell, path, len(points), len(polys), hmin, hmax))
        print(f"[level    ] ell={ell:3d}  (bg depth={ell}, centre depth={total_depth})  "
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
            ax.set_title(f"Square quadtree mesh, ell={ell}, k={args.k}")
            png_path = os.path.join(OUTDIR, f"centersquarefixed_graded_k{args.k}_N{ell}.png")
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
