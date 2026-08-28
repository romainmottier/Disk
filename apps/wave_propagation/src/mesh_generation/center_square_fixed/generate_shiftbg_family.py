#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPERIMENTAL prototype: "shift bands outward" construction, replacing
the failed ell/extra split (generate_slowbg_family.py: band0 becomes
huge -- many rings -- as soon as ell>0, since local refinement's reach
never scales with the background's resolution).

Instead of refining a small FIXED patch at the center while leaving the
background alone, this construction refines band0 (the CURRENT
boundary-touching ring) itself, at every step, and splits the result:
sub-cells that STILL touch the true domain boundary stay band0 (so
h_max genuinely shrinks every step); sub-cells that don't get PROMOTED
into a new band1, pushing every existing band index up by one (old
band1 -> band2, etc.). Two quadtree passes per step (not one) so the
promoted material forms a genuine 2-cell-thick radial shell -- a
1-cell-thick shell was found (empirically, via
MLTS_OVERLAP_RINGS-style adjacency dilation) to be fully swallowed by
the next-deeper band's protection margin, leaving it empty.

band0 is EXACTLY 1 cell thick at every step (h_max shrinks by 2x each
step): the fix requested after generate_slowbg_family.py showed a
"super wide yellow band" once ell>0. The innermost region (created at
N=0, never touched again by this operation) stays fixed size -- but
that's fine: for a globally smooth solution, the error is dominated by
h_max (the WORST/largest cell), not h_min, so h_max shrinking is what
actually matters for the convergence this family is meant to test.

Usage: python3 generate_shiftbg_family.py --levels 0 1 2 3 4 5 6
"""

import os
import argparse
from fractions import Fraction as F

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meshes")


def touches_boundary(x0, x1, y0, y1):
    return x0 == 0 or x1 == 1 or y0 == 0 or y1 == 1


def split_cell(x0, x1, y0, y1):
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    return [(x0, xm, y0, ym), (xm, x1, y0, ym), (x0, xm, ym, y1), (xm, x1, ym, y1)]


def shift_once(cells):
    """cells: list of (x0,x1,y0,y1,band). Refines band==0 cells by
    EXACTLY 2 quadtree passes, UNCONDITIONALLY (every band0 cell goes
    through both passes regardless of intermediate boundary status) --
    guarantees every resulting cell is the SAME size (no mixing of
    1-pass and 2-pass cells, which is what happened when pass 2 was
    only applied conditionally). Classification into new band0 (still
    touches the true domain boundary) vs promoted (doesn't) happens ONLY
    on the final, uniformly-2-passes-deep result. All other bands are
    shifted up by 1."""
    band0_cells = [(x0, x1, y0, y1) for (x0, x1, y0, y1, b) in cells if b == 0]
    other_cells = [(x0, x1, y0, y1, b + 1) for (x0, x1, y0, y1, b) in cells if b != 0]

    pass1 = [c for (x0, x1, y0, y1) in band0_cells for c in split_cell(x0, x1, y0, y1)]
    pass2 = [c for (x0, x1, y0, y1) in pass1 for c in split_cell(x0, x1, y0, y1)]

    new_band0, promoted = [], []
    for (x0, x1, y0, y1) in pass2:
        if touches_boundary(x0, x1, y0, y1):
            new_band0.append((x0, x1, y0, y1, 0))
        else:
            promoted.append((x0, x1, y0, y1, 1))

    return new_band0 + promoted + other_cells


def build_shiftbg_cells(n0, n_shifts):
    if n0 % 2 != 0:
        raise ValueError("n0 must be even")
    h0 = F(1, n0)
    cells = []
    for i in range(n0):
        for j in range(n0):
            x0, x1 = i * h0, (i + 1) * h0
            y0, y1 = j * h0, (j + 1) * h0
            band = 0 if touches_boundary(x0, x1, y0, y1) else 1
            cells.append((x0, x1, y0, y1, band))
    for _ in range(n_shifts):
        cells = shift_once(cells)
    return cells


def cells_to_polygons(cells):
    from collections import defaultdict
    horiz, vert = defaultdict(list), defaultdict(list)
    for (x0, x1, y0, y1, _b) in cells:
        horiz[y0].append((x0, x1)); horiz[y1].append((x0, x1))
        vert[x0].append((y0, y1)); vert[x1].append((y0, y1))

    def hang_h(y, lo, hi):
        pts = set()
        for a, b in horiz.get(y, []):
            if lo < a < hi: pts.add(a)
            if lo < b < hi: pts.add(b)
        return sorted(pts)

    def hang_v(x, lo, hi):
        pts = set()
        for a, b in vert.get(x, []):
            if lo < a < hi: pts.add(a)
            if lo < b < hi: pts.add(b)
        return sorted(pts)

    point_index, points = {}, []

    def gid(x, y):
        key = (x, y)
        if key not in point_index:
            point_index[key] = len(points)
            points.append((float(x), float(y)))
        return point_index[key]

    polygons, bands = [], []
    for (x0, x1, y0, y1, b) in cells:
        bottom = [(x, y0) for x in [x0] + hang_h(y0, x0, x1) + [x1]]
        right  = [(x1, y) for y in [y0] + hang_v(x1, y0, y1) + [y1]]
        top    = [(x, y1) for x in [x1] + hang_h(y1, x0, x1)[::-1] + [x0]]
        left   = [(x0, y) for y in [y1] + hang_v(x0, y0, y1)[::-1] + [y0]]
        ring = bottom[:-1] + right[:-1] + top[:-1] + left[:-1]
        polygons.append([gid(x, y) for x, y in ring])
        bands.append(b)

    bc_ids = sorted(pid for pid, (x, y) in enumerate(points) if x in (0.0, 1.0) or y in (0.0, 1.0))
    return points, polygons, bc_ids, bands


def write_poly2d(path, points, polys, bc_ids):
    with open(path, "w") as f:
        f.write(f"{len(points)} {len(polys)} 1\n")
        for x, y in points:
            f.write(f"{x:.15g} {y:.15g}\n")
        for poly in polys:
            f.write(f"{len(poly)} " + " ".join(str(v + 1) for v in poly) + "\n")
        f.write(" ".join(str(i + 1) for i in bc_ids) + "\n")


def h_min_max(points, polys):
    import math
    def diam(poly):
        pts = [points[v] for v in poly]
        return max(math.hypot(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1])
                   for i in range(len(pts)) for j in range(i+1, len(pts)))
    diams = [diam(p) for p in polys]
    return min(diams), max(diams)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n0", type=int, default=4)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    for N in args.levels:
        cells = build_shiftbg_cells(n0=args.n0, n_shifts=N)
        points, polys, bc_ids, bands = cells_to_polygons(cells)
        fname = f"centersquareshiftbg_graded_k{args.k}_N{N}.txt"
        path = os.path.join(OUTDIR, fname)
        write_poly2d(path, points, polys, bc_ids)
        hmin, hmax = h_min_max(points, polys)
        from collections import Counter
        c = Counter(bands)
        print(f"[N={N}]  cells={len(polys):6d}  h_min={hmin:.6g}  h_max={hmax:.6g}  "
              f"p~{round(hmax/hmin)}  L={max(bands)+1}  band_counts={dict(sorted(c.items()))}  -> {path}")


if __name__ == "__main__":
    main()
