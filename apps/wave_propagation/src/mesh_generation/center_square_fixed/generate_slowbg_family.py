#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPERIMENTAL prototype: unlike generate_fixedbg_family.py (ell=0 FIXED
forever -- confirmed empirically to plateau immediately, N=1,2,3 all
giving the same ~2.456e-6 error, since h_max never shrinks), this
family lets the BACKGROUND also refine, just much more slowly than the
local corner depth: ell = N//2 (background halves every 2 points, not
every 1), n0=4 (minimal -- only 1 ring of cells directly touches the
domain boundary at any given resolution, no spare/redundant coarse
shell), extra = 2*(N+1) (local depth keeps growing every point, same
schedule as before).

Every cell in the mesh, including the boundary-adjacent ring, shrinks
as N grows (ell increases) -- there is no permanently-frozen region,
addressing the plateau. Cost still grows much more slowly than the OLD
centersquareadaptive_graded family (ell=N, background refines at the
SAME rate as local depth): 40, 64, 136, 160, 376, 400, 1192 cells for
N=0..6, vs exponential blowup before.

Usage: python3 generate_slowbg_family.py --levels 0 1 2 3 4 5 6
"""

import os
import argparse
import generate_center_square_fixed_mesh as g

OUTDIR = g.OUTDIR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n0", type=int, default=4)
    ap.add_argument("--halo", type=float, default=0.01)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    for N in args.levels:
        ell = N // 2
        extra = 2 * (N + 1)
        pts, polys, bc_ids = g.generate_mesh(ell=ell, k=args.k, n0=args.n0, halo=args.halo, extra=extra)
        fname = f"centersquareslowbg_graded_k{args.k}_N{N}.txt"
        path = os.path.join(OUTDIR, fname)
        g.write_poly2d(path, pts, polys, bc_ids)
        hmin, hmax = g.h_min_max(pts, polys)
        p = round(hmax / hmin) if hmin > 0 else float("inf")
        print(f"[N={N}]  ell={ell}  extra={extra}  points={len(pts):6d}  cells={len(polys):6d}  "
              f"h_min={hmin:.6g}  h_max={hmax:.6g}  p~{p}  -> {path}")


if __name__ == "__main__":
    main()
