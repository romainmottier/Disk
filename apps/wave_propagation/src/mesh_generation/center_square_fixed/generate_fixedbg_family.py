#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPERIMENTAL, first-draft prototype (not wired into any driver/campaign
yet): a mesh family where the BACKGROUND is fixed (ell=0 always -- no
uniform whole-domain refinement, so band0 stays at a small, constant
cell count for every N), and only the LOCAL corner/center-only
refinement depth (`extra`) grows with N to reach deeper L each step
(extra=2*(N+1) -> L=2,3,4,5,... via the same natural threshold rule
used everywhere else in this investigation, since each band boundary
needs a 4x size gap = 2 quadtree halvings).

KNOWN, ACCEPTED TRADE-OFF (not yet resolved with the user, purely for
visual review at this stage): h_max (background cell size) never
shrinks since ell=0 for every N, so this family does NOT give true
global h-convergence -- it is meant only as a cheap way to exercise
deep multilevel L values without band0 exploding, not to replace the
validated final_tests/square/ convergence campaigns.

Usage: python3 generate_fixedbg_family.py --levels 0 1 2 3
"""

import os
import argparse
import generate_center_square_fixed_mesh as g

OUTDIR = g.OUTDIR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n0", type=int, default=6)
    ap.add_argument("--halo", type=float, default=0.01)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3])
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    for N in args.levels:
        extra = 2 * (N + 1)
        pts, polys, bc_ids = g.generate_mesh(ell=0, k=args.k, n0=args.n0, halo=args.halo, extra=extra)
        fname = f"centersquarefixedbg_graded_k{args.k}_N{N}.txt"
        path = os.path.join(OUTDIR, fname)
        g.write_poly2d(path, pts, polys, bc_ids)
        hmin, hmax = g.h_min_max(pts, polys)
        p = round(hmax / hmin) if hmin > 0 else float("inf")
        print(f"[N={N}]  ell=0 (fixed)  extra={extra}  points={len(pts):6d}  cells={len(polys):6d}  "
              f"h_min={hmin:.5g}  h_max={hmax:.5g}  p~{p}  -> {path}")


if __name__ == "__main__":
    main()
