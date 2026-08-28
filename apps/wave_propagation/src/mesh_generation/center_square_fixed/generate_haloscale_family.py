#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FINAL, accepted construction: `ell = N` (background genuinely halves
every N -- true global convergence, standard refinement schedule), and
a SMALL, FIXED local halo (0.01, same as generate_fixedbg_family.py)
for all `extra` local passes.

After 5 different attempts to also keep band0 thin (single big halo,
per-step re-tuned halo, outward "shift" peeling, ...), each failing a
different way (size inversion, empty bands, disproportionate band1,
cost explosion), this is the accepted trade-off: band0 = the ENTIRE
ell-refined background (large, but a SINGLE uniform size -- simple,
cheap to build, and correctly classified since the power-of-2
band-assignment fix). All local bands (band1..band(L-1)) come out
clean and evenly-sized (~24 cells each, 2-cell-thick radial shells),
exactly like generate_fixedbg_family.py's already-validated behavior,
because the small halo never has to "jump" -- it always operates on an
ordinary small local patch.

Usage: python3 generate_haloscale_family.py --levels 0 1 2 3 4 5 6
"""

import os
import argparse
import generate_center_square_fixed_mesh as g

OUTDIR = g.OUTDIR
SMALL_HALO = 0.01


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n0", type=int, default=4)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    for N in args.levels:
        ell = N
        extra = 2 * (N + 1)
        pts, polys, bc_ids = g.generate_mesh(ell=ell, k=args.k, n0=args.n0, halo=SMALL_HALO, extra=extra)
        fname = f"centersquarehaloscale_graded_k{args.k}_N{N}.txt"
        path = os.path.join(OUTDIR, fname)
        g.write_poly2d(path, pts, polys, bc_ids)
        hmin, hmax = g.h_min_max(pts, polys)
        p = round(hmax / hmin) if hmin > 0 else float("inf")
        print(f"[N={N}]  ell={ell}  extra={extra}  cells={len(polys):7d}  "
              f"h_min={hmin:.6g}  h_max={hmax:.6g}  p~{p}  -> {path}")


if __name__ == "__main__":
    main()
