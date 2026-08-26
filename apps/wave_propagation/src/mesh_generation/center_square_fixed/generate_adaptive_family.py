#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"centersquareadaptive_graded_k{k}_N{N}.txt" mesh family: n0=4, halo=0.0,
extra=2*N, so p_global = 4^N GROWS with N (1, 4, 16, 64, 256 for N=0..4).
Fed through the natural (unforced) threshold rule already used by
ERK4_MLTS_CenterSquare_FixedLevel_conv_test.hpp / ERK4_MLTS_CenterSquareCoarse_
PellSparse_conv_test.hpp, this p growth maps to L = 1, 1, 2, 3, 4 across
N=0..4 -- a genuine, monotonically-growing multilevel trace (unlike the
earlier centersquarecoarse_graded family, extra=N, which only reached L=2
at its very last point).

Usage: python3 generate_adaptive_family.py --levels 0 1 2 3 4
"""

import os
import argparse
import generate_center_square_fixed_mesh as g

OUTDIR = g.OUTDIR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    for N in args.levels:
        extra = 2 * N
        pts, polys, bc_ids = g.generate_mesh(ell=N, k=args.k, n0=4, halo=0.0, extra=extra)
        fname = f"centersquareadaptive_graded_k{args.k}_N{N}.txt"
        path = os.path.join(OUTDIR, fname)
        g.write_poly2d(path, pts, polys, bc_ids)
        hmin, hmax = g.h_min_max(pts, polys)
        p = round(hmax / hmin) if hmin > 0 else float("inf")
        print(f"[N={N}]  extra={extra}  points={len(pts):6d}  cells={len(polys):6d}  "
              f"h_min={hmin:.5g}  h_max={hmax:.5g}  p~{p}  -> {path}")

        if args.plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 6))
            for poly in polys:
                xs = [pts[v][0] for v in poly] + [pts[poly[0]][0]]
                ys = [pts[v][1] for v in poly] + [pts[poly[0]][1]]
                ax.plot(xs, ys, "-", color="#2a78d6", lw=0.5)
            ax.set_aspect("equal")
            ax.set_title(f"centersquareadaptive, N={N} (p~{p})", fontsize=11)
            png_path = os.path.join(OUTDIR, f"centersquareadaptive_graded_k{args.k}_N{N}.png")
            fig.savefig(png_path, dpi=150)
            plt.close(fig)


if __name__ == "__main__":
    main()
