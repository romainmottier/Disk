#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wrapper around generate_center_square_fixed_mesh.py's generate_mesh() to
produce the "centersq_pfix{pfam}_k{k}_N{N}.txt" mesh family: n0=4, halo=0.0,
extra=pfam (so p_global = 2^pfam is CONSTANT across the whole N=0..4 sweep).
This is the naming convention read directly by ERK4_MLTS_CenterSquare_
TwoLevel_conv_test.hpp, ERK4_MLTS_CenterSquare_ClassicalTwoLevel_conv_test.hpp
and ERK4_MLTS_CenterSquare_FixedLevel_conv_test.hpp (via MLTS_PFAM).

pfam=1..5 (p=2,4,8,16,32) already exist from an earlier ad-hoc run; this
script exists so deeper families (pfam=6,8,10 -> p=64,256,1024, needed for
genuine multilevel L=3,4,5 demonstrations) are reproducible.

Usage: python3 generate_pfix_family.py --pfam 6 8 10 --levels 0 1 2 3 4
"""

import os
import argparse
import generate_center_square_fixed_mesh as g

OUTDIR = g.OUTDIR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--pfam", type=int, nargs="+", required=True,
                     help="exponents: p_global = 2^pfam")
    ap.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    for pfam in args.pfam:
        for N in args.levels:
            pts, polys, bc_ids = g.generate_mesh(ell=N, k=args.k, n0=4, halo=0.0, extra=pfam)
            fname = f"centersq_pfix{pfam}_k{args.k}_N{N}.txt"
            path = os.path.join(OUTDIR, fname)
            g.write_poly2d(path, pts, polys, bc_ids)
            hmin, hmax = g.h_min_max(pts, polys)
            p = round(hmax / hmin) if hmin > 0 else float("inf")
            print(f"[pfam={pfam:2d} N={N}]  points={len(pts):6d}  cells={len(polys):6d}  "
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
                ax.set_title(f"centersq_pfix{pfam}, N={N} (p~{p})", fontsize=11)
                png_path = os.path.join(OUTDIR, f"centersq_pfix{pfam}_k{args.k}_N{N}.png")
                fig.savefig(png_path, dpi=150)
                plt.close(fig)


if __name__ == "__main__":
    main()
