#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generic viewer for DiSk++'s poly_2d mesh text format (as written by
write_poly2d() in the various generate_*_mesh.py scripts under
mesh_generation/): plots every cell's outline. Used to produce mesh
illustration images for the multilevel LTS-RK4 convergence experiments
under prototypes/LTS/mlts_2026/square/, saved alongside their .txt
result files in build/apps/wave_propagation/lshape/results/.
"""

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_poly2d(path):
    with open(path) as f:
        header = f.readline().split()
        npts, npoly = int(header[0]), int(header[1])
        points = []
        for _ in range(npts):
            x, y = map(float, f.readline().split())
            points.append((x, y))
        polys = []
        for _ in range(npoly):
            parts = f.readline().split()
            n = int(parts[0])
            verts = [int(v) - 1 for v in parts[1:1 + n]]
            polys.append(verts)
    return points, polys


def plot_mesh(path, out_png, title):
    points, polys = read_poly2d(path)
    fig, ax = plt.subplots(figsize=(6, 6))
    for poly in polys:
        xs = [points[v][0] for v in poly] + [points[poly[0]][0]]
        ys = [points[v][1] for v in poly] + [points[poly[0]][1]]
        ax.plot(xs, ys, "-", color="#2a78d6", lw=0.6)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  {os.path.basename(path)} ({len(polys)} cells) -> {out_png}")


if __name__ == "__main__":
    mesh_txt = sys.argv[1]
    out_png = sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(mesh_txt)
    plot_mesh(mesh_txt, out_png, title)
