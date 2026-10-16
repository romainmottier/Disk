#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Conforming graded quadrilateral meshes of the L-shaped domain (no hanging
node), corner at the origin:
    Omega = (-1/2, 1/2)^2 \\ [0, 1/2) x (0, 1/2],
for the singular solutions r^gamma sin(gamma theta), gamma = 2/3
(problem.lshape_corner = {0, 0}).
===============================================================================
A uniform grid of n x n squares of (-1/2, 1/2)^2 (n even, the quadrant
[0, 1/2) x (0, 1/2] removed: 3 n^2 / 4 cells) is mapped radially in the max
norm |x|_inf = max(|x_1|, |x_2|):
        x = xh (2 |xh|_inf)^(1/mu - 1),   i.e.   |x|_inf = (1/2) (2 |xh|_inf)^(1/mu),
which maps the square (and the two edges of the corner, the axes) onto
itself, and grades the cells towards the corner: a cell at max-norm distance
rho has the radial size ~ (h / mu) (2 rho)^(1 - mu) and the tangential size
~ h (2 rho)^(1 - mu), h = 1/n (aspect ratio <= 1/mu), the condition of
Babuska-Kellogg-Pitkaranta for the approximation of r^alpha at the optimal
rate h^(k+1) with mu < (alpha + 1)/(k + 1): mu = gamma/(k + 1) for the
gradient (k = 1: 1/3, k = 2: 2/9, k = 3: 1/6). The cells are convex
quadrilaterals (trapezoids near the diagonals), the mesh is conforming.

Output: the "polygons" format of wave_propagation (mesh.type = "polygons",
mesh/domain.hpp): "n_points n_cells", points "x y", cells "4 i_1 .. i_4 tag"
(0-based, counterclockwise, tag 1), and a summary line (cells, side min/max).
Usage: generate_lshape_radial_graded_mesh.py --n 16 --k 1 --output mesh.txt [--mu MU]
===============================================================================
"""
import argparse


def build(n, mu):
    """Vertices of the mapped grid (only those of the L) and the quadrilaterals."""
    def mapped(i, j):
        xh, yh = -0.5 + i / n, -0.5 + j / n
        s = max(abs(xh), abs(yh))
        if s == 0.0:
            return 0.0, 0.0
        g = (2.0 * s) ** (1.0 / mu - 1.0)
        return xh * g, yh * g

    def in_l(i, j):          # cell [i, i+1] x [j, j+1] of the grid, removed quadrant: x >= 0 and y >= 0
        return not (i >= n // 2 and j >= n // 2)

    ids, pts, cells = {}, [], []

    def vid(i, j):
        if (i, j) not in ids:
            ids[(i, j)] = len(pts)
            pts.append(mapped(i, j))
        return ids[(i, j)]
    for j in range(n):
        for i in range(n):
            if in_l(i, j):
                cells.append([vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)])
    return pts, cells


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, required=True, help="grid cells per unit length (even)")
    ap.add_argument("--k", type=int, default=1, help="polynomial degree (sets mu = gamma / (k + 1))")
    ap.add_argument("--mu", type=float, default=None, help="grading exponent (default (2/3) / (k + 1))")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    if a.n % 2:
        raise SystemExit("--n must be even")
    mu = a.mu if a.mu is not None else (2.0 / 3.0) / (a.k + 1)
    pts, cells = build(a.n, mu)
    with open(a.output, "w") as f:
        f.write("# L-shaped domain, corner at the origin, conforming radially graded quadrilaterals "
                "(generate_lshape_radial_graded_mesh.py)\n")
        f.write(f"{len(pts)} {len(cells)}\n")
        for (x, y) in pts:
            f.write(f"{x:.17g} {y:.17g}\n")
        for c in cells:
            f.write("4 " + " ".join(map(str, c)) + " 1\n")
    diam = []
    for c in cells:
        P = [pts[v] for v in c]
        diam.append(max(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5 for p in P for q in P))
    print(f"{a.output}: n = {a.n}, mu = {mu:.4f}: {len(cells)} cells, diameter {min(diam):.3e} .. {max(diam):.3e}, "
          f"ratio {max(diam) / min(diam):.0f}")


if __name__ == "__main__":
    main()
