#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
mu-graded quadtree meshes of the L-shaped domain
    Omega = (0,1)^2 \\ [1/2,1) x (1/2,1],   reentrant corner c = (1/2,1/2),
for the singular solutions r^gamma sin(gamma theta), gamma = 2/3.
===============================================================================
Grading (Babuska-Kellogg-Pitkaranta; Apel): a cell T of side s_T is split
while
        s_T > h * rho_T^(1 - mu),   rho_T = |c - x_T|  (x_T: centre of T),
from a uniform background of side h, then the quadtree is balanced 2:1 (at
most one hanging node per edge). A singular function r^alpha is approximated
at the optimal rate h^(k+1) in L^2 by polynomials of degree k when
mu < (alpha + 1) / (k + 1); for the first-order wave system p = u_t ~ r^gamma
and v = grad u ~ r^(gamma-1), the velocity gives the restrictive condition
        mu = gamma / (k + 1)   (k = 1: 1/3, k = 2: 2/9, k = 3: 1/6),
under which the smallest cell has side ~ h^(1/mu) and the number of cells
stays proportional to h^-2. --mu sets mu (default gamma / (k + 1)).

Output: the "polygons" format of wave_propagation (mesh.type = "polygons",
mesh/domain.hpp): "n_points n_cells", points "x y", cells
"nv i_1 .. i_nv tag" (0-based, counterclockwise, tag 1), hanging nodes as
polygon vertices; and a summary line (cells, h_min, h_max, levels).
Usage: generate_lshape_mu_graded_mesh.py --n 16 --k 3 [--corner-origin] --output mesh.txt
===============================================================================
"""
import argparse
from fractions import Fraction


def build(n, mu, max_depth=60):
    """Leaves (i, j, d): cell [i, i+1] x [j, j+1] * 2^-d / n (exact dyadic integers)."""
    h = 1.0 / n
    c = 0.5

    def in_domain(i, j, d):
        s = 1.0 / (n * 2 ** d)
        x0, y0 = i * s, j * s
        return not (x0 >= c - 1e-15 and y0 >= c - 1e-15)

    def needs_split(i, j, d):
        s = 1.0 / (n * 2 ** d)
        xc, yc = (i + 0.5) * s, (j + 0.5) * s
        rho = ((xc - c) ** 2 + (yc - c) ** 2) ** 0.5
        return s > h * rho ** (1.0 - mu) and d < max_depth

    leaves = {}
    stack = [(i, j, 0) for i in range(n) for j in range(n) if in_domain(i, j, 0)]
    while stack:
        i, j, d = stack.pop()
        if needs_split(i, j, d):
            stack += [(2 * i + a, 2 * j + b, d + 1) for a in (0, 1) for b in (0, 1)]
        else:
            leaves[(i, j, d)] = True
    return balance(leaves, n)


def balance(leaves, n):
    """2:1 balance: a leaf having a neighbour two levels deeper is split, until no change. A neighbour of depth >= d+2
    along a side of a leaf of depth d contains one of the points of that side at the multiples of a quarter side."""
    changed = True
    while changed:
        changed = False
        d_max = max(d for (_, _, d) in leaves)
        N = n * 2 ** d_max                       # integer coordinates at depth d_max
        C = N // 2                               # corner

        def inside(px, py):
            return 0 <= px < N and 0 <= py < N and not (px >= C and py >= C)

        def leaf_at_most(px, py, depth):
            for e in range(depth, -1, -1):
                if (px >> (d_max - e), py >> (d_max - e), e) in leaves:
                    return True
            return False

        for (i, j, d) in list(leaves):
            if d >= d_max - 1:
                continue
            sc = 2 ** (d_max - d)
            x0, y0 = i * sc, j * sc
            step = max(1, sc // 4)
            split = False
            for t in range(0, sc, step):
                for (px, py) in ((x0 - 1, y0 + t), (x0 + sc, y0 + t), (x0 + t, y0 - 1), (x0 + t, y0 + sc)):
                    if inside(px, py) and not leaf_at_most(px, py, d + 1):
                        split = True
                        break
                if split:
                    break
            if split:
                del leaves[(i, j, d)]
                for a in (0, 1):
                    for b in (0, 1):
                        leaves[(2 * i + a, 2 * j + b, d + 1)] = True
                changed = True
    return leaves


def write(path, leaves, n, shift=Fraction(0)):
    """Polygons with hanging nodes, exact dyadic coordinates (Fraction) printed with 17 digits."""
    d_max = max(d for (_, _, d) in leaves)
    D = 2 ** d_max

    def corners(i, j, d):
        s = 2 ** (d_max - d)
        return [(i * s, j * s), ((i + 1) * s, j * s), ((i + 1) * s, (j + 1) * s), (i * s, (j + 1) * s)], s

    vertex_set = set()
    for (i, j, d) in leaves:
        cs, _ = corners(i, j, d)
        vertex_set.update(cs)
    ids = {v: k for k, v in enumerate(sorted(vertex_set))}
    cells = []
    for (i, j, d) in sorted(leaves):
        cs, s = corners(i, j, d)
        ring = []
        for a in range(4):
            p, q = cs[a], cs[(a + 1) % 4]
            ring.append(p)
            if s % 2 == 0:
                m = ((p[0] + q[0]) // 2, (p[1] + q[1]) // 2)
                if m in ids:
                    ring.append(m)
        cells.append([ids[v] for v in ring])
    pts = sorted(ids, key=ids.get)
    with open(path, "w") as f:
        f.write("# L-shaped domain, mu-graded quadtree (generate_lshape_mu_graded_mesh.py)\n")
        f.write(f"{len(pts)} {len(cells)}\n")
        for (x, y) in pts:
            f.write(f"{float(Fraction(x, n * D) - shift):.17g} {float(Fraction(y, n * D) - shift):.17g}\n")
        for ring in cells:
            f.write(f"{len(ring)} " + " ".join(map(str, ring)) + " 1\n")
    sides = [1.0 / (n * 2 ** d) for (_, _, d) in leaves]
    return len(cells), min(sides), max(sides), d_max


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, required=True, help="background cells per unit length (even)")
    ap.add_argument("--k", type=int, default=3, help="polynomial degree (sets mu = gamma / (k + 1))")
    ap.add_argument("--mu", type=float, default=None, help="grading exponent (default (2/3) / (k + 1))")
    ap.add_argument("--output", required=True)
    ap.add_argument("--corner-origin", action="store_true",
                    help="coordinates minus (1/2, 1/2): corner at the origin, domain (-1/2,1/2)^2 \\ [0,1/2)x(0,1/2] "
                         "(full relative precision of the coordinates in the smallest cells; problem.lshape_corner = {0, 0})")
    a = ap.parse_args()
    mu = a.mu if a.mu is not None else (2.0 / 3.0) / (a.k + 1)
    leaves = build(a.n, mu)
    nc, smin, smax, dmax = write(a.output, leaves, a.n, Fraction(1, 2) if a.corner_origin else Fraction(0))
    print(f"{a.output}: n = {a.n}, mu = {mu:.4f}: {nc} cells, side {smin:.3e} .. {smax:.3e}, depth {dmax}, "
          f"ratio {smax / smin:.0f}")


if __name__ == "__main__":
    main()
