#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Meshes of the study of what sets the stable step of an HHO cell: the unit
square, n x n squares, one modification of the central square:
  thin     the central square split at x0 + eps h into a thin rectangle
           (eps h x h) and a wide one; the squares above and below get the
           new vertex (a hanging node: a face of length eps h);
  face     a vertex inserted on the top edge of the central square at
           x0 + eps h (shared with the square above): two square cells with
           a face of length eps h, nothing else changes;
  corner   the corner (x0, y0) of the central square cut off: a triangle of
           legs eps h and a pentagon; the left and lower squares get the new
           vertices;
  aniso    every cell a rectangle h x h / a (eps = 1/a), no other change;
  slab     a whole column split at eps h: thin rectangles in every row.
Output: "polygons" format (tag 1). Prints the smallest diameter, face,
4 |T| / |dT| and |T| / diam over the cells.
Usage: generate_cfl_cells.py --kind thin --eps 0.125 --n 16 --output mesh.txt
===============================================================================
"""
import argparse
import numpy as np


def square_mesh(n, ny=None):
    ny = ny or n
    pts = {}
    def pid(x, y):
        k = (round(x * 2 ** 40), round(y * 2 ** 40))
        if k not in pts:
            pts[k] = (len(pts), (x, y))
        return pts[k][0]
    return pts, pid


def build(kind, eps, n):
    h = 1.0 / n
    c = n // 2
    x0, y0 = c * h, c * h
    pts = {}
    order = []

    def pid(x, y):
        k = (round(x * 2 ** 40), round(y * 2 ** 40))
        if k not in pts:
            pts[k] = len(order)
            order.append((x, y))
        return pts[k]

    cells = []
    if kind == "aniso":
        a = 1.0 / eps
        ny = int(round(n * a))
        for j in range(ny):
            for i in range(n):
                X0, Y0, hy = i * h, j * h / a, h / a
                cells.append([pid(X0, Y0), pid(X0 + h, Y0), pid(X0 + h, Y0 + hy), pid(X0, Y0 + hy)])
        return order, cells
    xs = x0 + eps * h
    for j in range(n):
        for i in range(n):
            X0, Y0 = i * h, j * h
            ring = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0, Y0 + h)]
            if kind == "thin" and i == c and j == c:
                cells.append([pid(*p) for p in [(X0, Y0), (xs, Y0), (xs, Y0 + h), (X0, Y0 + h)]])
                cells.append([pid(*p) for p in [(xs, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (xs, Y0 + h)]])
                continue
            if kind == "slab" and i == c:
                cells.append([pid(*p) for p in [(X0, Y0), (xs, Y0), (xs, Y0 + h), (X0, Y0 + h)]])
                cells.append([pid(*p) for p in [(xs, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (xs, Y0 + h)]])
                continue
            if kind == "thin" and i == c and j == c + 1:      # above: hanging node on its bottom edge
                ring = [(X0, Y0), (xs, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0, Y0 + h)]
            if kind == "thin" and i == c and j == c - 1:      # below: on its top edge
                ring = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (xs, Y0 + h), (X0, Y0 + h)]
            if kind == "face" and i == c and j == c:
                ring = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (xs, Y0 + h), (X0, Y0 + h)]
            if kind == "face" and i == c and j == c + 1:
                ring = [(X0, Y0), (xs, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0, Y0 + h)]
            if kind == "corner":
                e = eps * h
                if i == c and j == c:
                    cells.append([pid(*p) for p in [(X0, Y0), (X0 + e, Y0), (X0, Y0 + e)]])
                    ring = [(X0 + e, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0, Y0 + h), (X0, Y0 + e)]
                if i == c - 1 and j == c:                     # left: vertex on its right edge
                    ring = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + e), (X0 + h, Y0 + h), (X0, Y0 + h)]
                if i == c and j == c - 1:                     # below: vertex on its top edge
                    ring = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0 + e, Y0 + h), (X0, Y0 + h)]
            cells.append([pid(*p) for p in ring])
    return order, cells


def measures(order, cells):
    P = np.array(order)
    out = dict(diam=np.inf, face=np.inf, ap=np.inf, ad=np.inf)
    for r in cells:
        q = P[r]
        A = 0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1])
        e = np.hypot(*(np.roll(q, -1, axis=0) - q).T)
        d = max(np.hypot(*(q[i] - q[j])) for i in range(len(q)) for j in range(i + 1, len(q)))
        out["diam"] = min(out["diam"], d); out["face"] = min(out["face"], e.min())
        out["ap"] = min(out["ap"], 4 * A / e.sum()); out["ad"] = min(out["ad"], A / d)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=["none", "thin", "face", "corner", "aniso", "slab"], required=True)
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    order, cells = build(a.kind, a.eps, a.n)
    with open(a.output, "w") as f:
        f.write(f"# generate_cfl_cells.py --kind {a.kind} --eps {a.eps} --n {a.n}\n")
        f.write(f"{len(order)} {len(cells)}\n")
        for (x, y) in order:
            f.write(f"{x:.15g} {y:.15g}\n")
        for r in cells:
            f.write(f"{len(r)} " + " ".join(map(str, r)) + " 1\n")
    m = measures(order, cells)
    print(f"{a.output}: {len(cells)} cells; min diameter {m['diam']:.5g}, face {m['face']:.5g}, 4|T|/|dT| {m['ap']:.5g}, "
          f"|T|/diam {m['ad']:.5g}")


if __name__ == "__main__":
    main()
