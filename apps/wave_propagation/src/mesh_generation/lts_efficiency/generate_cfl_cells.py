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
  slab     a whole column split at eps h: thin rectangles in every row;
  cut      the central square cut by the line at angle --theta (of its
           normal, in degrees) and signed distance --offset h from its
           centre: two convex pieces, the cut points becoming hanging nodes
           of the neighbours;
  circle   every square crossed by the circle (--cx, --cy, --radius) cut
           along the chord between its two crossing points (an interface).
Output: "polygons" format (tag 1). Prints the smallest diameter, face,
4 |T| / |dT| and |T| / diam over the cells.
Usage: generate_cfl_cells.py --kind thin --eps 0.125 --n 16 --output mesh.txt
       generate_cfl_cells.py --kind cut --theta 30 --offset 0.2 --output mesh.txt
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


def clip(ring, nrm, d):
    """Part of the convex polygon ring where nrm . x <= d, and the crossing points (Sutherland-Hodgman)."""
    out, cross = [], []
    for i in range(len(ring)):
        a, b = np.array(ring[i]), np.array(ring[(i + 1) % len(ring)])
        fa, fb = nrm @ a - d, nrm @ b - d
        if fa <= 0:
            out.append(tuple(a))
        if fa * fb < 0:
            p = tuple(a + fa / (fa - fb) * (b - a))
            out.append(p); cross.append(p)
    return out, cross


def build_cut(kind, n, theta=0.0, offset=0.0, cx=0.5, cy=0.5, radius=0.25):
    """Squares, the central one cut by a line (cut) or each square crossed by a circle cut along its chord (circle);
    then every vertex lying inside an edge of another cell is inserted into it (hanging nodes)."""
    h = 1.0 / n
    rings = []
    for j in range(n):
        for i in range(n):
            X0, Y0 = i * h, j * h
            sq = [(X0, Y0), (X0 + h, Y0), (X0 + h, Y0 + h), (X0, Y0 + h)]
            if kind == "cut" and i == n // 2 and j == n // 2:
                t = np.radians(theta)
                nrm = np.array([np.cos(t), np.sin(t)])
                d = nrm @ np.array([X0 + h / 2, Y0 + h / 2]) + offset * h
            elif kind == "circle":
                inside = [np.hypot(x - cx, y - cy) < radius for x, y in sq]
                if all(inside) or not any(inside):
                    rings.append(sq); continue
                # crossing points of the circle with the edges, then the chord through them
                pts = []
                for k in range(4):
                    a, b = np.array(sq[k]), np.array(sq[(k + 1) % 4])
                    if inside[k] != inside[(k + 1) % 4]:
                        u = b - a; f = a - (cx, cy)
                        A, B, C = u @ u, 2 * f @ u, f @ f - radius ** 2
                        s = [(-B + sg * np.sqrt(B * B - 4 * A * C)) / (2 * A) for sg in (1, -1)]
                        pts.append(a + min((x for x in s if -1e-12 <= x <= 1 + 1e-12), key=lambda x: abs(x - 0.5)) * u)
                if len(pts) != 2:
                    raise SystemExit("circle: a square crossed more than twice")
                t = pts[1] - pts[0]
                nrm = np.array([-t[1], t[0]]) / np.hypot(*t)
                d = nrm @ pts[0]
            else:
                rings.append(sq); continue
            lo, _ = clip(sq, nrm, d)
            hi, _ = clip(sq, -nrm, -d)
            rings += [lo, hi]
    pts = {}
    order = []
    def pid(x, y):
        k = (round(x * 2 ** 40), round(y * 2 ** 40))
        if k not in pts:
            pts[k] = len(order)
            order.append((x, y))
        return pts[k]
    cells = [[pid(*p) for p in r] for r in rings]
    P = np.array(order)
    full = []
    for r in cells:
        ring = []
        for k in range(len(r)):
            a, b = P[r[k]], P[r[(k + 1) % len(r)]]
            u = b - a; L2 = u @ u
            s = (P - a) @ u / L2
            off = np.abs((P - a) @ np.array([-u[1], u[0]])) / np.sqrt(L2)
            on = [i for i in np.argsort(s) if 1e-12 < s[i] < 1 - 1e-12 and off[i] < 1e-12 * h]
            ring += [r[k]] + on
        full.append(ring)
    return order, full


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
    ap.add_argument("--kind", choices=["none", "thin", "face", "corner", "aniso", "slab", "cut", "circle"], required=True)
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--theta", type=float, default=0.0, help="cut: angle of the normal of the line (degrees)")
    ap.add_argument("--offset", type=float, default=0.0, help="cut: signed distance of the line to the centre / h")
    ap.add_argument("--cx", type=float, default=0.5, help="circle: centre")
    ap.add_argument("--cy", type=float, default=0.5)
    ap.add_argument("--radius", type=float, default=0.25)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    if a.kind in ("cut", "circle"):
        order, cells = build_cut(a.kind, a.n, a.theta, a.offset, a.cx, a.cy, a.radius)
        args = (f"--theta {a.theta} --offset {a.offset}" if a.kind == "cut" else
                f"--cx {a.cx} --cy {a.cy} --radius {a.radius}")
    else:
        order, cells = build(a.kind, a.eps, a.n)
        args = f"--eps {a.eps}"
    with open(a.output, "w") as f:
        f.write(f"# generate_cfl_cells.py --kind {a.kind} {args} --n {a.n}\n")
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
