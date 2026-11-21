#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Removes the cells of zero area of a mesh in the "polygons" format of wave_propagation (2026-11-16).

generate_head_model_fitted.py may leave a flat "cell" whose vertices lie on one line: it joins a cell side holding a
hanging node (two small cells) to the side of a larger cell without this node. Such a cell has a zero mass matrix (NaN
in every run) and a zero width (MLTS with time.level_size = "width" loops forever). Each one is removed, and its
interior vertices are inserted, in order, into the side of the neighbouring cell that runs between its two extreme
vertices (hanging nodes of that cell). The output is checked: no cell of zero area, and every side of a cell is shared
by exactly two cells, except on the boundary of the domain.
Usage: python3 fix_degenerate_cells.py <input mesh> <output mesh>
"""
import sys
import numpy as np


def read(path):
    header, body = [], []
    for line in open(path):
        (header if line.startswith("#") else body).append(line)
    body = [l for l in body if l.strip()]
    nv, nc = map(int, body[0].split()[:2])
    pts = np.array([list(map(float, body[1 + i].split()[:2])) for i in range(nv)])
    cells = []
    for j in range(nc):
        v = list(map(int, body[1 + nv + j].split()))
        cells.append((v[1:1 + v[0]], v[1 + v[0]:]))          # ring, trailing values (tag)
    return header, pts, cells


def area(pts, ring):
    P = pts[ring]
    return 0.5 * abs(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))


def sides(ring):
    return [(ring[i], ring[(i + 1) % len(ring)]) for i in range(len(ring))]


def main(src, dst):
    header, pts, cells = read(src)
    scale = np.ptp(pts, axis=0).max()
    flat = [j for j, (ring, _) in enumerate(cells) if area(pts, ring) < 1e-14 * scale * scale]
    owner = {}
    for j, (ring, _) in enumerate(cells):
        if j in flat:
            continue
        for a, b in sides(ring):
            owner.setdefault((min(a, b), max(a, b)), []).append(j)
    for j in flat:
        ring = cells[j][0]
        P = pts[ring]
        d = P[np.argmax(np.linalg.norm(P - P[0], axis=1))] - P[0]
        order = sorted(ring, key=lambda v: float((pts[v] - P[0]) @ d))      # along the line
        e0, e1, inner = order[0], order[-1], order[1:-1]
        key = (min(e0, e1), max(e0, e1))
        if key not in owner or len(owner[key]) != 1:
            sys.exit(f"cell {j}: no single cell owns the side {key}, cannot fix")
        k = owner[key][0]
        rk = cells[k][0]
        i = next(i for i, (a, b) in enumerate(sides(rk)) if {a, b} == {e0, e1})
        a, b = rk[i], rk[(i + 1) % len(rk)]
        ins = inner if a == e0 else inner[::-1]                               # from a to b
        cells[k] = (rk[:i + 1] + ins + rk[i + 1:], cells[k][1])
        print(f"cell {j} ({len(ring)} vertices on a line) removed; {len(inner)} hanging node(s) inserted in cell {k}")
    cells = [c for j, c in enumerate(cells) if j not in set(flat)]
    # checks
    bad = [j for j, (ring, _) in enumerate(cells) if area(pts, ring) < 1e-14 * scale * scale]
    count = {}
    for ring, _ in cells:
        for a, b in sides(ring):
            key = (min(a, b), max(a, b))
            count[key] = count.get(key, 0) + 1
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    def on_boundary(key):
        p, q = pts[key[0]], pts[key[1]]
        return any(abs(p[c] - x) < 1e-12 and abs(q[c] - x) < 1e-12 for c in (0, 1) for x in (lo[c], hi[c]))
    wrong = [key for key, n in count.items() if n > 2 or (n == 1 and not on_boundary(key))]
    if bad or wrong:
        sys.exit(f"check failed: {len(bad)} cells of zero area, {len(wrong)} sides not shared by two cells")
    with open(dst, "w") as f:
        f.writelines(header)
        f.write(f"{len(pts)} {len(cells)}\n")
        for x, y in pts:
            f.write(f"{float(x)!r} {float(y)!r}\n")                      # shortest exact decimal of each double
        for ring, rest in cells:
            f.write(" ".join(map(str, [len(ring)] + ring + rest)) + "\n")
    print(f"{dst}: {len(cells)} cells, checks passed (no zero area, every inner side shared by two cells)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
