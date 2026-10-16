#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Locally refined quadtree meshes of the unit square for the efficiency study of
the local time-stepping schemes (LTS, MLTS) against ERK4.
===============================================================================
Background of n x n squares (side h = 1/n). Refined features (any number):
    box:x0,x1,y0,y1     rectangle (area refinement)
    disk:cx,cy,r        disk (area refinement)
    line:x0,y0,x1,y1    segment (refinement along a curve)
    point:cx,cy         point (refinement around a point)
Gradings:
    --grading two       exactly two cell sizes: the background squares meeting a
                        feature (with positive area for box and disk) are split
                        uniformly down to the depth L, the others are kept (abrupt
                        transition, up to 2^L - 1 hanging nodes per edge);
    --grading rings --rings W
                        a square of side s and depth d < L is split while its
                        distance to a feature is below W s: every level keeps a
                        band about W squares wide around the finer ones (W = 0:
                        only the squares meeting a feature are split); then
                        balanced 2:1 with --balance.
Output: the "polygons" format of wave_propagation (mesh.type = "polygons",
mesh/domain.hpp): "n_points n_cells", points "x y", cells "nv i_1 .. i_nv tag"
(0-based, counterclockwise, tag 1), all the hanging nodes as polygon vertices;
and a line of statistics (cells per depth, h_max / h_min).
Usage: generate_refined_quadtree.py --n 32 --depth 3 --feature disk:0.5,0.5,0.1 --grading two --output mesh.txt
===============================================================================
"""
import argparse
import bisect
import math


def parse_feature(text):
    kind, _, args = text.partition(":")
    v = [float(a) for a in args.split(",")] if args else []
    n = {"box": 4, "disk": 3, "line": 4, "point": 2}
    if kind not in n or len(v) != n[kind]:
        raise SystemExit(f"bad feature {text!r}: box:x0,x1,y0,y1 | disk:cx,cy,r | line:x0,y0,x1,y1 | point:cx,cy")
    return kind, v


def dist_point_box(px, py, x0, x1, y0, y1):
    dx = max(x0 - px, 0.0, px - x1)
    dy = max(y0 - py, 0.0, py - y1)
    return math.hypot(dx, dy)


def dist_point_segment(px, py, ax, ay, bx, by):
    vx, vy = bx - ax, by - ay
    L2 = vx * vx + vy * vy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - ax) * vx + (py - ay) * vy) / L2))
    return math.hypot(px - ax - t * vx, py - ay - t * vy)


def segments_intersect(a, b, c, d):
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    return o1 * o2 < 0 and o3 * o4 < 0


def distance(feature, x0, x1, y0, y1):
    """Distance from the square [x0,x1]x[y0,y1] to a feature (0 when they meet)."""
    kind, v = feature
    if kind == "point":
        return dist_point_box(v[0], v[1], x0, x1, y0, y1)
    if kind == "disk":
        return max(0.0, dist_point_box(v[0], v[1], x0, x1, y0, y1) - v[2])
    if kind == "box":
        dx = max(v[0] - x1, 0.0, x0 - v[1])
        dy = max(v[2] - y1, 0.0, y0 - v[3])
        return math.hypot(dx, dy)
    ax, ay, bx, by = v
    if dist_point_box(ax, ay, x0, x1, y0, y1) == 0.0 or dist_point_box(bx, by, x0, x1, y0, y1) == 0.0:
        return 0.0
    corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    for k in range(4):
        if segments_intersect((ax, ay), (bx, by), corners[k], corners[(k + 1) % 4]):
            return 0.0
    d = min(dist_point_segment(cx, cy, ax, ay, bx, by) for (cx, cy) in corners)
    return min(d, dist_point_box(ax, ay, x0, x1, y0, y1), dist_point_box(bx, by, x0, x1, y0, y1))


def meets_area(feature, x0, x1, y0, y1):
    """The square meets the feature with positive area (box, disk) or at all (line, point), contact excluded."""
    kind, v = feature
    e = 1e-9 * (x1 - x0)
    if kind == "box":
        return min(x1, v[1]) - max(x0, v[0]) > e and min(y1, v[3]) - max(y0, v[2]) > e
    return distance(feature, x0 + e, x1 - e, y0 + e, y1 - e) == 0.0


def build(n, depth, features, grading, rings):
    """Leaves (i, j, d): square [i, i+1] x [j, j+1] * 2^-d / n."""
    leaves = set()
    if grading == "two":
        for i in range(n):
            for j in range(n):
                s = 1.0 / n
                if any(meets_area(f, i * s, (i + 1) * s, j * s, (j + 1) * s) for f in features):
                    m = 2 ** depth
                    for a in range(m):
                        for b in range(m):
                            leaves.add((i * m + a, j * m + b, depth))
                else:
                    leaves.add((i, j, 0))
        return leaves
    stack = [(i, j, 0) for i in range(n) for j in range(n)]
    while stack:
        i, j, d = stack.pop()
        s = 1.0 / (n * 2 ** d)
        x0, y0 = i * s, j * s
        if d >= depth:
            split = False
        elif rings > 0:
            split = min(distance(f, x0, x0 + s, y0, y0 + s) for f in features) < rings * s
        else:
            split = any(meets_area(f, x0, x0 + s, y0, y0 + s) for f in features)
        if split:
            stack += [(2 * i + a, 2 * j + b, d + 1) for a in (0, 1) for b in (0, 1)]
        else:
            leaves.add((i, j, d))
    return leaves


def balance(leaves, n):
    """2:1 balance: a leaf having a neighbour two levels deeper is split, until no change."""
    leaves = set(leaves)
    changed = True
    while changed:
        changed = False
        d_max = max(d for (_, _, d) in leaves)
        N = n * 2 ** d_max

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
                    if 0 <= px < N and 0 <= py < N and not leaf_at_most(px, py, d + 1):
                        split = True
                        break
                if split:
                    break
            if split:
                leaves.remove((i, j, d))
                for a in (0, 1):
                    for b in (0, 1):
                        leaves.add((2 * i + a, 2 * j + b, d + 1))
                changed = True
    return leaves


def write(path, leaves, n, header):
    """Polygons with all their hanging nodes (integer coordinates on the finest grid)."""
    d_max = max(d for (_, _, d) in leaves)
    corners = {}
    for (i, j, d) in leaves:
        s = 2 ** (d_max - d)
        corners[(i, j, d)] = (i * s, j * s, s)
    verts = set()
    for (x, y, s) in corners.values():
        verts.update(((x, y), (x + s, y), (x + s, y + s), (x, y + s)))
    by_y, by_x = {}, {}
    for (x, y) in verts:
        by_y.setdefault(y, []).append(x)
        by_x.setdefault(x, []).append(y)
    for v in by_y.values():
        v.sort()
    for v in by_x.values():
        v.sort()
    ids = {v: k for k, v in enumerate(sorted(verts))}

    def on_h(y, xa, xb):        # vertices on the horizontal edge y, from xa to xb (exclusive of xb)
        row = by_y[y]
        lo, hi = min(xa, xb), max(xa, xb)
        pts = row[bisect.bisect_left(row, lo):bisect.bisect_right(row, hi)]
        pts = pts if xa < xb else pts[::-1]
        return [(x, y) for x in pts[:-1]]

    def on_v(x, ya, yb):
        col = by_x[x]
        lo, hi = min(ya, yb), max(ya, yb)
        pts = col[bisect.bisect_left(col, lo):bisect.bisect_right(col, hi)]
        pts = pts if ya < yb else pts[::-1]
        return [(x, y) for y in pts[:-1]]

    cells = []
    for key in sorted(leaves):
        x, y, s = corners[key]
        ring = on_h(y, x, x + s) + on_v(x + s, y, y + s) + on_h(y + s, x + s, x) + on_v(x, y + s, y)
        cells.append([ids[v] for v in ring])
    D = n * 2 ** d_max
    with open(path, "w") as f:
        f.write(f"# {header}\n")
        f.write(f"{len(ids)} {len(cells)}\n")
        for (x, y) in sorted(ids, key=ids.get):
            f.write(f"{x / D:.17g} {y / D:.17g}\n")
        for ring in cells:
            f.write(f"{len(ring)} " + " ".join(map(str, ring)) + " 1\n")
    per_depth = {}
    for (_, _, d) in leaves:
        per_depth[d] = per_depth.get(d, 0) + 1
    return len(cells), per_depth, max(len(r) for r in cells)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, required=True, help="background squares per side")
    ap.add_argument("--depth", type=int, required=True, help="refinement depth L of the features")
    ap.add_argument("--feature", action="append", required=True, help="box:.. | disk:.. | line:.. | point:..")
    ap.add_argument("--grading", choices=["two", "rings"], default="rings")
    ap.add_argument("--rings", type=float, default=1.0, help="band width W of each level (grading rings)")
    ap.add_argument("--balance", action="store_true", help="2:1 balance (grading rings)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    feats = [parse_feature(t) for t in a.feature]
    leaves = build(a.n, a.depth, feats, a.grading, a.rings)
    if a.balance:
        leaves = balance(leaves, a.n)
    header = (f"generate_refined_quadtree.py --n {a.n} --depth {a.depth} --grading {a.grading} --rings {a.rings}"
              f"{' --balance' if a.balance else ''} " + " ".join("--feature " + t for t in a.feature))
    nc, per_depth, nvmax = write(a.output, leaves, a.n, header)
    dmin, dmax = min(per_depth), max(per_depth)
    print(f"{a.output}: {nc} cells, per depth " + " ".join(f"{d}:{per_depth[d]}" for d in sorted(per_depth))
          + f", h_max / h_min = {2 ** (dmax - dmin)}, at most {nvmax} vertices per cell")


if __name__ == "__main__":
    main()
