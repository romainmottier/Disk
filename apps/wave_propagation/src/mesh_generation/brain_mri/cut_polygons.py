#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitted interfaces by cutting the cells of a polygonal mesh along the zero level of a function phi (phi > 0: inside).
The vertices with |phi| < delta x the local cell size are snapped onto the interface (phi = 0: the cut passes through
them, no sliver); along an edge whose ends have strictly opposite signs, the crossing point interpolates phi linearly
(both neighbours see the same points: the edges between consecutive vertices, hanging nodes included). A cell with
exactly two cut points (snapped vertices or crossings) not adjacent along its boundary is split along the segment
joining them into two convex polygons (a straight cut of a convex polygon); a cell whose vertices are all on one side
(or on the interface) keeps that side; the others (more cut points: features thinner than the cell) are kept whole and
reported.
"""
import numpy as np


def width(q):
    """Smallest width (altitude) of a convex polygon: the stable step of an HHO cell scales with it."""
    best = np.inf
    n = len(q)
    for i in range(n):
        t = q[(i + 1) % n] - q[i]
        L = np.hypot(*t)
        if L < 1e-14:
            continue
        best = min(best, np.max(np.abs((q - q[i]) @ (np.array([-t[1], t[0]]) / L))))
    return best


def cut_mesh(points, cells, phi_at, delta=0.2, min_width=0.0):
    """points: array (n, 2); cells: list of (ring (vertex indices, counterclockwise), data); phi_at(x, y) -> phi
    (vectorized). Returns (points, pieces, stats, edge_pt): pieces = list of (ring, data, side) with side +1 (phi > 0),
    -1, or 0 (kept whole: more than two cut points); edge_pt[(a, b)] (a < b): the crossing point of the edge (a, b)
    (the neighbours of a cell must insert it into their rings)."""
    P = np.asarray(points, dtype=float)
    phi = np.asarray(phi_at(P[:, 0], P[:, 1]), dtype=float)
    # local size: the smallest, over the cells around a vertex, of their longest edge (the side of the square a cell
    # comes from, also for a piece of a cell cut by a previous pass)
    h = np.full(len(P), np.inf)
    for ring, _ in cells:
        q = P[ring]
        e = np.max(np.hypot(*(np.roll(q, -1, axis=0) - q).T))
        h[ring] = np.minimum(h[ring], e)
    sgn = np.sign(phi)
    sgn[np.abs(phi) < delta * h] = 0                    # snapped onto the interface
    new_pts = list(map(tuple, P))
    edge_pt = {}

    def crossing(a, b):
        key = (min(a, b), max(a, b))
        if key not in edge_pt:
            t = phi[a] / (phi[a] - phi[b])
            x = P[a] + t * (P[b] - P[a])
            edge_pt[key] = len(new_pts)
            new_pts.append((x[0], x[1]))
        return edge_pt[key]

    pieces = []
    stats = dict(cut=0, whole_in=0, whole_out=0, multi=0)
    for ring, data in cells:
        s = sgn[ring]
        if np.all(s >= 0) and np.any(s > 0):
            pieces.append((list(ring), data, 1)); stats["whole_in"] += 1; continue
        if np.all(s <= 0) and np.any(s < 0):
            pieces.append((list(ring), data, -1)); stats["whole_out"] += 1; continue
        if not np.any(s):                                # every vertex on the interface: tiny cell, keep
            pieces.append((list(ring), data, 0)); stats["multi"] += 1; continue
        # augmented ring: vertices and crossing points of the edges with strictly opposite signs; a run of snapped
        # vertices is a cut point (its first vertex) only if the signs before and after it differ (else the interface
        # only touches the cell there)
        aug, cut = [], []
        n = len(ring)
        nz = [k for k in range(n) if s[k] != 0]
        for k in range(n):
            a, b = ring[k], ring[(k + 1) % n]
            aug.append(a)
            if s[k] == 0 and s[(k - 1) % n] != 0:
                # first vertex of a run of zeros: compare the nonzero signs around the run
                before = s[(k - 1) % n]
                j = k
                while s[j % n] == 0:
                    j += 1
                if before != s[j % n]:
                    cut.append(len(aug) - 1)
            if s[k] * s[(k + 1) % n] < 0:
                aug.append(crossing(a, b))
                cut.append(len(aug) - 1)
        m = len(aug)
        if len(cut) != 2 or (cut[1] - cut[0]) % m in (1, m - 1):
            # more cut points, or the two are neighbours along the boundary (no area on one side): keep whole, its
            # crossing points are inserted by insert_crossings
            pieces.append((list(ring), data, 0)); stats["multi"] += 1; continue
        i, j = cut
        poly_a = aug[i:j + 1]
        poly_b = aug[j:] + aug[:i + 1]
        if min_width > 0:
            # a piece thinner than min_width x the cell size (the cut segment, through vertices snapped within delta of
            # the curve, may pass close to another vertex): keep the cell whole
            Q = np.array(new_pts)
            hcell = np.max(np.hypot(*(np.roll(P[ring], -1, axis=0) - P[ring]).T))
            if min(width(Q[poly_a]), width(Q[poly_b])) < min_width * hcell:
                pieces.append((list(ring), data, 0)); stats["multi"] += 1; continue
        stats["cut"] += 1
        side_a = 1 if np.any(np.array([sgn[v] for v in poly_a if v < len(P)]) > 0) else -1
        pieces.append((poly_a, data, side_a))
        pieces.append((poly_b, data, -side_a))
    return np.array(new_pts), pieces, stats, edge_pt


def _valid(q, min_ratio, convex=False):
    """Simple polygon (no crossing edges), counterclockwise, (convex,) and 4|T|/|dT| >= min_ratio x its diameter."""
    n = len(q)
    A = 0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1])
    if A <= 0:
        return False
    if convex:
        d = np.roll(q, -1, axis=0) - q
        cr = d[:, 0] * np.roll(d, -1, axis=0)[:, 1] - d[:, 1] * np.roll(d, -1, axis=0)[:, 0]
        if np.any(cr < -1e-12 * np.max(np.abs(cr))):
            return False
    o = lambda p, r, s_: (r[0] - p[0]) * (s_[1] - p[1]) - (r[1] - p[1]) * (s_[0] - p[0])
    for i in range(n):
        a, b = q[i], q[(i + 1) % n]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            c, d = q[j], q[(j + 1) % n]
            if o(a, b, c) * o(a, b, d) < 0 and o(c, d, a) * o(c, d, b) < 0:
                return False
    per = np.sum(np.hypot(*(np.roll(q, -1, axis=0) - q).T))
    diam = max(np.hypot(*(q[i] - q[j])) for i in range(n) for j in range(i))
    return 4 * A / per >= min_ratio * diam


def move_to_level_set(points, cells, phi_at, delta, frozen=None, iters=4, min_ratio=0.15):
    """Move the vertices with |phi| < delta x the local cell size onto the zero level of phi (Newton steps along the
    gradient, by centred differences), so that a cut through them follows the curve exactly instead of passing within
    delta of it. A vertex is only moved if all its cells stay convex polygons (the cuts assume convex cells) with
    4|T|/|dT| >= min_ratio x their diameter (vertices taken by increasing |phi|); otherwise it stays (the cut then passes through it, within delta of
    the curve, as without moving). cells: list of rings (or of (ring, data)); frozen: boolean mask of the vertices that
    must not move (e.g. on another interface). The cells sharing a vertex move with it (conforming). Returns (points,
    number moved, number rejected, number of cells no longer convex)."""
    P = np.array(points, dtype=float)
    rings = [c[0] if isinstance(c, tuple) else c for c in cells]
    h = np.full(len(P), np.inf)
    around = {}
    for k, ring in enumerate(rings):
        q = P[ring]
        h[ring] = np.minimum(h[ring], np.max(np.hypot(*(np.roll(q, -1, axis=0) - q).T)))
        for v in ring:
            around.setdefault(v, []).append(k)
    f = np.asarray(phi_at(P[:, 0], P[:, 1]), dtype=float)
    sel = (np.abs(f) < delta * h) & np.isfinite(h)
    if frozen is not None:
        sel &= ~frozen
    idx = np.nonzero(sel)[0]
    X = P[idx].copy(); e = 1e-4 * h[idx]
    for _ in range(iters):
        fx = np.asarray(phi_at(X[:, 0], X[:, 1]), dtype=float)
        gx = (np.asarray(phi_at(X[:, 0] + e, X[:, 1])) - np.asarray(phi_at(X[:, 0] - e, X[:, 1]))) / (2 * e)
        gy = (np.asarray(phi_at(X[:, 0], X[:, 1] + e)) - np.asarray(phi_at(X[:, 0], X[:, 1] - e))) / (2 * e)
        g2 = np.maximum(gx * gx + gy * gy, 1e-12)
        X[:, 0] -= fx * gx / g2; X[:, 1] -= fx * gy / g2
    moved = rejected = 0
    for k in np.argsort(np.abs(f[idx])):
        v = idx[k]; old = P[v].copy(); P[v] = X[k]
        if all(_valid(P[rings[c]], min_ratio, convex=True) for c in around[v]):
            moved += 1
        else:
            P[v] = old; rejected += 1
    nonconvex = 0
    for ring in rings:
        q = P[ring]; d = np.roll(q, -1, axis=0) - q
        cr = d[:, 0] * np.roll(d, -1, axis=0)[:, 1] - d[:, 1] * np.roll(d, -1, axis=0)[:, 0]
        nonconvex += bool(np.any(cr < -1e-12 * np.max(np.abs(cr))))
    return P, moved, rejected, nonconvex

def write_polygons(path, points, pieces, header):
    """polygons format (mesh.type = "polygons"): pieces = list of (ring, tag)."""
    used = sorted({v for ring, _ in pieces for v in ring})
    idx = {v: k for k, v in enumerate(used)}
    with open(path, "w") as f:
        f.write("# " + header + "\n")
        f.write(f"{len(used)} {len(pieces)}\n")
        for v in used:
            f.write(f"{points[v][0]:.15g} {points[v][1]:.15g}\n")
        for ring, tag in pieces:
            f.write(f"{len(ring)} " + " ".join(str(idx[v]) for v in ring) + f" {tag}\n")
