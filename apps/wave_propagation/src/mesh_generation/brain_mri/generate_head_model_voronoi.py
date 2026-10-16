#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Polytopal (Voronoi) head model with fitted solid / fluid interfaces.
===============================================================================
The seeds are the centres of the cells of the quadtree of
generate_head_model_fitted.py (the same local sizes: fluids, bone, implant),
moved at random by up to --jitter x the cell size; the Voronoi cells of the
seeds (bounded by the seeds mirrored across the sides of the box) are convex
polygons with any number of vertices, conforming (shared Voronoi vertices).
They are then cut along the smoothed skull and the exact implant
(cut_polygons.py, two passes, vertices snapped within --delta x the cell size;
with --move, these vertices are moved onto the interfaces)
and tagged as in generate_head_model_fitted.py. No cell of the result is a
square: the test of genericity of the HHO discretization.
Usage:
    python3 generate_head_model_voronoi.py --output head.txt [--plate 50 66] --margin 40 --base 4 \\
        --size-soft 1 --size-bone 0.5 --interface 0.5 --size-implant 0.125 [--jitter 0.3] [--seed 1]
===============================================================================
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import Voronoi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_head_model as G                       # noqa: E402
import generate_head_model_fitted as F                # noqa: E402
import cut_polygons as C                             # noqa: E402

BONE = (G.CORTICAL, G.CANCELLOUS)


def voronoi_cells(seeds, W, H):
    """Convex Voronoi cells of the seeds in the box [0, W] x [0, H] (seeds mirrored across the four sides): points,
    counterclockwise rings (indices into points), one per seed."""
    n = len(seeds)
    mirrors = [seeds,
               np.column_stack([-seeds[:, 0], seeds[:, 1]]), np.column_stack([2 * W - seeds[:, 0], seeds[:, 1]]),
               np.column_stack([seeds[:, 0], -seeds[:, 1]]), np.column_stack([seeds[:, 0], 2 * H - seeds[:, 1]])]
    vor = Voronoi(np.vstack(mirrors))
    V = vor.vertices.copy()
    # snap the vertices on the sides of the box (round-off of the mirror construction)
    tol = 1e-9 * max(W, H)
    for c, v in ((0, 0.0), (0, W), (1, 0.0), (1, H)):
        V[np.abs(V[:, c] - v) < tol, c] = v
    rings = []
    for i in range(n):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or len(reg) < 3:
            raise RuntimeError(f"unbounded Voronoi cell of seed {i}")
        q = V[reg]
        c = q.mean(axis=0)
        ang = np.arctan2(q[:, 1] - c[1], q[:, 0] - c[0])
        rings.append([reg[k] for k in np.argsort(ang)])
    # keep the vertices used, renumbered
    used = sorted({v for r in rings for v in r})
    idx = {v: k for k, v in enumerate(used)}
    return V[used], [[idx[v] for v in r] for r in rings]


def merge_close(P, rings, tol):
    """Merge the vertices closer than tol (degenerate Voronoi edges): returns points, rings without repeated vertices."""
    from scipy.spatial import cKDTree
    tree = cKDTree(P)
    pairs = tree.query_pairs(tol)
    parent = list(range(len(P)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    out = []
    for r in rings:
        rr = []
        for v in r:
            f = find(v)
            if not rr or rr[-1] != f:
                rr.append(f)
        if len(rr) > 1 and rr[0] == rr[-1]:
            rr.pop()
        out.append(rr)
    return P, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name, typ, dflt in (("--data", str, os.path.expanduser("~/Downloads/mri_dataset/data")), ("--subject", int, 1),
                            ("--slice", int, 80), ("--pixel", float, 1.0), ("--margin", int, 12),
                            ("--head-threshold", float, 60.0), ("--diploe-threshold", float, 120.0), ("--base", float, 2.0),
                            ("--size-soft", float, 2.0), ("--size-bone", float, 1.0), ("--interface", float, 1.0),
                            ("--plate-thickness", float, 0.6), ("--screw-diameter", float, 1.5),
                            ("--screw-length", float, 4.0), ("--size-implant", float, 0.125), ("--smooth", float, 1.0),
                            ("--delta", float, 0.3), ("--jitter", float, 0.3), ("--seed", int, 1),
                            ("--min-width", float, 0.15)):
        ap.add_argument(name, type=typ, default=dflt)
    ap.add_argument("--plate", type=float, nargs=2, default=None)
    ap.add_argument("--move", action="store_true",
                    help="move the vertices within --delta x the cell size onto the interfaces (default: the cut only "
                         "passes through them, within delta of the curve)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    seg = np.load(os.path.join(a.data, "segmentation", f"sub-{a.subject:03d}_T1w.npy"))[a.slice]
    t1 = G.read_nifti(os.path.join(a.data, "mri", f"sub-{a.subject:03d}_T1w.nii.gz"))[a.slice]
    tag, _ = G.crop(G.tissues(seg, t1, a.head_threshold, a.diploe_threshold), a.margin, int(round(a.base / a.pixel)))
    H, W = tag.shape
    px = a.pixel
    smallest = min(a.interface, a.size_bone, a.size_implant if a.plate else px, px)
    scale = 2 ** max(0, int(math.ceil(math.log2(px / smallest - 1e-9))))
    unit = px / scale
    # level sets (as generate_head_model_fitted.py)
    bone = np.isin(tag, BONE).astype(float)
    phib_grid = ndi.gaussian_filter(bone, a.smooth) - 0.5
    coef = ndi.spline_filter(phib_grid, order=3)

    def smooth_bone(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        return ndi.map_coordinates(coef, [y / px - 0.5, x / px - 0.5], order=3, prefilter=False, mode="nearest")
    # exact signed distance to the zero level of the smoothed indicator: the level line on a 0.05 mm grid (contourpy),
    # densified, distances by a k-d tree, sign of the smoothed indicator (linear interpolation of phi along the edges
    # then places the crossing points where the curve is, whatever the distance of the vertices)
    import contourpy
    from scipy.spatial import cKDTree
    hc = 0.05
    xs = np.arange(0.0, W * px + hc / 2, hc); ys = np.arange(0.0, H * px + hc / 2, hc)
    XS, YS = np.meshgrid(xs, ys)
    lines = contourpy.contour_generator(xs, ys, smooth_bone(XS, YS)).lines(0.0)
    del XS, YS
    pts = []
    for ln in lines:
        for k in range(len(ln) - 1):
            a_, b_ = ln[k], ln[k + 1]
            m = max(1, int(np.ceil(np.hypot(*(b_ - a_)) / (0.2 * hc))))
            pts.append(a_ + np.linspace(0, 1, m, endpoint=False)[:, None] * (b_ - a_))
    tree = cKDTree(np.vstack(pts))

    def phi_bone(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        d, _ = tree.query(np.column_stack([x.ravel(), y.ravel()]))
        return (np.sign(smooth_bone(x, y)).ravel() * d).reshape(x.shape)
    phi_ti = None
    if a.plate:
        _, geo = G.implant_raster(tag, scale, px, a.plate, a.plate_thickness, a.screw_diameter, a.screw_length)
        xc, yc = geo["center"]; R = geo["radius"]; t = a.plate_thickness
        b0, b1 = (math.radians(v) for v in geo["angles"])
        sb = [b0 + 0.5 * a.screw_diameter / R + 1.0 / R, b1 - 0.5 * a.screw_diameter / R - 1.0 / R]

        def phi_ti(x, y, inner=R):
            x, y = np.asarray(x, float), np.asarray(y, float)
            r = np.hypot(x - xc, y - yc)
            th = np.arctan2(y - yc, x - xc)
            p = np.minimum.reduce([r - inner, R + t - r, R * (th - b0), R * (b1 - th)])
            for b in sb:
                along = (x - xc) * math.cos(b) + (y - yc) * math.sin(b)
                across = np.abs(-(x - xc) * math.sin(b) + (y - yc) * math.cos(b))
                p = np.maximum(p, np.minimum.reduce([along - (R - a.screw_length), R + t - along,
                                                     0.5 * a.screw_diameter - across]))
            return p

        def phi_solid(x, y):
            return np.maximum(phi_bone(x, y), phi_ti(x, y, inner=R - 1.5))
    else:
        phi_solid = phi_bone
    is_bone = np.isin(tag, BONE)
    _, ib = ndi.distance_transform_edt(~is_bone, return_indices=True)
    _, if_ = ndi.distance_transform_edt(is_bone, return_indices=True)
    bone_tag = tag[ib[0], ib[1]]
    fluid_tag = tag[if_[0], if_[1]]

    def tag_at(x, y, side):
        i = min(max(int(y / px), 0), H - 1); j = min(max(int(x / px), 0), W - 1)
        if side > 0:
            if phi_ti is not None and float(np.atleast_1d(phi_ti(np.atleast_1d(x), np.atleast_1d(y)))[0]) > 0:
                return G.TITANIUM
            return int(bone_tag[i, j])
        return int(fluid_tag[i, j])

    # quadtree of the smooth regions: the local sizes
    Hf, Wf = H * scale, W * scale
    X, Y = np.meshgrid((np.arange(Wf) + 0.5) * unit, (np.arange(Hf) + 0.5) * unit)
    solid = phi_solid(X, Y) > 0
    I = np.clip((Y / px).astype(int), 0, H - 1); J = np.clip((X / px).astype(int), 0, W - 1)
    fine = np.where(solid, bone_tag[I, J], fluid_tag[I, J]).astype(np.int8)
    if phi_ti is not None:
        fine[phi_ti(X, Y) > 0] = G.TITANIUM
    del X, Y, I, J
    f = lambda mm: max(1, int(round(mm / unit)))
    sizes = {t_: f(a.size_soft) for t_ in (G.WATER, G.SCALP, G.CSF, G.GRAY, G.WHITE)}
    sizes.update({G.CORTICAL: f(a.size_bone), G.CANCELLOUS: f(a.size_bone), G.TITANIUM: f(a.size_implant)})
    leaves, _ = G.build(fine, sizes, (f(a.interface), f(a.size_implant)), f(a.base))
    leaves = G.balance(leaves, fine.shape)
    # seeds: centres of the leaves, jittered
    rng = np.random.default_rng(a.seed)
    s = np.array([lf[2] for lf in leaves], float) * unit
    c = np.array([((lf[1] + 0.5 * lf[2]) * unit, (lf[0] + 0.5 * lf[2]) * unit) for lf in leaves])
    seeds = c + (rng.uniform(-1, 1, size=c.shape) * a.jitter) * s[:, None]
    seeds[:, 0] = np.clip(seeds[:, 0], 1e-3 * s, W * px - 1e-3 * s)
    seeds[:, 1] = np.clip(seeds[:, 1], 1e-3 * s, H * px - 1e-3 * s)
    P, rings = voronoi_cells(seeds, W * px, H * px)
    P, rings = merge_close(P, rings, 1e-6 * unit)
    print(f"Voronoi: {len(rings)} cells, vertices per cell {min(map(len, rings))} .. {max(map(len, rings))}")
    # cuts: solid / fluid, then titanium / bone in the solid
    if a.move:
        P, nm, nr, nc = C.move_to_level_set(P, rings, phi_solid, a.delta, min_ratio=a.min_width)
        print(f"solid / fluid: {nm} vertices moved onto the interface ({nr} kept), {nc} cells no longer convex")
    P, rings1, _, side1, st1 = F.cut_pass(P, rings, [None] * len(rings), phi_solid, a.delta, min_width=a.min_width)
    print("solid / fluid cut:", st1)
    final_rings, final_side = rings1, side1
    if phi_ti is not None:
        if a.move:
            # the vertices of the solid / fluid interface stay on it; every cell around a moved vertex is checked
            hv = np.full(len(P), np.inf)
            for ring in rings1:
                q = P[ring]; hv[ring] = np.minimum(hv[ring], np.max(np.hypot(*(np.roll(q, -1, axis=0) - q).T)))
            on_solid = np.abs(phi_solid(P[:, 0], P[:, 1])) < 1e-2 * hv
            P, nm, nr, nc = C.move_to_level_set(P, rings1,
                                                lambda x, y: phi_ti(x, y), a.delta, frozen=on_solid,
                                                min_ratio=a.min_width)
            print(f"titanium / bone: {nm} vertices moved onto the interface ({nr} kept), {nc} cells no longer convex")
        P, rings2, parent2, side2, st2 = F.cut_pass(P, rings1, side1, lambda x, y: phi_ti(x, y), a.delta, only_side=1,
                                                    min_width=a.min_width)
        print("titanium / bone cut:", st2)
        final_rings = rings2
        final_side = [1 if s_ is not None else side1[k] for k, s_ in zip(parent2, side2)]
    tags = []
    for ring, sd in zip(final_rings, final_side):
        q = P[ring]
        cx, cy = q[:, 0].mean(), q[:, 1].mean()
        if sd == 0 or sd is None:
            sd = 1 if float(np.atleast_1d(phi_solid(np.atleast_1d(cx), np.atleast_1d(cy)))[0]) > 0 else -1
        tags.append(tag_at(cx, cy, sd))
    used = sorted({v for r in final_rings for v in r})
    idx = {v: k for k, v in enumerate(used)}
    with open(a.output, "w") as fo:
        fo.write(f"# polytopal (Voronoi) head model with fitted solid interfaces, subject {a.subject}, slice {a.slice} "
                 f"(generate_head_model_voronoi.py, jitter {a.jitter}, seed {a.seed}), lengths in mm, tags: "
                 + ", ".join(f"{t_} {n}" for t_, n in G.NAMES.items()) + "\n")
        fo.write(f"{len(used)} {len(final_rings)}\n")
        for v in used:
            fo.write(f"{P[v][0]:.15g} {P[v][1]:.15g}\n")
        for ring, tg in zip(final_rings, tags):
            fo.write(f"{len(ring)} " + " ".join(str(idx[v]) for v in ring) + f" {tg}\n")
    np.save(a.output.rsplit(".", 1)[0] + "_tissues.npy", tag)
    tags = np.array(tags)
    nv = np.array([len(r) for r in final_rings])
    print(f"{a.output}: {len(used)} points, {len(final_rings)} cells, vertices per cell {nv.min()} .. {nv.max()} "
          f"(mean {nv.mean():.2f})")
    for t_ in sorted(set(tags)):
        print(f"   tag {t_} ({G.NAMES[t_]}): {np.sum(tags == t_)} cells")


if __name__ == "__main__":
    main()
