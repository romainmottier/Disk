#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Head model of generate_head_model.py with fitted solid / fluid interfaces.
===============================================================================
The quadtree of generate_head_model.py follows the voxels: the boundaries of
the skull are staircases, whose re-entrant corners carry singular fields that
the least-squares stabilization dissipates (energy loss of order h; tilted
water / bone / water slab, 0.5 mm steps: 9 % of the energy of a wave crossing
it, 0.25 mm: 3.6 %). Here the boundaries of the solids are smooth curves:
  - bone: zero level of the bone indicator of the 1 mm voxels smoothed by a
    Gaussian of --smooth voxels (minus 1/2), cubic interpolation;
  - titanium plate and screws: their exact geometry (arc of a circle of
    thickness --plate-thickness, rectangles), the plate resting on the bone
    (the gap between its inner face and the smoothed bone is bone);
and the cells of the quadtree (built on the raster of these smooth regions)
crossing them are cut along the zero level, linearly interpolated along each
edge (cut_polygons.py): each cut cell gives two convex polygons, the new
vertices are shared by the neighbours. Two passes: solid / fluid, then
titanium / bone inside the solid. Tags: titanium inside the implant, else the
bone type (cortical or cancellous) of the nearest bone voxel; in the fluid,
the tag of the nearest fluid voxel (water, scalp, CSF, gray or white matter).
The tags and the options are those of generate_head_model.py, plus --smooth
and --delta (|phi| >= delta x the local cell size at the vertices: no
sliver).
Usage:
    python3 generate_head_model_fitted.py --output head.txt [--plate 50 66] --margin 40 --base 4 \\
        --size-soft 1 --size-bone 0.5 --interface 0.5 --size-implant 0.125
===============================================================================
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_head_model as G           # noqa: E402
from cut_polygons import cut_mesh          # noqa: E402

BONE = (G.CORTICAL, G.CANCELLOUS)


def leaves_to_polygons(leaves, unit):
    """Points (mm) and counterclockwise rings of the quadtree leaves (i, j, s) (rows = y), hanging nodes included
    (2:1 balanced: at most one per edge)."""
    corners = set()
    for (i, j, s) in leaves:
        corners.update([(i, j), (i, j + s), (i + s, j), (i + s, j + s)])
    ids = {p: k for k, p in enumerate(sorted(corners))}
    rings = []
    for (i, j, s) in leaves:
        h = s // 2
        ring = [(i, j)]
        if s % 2 == 0 and (i, j + h) in ids: ring.append((i, j + h))
        ring.append((i, j + s))
        if s % 2 == 0 and (i + h, j + s) in ids: ring.append((i + h, j + s))
        ring.append((i + s, j + s))
        if s % 2 == 0 and (i + s, j + h) in ids: ring.append((i + s, j + h))
        ring.append((i + s, j))
        if s % 2 == 0 and (i + h, j) in ids: ring.append((i + h, j))
        rings.append([ids[p] for p in ring])
    pts = np.array([(c * unit, r * unit) for (r, c) in sorted(ids, key=ids.get)])
    return pts, rings


def insert_crossings(rings, edge_pt):
    """Insert into every ring the crossing points (edge_pt[(a, b)], a < b) lying on its edges: conforming neighbours."""
    out = []
    for ring in rings:
        r = []
        n = len(ring)
        for k in range(n):
            a, b = ring[k], ring[(k + 1) % n]
            r.append(a)
            c = edge_pt.get((min(a, b), max(a, b)))
            if c is not None and c != a and c != b and c not in ring:
                r.append(c)
        out.append(r)
    return out


def cut_pass(points, rings, sides, phi_at, delta, only_side=None, min_width=0.0):
    """One cut along phi = 0 of the rings (only those with side == only_side if given); returns points, rings, sides,
    the side of each piece w.r.t. phi (+1, -1, 0: kept whole) and the crossing points by edge."""
    sel = [k for k in range(len(rings)) if only_side is None or sides[k] == only_side]
    P, pieces, stats, edge_pt = cut_mesh(points, [(rings[k], k) for k in sel], phi_at, delta, min_width)
    new_rings, new_parent, new_side = [], [], []
    done = set(sel)
    for ring, k, s in pieces:
        new_rings.append(ring); new_parent.append(k); new_side.append(s)
    for k in range(len(rings)):
        if k not in done:
            new_rings.append(rings[k]); new_parent.append(k); new_side.append(None)
    new_rings = insert_crossings(new_rings, edge_pt)
    return P, new_rings, new_parent, new_side, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.expanduser("~/Downloads/mri_dataset/data"))
    ap.add_argument("--subject", type=int, default=1)
    ap.add_argument("--slice", type=int, default=80)
    ap.add_argument("--output", required=True)
    ap.add_argument("--pixel", type=float, default=1.0)
    ap.add_argument("--margin", type=int, default=12)
    ap.add_argument("--head-threshold", type=float, default=60.0)
    ap.add_argument("--diploe-threshold", type=float, default=120.0)
    ap.add_argument("--base", type=float, default=2.0)
    ap.add_argument("--size-soft", type=float, default=2.0)
    ap.add_argument("--size-bone", type=float, default=1.0)
    ap.add_argument("--interface", type=float, default=1.0)
    ap.add_argument("--plate", type=float, nargs=2, default=None, metavar=("A0", "A1"))
    ap.add_argument("--plate-thickness", type=float, default=0.6)
    ap.add_argument("--screw-diameter", type=float, default=1.5)
    ap.add_argument("--screw-length", type=float, default=4.0)
    ap.add_argument("--size-implant", type=float, default=0.125)
    ap.add_argument("--smooth", type=float, default=1.0, help="Gaussian smoothing of the bone indicator (voxels)")
    ap.add_argument("--delta", type=float, default=0.2, help="no cut closer than delta x the cell size to a vertex")
    a = ap.parse_args()

    seg = np.load(os.path.join(a.data, "segmentation", f"sub-{a.subject:03d}_T1w.npy"))[a.slice]
    t1 = G.read_nifti(os.path.join(a.data, "mri", f"sub-{a.subject:03d}_T1w.nii.gz"))[a.slice]
    tag, (r0, c0) = G.crop(G.tissues(seg, t1, a.head_threshold, a.diploe_threshold), a.margin, int(round(a.base / a.pixel)))
    H, W = tag.shape
    px = a.pixel
    smallest = min(a.interface, a.size_bone, a.size_implant if a.plate else px, px)
    scale = 2 ** max(0, int(math.ceil(math.log2(px / smallest - 1e-9))))
    unit = px / scale

    # level set of the bone: smoothed indicator of the bone voxels, cubic interpolation (phi > 0 inside)
    bone = np.isin(tag, BONE).astype(float)
    phib_grid = ndi.gaussian_filter(bone, a.smooth) - 0.5
    coef = ndi.spline_filter(phib_grid, order=3)

    # its gradient (per mm): phi / |grad phi| is a signed distance near the zero level (snapping and sizes in mm)
    gy, gx = np.gradient(phib_grid, px)
    cgx, cgy = ndi.spline_filter(gx, order=3), ndi.spline_filter(gy, order=3)

    def phi_bone(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        c = [y / px - 0.5, x / px - 0.5]
        v = ndi.map_coordinates(coef, c, order=3, prefilter=False, mode="nearest")
        g = np.hypot(ndi.map_coordinates(cgx, c, order=3, prefilter=False, mode="nearest"),
                     ndi.map_coordinates(cgy, c, order=3, prefilter=False, mode="nearest"))
        return v / np.maximum(g, 0.05)

    geo = None
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
        phi_ti = None
        phi_solid = phi_bone

    # nearest bone / fluid voxel tags (for the pieces and the raster)
    is_bone = np.isin(tag, BONE)
    _, ib = ndi.distance_transform_edt(~is_bone, return_indices=True)
    _, if_ = ndi.distance_transform_edt(is_bone, return_indices=True)
    bone_tag = tag[ib[0], ib[1]]
    fluid_tag = tag[if_[0], if_[1]]

    def tag_at(x, y, side):
        i = np.clip((np.asarray(y) / px).astype(int), 0, H - 1)
        j = np.clip((np.asarray(x) / px).astype(int), 0, W - 1)
        if side > 0:
            if phi_ti is not None and float(np.atleast_1d(phi_ti(np.atleast_1d(x), np.atleast_1d(y)))[0]) > 0:
                return G.TITANIUM
            return int(bone_tag[i, j])
        return int(fluid_tag[i, j])

    # raster of the smooth regions (quadtree refinement only)
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
    leaves, counts = G.build(fine, sizes, (f(a.interface), f(a.size_implant)), f(a.base))
    leaves = G.balance(leaves, fine.shape)
    pts, rings = leaves_to_polygons(leaves, unit)
    print(f"quadtree: {len(rings)} cells")

    # pass 1: solid / fluid
    P, rings1, _, side1, st1 = cut_pass(pts, rings, [None] * len(rings), phi_solid, a.delta)
    print("solid / fluid cut:", st1)
    final_rings, final_side = rings1, side1
    if phi_ti is not None:
        P, rings2, parent2, side2, st2 = cut_pass(P, rings1, side1, lambda x, y: phi_ti(x, y), a.delta, only_side=1)
        print("titanium / bone cut:", st2)
        # side w.r.t. the solid: the pieces of the second pass are solid
        final_rings = rings2
        final_side = [1 if s is not None else side1[k] for k, s in zip(parent2, side2)]
    tags = []
    for ring, s in zip(final_rings, final_side):
        q = P[ring]
        cx, cy = q[:, 0].mean(), q[:, 1].mean()
        if s == 0 or s is None:
            s = 1 if float(np.atleast_1d(phi_solid(np.atleast_1d(cx), np.atleast_1d(cy)))[0]) > 0 else -1
        tags.append(tag_at(cx, cy, s))
    # drop the unused points, write
    used = sorted({v for r in final_rings for v in r})
    idx = {v: k for k, v in enumerate(used)}
    header = (f"# head model with fitted solid interfaces, subject {a.subject}, slice {a.slice} "
              f"(generate_head_model_fitted.py), lengths in mm, tags: " + ", ".join(f"{t_} {n}" for t_, n in G.NAMES.items()))
    with open(a.output, "w") as fo:
        fo.write(header + "\n")
        fo.write(f"{len(used)} {len(final_rings)}\n")
        for v in used:
            fo.write(f"{P[v][0]:.15g} {P[v][1]:.15g}\n")
        for ring, tg in zip(final_rings, tags):
            fo.write(f"{len(ring)} " + " ".join(str(idx[v]) for v in ring) + f" {tg}\n")
    stem = a.output.rsplit(".", 1)[0]
    np.save(stem + "_tissues.npy", tag)
    tags = np.array(tags)
    print(f"{a.output}: {len(used)} points, {len(final_rings)} cells, domain {W * px:g} x {H * px:g} mm")
    for t_ in sorted(set(tags)):
        print(f"   tag {t_} ({G.NAMES[t_]}): {np.sum(tags == t_)} cells")
    nv = np.array([len(r) for r in final_rings])
    print(f"   vertices per cell: {nv.min()} .. {nv.max()}")
    if geo:
        print(f"   plate: circle centre ({geo['center'][0]:.2f}, {geo['center'][1]:.2f}), radius {geo['radius']:.2f} mm, "
              f"angles {geo['angles'][0]:.1f} .. {geo['angles'][1]:.1f} deg")


if __name__ == "__main__":
    main()
