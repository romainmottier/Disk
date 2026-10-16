#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Seven-tissue head model of an axial MRI slice and its quadtree mesh, for the
elasto-acoustic application of the LTS paper (transcranial ultrasound):
skull and implant = elastic solids, soft tissues and water = acoustic fluids.
===============================================================================
Data (~/Downloads/mri_dataset, create_dataset.py): the T1-weighted MRI
data/mri/sub-XXX_T1w.nii.gz (1 mm voxels) and its segmentation
data/segmentation/sub-XXX_T1w.npy (same voxel grid, array [z, y, x]; labels:
0 background, 1 outer skull surface, 2 skull, 3 inner skull surface (FSL BET),
4 CSF, 5 gray matter, 6 white matter (FSL FAST), 7 rest of the intracranial
mask). The segmentation has no scalp and a homogeneous skull; the T1 image
gives both:
  - head = T1 above --head-threshold (closed, holes filled, largest
    component) or segmented; scalp = head outside the skull;
  - skull = labels 1-3; its interior (one voxel away from the scalp and from
    the intracranial tissues) with a bright T1 signal (fatty marrow, median
    3 x 3 above --diploe-threshold) is the cancellous diploe, the rest (the
    outer and inner tables) cortical bone.
Tags (problem.zones of the Lua input):
    1 cortical bone   elastic       5 cerebrospinal fluid  acoustic
    2 cancellous bone elastic       6 gray matter          acoustic
    3 water           acoustic      7 white matter         acoustic
    4 scalp           acoustic      8 titanium (implant)   elastic
Implant (--plate): a cranial fixation plate (titanium, thickness --plate-
thickness) on the outer surface of the skull, an arc of a circle fitted to the
outer skull surface (through the centres of the surface voxels, moved out
by half a voxel) between two polar angles (degrees, around the centre of the
skull), fixed by two screws (diameter --screw-diameter, depth --screw-
length) at its ends, radial, through the plate into the bone.

Mesh: square cells of a uniform background (side --base mm), split while a
cell holds several tags and is larger than --interface (mm), or is larger
than the size of one of its tags (--size-soft for the fluids, --size-bone,
--size-implant), then balanced 2:1 (at most one hanging node per edge, the
hanging nodes being polygon vertices). The tag of a cell is the majority tag
of the raster (pixels subdivided --scale times). Output: the "polygons"
format of wave_propagation (mesh/domain.hpp, tagged_polygon_mesh_reader),
lengths in mm, y upwards; the tissue raster (.npy), a picture of the mesh and
of the tissues, and a summary (cells and sizes per tag).

Usage:
    python3 generate_head_model.py --data ~/Downloads/mri_dataset/data --subject 1 --slice 80 \\
        --output head.txt [--plate 40 75] [--base 2] [--interface 1] [--size-bone 1] [--size-implant 0.125]
===============================================================================
"""
import argparse
import gzip
import math
import os
import struct
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_brain_quadtree_mesh import balance  # noqa: E402  (2:1 balance of the first head meshes)

CORTICAL, CANCELLOUS, WATER, SCALP, CSF, GRAY, WHITE, TITANIUM = 1, 2, 3, 4, 5, 6, 7, 8
NAMES = {1: "cortical bone", 2: "cancellous bone", 3: "water", 4: "scalp", 5: "CSF", 6: "gray matter",
         7: "white matter", 8: "titanium"}
COLORS = {1: "#e8e3d3", 2: "#c9b98f", 3: "#cfe3f1", 4: "#e9b9a8", 5: "#7fb2d9", 6: "#9c9c9c", 7: "#efefef",
          8: "#3b3b3b"}
NIFTI_TYPES = {2: "u1", 4: "<i2", 8: "<i4", 16: "<f4", 64: "<f8", 256: "i1", 512: "<u2", 768: "<u4"}


def read_nifti(path):
    """Voxel array [z, y, x] of a little-endian NIfTI-1 file (.nii or .nii.gz), with its scaling."""
    raw = (gzip.open if path.endswith(".gz") else open)(path, "rb").read()
    if struct.unpack("<i", raw[:4])[0] != 348:
        raise ValueError(f"{path}: not a little-endian NIfTI-1 file")
    dim = struct.unpack("<8h", raw[40:56])
    dtype = NIFTI_TYPES[struct.unpack("<h", raw[70:72])[0]]
    offset = int(struct.unpack("<f", raw[108:112])[0])
    slope, inter = struct.unpack("<2f", raw[112:120])
    nx, ny, nz = dim[1], dim[2], dim[3]
    a = np.frombuffer(raw[offset:offset + nx * ny * nz * np.dtype(dtype).itemsize], dtype=dtype)
    a = a.reshape(nz, ny, nx).astype(np.float64)
    if slope != 0:
        a = a * slope + inter
    return a


def tissues(seg, t1, head_threshold, diploe_threshold):
    """Tissue tags of a slice (rows = y, columns = x) from the segmentation and the T1 image."""
    skull = np.isin(seg, [1, 2, 3])
    intracranial = np.isin(seg, [4, 5, 6, 7])
    head = (t1 > head_threshold) | skull | intracranial
    head = ndi.binary_closing(head, structure=np.ones((3, 3)), iterations=2)
    head = ndi.binary_fill_holes(head)
    lab, n = ndi.label(head)
    if n > 1:
        head = lab == (1 + np.argmax(ndi.sum(head, lab, range(1, n + 1))))
    tag = np.full(seg.shape, WATER, dtype=np.int8)
    tag[head] = SCALP
    tag[skull] = CORTICAL
    # diploe: skull voxels one voxel away from the other tissues, bright in T1 (fatty marrow)
    interior = ndi.binary_erosion(skull, structure=np.ones((3, 3)))
    bright = ndi.median_filter(t1, size=3) > diploe_threshold
    tag[interior & bright] = CANCELLOUS
    tag[np.isin(seg, [4, 7])] = CSF
    tag[seg == 5] = GRAY
    tag[seg == 6] = WHITE
    return tag


def crop(tag, margin, base):
    """Box of the head plus a water margin (pixels), sizes rounded up to multiples of base (pixels)."""
    rows, cols = np.nonzero(tag != WATER)
    r0, r1, c0, c1 = rows.min() - margin, rows.max() + 1 + margin, cols.min() - margin, cols.max() + 1 + margin
    H = int(math.ceil((r1 - r0) / base) * base)
    W = int(math.ceil((c1 - c0) / base) * base)
    out = np.full((H, W), WATER, dtype=np.int8)
    sr0, sc0 = max(r0, 0), max(c0, 0)
    sr1, sc1 = min(r0 + H, tag.shape[0]), min(c0 + W, tag.shape[1])
    out[sr0 - r0:sr1 - r0, sc0 - c0:sc1 - c0] = tag[sr0:sr1, sc0:sc1]
    return out, (r0, c0)


def implant_raster(tag, scale, pixel, plate, thickness, screw_d, screw_l):
    """Raster (scale x finer, rows = y upwards) of the titanium plate and screws; geometry summary."""
    H, W = tag.shape
    skull = np.isin(tag, [CORTICAL, CANCELLOUS])
    outer = skull & ~ndi.binary_erosion(skull) & ndi.binary_dilation(tag == SCALP)
    ys, xs = np.nonzero(outer)
    ys, xs = (ys + 0.5) * pixel, (xs + 0.5) * pixel
    cy, cx = np.nonzero(skull)
    cx0, cy0 = (cx.mean() + 0.5) * pixel, (cy.mean() + 0.5) * pixel
    a0, a1 = sorted(math.radians(a) for a in plate)
    ang = np.arctan2(ys - cy0, xs - cx0)
    sel = (ang > a0 - 0.1) & (ang < a1 + 0.1)
    # circle through the outer surface in the sector (algebraic least squares fit)
    A = np.column_stack([xs[sel], ys[sel], np.ones(sel.sum())])
    b = xs[sel] ** 2 + ys[sel] ** 2
    s = np.linalg.lstsq(A, b, rcond=None)[0]
    xc, yc = s[0] / 2, s[1] / 2
    R = math.sqrt(s[2] + xc ** 2 + yc ** 2) + 0.5 * pixel      # centres of the surface voxels -> outer surface
    Hf, Wf = H * scale, W * scale
    u = pixel / scale
    X, Y = np.meshgrid((np.arange(Wf) + 0.5) * u, (np.arange(Hf) + 0.5) * u)
    r = np.hypot(X - xc, Y - yc)
    th = np.arctan2(Y - yc, X - xc)
    # angles of the plate ends seen from the centre of the fitted circle
    def ang_c(a):
        return math.atan2(cy0 + R * math.sin(a) - yc, cx0 + R * math.cos(a) - xc)
    b0, b1 = sorted((ang_c(a0), ang_c(a1)))
    ti = (r >= R) & (r <= R + thickness) & (th >= b0) & (th <= b1)
    screws = []
    for b in (b0 + 0.5 * screw_d / R + 1.0 / R, b1 - 0.5 * screw_d / R - 1.0 / R):
        # screw axis: radial, from the outer face of the plate to screw_l inside the bone
        ex, ey = math.cos(b), math.sin(b)
        px, py = X - xc, Y - yc
        along = px * ex + py * ey
        across = np.abs(-px * ey + py * ex)
        ti |= (along <= R + thickness) & (along >= R - screw_l) & (across <= screw_d / 2)
        screws.append((xc + R * ex, yc + R * ey))
    geo = dict(center=(xc, yc), radius=R, angles=(math.degrees(b0), math.degrees(b1)), screws=screws,
               length=R * (b1 - b0))
    return ti, geo


def build(fine, sizes, interface, base):
    """Quadtree leaves (i, j, s) of the raster `fine` (rows i, columns j, side s in raster units)."""
    H, W = fine.shape
    present = sorted(int(t) for t in np.unique(fine))
    cum = {t: np.pad((fine == t).astype(np.int32), ((1, 0), (1, 0))).cumsum(0).cumsum(1) for t in present}

    def counts(i, j, s):
        return {t: int(c[i + s, j + s] - c[i, j + s] - c[i + s, j] + c[i, j]) for t, c in cum.items()}

    leaves = []
    stack = [(i, j, base) for i in range(0, H, base) for j in range(0, W, base)]
    while stack:
        i, j, s = stack.pop()
        cnt = counts(i, j, s)
        kinds = [t for t, n in cnt.items() if n > 0]
        target = min(sizes[t] for t in kinds)
        if len(kinds) > 1:
            target = min(target, interface[0] if TITANIUM not in kinds else interface[1])
        if s > target and s > 1:
            h = s // 2
            stack += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]
        else:
            leaves.append((i, j, s))
    return leaves, counts


def write_polygons(path, leaves, counts, unit, header):
    """Counterclockwise polygons (rows = y upwards here), hanging nodes as vertices; tag = majority tag."""
    corners = set()
    for (i, j, s) in leaves:
        corners.update([(i, j), (i, j + s), (i + s, j), (i + s, j + s)])
    ids = {p: k for k, p in enumerate(sorted(corners))}
    cells = []
    for (i, j, s) in leaves:
        # (row, col) with rows upwards: (i, j) bottom-left -> (i, j+s) bottom-right -> (i+s, j+s) -> (i+s, j)
        h = s // 2
        ring = [(i, j)]
        if s % 2 == 0 and (i, j + h) in ids: ring.append((i, j + h))
        ring.append((i, j + s))
        if s % 2 == 0 and (i + h, j + s) in ids: ring.append((i + h, j + s))
        ring.append((i + s, j + s))
        if s % 2 == 0 and (i + s, j + h) in ids: ring.append((i + s, j + h))
        ring.append((i + s, j))
        if s % 2 == 0 and (i + h, j) in ids: ring.append((i + h, j))
        cnt = counts(i, j, s)
        tag = max(cnt.items(), key=lambda kv: (kv[1], -kv[0]))[0]
        cells.append(([ids[p] for p in ring], tag, s))
    pts = sorted(ids, key=ids.get)
    with open(path, "w") as f:
        f.write(header)
        f.write(f"{len(pts)} {len(cells)}\n")
        for (r, c) in pts:
            f.write(f"{c * unit:.10g} {r * unit:.10g}\n")
        for ring, tag, _ in cells:
            f.write(f"{len(ring)} " + " ".join(map(str, ring)) + f" {tag}\n")
    return cells, pts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.expanduser("~/Downloads/mri_dataset/data"))
    ap.add_argument("--subject", type=int, default=1)
    ap.add_argument("--slice", type=int, default=80, help="index of the slice (first array axis)")
    ap.add_argument("--output", required=True, help="mesh file (.txt)")
    ap.add_argument("--pixel", type=float, default=1.0, help="voxel size in mm (default 1)")
    ap.add_argument("--margin", type=int, default=12, help="water around the head in pixels (default 12)")
    ap.add_argument("--head-threshold", type=float, default=60.0, help="T1 level of the head (default 60)")
    ap.add_argument("--diploe-threshold", type=float, default=120.0, help="T1 level of the diploe (default 120)")
    ap.add_argument("--base", type=float, default=2.0, help="background cell side in mm (default 2)")
    ap.add_argument("--size-soft", type=float, default=2.0, help="largest cell side in the fluids in mm (default 2)")
    ap.add_argument("--size-bone", type=float, default=1.0, help="largest cell side in the bone in mm (default 1)")
    ap.add_argument("--interface", type=float, default=1.0, help="cell side at the tissue interfaces in mm (default 1)")
    ap.add_argument("--plate", type=float, nargs=2, default=None, metavar=("A0", "A1"),
                    help="titanium plate between these polar angles (degrees, around the centre of the skull)")
    ap.add_argument("--plate-thickness", type=float, default=0.6, help="plate thickness in mm (default 0.6)")
    ap.add_argument("--screw-diameter", type=float, default=1.5, help="screw diameter in mm (default 1.5)")
    ap.add_argument("--screw-length", type=float, default=4.0, help="screw depth below the plate in mm (default 4)")
    ap.add_argument("--size-implant", type=float, default=0.125, help="cell side in the implant in mm (default 0.125)")
    a = ap.parse_args()

    seg = np.load(os.path.join(a.data, "segmentation", f"sub-{a.subject:03d}_T1w.npy"))[a.slice]
    t1 = read_nifti(os.path.join(a.data, "mri", f"sub-{a.subject:03d}_T1w.nii.gz"))[a.slice]
    tag, (r0, c0) = crop(tissues(seg, t1, a.head_threshold, a.diploe_threshold), a.margin,
                         int(round(a.base / a.pixel)))
    smallest = min(a.interface, a.size_bone, a.size_implant if a.plate else a.pixel, a.pixel)
    scale = 2 ** max(0, int(math.ceil(math.log2(a.pixel / smallest - 1e-9))))      # raster units per pixel
    unit = a.pixel / scale
    fine = np.kron(tag, np.ones((scale, scale), dtype=np.int8))
    geo = None
    if a.plate:
        ti, geo = implant_raster(tag, scale, a.pixel, a.plate, a.plate_thickness, a.screw_diameter, a.screw_length)
        fine[ti] = TITANIUM
    f = lambda mm: max(1, int(round(mm / unit)))
    sizes = {t: f(a.size_soft) for t in (WATER, SCALP, CSF, GRAY, WHITE)}
    sizes.update({CORTICAL: f(a.size_bone), CANCELLOUS: f(a.size_bone), TITANIUM: f(a.size_implant)})
    leaves, counts = build(fine, sizes, (f(a.interface), f(a.size_implant)), f(a.base))
    leaves = balance(leaves, fine.shape)
    header = (f"# head model, subject {a.subject}, slice {a.slice} (generate_head_model.py), lengths in mm, tags: "
              + ", ".join(f"{t} {n}" for t, n in NAMES.items()) + "\n")
    cells, pts = write_polygons(a.output, leaves, counts, unit, header)
    stem = a.output.rsplit(".", 1)[0]
    np.save(stem + "_tissues.npy", tag)

    sizes_mm = np.array([s * unit for (_, _, s) in cells])
    ctags = np.array([t for _, t, _ in cells])
    Hmm, Wmm = tag.shape[0] * a.pixel, tag.shape[1] * a.pixel
    print(f"{a.output}: {len(pts)} points, {len(cells)} cells, domain {Wmm:g} x {Hmm:g} mm "
          f"(slice pixels: rows {r0}.., columns {c0}..), side {sizes_mm.min():g} .. {sizes_mm.max():g} mm")
    for t in sorted(set(ctags)):
        sel = ctags == t
        print(f"   tag {t} ({NAMES[t]}): {sel.sum()} cells, side {sizes_mm[sel].min():g} .. {sizes_mm[sel].max():g} mm, "
              f"area {np.sum(sizes_mm[sel] ** 2):.0f} mm^2")
    if geo:
        print(f"   plate: circle centre ({geo['center'][0]:.2f}, {geo['center'][1]:.2f}), radius {geo['radius']:.2f} mm, "
              f"angles {geo['angles'][0]:.1f} .. {geo['angles'][1]:.1f} deg, length {geo['length']:.1f} mm, screws at "
              + ", ".join(f"({x:.1f}, {y:.1f})" for x, y in geo["screws"]))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        P = np.array([[c * unit, r * unit] for (r, c) in pts])
        fig, ax = plt.subplots(figsize=(Wmm / 20, Hmm / 20))
        ax.add_collection(PolyCollection([P[ring] for ring, _, _ in cells], facecolors=[COLORS[t] for t in ctags],
                                         edgecolors="k", linewidths=0.05))
        ax.set_xlim(0, Wmm); ax.set_ylim(0, Hmm); ax.set_aspect("equal")
        fig.tight_layout(); fig.savefig(stem + ".png", dpi=250)
        print(f"   picture {stem}.png")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
