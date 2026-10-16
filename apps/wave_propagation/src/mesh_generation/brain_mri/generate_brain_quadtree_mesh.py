#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Structured quadtree mesh of a head slice for the elasto-acoustic application
(skull = elastic solid, soft tissues and water = acoustic fluids), from a
labelled MRI segmentation slice of ~/Downloads/mri_dataset (create_dataset.py,
data/dataset/subject_XXXX_slice_YYYY.npy).
===============================================================================
Labels of the segmentation (create_dataset.py): 0 background, 1 outer cortical
bone, 2 cancellous bone, 3 inner cortical bone, 4 CSF, 5 gray matter, 6 white
matter, 7 CSF / scalp boundary. Mesh tags (problem.zones of the Lua input):
    1  cortical bone   (labels 1, 3)      elastic
    2  cancellous bone (label 2)          elastic
    3  water           (label 0, and the margin added around the head)
    4  soft tissues    (labels 4..7: CSF, gray and white matter, scalp)
       or, with --tissues separate: 4 CSF (4, 7), 5 gray matter, 6 white matter
    --skull homogeneous: all the bone labels (1, 2, 3) in tag 1 (homogeneous
    skull of the transcranial ultrasound benchmarks)

Mesh: square cells of a uniform background (side --base pixels) refined by
quadtree splitting
  - down to --interface pixels where a cell contains several tags,
  - down to --bone pixels in every cell containing bone and another tag (the
    thin cortical layers and the surfaces of the skull, fast medium: these
    cells set the CFL limit of a global explicit scheme, the setting of the
    local time stepping), and to --bone-interior pixels inside the cancellous
    bone (default --bone),
  - down to --roi-bone pixels in the bone cells meeting the box --roi (mm):
    a locally refined region of interest of the skull (e.g. the temporal
    window crossed by an ultrasound beam),
  - down to --inclusion-size pixels in the cells meeting the disk
    --inclusion X Y R (mm), whose cells (centre in the disk) take the tag 5
    (a small implant, e.g. the tip of a deep brain stimulation electrode),
then balanced 2:1 (a cell edge carries at most one hanging node). The tag of
a cell is the majority tag of its pixels (the pixel's tag below one pixel).
Hanging nodes are polygon vertices (cells with 5 to 8 vertices).

Output: the "polygons" format of wave_propagation (mesh.type = "polygons",
mesh/domain.hpp, tagged_polygon_mesh_reader): "n_points n_cells", points
"x y" (mm, y upwards), cells "nv i_1 .. i_nv tag" (0-based, counterclockwise);
plus a PNG of the mesh and a summary (cells per tag, h_min, h_max).

Usage:
    python3 generate_brain_quadtree_mesh.py SLICE.npy --output brain.txt \\
        [--pixel 1.0] [--margin 16] [--base 8] [--interface 2] [--bone 0.5] [--bone-interior 2]
===============================================================================
"""

import argparse
import math
import numpy as np

LABEL_TO_TAG = {0: 3, 1: 1, 2: 2, 3: 1, 4: 4, 5: 4, 6: 4, 7: 4}
LABEL_TO_TAG_SEPARATE = {0: 3, 1: 1, 2: 2, 3: 1, 4: 4, 5: 5, 6: 6, 7: 4}
HOMOGENEOUS_SKULL = {2: 1}
BONE_TAGS = (1, 2)


def pad_to(tags, margin, base):
    """Water margin around the slice, then sizes rounded up to multiples of base."""
    h, w = tags.shape
    H = int(math.ceil((h + 2 * margin) / base) * base)
    W = int(math.ceil((w + 2 * margin) / base) * base)
    out = np.full((H, W), 3, dtype=np.int64)
    r0, c0 = (H - h) // 2, (W - w) // 2
    out[r0:r0 + h, c0:c0 + w] = tags
    return out


def build_quadtree(tags, base, interface, bone, scale, bone_interior=None, roi=None, roi_bone=None, pixel=1.0,
                   inclusion=None, inclusion_size=None):
    """Leaves (i, j, s) on the fine grid of `scale` units per pixel: rows i, columns j, side s (fine units)."""
    H, W = tags.shape
    sb, si, sbo = int(base * scale), int(round(interface * scale)), max(1, int(round(bone * scale)))
    sbi = sbo if bone_interior is None else max(1, int(round(bone_interior * scale)))
    sroi = None if roi is None else max(1, int(round(roi_bone * scale)))
    Hf = H * scale

    def in_roi(i, j, s):
        # cell rows [i, i+s), columns [j, j+s) in fine units; roi = (x0, x1, y0, y1) in mm, y upwards
        if roi is None: return False
        u = pixel / scale
        x0, x1, y1, y0 = j * u, (j + s) * u, (Hf - i) * u, (Hf - i - s) * u
        return x0 < roi[1] and x1 > roi[0] and y0 < roi[3] and y1 > roi[2]

    sinc = None if inclusion is None else max(1, int(round(inclusion_size * scale)))

    def meets_inclusion(i, j, s):
        if inclusion is None: return False
        u = pixel / scale
        x0, x1, y0, y1 = j * u, (j + s) * u, (Hf - i - s) * u, (Hf - i) * u
        cx, cy, r = inclusion
        dx, dy = max(x0 - cx, 0.0, cx - x1), max(y0 - cy, 0.0, cy - y1)
        return dx * dx + dy * dy < r * r
    # pixel-level summed indicators of each tag (fast "which tags are in this cell" queries)
    present = sorted(set(np.unique(tags)))
    cum = {t: np.pad((tags == t).astype(np.int64), ((1, 0), (1, 0))).cumsum(0).cumsum(1) for t in present}

    def counts(i, j, s):
        # pixels covered by the fine-grid cell (i, j, s): rows [i/scale, (i+s)/scale)
        r0, r1 = i // scale, -(-(i + s) // scale)
        c0, c1 = j // scale, -(-(j + s) // scale)
        return {t: int(c[r1, c1] - c[r0, c1] - c[r1, c0] + c[r0, c0]) for t, c in cum.items()}

    leaves = []
    stack = [(i, j, sb) for i in range(0, H * scale, sb) for j in range(0, W * scale, sb)]
    while stack:
        i, j, s = stack.pop()
        cnt = counts(i, j, s)
        kinds = [t for t, n in cnt.items() if n > 0]
        target = s
        if len(kinds) > 1: target = si
        if any(t in BONE_TAGS for t in kinds):
            target = min(target, sbo if len(kinds) > 1 else sbi)
            if in_roi(i, j, s): target = min(target, sroi)
        if meets_inclusion(i, j, s): target = min(target, sinc)
        if s > target and s > 1:
            h = s // 2
            stack += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]
        else:
            leaves.append((i, j, s))
    return leaves, counts


def balance(leaves, shape_fine):
    """2:1 balance: a leaf with a neighbour smaller than half its side is split, until no change."""
    H, W = shape_fine
    while True:
        size = np.zeros((H, W), dtype=np.int64)
        for (i, j, s) in leaves:
            size[i:i + s, j:j + s] = s
        new, changed = [], False
        for (i, j, s) in leaves:
            if s >= 4:
                nb = []
                if i > 0: nb.append(size[i - 1, j:j + s].min())
                if i + s < H: nb.append(size[i + s, j:j + s].min())
                if j > 0: nb.append(size[i:i + s, j - 1].min())
                if j + s < W: nb.append(size[i:i + s, j + s].min())
                if nb and min(nb) < s // 2:
                    h = s // 2
                    new += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]
                    changed = True
                    continue
            new.append((i, j, s))
        leaves = new
        if not changed:
            return leaves


INCLUSION = None


def leaf_tag(tags, i, j, s, scale, counts):
    if INCLUSION is not None:
        cx, cy, r, u, Hf = INCLUSION
        x, y = (j + s / 2) * u, (Hf - i - s / 2) * u
        if (x - cx) ** 2 + (y - cy) ** 2 < r * r: return 5
    if s < scale:
        return int(tags[i // scale, j // scale])
    cnt = counts(i, j, s)
    return max(cnt.items(), key=lambda kv: (kv[1], -kv[0]))[0]


def write_polygons(path, leaves, tags, scale, counts, pixel, Hf):
    """Polygons with their hanging nodes, counterclockwise with y upwards (y = (Hf - row) * pixel / scale)."""
    unit = pixel / scale
    corners = set()
    for (i, j, s) in leaves:
        corners.update([(i, j), (i, j + s), (i + s, j), (i + s, j + s)])
    ids = {p: k for k, p in enumerate(sorted(corners))}
    cells = []
    for (i, j, s) in leaves:
        # corners in (row, col): bottom-left (i+s, j) -> bottom-right (i+s, j+s) -> top-right (i, j+s) -> top-left (i, j)
        # (rows downwards, y upwards: counterclockwise in (x, y)); hanging node = midpoint of an edge if it is a corner
        h = s // 2
        ring = [(i + s, j)]
        if s % 2 == 0 and (i + s, j + h) in ids: ring.append((i + s, j + h))
        ring.append((i + s, j + s))
        if s % 2 == 0 and (i + h, j + s) in ids: ring.append((i + h, j + s))
        ring.append((i, j + s))
        if s % 2 == 0 and (i, j + h) in ids: ring.append((i, j + h))
        ring.append((i, j))
        if s % 2 == 0 and (i + h, j) in ids: ring.append((i + h, j))
        cells.append(([ids[p] for p in ring], leaf_tag(tags, i, j, s, scale, counts)))
    pts = sorted(ids, key=ids.get)
    with open(path, "w") as f:
        f.write("# head slice quadtree mesh (generate_brain_quadtree_mesh.py), lengths in mm\n")
        f.write(f"{len(pts)} {len(cells)}\n")
        for (r, c) in pts:
            f.write(f"{c * unit:.10g} {(Hf - r) * unit:.10g}\n")
        for ring, tag in cells:
            f.write(f"{len(ring)} " + " ".join(map(str, ring)) + f" {tag}\n")
    return cells, pts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("slice", help="segmentation slice (.npy)")
    ap.add_argument("--output", required=True, help="mesh file (.txt)")
    ap.add_argument("--pixel", type=float, default=1.0, help="pixel size in mm (default 1)")
    ap.add_argument("--margin", type=int, default=16, help="water margin in pixels (default 16)")
    ap.add_argument("--base", type=int, default=8, help="background cell side in pixels, power of 2 (default 8)")
    ap.add_argument("--interface", type=float, default=2.0, help="cell side at tissue interfaces in pixels (default 2)")
    ap.add_argument("--bone", type=float, default=0.5, help="cell side at the bone interfaces in pixels (default 0.5)")
    ap.add_argument("--bone-interior", type=float, default=None,
                    help="cell side inside the cancellous bone in pixels (default: --bone)")
    ap.add_argument("--tissues", choices=["merged", "separate"], default="merged")
    ap.add_argument("--skull", choices=["layered", "homogeneous"], default="layered")
    ap.add_argument("--roi", type=float, nargs=4, default=None, metavar=("X0", "X1", "Y0", "Y1"),
                    help="box (mm) of the locally refined bone")
    ap.add_argument("--roi-bone", type=float, default=0.25, help="cell side of the bone in the box in pixels")
    ap.add_argument("--inclusion", type=float, nargs=3, default=None, metavar=("X", "Y", "R"),
                    help="disk (mm) of a small inclusion, tag 5")
    ap.add_argument("--inclusion-size", type=float, default=0.125, help="cell side in the inclusion in pixels")
    ap.add_argument("--png", default=None, help="picture of the mesh (default: output with .png)")
    a = ap.parse_args()

    seg = np.load(a.slice)
    lut = dict(LABEL_TO_TAG if a.tissues == "merged" else LABEL_TO_TAG_SEPARATE)
    if a.skull == "homogeneous":
        lut = {k: HOMOGENEOUS_SKULL.get(v, v) for k, v in lut.items()}
    tags = np.vectorize(lambda v: lut.get(int(v), 3))(seg)
    tags = pad_to(tags, a.margin, a.base)
    smallest = min(a.bone, a.interface, 1.0, a.roi_bone if a.roi else 1.0, a.inclusion_size if a.inclusion else 1.0)
    scale = max(1, int(round(1.0 / smallest)))                             # fine units per pixel
    H, W = tags.shape
    leaves, counts = build_quadtree(tags, a.base, a.interface, a.bone, scale, a.bone_interior, a.roi, a.roi_bone, a.pixel,
                                    a.inclusion, a.inclusion_size)
    global INCLUSION
    if a.inclusion: INCLUSION = (a.inclusion[0], a.inclusion[1], a.inclusion[2], a.pixel / scale, H * scale)
    leaves = balance(leaves, (H * scale, W * scale))
    cells, pts = write_polygons(a.output, leaves, tags, scale, counts, a.pixel, H * scale)

    unit = a.pixel / scale
    sizes = np.array([s * unit for (_, _, s) in leaves])
    ctags = np.array([t for _, t in cells])
    print(f"{a.output}: {len(pts)} points, {len(cells)} cells, domain {W * a.pixel:g} x {H * a.pixel:g} mm, "
          f"side {sizes.min():g} .. {sizes.max():g} mm")
    for t in sorted(set(ctags)):
        sel = ctags == t
        print(f"   tag {t}: {sel.sum()} cells, side {sizes[sel].min():g} .. {sizes[sel].max():g} mm")

    png = a.png or a.output.rsplit(".", 1)[0] + ".png"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        colors = {1: "#d95f02", 2: "#e6ab02", 3: "#a6cee3", 4: "#f4cae4", 5: "#1b1b1b", 6: "#fbe3f0"}
        P = np.array([[c * unit, (H * scale - r) * unit] for (r, c) in pts])
        polys = [P[ring] for ring, _ in cells]
        fig, ax = plt.subplots(figsize=(W / 25, H / 25))
        ax.add_collection(PolyCollection(polys, facecolors=[colors.get(t, "w") for t in ctags], edgecolors="k",
                                         linewidths=0.08))
        ax.set_xlim(0, W * a.pixel); ax.set_ylim(0, H * a.pixel); ax.set_aspect("equal")
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        fig.tight_layout(); fig.savefig(png, dpi=200)
        print(f"   picture {png}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
