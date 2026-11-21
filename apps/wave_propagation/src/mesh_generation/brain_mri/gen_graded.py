#!/usr/bin/env python3
"""
Diagnostic G, 2026-11-16 (lts_work/head_fine): where does the plate lose its energy (0.44 % over 80 us on fit2_plate,
0.38 % with the implant cells halved, fit2i)? Hypothesis: in the cells of the 2:1-balanced quadtree that surround the
corners of the titanium (plate ends, screws), where the fields are singular: the sizes there jump from 0.0625 mm to
0.5 mm in three one-cell layers, so the large cells sit at about one cell size from the corners.
Same mesh as fit2_plate (generate_head_model_fitted.py with the options of gen_meshes.sh), but the quadtree is graded
towards the titanium: a cell at distance d from it has a side at most max(size_implant, grade x d). The generator of
the repository is not modified: its function build() is replaced here, at run time, by the same function with this
extra rule (then the usual 2:1 balance and cuts).
Usage: python3 gen_graded.py <grade> <options of generate_head_model_fitted.py, with --output>
"""
import os
import sys

import numpy as np
from scipy import ndimage as ndi

BRAIN = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri"
sys.path.insert(0, BRAIN)
import generate_head_model as G            # noqa: E402
import generate_head_model_fitted as F     # noqa: E402

grade = float(sys.argv[1])


def build_graded(fine, sizes, interface, base):
    """generate_head_model.build plus: side <= max(implant size, grade x distance to the titanium), raster units"""
    H, W = fine.shape
    dist = ndi.distance_transform_edt(fine != G.TITANIUM)           # raster units, 0 in the titanium
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
            target = min(target, interface[0] if G.TITANIUM not in kinds else interface[1])
        if s > interface[1]:
            target = min(target, max(interface[1], int(grade * dist[i:i + s, j:j + s].min())))
        if s > target and s > 1:
            h = s // 2
            stack += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]
        else:
            leaves.append((i, j, s))
    print(f"graded build (side <= max(implant, {grade:g} x distance to the titanium)): {len(leaves)} leaves")
    return leaves, counts


G.build = build_graded
sys.argv = [os.path.join(BRAIN, "generate_head_model_fitted.py")] + sys.argv[2:]
F.main()
