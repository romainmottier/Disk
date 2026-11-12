#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data of the teaching animation Anim_IV.19_brain (course "Numerical methods for ODEs", Basel, session 26: Teaching
repository, BASEL/ODE/Handsout/Course/src/animations_chapter4.py, function anim_brain, which reads the .npz written
here): transcranial focused ultrasound through a titanium plate, test case sec:num-head of the LTS_HHO preprint
(main.tex).

The run is the job fit/plate_li of article/jobs.txt (input article/inputs/tfus_fit.lua, mesh.level = 2, LI, 6 threads,
about 15 min) with 200 snapshots, one every 7 steps (dt = 80 / 1394 us), i.e. the Lua statements
    mesh.level = 2
    output.snapshots = {}
    for j = 0, 199 do output.snapshots[j + 1] = 7 * j * 80.0 / 1394 end
appended to a copy of tfus_fit.lua, run in its own directory:
    OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 wave_propagation case.lua
The pressure at the target without the head and with the intact skull comes from the jobs fit/water_li and
fit/intact_li of article/jobs.txt.

Writes, on the box BOX with square pixels PX [mm]:
    field     frames of the pressure (fluids) and mean stress -tr(sigma)/2 (bone, titanium), int8, in signed square
              root: q = round(127 sign(v) sqrt(min(|v| / LIM, 1)))
    tissue    tags of the cells of the mesh (1 cortical bone, 2 cancellous bone, 3 water, 4 scalp, 5 CSF, 6 gray
              matter, 7 white matter, 8 titanium)
    t         times of the frames [us]
    sig_*     pressure at the target (geometric focus) without the head, with the intact skull and with the plate
    zoom_*    cells of the mesh around the plate (vertices, tags, implicit in LI: diameter below fine_size)
Usage: python3 make_tfus_movie.py <movie run> <fit runs> <output .npz>
"""
import os, re, sys
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
MESH = os.path.join(HERE, "..", "..", "..", "src", "mesh_generation", "brain_mri", "meshes", "head_s80_fit_plate.txt")
BOX, PX, LIM = (25.0, 236.0, 25.0, 265.0), 0.5, 3.0
ZOOM, FINE = (138.0, 164.0, 196.0, 222.0), 0.75       # box of Fig. tfus-model (right), time.fine_size of tfus_fit.lua


def pixel_map(polys, m, box, px):
    """For every pixel centre of box: the cell containing it (-1: none) and its sub-square (row, column) in the m x m
    samples of the bounding box of the cell (rows bottom to top). Axis-aligned rectangles are painted on their box,
    the other (cut) cells only on the pixels inside their polygon, after the rectangles."""
    X0, X1, Y0, Y1 = box
    nx, ny = int(round((X1 - X0) / px)), int(round((Y1 - Y0) / px))
    xc = X0 + (np.arange(nx) + 0.5) * px
    yc = Y0 + (np.arange(ny) + 0.5) * px
    cell = np.full((ny, nx), -1, dtype=np.int32)
    row = np.zeros((ny, nx), dtype=np.int8)
    col = np.zeros((ny, nx), dtype=np.int8)
    later = []
    for k, q in enumerate(polys):
        x0, x1, y0, y1 = q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()
        if x1 <= X0 or x0 >= X1 or y1 <= Y0 or y0 >= Y1:
            continue
        area = abs(0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1]))
        if abs(area - (x1 - x0) * (y1 - y0)) >= 1e-9 * (x1 - x0) * (y1 - y0):
            later.append(k)
            continue
        i0, i1 = np.searchsorted(yc, y0), np.searchsorted(yc, y1)
        j0, j1 = np.searchsorted(xc, x0), np.searchsorted(xc, x1)
        if i1 <= i0 or j1 <= j0:
            continue
        cell[i0:i1, j0:j1] = k
        row[i0:i1, j0:j1] = np.minimum(((yc[i0:i1] - y0) / (y1 - y0) * m).astype(int), m - 1)[:, None]
        col[i0:i1, j0:j1] = np.minimum(((xc[j0:j1] - x0) / (x1 - x0) * m).astype(int), m - 1)[None, :]
    for k in later:
        q = polys[k]
        x0, x1, y0, y1 = q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()
        i0, i1 = np.searchsorted(yc, y0), np.searchsorted(yc, y1)
        j0, j1 = np.searchsorted(xc, x0), np.searchsorted(xc, x1)
        if i1 <= i0 or j1 <= j0:
            continue
        XX, YY = np.meshgrid(xc[j0:j1], yc[i0:i1])
        inside = Path(q).contains_points(np.column_stack([XX.ravel(), YY.ravel()]), radius=1e-9).reshape(XX.shape)
        r = np.minimum(((yc[i0:i1] - y0) / (y1 - y0) * m).astype(int), m - 1)
        c = np.minimum(((xc[j0:j1] - x0) / (x1 - x0) * m).astype(int), m - 1)
        cell[i0:i1, j0:j1][inside] = k
        row[i0:i1, j0:j1][inside] = np.broadcast_to(r[:, None], XX.shape)[inside]
        col[i0:i1, j0:j1][inside] = np.broadcast_to(c[None, :], XX.shape)[inside]
    return cell, row, col


def read_cells(path):
    """vertices of the cells, in the order of the samples (<table>_cells.dat: medium, number of vertices, vertices)."""
    polys = []
    for line in open(path):
        if not line.startswith("#"):
            v = line.split()
            polys.append(np.array(v[2:2 + 2 * int(v[1])], dtype=float).reshape(-1, 2))
    return polys


def read_mesh(path):
    """vertices, vertex indices and tags of the cells of a polygonal mesh file (points, then cells: n, vertex
    indices, tag)."""
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    npt, nc = map(int, lines[0].split())
    pts = np.array([l.split() for l in lines[1:1 + npt]], dtype=float)
    polys, rings, tags = [], [], []
    for l in lines[1 + npt:1 + npt + nc]:
        v = list(map(int, l.split()))
        rings.append(v[1:1 + v[0]])
        polys.append(pts[rings[-1]])
        tags.append(v[1 + v[0]])
    return polys, rings, np.array(tags)


def fine_region(polys, rings, size, layers=1):
    """implicit cells of LI: diameter below time.fine_size, plus time.fine_layers layers of neighbours by vertex."""
    fine = np.array([np.max(np.linalg.norm(q[:, None, :] - q[None, :, :], axis=2)) < size for q in polys])
    for _ in range(layers):
        hit = set(v for k in np.nonzero(fine)[0] for v in rings[k])
        fine = fine | np.array([any(v in hit for v in r) for r in rings])
    return fine


def read_samples(path):
    """requested time and the m x m samples of every cell (rows bottom to top)."""
    with open(path) as f:
        head = f.readline()
        a = np.fromstring(f.read(), sep=" ")
    t = float(re.search(r"requested ([0-9.eE+-]+)", head).group(1).rstrip(":"))
    n = a.size
    for m in range(1, 9):
        if n % (4 + m * m) == 0 and (n // (4 + m * m)) * (4 + m * m) == n:
            D = a.reshape(-1, 4 + m * m)
            if np.all(D[:, 1] >= D[:, 0]) and np.all(D[:, 3] >= D[:, 2]):
                return t, D[:, 4:].reshape(-1, m, m)
    raise RuntimeError("cannot read " + path)


def read_sensor(path):
    D = np.loadtxt(path, delimiter=",", comments="#", skiprows=2)
    return D[:, 0], D[:, 1]


def main():
    movie, fit, out = sys.argv[1:4]
    files = sorted((int(m.group(1)), f) for f in os.listdir(movie)
                   for m in [re.match(r"c_snapshot_(\d+)_samples\.dat$", f)] if m)
    polys = read_cells(os.path.join(movie, "c_cells.dat"))
    t0, v0 = read_samples(os.path.join(movie, files[0][1]))
    cell, row, col = pixel_map(polys, v0.shape[1], BOX, PX)
    inside = cell >= 0
    frames, times = [], []
    for i, f in files:
        t, v = read_samples(os.path.join(movie, f))
        img = np.zeros(cell.shape)
        img[inside] = v[cell[inside], row[inside], col[inside]]
        frames.append(np.round(127 * np.sign(img) * np.sqrt(np.minimum(np.abs(img) / LIM, 1.0))).astype(np.int8))
        times.append(t)
        print(f"\r{len(times)}/{len(files)} t = {t:.2f}", end="", flush=True)
    print()
    mpolys, rings, tags = read_mesh(MESH)
    c, _, _ = pixel_map(mpolys, 1, BOX, PX)
    tissue = np.where(c >= 0, tags[np.maximum(c, 0)], 3).astype(np.uint8)
    # cells around the plate, and the implicit region of LI
    fine = fine_region(mpolys, rings, FINE)
    X0, X1, Y0, Y1 = ZOOM
    sel = [k for k, q in enumerate(mpolys) if q[:, 0].max() > X0 and q[:, 0].min() < X1 and q[:, 1].max() > Y0
           and q[:, 1].min() < Y1]
    zv = np.concatenate([mpolys[k] for k in sel]).astype(np.float32)
    zn = np.array([len(mpolys[k]) for k in sel], dtype=np.int32)
    sig = {}
    for name, d in (("water", os.path.join(fit, "water_li")), ("intact", os.path.join(fit, "intact_li")),
                    ("plate", movie)):
        ts, p = read_sensor(os.path.join(d, "c_focus.csv"))
        sig["sig_t_" + name], sig["sig_" + name] = ts.astype(np.float32), p.astype(np.float32)
    log = open(os.path.join(movie, "log.txt")).read()
    dofs = float(re.findall(r"\s([0-9.]+e\+0[0-9])\s+--\s+[0-9]+\s+--", log)[-1])
    np.savez_compressed(out, field=np.array(frames), t=np.array(times, dtype=np.float32), lim=LIM,
                        box=np.array(BOX), px=PX, tissue=tissue, zoom_box=np.array(ZOOM), zoom_vertices=zv,
                        zoom_nv=zn, zoom_tags=tags[sel].astype(np.uint8), zoom_fine=fine[sel],
                        n_cells=len(mpolys), n_fine=int(np.sum(fine)), dofs=dofs, **sig)
    print(f"{out}: {len(times)} frames {cell.shape}, {os.path.getsize(out) / 2**20:.1f} MB, "
          f"cells {len(mpolys)}, implicit in LI {int(np.sum(fine))}, unknowns {dofs:.4g}")


if __name__ == "__main__":
    main()
