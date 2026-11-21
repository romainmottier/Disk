#!/usr/bin/env python3
"""
Numbers of Section 5.4 of the LTS_HHO preprint on the head meshes twice finer (2026-11-16, lts_work/head_fine):
mesh sizes (Table tfus-sizes), split cells and smallest width (Meshes), parts of MLTS and LI (Table tfus-levels), time
steps, accuracy in time (pressure at the geometric focus: LI at dt and dt/2, MLTS), energy losses (fitted, staircase,
Voronoi), energy held by the solids, focal metrics, costs (Table tfus-cost, when the timed runs exist).
Usage: python3 analyze_fine.py [<campaign dir>] [--plate P]   (default: the directory of this script, P = fit2)
--plate P (2026-11-17, plan B: P = fit2g, quadtree graded towards the titanium): the runs with the plate, the ERK(4)
stable step and the plate mesh come from out/P and meshes/head_s80_P_plate.txt (names "P/..." below); the intact and
water runs stay those of out/fit2.
Read only: prints the numbers, writes nothing.
"""
import os, re, sys
import numpy as np

ARGV = sys.argv[1:]
PLATE = "fit2"
if "--plate" in ARGV:
    _k = ARGV.index("--plate")
    PLATE, ARGV = ARGV[_k + 1], ARGV[:_k] + ARGV[_k + 2:]
W = os.path.abspath(ARGV[0] if ARGV else os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(W, "out")
MESH = os.path.join(W, "meshes")


def R(d):
    """run directory of a name: "P/<run>" is a run with the plate, in out/PLATE"""
    return PLATE + d[1:] if d.startswith("P/") else d


def log(d):
    f = os.path.join(OUT, R(d), "log.txt")
    return re.sub(r"\x1b\[[0-9;]*m", "", open(f).read()) if os.path.exists(f) else ""


def done(d):
    return os.path.exists(os.path.join(OUT, R(d), "time.txt"))


def csv(p):
    L = [l for l in open(p) if not l.startswith("#")]
    return np.array([[float(x) for x in l.split(",")] for l in L[1:] if l.strip()])


def total_unknowns(d):
    m = re.search(r"LI work per macro step \(unknowns x products\): 4 x ([0-9.e+]+)", log(d))
    if m:
        return float(m.group(1))
    m = re.search(r"MLTS work: .* per macro step = .* ERK4 steps \(([0-9]+) unknowns\)", log(d))
    return float(m.group(1)) if m else None


def cells(d):
    m = re.search(r"matrix-free operator: ([0-9]+) cells", log(d))
    return int(m.group(1)) if m else None


def read_mesh(path):
    body = [l for l in open(path) if not l.startswith("#") and l.strip()]
    nv, nc = map(int, body[0].split()[:2])
    pts = np.array([list(map(float, body[1 + i].split()[:2])) for i in range(nv)])
    cl = []
    for j in range(nc):
        v = list(map(int, body[1 + nv + j].split()))
        cl.append((v[1:1 + v[0]], v[1 + v[0]] if len(v) > 1 + v[0] else -1))
    return pts, cl


def width(P):
    w = np.inf
    for i in range(len(P)):
        a, b = P[i], P[(i + 1) % len(P)]
        e = b - a
        L = np.hypot(*e)
        if L > 0:
            w = min(w, np.max(np.abs((P - a) @ (np.array([-e[1], e[0]]) / L))))
    return w


def energy_loss(d):
    f = os.path.join(OUT, R(d), "c_energy.dat")
    if not os.path.exists(f):
        return None
    E = np.loadtxt(f, comments="#")
    return (1 - E[-1, 1] / E[0, 1]) * 100, E


def compare(a, b, interp=False):
    """max over the common times of |p_a - p_b| at the geometric focus, in % of the peak of b (interp: b linearly
    interpolated at the times of a, for runs with different time steps)"""
    X, Y = csv(os.path.join(OUT, R(a), "c_focus.csv")), csv(os.path.join(OUT, R(b), "c_focus.csv"))
    if interp:
        return np.max(np.abs(X[:, 1] - np.interp(X[:, 0], Y[:, 0], Y[:, 1]))) / np.max(np.abs(Y[:, 1])) * 100
    tX, tY = np.round(X[:, 0], 6), np.round(Y[:, 0], 6)
    _, iX, iY = np.intersect1d(tX, tY, return_indices=True)
    return np.max(np.abs(X[iX, 1] - Y[iY, 1])) / np.max(np.abs(Y[:, 1])) * 100


print(f"campaign {W}, runs with the plate: out/{PLATE} (P)\n")
# ------------------------------------------------------------------ meshes
print("== meshes (Table tfus-sizes): cells, unknowns (k = 3, setting 1E)")
for lab, d in (("quadtree split, plate", "P/plate_li"), ("quadtree split, intact", "fit2/intact_li"),
               ("water alone (intact mesh)", "fit2/water_li"), ("staircase, plate", "tfus2/plate_li"),
               ("staircase, intact", "tfus2/intact_li"), ("Voronoi split, plate", "vor2/plate_mlts")):
    n = total_unknowns(d)
    print(f"  {lab:28s} cells {cells(d)}  unknowns {n:.4g}" if n else f"  {lab:28s} (missing {d})")
for name in (f"{PLATE}_plate", "fit2_intact", "vor2_plate"):
    f = os.path.join(MESH, name + ".log")
    if os.path.exists(f):
        t = open(f).read()
        cuts = [int(x) for x in re.findall(r"'cut': ([0-9]+)", t)]
        nc = int(re.search(r"([0-9]+) points, ([0-9]+) cells", t).group(2))
        print(f"  {name}: cut cells {cuts} -> split pieces {2 * sum(cuts)} = {2 * sum(cuts) / nc * 100:.1f} % of {nc} cells")
pts, cl = read_mesh(os.path.join(MESH, f"head_s80_{PLATE}_plate.txt"))
ws = np.array([width(pts[r]) for r, _ in cl])
print(f"  {PLATE}_plate smallest width {ws.min():.4f} mm (titanium {min(w for w, (_, t) in zip(ws, cl) if t == 8):.4f})")

# ------------------------------------------------------------------ time steps, parts of MLTS and LI
print("\n== time integration")
dt = 80.0 / 2788
dterk = None
if os.path.exists(os.path.join(OUT, R("P/cfl_erk"), "c.dat")):
    c = np.loadtxt(os.path.join(OUT, R("P/cfl_erk"), "c.dat"), comments="#")
    dterk = 0.9 * c[2]
    print(f"  ERK(4): stable dt {c[2]:.5g} us (N {int(c[3])} stable / {int(c[4])} unstable over 0.5 us), production "
          f"0.9 x = {dterk:.4g} us, {80 / dterk:.0f} steps over 80 us; MLTS/LI dt = 80/2788 = {dt:.5g} us "
          f"= {dt / dterk:.1f} ERK(4) steps")
else:
    print("  ERK(4): stable step missing (P/cfl_erk)")
t = log("P/plate_mlts_w")
m = re.search(r"MLTS: ([0-9]+) levels, cells per level \(coarsest first\): ([0-9 ]+), ratios ([0-9 ]+)", t)
if m:
    nl, cpl, rat = int(m.group(1)), list(map(int, m.group(2).split())), list(map(int, m.group(3).split()))
    bands = [float(x) for x in re.findall(r"band ([0-9.e+]+)", t)]
    tot = total_unknowns("P/plate_mlts_w") or sum(bands)
    w = re.search(r"MLTS work: ([0-9.e+]+) per macro step = ([0-9.]+) ERK4 steps .*terminal ([0-9.]+) %", t)
    print(f"  MLTS: {nl} levels, total unknowns {tot:.6g}; work {w.group(2)} ERK4 steps per macro step, terminal "
          f"{w.group(3)} %" if w else "")
    div = 1
    for l in range(nl):
        if l > 0:
            div *= rat[l - 1]
        r_ = f"{dt / div / dterk:.3g}" if dterk else "?"
        print(f"    level {l}: cells {cpl[l]}  unknowns {bands[l]:.0f} ({bands[l] / tot * 100:.2g} %)  step dt/{div}  "
              f"ratio to ERK(4) {r_}")
else:
    print("  MLTS levels: missing (P/plate_mlts_w)")
t = log("P/plate_li")
m = re.search(r"LI fine system: ([0-9]+) unknowns \(([0-9]+) cell, ([0-9]+) face\).*condensed on the faces: ([0-9]+) "
              r"unknowns, nnz\(L \+ U\) ([0-9]+)", t)
fc = re.search(r"fine cells = ([0-9]+)/([0-9]+)", t)
if m and fc:
    tot = total_unknowns("P/plate_li")
    imp, sysn, nnz = int(m.group(1)), int(m.group(4)), int(m.group(5))
    ci, ct = int(fc.group(1)), int(fc.group(2))
    print(f"  LI: explicit cells {ct - ci} unknowns {tot - imp:.0f} ({(tot - imp) / tot * 100:.3g} %); implicit cells {ci} "
          f"({ci / ct * 100:.3g} % of cells) unknowns {imp} ({imp / tot * 100:.3g} %); linear system {sysn} "
          f"({sysn / tot * 100:.2g} %), nnz(L+U) {nnz:.3g}")

# ------------------------------------------------------------------ accuracy in time
print("\n== accuracy in time (pressure at the geometric focus, plate, max difference in % of the peak)")
for a, b in (("P/plate_li", "P/plate_li_dt2"), ("P/plate_mlts_w", "P/plate_li_dt2"),
             ("P/plate_mlts_w", "P/plate_li")):
    if done(a) and done(b):
        print(f"  {a} vs {b}: {compare(a, b):.4f} %")
    else:
        print(f"  {a} vs {b}: missing")

# ------------------------------------------------------------------ energy
print("\n== energy loss over 80 us (%) and energy of the solids (% of E(0))")
for d in ("fit2/water_li", "fit2/intact_li", "P/plate_li", "P/plate_li_dt2", "P/plate_mlts_w", "fit2/plate_li",
          "fit2i/plate_li", "fit2g/plate_li", "tfus2/intact_li", "tfus2/plate_li", "vor2/plate_mlts"):
    if d != R(d) and R(d) in ("fit2/plate_li", "fit2g/plate_li"):
        continue   # printed below under its own name
    r = energy_loss(d) if done(d) else None
    if r is None:
        print(f"  {d:20s} missing")
        continue
    loss, E = r
    s = f"  {d:20s} loss {loss:.3f} %"
    if E.shape[1] > 3:
        tt = E[:, 0]
        k = np.argmax(E[:, 3])
        s += f"; solids: max {E[k, 3] / E[0, 1] * 100:.1f} % at {tt[k]:.1f} us"
        for x in (44.0, 80.0):
            j = min(np.searchsorted(tt, x), len(tt) - 1)
            s += f", {E[j, 3] / E[0, 1] * 100:.1f} % at {x:.0f} us"
    print(s)

# ------------------------------------------------------------------ focus
print("\n== focus: peak |p| at the geometric focus (sensor) and ratios")
pk = {}
for d in ("fit2/water_li", "fit2/intact_li", "P/plate_li", "vor2/plate_mlts"):
    if done(d):
        A = csv(os.path.join(OUT, R(d), "c_focus.csv"))
        pk[d] = np.max(np.abs(A[:, 1]))
        print(f"  {d:18s} {pk[d]:.4f}")
if "fit2/intact_li" in pk and "P/plate_li" in pk:
    r = pk["P/plate_li"] / pk["fit2/intact_li"]
    print(f"  plate / intact at the target: {r:.4f} = {20 * np.log10(r):.2f} dB, intensity x{r * r:.3f} "
          f"({(1 - r * r) * 100:.1f} % less)")
if "vor2/plate_mlts" in pk and "P/plate_li" in pk:
    print(f"  Voronoi vs quadtree, peak at the target: {(pk['vor2/plate_mlts'] / pk['P/plate_li'] - 1) * 100:+.2f} %, "
          f"signals: {compare('vor2/plate_mlts', 'P/plate_li', True):.2f} % of the peak")
print("  (maps, shifts and dB re the water focal peak: make_tfus.py --focus-row, tfus_focus.txt)")
