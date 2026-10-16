#!/usr/bin/env python3
"""What sets the stable step: for each mesh, the smallest diameter, face, 4|T|/|dT|, |T|/diam and minimal width
(smallest altitude of a convex polygon over its edges) of its cells; stable steps of runs/k<k>_<mesh>/c.dat.
Usage: analyze_cfl.py <out of run_cfl_cells.sh> (writes <out>/cfl.json)"""
import glob, os, re, json, sys
import numpy as np

os.chdir(sys.argv[1] if len(sys.argv) > 1 else ".")

def read(path):
    L = [l for l in open(path) if not l.startswith("#")]
    npt, nc = map(int, L[0].split())
    P = np.array([list(map(float, l.split())) for l in L[1:1 + npt]])
    cells = [list(map(int, l.split()))[1:-1] for l in L[1 + npt:1 + npt + nc]]
    return P, cells

def width(q):
    best = np.inf
    n = len(q)
    for i in range(n):
        a, b = q[i], q[(i + 1) % n]
        t = b - a; L = np.hypot(*t)
        if L < 1e-14: continue
        nrm = np.array([-t[1], t[0]]) / L
        best = min(best, np.max(np.abs((q - a) @ nrm)))
    return best

out = []
for m in sorted(glob.glob("meshes/*.txt")):
    b = os.path.basename(m)[:-4]
    P, cells = read(m)
    mm = dict(diam=np.inf, face=np.inf, ap=np.inf, ad=np.inf, w=np.inf)
    for r in cells:
        q = P[r]
        A = 0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1])
        e = np.hypot(*(np.roll(q, -1, axis=0) - q).T)
        d = max(np.hypot(*(q[i] - q[j])) for i in range(len(q)) for j in range(i + 1, len(q)))
        mm["diam"] = min(mm["diam"], d); mm["face"] = min(mm["face"], e[e > 1e-14].min())
        mm["ap"] = min(mm["ap"], 4 * A / e.sum()); mm["ad"] = min(mm["ad"], A / d); mm["w"] = min(mm["w"], width(q))
    for k in (3, 1):
        f = f"runs/k{k}_{b}/c.dat"
        if os.path.exists(f):
            dt = float(open(f).read().split("\n")[1].split()[2])
            out.append(dict(mesh=b, k=k, dt=dt, **mm))
json.dump(out, open("cfl.json", "w"), indent=1)
print(f"{'mesh':18} k {'dt':>9} | dt/face  dt/diam  dt/(4A/P)  dt/(A/d)  dt/width")
for r in out:
    print(f"{r['mesh']:18} {r['k']} {r['dt']:9.3e} | {r['dt']/r['face']:7.3f} {r['dt']/r['diam']:8.4f} {r['dt']/r['ap']:9.3f} {r['dt']/r['ad']:9.3f} {r['dt']/r['w']:8.3f}")
