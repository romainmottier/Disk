#!/usr/bin/env python3
"""E2 table: largest stable step of LTS / MLTS against the overlap, relative to the uniform 16 x 16 mesh (7.2993e-3)."""
import os, re, glob
R = "runs/e2"
ref = 1.0 / 137
rows = {}
for d in sorted(glob.glob(f"{R}/*")):
    name = os.path.basename(d)
    m = re.match(r"(two|grad)_L(\d)_(erk|lts|mlts4|mlts)(?:_o(\d))?(_T10)?$", name)
    if not m:
        continue
    fam, L, sch, o, T = m.group(1), int(m.group(2)), m.group(3), m.group(4), m.group(5) or ""
    f = f"{d}/c.dat"
    if not os.path.exists(f):
        rows[(fam, sch, o, T, L)] = None
        continue
    v = open(f).read().split("\n")[1].split()
    dt = float(v[2])
    rows[(fam, sch, o, T, L)] = dt
keys = sorted({k[:4] for k in rows}, key=lambda k: (k[0], k[3], k[1], k[2] or ""))
print(f"{'family':6} {'scheme':6} {'ov':3} {'tf':4} " + " ".join(f"{'p=' + str(2 ** L):>8}" for L in range(1, 6)))
for k in keys:
    vals = []
    for L in range(1, 6):
        v = rows.get(k + (L,), "")
        vals.append(f"{'':>8}" if v == "" else f"{'fail':>8}" if v is None else f"{v / ref:8.3f}")
    print(f"{k[0]:6} {k[1]:6} {k[2] or '-':3} {('10' if k[3] else '1'):4} " + " ".join(vals))
