#!/usr/bin/env python3
"""Per-call overhead of the products: T_step = kappa (W + n0 C), C = products (kernel calls) per macro step; kappa from
ERK4 on the same mesh and thread count; n0 fitted (least squares on log T) per scheme family, 1 and 6 threads."""
import glob, json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import read_run

def calls_mlts(r, reuse):
    c = 0
    for lv, band, w, mult in r["levels_work"]:
        c += (4 if lv == "terminal" else (8 if reuse else 7)) * float(mult)
    return c

recs = []
for exp in ("e4", "e5t"):
    groups = {}
    for d in glob.glob(f"runs/{exp}/*"):
        m = re.match(r"(.+)_t(\d+)_(\w+)$", os.path.basename(d))
        if m: groups.setdefault((m.group(1), int(m.group(2))), {})[m.group(3)] = d
    for (mesh, th), g in groups.items():
        if "erk" not in g: continue
        e = read_run(g["erk"])
        for s, d in g.items():
            if s == "erk": continue
            r = read_run(d)
            if r is None or "W" not in r: continue
            M = r["steps"]
            if s == "lts":
                p = round(r["h_max"] / r["h_min"]); t = (r["t_coarse"] + r["t_fine"]) / M; C = 8 + 4 * p
            elif s.startswith("mlts") or s.startswith("all2") or s.startswith("best"):
                t = (sum(r["t_coef"]) + r["t_term"] + r["t_close"]) / M; C = calls_mlts(r, s.endswith("r"))
            else:
                continue
            P = 2 ** (len(r["cells_per_level"]) - 1) if "cells_per_level" in r and s != "lts" else None
            ne = e["steps"]
            kappa = e["tint"] / (ne * 4 * r["N"])
            recs.append(dict(exp=exp, mesh=mesh, th=th, s=s, t=t, W=r["W"], C=C, kappa=kappa, N=r["N"],
                             erk_step=e["tint"] / ne, ratio=ne / M))
for th in (1, 6):
    for fam in ("lts", "mlts"):
        R = [x for x in recs if x["th"] == th and ((x["s"] == "lts") == (fam == "lts"))]
        y = np.array([x["t"] / x["kappa"] for x in R]); W = np.array([x["W"] for x in R]); C = np.array([x["C"] for x in R])
        # fit y = a W + n0 C (a: efficiency of the products, n0: per call), least squares relative
        A = np.column_stack([W / y, C / y]); sol, *_ = np.linalg.lstsq(A, np.ones_like(y), rcond=None)
        pred = sol[0] * W + sol[1] * C
        err = np.abs(pred / y - 1)
        print(f"{th} thread(s) {fam:5}: T = kappa ({sol[0]:.2f} W + {sol[1]:.0f} C), {len(R)} runs, error median {np.median(err):.0%}, max {err.max():.0%}")
json.dump(recs, open("overhead.json", "w"), indent=1)
