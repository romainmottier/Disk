#!/usr/bin/env python3
"""E5 timing: measured gains over ERK4 on the graded quadtrees (runs/e5t), times of the time loop (the construction of
the levels and the factorization of LI excluded), and ideal gain G_id = 2^D sum N_l / sum 2^l N_l from the cells per
level of MLTS with all the ratios two, and peak memories. Writes e5t.json for the figures (multilevel_measured)."""
import glob, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import read_run


def loop_time(s, r):
    if s == "erk":
        return r["tint"]
    if s == "lts":
        return r["t_coarse"] + r["t_fine"]
    if s == "li":
        return r["t_coarse"] + r["t_solve"] + r["t_stage"]
    return sum(r["t_coef"]) + r["t_term"] + r["t_close"]


groups = {}
for d in sorted(glob.glob("runs/e5t/*")):
    m = re.match(r"(.+)_t(\d+)_(\w+)$", os.path.basename(d))
    if m and os.path.exists(f"{d}/time.txt"):
        groups.setdefault((m.group(1), int(m.group(2))), {})[m.group(3)] = read_run(d)
out = []
print(f"{'mesh':12} {'th':>2} | {'all2':>6} {'all2r':>6} {'bestr':>6} {'lts':>6} {'li':>6} | {'G_id':>6}")
for (mesh, th), g in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if "erk" not in g:
        continue
    te = loop_time("erk", g["erk"])
    G = {s: te / loop_time(s, r) for s, r in g.items() if s != "erk" and r is not None}
    c = g["all2"]["cells_per_level"]
    D = len(c) - 1
    gid = 2 ** D * sum(c) / sum(2 ** l * n for l, n in enumerate(c))
    mem = {s: r["mem"] for s, r in g.items() if r is not None and r.get("mem")}        # peak memory (MiB)
    out.append(dict(mesh=mesh, threads=th, G=G, Gid=gid, cells_per_level=c, mem=mem))
    f = lambda k: f"{G[k]:6.2f}" if k in G else f"{'-':>6}"
    print(f"{mesh:12} {th:2d} | {f('all2')} {f('all2r')} {f('bestr')} {f('lts')} {f('li')} | {gid:6.2f}")
json.dump(out, open("e5t.json", "w"), indent=1)
