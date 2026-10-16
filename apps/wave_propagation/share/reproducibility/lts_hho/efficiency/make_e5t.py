#!/usr/bin/env python3
"""E5 timing: ERK4, MLTS all ratios 2 (direct, reuse), the best composition of the work model, two levels (LTS, ratio
2^D) and LI, overlap 0 (graded 2:1 meshes: stable), 1 and 6 threads; coarse step 0.8 x the limit of the background."""
import json, os, math
M_ = os.path.abspath("meshes")
best = {"pt_W1_D8": "{2, 4, 4, 8}", "pt_W2_D8": "{4}", "pt_W4_D8": "{4}", "pt_W1_D6": "{2, 4, 8}",
        "line_W1_D5": "{4, 8}", "disk_W1_D4": "{2, 8}"}
W = json.load(open("e5w.json"))
jobs = ["# E5 timing"]
for mesh, br in best.items():
    D = int(mesh.split("_D")[1]); P = 2 ** D
    N = W[f"{mesh}_o0"]["N"]
    Mst = max(10, min(4000, round(2e9 / (4 * N * P))))
    h = math.sqrt(2) / 32
    tf = Mst * 0.8 * 0.0826 * h
    mf = f"mesh.files = {{ '{M_}/{mesh}.txt' }}"
    c = f"time.final_time = {tf:.10g}"
    for th in (1, 6):
        b = f"e5t/{mesh}_t{th}"
        jobs += [f"{b}_erk|{th}|run_poly.lua|{mf}|time.scheme = 'erk'; time.n_steps = {Mst * P}; {c}",
                 f"{b}_all2|{th}|run_poly.lua|{mf}|time.scheme = 'mlts'; time.overlap = 0; time.n_steps = {Mst}; {c}",
                 f"{b}_all2r|{th}|run_poly.lua|{mf}|time.scheme = 'mlts'; time.overlap = 0; time.mlts_reuse = true; time.n_steps = {Mst}; {c}",
                 f"{b}_bestr|{th}|run_poly.lua|{mf}|time.scheme = 'mlts'; time.overlap = 0; time.mlts_reuse = true; time.level_ratios = {br}; time.n_steps = {Mst}; {c}",
                 f"{b}_lts|{th}|run_poly.lua|{mf}|time.scheme = 'lts'; time.fine_layers = 0; time.n_steps = {Mst}; {c}",
                 f"{b}_li|{th}|run_poly.lua|{mf}|time.scheme = 'li'; time.fine_layers = 0; time.li_system = 'constraint'; time.li_coupling = 'stage'; time.li_condensation = true; time.n_steps = {Mst}; {c}"]
open("jobs_e5t.txt", "w").write("\n".join(jobs) + "\n")
print(len(jobs) - 1, "jobs")
