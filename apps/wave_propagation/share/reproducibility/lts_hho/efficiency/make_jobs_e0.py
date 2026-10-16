#!/usr/bin/env python3
"""E0: LTS = MLTS with two levels (same partition, results equal up to round-off; reuse of the chain equal to the direct
MLTS bit for bit): two-size mesh, graded box, L-shaped domain with source and Dirichlet data, elasto-acoustic medium."""
import os
M = os.path.abspath("meshes")
two = f"mesh.files = {{ '{M}/two_n16_m2_L3.txt' }}"
J = ["# E0 LTS = MLTS with two levels"]
for s, extra in (("lts", "time.scheme = 'lts'"), ("mlts", "time.scheme = 'mlts'; time.level_ratios = {8}"),
                 ("mltsr", "time.scheme = 'mlts'; time.level_ratios = {8}; time.mlts_reuse = true")):
    J.append(f"e0/two_L3_{s}|6|run_poly.lua|{two}|{extra}; time.final_time = 0.5; time.n_steps = 100")
    J.append(f"e0/grad_L3_{s}|6|cart.lua|mesh.level = 1; mesh.refine_levels = 3|{extra}; time.n_steps = 160|output.profile = true")
for s, extra in (("lts", "time.scheme = 'lts'; time.substeps = 1024"), ("mlts", "time.scheme = 'mlts'; time.level_ratios = {1024}"),
                 ("mltsr", "time.scheme = 'mlts'; time.level_ratios = {1024}; time.mlts_reuse = true")):
    J.append(f"e0/ls_{s}|6|lshape.lua|mesh.level = 1|{extra}; time.n_steps = 40|output.profile = true")
J.append("e0/ea_lts|6|ea.lua|time.scheme = 'lts'")
J.append("e0/ea_mlts|6|ea.lua|time.scheme = 'mlts'; time.level_ratios = {2}")
open("jobs_e0.txt", "w").write("\n".join(J) + "\n")
print(len(J) - 1, "jobs")
