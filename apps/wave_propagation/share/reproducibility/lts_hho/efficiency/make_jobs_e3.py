#!/usr/bin/env python3
"""E3: time error as a function of the overlap (16 x 16 mesh, central box refined 3 times, smooth solution), against
the same scheme at a 4 times smaller step; MLTS direct and reuse (same scheme: equal up to round-off)."""
lines = ["# E3 time error vs overlap"]
base = "cart.lua|mesh.level = 1; mesh.refine_levels = 3"
study = "study = { type = 'time', time_levels = 4, reference = 'numerical', reference_factor = 4 }"
for o in range(3):
    lines.append(f"e3/lts_o{o}|6|{base}|time.scheme = 'lts'; time.fine_layers = {o}; time.n_steps = 160|{study}")
    lines.append(f"e3/mlts_o{o}|6|{base}|time.scheme = 'mlts'; time.overlap = {o}; time.n_steps = 160|{study}")
lines.append(f"e3/mltsr_o1|6|{base}|time.scheme = 'mlts'; time.mlts_reuse = true; time.n_steps = 160|{study}")
lines.append(f"e3/erk|6|{base}|time.scheme = 'erk'; time.n_steps = 1280|{study}")
open("jobs_e3.txt", "w").write("\n".join(lines) + "\n")
print(len(lines) - 1, "jobs")
