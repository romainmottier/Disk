#!/usr/bin/env python3
"""E2: stability of LTS / MLTS as a function of the overlap (cfl studies), two-size and graded meshes."""
import os
M = os.path.abspath("meshes")
lines = ["# E2 stability vs overlap: dir|threads|input|lua lines"]
for L in range(1, 6):
    p = 2 ** L
    mesh = f"mesh.files = {{ '{M}/two_n16_m2_L{L}.txt' }}"
    lines.append(f"e2/two_L{L}_erk|6|cfl_poly.lua|{mesh}|time.scheme = 'erk'; time.n_steps = {160 * p}")
    for o in range(3):
        lines.append(f"e2/two_L{L}_lts_o{o}|6|cfl_poly.lua|{mesh}|time.scheme = 'lts'; time.fine_layers = {o}")
        lines.append(f"e2/two_L{L}_mlts_o{o}|6|cfl_poly.lua|{mesh}|time.scheme = 'mlts'; time.level_ratios = {{{p}}}; time.overlap = {o}")
    for o in range(3):
        lines.append(f"e2/grad_L{L}_lts_o{o}|6|cfl_cart.lua|mesh.refine_levels = {L}|time.scheme = 'lts'; time.fine_layers = {o}")
        lines.append(f"e2/grad_L{L}_mlts_o{o}|6|cfl_cart.lua|mesh.refine_levels = {L}|time.scheme = 'mlts'; time.overlap = {o}")
    if L in (2, 4):
        for o in range(2):
            lines.append(f"e2/grad_L{L}_mlts4_o{o}|6|cfl_cart.lua|mesh.refine_levels = {L}|time.scheme = 'mlts'; "
                         f"time.level_ratios = {{4}}; time.overlap = {o}")
# long horizon (t_f = 10, about 1400 coarse steps): weak instabilities
for L in range(1, 4):
    p = 2 ** L
    mesh = f"mesh.files = {{ '{M}/two_n16_m2_L{L}.txt' }}"
    for o in range(2):
        lines.append(f"e2/two_L{L}_lts_o{o}_T10|6|cfl_poly.lua|{mesh}|time.scheme = 'lts'; time.fine_layers = {o}; "
                     f"time.final_time = 10.0; time.n_steps = 1600")
        lines.append(f"e2/two_L{L}_mlts_o{o}_T10|6|cfl_poly.lua|{mesh}|time.scheme = 'mlts'; time.level_ratios = {{{p}}}; "
                     f"time.overlap = {o}; time.final_time = 10.0; time.n_steps = 1600")
        lines.append(f"e2/grad_L{L}_lts_o{o}_T10|6|cfl_cart.lua|mesh.refine_levels = {L}|time.scheme = 'lts'; "
                     f"time.fine_layers = {o}; time.final_time = 10.0; time.n_steps = 1600")
        lines.append(f"e2/grad_L{L}_mlts_o{o}_T10|6|cfl_cart.lua|mesh.refine_levels = {L}|time.scheme = 'mlts'; "
                     f"time.overlap = {o}; time.final_time = 10.0; time.n_steps = 1600")
open("jobs_e2.txt", "w").write("\n".join(lines) + "\n")
print(len(lines) - 1, "jobs")
