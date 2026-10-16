#!/usr/bin/env python3
"""E4: cost of the two-level schemes against ERK4 on two-size meshes (fine fraction theta, ratio p, geometry).
Writes make_meshes_e4.sh and jobs_e4.txt (dir|threads|input|lua lines). Time steps: dt_L = 0.8 x the coarse limit
(nu = c dt / h_max = 0.0826 on these meshes, h_max the diameter of a background square), dt_E = dt_L / p; M macro steps
so that ERK4 does about 2e9 unknown-products."""
import os, math
G = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/generate_refined_quadtree.py"
M = os.path.abspath("meshes")
# (p, n, m): central block of m x m background squares split to depth log2 p
cases = {2: [(64, 2), (64, 4), (64, 8), (64, 16), (64, 24), (64, 32), (64, 48)],
         4: [(64, 1), (64, 2), (64, 4), (64, 8), (64, 12), (64, 16), (64, 24)],
         8: [(128, 1), (64, 1), (64, 2), (64, 4), (64, 6), (64, 8), (64, 12)],
         16: [(256, 1), (128, 1), (64, 1), (64, 2), (64, 3), (64, 4), (64, 6)]}
meshes, jobs = [], ["# E4 two-level cost: dir|threads|input|lua lines"]

def unknowns(nc_coarse, nc_fine):
    # 24 cell unknowns (k = 3, acoustic), 4 per face, about 2 faces per cell
    return 32 * (nc_coarse + nc_fine)

def add(name, gen_args, p, nc_coarse, nc_fine, n):
    path = f"{M}/{name}.txt"
    meshes.append(f"python3 {G} {gen_args} --output {path}")
    N = unknowns(nc_coarse, nc_fine)
    Mst = max(10, min(4000, round(2e9 / (4 * p * N))))
    h = math.sqrt(2) / n
    tf = Mst * 0.8 * 0.0826 * h
    mesh = f"mesh.files = {{ '{path}' }}"
    common = f"time.final_time = {tf:.10g}"
    for th in (1, 6):
        jobs.append(f"e4/{name}_t{th}_erk|{th}|run_poly.lua|{mesh}|time.scheme = 'erk'; time.n_steps = {Mst * p}; {common}")
        jobs.append(f"e4/{name}_t{th}_lts|{th}|run_poly.lua|{mesh}|time.scheme = 'lts'; time.n_steps = {Mst}; {common}")
        jobs.append(f"e4/{name}_t{th}_mlts|{th}|run_poly.lua|{mesh}|time.scheme = 'mlts'; time.level_ratios = {{{p}}}; "
                    f"time.n_steps = {Mst}; {common}")
        jobs.append(f"e4/{name}_t{th}_mltsr|{th}|run_poly.lua|{mesh}|time.scheme = 'mlts'; time.level_ratios = {{{p}}}; "
                    f"time.mlts_reuse = true; time.n_steps = {Mst}; {common}")
        jobs.append(f"e4/{name}_t{th}_li|{th}|run_poly.lua|{mesh}|time.scheme = 'li'; time.li_system = 'constraint'; "
                    f"time.li_coupling = 'stage'; time.li_condensation = true; time.n_steps = {Mst}; {common}")

for p, lst in cases.items():
    L = int(round(math.log2(p)))
    for (n, m) in lst:
        a, b = (n - m) // 2, (n + m) // 2
        box = f"box:{a / n!r},{b / n!r},{a / n!r},{b / n!r}"
        add(f"two_n{n}_m{m}_L{L}", f"--n {n} --depth {L} --grading two --feature {box}", p, n * n - m * m, m * m * 4 ** L, n)
# geometry: a strip of one square across the domain, 16 scattered single squares (p = 4 and 16), n = 64
for p in (4, 16):
    L = int(round(math.log2(p)))
    n = 64
    y = 31.5 / 64
    add(f"strip_n64_L{L}", f"--n 64 --depth {L} --grading two --feature line:0.0,{y!r},1.0,{y!r}", p, n * n - n, n * 4 ** L, n)
    pts = " ".join(f"--feature point:{(8 * i + 4.5) / 64!r},{(8 * j + 4.5) / 64!r}" for i in range(0, 8, 2) for j in range(0, 8, 2))
    add(f"dots_n64_L{L}", f"--n 64 --depth {L} --grading two {pts}", p, n * n - 16, 16 * 4 ** L, n)
open("make_meshes_e4.sh", "w").write("#!/bin/bash\n" + "\n".join(meshes) + "\n")
open("jobs_e4.txt", "w").write("\n".join(jobs) + "\n")
print(len(meshes), "meshes,", len(jobs) - 1, "jobs")
