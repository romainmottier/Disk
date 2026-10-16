#!/usr/bin/env python3
"""E5: multilevel. Graded quadtrees (band of W squares per level, 2:1 balanced) refined around a point, a segment or a
disk down to the depth D; for each mesh, the MLTS work per macro step (output.profile, one macro step) of every
composition of D into level ratios 2^a (2^(D-1) choices: all the ways of merging consecutive levels), direct and
reuse. Usage: make_e5.py OVERLAP  -> make_meshes_e5.sh, jobs_e5w.txt"""
import os, sys, itertools
OV = int(sys.argv[1]) if len(sys.argv) > 1 else 1
G = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/generate_refined_quadtree.py"
M = os.path.abspath("meshes")
specs = []
for W in (1, 2, 4):
    for D in (4, 6, 8):
        specs.append((f"pt_W{W}_D{D}", f"--n 32 --depth {D} --grading rings --rings {W} --balance --feature point:0.5,0.5", D))
for W in (1, 2):
    for D in (3, 5):
        specs.append((f"line_W{W}_D{D}", f"--n 32 --depth {D} --grading rings --rings {W} --balance "
                                         f"--feature line:0.125,0.5,0.875,0.5", D))
        specs.append((f"disk_W{W}_D{D - 1}", f"--n 32 --depth {D - 1} --grading rings --rings {W} --balance "
                                             f"--feature disk:0.5,0.5,0.125", D - 1))

def compositions(D):
    for cuts in itertools.product((0, 1), repeat=D - 1):
        parts, a = [], 1
        for c in cuts:
            if c: parts.append(a); a = 1
            else: a += 1
        parts.append(a)
        yield parts

meshes, jobs = [], [f"# E5 work of the level compositions (overlap {OV})"]
for name, args, D in specs:
    path = f"{M}/{name}.txt"
    meshes.append(f"python3 {G} {args} --output {path}")
    mesh = f"mesh.files = {{ '{path}' }}"
    for parts in compositions(D):
        r = ", ".join(str(2 ** a) for a in parts)
        tag = "-".join(str(a) for a in parts)
        for reuse in (False, True):
            jobs.append(f"e5w/{name}_o{OV}_{tag}{'_r' if reuse else ''}|6|run_poly.lua|{mesh}|time.scheme = 'mlts'; "
                        f"time.level_ratios = {{{r}}}; time.overlap = {OV}; time.mlts_reuse = {str(reuse).lower()}; "
                        f"time.n_steps = 1; time.final_time = 1e-4")
open("make_meshes_e5.sh", "w").write("#!/bin/bash\n" + "\n".join(meshes) + "\n")
open(f"jobs_e5w_o{OV}.txt", "w").write("\n".join(jobs) + "\n")
print(len(meshes), "meshes,", len(jobs) - 1, "jobs")
