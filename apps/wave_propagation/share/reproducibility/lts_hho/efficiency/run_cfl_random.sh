#!/bin/bash
# main.tex, what sets the CFL limit of a cell, random cut cells (complement of run_cfl_cells.sh): 16 x 16 squares with
#   cut_<i>     the central square cut by a random line (40 meshes, seed 1: angle uniform, distance of the line to the
#               farthest corner log-uniform between h/150 and the half-diagonal, so that the smaller piece ranges from a
#               sliver to a half square),
#   circle_<i>  every square crossed by a random circle cut along its chord (20 meshes, seed 2: centre in
#               [0.4, 0.6]^2, radius in [0.2, 0.3], as an interface),
# largest stable step of ERK4 for k = 3 (and k = 1 on the circle meshes, 2026-11-06), then cfl.json (analyze_cfl.py).
# Resumable; no timing involved.
#   usage: run_cfl_random.sh <out> [wave_propagation executable]
here=$(cd "$(dirname "$0")" && pwd)
gen=$here/../../../../src/mesh_generation/lts_efficiency/generate_cfl_cells.py
out=${1:?usage: run_cfl_random.sh <out> [executable]}; exe=${2:-wave_propagation}
mkdir -p "$out/meshes" "$out/runs"; out=$(cd "$out" && pwd)
python3 - "$gen" "$out/meshes" <<'EOF'
import subprocess, sys
import numpy as np
gen, dst = sys.argv[1], sys.argv[2]
r = np.random.default_rng(1)
for i in range(40):
    th = r.uniform(0, 180)
    dmax = (abs(np.cos(np.radians(th))) + abs(np.sin(np.radians(th)))) / 2
    off = (dmax - np.exp(r.uniform(np.log(1 / 150), np.log(dmax)))) * r.choice((-1, 1))
    subprocess.run(["python3", gen, "--kind", "cut", "--theta", f"{th:.6f}", "--offset", f"{off:.6f}",
                    "--output", f"{dst}/cut_{i:02d}.txt"], check=True, stdout=subprocess.DEVNULL)
r = np.random.default_rng(2)
i = 0
while i < 20:
    cx, cy, rad = r.uniform(0.4, 0.6), r.uniform(0.4, 0.6), r.uniform(0.2, 0.3)
    p = subprocess.run(["python3", gen, "--kind", "circle", "--cx", f"{cx:.6f}", "--cy", f"{cy:.6f}", "--radius",
                        f"{rad:.6f}", "--output", f"{dst}/circle_{i:02d}.txt"], stdout=subprocess.DEVNULL)
    i += p.returncode == 0
EOF
for m in "$out"/meshes/*.txt; do
  b=$(basename "$m" .txt)
  case $b in circle_*) ks="3 1";; *) ks=3;; esac
  for k in $ks; do
    d="$out/runs/k${k}_$b"
    [ -f "$d/c.dat" ] && continue
    rm -rf "$d"; mkdir -p "$d"; cp "$here/inputs/cfl_cells.lua" "$d/case.lua"
    echo "mesh.files = { '$m' }" >> "$d/case.lua"; echo "space.k = $k" >> "$d/case.lua"
    ( cd "$d" && env OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 "$exe" case.lua > log.txt 2>&1 )
    echo "$b k = $k: $(sed -n 2p "$d/c.dat" | awk '{print $3}')"
  done
done
python3 "$here/analyze_cfl.py" "$out"
