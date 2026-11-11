#!/bin/bash
# main.tex, what sets the CFL limit of a cell, setting (2E): the meshes of run_cfl_cells.sh and run_cfl_random.sh
# (their meshes/ directories), largest stable step of ERK4 for k = 2 (mixed order, O(1/h) stabilization, weight 3 gamma*
# computed on each mesh: inputs/cfl_cells2.lua), then cfl.json (analyze_cfl.py, with gamma*). Resumable; no timing.
#   usage: run_cfl_2e.sh <out> <executable> <meshes directory>...
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_cfl_2e.sh <out> <executable> <meshes directory>...}; exe=${2:?}; shift 2
mkdir -p "$out/meshes" "$out/runs"; out=$(cd "$out" && pwd)
for d in "$@"; do for m in "$d"/*.txt; do ln -sf "$(realpath "$m")" "$out/meshes/$(basename "$m")"; done; done
for m in "$out"/meshes/*.txt; do
  b=$(basename "$m" .txt)
  d="$out/runs/k2_$b"
  [ -f "$d/c.dat" ] && continue
  rm -rf "$d"; mkdir -p "$d"; cp "$here/inputs/cfl_cells2.lua" "$d/case.lua"
  echo "mesh.files = { '$m' }" >> "$d/case.lua"
  ( cd "$d" && env OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 "$exe" case.lua > log.txt 2>&1 )
  echo "$b: $(sed -n 2p "$d/c.dat" | awk '{print $3}')"
done
python3 "$here/analyze_cfl.py" "$out"
