#!/bin/bash
# main.tex 6.3 (what sets the stable step of a cell): meshes of generate_cfl_cells.py (16 x 16 squares, one
# modified cell, or anisotropic cells), largest stable step of ERK4 for k = 3 (and k = 1 for thin, face, none), then
# cfl.json (analyze_cfl.py) for ../make_cfl_cells.py. Resumable; no timing involved.
#   usage: run_cfl_cells.sh <out> [wave_propagation executable]
here=$(cd "$(dirname "$0")" && pwd)
gen=$here/../../../../src/mesh_generation/lts_efficiency/generate_cfl_cells.py
out=${1:?usage: run_cfl_cells.sh <out> [executable]}; exe=${2:-wave_propagation}
mkdir -p "$out/meshes" "$out/runs"; out=$(cd "$out" && pwd)
python3 "$gen" --kind none --output "$out/meshes/none.txt" > /dev/null
for e in 0.5 0.25 0.125 0.0625 0.03125 0.015625; do
  for k in thin face corner; do python3 "$gen" --kind $k --eps $e --output "$out/meshes/${k}_$e.txt" > /dev/null; done
done
for e in 0.5 0.25 0.125; do python3 "$gen" --kind aniso --eps $e --output "$out/meshes/aniso_$e.txt" > /dev/null; done
for e in 0.25 0.0625; do python3 "$gen" --kind slab --eps $e --output "$out/meshes/slab_$e.txt" > /dev/null; done
for m in "$out"/meshes/*.txt; do
  b=$(basename "$m" .txt)
  for k in 3 1; do
    if [ $k = 1 ]; then case $b in thin_*|face_*|none) ;; *) continue;; esac; fi
    d="$out/runs/k${k}_$b"
    [ -f "$d/c.dat" ] && continue
    rm -rf "$d"; mkdir -p "$d"; cp "$here/inputs/cfl_cells.lua" "$d/case.lua"
    echo "mesh.files = { '$m' }" >> "$d/case.lua"; echo "space.k = $k" >> "$d/case.lua"
    ( cd "$d" && env OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 "$exe" case.lua > log.txt 2>&1 )
    echo "$b k = $k: $(sed -n 2p "$d/c.dat" | awk '{print $3}')"
  done
done
python3 "$here/analyze_cfl.py" "$out"
