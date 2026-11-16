#!/bin/bash
# Campaign of 2026-11-15 of the LTS_HHO preprint (main.tex): the timings of Sections 5.2.2 (L-shaped domain), 5.3
# (efficiency) and 5.4 (head model) under the measurement rules of the preprint:
#   - Sections 5.2.2 and 5.3 on ONE thread (work of the schemes), the head model on six threads (practice) and one;
#   - LI solves its fine system with MKL Pardiso in both settings, with the same options (time.li_solver = 'pardiso'
#     in the setting (1E), time.li_refinement = 0 in both: nested dissection, parallel substitutions, no iterative
#     refinement);
#   - every timed run alone on the machine (NOTHING else running, no numpy/BLAS, no LaTeX build), three times
#     (<out>/out, <out>/out_r2, <out>/out_r3), best of three; then the untimed full runs of the head model
#     (<out>/out_full: signals, snapshots, energies) and the algebraic memory (<out>/out_mem, 3 steps each).
# Resumable: a job is done when its time.txt exists (written only on success). Job line:
#   dir|threads|input (../inputs)|Lua statement|Lua statement|...
# Needs a build with Intel MKL (-DHAVE_INTEL_MKL). About 6 h 30 on the machine of the preprint.
#   usage: run_pardiso.sh <out> <wave_propagation executable>
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_pardiso.sh <out> <executable>}
exe=$(readlink -f "${2:?usage: run_pardiso.sh <out> <executable>}")
inputs="$here/../inputs"
export LC_ALL=C

run_set() {   # run_set <jobs file> <out>
  mkdir -p "$2"
  while IFS= read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    IFS='|' read -r -a f <<< "$line"
    d="$2/${f[0]}"; th=${f[1]}; in="$inputs/${f[2]}"
    [ -f "$d/time.txt" ] && continue
    rm -rf "$d"; mkdir -p "$d"; cp "$in" "$d/case.lua"
    for l in "${f[@]:3}"; do echo "$l" >> "$d/case.lua"; done
    echo "$(date +%T) ${f[0]}"
    ( cd "$d" && /usr/bin/time -f "%e s %M KB" -o time.txt.part env OMP_NUM_THREADS=$th MKL_NUM_THREADS=$th "$exe" case.lua \
        > log.txt 2>&1 && mv time.txt.part time.txt ) || echo "$(date +%T) FAILED ${f[0]}"
  done < "$1"
}

for o in out out_r2 out_r3; do
  for j in jobs_cost.txt jobs_ls1.txt jobs_ls2.txt jobs_head.txt; do
    run_set "$here/$j" "$out/$o"
  done
done
run_set "$here/jobs_full.txt" "$out/out_full"
run_set "$here/jobs_mem.txt" "$out/out_mem"
echo "$(date +%T) finished"
