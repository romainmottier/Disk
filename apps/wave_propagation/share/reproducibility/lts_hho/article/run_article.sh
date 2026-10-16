#!/bin/bash
# Campaign of Sections 4 to 6 of the LTS_HHO preprint (main.tex): runs the jobs of jobs.txt one at a time
# (clean timings: nothing else may run on the machine meanwhile), resumable (a job is done when <out>/<dir>/time.txt
# exists). Job line: dir|threads|input (inputs/)|Lua statement|Lua statement|...
#   usage: run_article.sh <out> [wave_propagation executable]
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_article.sh <out> [executable]}
exe=${2:-wave_propagation}
mkdir -p "$out"
while IFS= read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  IFS='|' read -r -a f <<< "$line"
  d="$out/${f[0]}"; th=${f[1]}; in="$here/inputs/${f[2]}"
  [ -f "$d/time.txt" ] && continue
  rm -rf "$d"; mkdir -p "$d"; cp "$in" "$d/case.lua"
  for l in "${f[@]:3}"; do echo "$l" >> "$d/case.lua"; done
  echo "$(date +%T) ${f[0]}"
  ( cd "$d" && /usr/bin/time -f "%e s %M KB" -o time.txt env OMP_NUM_THREADS=$th MKL_NUM_THREADS=$th "$exe" case.lua > log.txt 2>&1 )
done < "$here/jobs.txt"
