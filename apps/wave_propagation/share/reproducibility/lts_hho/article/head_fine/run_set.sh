#!/bin/bash
# Pardiso campaign (2026-11-15): runs the jobs of a jobs file one at a time, resumable (a job is done when its time.txt
# exists; time.txt written only on success), records the load of the machine before and after each job (load.txt)
#   usage: run_set.sh <jobs file> <out> <executable> <inputs dir>
export LC_ALL=C
jobs=${1:?}; out=${2:?}; exe=$(readlink -f "${3:?}"); inputs=${4:?}
mkdir -p "$out"
while IFS= read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  IFS='|' read -r -a f <<< "$line"
  d="$out/${f[0]}"; th=${f[1]}; in="$inputs/${f[2]}"
  [ -f "$d/time.txt" ] && continue
  [ -f "$in" ] || { echo "$(date +%T) ERROR missing input $in"; continue; }
  rm -rf "$d"; mkdir -p "$d"; cp "$in" "$d/case.lua"
  for l in "${f[@]:3}"; do echo "$l" >> "$d/case.lua"; done
  echo "$(date +%T) ${f[0]}"
  { echo "before $(date +%s) $(cat /proc/loadavg)"; top -b -n1 -o %CPU -w 160 | sed -n 7,12p; } > "$d/load.txt"
  ( cd "$d" && /usr/bin/time -f "%e s %M KB" -o time.txt.part env OMP_NUM_THREADS=$th MKL_NUM_THREADS=$th "$exe" case.lua > log.txt 2>&1 \
    && mv time.txt.part time.txt ) || echo "$(date +%T) FAILED ${f[0]}"
  { echo "after $(date +%s) $(cat /proc/loadavg)"; top -b -n1 -o %CPU -w 160 | sed -n 7,12p; } >> "$d/load.txt"
done < "$jobs"
echo "$(date +%T) done $jobs"
