#!/bin/bash
# Voronoi (fine), part 2 (2026-11-16, alone): full MLTS run and memory probe of jobs_vor_2.txt (made by chain_V1.sh)
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
[ -f jobs_vor_2.txt ] || { echo "$(date +%T) missing jobs_vor_2.txt (chain_V1.sh)"; exit 1; }
./run_set.sh jobs_vor_2.txt out ./wave_propagation_hf $I
