#!/bin/bash
# 2026-11-17, plan B (if the plate mesh graded towards the titanium, fit2g, cuts the energy loss of the plate): plate
# runs redone on head_s80_fit2g_plate.txt; the intact and water runs (fit2) and the staircase runs (tfus2) are unchanged.
# After chains C and G, alone on the machine, 45G. Order: ERK(4) stable step with the short LI runs (untimed), timed
# quadtree runs at night (Tq, nothing else may run), LI at dt/2 with the probes and the full MLTS run (untimed), Voronoi
# stable step (V1) and full run (V2, untimed), timed Voronoi runs (Tv, nothing else may run).
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
./run_set.sh jobs_g_1a.txt out ./wave_propagation_hf $I &
./run_set.sh jobs_g_1b.txt out ./wave_propagation_hf $I &
wait
./chain_T.sh q fit2g || exit 1
echo "$(date +%T) TIMED QUADTREE RUNS DONE"
./run_set.sh jobs_g_2a.txt out ./wave_propagation_hf $I &
./run_set.sh jobs_g_2b.txt out ./wave_propagation_hf $I &
wait
echo "$(date +%T) UNTIMED PLATE RUNS DONE"
./chain_V1.sh || exit 1
./chain_V2.sh
echo "$(date +%T) VORONOI FULL RUN DONE, TIMED VORONOI RUNS START"
./chain_T.sh v || exit 1
echo "$(date +%T) chain_after_G done"
