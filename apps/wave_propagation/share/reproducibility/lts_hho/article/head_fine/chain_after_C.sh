#!/bin/bash
# 2026-11-16: after chains C and D, alone on the machine, 45G. Order chosen so that the paper can be updated during
# the day: timed quadtree runs at night (Tq, nothing else may run), Voronoi stable step (V1) and full run (V2), untimed
# (compiling and plotting allowed), then the timed Voronoi runs (Tv, nothing else may run)
cd "$(dirname "$0")"
./chain_T.sh q || exit 1
echo "$(date +%T) TIMED QUADTREE RUNS DONE"
./chain_V1.sh || exit 1
./chain_V2.sh
echo "$(date +%T) VORONOI FULL RUN DONE, TIMED VORONOI RUNS START"
./chain_T.sh v || exit 1
echo "$(date +%T) chain_after_C done"
