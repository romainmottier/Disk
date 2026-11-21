#!/bin/bash
# 2026-11-17: Voronoi part of chain_after_G.sh (V1: stable step, V2: full MLTS run and memory probes, Tv: timed runs),
# relaunched alone at MemoryMax 60G after the OOM of vor2/cfl_mlts at 45G. The fine Voronoi mesh shares no block, and
# its setup peak is ~53 GB (coarse Voronoi probe: 13.6 GB max RSS for 4.9 GB of algebraic memory, x3.9 cells). The user
# allowed this exception to the 50 GB rule for these runs only, one job at a time. vor2/mem_li is expected to exceed 60G
# (LI factors >= 80 GB): it runs last of V2, and its log may still give nnz(L + U) before the OOM.
cd "$(dirname "$0")"
echo "$(date +%T) chain_after_V start"
./chain_V1.sh || exit 1
./chain_V2.sh
echo "$(date +%T) VORONOI FULL RUN DONE, TIMED VORONOI RUNS START"
./chain_T.sh v || exit 1
echo "$(date +%T) chain_after_V done"
