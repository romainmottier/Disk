#!/bin/bash
# step 2 of the fine-mesh campaign: chains A and B in parallel (untimed runs, six threads each)
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
./run_set.sh jobs_full_A.txt out ./wave_propagation_hf $I > chain_A.log 2>&1 &
./chain_B.sh > chain_B.log 2>&1 &
wait
echo "$(date +%T) step 2 done"
