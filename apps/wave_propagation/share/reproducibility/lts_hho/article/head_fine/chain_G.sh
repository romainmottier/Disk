#!/bin/bash
# 2026-11-16 diagnostic G: plate mesh graded towards the titanium (gen_graded.py: side <= max(0.0625, d/4)), then LI as
# fit2/plate_li (jobs_diag_G.txt -> out/fit2g/plate_li). Compare its energy loss with fit2 (0.440 %) and fit2i (0.384 %).
cd "$(dirname "$0")"
PY=/home/mottie0000/Github/env_py/bin/python
B=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
S="--margin 40 --base 4 --size-soft 0.5 --size-bone 0.25 --interface 0.25 --size-implant 0.0625"
echo "$(date +%T) mesh fit2g_plate"
"$PY" gen_graded.py 0.25 --output $M/head_s80_fit2g_plate.txt $S --delta 0.3 --plate 50 66 > $M/fit2g_plate.log 2>&1 \
  || { echo "$(date +%T) mesh FAILED"; exit 1; }
"$PY" $B/fix_degenerate_cells.py $M/head_s80_fit2g_plate.txt $M/head_s80_fit2g_plate_fixed.txt >> $M/fit2g_plate.log 2>&1 \
  && mv $M/head_s80_fit2g_plate_fixed.txt $M/head_s80_fit2g_plate.txt || { echo "$(date +%T) fix FAILED"; exit 1; }
echo "$(date +%T) mesh done: $(grep -E 'graded build|quadtree:|points,' $M/fit2g_plate.log | tr '\n' ' ')"
./run_set.sh jobs_diag_G.txt out ./wave_propagation_hf $I
