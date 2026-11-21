#!/bin/bash
# 2026-11-16: head meshes twice finer than those of the preprint (all cell sizes halved), same generators and options
# (LISEZMOI: --margin 40 --base 4 --size-soft 1 --size-bone 0.5 --interface 0.5 --size-implant 0.125 [--delta 0.3])
G=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
PY=/home/mottie0000/Github/env_py/bin/python
S="--margin 40 --base 4 --size-soft 0.5 --size-bone 0.25 --interface 0.25 --size-implant 0.0625"
cd $G
run() { local t0=$(date +%s); echo "$(date +%T) start $1"; "$PY" "${@:2}" > "$M/$1.log" 2>&1 && echo "$(date +%T) done $1 in $(( $(date +%s) - t0 )) s" || echo "$(date +%T) FAILED $1"; }
run fit2_plate  generate_head_model_fitted.py  --output $M/head_s80_fit2_plate.txt  $S --delta 0.3 --plate 50 66
"$PY" fix_degenerate_cells.py $M/head_s80_fit2_plate.txt $M/head_s80_fit2_plate_fixed.txt && mv $M/head_s80_fit2_plate_fixed.txt $M/head_s80_fit2_plate.txt   # 2026-11-16: one flat cell (NaN)
run fit2_intact generate_head_model_fitted.py  --output $M/head_s80_fit2_intact.txt $S --delta 0.3
run hf2_plate   generate_head_model.py         --output $M/head_s80_hf2_plate.txt   $S --plate 50 66 --fill-gap   # 2026-11-16: plate on bone (was a 5.2 mm^2 scalp gap)
run hf2_intact  generate_head_model.py         --output $M/head_s80_hf2_intact.txt  $S
run vor2_plate  generate_head_model_voronoi.py --output $M/head_s80_vor2_plate.txt  $S --plate 50 66
# 2026-11-17 (plan B, Section 5.4): plate mesh graded towards the titanium, used by every run with the plate (made by
# chain_G.sh): gen_graded.py replaces build() of generate_head_model.py at run time, side <= max(0.0625, 0.25 x distance)
run fit2g_plate gen_graded.py 0.25 --output $M/head_s80_fit2g_plate.txt $S --delta 0.3 --plate 50 66
"$PY" fix_degenerate_cells.py $M/head_s80_fit2g_plate.txt $M/head_s80_fit2g_plate_fixed.txt && mv $M/head_s80_fit2g_plate_fixed.txt $M/head_s80_fit2g_plate.txt
echo "$(date +%T) all done"
