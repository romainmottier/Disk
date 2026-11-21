#!/bin/bash
# chain B: jobs_full_B1.txt, then the Voronoi MLTS full run at 0.9 of the stable step found by vor2/cfl_mlts
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
./run_set.sh jobs_full_B1.txt out ./wave_propagation_hf $I
dt=$(awk '!/^#/ {print $3; exit}' out/vor2/cfl_mlts/c.dat 2>/dev/null)
if [ -z "$dt" ]; then echo "$(date +%T) no stable step for the Voronoi mesh"; exit 1; fi
n=$(awk -v dt="$dt" 'BEGIN { n = 80.0 / (0.9 * dt); printf "%d", (n == int(n)) ? n : int(n) + 1 }')
echo "$(date +%T) Voronoi MLTS: stable dt $dt, production n_steps $n"
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
cat > jobs_full_B2.txt <<EOT
# made by chain_B.sh: MLTS on the fine Voronoi mesh, 80 us at 0.9 of the stable step ($dt us): $n steps
vor2/plate_mlts|6|tfus_vor.lua|mesh.files = { "$M/head_s80_vor2_plate.txt" }|time.scheme = 'mlts'; time.level_size = 'width'; time.mlts_reuse = true; time.level_ratios = {2, 2, 4, 8}; time.n_steps = $n
EOT
./run_set.sh jobs_full_B2.txt out ./wave_propagation_hf $I
