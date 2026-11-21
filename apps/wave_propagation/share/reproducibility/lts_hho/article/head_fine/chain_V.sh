#!/bin/bash
# Voronoi (fine) alone (its matrix-free operator has no shared block: > 27 GB): stable step of MLTS (random data, 1 us),
# then the MLTS full run at 0.9 of it, and the memory probe (3 steps)
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
V2="mesh.files = { \"$M/head_s80_vor2_plate.txt\" }"
MLTS="time.scheme = 'mlts'; time.level_size = 'width'; time.mlts_reuse = true; time.level_ratios = {2, 2, 4, 8}"
NOOUT="output.snapshots = {}; output.sensors = {}; output.energy = false; output.peak_samples = 0; output.profile = false"
rm -rf out/vor2/cfl_mlts
cat > jobs_vor_1.txt <<EOT
vor2/cfl_mlts|6|tfus_vor.lua|$V2|$MLTS|study = { type = 'cfl', cfl_random = true, cfl_decrement = 0.01 }|$NOOUT|time.final_time = 1.0; time.n_steps = 29
EOT
./run_set.sh jobs_vor_1.txt out ./wave_propagation_hf $I
dt=$(awk '!/^#/ {print $3; exit}' out/vor2/cfl_mlts/c.dat 2>/dev/null)
if [ -z "$dt" ]; then echo "$(date +%T) no stable step for the Voronoi mesh"; exit 1; fi
n=$(awk -v dt="$dt" 'BEGIN { n = 80.0 / (0.9 * dt); printf "%d", (n == int(n)) ? n : int(n) + 1 }')
echo "$(date +%T) Voronoi MLTS: stable dt $dt, production n_steps $n"
MEM="local dt_mem = time.final_time / time.n_steps; time.final_time = 3 * dt_mem; time.n_steps = 3; output.algebraic_memory = true; output.snapshots = {}; output.sensors = {}; output.energy = false; output.peak_samples = 0"
cat > jobs_vor_2.txt <<EOT
# made by chain_V.sh: MLTS on the fine Voronoi mesh, 80 us at 0.9 of the stable step ($dt us): $n steps
vor2/plate_mlts|6|tfus_vor.lua|$V2|$MLTS; time.n_steps = $n
vor2/mem_mlts|6|tfus_vor.lua|$V2|$MLTS; time.n_steps = $n|$MEM
EOT
./run_set.sh jobs_vor_2.txt out ./wave_propagation_hf $I
