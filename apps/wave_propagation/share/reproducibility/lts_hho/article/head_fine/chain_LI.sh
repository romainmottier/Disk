#!/bin/bash
# 2026-11-18 (user: "relance le voronoi avec plus de memoire"): LI on the fine Voronoi mesh, alone, MemoryMax 64G, with
# the binary wave_propagation_hf2 (option solver.release_assembly, off by default: the assembled LHS and COUPLING are
# freed and malloc_trim returns the memory of the setup once the matrix-free operator is built; results unchanged, checked
# on the coarse Voronoi mesh in probe_malloc/: peak 19.2 -> 14.6 GB there). The first attempt (vor2/mem_li, old binary)
# was killed at 60G during the factorization (kept in out/vor2/mem_li_oom_60G_old_binary).
#   L1: stable step of LI (random data, 1 us, energy growth < 5 %), as vor2/cfl_mlts for MLTS
#   L2: memory probe (3 steps, algebraic memory) at 0.9 of the stable step (80 us)
#   T : timed runs over the first tenth of the steps (8 us), three times (out_t, out_t_r2, out_t_r3), six threads, as
#       vor/plate_mlts_8us (nothing else on the machine during T)
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
V2="mesh.files = { \"$M/head_s80_vor2_plate.txt\" }"
LI="time.li_solver = 'pardiso'; time.li_refinement = 0; time.fine_size = 0.375; solver.release_assembly = true"
NOOUT="output.snapshots = {}; output.sensors = {}; output.energy = false; output.peak_samples = 0"
EXE=./wave_propagation_hf2
echo "$(date +%T) chain_LI start"
cat > jobs_vor_li1.txt <<EOT
# made by chain_LI.sh: stable step of LI on the fine Voronoi mesh
vor2/cfl_li|6|tfus_vor.lua|$V2|$LI|study = { type = 'cfl', cfl_random = true, cfl_decrement = 0.01 }|$NOOUT; output.profile = false|time.final_time = 1.0; time.n_steps = 29
EOT
./run_set.sh jobs_vor_li1.txt out $EXE $I
dt=$(awk '!/^#/ {print $3; exit}' out/vor2/cfl_li/c.dat 2>/dev/null)
if [ -z "$dt" ]; then echo "$(date +%T) no stable step for LI on the Voronoi mesh"; exit 1; fi
n=$(awk -v dt="$dt" 'BEGIN { n = 80.0 / (0.9 * dt); printf "%d", (n == int(n)) ? n : int(n) + 1 }')
k=$(( (n + 5) / 10 ))
echo "$(date +%T) Voronoi LI: stable dt $dt, production n_steps $n, timed runs $k steps"
MEM="local dt_mem = time.final_time / time.n_steps; time.final_time = 3 * dt_mem; time.n_steps = 3; output.algebraic_memory = true; $NOOUT"
cat > jobs_vor_li2.txt <<EOT
# made by chain_LI.sh: LI on the fine Voronoi mesh, 80 us at 0.9 of the stable step ($dt us): $n steps
vor2/mem_li|6|tfus_vor.lua|$V2|$LI; time.n_steps = $n|$MEM
EOT
if [ -d out/vor2/mem_li ] && [ ! -f out/vor2/mem_li/time.txt ] && [ ! -d out/vor2/mem_li_oom_60G_old_binary ]; then
  mv out/vor2/mem_li out/vor2/mem_li_oom_60G_old_binary
fi
./run_set.sh jobs_vor_li2.txt out $EXE $I
echo "$(date +%T) LI MEMORY PROBE DONE, TIMED LI RUNS START"
cat > jobs_timed_vli.txt <<EOT
# made by chain_LI.sh: Voronoi LI, first $k of the $n steps (8 us), six threads
vor/plate_li_8us|6|tfus_vor.lua|$V2|$LI; time.n_steps = $n|local dt_x = time.final_time / time.n_steps; time.final_time = $k * dt_x; time.n_steps = $k|$NOOUT
EOT
for o in out_t out_t_r2 out_t_r3; do
  ./run_set.sh jobs_timed_vli.txt $o $EXE $I
done
echo "$(date +%T) chain_LI done"
