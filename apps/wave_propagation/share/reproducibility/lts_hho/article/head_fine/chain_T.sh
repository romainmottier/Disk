#!/bin/bash
# Timed runs of the fine head model (2026-11-16, Section 5.4 on the meshes twice finer), under the rules of the
# preprint: every run ALONE on the machine (nothing else: no other unit, no numpy, no LaTeX build), three times
# (out_t, out_t_r2, out_t_r3), best of three; six threads (practice) and one thread (work).
#   - ERK(4) over 2 us (user, 2026-11-16) at 0.9 of its stable step (out/P/cfl_erk), extrapolated x 40;
#   - MLTS and LI over the first 278 of the 2788 macro steps of the full runs (8 us), extrapolated x 2788/278, the
#     factorization of LI ("fine system") counted once;
#   - Voronoi mesh: MLTS only (LI does not fit in memory), first tenth of the steps of vor2/plate_mlts, six threads.
# Then the layout of make_tfus.py --short in runs_view (fit/, t1/, vor/, with _r2 and _r3).
#   usage: chain_T.sh q [P] (quadtree: ERK(4), MLTS, LI on six and one threads, plate mesh head_s80_P_plate.txt,
#                            needs out/P/cfl_erk; P = fit2 by default, fit2g = graded towards the titanium)
#          chain_T.sh v   (Voronoi: MLTS on six threads, needs jobs_vor_2.txt of chain_V1.sh)
part=${1:?usage: chain_T.sh q|v [P]}
P=${2:-fit2}
cd "$(dirname "$0")"
I=/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/share/reproducibility/lts_hho/article/inputs
M=/home/mottie0000/Github/Diskpp/lts_work/head_fine/meshes
F2="mesh.files = { \"$M/head_s80_fit2_intact.txt\", \"$M/head_s80_${P}_plate.txt\" }; mesh.level = 2"
V2="mesh.files = { \"$M/head_s80_vor2_plate.txt\" }"
MLTS="time.scheme = 'mlts'; time.level_size = 'width'; time.mlts_reuse = true; time.level_ratios = {2, 2, 4, 8}"
LI="time.li_solver = 'pardiso'; time.li_refinement = 0; time.fine_size = 0.375"
NOOUT="output.snapshots = {}; output.sensors = {}; output.energy = false; output.peak_samples = 0"
if [ "$part" = q ]; then
dt=$(awk '!/^#/ {print $3; exit}' out/$P/cfl_erk/c.dat 2>/dev/null)
[ -n "$dt" ] || { echo "$(date +%T) no stable step for ERK(4) (out/$P/cfl_erk)"; exit 1; }
nerk=$(awk -v dt="$dt" 'BEGIN { n = 2.0 / (0.9 * dt); printf "%d", (n == int(n)) ? n : int(n) + 1 }')
S8="time.final_time = 278 * 80.0 / 2788; time.n_steps = 278"
cat > jobs_timed_q.txt <<EOT
# made by chain_T.sh q $P: ERK(4) stable step $dt us -> $nerk steps over 2 us (plate mesh head_s80_${P}_plate.txt)
fit/plate_erk_2us|6|tfus_fit.lua|$F2|time.scheme = 'erk'; time.final_time = 2.0; time.n_steps = $nerk|$NOOUT
fit/plate_mlts_8us|6|tfus_fit.lua|$F2|$MLTS; $S8|$NOOUT
fit/plate_li_8us|6|tfus_fit.lua|$F2|$LI; $S8|$NOOUT
t1/erk_2us|1|tfus_fit.lua|$F2|time.scheme = 'erk'; time.final_time = 2.0; time.n_steps = $nerk|$NOOUT
t1/mlts_8us|1|tfus_fit.lua|$F2|$MLTS; $S8|$NOOUT
t1/li_8us|1|tfus_fit.lua|$F2|$LI; $S8|$NOOUT
EOT
jobs=jobs_timed_q.txt
else
nv=$(sed -n 's/.*vor2\/plate_mlts|.*time.n_steps = \([0-9]*\).*/\1/p' jobs_vor_2.txt 2>/dev/null | head -1)
[ -n "$nv" ] || { echo "$(date +%T) no step count for the Voronoi mesh (jobs_vor_2.txt)"; exit 1; }
kv=$(( (nv + 5) / 10 ))
cat > jobs_timed_v.txt <<EOT
# made by chain_T.sh v: Voronoi MLTS, first $kv of the $nv steps of vor2/plate_mlts
vor/plate_mlts_8us|6|tfus_vor.lua|$V2|$MLTS; time.n_steps = $nv|local dt_x = time.final_time / time.n_steps; time.final_time = $kv * dt_x; time.n_steps = $kv|$NOOUT
EOT
jobs=jobs_timed_v.txt
fi
for o in out_t out_t_r2 out_t_r3; do
  ./run_set.sh $jobs $o ./wave_propagation_hf $I
done
echo "$(date +%T) timed runs ($jobs) done"
