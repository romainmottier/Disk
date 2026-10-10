#!/bin/bash
#
# Numerical results of Section 5 of "Unfitted hybrid high-order methods
# stabilized by polynomial extension for elliptic interface problems"
# (Burman, Ern, Mottier): stabilization by polynomial extension (-p0) and,
# for comparison, by cell agglomeration (-p1, directories *_agglo).
#
# Usage: run_paper.sh <path/to/unfitted_HHO> [output_dir]
# Then:  python3 plot_paper_results.py [output_dir]
#
# Each run is performed in its own directory, which contains the log, the
# errors (solution_error_file.txt) and the other outputs of the run.

set -u

# decimal point independent of the user's locale (awk, printf)
export LC_ALL=C

if [ $# -lt 1 ]; then
    echo "Usage: $0 <path/to/unfitted_HHO> [output_dir]"
    exit 1
fi

BIN=$(realpath "$1")
OUT=${2:-paper_results}
mkdir -p "$OUT"
cd "$OUT" || exit 1

# The paper uses the parameters of ProtoN, which are not the defaults of the
# app: eta_N = (k+1)^2 (-e -1), the scaled monomial bases centered on the
# subcells (-b 1), the choice of the neighbors of the ill-cut cells of ProtoN
# (-P 0), the full system without static condensation, solved by
# SimplicialLDLT (-c 0 -s 0), the polyline interface of 2^r segments (-g 0)
# and local operators computed in every cell (-u 0). A later -e, -b, -P, -c,
# -s, -g or -u in <options> overrides them.

# run <dir> <options>: polynomial extension in <dir>
run() {
    local dir=$1
    shift
    mkdir -p "$dir"
    echo "  $dir"
    (cd "$dir" && "$BIN" -p0 -e -1 -b 1 -P 0 -c 0 -s 0 -g 0 -u 0 "$@" > log.txt 2>&1) || echo "    FAILED (see $dir/log.txt)"
}

# run_agglo <dir> <options>: cell agglomeration in <dir>_agglo
run_agglo() {
    local dir=${1}_agglo
    shift
    mkdir -p "$dir"
    echo "  $dir"
    (cd "$dir" && "$BIN" -p1 -b 1 -P 0 -c 0 -s 0 -g 0 -u 0 "$@" > log.txt 2>&1) || echo "    FAILED (see $dir/log.txt)"
}

echo "Fig. 5: coarsest meshes, pairing and agglomeration (silo, okc and cells_l0.txt files)"
for i in 0 1; do
    name=fig5_circle
    [ $i = 1 ] && name=fig5_flower
    run       $name -i$i -t0 -k0 -l0 -r8 -f1
    run_agglo $name -i$i -t0 -k0 -l0 -r8 -f1
done

echo "Fig. 6: convergence, solution (56)"
for i in 0 1; do
    name=fig6_circle
    [ $i = 1 ] && name=fig6_flower
    run       $name -i$i -t0 -k3 -l4 -r8
    run_agglo $name -i$i -t0 -k3 -l4 -r8
done

echo "Fig. 6 with eta_N = 20 (value given in Section 5.1 of the paper)"
run fig6_circle_etaN20 -i0 -t0 -k3 -l4 -r8 -e 20
run fig6_flower_etaN20 -i1 -t0 -k3 -l4 -r8 -e 20

echo "Fig. 7: flagging parameter and centering of the bases, k = 3"
for a in 0 0.1 0.2 0.3; do
    for b in 1 0; do
        name=fig7_a${a}_subcell
        [ $b = 0 ] && name=fig7_a${a}_cell
        run       $name -i0 -t0 -k3 -K3 -l4 -r8 -a $a -b $b
        run_agglo $name -i0 -t0 -k3 -K3 -l4 -r8 -a $a -b $b
    done
done

echo "Fig. 8: block sparsity pattern, k = 0, coarsest mesh"
run       fig8 -i0 -t0 -k0 -l0 -r8 -z1
run_agglo fig8 -i0 -t0 -k0 -l0 -r8 -z1

echo "Fig. 9 left: solution (57), contrast kappa_2 = 10^m, finest mesh"
for m in 0 1 2 3 4; do
    run fig9_left_m$m -i0 -t1 -k3 -l4 -L4 -r8 -m $m
done

echo "Fig. 9 right: solution (57), sub-triangulation parameter r, k = 3"
for r in 2 4 6 8; do
    for m in 0 1 2 3 4; do
        run fig9_right_r${r}_m$m -i0 -t1 -k3 -K3 -l4 -r $r -m $m
    done
done

echo "Fig. 10: solutions (58a) and (58b), kappa_2 = 10^4"
run fig10_left  -i0 -t2 -k3 -l4 -r8 -m4
run fig10_right -i0 -t3 -k3 -l4 -r8 -m4

echo "Fig. 11: solution (59)"
run fig11_left -i0 -t4 -k3 -l4 -r10
for r in 4 6 8 10 11; do
    run fig11_right_r$r -i0 -t4 -k3 -K3 -l4 -r $r
done

echo "Fig. 12 left: condition number, circle of radius 1/3 + i/32"
for i in -4 -3 -2 -1 0 1 2 3 4; do
    R=$(awk -v i=$i 'BEGIN { printf "%.17g", 1.0/3.0 + i/32.0 }')
    run       fig12_left_i$i -i0 -t0 -k3 -l0 -r8 -R $R -q1
    run_agglo fig12_left_i$i -i0 -t0 -k3 -l0 -r8 -R $R -q1
done

echo "Fig. 12 right: condition number, square interface shifted by 0.5*10^-p"
for p in 2 3 4 5 6 7 8 9; do
    run       fig12_right_p$p -i2 -t0 -k3 -l0 -r8 -d $p -q1
    run_agglo fig12_right_p$p -i2 -t0 -k3 -l0 -r8 -d $p -q1
done
