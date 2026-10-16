#!/usr/bin/env python3
"""Numerical experiments of the boundary value correction (unfitted_HHO -p2).

    python3 run_bvc.py <path/to/unfitted_HHO> <outdir> [--set 2d|3d|all]
                       [--prototype 2|3] [--mem GB] [--threads N] [--list]

--prototype 3 runs the same experiments with the stabilization by polynomial
extension (-p3) instead of cell agglomeration (-p2).

Each run is done in <outdir>/<name> (bvc_convergence.txt, bvc_geometry.txt,
solution_error_file.txt, out.txt). The runs are scheduled within a memory
budget (estimated from the number of face unknowns: the 3D runs with k = 2, 3
on the finest meshes need 10-20 GB each) and a thread budget; the finished
runs (exit 0 in <outdir>/done.txt) are skipped when the script is restarted.

2D (quadrangles, ProtoN pairing -P0: the errors are those of ProtoN's
cuthho_bvc_results.txt): Tables 1-10 of the paper, i.e. the full method on
the circle (kappa_2 = 1e4, h = 1/8..1/128) and the flower (h = 1/16..1/128),
the corrections one at a time, the consistency benchmark, the interface
refinement without correction, the contrast, the Taylor order, the patch
test, the quadrature check and the geometry; full BVC on triangles.

3D (hexahedra, tetrahedra): sphere (kappa_2 = 1e4) and 3D flower
(discontinuous solution), on two series of meshes (N = 4, 8, 16, 32 and
N = 5, 10, 20 for the sphere, N = 8, 16, 32 and N = 6, 12, 24 for the flower)
to sample the curves (the finest meshes for k = 3 on tetrahedra need 22 and
46 GB: they are skipped if --mem is smaller); full method and no correction (and the corrections one
at a time on the flower, hexahedra), contrast, patch test, geometry.
"""
import argparse, os, subprocess, sys, time


def fbs(k):
    return (k + 1) * (k + 2) // 2


def memory(M, k, N):
    """estimated memory (GB) of a run on the meshes N^3 (Pardiso LU)"""
    if M in (0, 1):
        return 1.0
    if M == 3:   # hexahedra: k = 2, N = 32, 585e3 face unknowns, 13.8 GB
        dofs = fbs(k) * 3 * N**3
        return 13.8 * (dofs / 585e3)**(4/3) * (fbs(k) / 6)**(2/3) + 0.3
    dofs = fbs(k) * 12 * N**3   # tetrahedra: k = 1, N = 32, 1.17e6, 13 GB
    return 13.0 * (dofs / 1.17e6)**(4/3) * (fbs(k) / 3)**(2/3) + 0.3


def jobs_2d():
    J = []
    def add(name, args):
        J.append((name, '-p2 -x0 ' + args, 1.0, 2))
    C = '-M0 -P0 -i0 -t0 -m4 -N8 -K0 -k3 -l4'      # circle, N = 8..128
    F = '-M0 -P0 -i7 -t1 -m0 -N16 -K0 -k3 -l3'     # flower, N = 16..128
    for v in ['full', 'naive', 'taylor', 'gam', 'ext']:
        add(f'circle_{v}', f'{C} -v {v}')
        add(f'flower_{v}', f'{F} -v {v}')
    for n in [2, 3, 4]:
        add(f'circle_naive_n{n}', f"{C.replace('-K0', '-K2')} -v naive -n{n}")
        add(f'flower_naive_n{n}', f"{F.replace('-K0', '-K2')} -v naive -n{n}")
    for m in [0, 2, 4, 6]:
        add(f'contrast_m{m}', f'-M0 -P0 -i0 -t0 -m{m} -N8 -K0 -k3 -L3 -l3 -v full')
    for T in [0, 1, 2]:
        add(f'circle_T{T}', f"{C.replace('-K0', '-K1')} -v gam -T{T}")
        add(f'flower_T{T}', f"{F.replace('-K0', '-K1')} -v gam -T{T}")
    add('circle_full_P1', f"{C.replace('-P0', '-P1')} -v full")
    add('flower_full_P1', f"{F.replace('-P0', '-P1')} -v full")
    add('circle_full_Q4', f'{C} -v full -Q8')
    add('flower_full_Q4', f'{F} -v full -Q8')
    add('patch2d', '-M0 -P0 -i8 -t2 -N4 -K0 -k3 -l2 -v full')
    add('patch2d_tri', '-M1 -i8 -t2 -N4 -K0 -k3 -l2 -v full')
    add('geo_circle', '-M0 -P0 -i0 -t3 -N8 -K0 -k1 -l4')
    add('geo_flower', '-M0 -P0 -i7 -t3 -N8 -K0 -k1 -l4')
    add('circle_full_tri', '-M1 -i0 -t0 -m4 -N8 -K0 -k3 -l4 -v full')
    add('flower_full_tri', '-M1 -i7 -t1 -m0 -N16 -K0 -k3 -l3 -v full')
    add('flower_naive_tri', '-M1 -i7 -t1 -m0 -N16 -K0 -k3 -l3 -v naive')
    return J


def jobs_3d():
    J = []
    def add(name, M, prob, var, k, N0, L):
        m = memory(M, k, N0 * 2**L)
        J.append((f'{name}_M{M}_{var}_k{k}_N{N0}',
                  f'-p2 -M{M} -x0 {prob} -N{N0} -K{k} -k{k} -l{L} -v {var}', m, 12 if m > 6 else 3))
    S = '-i0 -t0 -m4'    # sphere, u_i = rho^6/kappa_i
    F = '-i7 -t1 -m0'    # 3D flower, discontinuous solution
    for var in ['full', 'naive']:
        for k in range(4):
            add('sphere', 3, S, var, k, 4, 2 if k == 3 else 3)
            add('sphere', 3, S, var, k, 5, 2)
            add('sphere', 4, S, var, k, 4, 3 if k <= 1 else 2)
            add('sphere', 4, S, var, k, 5, 2 if k <= 1 else 1)
    for var in ['full', 'naive', 'taylor', 'gam']:
        for k in range(4):
            add('flower', 3, F, var, k, 8, 1 if k == 3 else 2)
            add('flower', 3, F, var, k, 6, 2)
    for var in ['full', 'naive']:
        for k in range(4):
            add('flower', 4, F, var, k, 8, 2 if k <= 1 else 1)
            add('flower', 4, F, var, k, 6, 2 if k <= 1 else 1)
            if k == 2:
                add('flower', 4, F, var, k, 12, 1)
        # finest meshes for k = 3 on tetrahedra (about 46 GB and 22 GB)
        add('flower', 4, F, var, 3, 24, 0)
        add('sphere', 4, S, var, 3, 20, 0)
    for m in [0, 2, 4, 6]:
        J.append((f'contrast_M3_m{m}', f'-p2 -M3 -x0 -i0 -t0 -m{m} -N16 -K0 -k3 -l0 -v full', memory(3, 3, 16), 3))
        J.append((f'contrast_M4_m{m}', f'-p2 -M4 -x0 -i0 -t0 -m{m} -N8 -K0 -k3 -l0 -v full', memory(4, 3, 8), 3))
    for M in [3, 4]:
        J.append((f'patch_M{M}', f'-p2 -M{M} -x0 -i8 -t2 -N4 -K0 -k3 -l2', memory(M, 3, 16), 12))
        J.append((f'geo_sphere_M{M}', f'-p2 -M{M} -x0 -i0 -t3 -N4 -K0 -k1 -l3', 2, 2))
        J.append((f'geo_flower_M{M}', f'-p2 -M{M} -x0 -i7 -t3 -N8 -K0 -k1 -l2', 2, 2))
        J.append((f'quad_sphere_M{M}', f'-p2 -M{M} -x0 -i0 -t0 -m4 -N4 -K0 -k3 -l1 -v full -Q4', memory(M, 3, 8), 3))
        J.append((f'quad0_sphere_M{M}', f'-p2 -M{M} -x0 -i0 -t0 -m4 -N4 -K0 -k3 -l1 -v full', memory(M, 3, 8), 3))
    return J


def run(jobs, binary, outdir, mem_budget, thr_budget):
    os.makedirs(outdir, exist_ok=True)
    done_path = os.path.join(outdir, 'done.txt')
    done = set()
    if os.path.exists(done_path):
        for l in open(done_path):
            t = l.split()
            if len(t) >= 3 and t[1] == 'exit' and t[2] == '0':
                done.add(t[0])
    todo = sorted([j for j in jobs if j[0] not in done], key=lambda j: -j[2])
    for j in [j for j in todo if j[2] > mem_budget]:
        print(f'{j[0]}: skipped, needs about {j[2]:.0f} GB (--mem {mem_budget:g})', flush=True)
    todo = [j for j in todo if j[2] <= mem_budget]
    running, used_mem, used_thr = {}, 0.0, 0
    while todo or running:
        for j in list(todo):
            name, args, mem, thr = j
            if running and (used_mem + mem > mem_budget or used_thr + thr > thr_budget):
                continue
            d = os.path.join(outdir, name)
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, 'args.txt'), 'w').write(args + '\n')
            env = dict(os.environ, LC_ALL='C', OMP_NUM_THREADS=str(thr), MKL_NUM_THREADS=str(thr))
            p = subprocess.Popen([binary] + args.split(), cwd=d, env=env,
                                 stdout=open(os.path.join(d, 'out.txt'), 'w'), stderr=subprocess.STDOUT)
            running[p] = (j, time.time())
            used_mem += mem
            used_thr += thr
            todo.remove(j)
        time.sleep(2)
        for p in list(running):
            if p.poll() is None:
                continue
            (name, args, mem, thr), t0 = running.pop(p)
            used_mem -= mem
            used_thr -= thr
            with open(done_path, 'a') as f:
                f.write(f'{name} exit {p.returncode} {time.time() - t0:.0f}s\n')
            print(f'{name}: exit {p.returncode} ({time.time() - t0:.0f} s), {len(todo)} to start', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('binary')
    ap.add_argument('outdir')
    ap.add_argument('--set', choices=['2d', '3d', 'all'], default='all')
    ap.add_argument('--prototype', type=int, choices=[2, 3], default=2,
                    help='2: cell agglomeration, 3: polynomial extension')
    ap.add_argument('--mem', type=float, default=40, help='memory budget (GB)')
    ap.add_argument('--threads', type=int, default=os.cpu_count(), help='thread budget')
    ap.add_argument('--list', action='store_true', help='print the runs and exit')
    a = ap.parse_args()
    J = []
    if a.set in ('2d', 'all'):
        J += [(n, x, m, t) for n, x, m, t in jobs_2d()]
    if a.set in ('3d', 'all'):
        J += [('3d_' + n, x, m, t) for n, x, m, t in jobs_3d()]
    if a.prototype == 3:
        # default pairing (largest subcell) for all the runs; the full method
        # is also run with the pairing of ProtoN (*_full_P0)
        K = []
        for n, x, m, t in J:
            if n.endswith('_P1'):
                continue
            x = x.replace('-p2', '-p3')
            if n in ('circle_full', 'flower_full'):
                K.append((n + '_P0', x, m * 1.3, t))
            K.append((n, x.replace(' -P0', ''), m * 1.3, t))
        J = K
    if a.list:
        for n, x, m, t in J:
            print(f'{n:28s} {m:5.1f} GB {t:2d} threads   {x}')
        sys.exit(0)
    run(J, os.path.abspath(a.binary), a.outdir, a.mem, a.threads)
