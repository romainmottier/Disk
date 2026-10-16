#!/usr/bin/env python3
"""Cost of the unfitted HHO methods with quadratures on the curved geometry
(-p0 polynomial extension, -p1 cell agglomeration, exact interface) against
the boundary value correction on the flat interface (-p2 agglomeration, -p3
polynomial extension); in 2D, also -p0, -p1 with the polyline of 2^8 segments
per cut cell of ProtoN and of the paper (-g0 -r8).

    python3 compare_bvc_speed.py <unfitted_HHO> <outdir> [--threads N] [--run] [--plot]

--run runs the experiments one at a time (<outdir>/<case>_p<p>_k<k>), --plot
writes <outdir>/bvc_speed.md and <outdir>/figures/bvc_speed_*.png: wall times
of the geometry of the cut (mesh), of the assembly with static condensation,
of the solve, and the H1 error (own unknowns of the subcells) against the
time of mesh + assembly + solve. Test cases common to the four prototypes:
contrast rho^6 with kappa_2/kappa_1 = 1e4 (circle, sphere; -t1 for -p0/-p1,
-t0 for -p2/-p3) and the jumps (59) on the flower of the BVC paper (-t4 -i7).
"""
import argparse, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CASES = {   # name: (mesh, base N, levels per k, problem)
    '2d_circle':     (0, 8,  {1: 4, 2: 4, 3: 4}, 'contrast'),
    '2d_flower':     (0, 16, {1: 3, 2: 3, 3: 3}, 'jumps'),
    '3d_sphere_hex': (3, 4,  {1: 3, 2: 3, 3: 2}, 'contrast'),
    '3d_sphere_tet': (4, 4,  {1: 3, 2: 2, 3: 2}, 'contrast'),
    '3d_flower_hex': (3, 8,  {1: 2, 2: 2, 3: 1}, 'jumps'),
}
METHODS = {0: 'extension, curved', 1: 'agglomeration, curved', 2: 'BVC + agglomeration', 3: 'BVC + extension',
           4: 'extension, polyline 2^8', 5: 'agglomeration, polyline 2^8'}


def run_dir(outdir, case, p, k):
    """-p4, -p5: -p0, -p1 with the polyline of 2^8 segments of ProtoN (2D)"""
    if p >= 4:
        return os.path.join(outdir, f'{case}_r8_p{p-4}_k{k}')
    return os.path.join(outdir, f'{case}_p{p}_k{k}')


def methods_of(case):
    return [0, 1, 2, 3, 4, 5] if case.startswith('2d') else [0, 1, 2, 3]


def args_of(case, p, k):
    M, N0, L, prob = CASES[case]
    q = p - 4 if p >= 4 else p
    if prob == 'contrast':
        t = ('-t1' if q in (0, 1) else '-t0') + ' -m4 -i0'
    else:
        t = '-t4 -i7'
    extra = ' -g0 -r8' if p >= 4 else ''
    return f'-p{q} -M{M} -x0 {t} -N{N0} -K{k} -k{k} -l{L[k]}{extra}'


ANSI = re.compile(r'\x1b\[[0-9;]*m')


def parse(d):
    """levels: list of dict(l, h, ndof, H1, L2, t_mesh, t_asm, t_solve, t_err)"""
    out = os.path.join(d, 'out.txt')
    if not os.path.exists(out):
        return []
    levels, cur = [], None
    for line in open(out, errors='replace'):
        line = ANSI.sub('', line)
        m = re.search(r'Space refinment level -l : (\d+)', line)
        if m:
            cur = dict(l=int(m.group(1)))
            levels.append(cur)
            continue
        if cur is None:
            continue
        m = re.search(r'(\d+) face unknowns', line)
        if m:
            cur['ndof'] = int(m.group(1))
        m = re.search(r'Times \(s\): mesh ([\d.e+-]+), assembly ([\d.e+-]+)(?: \(cut cells [\d.e+-]+\))?, solve ([\d.e+-]+), errors ([\d.e+-]+)', line)
        if m:
            cur['t_mesh'], cur['t_asm'], cur['t_solve'], cur['t_err'] = map(float, m.groups())
    # errors and h from solution_error_file.txt (own unknowns of the subcells)
    ef = os.path.join(d, 'solution_error_file.txt')
    if os.path.exists(ef):
        lvl = None
        vals = {}
        for line in open(ef):
            m = re.match(r'Space refinment level -l : (\d+)', line)
            if m:
                lvl = int(m.group(1))
                vals.setdefault(lvl, {})
                continue
            if lvl is None or '=' not in line:
                continue
            key, v = line.split('=', 1)
            try:
                vals[lvl][key.strip()] = float(v)
            except ValueError:
                pass
        for lv in levels:
            v = vals.get(lv['l'], {})
            lv['h'] = v.get('Characteristic h size', float('nan'))
            lv['H1'] = v.get('H1-norm error with own unknowns on ill-cut subcells', v.get('H1-norm error', float('nan')))
            lv['L2'] = v.get('L2-norm error with own unknowns on ill-cut subcells', v.get('L2-norm error', float('nan')))
    return [lv for lv in levels if 't_solve' in lv]


def run(binary, outdir, threads):
    env = dict(os.environ, LC_ALL='C', OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    for case in CASES:
        for k in [1, 2, 3]:
            for p in methods_of(case):
                d = run_dir(outdir, case, p, k)
                if parse(d):
                    continue
                os.makedirs(d, exist_ok=True)
                a = args_of(case, p, k)
                open(os.path.join(d, 'args.txt'), 'w').write(a + '\n')
                print(case, p, k, a, flush=True)
                subprocess.run([binary] + a.split(), cwd=d, env=env,
                               stdout=open(os.path.join(d, 'out.txt'), 'w'), stderr=subprocess.STDOUT)


def plot(outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
    MARK = ['o', 's', '^', 'D', 'v', 'P']
    INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
    plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
                         'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                         'text.color': INK, 'font.size': 9, 'axes.titlesize': 9.5, 'legend.fontsize': 7.5,
                         'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
                         'lines.linewidth': 1.5, 'lines.markersize': 5, 'legend.frameon': False})
    figs = os.path.join(outdir, 'figures')
    os.makedirs(figs, exist_ok=True)
    md = ['# Cost: quadratures on the curved geometry against the boundary value correction\n',
          'Wall times (s) per mesh: geometry of the cut (mesh), assembly with static condensation, solve; '
          'H1 error with the own unknowns of the subcells. -p0/-p1: exact (curved) interface, '
          '-p2/-p3: flat interface with boundary value correction.\n']
    data = {}
    for case in CASES:
        for k in [1, 2, 3]:
            for p in methods_of(case):
                data[(case, k, p)] = parse(run_dir(outdir, case, p, k))
    for case in CASES:
        md.append(f'## {case}\n')
        for k in [1, 2, 3]:
            if not any(data[(case, k, p)] for p in methods_of(case)):
                continue
            md.append(f'### k = {k}\n')
            md.append('| method | h | ndof | geometry | assembly | solve | total | H1 |')
            md.append('|---|---|---|---|---|---|---|---|')
            for p in methods_of(case):
                q = p - 4 if p >= 4 else p
                for lv in data[(case, k, p)]:
                    tot = lv['t_mesh'] + lv['t_asm'] + lv['t_solve']
                    md.append(f"| -p{q} {METHODS[p]} | {lv.get('h', float('nan')):.3e} | {lv.get('ndof', 0)} | "
                              f"{lv['t_mesh']:.3g} | {lv['t_asm']:.3g} | {lv['t_solve']:.3g} | {tot:.3g} | "
                              f"{lv.get('H1', float('nan')):.3e} |")
            md.append('')
        # figure: H1 error against time, one panel per k
        fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.6))
        for ax, k in zip(axs, [1, 2, 3]):
            for p in methods_of(case):
                lv = data[(case, k, p)]
                if not lv:
                    continue
                q = p - 4 if p >= 4 else p
                t = [x['t_mesh'] + x['t_asm'] + x['t_solve'] for x in lv]
                e = [x.get('H1', float('nan')) for x in lv]
                ax.loglog(t, e, color=SLOTS[p], marker=MARK[p], linestyle='-' if p in (2, 3) else '--',
                          label=f'-p{q} {METHODS[p]}', markeredgecolor=SURF, markeredgewidth=0.6)
            ax.set_title(f'{case}, k = {k}', loc='left')
            ax.set_xlabel('wall time: geometry + assembly + solve (s)')
            ax.set_ylabel(r'$|u - u_h|_{H^1}$')
        axs[0].legend(loc='lower left')
        fig.tight_layout()
        fig.savefig(os.path.join(figs, f'bvc_speed_{case}.png'), dpi=160)
        fig.savefig(os.path.join(figs, f'bvc_speed_{case}.pdf'))
        plt.close(fig)
        md.append(f'![{case}](figures/bvc_speed_{case}.png)\n')
    open(os.path.join(outdir, 'bvc_speed.md'), 'w').write('\n'.join(md) + '\n')
    print('written', os.path.join(outdir, 'bvc_speed.md'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('binary')
    ap.add_argument('outdir')
    ap.add_argument('--threads', type=int, default=12)
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--plot', action='store_true')
    a = ap.parse_args()
    if a.run:
        run(os.path.abspath(a.binary), a.outdir, a.threads)
    if a.plot or not a.run:
        plot(a.outdir)
