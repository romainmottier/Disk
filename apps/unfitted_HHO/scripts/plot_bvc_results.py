#!/usr/bin/env python3
"""Tables and figures of the boundary value correction experiments.

    python3 plot_bvc_results.py <results dir of run_bvc.py> [--proton cuthho_bvc_results.txt]
                                [--compare <other results dir> [--labels a,b]]

writes <dir>/bvc_tables.md (tables, least-squares rates, comparison with
ProtoN's results if given) and <dir>/figures/bvc*.png/pdf. With --compare
(e.g. the results of -p2 and of -p3), the full method of both directories is
compared (figures bvc*_compare, section of the tables).
"""
import argparse, math, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

HEAD = re.compile(r'#\s*(.*), k = (\d+), variant (\w+) \(m = (\d+), gamma (\d)\)')


def parse_table(path):
    """bvc_convergence.txt (or ProtoN's table) -> {(k, variant): [rows]}"""
    res, cur = {}, None
    if not os.path.exists(path):
        return res
    for line in open(path):
        line = line.strip()
        m = HEAD.match(line)
        if m:
            cur = (int(m.group(2)), m.group(3))
            res[cur] = []
            continue
        if not line or line.startswith('#') or cur is None:
            continue
        t = line.split()
        f = lambda x: float(x)
        if len(t) == 9:      # ProtoN: N h ndof H1 r L2 r EN r
            res[cur].append(dict(N=int(t[0]), h=f(t[1]), ndof=int(t[2]), H1=f(t[3]), L2=f(t[5]), EN=f(t[7])))
        elif len(t) == 10:   # diskpp: N l h ndof H1 r L2 r EN r
            res[cur].append(dict(N=int(t[0]), h=f(t[2]), ndof=int(t[3]), H1=f(t[4]), L2=f(t[6]), EN=f(t[8])))
    return res


def parse_proton(path):
    blocks, cur = {}, None
    for l in open(path):
        if l.startswith('# ') and ', k = ' in l:
            cur = l[2:].split('(')[0]
        if cur is not None:
            blocks.setdefault(cur, []).append(l)
    out = {}
    for p, lines in blocks.items():
        tmp = os.path.join('/dev/shm' if os.path.isdir('/dev/shm') else '.', '_bvc_proton.txt')
        open(tmp, 'w').write(''.join(lines))
        out[p] = parse_table(tmp)
        os.remove(tmp)
    return out


def ls_rate(rows, q):
    pts = [(r['h'], r[q]) for r in rows if r[q] > 0]
    if len(pts) < 2:
        return float('nan')
    h, e = np.array(pts).T
    return np.polyfit(np.log(h), np.log(e), 1)[0]


# ---------------------------------------------------------------- figures
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'text.color': INK, 'font.size': 9, 'axes.titlesize': 9.5, 'legend.fontsize': 8,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
    'lines.linewidth': 1.5, 'lines.markersize': 5, 'legend.frameon': False,
})
QLABEL = {'H1': r'$|u - u_h|_{H^1}$', 'L2': r'$\|u - u_h\|_{L^2}$', 'EN': 'energy error'}
QOPT = {'H1': 1, 'L2': 2, 'EN': 1}


def series(ax, rows, q, i, label, dashed=False):
    rows = [r for r in rows if r[q] > 0]
    if not rows:
        return
    ax.loglog([r['h'] for r in rows], [r[q] for r in rows], color=SLOTS[i], marker=MARKERS[i], label=label,
              linestyle='--' if dashed else '-', markeredgecolor=SURF, markeredgewidth=0.6)


def slope(ax, rows, q, p):
    """reference h^p below the last point"""
    rows = [r for r in rows if r[q] > 0]
    if len(rows) < 2:
        return
    h = np.array([rows[0]['h'], rows[-1]['h']])
    e = rows[-1][q] * 0.35
    ax.loglog(h, e * (h / h[-1])**p, color=INK2, linestyle=':', linewidth=1.0)
    ax.text(h[-1], e, f'  $h^{{{p}}}$', color=INK2, fontsize=7.5, va='center', ha='left')


def finish(ax, title, q):
    ax.set_title(title, loc='left')
    ax.set_xlabel('$h$ (largest cell diameter)')
    ax.set_ylabel(QLABEL[q])
    lines = [l for l in ax.get_lines() if len(l.get_xdata()) > 0]
    if lines:
        lo = min(min(l.get_xdata()) for l in lines)
        hi = max(max(l.get_xdata()) for l in lines)
        ticks = [t for t in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5] if lo / 1.3 <= t <= hi * 1.3]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.margins(x=0.12)
    ax.invert_xaxis()


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, name + '.png'), dpi=160)
    fig.savefig(os.path.join(outdir, name + '.pdf'))
    plt.close(fig)


M_OF_K = [0, 1, 1, 2, 2]


# ---------------------------------------------------------------- data access
class Results:
    def __init__(self, d):
        self.d = d

    def t2(self, name):
        return parse_table(os.path.join(self.d, name, 'bvc_convergence.txt'))

    def t3(self, prob, M, var, k, series=None):
        """3D: rows of the meshes of all the series, sorted by decreasing h"""
        rows = []
        for n in sorted(os.listdir(self.d)):
            m = re.match(rf'3d_{prob}_M{M}_{var}_k{k}_N(\d+)$', n)
            if m and (series is None or int(m.group(1)) in series):
                rows += parse_table(os.path.join(self.d, n, 'bvc_convergence.txt')).get((k, var), [])
        uniq = {}
        for r in rows:
            uniq[r['N']] = r
        return sorted(uniq.values(), key=lambda r: -r['h'])


def fmt_rows(rows, qs=('H1', 'L2', 'EN')):
    out = ['| N | h | ndof | ' + ' | '.join(f'{q} | rate' for q in qs) + ' |',
           '|---|---|---|' + '---|---|' * len(qs)]
    for i, r in enumerate(rows):
        cells = [str(r['N']), f"{r['h']:.3e}", str(r['ndof'])]
        for q in qs:
            cells.append(f"{r[q]:.3e}" if r[q] > 0 else '-')
            if i == 0 or r[q] <= 0 or rows[i-1][q] <= 0:
                cells.append('-')
            else:
                cells.append(f"{math.log(rows[i-1][q] / r[q]) / math.log(rows[i-1]['h'] / r['h']):.2f}")
        out.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('--proton', default=None, help="ProtoN's cuthho_bvc_results.txt")
    ap.add_argument('--compare', default=None, help='second results directory')
    ap.add_argument('--labels', default='cell agglomeration,polynomial extension',
                    help='labels of the two directories')
    a = ap.parse_args()
    R = Results(a.dir)
    figs = os.path.join(a.dir, 'figures')
    md = []
    P = md.append
    have2d = os.path.exists(os.path.join(a.dir, 'circle_full'))
    have3d = any(n.startswith('3d_') for n in os.listdir(a.dir))

    if have2d:
        P('# 2D: quadrangles, one chord per cut cell\n')
        if a.proton:
            ref = parse_proton(a.proton)
            P("## Comparison with ProtoN (cuthho_bvc_results.txt)\n")
            P('| problem | variant | k | meshes | ndof identical | max rel. diff. H1 | L2 | energy |')
            P('|---|---|---|---|---|---|---|---|')
            for prob in ['circle', 'flower']:
                for var in ['full', 'naive']:
                    mine = R.t2(f'{prob}_{var}')
                    for k in range(4):
                        r, m = ref[prob][(k, var)], mine[(k, var)]
                        nd = len(r) == len(m) and all(x['ndof'] == y['ndof'] for x, y in zip(r, m))
                        dif = [max(abs(x[q] - y[q]) / x[q] for x, y in zip(r, m)) for q in ['H1', 'L2', 'EN']]
                        P(f"| {prob} | {var} | {k} | {','.join(str(x['N']) for x in r)} | {nd} | "
                          + ' | '.join(f'{x:.1e}' for x in dif) + ' |')
            P('')
        P('## Least-squares rates H1 / L2 / energy (whole range)\n')
        P('| run | k = 0 | k = 1 | k = 2 | k = 3 |')
        P('|---|---|---|---|---|')
        for name, var in [(f'{p}_{v}', v) for p in ['circle', 'flower'] for v in ['naive', 'gam', 'taylor', 'full', 'ext']] + \
                         [('circle_full_P1', 'full'), ('flower_full_P1', 'full'), ('circle_full_tri', 'full'),
                          ('flower_full_tri', 'full'), ('flower_naive_tri', 'naive')]:
            d = R.t2(name)
            if not d:
                continue
            P(f'| {name} | ' + ' | '.join('/'.join(f'{ls_rate(d[(k, var)], q):.2f}' for q in ['H1', 'L2', 'EN'])
                                          for k in range(4)) + ' |')
        P('')
        for prob in ['circle', 'flower']:
            d = R.t2(f'{prob}_full')
            for k in range(4):
                P(f'### {prob}, full BVC, k = {k} (m = {M_OF_K[k]})\n')
                P(fmt_rows(d[(k, 'full')]) + '\n')
        P('## Interface refinement without correction against BVC with one chord (H1, k = 2, 3)\n')
        P('| problem | k | run | errors | LS rate |')
        P('|---|---|---|---|---|')
        for prob in ['circle', 'flower']:
            for k in [2, 3]:
                rows = [('BVC, 1 chord', R.t2(f'{prob}_full')[(k, 'full')]),
                        ('no BVC, 1 chord', R.t2(f'{prob}_naive')[(k, 'naive')])]
                rows += [(f'no BVC, {2**n} chords', R.t2(f'{prob}_naive_n{n}')[(k, 'naive')]) for n in [2, 3, 4]]
                for lab, r in rows:
                    errs = ' '.join('%.2e' % x['H1'] for x in r)
                    P(f"| {prob} | {k} | {lab} | {errs} | {ls_rate(r, 'H1'):.2f} |")
        P('')
        P('## Contrast (circle, h = 1/64, full BVC): H1 / L2 / energy\n')
        P('| k | kappa2 = 1 | 1e2 | 1e4 | 1e6 |')
        P('|---|---|---|---|---|')
        for k in range(4):
            P(f'| {k} | ' + ' | '.join(
                ' / '.join(f"{R.t2(f'contrast_m{m}')[(k, 'full')][-1][q]:.3e}" for q in ['H1', 'L2', 'EN'])
                for m in [0, 2, 4, 6]) + ' |')
        P('')
        P('## Taylor order m (gamma_T on): LS rates H1 / L2\n')
        P('| problem | k | m = 0 | m = 1 | m = 2 |')
        P('|---|---|---|---|---|')
        for prob in ['circle', 'flower']:
            for k in [1, 2, 3]:
                P(f'| {prob} | {k} | ' + ' | '.join(
                    f"{ls_rate(R.t2(f'{prob}_T{T}')[(k, 'gam')], 'H1'):.2f} / {ls_rate(R.t2(f'{prob}_T{T}')[(k, 'gam')], 'L2'):.2f}"
                    for T in [0, 1, 2]) + ' |')
        P('')
        P('## Patch test and quadrature check\n')
        for name in ['patch2d', 'patch2d_tri']:
            for (k, var), r in sorted(R.t2(name).items()):
                P(f"- {name}, k = {k}: max H1 {max(x['H1'] for x in r):.1e}, L2 {max(x['L2'] for x in r):.1e}, "
                  f"energy {max(x['EN'] for x in r):.1e}")
        for prob in ['circle', 'flower']:
            A, B = R.t2(f'{prob}_full'), R.t2(f'{prob}_full_Q4')
            if B:
                dif = {q: max(abs(x[q] - y[q]) / x[q] for k in range(4) for x, y in zip(A[(k, 'full')], B[(k, 'full')]))
                       for q in ['H1', 'L2', 'EN']}
                P(f"- {prob}, all quadrature degrees + 8: max relative change H1 {dif['H1']:.1e}, "
                  f"L2 {dif['L2']:.1e}, energy {dif['EN']:.1e}")
        P('')
        P('## Geometry (Table 1 of the paper)\n')
        for g in ['geo_circle', 'geo_flower']:
            p = os.path.join(a.dir, g, 'bvc_geometry.txt')
            if os.path.exists(p):
                P(f'{g}:\n```\n' + open(p).read() + '```\n')

        # figures 2D
        fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.4))
        for j, (prob, title) in enumerate([('circle', r'circle, $\kappa_2/\kappa_1 = 10^4$'),
                                           ('flower', 'flower, discontinuous solution')]):
            d = R.t2(f'{prob}_full')
            for i, q in enumerate(['H1', 'L2']):
                for k in range(4):
                    series(axs[i, j], d[(k, 'full')], q, k, f'k = {k} (m = {M_OF_K[k]})')
                    slope(axs[i, j], d[(k, 'full')], q, k + QOPT[q])
                finish(axs[i, j], f'2D {title}', q)
        axs[0, 0].legend(loc='lower left')
        save(fig, figs, 'bvc2d_full')

        VAR = [('naive', 'no correction'), ('gam', r'$\gamma_T$ only'), ('taylor', 'Taylor only'), ('full', 'full BVC')]
        fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.5))
        for ax, (prob, k) in zip(axs, [('circle', 2), ('flower', 2), ('flower', 3)]):
            for i, (var, lab) in enumerate(VAR):
                series(ax, R.t2(f'{prob}_{var}')[(k, var)], 'H1', i, lab)
            slope(ax, R.t2(f'{prob}_full')[(k, 'full')], 'H1', k + 1)
            finish(ax, f'2D {prob}, k = {k}', 'H1')
        axs[0].legend(loc='lower left')
        save(fig, figs, 'bvc2d_corrections')

        fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.6))
        for ax, prob in zip(axs, ['flower', 'circle']):
            series(ax, R.t2(f'{prob}_full')[(3, 'full')], 'H1', 0, 'BVC, 1 chord')
            series(ax, R.t2(f'{prob}_naive')[(3, 'naive')], 'H1', 1, 'no BVC, 1 chord')
            for i, n in enumerate([2, 3, 4]):
                series(ax, R.t2(f'{prob}_naive_n{n}')[(3, 'naive')], 'H1', i + 2, f'no BVC, {2**n} chords')
            slope(ax, R.t2(f'{prob}_full')[(3, 'full')], 'H1', 4)
            finish(ax, f'2D {prob}, k = 3', 'H1')
        axs[0].legend(loc='lower left')
        save(fig, figs, 'bvc2d_chords')

        fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.4))
        for ax, q in zip(axs, ['H1', 'EN']):
            for k in range(4):
                es = [R.t2(f'contrast_m{m}')[(k, 'full')][-1][q] for m in [0, 2, 4, 6]]
                ax.loglog([1, 1e2, 1e4, 1e6], es, color=SLOTS[k], marker=MARKERS[k], label=f'k = {k}',
                          markeredgecolor=SURF, markeredgewidth=0.6)
            ax.set_xlabel(r'$\kappa_2/\kappa_1$')
            ax.set_ylabel(QLABEL[q])
            ax.set_title('2D circle, $h = 1/64$, full BVC', loc='left')
        axs[0].legend(loc='upper right', ncol=2)
        save(fig, figs, 'bvc2d_contrast')

    if have3d:
        P('# 3D: hexahedra (M3) and tetrahedra (M4), flat interface\n')
        MESH = {3: 'hexahedra', 4: 'tetrahedra'}
        PROB = {'sphere': r'sphere, $\kappa_2/\kappa_1 = 10^4$', 'flower': 'flower, discontinuous solution'}
        P('## Least-squares rates H1 / L2 / energy (all the meshes; [finest three])\n')
        P('| problem | mesh | variant | k = 0 | k = 1 | k = 2 | k = 3 |')
        P('|---|---|---|---|---|---|---|')
        for prob in ['sphere', 'flower']:
            for M in [3, 4]:
                for var in ['full', 'naive', 'taylor', 'gam']:
                    cells = []
                    for k in range(4):
                        r = R.t3(prob, M, var, k)
                        if not r:
                            cells.append('')
                            continue
                        cells.append('/'.join(f'{ls_rate(r, q):.2f}' for q in ['H1', 'L2', 'EN']) + ' [' +
                                     '/'.join(f'{ls_rate(r[-3:], q):.2f}' for q in ['H1', 'L2', 'EN']) + ']')
                    if any(cells):
                        P(f'| {prob} | {MESH[M]} | {var} | ' + ' | '.join(cells) + ' |')
        P('')
        for prob in ['sphere', 'flower']:
            for M in [3, 4]:
                for var in ['full', 'naive']:
                    for k in range(4):
                        r = R.t3(prob, M, var, k)
                        if r:
                            P(f'### {prob}, {MESH[M]}, {var}, k = {k}\n')
                            P(fmt_rows(r) + '\n')
        P('## Contrast (sphere, full BVC): H1 / L2 / energy\n')
        P('| mesh | k | kappa2 = 1 | 1e2 | 1e4 | 1e6 |')
        P('|---|---|---|---|---|---|')
        for M in [3, 4]:
            for k in range(4):
                cells = []
                for m in [0, 2, 4, 6]:
                    t = parse_table(os.path.join(a.dir, f'3d_contrast_M{M}_m{m}', 'bvc_convergence.txt'))
                    cells.append(' / '.join(f"{t[(k, 'full')][-1][q]:.3e}" for q in ['H1', 'L2', 'EN']) if t else '')
                P(f'| {MESH[M]} | {k} | ' + ' | '.join(cells) + ' |')
        P('')
        P('## Patch test, quadrature, geometry\n')
        for M in [3, 4]:
            t = parse_table(os.path.join(a.dir, f'3d_patch_M{M}', 'bvc_convergence.txt'))
            for (k, var), r in sorted(t.items()):
                P(f"- patch test, {MESH[M]}, k = {k}: max H1 {max(x['H1'] for x in r):.1e}, "
                  f"L2 {max(x['L2'] for x in r):.1e}, energy {max(x['EN'] for x in r):.1e}")
            A = parse_table(os.path.join(a.dir, f'3d_quad0_sphere_M{M}', 'bvc_convergence.txt'))
            B = parse_table(os.path.join(a.dir, f'3d_quad_sphere_M{M}', 'bvc_convergence.txt'))
            if A and B:
                dif = {q: max(abs(x[q] - y[q]) / x[q] for k in range(4) for x, y in zip(A[(k, 'full')], B[(k, 'full')]))
                       for q in ['H1', 'L2', 'EN']}
                P(f"- sphere, {MESH[M]}, quadrature degrees + 4: max relative change H1 {dif['H1']:.1e}, "
                  f"L2 {dif['L2']:.1e}, energy {dif['EN']:.1e}")
        for g in ['sphere_M3', 'sphere_M4', 'flower_M3', 'flower_M4']:
            p = os.path.join(a.dir, f'3d_geo_{g}', 'bvc_geometry.txt')
            if os.path.exists(p):
                P(f'\n{g}:\n```\n' + open(p).read() + '```')
        P('')

        # figures 3D: full method, H1 and L2, the four configurations
        for q in ['H1', 'L2']:
            fig, axs = plt.subplots(1, 4, figsize=(14.0, 3.6))
            for ax, (prob, M) in zip(axs, [('sphere', 3), ('sphere', 4), ('flower', 3), ('flower', 4)]):
                for k in range(4):
                    r = R.t3(prob, M, 'full', k)
                    if r:
                        series(ax, r, q, k, f'k = {k} (m = {M_OF_K[k]})')
                        slope(ax, r, q, k + QOPT[q])
                finish(ax, f'3D {prob}, {MESH[M]}', q)
            axs[0].legend(loc='lower left')
            save(fig, figs, f'bvc3d_full_{q}')
        # full against no correction
        fig, axs = plt.subplots(1, 4, figsize=(14.0, 3.6))
        for ax, (prob, M) in zip(axs, [('sphere', 3), ('sphere', 4), ('flower', 3), ('flower', 4)]):
            q = 'L2' if prob == 'sphere' else 'H1'
            for i, k in enumerate([1, 2, 3]):
                rf, rn = R.t3(prob, M, 'full', k), R.t3(prob, M, 'naive', k)
                if rf:
                    series(ax, rf, q, k, f'k = {k}, full BVC')
                if rn:
                    series(ax, rn, q, k, f'k = {k}, no correction', dashed=True)
            finish(ax, f'3D {prob}, {MESH[M]}', q)
        axs[0].legend(loc='lower left', fontsize=7)
        save(fig, figs, 'bvc3d_full_vs_naive')
        # corrections one at a time (flower, hexahedra)
        VAR = [('naive', 'no correction'), ('gam', r'$\gamma_T$ only'), ('taylor', 'Taylor only'), ('full', 'full BVC')]
        fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.6))
        for ax, k in zip(axs, [2, 3]):
            for i, (var, lab) in enumerate(VAR):
                r = R.t3('flower', 3, var, k)
                if r:
                    series(ax, r, 'H1', i, lab)
            r = R.t3('flower', 3, 'full', k)
            if r:
                slope(ax, r, 'H1', k + 1)
            finish(ax, f'3D flower, hexahedra, k = {k}', 'H1')
        axs[0].legend(loc='lower left')
        save(fig, figs, 'bvc3d_corrections')

    if a.compare:
        R2 = Results(a.compare)
        lab = a.labels.split(',')
        P(f'# Full BVC: {lab[0]} against {lab[1]}\n')
        P('Least-squares rates H1 / L2 / energy (2D: all the meshes; 3D: finest three) and H1 error on the finest mesh.\n')
        P(f'| case | k | {lab[0]} | {lab[1]} |')
        P('|---|---|---|---|')
        def cell(rows, three):
            if not rows:
                return ''
            r = rows[-3:] if three else rows
            return '/'.join(f'{ls_rate(r, q):.2f}' for q in ['H1', 'L2', 'EN']) + f" ({rows[-1]['H1']:.2e})"
        for prob in ['circle', 'flower']:
            for k in range(4):
                c1 = R.t2(f'{prob}_full').get((k, 'full'), [])
                c2 = R2.t2(f'{prob}_full').get((k, 'full'), [])
                if c1 or c2:
                    P(f'| 2D {prob} | {k} | {cell(c1, False)} | {cell(c2, False)} |')
        for prob in ['sphere', 'flower']:
            for M, mesh in [(3, 'hexahedra'), (4, 'tetrahedra')]:
                for k in range(4):
                    c1, c2 = R.t3(prob, M, 'full', k), R2.t3(prob, M, 'full', k)
                    if c1 or c2:
                        P(f'| 3D {prob}, {mesh} | {k} | {cell(c1, True)} | {cell(c2, True)} |')
        P('')
        def compare_panel(ax, r1, r2, title):
            for k in range(4):
                if r1.get(k):
                    series(ax, r1[k], 'H1', k, f'k = {k}, {lab[0]}')
                if r2.get(k):
                    series(ax, r2[k], 'H1', k, f'k = {k}, {lab[1]}', dashed=True)
            finish(ax, title, 'H1')
        fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.8))
        for ax, prob in zip(axs, ['circle', 'flower']):
            r1 = {k: R.t2(f'{prob}_full').get((k, 'full'), []) for k in range(4)}
            r2 = {k: R2.t2(f'{prob}_full').get((k, 'full'), []) for k in range(4)}
            compare_panel(ax, r1, r2, f'2D {prob}, full BVC')
        axs[0].legend(loc='lower left', fontsize=6.5)
        save(fig, figs, 'bvc2d_compare')
        if have3d:
            fig, axs = plt.subplots(1, 4, figsize=(14.0, 3.8))
            for ax, (prob, M, mesh) in zip(axs, [('sphere', 3, 'hexahedra'), ('sphere', 4, 'tetrahedra'),
                                                 ('flower', 3, 'hexahedra'), ('flower', 4, 'tetrahedra')]):
                r1 = {k: R.t3(prob, M, 'full', k) for k in range(4)}
                r2 = {k: R2.t3(prob, M, 'full', k) for k in range(4)}
                compare_panel(ax, r1, r2, f'3D {prob}, {mesh}')
            axs[0].legend(loc='lower left', fontsize=6.5)
            save(fig, figs, 'bvc3d_compare')

    open(os.path.join(a.dir, 'bvc_tables.md'), 'w').write('\n'.join(md) + '\n')
    print('written', os.path.join(a.dir, 'bvc_tables.md'), 'and', figs)


if __name__ == '__main__':
    main()
