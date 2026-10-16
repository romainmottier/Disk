#!/usr/bin/env python3
"""Numerical results of the paper on the unfitted HHO method with boundary
value correction (prototypes -p2: cell agglomeration, -p3: polynomial
extension) in 2D and 3D.

    python3 run_bvc_paper.py <path/to/unfitted_HHO> <outdir> [--run] [--plot]
                             [--paper DIR] [--mem GB] [--threads N]
                             [--only group,...] [--list]

--run runs the experiments (finished runs are skipped when restarted):
first all the groups but "efficiency" within the memory and thread budgets,
then the group "efficiency" one run at a time on --threads threads (wall
times). --plot writes the figures <outdir>/paper/figures/*.pdf (curves and
bar charts) and <outdir>/paper/report_figures/*.pdf (details of the
validation), and --paper DIR copies the former into DIR/figures.

Groups:
  verif        patch tests (quadrangles, triangles, hexahedra, tetrahedra),
               geometric checks (-t3), geometry of the cut meshes (-f)
  conv2d       circle (kappa_2/kappa_1 = 1e4) and flower (discontinuous
               solution), quadrangles and triangles, k = 0..3
  corrections  the two corrections one at a time (2D circle and flower, 3D
               flower on hexahedra), interface refinement without correction
  contrast     robustness with respect to kappa_2/kappa_1 (2D, 3D)
  taylor       order m of the Taylor correction
  stab         weight of the stabilization s°
  conv3d       sphere (kappa_2/kappa_1 = 1e4) and 3D flower (discontinuous
               solution), hexahedra and tetrahedra, k = 0..3
  geometries   other interfaces (ellipse/ellipsoid, rounded square/cube,
               annulus/torus) with the jumps (59), local refinement -x
  efficiency   error against wall time: BVC against the quadratures on the
               curved geometry (-p0, -p1)

BVC settings: weight (k+1) kappa_i/h_T of s° (-A -1), GMRES preconditioned by
the Cholesky factorization of the symmetric part (-S 2), default pairing,
threshold 0.3, eta = 1, eta_N = 20 (k+1)^2, Taylor order m(k).
"""
import argparse, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_bvc import run as run_jobs, memory  # noqa: E402

BVC = '-A -1 -S 2'


def fbs(k):
    return (k + 1) * (k + 2) // 2


def mem3d(M, k, N, p):
    """same memory for both stabilizations (measured)"""
    return memory(M, k, N)


def jobs():
    """the runs of the paper and of the validation; every convergence
    curve has at least five meshes, the same for all the degrees"""
    J = []   # (group, name, args, mem, threads)

    def add(group, name, args, mem=1.5, thr=2):
        J.append((group, name, args, mem, thr if mem <= 6 else 12))

    # validation: the 2D test cases of [Bur_BVC] with its settings (cell
    # agglomeration with the pairing of ProtoN, s° weighted by kappa_i/h_T, LU)
    B = '-p2 -M0 -x0 -P0 -A 1 -S 0'
    BC = f'{B} -i0 -t0 -m4 -N8 -k3 -l4'      # N = 8 ... 128
    BF = f'{B} -i7 -t1 -m0 -N16 -k3 -l4'     # N = 16 ... 256
    for prob, a in [('circle', BC), ('flower', BF)]:
        for var in ['full', 'naive', 'gam', 'taylor']:
            add('burman', f'burman_{prob}_{var}', f'{a} -K0 -v {var}')
        for n in [2, 3, 4]:
            add('burman', f'burman_{prob}_naive_n{n}', f'{a} -K2 -v naive -n{n}')
        for T in [0, 1, 2]:
            add('burman', f'burman_{prob}_T{T}', f'{a} -K1 -v gam -T{T}')
    for m in [0, 2, 4, 6]:
        add('burman', f'burman_contrast_m{m}', f'{B} -i0 -t0 -m{m} -N8 -K0 -k3 -L3 -l3 -v full')
    for name, a in [('circle', '-i0'), ('flower', '-i7')]:
        add('burman', f'burman_geo_{name}', f'{B} {a} -t3 -N8 -K0 -k1 -l4')

    # verification (report) and meshes (paper)
    for M, mesh in [(0, 'quad'), (1, 'tri'), (3, 'hex'), (4, 'tet')]:
        for p in [2, 3]:
            add('verif', f'verif_patch_{mesh}_p{p}', f'-p{p} -M{M} -x0 -i8 -t2 -N4 -K0 -k3 -l2 {BVC}',
                mem3d(M, 3, 16, p))
    for name, a in [('sphere_hex', '-M3 -i0 -N4 -l3'), ('sphere_tet', '-M4 -i0 -N4 -l3'),
                    ('flower_hex', '-M3 -i7 -N8 -l2'), ('flower_tet', '-M4 -i7 -N8 -l2')]:
        add('verif', f'verif_geo_{name}', f'-p3 -x0 -t3 -K0 -k1 {a}', 2)
    for name, a in [('circle_p2', '-p2 -M0 -i0 -t0 -m4'), ('circle_p3', '-p3 -M0 -i0 -t0 -m4'),
                    ('flower_p2', '-p2 -M0 -i7 -t1'), ('flower_p3', '-p3 -M0 -i7 -t1')]:
        add('verif', f'verif_mesh_{name}', f'{a} -x0 -N16 -K1 -k1 -l0 -f1 {BVC}')
    for name, a in [('sphere', '-i0 -t0 -m4 -N8'), ('flower', '-i7 -t1 -N8')]:
        add('verif', f'verif_mesh3d_{name}', f'-p2 -M3 -x0 {a} -K1 -k1 -l0 -f1 {BVC}')
    # thumbnails of the interfaces (flat interface of fine meshes)
    for i in [0, 7, 4, 5, 6]:
        t = '-t0 -m4' if i == 0 else ('-t1' if i == 7 else '-t4')
        add('verif', f'thumb2d_i{i}', f'-p2 -M0 -x3 -i{i} {t} -N128 -K0 -k0 -l0 -f1 {BVC}')
        add('verif', f'thumb3d_i{i}', f'-p2 -M3 -x3 -i{i} {t} -N32 -K0 -k0 -l0 -f1 {BVC}', 2)

    # 2D, proposed settings: extension against agglomeration
    C = '-M0 -x0 -i0 -t0 -m4 -N8 -K0 -k3 -l4'    # N = 8 ... 128
    F = '-M0 -x0 -i7 -t1 -N16 -K0 -k3 -l4'       # N = 16 ... 256
    for p in [2, 3]:
        add('conv2d', f'conv2d_circle_p{p}', f'-p{p} {C} {BVC}')
        add('conv2d', f'conv2d_flower_p{p}', f'-p{p} {F} {BVC}')

    # contrast: kappa_2/kappa_1 = 1 ... 10^8 (the method assumes kappa_1 <= kappa_2)
    for m in range(9):
        add('contrast', f'contrast2d_m{m}', f'-p3 -M0 -x0 -i0 -t0 -m{m} -N8 -K0 -k3 -L3 -l3 {BVC}')
        add('contrast', f'contrast3d_m{m}', f'-p3 -M3 -x0 -i0 -t0 -m{m} -N16 -K0 -k3 -l0 {BVC}',
            mem3d(3, 3, 16, 3))

    # weight of s°: kappa_i/h_T of [Bur_BVC] against (k+1) kappa_i/h_T
    for a, tag in [('1', 'A1'), ('-1', 'Ak')]:
        b = BVC.replace('-A -1', f'-A {a}')
        add('stab', f'stab_flower_{tag}', f'-p3 {F} {b}')
        for k in [1, 2, 3]:   # sphere, hexahedra: N = 4, 5, 8, 10, 16, 20
            add('stab', f'stab_sphere_hex_{tag}_k{k}', f'-p3 -M3 -x0 -i0 -t0 -m4 -N4 -K{k} -k{k} -l2 {b}',
                mem3d(3, k, 16, 3))
            add('stab', f'stab_sphere_hex_{tag}_k{k}_N5', f'-p3 -M3 -x0 -i0 -t0 -m4 -N5 -K{k} -k{k} -l2 {b}',
                mem3d(3, k, 20, 3))

    # 3D convergence: two families of meshes, merged (same meshes for all the
    # degrees): sphere, hexahedra N = 4, 5, 8, 10, 16, 20, 32, tetrahedra
    # N = 4, 5, 8, 10, 16, 20; flower, hexahedra N = 6, 8, 12, 16, 24, 32,
    # tetrahedra N = 6, 8, 12, 16, 24
    S = '-i0 -t0 -m4'
    Fl = '-i7 -t1 -m0'
    for p in [2, 3]:
        for k in range(4):
            for prob, M, N0, L in [('sphere', 3, 4, 2 if k == 3 else 3), ('sphere', 3, 5, 2),
                                   ('sphere', 4, 4, 3 if k <= 1 else 2), ('sphere', 4, 5, 2 if k <= 1 else 1),
                                   ('flower', 3, 8, 1 if k == 3 else 2), ('flower', 3, 6, 2),
                                   ('flower', 4, 8, 2 if k <= 1 else 1), ('flower', 4, 6, 2 if k <= 1 else 1)]:
                a = S if prob == 'sphere' else Fl
                add('conv3d', f'conv3d_{prob}_M{M}_p{p}_full_k{k}_N{N0}',
                    f'-p{p} -M{M} -x0 {a} -N{N0} -K{k} -k{k} -l{L} -v full {BVC}', mem3d(M, k, N0 * 2**L, p))
            if k == 3:   # N = 32 on hexahedra (39 GB)
                add('conv3d', f'conv3d_sphere_M3_p{p}_full_k3_N32',
                    f'-p{p} -M3 -x0 {S} -N32 -K3 -k3 -l0 -v full {BVC}', mem3d(3, 3, 32, p))
                add('conv3d', f'conv3d_flower_M3_p{p}_full_k3_N32',
                    f'-p{p} -M3 -x0 {Fl} -N32 -K3 -k3 -l0 -v full {BVC}', mem3d(3, 3, 32, p))
            if k >= 2:   # N = 20 (sphere) and 24 (flower) on tetrahedra (up to 46 GB)
                add('conv3d', f'conv3d_sphere_M4_p{p}_full_k{k}_N20',
                    f'-p{p} -M4 -x0 {S} -N20 -K{k} -k{k} -l0 -v full {BVC}', mem3d(4, k, 20, p))
                add('conv3d', f'conv3d_flower_M4_p{p}_full_k{k}_N24',
                    f'-p{p} -M4 -x0 {Fl} -N24 -K{k} -k{k} -l0 -v full {BVC}', mem3d(4, k, 24, p))

    # other interfaces, jumps (59), local refinement of the cells not cut once:
    # 2D N = 16 ... 256, 3D N = 8, 12, 16, 24, 32
    for i in [4, 5, 6]:
        add('geometries', f'geom2d_i{i}', f'-p3 -M0 -x3 -i{i} -t4 -N16 -K1 -k3 -l4 {BVC}')
        for k in [1, 2, 3]:
            L = 2 if k <= 2 else 1
            add('geometries', f'geom3d_i{i}_k{k}', f'-p3 -M3 -x3 -i{i} -t4 -N8 -K{k} -k{k} -l{L} {BVC}',
                mem3d(3, k, 8 * 2**L, 3))
            add('geometries', f'geom3d_i{i}_k{k}_N12', f'-p3 -M3 -x3 -i{i} -t4 -N12 -K{k} -k{k} -l1 {BVC}',
                mem3d(3, k, 24, 3))
        add('geometries', f'geom3d_i{i}_k3_N32', f'-p3 -M3 -x3 -i{i} -t4 -N32 -K3 -k3 -l0 {BVC}', mem3d(3, 3, 32, 3))
    return J


EFF_CASES = {   # name: (mesh, problem, k, [(N0, levels), ...] of the two families of meshes)
    '2d_flower':     (0, 'jumps', 3, [(16, 4)]),
    '3d_sphere_hex': (3, 'contrast', 2, [(4, 3), (5, 2)]),
    '3d_sphere_tet': (4, 'contrast', 3, [(4, 2), (5, 2)]),
    '3d_flower_hex': (3, 'jumps', 2, [(8, 2), (6, 2)]),
}
EFF_METHODS = ['p3', 'p0', 'g0r8', 'g0r4']   # BVC, exact curved interface, refined interface r = 8, 4


def eff_name(case, m, k, N0, first):
    return f'eff_{case}_{m}_k{k}' if first else f'eff_{case}_{m}_k{k}_N{N0}'


def eff_jobs():
    """efficiency: one run at a time; BVC (extension, proposed settings), and
    the quadratures on the exact curved interface (-p0 -g -1) or on a refined
    interface (-p0 -g 0 -r r), polynomial extension"""
    J = []
    for case, (M, prob, k, fams) in EFF_CASES.items():
        t = '-t1 -m4 -i0' if prob == 'contrast' else '-t4 -i7'
        for j, (N0, L) in enumerate(fams):
            for m in EFF_METHODS:
                if m == 'p3':
                    a = f'-p3 -M{M} -x0 {t.replace("-t1", "-t0")} -N{N0} -K{k} -k{k} -l{L} {BVC}'
                elif m == 'p0':
                    a = f'-p0 -M{M} -x0 {t} -N{N0} -K{k} -k{k} -l{L}'
                else:
                    # the refined interface cannot be built on the mesh N = 8 of the 3D flower
                    lmin = ' -L1' if (case == '3d_flower_hex' and N0 == 8) else ''
                    a = f'-p0 -M{M} -x0 {t} -N{N0} -K{k} -k{k}{lmin} -l{L} -g 0 -r {m[3:]}'
                J.append(('efficiency', eff_name(case, m, k, N0, j == 0), a, mem3d(M, k, N0 * 2**L, 3), 0))
    return J


# --------------------------------------------------------------------------- figures (style)
# Style of the figures of the LTS-HHO preprint: Computer Modern fonts, dotted
# grid, observed order on each segment of the convergence curves, dashed
# reference slope, annotation box in the upper left corner, legend in the
# lower right corner. No tables: every result is a curve or a bar chart.
RED, BLUE, BLACK, GREEN, GRAY = '#ac1117', '#1764ab', '#000000', '#238b45', '#666666'


def paper_style(plt):
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
        'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False,
        'font.size': 8, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7.5,
        'axes.titlesize': 8.5,
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'axes.edgecolor': 'black', 'axes.linewidth': 0.8, 'axes.labelcolor': 'black',
        'xtick.color': 'black', 'ytick.color': 'black', 'text.color': 'black',
        'axes.grid': True, 'axes.grid.which': 'major', 'grid.linestyle': ':', 'grid.color': '#999999',
        'grid.linewidth': 0.7, 'axes.axisbelow': True,
        'lines.linewidth': 1.4, 'lines.markersize': 4.5,
        'legend.frameon': True, 'legend.fancybox': False, 'legend.edgecolor': '#999999', 'legend.framealpha': 1.0,
        'legend.borderpad': 0.4, 'legend.handlelength': 2.2,
        'axes.titlelocation': 'left',
    })


def family(N):
    """meshes N0 2^l of one family share the odd part of N"""
    while N % 2 == 0:
        N //= 2
    return N


def split_families(rows):
    """main family (the one with the finest mesh, then the most meshes) and the other points"""
    fams = {}
    for r in rows:
        fams.setdefault(family(r['N']), []).append(r)
    if not fams:
        return [], []
    main = max(fams, key=lambda f: (max(r['N'] for r in fams[f]), len(fams[f])))
    other = [r for f, rr in fams.items() if f != main for r in rr]
    return sorted(fams[main], key=lambda r: r['h']), sorted(other, key=lambda r: r['h'])


def curve(ax, rows, q, color, marker, ls='-', label=None, orders=None, open_marker=False, families=False):
    """error q against h; orders: 'above' / 'below' / None writes the observed
    order on each segment; families: only the main family of meshes is joined
    (the other meshes are drawn as open markers)"""
    rows = sorted([r for r in rows if r.get(q, 0) > 0], key=lambda r: r['h'])
    if not rows:
        return []
    main, other = split_families(rows) if families else (rows, [])
    h = [r['h'] for r in main]
    e = [r[q] for r in main]
    ax.loglog(h, e, color=color, marker=marker, linestyle=ls, label=label, zorder=3,
              markerfacecolor='white' if open_marker else color, markeredgecolor=color, markeredgewidth=1.0)
    if other:
        ax.loglog([r['h'] for r in other], [r[q] for r in other], color=color, marker=marker, linestyle='none',
                  zorder=3, markerfacecolor='white', markeredgecolor=color, markeredgewidth=1.0)
    if orders:
        up = orders == 'above'
        for i in range(len(h) - 1):
            p = math.log(e[i + 1] / e[i]) / math.log(h[i + 1] / h[i])
            ax.annotate(f'{p:.2f}', (math.sqrt(h[i] * h[i + 1]), math.sqrt(e[i] * e[i + 1])),
                        xytext=(-5, 4) if up else (5, -4), textcoords='offset points',
                        ha='right' if up else 'left', va='bottom' if up else 'top', color=color, fontsize=6.5,
                        zorder=4)
    return rows


def ref_slope(ax, rows_list, q, p, factor=0.3):
    """dashed reference line of slope p below the curves"""
    pts = [(r['h'], r[q]) for rows in rows_list for r in rows if r.get(q, 0) > 0]
    if not pts:
        return
    c = min(e / h**p for h, e in pts) * factor
    lo, hi = min(h for h, _ in pts), max(h for h, _ in pts)
    hs = [lo / 1.15, hi * 1.15]
    ax.loglog(hs, [c * x**p for x in hs], linestyle='--', color=GRAY, linewidth=1.0, zorder=2)


def klabel(ax, rows_list, q, text):
    """label of a group of curves (e.g. k = 2), left of their finest points"""
    fin = [min(rows, key=lambda r: r['h']) for rows in rows_list if rows]
    fin = [r for r in fin if r.get(q, 0) > 0]
    if not fin:
        return
    h = min(r['h'] for r in fin)
    e = math.exp(sum(math.log(r[q]) for r in fin) / len(fin))
    ax.annotate(text, (h, e), xytext=(-7, 0), textcoords='offset points', ha='right', va='center', fontsize=7.5)


def tag(ax, text, loc='upper left'):
    x, ha = (0.04, 'left') if 'left' in loc else (0.96, 'right')
    y, va = (0.95, 'top') if 'upper' in loc else (0.05, 'bottom')
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7.5, zorder=6,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#999999', linewidth=0.8))


def frame_h(ax, xlabel='$h$', left=2.0):
    """x limits with room for the labels of the curves on the left"""
    data = [l for l in ax.get_lines() if l.get_linestyle() != '--' or l.get_marker() not in ('None', None, '')]
    xs = [x for l in data for x in l.get_xdata() if x > 0]
    ys = [y for l in data for y in l.get_ydata() if y > 0]
    if xs:
        ax.set_xlim(min(xs) / left, max(xs) * 1.3)
    if ys:
        lo, hi = min(ys), max(ys)
        r = max(hi / lo, 10.0)
        ax.set_ylim(lo / r**0.06, hi * r**0.18)     # room for the annotation box
    ax.set_xlabel(xlabel)


def log_minor_labels(fig):
    """labels 2 10^e and 5 10^e on the log axes spanning less than 1.5 decades
    (the default labels use a glyph missing from the Computer Modern text font)"""
    from matplotlib.ticker import FuncFormatter

    def make(axis_lim, ms):
        def fmt(x, pos):
            lo, hi = sorted(axis_lim())
            if lo <= 0 or math.log10(hi / lo) > 1.5:
                return ''
            e = math.floor(math.log10(x) + 1e-9)
            m = round(x / 10**e)
            return rf'${m}{{\cdot}}10^{{{e}}}$' if m in ms else ''
        return fmt
    for ax in fig.axes:
        if getattr(ax, 'name', '') == '3d':
            continue
        # narrow panels: 5 10^e only
        ms = (2, 5) if ax.get_position().width * fig.get_figwidth() > 2.6 else (5,)
        if ax.get_xscale() == 'log':
            ax.xaxis.set_minor_formatter(FuncFormatter(make(ax.get_xlim, ms)))
        if ax.get_yscale() == 'log':
            ax.yaxis.set_minor_formatter(FuncFormatter(make(ax.get_ylim, ms)))


QLAB = {'H1': r'$|u-u_h|_{H^1}$', 'L2': r'$\|u-u_h\|_{L^2}$', 'EN': 'energy error'}
QOPT = {'H1': 1, 'L2': 2, 'EN': 1}


# --------------------------------------------------------------------------- 3D views of the cut meshes
class OrthoView:
    """orthographic projection and painter's algorithm (vector output)"""

    def __init__(self, azim=-60.0, elev=28.0):
        import numpy as np
        a, e = math.radians(azim), math.radians(elev)
        self.d = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])  # to the viewer
        self.r = np.array([-math.sin(a), math.cos(a), 0.0])
        self.u = np.cross(self.d, self.r)
        l = -0.45 * self.r + 0.55 * self.u + 0.70 * self.d
        self.light = l / np.linalg.norm(l)
        self.prims = []   # (depth, polygon 2D, facecolor, edgecolor, linewidth)

    def proj(self, P):
        import numpy as np
        P = np.atleast_2d(np.asarray(P, dtype=float))
        return np.stack([P @ self.r, P @ self.u], axis=1), P @ self.d

    def segment(self, p, q, color, lw, subdiv=6, bias=0.0):
        import numpy as np
        p, q = np.asarray(p, dtype=float), np.asarray(q, dtype=float)
        for i in range(subdiv):
            xy, z = self.proj([p + (q - p) * (i / subdiv), p + (q - p) * ((i + 1) / subdiv)])
            self.prims.append((z.mean() + bias, xy, 'none', color, lw))

    def polygon(self, P, base, edge='none', lw=0.0, ambient=0.5, alpha=1.0, bias=0.0, flat=False):
        import numpy as np
        P = np.asarray(P, dtype=float)
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n /= max(np.linalg.norm(n), 1e-300)
        c = np.asarray(base, dtype=float)
        if not flat:
            c = c * (ambient + (1 - ambient) * abs(n @ self.light))
        xy, z = self.proj(P)
        self.prims.append((z.mean() + bias, xy, tuple(np.clip(c, 0, 1)) + (alpha,), edge, lw))

    def draw(self, ax, margin=0.02):
        import numpy as np
        from matplotlib.collections import PolyCollection
        prims = sorted(self.prims, key=lambda p: p[0])
        ax.add_collection(PolyCollection([p[1] for p in prims], facecolors=[p[2] for p in prims],
                                         edgecolors=[p[3] for p in prims], linewidths=[p[4] for p in prims],
                                         closed=True, joinstyle='round', capstyle='round'))
        xy = np.concatenate([p[1] for p in prims])
        ax.set_xlim(xy[:, 0].min() - margin, xy[:, 0].max() + margin)
        ax.set_ylim(xy[:, 1].min() - margin, xy[:, 1].max() + margin)
        ax.set_aspect('equal')
        ax.set_axis_off()


ITF_BLUE = (0x8f / 255, 0xb8 / 255, 0xe6 / 255)


def read_cut_mesh(cells_file, tri_file):
    import numpy as np
    cells = {}
    for line in open(cells_file):
        t = line.split()
        if t[0] == 'cell':
            cells[int(t[1])] = (np.array(list(map(float, t[2:5]))), np.array(list(map(float, t[5:8]))),
                                int(t[8]), int(t[9]))
    return cells, np.loadtxt(tri_file)


def box_corners(lo, hi):
    import numpy as np
    return [np.array([hi[0] if i & 1 else lo[0], hi[1] if i & 2 else lo[1], hi[2] if i & 4 else lo[2]])
            for i in range(8)]


BOX_EDGES = [(i, i | b) for i in range(8) for b in (1, 2, 4) if (i | b) != i]
BOX_FACES = [(0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4), (2, 3, 7, 6), (0, 2, 6, 4), (1, 3, 7, 5)]


def render_cut_mesh(ax, cells_file, tri_file, N, azim=-60.0, elev=28.0):
    """the background mesh on the three back walls of the domain, the edges of
    the cut cells, and the flat interface Gamma_h (shaded, with its triangles)"""
    import numpy as np
    V = OrthoView(azim, elev)
    cells, tris = read_cut_mesh(cells_file, tri_file)
    # back walls: the faces of (0,1)^3 whose outward normal points away from the viewer
    for axis in range(3):
        for side in (0.0, 1.0):
            normal = np.zeros(3)
            normal[axis] = 1.0 if side == 1.0 else -1.0
            if normal @ V.d >= 0:
                continue
            o1, o2 = [a for a in range(3) if a != axis]
            for i in range(N):
                for j in range(N):
                    P = []
                    for (a_, b_) in [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]:
                        x = np.zeros(3)
                        x[axis], x[o1], x[o2] = side, a_ / N, b_ / N
                        P.append(x)
                    V.polygon(P, (0.965, 0.965, 0.965), edge='#c8c8c8', lw=0.3, flat=True, bias=-10)
    # edges of the cut cells
    edges = set()
    for lo, hi, loc, ct in cells.values():
        if loc != 2:
            continue
        c = box_corners(lo, hi)
        for i, j in BOX_EDGES:
            edges.add((tuple(np.round(c[i], 9)), tuple(np.round(c[j], 9))))
    for p, q in edges:
        V.segment(p, q, '#a8a8a8', 0.3)
    # flat interface
    for row in tris:
        V.polygon(row[:9].reshape(3, 3), ITF_BLUE, edge='#2f5f94', lw=0.15)
    V.draw(ax)


def render_cut_cell(ax, cells_file, tri_file, cell=None, azim=None, elev=28.0):
    """one cut cell: its faces, the chords of its cut faces (red), which form a
    closed loop, and the triangles of Gamma_h in the cell"""
    import numpy as np
    V = OrthoView(azim if azim is not None else -60.0, elev)
    cells, tris = read_cut_mesh(cells_file, tri_file)

    def loop_size(cid):
        lo_, hi_ = cells[cid][0], cells[cid][1]
        h_ = hi_[0] - lo_[0]
        pts = set()
        for row in tris[tris[:, 9].astype(int) == cid]:
            for P in row[:9].reshape(3, 3):
                if sum(min(abs(P[k_] - lo_[k_]), abs(P[k_] - hi_[k_])) < 1e-9 * h_ for k_ in range(3)) >= 2:
                    pts.add(tuple(np.round(P, 8)))
        return len(pts)
    if cell is None:   # a cut cell whose loop has the most chords
        ids = np.unique(tris[:, 9].astype(int))
        cell = int(max(ids, key=lambda c_: (loop_size(c_), -c_)))
    lo, hi, loc, ct = cells[cell]
    if azim is None:   # view along the mean normal of the patch, slightly tilted
        mine0 = tris[tris[:, 9].astype(int) == cell]
        n = np.zeros(3)
        for row in mine0:
            T3 = row[:9].reshape(3, 3)
            n += np.cross(T3[1] - T3[0], T3[2] - T3[0])
        n /= np.linalg.norm(n)
        if n @ np.array([1.0, -1.0, 1.0]) < 0:
            n = -n
        V = OrthoView(math.degrees(math.atan2(n[1], n[0])) + 25, math.degrees(math.asin(n[2])) + 12)
    c = box_corners(lo, hi)
    h = hi[0] - lo[0]
    for f in BOX_FACES:
        V.polygon([c[i] for i in f], (0.93, 0.93, 0.93), alpha=0.35, flat=True)
    for i, j in BOX_EDGES:
        V.segment(c[i], c[j], '#000000', 0.8)
    mine = tris[tris[:, 9].astype(int) == cell]
    for row in mine:
        T3 = row[:9].reshape(3, 3)
        V.polygon(T3, ITF_BLUE, edge='#2f5f94', lw=0.6)
        for a_, b_ in [(0, 1), (1, 2), (2, 0)]:
            p, q = T3[a_], T3[b_]
            on_face = any(abs(p[d_] - v) < 1e-9 * h and abs(q[d_] - v) < 1e-9 * h
                          for d_ in range(3) for v in (lo[d_], hi[d_]))
            if on_face:
                V.segment(p, q, '#ac1117', 1.6, subdiv=1, bias=1e-6)
    for row in mine:   # vertices of the chords (on the edges of the cell)
        for P in row[:9].reshape(3, 3):
            on_edge = sum(min(abs(P[d_] - lo[d_]), abs(P[d_] - hi[d_])) < 1e-9 * h for d_ in range(3)) >= 2
            if on_edge:
                xy, z = V.proj(P)
                ax.plot(xy[0, 0], xy[0, 1], 'o', color='#ac1117', markersize=2.5, zorder=5)
    V.draw(ax, margin=0.05 * h)
    return cell


# figures of the paper, in the order of the paper
PAPER_FIGURES = ['interface_circle', 'interface_flower', 'interface_ellipse', 'interface_rounded',
                 'interface_annulus', 'meshes2d', 'meshes3d', 'validation2d', 'weight', 'extension2d', 'conv3d_hex',
                 'conv3d_tet', 'contrast', 'geometries', 'efficiency', 'efficiency_times']


def energy_norm_contrast(dim, kap2, R=1.0 / 3.0, kap1=1.0):
    """|| kappa^{1/2} grad u || for u_i = rho^6/kappa_i (test case 1) on (0,1)^dim"""
    import numpy as np
    x, w = np.polynomial.legendre.leggauss(8)
    x, w = 0.5 * (x + 1), 0.5 * w
    if dim == 2:
        X, Y = np.meshgrid(x, x, indexing='ij')
        W = np.outer(w, w)
        box = np.sum(W * 36 * ((X - .5)**2 + (Y - .5)**2)**5)
        ball = 6 * math.pi * R**12
    else:
        X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
        W = w[:, None, None] * w[None, :, None] * w[None, None, :]
        box = np.sum(W * 36 * ((X - .5)**2 + (Y - .5)**2 + (Z - .5)**2)**5)
        ball = 36 * 4 * math.pi * R**13 / 13
    return math.sqrt(ball / kap1 + (box - ball) / kap2)


def plot(outdir):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import plot_bvc_results as P
    import compare_bvc_speed as CS
    paper_style(plt)

    figs = os.path.join(outdir, 'paper', 'figures')          # paper
    rfigs = os.path.join(outdir, 'paper', 'report_figures')  # report (validation details)
    os.makedirs(figs, exist_ok=True)
    os.makedirs(rfigs, exist_ok=True)

    def T(name):
        return P.parse_table(os.path.join(outdir, name, 'bvc_convergence.txt'))

    def merged(prefix, k, var, fam=None):
        """rows of the runs <prefix>_k<k>_N<N0>, of the family fam of meshes"""
        rows = []
        for n in sorted(os.listdir(outdir)):
            if re.match(rf'{prefix}_k{k}_N\d+$', n):
                rows += T(n).get((k, var), [])
        uniq = {r['N']: r for r in rows if fam is None or family(r['N']) == fam}
        return sorted(uniq.values(), key=lambda r: -r['h'])

    saved = []

    def save(fig, name, rect=(0, 0, 1, 1), report=False):
        saved.append(name)
        log_minor_labels(fig)
        fig.tight_layout(rect=rect)
        d = rfigs if report else figs
        fig.savefig(os.path.join(d, name + '.pdf'))
        fig.savefig(os.path.join(d, name + '.png'), dpi=150)
        plt.close(fig)

    KSTY = [(BLACK, 'o'), (RED, '^'), (BLUE, 'D'), (GREEN, 's')]

    def legend_k(ax, ks, styles=(), loc='lower right', fontsize=6):
        """legend: one color per k, then the line styles of the compared methods"""
        h = [Line2D([], [], color=KSTY[k][0], marker=KSTY[k][1], label=f'$k={k}$') for k in ks]
        h += [Line2D([], [], color='#444444', marker=m, linestyle=ls, markerfacecolor=fc, label=lab)
              for lab, ls, m, fc in styles]
        ax.legend(handles=h, loc=loc, fontsize=fontsize, handlelength=1.8, borderpad=0.3, labelspacing=0.25)

    def common(lists):
        """the rows of the meshes (N) common to all the lists"""
        Ns = None
        for rows in lists:
            n = {r['N'] for r in rows}
            Ns = n if Ns is None else Ns & n
        return [sorted([r for r in rows if r['N'] in (Ns or set())], key=lambda r: -r['h']) for rows in lists]

    def kpanel(ax, get_main, get_cmp, q, ks, title, styles, cmp_markers_only=False):
        """all the degrees on one graph, on the meshes common to all the curves: main method
        (solid, observed orders), compared method (dashed, open markers; markers only, at the
        available meshes, if cmp_markers_only); legend in the panel"""
        ks = list(ks)
        mains = [get_main(k) for k in ks]
        cmps = [get_cmp(k) if get_cmp is not None else [] for k in ks]
        if get_cmp is not None and not cmp_markers_only:
            both = common(mains + cmps)
            mains, cmps = both[:len(ks)], both[len(ks):]
        else:
            mains = common(mains)
        for k, r, c in zip(ks, mains, cmps):
            col, mk = KSTY[k]
            curve(ax, r, q, col, mk, '-', orders='above')
            ref_slope(ax, [r], q, k + QOPT[q], factor=0.4)
            if get_cmp is not None:
                c = sorted([x for x in c if x.get(q, 0) > 0], key=lambda x: x['h'])
                if cmp_markers_only:
                    ax.loglog([x['h'] for x in c], [x[q] for x in c], linestyle='none', marker='o', markersize=9,
                              markerfacecolor='none', markeredgecolor=col, markeredgewidth=0.9, zorder=4)
                else:
                    curve(ax, c, q, col, mk, '--', open_marker=True)
        ax.set_ylabel(QLAB[q])
        frame_h(ax, left=1.3)
        tag(ax, title)
        legend_k(ax, ks, styles)

    def union(*lists):
        """rows of several runs, merged by mesh (N)"""
        u = {}
        for rows in lists:
            for r in rows:
                u[r['N']] = r
        return sorted(u.values(), key=lambda r: -r['h'])

    # ---------------- meshes and flat interfaces
    def cells_file(d):
        f = os.path.join(outdir, d, 'cells_l0.txt')
        return f if os.path.exists(f) else None
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.3))
    for ax, (d, title) in zip(axs, [('verif_mesh_circle_p2', 'circle, agglomeration'),
                                    ('verif_mesh_circle_p3', 'circle, extension'),
                                    ('verif_mesh_flower_p2', 'flower, agglomeration'),
                                    ('verif_mesh_flower_p3', 'flower, extension')]):
        f = cells_file(d)
        if not f:
            continue
        cells, itf = {}, []
        for line in open(f):
            t = line.split()
            if t[0] == 'cell':
                cells[int(t[1])] = dict(box=list(map(float, t[2:6])), loc=int(t[6]), cut=int(t[7]),
                                        stab=int(t[8]), group=int(t[9]))
            elif t[0] == 'interface':
                pts = list(map(float, t[2:]))
                itf.append((pts[0::2], pts[1::2]))
        for c in cells.values():
            x0, y0, x1, y1 = c['box']
            fc = '#ffffff'
            if c['group'] >= 0:
                fc = '#cfcfcf'
            elif c['cut'] == 1:
                fc = '#a9c6e8'
            elif c['cut'] == 2:
                fc = '#f2b5b0'
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=fc, edgecolor='#8c8c8c', linewidth=0.4))
        for xs, ys in itf:
            ax.plot(xs, ys, color='black', linewidth=1.0)
        for i, c in cells.items():
            if c['stab'] >= 0 and 'extension' in title:
                a, b = cells[i]['box'], cells[c['stab']]['box']
                ax.annotate('', xy=((b[0] + b[2]) / 2, (b[1] + b[3]) / 2), xytext=((a[0] + a[2]) / 2, (a[1] + a[3]) / 2),
                            arrowprops=dict(arrowstyle='->', color=RED, lw=0.8))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    save(fig, 'meshes2d')

    # 3D cut meshes (hexahedra, N = 8): background mesh on the back walls, cut
    # cells, flat interface; one cut cell with the loop of chords of its cut faces
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.7])
    for j, (d, title, el) in enumerate([('verif_mesh3d_sphere', 'sphere', 28), ('verif_mesh3d_flower', '3D flower', 40)]):
        cf, tf = os.path.join(outdir, d, 'cells_l0.txt'), os.path.join(outdir, d, 'flat_interface_l0.txt')
        if os.path.exists(cf) and os.path.exists(tf):
            ax = fig.add_subplot(gs[0, j])
            render_cut_mesh(ax, cf, tf, 8, elev=el)
            ax.set_title(f'{title}, $N = 8$', loc='center', fontsize=8)
    cf = os.path.join(outdir, 'verif_mesh3d_sphere', 'cells_l0.txt')
    tf = os.path.join(outdir, 'verif_mesh3d_sphere', 'flat_interface_l0.txt')
    if os.path.exists(cf) and os.path.exists(tf):
        ax = fig.add_subplot(gs[0, 2])
        render_cut_cell(ax, cf, tf)
        ax.set_title('one cut cell of the sphere', loc='center', fontsize=8)
    save(fig, 'meshes3d')

    # thumbnails of the interfaces, one figure per interface: 2D (left), 3D (right)
    import numpy as np
    TH = [(0, 'circle'), (7, 'flower'), (4, 'ellipse'), (5, 'rounded'), (6, 'annulus')]
    for i, tag_ in TH:
        fig, axs = plt.subplots(1, 2, figsize=(2.4, 1.15), gridspec_kw=dict(wspace=0.05))
        ax = axs[0]
        f = os.path.join(outdir, f'thumb2d_i{i}', 'cells_l0.txt')
        if os.path.exists(f):
            for line in open(f):
                t = line.split()
                if t[0] == 'interface':
                    pts = list(map(float, t[2:]))
                    ax.plot(pts[0::2], pts[1::2], color=BLUE, linewidth=1.1, solid_capstyle='round')
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor='#666666', linewidth=0.5))
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_aspect('equal')
        ax.set_axis_off()
        ax = axs[1]
        f = os.path.join(outdir, f'thumb3d_i{i}', 'flat_interface_l0.txt')
        if os.path.exists(f):
            V = OrthoView(-60.0, 40.0 if i == 7 else 28.0)
            for row in np.loadtxt(f):
                V.polygon(row[:9].reshape(3, 3), ITF_BLUE, edge='none', ambient=0.35)
            for a_, b_ in BOX_EDGES:   # the domain (0,1)^3
                ca = np.array([a_ & 1, (a_ >> 1) & 1, (a_ >> 2) & 1], dtype=float)
                cb = np.array([b_ & 1, (b_ >> 1) & 1, (b_ >> 2) & 1], dtype=float)
                V.segment(ca, cb, '#999999', 0.4)
            V.draw(ax, margin=0.02)
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.98)
        saved.append(f'interface_{tag_}')
        fig.savefig(os.path.join(figs, f'interface_{tag_}.pdf'))
        fig.savefig(os.path.join(figs, f'interface_{tag_}.png'), dpi=150)
        plt.close(fig)

    # ================= A. validation: the 2D test cases of [Bur_BVC] with its settings
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2))
    for i, (prob, title) in enumerate([('circle', r'circle, $\kappa_2/\kappa_1 = 10^4$'),
                                       ('flower', 'flower, discontinuous')]):
        for j, q in enumerate(['H1', 'L2']):
            kpanel(axs[i, j], lambda k, prob=prob: T(f'burman_{prob}_full').get((k, 'full'), []), None,
                   q, range(4), title, [])
    save(fig, 'validation2d')

    # ================= B. improvements
    # B1: weight of s°: kappa_i/h_T of [Bur_BVC] against (k+1) kappa_i/h_T
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    WST = [(r'$(k{+}1)\kappa_i/h_T$', '-', 'o', BLACK), (r'$\kappa_i/h_T$', '--', 'o', 'white')]
    kpanel(axs[0], lambda k: T('stab_flower_Ak').get((k, 'full'), []),
           lambda k: T('stab_flower_A1').get((k, 'full'), []), 'H1', [1, 2, 3], '2D flower', WST)
    kpanel(axs[1], lambda k: union(T(f'stab_sphere_hex_Ak_k{k}').get((k, 'full'), []),
                                   T(f'stab_sphere_hex_Ak_k{k}_N5').get((k, 'full'), [])),
           lambda k: union(T(f'stab_sphere_hex_A1_k{k}').get((k, 'full'), []),
                           T(f'stab_sphere_hex_A1_k{k}_N5').get((k, 'full'), [])),
           'H1', [1, 2, 3], 'sphere, hexahedra', WST)
    save(fig, 'weight')

    # B2: polynomial extension against cell agglomeration (2D, same weight)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2))
    for i, (prob, title) in enumerate([('circle', r'circle, $\kappa_2/\kappa_1 = 10^4$'),
                                       ('flower', 'flower, discontinuous')]):
        d3, d2 = T(f'conv2d_{prob}_p3'), T(f'conv2d_{prob}_p2')
        for j, q in enumerate(['H1', 'L2']):
            kpanel(axs[i, j], lambda k, d3=d3: d3.get((k, 'full'), []), lambda k, d2=d2: d2.get((k, 'full'), []),
                   q, range(4), title, [('extension', '-', 'o', BLACK), ('agglomeration', '--', 'o', 'white')])
    save(fig, 'extension2d')

    # ================= C. new test cases
    # C1, C2: 3D, same meshes for all the degrees (one family of meshes per panel)
    for M, mesh in [(3, 'hex'), (4, 'tet')]:
        fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2))
        for i, (prob, title) in enumerate([('sphere', r'sphere, $\kappa_2/\kappa_1 = 10^4$'),
                                           ('flower', '3D flower, discontinuous')]):
            for j, q in enumerate(['H1', 'L2']):
                kpanel(axs[i, j],
                       lambda k, prob=prob: merged(f'conv3d_{prob}_M{M}_p3_full', k, 'full'),
                       lambda k, prob=prob: merged(f'conv3d_{prob}_M{M}_p2_full', k, 'full'),
                       q, range(4), title, [('extension', '-', 'o', BLACK), ('agglomeration', '--', 'o', 'white')])
        save(fig, f'conv3d_{mesh}')

    # C3: contrast, relative energy error, kappa_2/kappa_1 = 10^-8 ... 10^8
    MS = list(range(9))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    for ax, (pre, dim, title) in zip(axs, [('contrast2d', 2, 'circle, $N = 64$'),
                                           ('contrast3d', 3, 'sphere, hexahedra, $N = 16$')]):
        for k in range(4):
            xs, es = [], []
            for m in MS:
                r = T(f'{pre}_m{m}').get((k, 'full'), [])
                if r and r[-1]['EN'] > 0:
                    xs.append(10.0**m)
                    es.append(r[-1]['EN'] / energy_norm_contrast(dim, 10.0**m))
            col, mk = KSTY[k]
            ax.loglog(xs, es, color=col, marker=mk, zorder=3, label=f'$k={k}$')
        ax.set_xticks([10.0**m for m in MS[::2]])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlabel(r'$\kappa_2/\kappa_1$')
        ax.set_ylabel('relative energy error')
        ys = [y for l in ax.get_lines() for y in l.get_ydata() if y > 0]
        if ys:
            ax.set_ylim(min(ys) / 300, max(ys) * 30)
        tag(ax, title)
        ax.legend(loc='lower right', ncol=4, fontsize=6, handlelength=1.8, columnspacing=1.0)
    save(fig, 'contrast')

    # C4: other interfaces, same meshes for all the degrees
    names = {4: ('ellipse', 'ellipsoid'), 5: ('rounded square', 'rounded cube'), 6: ('annulus', 'torus')}
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.2))
    for j, i in enumerate([4, 5, 6]):
        d = T(f'geom2d_i{i}')
        for k in [1, 2, 3]:
            col, mk = KSTY[k]
            r = d.get((k, 'full'), [])
            curve(axs[0, j], r, 'H1', col, mk, '-', orders='above')
            ref_slope(axs[0, j], [r], 'H1', k + 1, factor=0.4)
        rows3 = common([union(T(f'geom3d_i{i}_k{k}').get((k, 'full'), []),
                              T(f'geom3d_i{i}_k{k}_N12').get((k, 'full'), []),
                              T(f'geom3d_i{i}_k{k}_N32').get((k, 'full'), [])) for k in [1, 2, 3]])
        for k, r in zip([1, 2, 3], rows3):
            col, mk = KSTY[k]
            curve(axs[1, j], r, 'H1', col, mk, '-', orders='above')
            ref_slope(axs[1, j], [r], 'H1', k + 1, factor=0.4)
        for row in range(2):
            axs[row, j].set_ylabel(QLAB['H1'])
            frame_h(axs[row, j], left=1.3)
            tag(axs[row, j], names[i][row])
            legend_k(axs[row, j], [1, 2, 3])
    save(fig, 'geometries')

    # C5: efficiency: boundary value correction (flat interface) against the quadratures on a
    # refined interface (polyline of 2^r segments per cut face/cell) and on the exact curved interface
    METH = [('p3', 'flat, BVC', RED, '^', '-', False), ('p0', 'curved, exact', BLUE, 'D', '--', True),
            ('g0r8', 'polyline $r = 8$', BLACK, 'o', '--', True), ('g0r4', 'polyline $r = 4$', GRAY, 's', '--', True)]
    EFF = [('2d_flower', 3, '2D flower'), ('3d_sphere_hex', 2, 'sphere, hexahedra'),
           ('3d_sphere_tet', 3, 'sphere, tetrahedra'), ('3d_flower_hex', 2, '3D flower, hexahedra')]

    def eff_levels(case, m, k):
        u = {}
        for j, (N0, L) in enumerate(EFF_CASES[case][3]):
            for lv in CS.parse(os.path.join(outdir, eff_name(case, m, k, N0, j == 0))):
                if 'h' in lv:
                    u[round(lv['h'], 12)] = lv
        return [u[h] for h in sorted(u, reverse=True)]
    eff = {(case, m): eff_levels(case, m, k) for case, k, _ in EFF for m, *_ in METH}
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0))
    for ax, (case, k, title) in zip(axs.flat, EFF):
        for m, lab, col, mk, ls, op in METH:
            lv = eff[(case, m)]
            if not lv:
                continue
            ax.loglog([x['t_mesh'] + x['t_asm'] + x['t_solve'] for x in lv], [x.get('H1', float('nan')) for x in lv],
                      color=col, marker=mk, linestyle=ls, label=lab, zorder=3,
                      markerfacecolor='white' if op else col, markeredgecolor=col, markeredgewidth=1.0)
        ax.set_xlabel('wall time [s]')
        ax.set_ylabel(QLAB['H1'])
        tag(ax, f'{title}, $k={k}$', loc='upper right')
        ax.legend(loc='lower left', fontsize=6, handlelength=2.2, borderpad=0.3, labelspacing=0.25)
    save(fig, 'efficiency')

    fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.9))
    PARTS = [('t_mesh', 'geometry', '#d9d9d9'), ('t_asm', 'assembly', '#8c8c8c'), ('t_solve', 'solution', '#262626')]
    for ax, (case, k, title) in zip(axs, EFF):
        for j, (m, lab, col, mk, ls, op) in enumerate(METH):
            lv = eff[(case, m)]
            if not lv:
                continue
            x = lv[-1]
            bottom = 0
            for key, plab, pc in PARTS:
                ax.bar(j, x[key], bottom=bottom, width=0.7, color=pc, edgecolor=col, linewidth=1.2,
                       label=plab if j == 0 else None, zorder=3)
                bottom += x[key]
            ax.annotate(f"{x.get('H1', float('nan')):.1e}".replace('e-0', 'e-'), (j, bottom), xytext=(0, 2),
                        textcoords='offset points', ha='center', va='bottom', fontsize=5.5, rotation=90)
        ax.set_xticks(range(len(METH)))
        ax.set_xticklabels(['BVC', 'exact', '$r=8$', '$r=4$'], fontsize=6.5)
        ax.grid(axis='x', visible=False)
        ax.set_title(f'{title}, $k={k}$', fontsize=7.5)
        top = max(sum(lv[-1][key] for key, _, _ in PARTS) for lv in (eff[(case, m)] for m, *_ in METH) if lv)
        ax.set_ylim(0, 1.75 * top)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(facecolor=pc, edgecolor='black', linewidth=0.6, label=plab) for _, plab, pc in PARTS],
                  loc='upper left', fontsize=5.5, handlelength=1.2, borderpad=0.3, labelspacing=0.2)
    axs[0].set_ylabel('wall time [s]')
    save(fig, 'efficiency_times')

    # ================= report: details of the validation of [Bur_BVC] (its settings)
    # corrections one at a time
    VAR = [('naive', 'no correction', BLACK, 'o', '--', True), ('gam', r'$\gamma_T$ only', BLUE, 'D', '-', False),
           ('taylor', 'Taylor only', GREEN, 's', '-', False), ('full', 'full method', RED, '^', '-', False)]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0))
    for row, prob in enumerate(['circle', 'flower']):
        for col_, k in enumerate([2, 3]):
            ax = axs[row, col_]
            for v, lab, col, mk, ls, op in VAR:
                curve(ax, T(f'burman_{prob}_{v}').get((k, v), []), 'H1', col, mk, ls, label=lab, open_marker=op,
                      orders={'full': 'above', 'naive': 'below'}.get(v))
            ax.set_ylabel(QLAB['H1'])
            frame_h(ax, left=1.3)
            tag(ax, f'{prob}, $k={k}$')
            ax.legend(loc='lower right', fontsize=6)
    save(fig, 'burman_corrections', report=True)
    # interface refinement without correction
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, prob in zip(axs, ['circle', 'flower']):
        curve(ax, T(f'burman_{prob}_full').get((3, 'full'), []), 'H1', RED, '^', '-', label='full method, 1 chord',
              orders='above')
        curve(ax, T(f'burman_{prob}_naive').get((3, 'naive'), []), 'H1', BLACK, 'o', '--',
              label='no correction, 1 chord', open_marker=True)
        for n, c, col in [(2, 4, '#4d4d4d'), (3, 8, '#7f7f7f'), (4, 16, '#ababab')]:
            curve(ax, T(f'burman_{prob}_naive_n{n}').get((3, 'naive'), []), 'H1', col, 'o', '--',
                  label=f'no correction, {c} chords')
        ax.set_ylabel(QLAB['H1'])
        frame_h(ax, left=1.3)
        tag(ax, f'{prob}, $k=3$')
        ax.legend(loc='lower right', fontsize=6)
    save(fig, 'burman_chords', report=True)
    # Taylor order
    TM = [(0, BLACK, 'o', '--', True), (1, BLUE, 'D', '-', False), (2, RED, '^', '-', False)]
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.2))
    for i, prob in enumerate(['circle', 'flower']):
        for j, k in enumerate([1, 2, 3]):
            ax = axs[i, j]
            for m, col, mk, ls, op in TM:
                curve(ax, T(f'burman_{prob}_T{m}').get((k, 'gam'), []), 'H1', col, mk, ls, label=f'$m={m}$',
                      open_marker=op, orders='above' if m == [0, 1, 1, 2][k] else None)
            ax.set_ylabel(QLAB['H1'])
            frame_h(ax, left=1.3)
            tag(ax, f'{prob}, $k={k}$')
            ax.legend(loc='lower right', fontsize=6)
    save(fig, 'burman_taylor', report=True)
    # contrast, settings of [Bur_BVC]
    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    for k in range(4):
        es = []
        for m in [0, 2, 4, 6]:
            r = T(f'burman_contrast_m{m}').get((k, 'full'), [])
            es.append(r[-1]['H1'] if r else float('nan'))
        ax.loglog([1, 1e2, 1e4, 1e6], es, color=KSTY[k][0], marker=KSTY[k][1], label=f'$k={k}$')
    ax.set_xlabel(r'$\kappa_2/\kappa_1$')
    ax.set_ylabel(QLAB['H1'])
    tag(ax, 'circle, $N = 64$')
    ax.legend(loc='upper right', ncol=2, fontsize=6.5)
    save(fig, 'burman_contrast', report=True)
    # geometric quantities, normalized, same meshes in each panel
    def geo(name, Ns=None):
        f = os.path.join(outdir, name, 'bvc_geometry.txt')
        rows = []
        if os.path.exists(f):
            for line in open(f):
                if line.startswith('#') or not line.strip():
                    continue
                t = line.split()
                h = float(t[1])
                rows.append(dict(h=h, dh2=float(t[10]), th=float(t[12]), rrho=float(t[13]), rth=float(t[14])))
        return rows
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6))
    G2 = [('burman_geo_circle', 'circle', BLACK, 'o'), ('burman_geo_flower', 'flower', RED, '^')]
    G3 = [('verif_geo_sphere_hex', 'sphere, hex.', BLUE, 'D'), ('verif_geo_sphere_tet', 'sphere, tet.', BLUE, 'd'),
          ('verif_geo_flower_hex', '3D flower, hex.', GREEN, 's'), ('verif_geo_flower_tet', '3D flower, tet.', GREEN, 'v')]
    for row, (G, hmin) in enumerate([(G2, 0), (G3, 0.05)]):
        for c, (qmax, qrms, ylab) in enumerate([('dh2', 'rrho', r'$\max|\varrho_h|/h^2$ (solid), rms $/h^2$ (dashed)'),
                                                ('th', 'rth', r'$\max|\mathbf{n}_h-\mathbf{n}^*|/h$ (solid), rms $/h$ (dashed)')]):
            ax = axs[row, c]
            for name, lab, col, mk in G:
                r = [x for x in geo(name) if x['h'] >= hmin and x['h'] <= 0.25]
                ax.semilogx([x['h'] for x in r], [x[qmax] for x in r], color=col, marker=mk, label=lab)
                ax.semilogx([x['h'] for x in r], [x[qrms] for x in r], color=col, marker=mk, linestyle='--',
                            markerfacecolor='white')
            ax.set_xlabel('$h$')
            ax.set_ylabel(ylab, fontsize=7.5)
            ax.set_ylim(bottom=0)
            tag(ax, '2D' if row == 0 else '3D, $N = 8, 16, 32$', loc='upper right')
            ax.legend(loc='upper left', fontsize=6)
    save(fig, 'geometry', report=True)
    # patch test
    MESHES = [('quad', 'quadrangles', BLACK), ('tri', 'triangles', RED), ('hex', 'hexahedra', BLUE),
              ('tet', 'tetrahedra', GREEN)]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), sharey=True)
    for ax, (p, lab) in zip(axs, [(2, 'cell agglomeration'), (3, 'polynomial extension')]):
        w = 0.19
        for j, (M, mesh, col) in enumerate(MESHES):
            d = T(f'verif_patch_{M}_p{p}')
            vals = [max((x['H1'] for x in d.get((k, 'full'), [])), default=float('nan')) for k in range(4)]
            ax.bar([k + (j - 1.5) * w for k in range(4)], vals, width=w, color=col, edgecolor='black',
                   linewidth=0.4, label=mesh, zorder=3)
        ax.set_yscale('log')
        ax.set_ylim(1e-16, 1e-7)
        ax.set_xticks(range(4))
        ax.set_xticklabels([f'$k={k}$' for k in range(4)])
        ax.grid(axis='x', visible=False)
        tag(ax, lab)
        ax.legend(loc='upper right', ncol=2, fontsize=6)
    axs[0].set_ylabel(r'$\max_{N}\, |u-u_h|_{H^1}$')
    save(fig, 'patch', report=True)

    print('figures in', figs, 'and', rfigs)
    return [n for n in saved if n in PAPER_FIGURES]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('binary')
    ap.add_argument('outdir')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--mem', type=float, default=60, help='memory budget (GB)')
    ap.add_argument('--threads', type=int, default=os.cpu_count(), help='thread budget (parallel phase)')
    ap.add_argument('--eff-threads', type=int, default=12, help='threads of the efficiency runs')
    ap.add_argument('--only', default=None, help='comma-separated groups')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--paper', default=None, help='LaTeX folder: copy the figures (pdf) into DIR/figures')
    a = ap.parse_args()
    groups = a.only.split(',') if a.only else None
    J = [j for j in jobs() if groups is None or j[0] in groups]
    E = [j for j in eff_jobs() if groups is None or 'efficiency' in groups]
    if a.list:
        for g, n, x, m, t in J + E:
            print(f'{g:12s} {n:40s} {m:5.1f} GB   {x}')
        print(len(J) + len(E), 'runs')
        sys.exit(0)
    if a.run:
        os.makedirs(a.outdir, exist_ok=True)
        # a run whose arguments changed is run again
        done_path = os.path.join(a.outdir, 'done.txt')
        if os.path.exists(done_path):
            changed = set()
            for g, n, x, m, t in J + E:
                f = os.path.join(a.outdir, n, 'args.txt')
                if os.path.exists(f) and open(f).read().strip() != x.strip():
                    changed.add(n)
            if changed:
                lines = [l for l in open(done_path) if l.split()[0] not in changed]
                open(done_path, 'w').writelines(lines)
                print('arguments changed, run again:', ' '.join(sorted(changed)), flush=True)
        run_jobs([(n, x, m, t) for g, n, x, m, t in J], os.path.abspath(a.binary), a.outdir, a.mem, a.threads)
        # efficiency: one run at a time, fixed number of threads
        run_jobs([(n, x, m, a.eff_threads) for g, n, x, m, t in E], os.path.abspath(a.binary), a.outdir,
                 1e9, a.eff_threads)
    saved = plot(a.outdir) if a.plot else None
    if a.paper:
        import glob, shutil
        src = os.path.join(a.outdir, 'paper', 'figures')
        files = ([os.path.join(src, n + '.pdf') for n in saved] if saved is not None
                 else glob.glob(os.path.join(src, '*.pdf')))
        os.makedirs(os.path.join(a.paper, 'figures'), exist_ok=True)
        for f in files:
            shutil.copy(f, os.path.join(a.paper, 'figures'))
        print(len(files), 'figures copied into', os.path.join(a.paper, 'figures'))
