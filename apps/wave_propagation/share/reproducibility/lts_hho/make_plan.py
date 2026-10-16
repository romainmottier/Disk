#!/usr/bin/env python3
"""Tables (LaTeX) and figures of numerical_results_plan.tex (preprint LTS_HHO) from the runs of run_plan.sh and the
reduced models.
    usage: make_plan.py <out of run_plan.sh> <figure dir> [<modele_final dir>] [<dir of plan_extra.npz and figures>]
The head meshes (pictures) are read in <src>/mesh_generation/brain_mri/meshes."""
import os, re, glob, math, shutil, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = os.path.dirname(os.path.abspath(__file__))
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(P, "runs")
OUT = sys.argv[2] if len(sys.argv) > 2 else "plan_figures"
TOY = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser("~/Github/LTS_RK4/modele_final")
EXTRA = sys.argv[4] if len(sys.argv) > 4 else os.path.join(P, "reduced_models")
MESHES = os.path.join(P, "..", "..", "..", "src", "mesh_generation", "brain_mri", "meshes")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9})


def dat(d):
    f = os.path.join(R, d, "c.dat")
    if not os.path.exists(f): return None, []
    head, rows = None, []
    for l in open(f):
        if l.startswith("#"): head = l[1:].split(); continue
        if l.strip(): rows.append([float(x) for x in l.split()])
    return head, rows


def col(d, name, row=-1):
    h, r = dat(d)
    if h is None or not r or name not in h: return None
    return r[row][h.index(name)]


def cols(d, name):
    h, r = dat(d)
    if h is None or name not in h: return []
    return [x[h.index(name)] for x in r]


def tm(d):
    f = os.path.join(R, d, "time.txt")
    if not os.path.exists(f): return None, None
    s = open(f).read().split()
    try:
        i = s.index("s")
        return float(s[i - 1]), float(s[i + 1]) / 1024.0
    except (ValueError, IndexError):
        return None, None


def log(d):
    f = os.path.join(R, d, "log.txt")
    return re.sub(r"\x1b\[[0-9;]*m", "", open(f).read()) if os.path.exists(f) else ""


def fmt(x, f="{:.2e}"):
    return "--" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f.format(x)


def orders(h, e):
    return [math.log(e[i] / e[i + 1]) / math.log(h[i] / h[i + 1]) if e[i] > 0 and e[i + 1] > 0 else float("nan")
            for i in range(len(e) - 1)]


LABEL = {"erk": "ERK4", "lts": "LTS", "mlts": "MLTS", "li": "LI", "mlts_substep": "MLTS (per sub-window)",
         "lts_direct": "LTS (direct faces)"}


def write(name, s):
    open(os.path.join(OUT, name), "w").write(s)


# ------------------------------------------------------------------------------------------------- reduced models
def toy():
    for src, dst in [("comparaison/figures/rho_ordre2.png", "toy_rho_ordre2.png"), ("comparaison/figures/rho_ordre1.png", "toy_rho_ordre1.png"),
                     ("comparaison/figures/dtc_ordre2.png", "toy_dtc_ordre2.png"), ("comparaison/figures/dtc_ordre1.png", "toy_dtc_ordre1.png"),
                     ("mlts/figures/rho_order2_mlts.png", "toy_rho_order2_mlts.png"),
                     ("mlts/figures/rho_order1_buffer_mlts.png", "toy_rho_order1_buffer_mlts.png"),
                     ("mlts/figures/depth_vs_buffer.png", "toy_depth_vs_buffer.png")]:
        if os.path.exists(os.path.join(TOY, src)): shutil.copy(os.path.join(TOY, src), os.path.join(OUT, dst))
    for f in ("toy_li_dtc.png", "toy_mlts_alg4_am.png", "toy_order.png"):
        if os.path.exists(os.path.join(EXTRA, f)): shutil.copy(os.path.join(EXTRA, f), os.path.join(OUT, f))
    z = np.load(os.path.join(EXTRA, "plan_extra.npz"))
    li = {k: (z[f"li_{k}_rho"].max(), z[f"li_{k}_dtc"].min(), z[f"li_{k}_dtc"].max()) for k in ("order2", "order1_nobuf", "order1_buf")}
    ml = {(k, n): (z[f"mlts_{k}_{n}_rho"].max(), z[f"mlts_{k}_{n}_dtc"].min(), z[f"mlts_{k}_{n}_dtc"].max())
          for k in ("order2", "order1_buf") for n in ("alg4", "am")}
    names, errs = list(z["order_names"]), z["order_errs"]
    rows = []
    for n, e in zip(names, errs):
        o = orders([1, 0.5], [e[-2], e[-1]])[0]
        rows.append(f"{n} & {e[0]:.1e} & {e[-1]:.1e} & {o:.2f}\\\\")
    s = r"""\begin{itemize}
\item \emph{LI-SDIRK(3,4)} (Figure~\ref{fig:toy-li}), $p\leq2048$, $\varepsilon\in[0.01,0.5]$: order-2 model, $\max\rho(C_{\rm LI}(\Delta t_{\rm ref})) = %.4f$
and $\Delta t_c/\Delta t_{\rm ref}\in[%.2f,%.2f]$; order-1 model with overlap: $%.4f$, $[%.2f,%.2f]$; order-1 model
\emph{without} overlap: $\max\rho = %.1f$, $\Delta t_c/\Delta t_{\rm ref}$ down to $%.2f$. The implicit treatment of the
fine region does not remove the need for an overlap in the first-order form: the instability comes from the explicit
coarse polynomial (the fine coupling enters $\rm{w}^{(2)}$), which LI shares with LTS.
\item \emph{MLTS, Algorithm~4 vs Almquist--Mehlin} (Figure~\ref{fig:toy-new}, left), $p_1=p_2=p\leq64$: order 2,
$\max\rho = %.6f$ (Alg.~4) and $%.6f$ (A--M); order 1 with buffer, $\max\rho = %.6f$ and $%.6f$,
$\Delta t_c/\Delta t_{\rm ref}\in[%.2f,%.2f]$ and $[%.2f,%.2f]$. Same stability; the amplification matrices differ by
$\cal{O}(10^{-7})$ at $\Delta t = 0.3$ (two different order-4 schemes) and coincide for two levels.
\item \emph{Orders} (Figure~\ref{fig:toy-new}, right; relative error at $T=4$, $\Delta t = 0.4\cdot2^{-j}$, $j=0..5$):
\end{itemize}
\begin{center}\small\begin{tabular}{lccc}\toprule
model: scheme & error ($\Delta t=0.4$) & error ($\Delta t=0.0125$) & last order\\\midrule
%s
\bottomrule\end{tabular}\end{center}
LI reaches order 4 only once $\omega_{\rm fine}\Delta t$ is small (orders $3.95$--$3.99$ for $\Delta t\leq 3\times10^{-3}$;
at the coarse CFL limit, $\omega_{\rm fine}\Delta t = p\,\omega_{\rm coarse}\Delta t$ and the error of the fine modes is of
order one): LI is accurate when the solution has no high-frequency content in the fine region, which is the case of
the smooth and of the L-shape solutions below.
""" % (li["order2"] + li["order1_buf"] + (li["order1_nobuf"][0], li["order1_nobuf"][1]) + (ml[("order2", "alg4")][0], ml[("order2", "am")][0],
       ml[("order1_buf", "alg4")][0], ml[("order1_buf", "am")][0], ml[("order1_buf", "alg4")][1], ml[("order1_buf", "alg4")][2],
       ml[("order1_buf", "am")][1], ml[("order1_buf", "am")][2], "\n".join(rows)))
    write("toy_summary.tex", s)


# ------------------------------------------------------------------------------------------------- algebraic comparison
def alg():
    rows = []
    for o, lab in (("o1", "first order, $k=3$"), ("o2", "second order, $k=2$ mixed")):
        for l in range(4):
            d = f"alg/{o}_l{l}"
            t, mb = tm(d)
            lg = log(d)
            h = col(d, "h")
            if o == "o1":
                e1, e2 = col(d, "acoustic_p"), col(d, "acoustic_q")
                cpu, dofs, fdofs = col(d, "cpu"), col(d, "dofs"), col(d, "face_dofs")
                cells = int(dofs - fdofs) if dofs is not None and fdofs is not None else None
                faces = int(fdofs) if fdofs is not None else None
            else:
                e1, e2 = col(d, "acoustic_phi"), col(d, "acoustic_grad_phi")
                cpu = t
                m = re.search(r"unknowns: cells (\d+), dynamic faces (\d+), algebraic faces (\d+)", lg)
                cells, faces = (int(m.group(1)), int(m.group(3))) if m else (None, None)
            nt = 160 * 2 ** l
            rows.append(f"{lab if l == 0 else ''} & ${8 * 2 ** l}^2$ & {fmt(cells, '{:d}')} & {fmt(faces, '{:d}')} & "
                        f"{fmt(cpu / nt * 1e3 if cpu else None, '{:.3g}')} & {fmt(mb, '{:.0f}')} & {fmt(e1)} & {fmt(e2)}\\\\")
        if o == "o1": rows.append(r"\midrule")
    cfl = []
    for o in ("o1", "o2"):
        d = f"alg/{o}_cfl_l1"
        cfl.append((col(d, "stable_steps"), col(d, "cfl_h_max")))
    om = re.search(r"omega_max ~ ([0-9.e+-]+)", log("alg/o2_l1"))
    s = r"""\begin{table}[H]\centering\small
\caption{The two formulations on the unit square (acoustic, $t_f=1$, ERK4, $\Delta t = 1/(160\cdot2^l)$ on the
mesh $8\cdot2^l$ squared, fast modes, 6 threads). Unknowns: cell unknowns (first order: $(\bd{v},p)$; second order:
the state $(p,\phi)$) and face unknowns; time per step: CPU of the run (assembly included) divided by the number of
steps; memory: maximum resident size; errors: $L^2$ errors of $p$, $\bd{v}$ (first order) and of $u$, $\nabla u$
(second order) at $t_f$. Largest stable step on the $16^2$ mesh (CFL study, random data): first order $N=%s$ steps,
$c\Delta t/h = %s$; second order $N=%s$, $c\Delta t/h = %s$%s.}\label{tab:alg}
\begin{tabular}{llrrrrrr}\toprule
formulation & mesh & cell unk. & face unk. & ms/step & MB & error 1 & error 2\\\midrule
%s
\bottomrule\end{tabular}\end{table}
""" % (fmt(cfl[0][0], "{:.0f}"), fmt(cfl[0][1], "{:.3f}"), fmt(cfl[1][0], "{:.0f}"), fmt(cfl[1][1], "{:.3f}"),
       (r" ($\omega_{\max} = %s$, $2\sqrt2/\omega_{\max}$ on the $16^2$ mesh)" % om.group(1)) if om else "", "\n".join(rows))
    write("tab_alg.tex", s)


# ------------------------------------------------------------------------------------------------- stability
def stab():
    sch1 = [("erk", "ERK4"), ("lts", "LTS"), ("mlts", "MLTS"), ("mlts_o0", "MLTS, no ring"), ("li", "LI")]
    sch2 = [("erk", "ERK4"), ("lts", "LTS"), ("mlts", "MLTS"), ("li", "LI")]
    lines = []
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    for k, (o, sch, Ls) in enumerate((("stab1", sch1, range(5)), ("stab2", sch2, range(4)))):
        for s, lab in sch:
            hs, dts = [], []
            for L in Ls:
                d = f"{o}/L{L}_{s}"
                hmin, dt = col(d, "h_min"), col(d, "dt")
                st = col(d, "stable_steps")
                if hmin is not None and dt is not None and st is not None and st > 1:
                    hs.append(hmin); dts.append(dt)
            if hs: ax[k].loglog(hs, dts, "o-", ms=3, label=lab)
        ax[k].set_xlabel(r"$h_{\min}$"); ax[k].set_title("first order, $k=3$" if k == 0 else "second order, $k=2$ mixed")
        ax[k].legend(fontsize=7)
    ax[0].set_ylabel(r"largest stable $\Delta t$")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "stab_dtc.png"), dpi=150); plt.close(fig)
    for o, sch, Ls, lab in (("stab1", sch1, range(5), "first order"), ("stab2", sch2, range(4), "second order")):
        head = " & ".join(["$L$", "$h_{\\min}$"] + [x[1] for x in sch])
        rows = []
        for L in Ls:
            hmin = col(f"{o}/L{L}_erk", "h_min")
            vals = []
            for s, _ in sch:
                d = f"{o}/L{L}_{s}"
                dt, st = col(d, "dt"), col(d, "stable_steps")
                vals.append("$\\infty$" if st == 1 else fmt(dt, "{:.3e}"))
            rows.append(" & ".join([str(L), fmt(hmin, "{:.2e}")] + vals) + "\\\\")
        lines.append(r"\multicolumn{%d}{l}{\emph{%s}}\\ %s\\\midrule %s" % (len(sch) + 2, lab, head, "\n".join(rows)))
    s = r"""\begin{table}[H]\centering\small
\caption{Largest stable $\Delta t$ (CFL study, $t_f=1$, random data) on the $18\times18$ mesh with the box refined $L$
times. $\infty$: stable down to one step (LI at $L=0$: every cell is fine, the scheme is SDIRK(3,4)).}\label{tab:stab}
\begin{tabular}{%s}\toprule
%s
\bottomrule\end{tabular}\end{table}
""" % ("l" * 7, "\n\\midrule\n".join(lines))
    write("tab_stab.tex", s)


# ------------------------------------------------------------------------------------------------- smooth convergence
def conv():
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    out = []
    for k, (o, e1, e2, lab, schemes) in enumerate((("conv1", "acoustic_p", "acoustic_q", "first order", ["lts", "mlts", "li"]),
                                                  ("conv2", "acoustic_phi", "acoustic_grad_phi", "second order", ["lts", "mlts", "li"]))):
        rows = []
        for s in schemes:
            d = f"{o}/{s}"
            h, a, b = cols(d, "h"), cols(d, e1), cols(d, e2)
            if not h: continue
            ax[k].loglog(h, a, "o-", ms=3, label=f"{s} ({e1.split('_')[-1]})")
            ax[k].loglog(h, b, "s--", ms=3, label=f"{s} ({e2.split('_', 1)[-1]})")
            oa, ob = orders(h, a), orders(h, b)
            t, mb = tm(d)
            rows.append(f"{lab} & {LABEL.get(s, s)} & " + " ".join(f"{x:.2e}" for x in a) + " & " + " ".join(f"{x:.2f}" for x in oa)
                        + " & " + " ".join(f"{x:.2f}" for x in ob) + f" & {fmt(t, '{:.0f}')}\\\\")
        he, ae, be = [], [], []
        for l in range(1, 5):
            d = f"{o}/erk_l{l}"
            if col(d, "h") is not None: he.append(col(d, "h")); ae.append(col(d, e1)); be.append(col(d, e2))
        if he:
            ax[k].loglog(he, ae, "k^:", ms=4, label="ERK4 (p)"); ax[k].loglog(he, be, "kv:", ms=4, label="ERK4 (grad)")
            rows.append(f"{lab} & ERK4 & " + " ".join(f"{x:.2e}" for x in ae) + " & " + " ".join(f"{x:.2f}" for x in orders(he, ae))
                        + " & " + " ".join(f"{x:.2f}" for x in orders(he, be)) + " & \\\\")
        if he or rows:
            hh = np.array(cols(f"{o}/mlts", "h") or he)
            if len(hh): ax[k].loglog(hh, (ae[0] if ae else 1e-4) * (hh / hh[0]) ** 4, "k-", lw=0.5, label="slope 4")
        ax[k].set_xlabel(r"$h_{\max}$"); ax[k].set_title(lab); ax[k].legend(fontsize=6, ncol=2)
        out += rows
    ax[0].set_ylabel(r"$L^2$ error at $t_f$")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "conv_space.png"), dpi=150); plt.close(fig)
    s = r"""\begin{table}[H]\centering\scriptsize
\caption{Space convergence (haloscale meshes, $\Delta t = 0.05h_{\max}$; ERK4: $\Delta t = 0.05 h_{\min}$): $L^2$ errors
of $p$ (first order) or $u$ (second order) on the successive meshes, orders of these errors and of the second field
($\bd{v}$, $\nabla u$), CPU time (s) of the whole study (6 threads).}\label{tab:conv}
\begin{tabular}{llp{5.6cm}p{2.3cm}p{2.3cm}r}\toprule
formulation & scheme & errors & orders & orders (2nd field) & CPU\\\midrule
%s
\bottomrule\end{tabular}\end{table}
""" % "\n".join(out)
    write("tab_conv.tex", s)
    rows = []
    for o, lab in (("tconv1", "first order"), ("tconv2", "second order")):
        for s in ("lts", "lts_direct", "mlts", "mlts_substep", "li"):
            d = f"{o}/{s}"
            dt, e = cols(d, "dt"), cols(d, "mass_norm")
            if not dt: continue
            rows.append(f"{lab} & {LABEL.get(s, s)} & " + " ".join(f"{x:.2e}" for x in e) + " & " + " ".join(f"{x:.2f}" for x in orders(dt, e)) + "\\\\")
    s = r"""\begin{table}[H]\centering\small
\caption{Convergence in time (mesh $N=2$, $t_f=1$, $\Delta t = 1/256,\ldots,1/2048$, against the same scheme at
$\Delta t/4$ of the finest step): mass norm of the difference and orders.}\label{tab:tconv}
\begin{tabular}{llp{6.2cm}l}\toprule
formulation & scheme & difference & orders\\\midrule
%s
\bottomrule\end{tabular}\end{table}
""" % "\n".join(rows)
    write("tab_tconv.tex", s)


# ------------------------------------------------------------------------------------------------- L-shape
def lshape():
    rows = []
    data = {}
    for s in ("mlts", "li", "lts", "erk", "li_stable", "li_ref"):
        for l in range(1, 5):
            d = f"lshape/{s}_l{l}"
            h, ep, ev, cpu = col(d, "h"), col(d, "acoustic_p"), col(d, "acoustic_q"), col(d, "cpu")
            t, mb = tm(d)
            lg = log(d)
            m = re.search(r"LI fine system: (\d+) unknowns, nnz (\d+), nnz\(L \+ U\) (\d+) \(~([0-9.]+) MB\)", lg)
            lv = re.search(r"MLTS: (\d+) levels", lg)
            data[(s, l)] = (h, ep, ev, cpu, mb, m, lv)
    for l in range(1, 5):
        for s, lab in (("erk", "ERK4"), ("lts", "LTS"), ("mlts", "MLTS"), ("li", "LI"),
                       ("li_stable", r"LI, $\Delta t = 1.17\times10^{-4}$"), ("li_ref", r"LI ref., $\Delta t = 1.04\times10^{-3}$")):
            h, ep, ev, cpu, mb, m, lv = data[(s, l)]
            if h is None and cpu is None: continue
            extra = f"{m.group(1)} unk., {float(m.group(4)):.0f} MB" if m else (f"{lv.group(1)} levels" if lv else "")
            rows.append(f"{l - 1} & {lab} & {fmt(h, '{:.3f}')} & {fmt(ep)} & {fmt(ev)} & {fmt(cpu, '{:.1f}')} & {fmt(mb, '{:.0f}')} & {extra}\\\\")
        rows.append(r"\midrule")
    s = r"""\begin{table}[H]\centering\small
\caption{L-shaped domain: $L^2$ errors at $t_f=0.5$, CPU time of the run (s, 6 threads, assembly included), maximum
resident memory (MB), and the number of levels (MLTS) or the size and memory of the LU factorization of the fine
system (LI).}\label{tab:lshape}
\begin{tabular}{llrrrrrl}\toprule
$N$ & scheme & $h_{\max}$ & error $p$ & error $\bd{v}$ & CPU & MB & \\\midrule
%s
\bottomrule\end{tabular}\end{table}
""" % "\n".join(rows[:-1])
    write("tab_lshape.tex", s)
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for s, mk in (("mlts", "o-"), ("li", "s-"), ("lts", "^-"), ("erk", "v-")):
        # LI on N = 3: the run at its stability limit (the run at 0.05 h_max blows up)
        src = {l: (("li_stable", l) if (s == "li" and l == 4) else (s, l)) for l in range(1, 5)}
        pts = [(data[src[l]][0], data[src[l]][1], data[src[l]][2], data[src[l]][3]) for l in range(1, 5)
               if data[src[l]][1] is not None and data[src[l]][1] < 1.0]
        if not pts: continue
        h, ep, ev, cpu = map(np.array, zip(*pts))
        ax[0].loglog(h, ep, mk, ms=3, label=f"{s} (p)"); ax[0].loglog(h, ev, mk.replace("-", "--"), ms=3, label=f"{s} (v)")
        ax[1].loglog(cpu, ep, mk, ms=3, label=f"{s} (p)")
    ax[0].set_xlabel(r"$h_{\max}$"); ax[0].set_ylabel(r"$L^2$ error"); ax[0].legend(fontsize=6, ncol=2)
    ax[1].set_xlabel("CPU (s)"); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "lshape.png"), dpi=150); plt.close(fig)


# ------------------------------------------------------------------------------------------------- head application
def brain():
    rows = []
    cfgs = ((1, "b8: layered skull"), (2, "w8: temporal window"), (3, "i8: implant"))
    for lv, lab in cfgs:
        erk1 = None
        for s, sl in (("erk", "ERK4"), ("lts", "LTS (2 levels)"), ("mlts", "MLTS"), ("li", "LI")):
            c = f"brain/cfl_l{lv}_{s}"
            dt = col(c, "dt")
            w6, m6 = tm(f"brain/t6_l{lv}_{s}")
            h6, _ = tm(f"brain/h6_l{lv}_{s}")
            w1, _ = tm(f"brain/t1_l{lv}_{s}")
            h1, _ = tm(f"brain/h1_l{lv}_{s}")
            i6 = 2 * (w6 - h6) if w6 and h6 else None
            i1 = 2 * (w1 - h1) if w1 and h1 else None
            if s == "erk": erk1, erk6 = i1, i6
            lg = log(f"brain/t6_l{lv}_{s}") or log(c)
            m = re.search(r"MLTS: (\d+) levels, cells per level \(coarsest first\):([0-9 ]+)", lg)
            li = re.search(r"LI fine system: (\d+) unknowns.*?\(~([0-9.]+) MB\)", lg)
            info = (f"{m.group(1)} lev.: {m.group(2).strip()}" if m else "") + (f"{li.group(1)} unk., LU {float(li.group(2)):.0f} MB" if li else "")
            sp6 = f"{erk6 / i6:.2f}" if (i6 and erk6) else "--"
            sp1 = f"{erk1 / i1:.2f}" if (i1 and erk1) else "--"
            rows.append(f"{lab if s == 'erk' else ''} & {sl} & {fmt(dt, '{:.4f}')} & {fmt(w6, '{:.1f}')} & {fmt(i6, '{:.1f}')} ({sp6}) & "
                        f"{fmt(i1, '{:.1f}')} ({sp1}) & {fmt(m6, '{:.0f}')} & {info}\\\\")
        rows.append(r"\midrule")
    write("tab_brain.tex", r"""\begin{table}[H]\centering\scriptsize
\caption{Head slice: largest stable $\Delta t$ ($\mu$s, CFL study, random data); run to $t_f=100\,\mu$s at $90\,\%%$ of the
limit of each scheme without output: wall time (s, 6 threads, setup included), integration time
$2(T_{100}-T_{50})$ (s) with 6 threads and with 1 thread (in brackets: speed-up over ERK4), maximum memory (MB), levels
and cells per level, or size and memory of the LI factorization.}\label{tab:brain}
\begin{tabular}{llrrrrrp{4.6cm}}\toprule
mesh & scheme & $\Delta t_{\rm crit}$ & wall & integr. 6 thr. & integr. 1 thr. & MB & \\\midrule
%s
\bottomrule\end{tabular}\end{table}
""" % "\n".join(rows[:-1]))
    # sensors
    fig, ax = plt.subplots(3, 2, figsize=(10, 8))
    for k, (lv, name, fld) in enumerate(((1, "brain_center", "acoustic_pressure"), (1, "water_back", "acoustic_pressure"),
                                         (2, "brain_center", "acoustic_pressure"), (2, "skull_window", "elastic_velocity"),
                                         (3, "near_implant", "acoustic_pressure"), (3, "water_back", "acoustic_pressure"))):
        a = ax[k // 2][k % 2]
        for s, st in (("erk", "k-"), ("lts", "b--"), ("mlts", "r:"), ("li", "g-.")):
            f = os.path.join(R, f"brain/o_l{lv}_{s}", f"c_{name}.csv")
            if not os.path.exists(f): continue
            lines = [l for l in open(f) if not l.startswith("#") and l.strip()]
            A = np.array([[float(x) for x in l.split(",")] for l in lines[1:]])
            if A.size == 0: continue
            t = A[:, 0]
            y = A[:, 1] if A.shape[1] == 2 else np.sqrt((A[:, 1:] ** 2).sum(axis=1))     # |v| for the vector fields
            a.plot(t, y, st, lw=0.8, label=s)
        a.set_title(f"{['b8', 'w8', 'i8'][lv - 1]}: {name} ({fld})", fontsize=8); a.legend(fontsize=6)
        a.set_xlabel(r"$t$ ($\mu$s)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "brain_sensors.png"), dpi=150); plt.close(fig)
    # meshes and snapshots (MLTS runs with observers)
    for f, g in (("b8.png", "head_s80_layered_b8.png"), ("w8.png", "head_s80_window_w8.png"), ("i8.png", "head_s80_implant_i8.png")):
        src = os.path.join(MESHES, g)
        if os.path.exists(src): shutil.copy(src, os.path.join(OUT, "brain_mesh_" + f))
    from matplotlib.collections import PolyCollection
    for lv in (1, 2, 3):
        d = os.path.join(R, f"brain/o_l{lv}_mlts")
        cf = os.path.join(d, "c_cells.dat")
        if not os.path.exists(cf): continue
        cells = []
        for l in open(cf):
            if l.startswith("#"): continue
            v = l.split(); n = int(v[1])
            cells.append((int(v[0]), np.array(v[2:2 + 2 * n], float).reshape(n, 2)))
        snaps = sorted(glob.glob(os.path.join(d, "c_snapshot_*.dat")))
        if not snaps: continue
        fig, ax = plt.subplots(1, len(snaps), figsize=(3 * len(snaps), 3.6))
        for a, sf in zip(np.atleast_1d(ax), snaps):
            t = re.search(r"t = ([0-9.e+-]+)", open(sf).readline()).group(1)
            vals = [float(l.split()[0]) for l in open(sf) if not l.startswith("#")]
            vals = np.array(vals)
            ac = np.array([c[0] for c in cells]) == 1
            vmax = np.percentile(np.abs(vals[ac]), 99.5) or 1.0
            pc = PolyCollection([c[1] for c in cells], array=np.where(ac, vals, np.nan), cmap="RdBu_r", edgecolors="none")
            pc.set_clim(-vmax, vmax)
            a.add_collection(pc)
            pe = PolyCollection([c[1] for c, e in zip(cells, ~ac) if e], facecolors="none", edgecolors="k", linewidths=0.05)
            a.add_collection(pe)
            a.autoscale(); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
            a.set_title(f"t = {float(t):.0f} us", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f"brain_snapshots_l{lv}.png"), dpi=170); plt.close(fig)


if __name__ == "__main__":
    for f in (toy, alg, stab, conv, lshape, brain):
        try:
            f()
        except Exception as e:
            print(f.__name__, "failed:", repr(e))
    print("written to", OUT)
