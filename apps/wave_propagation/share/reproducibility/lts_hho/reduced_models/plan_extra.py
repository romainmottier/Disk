#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reduced models (no spatial discretization) for the plan of the numerical results of the LTS_HHO preprint:
  1. LI-SDIRK(3,4) (Section 5): amplification matrix C_LI(dt) of one step on the order-2 model (4x4), the order-1 model
     without overlap (4x4) and with overlap (6x6) of comparaison/lts_rk4_comparaison.py; rho(C_LI(dt_ref)) and the
     critical step dt_c(p) for p up to 2048.
  2. MLTS-RK4: Algorithm 4 of the preprint (own coefficients from powers of the full B, accumulated and Taylor-shifted)
     against the Almquist-Mehlin recursion of mlts/mlts_engine.py (powers of B P_{l-1} and the re-centred parent
     polynomial), on the 3-level models of mlts/compute_mlts.py: rho(C(dt_ref)) and dt_c for p1 = p2 = p.
  3. Order of accuracy: error at T of ERK4, LTS-RK4, MLTS (both forms) and LI-SDIRK against expm(B T) y0.
Results: plan_extra.npz and the figures in --out.
"""
import argparse, os, sys, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# engines of ~/Github/LTS_RK4/modele_final (comparaison/lts_rk4_comparaison.py, mlts/mlts_engine.py, mlts/compute_mlts.py)
MODELE_FINAL = os.environ.get("MODELE_FINAL", os.path.expanduser("~/Github/LTS_RK4/modele_final"))
sys.path.insert(0, os.path.join(MODELE_FINAL, "comparaison"))
sys.path.insert(0, os.path.join(MODELE_FINAL, "mlts"))
import lts_rk4_comparaison as two                                   # noqa: E402
import mlts_engine as am                                            # noqa: E402
from compute_mlts import (build_B_order2_3tier, build_B_order1_3tier, PS_3, PS_4)   # noqa: E402

EPS_LIST = [0.01, 0.05, 0.1, 0.25, 0.5]
P_GRID = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]


# ----------------------------------------------------------------------------------------------- LI-SDIRK(3,4)
def sdirk34():
    g = math.cos(math.pi / 18.0) / math.sqrt(3.0) + 0.5
    xi = 1.0 / (6.0 * (2.0 * g - 1.0) ** 2)
    a = np.array([[g, 0, 0], [0.5 - g, g, 0], [2 * g, 1 - 4 * g, g]])
    b = np.array([xi, 1 - 2 * xi, xi])
    c = np.array([g, 0.5, 1 - g])
    return a, b, c


def li_amplification(B, P, dts, coupling="taylor"):
    """One LI step: coarse Taylor polynomial w^(j) = B (I-P) B^j y (ERK4, degree 3), SDIRK(3,4) stages
    (I - g dt B P) Z_i = y + dt sum_{l<i} a_il k_l + g dt w_i, k_i = B P Z_i + w_i, y += dt sum b_i k_i, with
    w_i = w(c_i dt) (coupling "taylor", the original LI) or w_i = sum_j beta_ij w^(j), beta_i0 = 1,
    beta_ij = dt sum_{l<=i} a_il (c_l dt)^(j-1) / (j-1)! (coupling "stage": coarse stage values of the SDIRK tableau,
    time.li_coupling = "stage" of wave_propagation)."""
    a, b, c = sdirk34()
    n = B.shape[0]
    I = np.eye(n, dtype=B.dtype)
    Bc, Bf = B @ (I - P), B @ P
    W = [Bc, Bc @ B, Bc @ B @ B, Bc @ B @ B @ B]
    out = []
    for dt in np.atleast_1d(dts):
        A = I - a[0, 0] * dt * Bf
        K = []
        with np.errstate(over="ignore", invalid="ignore"):
            for i in range(3):
                tau = c[i] * dt
                if coupling == "stage":
                    be = [1.0] + [dt * sum(a[i, l] * (c[l] * dt) ** (j - 1) / math.factorial(j - 1) for l in range(i + 1))
                                  for j in (1, 2, 3)]
                    w = be[0] * W[0] + be[1] * W[1] + be[2] * W[2] + be[3] * W[3]
                else:
                    w = W[0] + tau * W[1] + tau ** 2 / 2 * W[2] + tau ** 3 / 6 * W[3]
                R = I + sum(dt * a[i, l] * K[l] for l in range(i)) if i else I.copy()
                Z = np.linalg.solve(A, R + a[i, i] * dt * w)
                K.append(Bf @ Z + w)
            out.append(I + dt * sum(b[i] * K[i] for i in range(3)))
    return np.array(out)


# ----------------------------------------------------------------------------------------------- MLTS, Algorithm 4
def taylor_shift(w, s):
    """Coefficients of sum_k t^k/k! w_k re-expanded at t = s."""
    return [sum(s ** (q - j) / math.factorial(q - j) * w[q] for q in range(j, 4)) for j in range(4)]


def mlts_alg4_amplification(B, Ps, ps, dts):
    """Algorithm 4 (window variant): Ps = [P_1..P_Lmax] (toy notation, P_l selects the levels >= l), ps = ratios.
    Level l < Lmax: own coefficients w_l^(i) = B (P_l - P_{l+1}) B^i Y (P_0 = I) from the current state, added to the
    inherited ones, Taylor-shifted to each sub-window; level Lmax: RK4 of y' = B P_Lmax y + w(tau)."""
    n = B.shape[0]
    I = np.eye(n, dtype=B.dtype)
    Lmax = len(Ps)
    proj = [I] + list(Ps) + [np.zeros_like(I)]

    def step(level, Y, h, w):
        if level == Lmax:
            Bf = B @ proj[Lmax]
            pw = lambda t: w[0] + t * w[1] + t * t / 2 * w[2] + t ** 3 / 6 * w[3]
            k1 = pw(0.0) + Bf @ Y
            k2 = pw(h / 2) + Bf @ (Y + h / 2 * k1)
            k3 = pw(h / 2) + Bf @ (Y + h / 2 * k2)
            k4 = pw(h) + Bf @ (Y + h * k3)
            return Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Q = proj[level] - proj[level + 1]
        g = [Y]
        for _ in range(3): g.append(B @ g[-1])
        w = [w[i] + B @ (Q @ g[i]) for i in range(4)]
        p = ps[level]
        hs = h / p
        for m in range(p):
            Y = step(level + 1, Y, hs, taylor_shift(w, m * hs))
        return Y

    out = []
    with np.errstate(over="ignore", invalid="ignore"):
        for dt in np.atleast_1d(dts):
            zero = np.zeros_like(I)
            out.append(step(0, I.copy(), dt, [zero] * 4))
    return np.array(out)


def rho(M):
    return two.rho_batch(M)


def dt_crit(amp, dt_max=4.0, n=200, tol=1e-5):
    dts = np.linspace(dt_max / n, dt_max, n)
    r = rho(amp(dts))
    bad = np.where(r > 1 + tol)[0]
    return dt_max if len(bad) == 0 else (0.0 if bad[0] == 0 else dts[bad[0] - 1])


# ----------------------------------------------------------------------------------------------- models
P_O1_NOBUF = np.diag([0.0, 1.0, 0.0, 1.0])


def build_B_order1_nobuf(eps, p):
    G = np.array([[1.0, eps], [eps, p * 1.0]], dtype=complex)
    B = np.zeros((4, 4), dtype=complex)
    B[0:2, 2:4] = G
    B[2:4, 0:2] = -G
    return B


LI_MODELS = {
    "order2": (two.build_B_order2, two.P_ORDER2, "order 2 (4x4)"),
    "order1_nobuf": (build_B_order1_nobuf, P_O1_NOBUF, "order 1, no overlap (4x4)"),
    "order1_buf": (two.build_B_order1, two.P_ORDER1, "order 1, overlap (6x6)"),
}
MLTS_MODELS = {
    "order2": (build_B_order2_3tier, PS_3, "order 2, 3 levels"),
    "order1_buf": (build_B_order1_3tier, PS_4, "order 1, 3 levels, overlap"),
}


def erk4_matrix(B, dt):
    I = np.eye(B.shape[0], dtype=B.dtype)
    z = dt * B
    return I + z + z @ z / 2 + z @ z @ z / 6 + z @ z @ z @ z / 24


def expm(M):
    w, V = np.linalg.eig(M)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    res = {}

    # 1. LI: rho at the p = 1 critical step of ERK4 / LTS (dt_ref) and dt_c(p) / dt_ref
    print("== LI-SDIRK(3,4) on the reduced models")
    for key, (build, P, label) in LI_MODELS.items():
        mod = dict(build_B=build, P=P)
        rr, dc, dref = [], [], []
        for e in EPS_LIST:
            dr = two.dt_ref(mod, e)
            dref.append(dr)
            rr.append([rho(li_amplification(build(e, p), P, [dr]))[0] for p in P_GRID])
            dc.append([dt_crit(lambda d, p=p: li_amplification(build(e, p), P, d)) / dr for p in P_GRID])
        res[f"li_{key}_rho"], res[f"li_{key}_dtc"], res[f"li_{key}_dtref"] = np.array(rr), np.array(dc), np.array(dref)
        print(f"  {label}: max rho(dt_ref) = {np.max(rr):.6f}, dt_c/dt_ref in [{np.min(dc):.3f}, {np.max(dc):.3f}]")

    # 2. MLTS Algorithm 4 vs Almquist-Mehlin, p1 = p2 = p
    print("== MLTS: Algorithm 4 vs Almquist-Mehlin (3 levels)")
    P3 = [1, 2, 4, 8, 16, 32, 64]
    for key, (build, Ps, label) in MLTS_MODELS.items():
        for name, amp in (("alg4", mlts_alg4_amplification), ("am", am.mlts_amplification_batch)):
            rr, dc = [], []
            for e in EPS_LIST:
                ones = lambda dts, e=e: amp(build([1, 1], e), Ps, [1, 1], dts)
                dr = dt_crit(ones, n=400)
                rr.append([rho(amp(build([p, p], e), Ps, [p, p], np.array([dr])))[0] for p in P3])
                dc.append([dt_crit(lambda d, p=p: amp(build([p, p], e), Ps, [p, p], d)) / dr for p in P3])
            res[f"mlts_{key}_{name}_rho"], res[f"mlts_{key}_{name}_dtc"] = np.array(rr), np.array(dc)
            print(f"  {label}, {name}: max rho(dt_ref) = {np.max(rr):.6f}, dt_c/dt_ref in [{np.min(dc):.3f}, {np.max(dc):.3f}]")
    # the two MLTS forms coincide for 2 levels
    B2 = two.build_B_order2(0.1, 8)
    d = np.abs(mlts_alg4_amplification(B2, [two.P_ORDER2], [8], [0.5]) - two.amplification_batch(B2, two.P_ORDER2, 8, [0.5])).max()
    print(f"  check: Algorithm 4 with 2 levels = LTS-RK4 (Algorithm 3): max |C_alg4 - C_lts| = {d:.2e}")
    res["check_alg4_lts"] = d

    # 3. order of accuracy, eps = 0.1, T = 4
    print("== order of accuracy (error at T = 4 against expm)")
    T, e = 4.0, 0.1
    rng = np.random.default_rng(1)
    cases = []
    B = two.build_B_order2(e, 8); y0 = rng.standard_normal(4)
    cases += [("order2", "ERK4 (dt/8)", B, lambda d, B=B: np.linalg.matrix_power(erk4_matrix(B, d / 8), 8), y0),
              ("order2", "LTS p=8", B, lambda d, B=B: two.amplification_batch(B, two.P_ORDER2, 8, [d])[0], y0),
              ("order2", "LI p=8", B, lambda d, B=B: li_amplification(B, two.P_ORDER2, [d])[0], y0)]
    B = build_B_order2_3tier([4, 4], e); y0 = rng.standard_normal(6)
    cases += [("order2", "MLTS alg4 4x4", B, lambda d, B=B: mlts_alg4_amplification(B, PS_3, [4, 4], [d])[0], y0),
              ("order2", "MLTS AM 4x4", B, lambda d, B=B: am.mlts_amplification_batch(B, PS_3, [4, 4], np.array([d]))[0], y0)]
    B = two.build_B_order1(e, 8); y0 = rng.standard_normal(6)
    cases += [("order1", "LTS p=8", B, lambda d, B=B: two.amplification_batch(B, two.P_ORDER1, 8, [d])[0], y0),
              ("order1", "LI p=8", B, lambda d, B=B: li_amplification(B, two.P_ORDER1, [d])[0], y0)]
    B = build_B_order1_3tier([4, 4], e); y0 = rng.standard_normal(8)
    cases += [("order1", "MLTS alg4 4x4", B, lambda d, B=B: mlts_alg4_amplification(B, PS_4, [4, 4], [d])[0], y0),
              ("order1", "MLTS AM 4x4", B, lambda d, B=B: am.mlts_amplification_batch(B, PS_4, [4, 4], np.array([d]))[0], y0)]
    dts = [0.4 / 2 ** k for k in range(6)]
    order_rows = []
    for model, name, B, C, y0 in cases:
        ex = expm(B * T) @ y0
        errs = []
        for d in dts:
            N = int(round(T / d))
            y = np.linalg.matrix_power(C(d), N) @ y0
            errs.append(np.linalg.norm(y - ex) / np.linalg.norm(ex))
        orders = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
        order_rows.append((model, name, errs, orders))
        print(f"  {model:7s} {name:14s} " + " ".join(f"{x:.2e}" for x in errs) + "  orders " + " ".join(f"{o:.2f}" for o in orders))
    res["order_dts"] = np.array(dts)
    res["order_errs"] = np.array([r[2] for r in order_rows])
    res["order_names"] = np.array([f"{r[0]}: {r[1]}" for r in order_rows])
    np.savez(os.path.join(a.out, "plan_extra.npz"), **res)

    # figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(EPS_LIST)))
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for k, key in enumerate(LI_MODELS):
        for ie, e in enumerate(EPS_LIST):
            ax[k].semilogx(P_GRID, res[f"li_{key}_dtc"][ie], "o-", ms=3, color=cols[ie], label=fr"$\varepsilon={e}$")
        ax[k].set_title("LI-SDIRK(3,4): " + LI_MODELS[key][2]); ax[k].set_xlabel("$p$")
        ax[k].set_ylim(0, 1.6); ax[k].axhline(1, color="k", lw=0.5)
    ax[0].set_ylabel(r"$\Delta t_c(p)/\Delta t_{\rm ref}$"); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "toy_li_dtc.png"), dpi=150)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for k, key in enumerate(MLTS_MODELS):
        for name, ls in (("alg4", "-"), ("am", "--")):
            for ie, e in enumerate(EPS_LIST):
                ax[k].semilogx(P3, res[f"mlts_{key}_{name}_dtc"][ie], "o" + ls, ms=3, color=cols[ie],
                               label=(fr"$\varepsilon={e}$" if name == "alg4" else None))
        ax[k].set_title("MLTS, " + MLTS_MODELS[key][2] + " (solid: Alg. 4, dashed: A-M)", fontsize=9)
        ax[k].set_xlabel("$p_1=p_2=p$"); ax[k].set_ylim(0, 1.6); ax[k].axhline(1, color="k", lw=0.5)
    ax[0].set_ylabel(r"$\Delta t_c/\Delta t_{\rm ref}$"); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "toy_mlts_alg4_am.png"), dpi=150)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for k, model in enumerate(("order2", "order1")):
        for (m, name, errs, _) in order_rows:
            if m == model: ax[k].loglog(dts, errs, "o-", ms=3, label=name)
        ax[k].loglog(dts, [errs[0] * (d / dts[0]) ** 4 for d in dts], "k:", label="slope 4")
        ax[k].set_title(f"{model} reduced model, error at T = {T:g}"); ax[k].set_xlabel(r"$\Delta t$")
        ax[k].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "toy_order.png"), dpi=150)
    print("figures:", a.out)


if __name__ == "__main__":
    main()
