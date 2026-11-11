#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reduced models (no spatial discretization) of Section 4 of the LTS_HHO preprint: a chain of four tiers, one per level,
the ratio p between the time steps (and the cell sizes) of two consecutive levels, the total ratio p^3, coupling eps
between neighbouring unknowns:
  order2        second-order formulation, z'' + A z = 0, stiffnesses p^(2l), l = 0..3 (8 x 8)
  order1_nobuf  first-order formulation, y' = [0 G; -G 0] y, gradients p^l (8 x 8)
  order1_buf    first-order formulation with overlap: each level l >= 1 holds an unknown of the size of level l-1
                (one ring of overlap) and one of its own size, chain 1 - 1 - p - p - p^2 - p^2 - p^3 (14 x 14)
Schemes: ERK4 with the global step, MLTS with four levels (Algorithm 2 of the preprint, window variant, ratios p, p,
p) and LI = ERK(4)-SDIRK(3,4) with the stage values of the coarse polynomial (time.li_coupling = "stage"), the fine
region being the levels >= 1. dt_ref(eps): critical step of ERK4 at p = 1. Results: rho(C_p(dt_ref)) and the
critical step dt_c(p) / dt_ref (largest step below which every step is stable, grid of 200 steps up to 4).
Output: reduced_mlts4_li.npz (keys <model>_<scheme>_rho, <model>_<scheme>_dtc, arrays eps x p, and p, eps).
Usage: python3 reduced_mlts4_li.py [--out DIR] [--pmax 16]   (engines of plan_extra.py)
"""
import argparse
import math
import os
import numpy as np

import plan_extra as pe          # li_amplification, erk4_matrix, rho, EPS_LIST

L = 4


def chain(diag, eps):
    n = len(diag)
    G = np.diag(np.array(diag, dtype=float))
    for i in range(n - 1):
        G[i, i + 1] = G[i + 1, i] = eps
    return G


def build_order2(eps, p):
    A = chain([float(p) ** (2 * l) for l in range(L)], eps)
    B = np.zeros((2 * L, 2 * L))
    B[:L, L:] = np.eye(L)
    B[L:, :L] = -A
    return B, list(range(L)) * 2


def build_order1(G, levels):
    n = G.shape[0]
    B = np.zeros((2 * n, 2 * n))
    B[:n, n:] = G
    B[n:, :n] = -G
    return B, levels * 2


def build_order1_nobuf(eps, p):
    return build_order1(chain([float(p) ** l for l in range(L)], eps), list(range(L)))


def build_order1_buf(eps, p):
    sizes, levels = [1.0], [0]
    for l in range(1, L):
        sizes += [float(p) ** (l - 1), float(p) ** l]
        levels += [l, l]
    return build_order1(chain(sizes, eps), levels)


MODELS = {"order2": build_order2, "order1_nobuf": build_order1_nobuf, "order1_buf": build_order1_buf}


def projectors(levels):
    """Toy notation of plan_extra: P_l selects the levels >= l, l = 1..L-1."""
    lv = np.array(levels)
    return [np.diag((lv >= l).astype(float)) for l in range(1, L)]


def taylor_shift(w, s):
    return [sum(s ** (q - j) / math.factorial(q - j) * w[q] for q in range(j, 4)) for j in range(4)]


def mlts_amplification(B, Ps, ps, dts):
    """plan_extra.mlts_alg4_amplification vectorized over the steps dts (same operations, batched)."""
    n = B.shape[0]
    dts = np.atleast_1d(np.asarray(dts, dtype=float))
    I = np.broadcast_to(np.eye(n), (len(dts), n, n)).copy()
    Lmax = len(Ps)
    proj = [np.eye(n)] + list(Ps) + [np.zeros((n, n))]
    Bf = B @ proj[Lmax]
    BQ = [B @ (proj[l] - proj[l + 1]) for l in range(Lmax)]

    def step(level, Y, h, w):
        if level == Lmax:
            pw = lambda t: w[0] + t * w[1] + t * t / 2 * w[2] + t ** 3 / 6 * w[3]
            k1 = pw(0.0) + Bf @ Y
            k2 = pw(h / 2) + Bf @ (Y + h / 2 * k1)
            k3 = pw(h / 2) + Bf @ (Y + h / 2 * k2)
            k4 = pw(h) + Bf @ (Y + h * k3)
            return Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        g = [Y]
        for _ in range(3):
            g.append(B @ g[-1])
        w = [w[i] + BQ[level] @ g[i] for i in range(4)]
        p = ps[level]
        hs = h / p
        for m in range(p):
            Y = step(level + 1, Y, hs, taylor_shift(w, m * hs))
        return Y

    with np.errstate(over="ignore", invalid="ignore"):
        zero = np.zeros_like(I)
        return step(0, I, dts[:, None, None], [zero] * 4)


def erk4_dtc(B, tol=1e-5):
    stable = lambda d: pe.rho(np.array([pe.erk4_matrix(B, d)]))[0] <= 1 + tol
    lo, hi = 0.0, 4.0
    if stable(hi):
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if stable(mid) else (lo, mid)
    return lo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=pe.HERE)
    ap.add_argument("--pmax", type=int, default=16)
    a = ap.parse_args()
    PS = list(range(1, a.pmax + 1))
    res = {"p": np.array(PS), "eps": np.array(pe.EPS_LIST)}
    # check: the batched engine against plan_extra.mlts_alg4_amplification
    B1, lv1 = build_order1_buf(0.1, 3)
    d = np.abs(mlts_amplification(B1, projectors(lv1), [3, 3, 3], [0.7])[0]
               - pe.mlts_alg4_amplification(B1, projectors(lv1), [3, 3, 3], [0.7])[0]).max()
    print(f"check: batched MLTS = plan_extra.mlts_alg4_amplification: {d:.2e}")
    res["check_mlts_batched"] = d
    for key, build in MODELS.items():
        out = {s: ([], []) for s in ("erk", "mlts4", "li")}
        for e in pe.EPS_LIST:
            B, lv = build(e, 1)
            dr = erk4_dtc(B)
            rows = {s: ([], []) for s in out}
            for p in PS:
                B, lv = build(e, p)
                Ps = projectors(lv)
                P = Ps[0]
                amps = {
                    "erk": lambda d: np.array([pe.erk4_matrix(B, x) for x in np.atleast_1d(d)]),
                    "mlts4": lambda d: mlts_amplification(B, Ps, [p] * (L - 1), d),
                    "li": lambda d: pe.li_amplification(B, P, d, coupling="stage"),
                }
                for s, amp in amps.items():
                    rows[s][0].append(pe.rho(amp([dr]))[0])
                    rows[s][1].append((erk4_dtc(B) if s == "erk" else pe.dt_crit(amp)) / dr)
            for s in out:
                out[s][0].append(rows[s][0]); out[s][1].append(rows[s][1])
            print(f"  {key} eps = {e}: dt_ref = {dr:.4f}, max rho MLTS {max(rows['mlts4'][0]):.6f}, "
                  f"LI {max(rows['li'][0]):.6f}; min dt_c/dt_ref MLTS {min(rows['mlts4'][1]):.3f}, "
                  f"LI {min(rows['li'][1]):.3f}", flush=True)
        for s in out:
            res[f"{key}_{s}_rho"], res[f"{key}_{s}_dtc"] = np.array(out[s][0]), np.array(out[s][1])
    np.savez(os.path.join(a.out, "reduced_mlts4_li.npz"), **res)
    print("saved", os.path.join(a.out, "reduced_mlts4_li.npz"))


if __name__ == "__main__":
    main()
