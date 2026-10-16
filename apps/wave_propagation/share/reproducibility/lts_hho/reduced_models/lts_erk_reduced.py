#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reduced models (no spatial discretization), complement of plan_extra.py for the stability section of the LTS_HHO
preprint: ERK4 with the global step and the two-level LTS-RK4 (comparaison/lts_rk4_comparaison.py) on the three models
of the LI study of plan_extra.py (order 2, 4x4; order 1 without overlap, 4x4; order 1 with overlap, 6x6), for the same
eps and p: rho(C(dt_ref)) and the critical step dt_c(p) / dt_ref (dt_ref: critical step of RK4 at p = 1); ERK4 with
the global step at dt_ref: rho and dt_c.
Results: reduced_lts_erk.npz (keys erk_<model>_dtc, erk_<model>_rho, lts_<model>_rho, lts_<model>_dtc).
Usage: python3 lts_erk_reduced.py [--out DIR]   (MODELE_FINAL=~/Github/LTS_RK4/modele_final by default)
"""
import argparse
import os
import numpy as np

import plan_extra as pe          # models, dt_crit, rho, erk4_matrix, EPS_LIST, P_GRID (and the modele_final engines)

two = pe.two


def erk4_dtc(B, tol=1e-5):
    """Critical step of ERK4 with the global step (bisection: the stability interval of RK4 on the spectrum of B,
    imaginary or with a small real part, is [0, dt_c])."""
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
    a = ap.parse_args()
    res = {}
    for key, (build, P, label) in pe.LI_MODELS.items():
        mod = dict(build_B=build, P=P)
        rr, dc, de, re = [], [], [], []
        for e in pe.EPS_LIST:
            dr = two.dt_ref(mod, e)
            rr.append([pe.rho(two.amplification_batch(build(e, p), P, p, [dr]))[0] for p in pe.P_GRID])
            dc.append([pe.dt_crit(lambda d, p=p: two.amplification_batch(build(e, p), P, p, d)) / dr for p in pe.P_GRID])
            de.append([erk4_dtc(build(e, p)) / dr for p in pe.P_GRID])
            re.append([pe.rho(np.array([pe.erk4_matrix(build(e, p), dr)]))[0] for p in pe.P_GRID])
        res[f"lts_{key}_rho"], res[f"lts_{key}_dtc"], res[f"erk_{key}_dtc"] = np.array(rr), np.array(dc), np.array(de)
        res[f"erk_{key}_rho"] = np.array(re)
        print(f"{label}: LTS max rho(dt_ref) = {np.max(rr):.6f}, dt_c/dt_ref in [{np.min(dc):.3f}, {np.max(dc):.3f}]; "
              f"ERK4 dt_c/dt_ref at p = {pe.P_GRID[-1]}: {np.array(de)[:, -1].min():.2e} .. {np.array(de)[:, -1].max():.2e}")
    np.savez(os.path.join(a.out, "reduced_lts_erk.npz"), **res)


if __name__ == "__main__":
    main()
