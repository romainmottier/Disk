#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reduced models, complement of plan_extra.py: the corrected LI-SDIRK(3,4) (coarse stage values of the SDIRK tableau,
time.li_coupling = "stage" of wave_propagation) on the three models of the LI study (order 2, 4x4; order 1 without
overlap, 4x4; order 1 with overlap, 6x6), for the same eps and p: rho(C(dt_ref)) and dt_c(p) / dt_ref; and the error at
T = 4 against expm(B T) y0 of both couplings (eps = 0.1, p = 8).
Results: reduced_li_stage.npz (keys li_stage_<model>_rho, li_stage_<model>_dtc, order_dts, order_errs, order_names).
Usage: python3 li_stage_reduced.py [--out DIR]
"""
import argparse
import os
import numpy as np

import plan_extra as pe

two = pe.two


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=pe.HERE)
    a = ap.parse_args()
    res = {}
    for key, (build, P, label) in pe.LI_MODELS.items():
        mod = dict(build_B=build, P=P)
        rr, dc = [], []
        for e in pe.EPS_LIST:
            dr = two.dt_ref(mod, e)
            rr.append([pe.rho(pe.li_amplification(build(e, p), P, [dr], "stage"))[0] for p in pe.P_GRID])
            dc.append([pe.dt_crit(lambda d, p=p: pe.li_amplification(build(e, p), P, d, "stage")) / dr for p in pe.P_GRID])
        res[f"li_stage_{key}_rho"], res[f"li_stage_{key}_dtc"] = np.array(rr), np.array(dc)
        print(f"{label}: LI (stage) max rho(dt_ref) = {np.max(rr):.6f}, dt_c/dt_ref in [{np.min(dc):.3f}, {np.max(dc):.3f}]")
    T, e = 4.0, 0.1
    rng = np.random.default_rng(1)
    dts = [0.4 / 2 ** k for k in range(6)]
    rows = []
    for model, build, P, n in (("order2", two.build_B_order2, two.P_ORDER2, 4), ("order1", two.build_B_order1, two.P_ORDER1, 6)):
        B = build(e, 8); y0 = rng.standard_normal(n)
        ex = pe.expm(B * T) @ y0
        for cp in ("taylor", "stage"):
            errs = []
            for d in dts:
                N = int(round(T / d))
                y = np.linalg.matrix_power(pe.li_amplification(B, P, [d], cp)[0], N) @ y0
                errs.append(np.linalg.norm(y - ex) / np.linalg.norm(ex))
            rows.append((f"{model}: LI ({cp}) p=8", errs))
            print(f"  {model} LI {cp:6s} " + " ".join(f"{x:.2e}" for x in errs) + "  orders "
                  + " ".join(f"{np.log2(errs[i] / errs[i + 1]):.2f}" for i in range(len(errs) - 1)))
    res["order_dts"] = np.array(dts)
    res["order_errs"] = np.array([r[1] for r in rows])
    res["order_names"] = np.array([r[0] for r in rows])
    np.savez(os.path.join(a.out, "reduced_li_stage.npz"), **res)


if __name__ == "__main__":
    main()
