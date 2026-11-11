#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check of the cost model (eq. gains) of the LTS_HHO preprint on the two meshes of Fig. cost_examples
(make_cost_examples.py): measured gains over ERK4 of MLTS with two levels, of MLTS with all its levels and of LI, in the
settings (1E) and (2E), against the model, and the cost rho of a fine solve of LI in products of each setting.
Runs: article/jobs.txt, cost/ (one thread, timed alone): each scheme for n and 2n coarse steps (dt0 = 1/400, ERK4 with
dt0 / P), so that the time of n coarse steps is the difference of the two wall times (the setup and the factorization of
LI cancel); best of the repetitions given. rho is read from LI: t_LI / t_ERK4 = (4 N + 3 (1 + rho) |S|) / (4 P N).
Usage: python3 analyze_cost_check.py <cost> [<cost_r2> <cost_r3>]
"""
import os, re, sys
sys.argv, ARGS = sys.argv[:1], sys.argv[1:]
from make_cost_examples import EXAMPLES, RHO, build, gains

TAG = {0: "a", 1: "b"}


def walls(run):
    out = []
    for d in ARGS:
        f = os.path.join(d, run, "time.txt")
        if os.path.exists(f):
            m = re.findall(r"([0-9.]+) s ([0-9]+) KB", open(f).read())
            if m:
                out.append(float(m[-1][0]))
    return out


def n_coarse(run):
    for d in ARGS:
        f = os.path.join(d, run, "case.lua")
        if os.path.exists(f):
            return round(float(re.findall(r"time.final_time = ([0-9.e-]+)", open(f).read())[-1]) * 400)


for k, E in enumerate(EXAMPLES):
    m = TAG[k]
    leaves = build(E["depth"], E["features"])
    G, N, L, band = gains(leaves, RHO)
    P = 2 ** (L - 1)
    names = {"mlts2": "MLTS, 2 levels", "mlts": f"MLTS, {L} levels", "li": "LI"}
    for s, key, kring in (("o1", "1E", 1), ("o2", "2E", 2)):
        t = {}
        for sc in ("erk", "mlts", "mlts2", "li"):
            a, b = walls(f"{s}_{m}_{sc}_n"), walls(f"{s}_{m}_{sc}_2n")
            t[sc] = (min(b) - min(a)) / n_coarse(f"{s}_{m}_{sc}_n") if a and b else float("nan")
        # |S|: the gains of the model give it back from LI (G_LI = P N / (N + 3/4 (1 + rho) |S|))
        S = (P * N / G[key]["LI"] - N) / (0.75 * (1 + RHO[key]))
        rho = (t["li"] / t["erk"] * 4 * P * N - 4 * N) / (3 * S) - 1
        print(f"({m}) setting ({key}): ERK4 {t['erk'] * 1e3:.2f} ms per coarse step, |S| = {S:.0f} cells, "
              f"rho of LI = {rho:.1f} (figure: {RHO[key]})")
        for sc, nm in names.items():
            print(f"    {nm:16} measured {t['erk'] / t[sc]:7.2f}   model {G[key][nm]:7.2f}")
