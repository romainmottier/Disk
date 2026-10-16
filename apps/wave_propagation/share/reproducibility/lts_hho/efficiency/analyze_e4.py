#!/usr/bin/env python3
"""E4: two-level schemes against ERK4 on two-size meshes. Measured gains (integration wall time of ERK4 at its step
over that of the scheme at the coarse step, same final time) against the cost models; unit costs kappa (ERK4),
kappa_S (fine sub-steps of LTS), tau_LU (LI solves). Writes e4.json for the figures."""
import glob, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import read_run

R = "runs/e4"
rows = {}
for d in sorted(glob.glob(f"{R}/*")):
    m = re.match(r"(.+)_t(\d+)_(erk|lts|mltsr|mlts|li)$", os.path.basename(d))
    if not m or not os.path.exists(f"{d}/time.txt"):
        continue
    r = read_run(d)
    if r is None or "tint" not in r:
        continue
    rows.setdefault((m.group(1), int(m.group(2))), {})[m.group(3)] = r

out = []
print(f"{'mesh':18} {'th':>2} {'p':>3} {'N':>8} {'theta':>6} {'sigma':>6} | {'G_LTS':>6} {'mod':>6} | {'G_ML':>6} {'mod':>6} |"
      f" {'G_MLr':>6} {'mod':>6} | {'G_LI':>6} {'rho':>6} | {'kS/k':>5}")
for (mesh, th), s in sorted(rows.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if "erk" not in s or "lts" not in s:
        continue
    e, l = s["erk"], s["lts"]
    p = round(l["h_max"] / l["h_min"])
    N, S, F = l["N"], l["S"], l["fine"]
    M = l["steps"]                                   # macro steps; ERK4: M p steps
    kappa = e["tint"] / (M * p * 4 * N)              # s per unknown-product, ERK4
    rec = dict(mesh=mesh, threads=th, p=p, N=N, S=S, fine=F, theta=F / N, sigma=S / N, M=M, kappa=kappa,
               t_erk=e["tint"])
    rec["G_lts"] = e["tint"] / l["tint"]
    rec["G_lts_model"] = p / (1 + (p + 1) * S / N)
    # efficiency of the fine sub-steps per unknown-product against ERK4
    rec["kS_over_k"] = l["t_fine"] / (M * 4 * p * S) / kappa if "t_fine" in l else None
    rec["kC_over_k"] = l["t_coarse"] / (M * (4 * N + 4 * S)) / kappa if "t_coarse" in l else None
    for key in ("mlts", "mltsr"):
        if key in s and "W" in s[key]:
            rec[f"G_{key}"] = e["tint"] / s[key]["tint"]
            rec[f"G_{key}_model"] = 4 * N * p / s[key]["W"]
            rec[f"W_{key}"] = s[key]["W"]
    if "li" in s and "t_solve" in s["li"]:
        li = s["li"]
        rec["G_li"] = e["tint"] / li["tint"]
        tau = li["t_solve"] / (3 * M)                # one solve
        rec["tau_lu"] = tau
        rec["rho"] = tau / (kappa * S)
        rec["nnz_lu"] = li.get("nnz_lu")
        rec["rows"] = li.get("rows")
        rec["t_fact"] = li["t_fact"]
        rec["G_li_model"] = p * 4 * N / (4 * N + 7 * S + 3 * rec["rho"] * S + li["t_fact"] / (M * kappa))
    # per-step gains (timers of the loops: the setup of the sets and the factorization excluded)
    if "t_coarse" in l: rec["Gs_lts"] = e["tint"] / (l["t_coarse"] + l["t_fine"])
    for key in ("mlts", "mltsr"):
        if key in s and "t_coef" in s[key]:
            m_ = s[key]; rec[f"Gs_{key}"] = e["tint"] / (sum(m_["t_coef"]) + m_["t_term"] + m_["t_close"])
            rec[f"setup_{key}"] = m_["tint"] - (sum(m_["t_coef"]) + m_["t_term"] + m_["t_close"])
    if "li" in s and "t_solve" in s["li"]:
        li = s["li"]; rec["Gs_li"] = e["tint"] / (li["t_coarse"] + li["t_solve"] + li["t_stage"])
    out.append(rec)
    f = lambda k: f"{rec[k]:6.2f}" if rec.get(k) is not None else f"{'-':>6}"
    print(f"{mesh:18} {th:2d} {p:3d} {int(N):8d} {rec['theta']:6.3f} {rec['sigma']:6.3f} | {f('G_lts')} {f('G_lts_model')} |"
          f" {f('G_mlts')} {f('G_mlts_model')} | {f('G_mltsr')} {f('G_mltsr_model')} | {f('G_li')} {f('rho')} |"
          f" {rec['kS_over_k'] if rec['kS_over_k'] is None else round(rec['kS_over_k'], 2):>5}")
json.dump(out, open("e4.json", "w"), indent=1)
print()
print(f"{'mesh':18} {'th':>2} {'p':>3} {'sigma':>6} | per step: {'LTS':>6} {'MLTS':>6} {'MLTSr':>6} {'LI':>6} | models: {'LTS':>6} {'MLTS':>6} {'MLTSr':>6} | setup MLTS")
for r in sorted(out, key=lambda r: (r['threads'], r['p'], r['sigma'])):
    g = lambda k: f"{r[k]:6.2f}" if r.get(k) is not None else f"{'-':>6}"
    print(f"{r['mesh']:18} {r['threads']:2d} {r['p']:3d} {r['sigma']:6.3f} |           {g('Gs_lts')} {g('Gs_mlts')} {g('Gs_mltsr')} {g('Gs_li')} |"
          f"         {g('G_lts_model')} {g('G_mlts_model')} {g('G_mltsr_model')} | {r.get('setup_mlts', 0):.2f} s")
