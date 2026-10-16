#!/usr/bin/env python3
"""E5: work per macro step (output.profile) of every composition of the depth D into level ratios 2^a, direct and
reuse, on graded meshes. For each mesh: all ratios 2 (D+1 levels), two levels (ratio 2^D), the best composition, the
composition chosen by the merging rule r_l / b_l < p_l m_{l+1} - m_l applied greedily from the finest level, and the
model gain 4 N 2^D / W against ERK4 (whose step is 2^D times smaller)."""
import glob, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import read_run

R = "runs/e5w"
data = {}
for d in sorted(glob.glob(f"{R}/*")):
    m = re.match(r"(.+)_o(\d)_([0-9-]+)(_r)?$", os.path.basename(d))
    if not m:
        continue
    r = read_run(d)
    if r is None or "W" not in r:
        continue
    mesh, ov, tag, reuse = m.group(1), int(m.group(2)), m.group(3), bool(m.group(4))
    data.setdefault((mesh, ov), {})[(tag, reuse)] = r

out = {}
for (mesh, ov), runs in sorted(data.items()):
    D = int(re.search(r"_D(\d)", mesh).group(1))
    any_r = next(iter(runs.values()))
    N = any_r["N"]
    res = {}
    for reuse in (False, True):
        cand = {tag: r["W"] for (tag, rr), r in runs.items() if rr == reuse}
        if not cand:
            continue
        all2 = "-".join(["1"] * D)
        best = min(cand, key=cand.get)
        res["reuse" if reuse else "direct"] = dict(
            all2=cand.get(all2), two=cand.get(str(D)), best=best, W_best=cand[best],
            G_all2=4 * N * 2 ** D / cand[all2] if all2 in cand else None,
            G_two=4 * N * 2 ** D / cand[str(D)] if str(D) in cand else None,
            G_best=4 * N * 2 ** D / cand[best], n=len(cand),
            ranking=sorted(cand.items(), key=lambda kv: kv[1])[:5])
    out[f"{mesh}_o{ov}"] = dict(N=N, D=D, cells=any_r.get("cells_per_level"), **res)
    for k in ("direct", "reuse"):
        if k in res:
            x = res[k]
            print(f"{mesh:12} o{ov} {k:6} N={int(N):7d}: G all-2 {x['G_all2']:7.2f}  two-level {x['G_two']:7.2f}  "
                  f"best {x['G_best']:7.2f} [{x['best']}]  (W_best/W_all2 = {x['W_best'] / x['all2']:.3f})")
json.dump(out, open("e5w.json", "w"), indent=1)
