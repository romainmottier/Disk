#!/usr/bin/env python3
"""
Analysis of the campaign run_pardiso.sh (2026-11-15) of the LTS_HHO preprint:
  1. L-shaped domain, one thread: CPU time of each run of Figs. ls, ls2, ls-form, written to <campaign>/times.json for
     make_article.py --times: best of three of the cpu column (construction and integration) for the full runs; for the
     longest MLTS runs, two runs of n and 2n macro steps (best of three each), cpu(N) = setup + N t with
     t = (cpu(2n) - cpu(n)) / n and setup = 2 cpu(n) - cpu(2n); checked on ls/mlts_l3 (1E) and ls2/mltsp_l2 (2E), which
     have both.
  2. <campaign>/view_article: the runs of Figs. ls, ls2, ls-form for make_article.py (errors from their c.dat): the
     full runs of the campaign, and the full runs of run_article.sh (<article>) for the longest MLTS runs (same code
     path: same errors).
  3. <campaign>/view_head: the layout of make_tfus.py --short (fit/, t1/, vor/: timed 8 us runs as <run>, <run>_r2,
     <run>_r3; full LI runs of out_full; MLTS signal fit/plate_mlts_w and staircase runs tfus/ of <article>).
  4. Head model, six threads: LI over 8 us extrapolated (factorization counted once) on the quadtree and Voronoi
     meshes, nonzeros of the factors and profile of LI.
Usage: python3 analyze_pardiso.py <campaign> [<article>]
  <campaign>: output directory of run_pardiso.sh (out, out_r2, out_r3, out_full, out_mem)
  <article>:  output directory of run_article.sh (published full runs: ls/mlts_l4-5, ls2/mltsp_l3-5,
              fit/plate_mlts_w, tfus/)
Then (figures of the preprint, in <figures>):
  python3 ../../make_article.py <campaign>/view_article <figures> --times <campaign>/times.json --memory <campaign>/out_mem
      (lshape.pdf, lshape2.pdf, lshape_form.pdf)
  python3 ../../make_cost_examples.py <figures> --runs <campaign>/out/cost <campaign>/out_r2/cost <campaign>/out_r3/cost \\
      --formulations --memory <campaign>/out_mem      (cost_formulations.pdf, cost_gains.pdf, cost_meshes.pdf)
  python3 ../../make_tfus.py <campaign>/view_head <figures> --memory <campaign>/out_mem --short --three-snapshots \\
      --skull-strip --focus-row --signals-row
      (tfus_*.pdf; since 2026-11-16: --three-snapshots, Fig. tfus-snap at t = 4, 20 and 44 us; --skull-strip, Fig.
      tfus-skull as strips in the frame of the plate, shear stress on a common scale; --focus-row, Fig. tfus-focus in
      one row with the axial profile in dB from the geometric focus; --signals-row, Fig. tfus-signals with the accuracy
      check and the energy only, the costs being in Table tfus-cost)
"""
import json, os, re, sys

if len(sys.argv) < 2:
    sys.exit(__doc__)
W = os.path.abspath(sys.argv[1])
ART = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else None
REPS = ("out", "out_r2", "out_r3")


def text(d):
    f = os.path.join(d, "log.txt")
    return re.sub(r"\x1b\[[0-9;]*m", "", open(f).read()) if os.path.exists(f) else ""


def cpu(d):
    """cpu column of the last row of c.dat (construction and integration), None if the run is missing or failed"""
    f = os.path.join(d, "c.dat")
    if not (os.path.exists(f) and os.path.exists(os.path.join(d, "time.txt"))):
        return None
    L = [l.split() for l in open(f) if l.strip()]
    names = L[0][1:]
    return float(L[-1][names.index("cpu")])


def best(run, f=cpu):
    v = [f(os.path.join(W, r, run)) for r in REPS]
    v = [x for x in v if x is not None]
    return (min(v), len(v)) if v else (None, 0)


def steps(run):
    m = re.search(r"final time [0-9.eE+-]+, ([0-9]+) steps", text(os.path.join(W, "out", run)))
    return int(m.group(1)) if m else None


def lshape_times():
    T, notes = {}, []
    full = ["ls/erk_l1"] + [f"ls/li2_l{l}" for l in range(1, 6)] + [f"ls/mlts_l{l}" for l in range(1, 4)] + \
           [f"ls2/li4_l{l}" for l in range(1, 6)] + [f"ls2/mltsp_l{l}" for l in range(1, 3)]
    for r in full:
        t, n = best(r)
        if t is not None:
            T[r] = t
            notes.append(f"{r:16s} full     {t:10.2f} s  (best of {n})")
    N = {"ls/mlts_l3": 80, "ls/mlts_l4": 120, "ls/mlts_l5": 160,
         "ls2/mltsp_l2": 75, "ls2/mltsp_l3": 100, "ls2/mltsp_l4": 150, "ls2/mltsp_l5": 200}
    for r, Nf in N.items():
        a, na = best(r + "_n"); b, nb = best(r + "_2n")
        if a is None or b is None:
            continue
        n = steps(r + "_n")
        t = (b - a) / n
        setup = 2 * a - b
        X = setup + Nf * t
        key = r if r not in T else r + "_extrapolated"
        T[key] = X
        notes.append(f"{key:24s} n = {n:3d}: cpu(n) {a:9.2f} cpu(2n) {b:9.2f} (best of {na}, {nb}), step {t:8.3f} s, "
                     f"setup {setup:7.2f} s -> {X:10.2f} s")
    for r in ("ls/mlts_l3", "ls2/mltsp_l2"):
        if r in T and r + "_extrapolated" in T:
            notes.append(f"check {r}: full {T[r]:.2f} s, extrapolated {T[r + '_extrapolated']:.2f} s "
                         f"({100 * (T[r + '_extrapolated'] / T[r] - 1):+.1f} %)")
    json.dump({k: v for k, v in T.items() if not k.endswith("_extrapolated")}, open(os.path.join(W, "times.json"), "w"),
              indent=1, sort_keys=True)
    return notes


def link(src, dst):
    if not os.path.exists(src):
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(src, dst)


def view_article():
    V = os.path.join(W, "view_article")
    new = ["ls/erk_l1"] + [f"ls/li2_l{l}" for l in range(1, 6)] + [f"ls/mlts_l{l}" for l in range(1, 4)] + \
          [f"ls2/li4_l{l}" for l in range(1, 6)] + [f"ls2/mltsp_l{l}" for l in range(1, 3)]
    for r in new:
        link(os.path.join(W, "out", r), os.path.join(V, r))
    if ART:
        for r in ("ls/mlts_l4", "ls/mlts_l5", "ls2/mltsp_l3", "ls2/mltsp_l4", "ls2/mltsp_l5"):
            link(os.path.join(ART, r), os.path.join(V, r))
    return V


def view_head():
    V = os.path.join(W, "view_head")
    for kind, runs in (("fit", ("plate_erk_8us", "plate_mlts_8us", "plate_li_8us")), ("t1", ("erk_8us", "mlts_8us", "li_8us")),
                       ("vor", ("plate_li_8us",))):
        for r in runs:
            for rep, suf in zip(REPS, ("", "_r2", "_r3")):
                link(os.path.join(W, rep, kind, r), os.path.join(V, kind, r + suf))
    for kind, r in (("fit", "plate_li"), ("fit", "plate_li_dt2"), ("fit", "intact_li"), ("fit", "water_li"), ("vor", "plate_li"),
                    ("fit", "plate_li_early")):                    # plate_li_early: jobs_snap.txt (2026-11-16)
        link(os.path.join(W, "out_full", kind, r), os.path.join(V, kind, r))
    if ART:
        link(os.path.join(ART, "fit", "plate_mlts_w"), os.path.join(V, "fit", "plate_mlts_w"))   # MLTS signal, energy
        if os.path.isdir(os.path.join(ART, "tfus")):                                         # staircase (MLTS)
            for r in os.listdir(os.path.join(ART, "tfus")):
                link(os.path.join(ART, "tfus", r), os.path.join(V, "tfus", r))
    return V


def integration(d):
    m = re.search(r"integration wall time \(s\): ([0-9.eE+-]+)", text(d))
    return float(m.group(1)) if m else None


def li_profile(d):
    t = text(d)
    m = re.search(r"LI time \(s\): fine system ([0-9.eE+-]+), coarse coefficients ([0-9.eE+-]+), stage solves "
                  r"([0-9.eE+-]+), stage products and updates ([0-9.eE+-]+)", t)
    z = re.search(r"nnz\(L \+ U\) ([0-9]+)", t)
    return (tuple(float(x) for x in m.groups()) if m else None), (int(z.group(1)) if z else None)


def li_extrapolated(d, n, N):
    (p, _), T = li_profile(d), integration(d)
    return None if p is None or T is None else (T - p[0]) * N / n + p[0]


def head_li():
    out = []
    q = [li_extrapolated(os.path.join(W, r, "fit", "plate_li_8us"), 139, 1394) for r in REPS]
    v = [li_extrapolated(os.path.join(W, r, "vor", "plate_li_8us"), 129, 1290) for r in REPS]
    q = [x for x in q if x]; v = [x for x in v if x]
    if q and v:
        out.append(f"LI six threads, extrapolated: quadtree {min(q) / 60:.1f} min, Voronoi {min(v) / 60:.1f} min "
                   f"(best of {len(q)}, {len(v)})")
    for k, d in (("quadtree", os.path.join(W, "out", "fit", "plate_li_8us")), ("Voronoi", os.path.join(W, "out", "vor", "plate_li_8us"))):
        p, z = li_profile(d)
        out.append(f"{k}: nnz(L + U) {z}, LI profile over 8 us (fine system, coarse coefficients, stage solves, "
                   f"stage products and updates; s) {p}")
    return out


if __name__ == "__main__":
    for l in lshape_times():
        print(l)
    print("view:", view_article())
    print("view:", view_head())
    for l in head_li():
        print(l)
