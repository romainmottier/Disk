#!/usr/bin/env python3
"""Parsing of the run logs of the efficiency study (output.profile)."""
import os, re

def num(s):
    return float(s)

def read_run(d):
    """dict: tint (integration wall time), N (unknowns), steps, W (work per macro step), erk_eq (W / 4N), S, fine,
    t_coarse, t_fine (LTS), t_coef (list), t_term, t_close (MLTS), cells per level (MLTS), dt, wall (time.txt)."""
    r = {"dir": d}
    log = os.path.join(d, "log.txt")
    if not os.path.exists(log):
        return None
    txt = open(log, errors="replace").read()
    txt = re.sub(r"\x1b\[[0-9;]*m", "", txt)
    m = re.search(r"integration wall time \(s\): ([0-9.eE+-]+)", txt)
    if m: r["tint"] = num(m.group(1))
    m = re.search(r"LTS work per macro step \(unknowns x products\): 4 x ([0-9.e+]+) \+ 4 x ([0-9.e+]+) \+ (\d+) x "
                  r"[0-9.e+]+ = ([0-9.e+]+) = ([0-9.e+-]+) ERK4 steps \(fine unknowns (\d+), region S", txt)
    if m:
        r.update(N=num(m.group(1)), S=num(m.group(2)), W=num(m.group(4)), erk_eq=num(m.group(5)), fine=num(m.group(6)))
    m = re.search(r"LTS time \(s\): coarse coefficients ([0-9.e+-]+), fine substeps ([0-9.e+-]+)", txt)
    if m: r.update(t_coarse=num(m.group(1)), t_fine=num(m.group(2)))
    m = re.search(r"LI work per macro step \(unknowns x products\): 4 x ([0-9.e+]+) \+ 4 x ([0-9.e+]+) \+ \d+ x "
                  r"[0-9.e+]+ = ([0-9.e+]+) = ([0-9.e+-]+) ERK4 steps, plus (\d+) solves \(fine unknowns (\d+), region S "
                  r"[0-9.e+]+, solved rows (\d+), nnz\(L \+ U\) (\d+)\)", txt)
    if m:
        r.update(N=num(m.group(1)), S=num(m.group(2)), W=num(m.group(3)), erk_eq=num(m.group(4)), fine=num(m.group(6)),
                 rows=num(m.group(7)), nnz_lu=num(m.group(8)))
    m = re.search(r"LI time \(s\): fine system ([0-9.e+-]+), coarse coefficients ([0-9.e+-]+), stage solves "
                  r"([0-9.e+-]+), stage products and updates ([0-9.e+-]+)", txt)
    if m: r.update(t_fact=num(m.group(1)), t_coarse=num(m.group(2)), t_solve=num(m.group(3)), t_stage=num(m.group(4)))
    m = re.search(r"MLTS work: ([0-9.e+]+) per macro step = ([0-9.e+-]+) ERK4 steps \((\d+) unknowns\), terminal "
                  r"([0-9.e+-]+) %", txt)
    if m: r.update(W=num(m.group(1)), erk_eq=num(m.group(2)), N=num(m.group(3)), term_pct=num(m.group(4)))
    m = re.search(r"MLTS work per macro step \(unknowns x products x windows\):(.*)", txt)
    if m:
        r["levels_work"] = re.findall(r"\[(\d+|terminal): band (\d+), ([0-9.e+]+) x ([0-9.e+]+)\]", m.group(1))
    m = re.search(r"MLTS time \(s\): coefficients([0-9.e+\- ]*), terminal ([0-9.e+-]+), closed form ([0-9.e+-]+)", txt)
    if m:
        r.update(t_coef=[num(x) for x in m.group(1).split()], t_term=num(m.group(2)), t_close=num(m.group(3)))
    m = re.search(r"MLTS: (\d+) levels, cells per level \(coarsest first\):([0-9 ]+)", txt)
    if m: r["cells_per_level"] = [int(x) for x in m.group(2).split()]
    m = re.search(r"nt = (\d+)", txt)
    if m: r["steps"] = int(m.group(1))
    m = re.search(r"h_max = ([0-9.e+-]+)\s+h_min = ([0-9.e+-]+)", txt)
    if m: r.update(h_max=num(m.group(1)), h_min=num(m.group(2)))
    m = re.search(r"fine cells = (\d+)/(\d+)", txt)
    if m: r.update(fine_cells=int(m.group(1)), cells=int(m.group(2)))
    dat = os.path.join(d, "c.dat")
    if os.path.exists(dat):
        lines = open(dat).read().split("\n")
        names = lines[0].lstrip("# ").split()
        vals = [float(x) for x in lines[1].split()]
        r.update(dict(zip(names, vals)))
    t = os.path.join(d, "time.txt")
    if os.path.exists(t):
        s = open(t).read().split()
        r["wall"] = float(s[0]) if s and s[0] != "failed" else None
    return r
