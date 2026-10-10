#!/usr/bin/env python3
#
# Plots of the numerical results of Section 5 of "Unfitted hybrid high-order
# methods stabilized by polynomial extension for elliptic interface problems",
# from the outputs of run_paper.sh (the results obtained with cell
# agglomeration, in the directories *_agglo, are plotted with dashed lines).
#
# Usage: python3 plot_paper_results.py [results_dir]

import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "paper_results"
COLORS = {0: "darkgreen", 1: "orange", 2: "darkred", 3: "darkblue"}
MARKERS = {0: "^", 1: "s", 2: "o", 3: "D"}


def parse_errors(path):
    """Return {k: {"h": [...], "H1": [...], "L2": [...], "cond": [...]}}"""
    res, k = {}, None
    if not os.path.exists(path):
        return res
    for line in open(path):
        line = line.strip()
        m = re.match(r"Polynomial degree k : (\d+)", line)
        if m:
            k = int(m.group(1))
            res[k] = {"h": [], "H1": [], "L2": [], "cond": []}
            continue
        for key, pat in (("h", r"Characteristic h size = (\S+)"),
                         ("H1", r"H1-norm error = (\S+)"),
                         ("L2", r"L2-norm error = (\S+)"),
                         ("cond", r"condition number: (\S+)")):
            m = re.match(pat, line)
            if m and k is not None:
                res[k][key].append(float(m.group(1)))
    return res


def errors(name):
    return parse_errors(os.path.join(RES, name, "solution_error_file.txt"))


def exists(name):
    return os.path.exists(os.path.join(RES, name, "solution_error_file.txt"))


def slope_triangle(ax, h0, e0, order, color="k"):
    h1 = h0 * 2.0
    ax.plot([h0, h1, h1, h0], [e0, e0, e0 * 2.0**order, e0], color=color, lw=0.8)


def loglog_errors(ax, name, ks=range(4), ls="-", label_suffix=""):
    res = errors(name)
    for k in ks:
        if k not in res or not res[k]["h"]:
            continue
        ax.loglog(res[k]["h"], res[k]["H1"], ls, color=COLORS[k], marker=MARKERS[k], ms=4,
                  label=("k = %d" % k) + label_suffix)


def finish(ax, xlabel="h", ylabel="H1-error"):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls=":", lw=0.5)
    ax.legend(fontsize=7)


def read_cells(path):
    cells, itfs = {}, []
    for line in open(path):
        t = line.split()
        if t[0] == "cell":
            cells[int(t[1])] = [float(x) for x in t[2:6]] + [int(x) for x in t[6:10]]
        elif t[0] == "interface":
            v = [float(x) for x in t[2:]]
            itfs.append((v[0::2], v[1::2]))
    return cells, itfs


def fig5():
    fig, axs = plt.subplots(2, 2, figsize=(9, 9))
    for row, shape in enumerate(("circle", "flower")):
        for col, (suffix, title) in enumerate((("", "polynomial extension"), ("_agglo", "cell agglomeration"))):
            ax = axs[row][col]
            path = os.path.join(RES, "fig5_%s%s" % (shape, suffix), "cells_l0.txt")
            if not os.path.exists(path):
                ax.set_visible(False)
                continue
            cells, itfs = read_cells(path)
            for i, (x0, y0, x1, y1, loc, cut, stab, group) in cells.items():
                color = "white"
                if suffix == "":
                    if cut == 1:
                        color = "yellowgreen"   # ill-cut on the negative side (T^{KO,1})
                    elif cut == 2:
                        color = "slateblue"     # ill-cut on the positive side (T^{KO,2})
                elif group >= 0:
                    color = "dimgray"           # agglomerated cells
                ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc=color, ec="k", lw=0.5))
            if suffix == "":
                for i, (x0, y0, x1, y1, loc, cut, stab, group) in cells.items():
                    if stab >= 0:
                        s = cells[stab]
                        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
                        sx, sy = 0.5 * (s[0] + s[2]), 0.5 * (s[1] + s[3])
                        c = "darkgreen" if cut == 1 else "darkblue"
                        ax.annotate("", xy=(cx + 0.4 * (sx - cx), cy + 0.4 * (sy - cy)), xytext=(cx, cy),
                                    arrowprops=dict(arrowstyle="->", color=c, lw=0.8))
            for xs, ys in itfs:
                ax.plot(xs, ys, "r-", lw=1)
            ax.set_aspect("equal")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_title("Fig. 5: %s, %s" % (shape, title), fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig5.png"), dpi=150)
    plt.close(fig)


def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, name, title in ((axs[0], "fig6_circle", "circle"), (axs[1], "fig6_flower", "flower")):
        loglog_errors(ax, name)
        if exists(name + "_agglo"):
            loglog_errors(ax, name + "_agglo", ls="--", label_suffix=" (agglo)")
        ax.set_title("Fig. 6: solution (56), %s interface" % title)
        finish(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig6.png"), dpi=150)
    plt.close(fig)


def fig6_etaN20():
    """Fig. 6 with eta_N = 20 (value given in Section 5.1 of the paper)."""
    if not exists("fig6_circle_etaN20"):
        return
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, name, title in ((axs[0], "fig6_circle", "circle"), (axs[1], "fig6_flower", "flower")):
        loglog_errors(ax, name + "_etaN20", label_suffix=" (eta_N = 20)")
        if exists(name + "_agglo"):
            loglog_errors(ax, name + "_agglo", ls="--", label_suffix=" (agglo)")
        ax.set_title("Fig. 6 with eta_N = 20: %s interface" % title)
        finish(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig6_etaN20.png"), dpi=150)
    plt.close(fig)


def fig7():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    cols = {"0": "darkgreen", "0.1": "orange", "0.2": "darkred", "0.3": "darkblue"}
    for ax, suffix, title in ((axs[0], "", "polynomial extension"), (axs[1], "_agglo", "cell agglomeration")):
        for a, col in cols.items():
            for center, ls in (("subcell", "-"), ("cell", "--")):
                res = errors("fig7_a%s_%s%s" % (a, center, suffix))
                if 3 in res:
                    ax.loglog(res[3]["h"], res[3]["H1"], ls, color=col, marker="o", ms=3,
                              label="theta = %s (%s-centered bases)" % (a, center))
        ax.set_title("Fig. 7: k = 3, %s" % title)
        finish(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig7.png"), dpi=150)
    plt.close(fig)


def fig8():
    names = [("fig8", "polynomial extension"), ("fig8_agglo", "cell agglomeration")]
    names = [n for n in names if os.path.exists(os.path.join(RES, n[0], "sparsity_k0_l0.txt"))]
    if not names:
        return
    fig, axs = plt.subplots(1, len(names), figsize=(5 * len(names), 5), squeeze=False)
    for ax, (name, title) in zip(axs[0], names):
        path = os.path.join(RES, name, "sparsity_k0_l0.txt")
        ncells = 0
        loc_i, loc_j, ext_i, ext_j = [], [], [], []
        for line in open(path):
            if line.startswith("#"):
                ncells = int(line.split()[2])
                continue
            i, j, flag = (int(x) for x in line.split())
            if flag:
                ext_i.append(i); ext_j.append(j)
            else:
                loc_i.append(i); loc_j.append(j)
        ax.plot(loc_j, loc_i, "s", color="navy", ms=1)
        ax.plot(ext_j, ext_i, "s", color="darkred", ms=1)
        ax.axvline(ncells - 0.5, color="k", ls="--", lw=0.8)
        ax.axhline(ncells - 0.5, color="k", ls="--", lw=0.8)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_title("Fig. 8: %s" % title)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig8.png"), dpi=150)
    plt.close(fig)


def fig9():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    ax = axs[0]
    for k in range(4):
        xs, ys = [], []
        for m in range(5):
            res = errors("fig9_left_m%d" % m)
            if k in res and res[k]["H1"]:
                xs.append(10.0**m)
                ys.append(res[k]["H1"][-1])
        if xs:
            ax.loglog(xs, ys, "-", color=COLORS[k], marker=MARKERS[k], ms=4, label="k = %d" % k)
    ax.set_title("Fig. 9 left: solution (57), finest mesh")
    finish(ax, xlabel="contrast kappa_2/kappa_1")
    ax = axs[1]
    rcols = {2: "darkgreen", 4: "orange", 6: "darkred", 8: "darkblue"}
    lss = {0: "-", 1: "--", 2: "-.", 3: ":", 4: (0, (3, 1, 1, 1))}
    for r, col in rcols.items():
        for m, ls in lss.items():
            res = errors("fig9_right_r%d_m%d" % (r, m))
            if 3 in res:
                ax.loglog(res[3]["h"], res[3]["H1"], ls=ls, color=col, marker="o", ms=3,
                          label="r = %d" % r if m == 0 else None)
    ax.set_title("Fig. 9 right: k = 3, kappa_2 = 10^m kappa_1, m = 0..4")
    finish(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig9.png"), dpi=150)
    plt.close(fig)


def fig10():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, name, title in ((axs[0], "fig10_left", "(58a)"), (axs[1], "fig10_right", "(58b)")):
        loglog_errors(ax, name)
        ax.set_title("Fig. 10: solution %s, kappa_2 = 10^4 kappa_1" % title)
        finish(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig10.png"), dpi=150)
    plt.close(fig)


def fig11():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    loglog_errors(axs[0], "fig11_left")
    axs[0].set_title("Fig. 11 left: solution (59), r = 10")
    finish(axs[0])
    rcols = {4: "black", 6: "orange", 8: "darkred", 10: "darkblue", 11: "gray"}
    for r, col in rcols.items():
        res = errors("fig11_right_r%d" % r)
        if 3 in res:
            axs[1].loglog(res[3]["h"], res[3]["H1"], "-", color=col, marker="o", ms=3, label="r = %d" % r)
    axs[1].set_title("Fig. 11 right: solution (59), k = 3")
    finish(axs[1])
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig11.png"), dpi=150)
    plt.close(fig)


def fig12():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    for suffix, ls in (("", "-"), ("_agglo", "--")):
        for k in range(4):
            xs, ys = [], []
            for i in range(-4, 5):
                res = errors("fig12_left_i%d%s" % (i, suffix))
                if k in res and res[k]["cond"]:
                    xs.append(1.0 / 3.0 + i / 32.0)
                    ys.append(res[k]["cond"][0])
            if xs:
                axs[0].semilogy(xs, ys, ls, color=COLORS[k], marker=MARKERS[k], ms=4,
                                label="k = %d%s" % (k, " (agglo)" if suffix else ""))
            xs, ys = [], []
            for p in range(2, 10):
                res = errors("fig12_right_p%d%s" % (p, suffix))
                if k in res and res[k]["cond"]:
                    xs.append(0.5 * 10.0**(-p))
                    ys.append(res[k]["cond"][0])
            if xs:
                axs[1].loglog(xs, ys, ls, color=COLORS[k], marker=MARKERS[k], ms=4,
                              label="k = %d%s" % (k, " (agglo)" if suffix else ""))
    axs[0].set_title("Fig. 12 left: circle of radius R")
    finish(axs[0], xlabel="R", ylabel="condition number")
    axs[1].set_title("Fig. 12 right: square interface")
    finish(axs[1], xlabel="delta", ylabel="condition number")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig12.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig5, fig6, fig6_etaN20, fig7, fig8, fig9, fig10, fig11, fig12):
        f()
    print("Figures written in", RES)
