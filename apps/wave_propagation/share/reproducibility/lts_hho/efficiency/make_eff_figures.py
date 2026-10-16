#!/usr/bin/env python3
"""Figures of efficiency_study.tex from the runs of the efficiency study (runs/e2, e3, e4, e4b, e5w, e5t)."""
import glob, json, math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import read_run

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Github/Preprints/LTS_HHO/figures_efficiency")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.5,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "mathtext.fontset": "cm", "savefig.dpi": 300})
GRID = dict(color="#c8c8c4", lw=0.4, ls=":")
PC = {2: "#2b6cb0", 4: "#2f855a", 8: "#c05621", 16: "#9b2c2c", 32: "#553c9a"}
SC = {"erk": ("ERK4", "#1a1a1a", "o"), "lts": ("LTS", "#dd6b20", "s"), "mlts": ("MLTS (Alg. 4)", "#b7202e", "^"),
      "mltsr": ("MLTS, chain reused", "#e07a85", "v"), "li": ("LI", "#2b6cb0", "D")}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("ok", name)


def style(ax, xl, yl, logx=True, logy=True):
    if logx: ax.set_xscale("log")
    if logy: ax.set_yscale("log")
    ax.set_xlabel(xl, labelpad=1); ax.set_ylabel(yl, labelpad=1)
    ax.grid(True, which="major", **GRID)
    ax.tick_params(which="major", length=2.5, width=0.5, pad=1.5)


# ------------------------------------------------------------------------------------------- E2 stability, overlap
def stability():
    ref = 1.0 / 137
    rows = {}
    for d in glob.glob("runs/e2/*"):
        m = re.match(r"(two|grad)_L(\d)_(erk|lts|mlts)(?:_o(\d))?$", os.path.basename(d))
        if not m or not os.path.exists(f"{d}/c.dat"):
            continue
        dt = float(open(f"{d}/c.dat").read().split("\n")[1].split()[2])
        rows[(m.group(1), m.group(3), m.group(4), int(m.group(2)))] = dt / ref
    P = np.array([2, 4, 8, 16, 32])
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.3), sharey=True)
    for ax, fam, title in ((axs[0], "two", "two sizes (abrupt transition)"), (axs[1], "grad", "graded 2:1 transition")):
        for o, ls in (("0", "-"), ("1", "--"), ("2", ":")):
            for s in ("lts", "mlts"):
                v = [rows.get((fam, s, o, L), np.nan) for L in range(1, 6)]
                lab, c, mk = SC[s]
                ax.plot(P, v, ls + mk, color=c, ms=3.2 if s == "lts" else 2.2, mfc=c if s == "lts" else "white",
                        label=f"{lab}, overlap {o}")
        if fam == "two":
            v = [rows.get(("two", "erk", None, L), np.nan) for L in range(1, 6)]
            ax.plot(P, v, "-o", color="#1a1a1a", ms=3, label="ERK4")
        ax.set_title(title, fontsize=7, pad=3)
        style(ax, r"$p=h_{\max}/h_{\min}$", r"$\Delta t_{\rm crit}/\Delta t_c$" if fam == "two" else "")
        ax.set_xscale("log", base=2)
    h, l = axs[0].get_legend_handles_labels()
    axs[1].legend(h, l, loc="lower left", fontsize=5.5, frameon=False, ncol=2)
    axs[1].set_ylim(axs[0].get_ylim())
    fig.tight_layout(pad=0.2, w_pad=0.8)
    save(fig, "stab_overlap.pdf")


# --------------------------------------------------------------------------------------------- E4 two levels
def two_levels():
    E = json.load(open("e4.json"))
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.35))
    s = np.linspace(0.002, 1.0, 400)
    for ax, th in ((axs[0], 1), (axs[1], 6)):
        for p in (2, 4, 8, 16):
            R = [r for r in E if r["threads"] == th and r["p"] == p and "Gs_lts" in r]
            ax.plot([r["sigma"] for r in R], [r["Gs_lts"] for r in R], "o", color=PC[p], ms=3, label=rf"$p={p}$")
            ax.plot(s, p / (1 + (p + 1) * s), "-", color=PC[p], lw=0.8)
            if th == 6:
                # fine region on one thread: products s_6 times more expensive (s_6: parallel speed-up of ERK4)
                t1 = {r["mesh"]: r["t_erk"] for r in E if r["threads"] == 1}
                s6 = np.median([t1[r["mesh"]] / r["t_erk"] for r in E if r["threads"] == 6 and r["mesh"] in t1])
                ax.plot(s, p / (1 + (p + 1) * s * s6), ":", color=PC[p], lw=0.8)
        ax.axhline(1, color="#6a6a66", lw=0.6, ls="--")
        style(ax, r"$\sigma=|S|/N$", "gain over ERK4 (per step)")
        ax.set_title(f"LTS, {th} thread" + ("s" if th > 1 else ""), fontsize=7, pad=3)
        ax.set_ylim(0.3, 25)
    axs[0].legend(loc="lower left", frameon=False)
    ax = axs[2]
    for key, (lab, c, mk) in (("lts", SC["lts"]), ("mlts", SC["mlts"]), ("mltsr", SC["mltsr"])):
        R = [r for r in E if r["threads"] == 1 and f"Gs_{key}" in r and f"G_{key}_model" in r]
        ax.plot([r[f"G_{key}_model"] for r in R], [r[f"Gs_{key}"] for r in R], mk, color=c, ms=3, label=lab, mfc="white" if key == "mltsr" else c)
    g = np.array([0.3, 30])
    ax.plot(g, g, "-", color="#6a6a66", lw=0.6)
    style(ax, "model gain (work only)", "measured gain (1 thread)")
    ax.set_xlim(0.3, 30); ax.set_ylim(0.3, 30)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout(pad=0.2, w_pad=0.6)
    save(fig, "two_levels.pdf")


def li_vs_lts():
    E = json.load(open("e4.json"))
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.35))
    ax = axs[0]
    for th, mk in ((1, "o"), (6, "s")):
        R = [r for r in E if r["threads"] == th and "rho" in r]
        ax.plot([r["fine"] for r in R], [r["rho"] for r in R], mk, color="#2b6cb0" if th == 1 else "#90cdf4", ms=3,
                label=f"{th} thread" + ("s" if th > 1 else ""))
    style(ax, "fine unknowns", r"$\rho=\tau_{\rm LU}/(\kappa|S|)$")
    ax.legend(loc="upper left", frameon=False)
    ax = axs[1]
    for th, mk in ((1, "o"), (6, "s")):
        R = [r for r in E if r["threads"] == th and "Gs_li" in r and "Gs_lts" in r]
        x = [r["p"] / (0.75 * (1 + r["rho"])) for r in R]
        y = [r["Gs_li"] / r["Gs_lts"] for r in R]
        ax.plot(x, y, mk, color="#2b6cb0" if th == 1 else "#90cdf4", ms=3, label=f"{th} thread" + ("s" if th > 1 else ""))
    ax.axhline(1, color="#6a6a66", lw=0.6, ls="--"); ax.axvline(1, color="#6a6a66", lw=0.6, ls="--")
    style(ax, r"$p/p^\star$, $p^\star=0.75\,(1+\rho)$", r"time LTS / time LI")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout(pad=0.2, w_pad=1.0)
    save(fig, "li_vs_lts.pdf")


# ------------------------------------------------------------------------------------------------ E5 multilevel
def multilevel():
    W = json.load(open("e5w.json"))
    names = [k for k in W if k.endswith("_o0")]
    order = sorted(names, key=lambda k: (k.split("_")[0] != "pt", k.split("_")[0], k))
    fig, ax = plt.subplots(figsize=(6.3, 2.3))
    x = np.arange(len(order))
    for j, (lab, key, c) in enumerate((("all ratios 2", "G_all2", "#e07a85"), ("two levels", "G_two", "#dd6b20"),
                                       ("best composition", "G_best", "#b7202e"))):
        ax.bar(x + (j - 1) * 0.27, [W[k]["reuse"][key] for k in order], width=0.26, color=c, label=lab,
               edgecolor="#1a1a1a", linewidth=0.3)
    ax.set_xticks(x); ax.set_xticklabels([k.replace("_o0", "").replace("_", " ") for k in order], rotation=40, fontsize=5.5,
                                         ha="right")
    ax.set_yscale("log"); ax.axhline(1, color="#6a6a66", lw=0.6, ls="--")
    ax.set_ylabel("model gain over ERK4", labelpad=1); ax.grid(True, axis="y", **GRID)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout(pad=0.2)
    save(fig, "multilevel_work.pdf")


def multilevel_measured():
    T = json.load(open("e5t.json"))
    O = json.load(open("overhead.json"))
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.4), gridspec_kw=dict(width_ratios=[1.6, 1]))
    meshes = ["pt_W1_D6", "pt_W1_D8", "pt_W2_D8", "pt_W4_D8", "line_W1_D5", "disk_W1_D4"]
    labels = {"pt_W1_D6": "point, W=1\n$p=64$", "pt_W1_D8": "point, W=1\n$p=256$", "pt_W2_D8": "point, W=2\n$p=256$",
              "pt_W4_D8": "point, W=4\n$p=256$", "line_W1_D5": "curve, W=1\n$p=32$", "disk_W1_D4": "area, W=1\n$p=16$"}
    ax = axs[0]
    x = np.arange(len(meshes))
    keys = (("all2", "MLTS, ratios 2", "#e07a85"), ("bestr", "MLTS, best levels", "#b7202e"), ("lts", "two levels (LTS)", "#dd6b20"),
            ("li", "LI", "#2b6cb0"))
    for j, (k, lab, c) in enumerate(keys):
        v = [next(r["G"][k] for r in T if r["mesh"] == m and r["threads"] == 1) for m in meshes]
        ax.bar(x + (j - 1.5) * 0.2, v, width=0.19, color=c, label=lab, edgecolor="#1a1a1a", linewidth=0.3)
    gid = [next(r["Gid"] for r in T if r["mesh"] == m and r["threads"] == 1) for m in meshes]
    ax.plot(x, gid, "_", color="#1a1a1a", ms=14, mew=1.2, label=r"ideal $G_{\rm id}$")
    ax.set_xticks(x); ax.set_xticklabels([labels[m] for m in meshes], fontsize=5.5)
    ax.set_yscale("log"); ax.axhline(1, color="#6a6a66", lw=0.6, ls="--"); ax.grid(True, axis="y", **GRID)
    ax.set_ylabel("measured gain over ERK4 (1 thread)", labelpad=1)
    ax.legend(loc="upper right", frameon=False, ncol=2, fontsize=5.5)
    ax = axs[1]
    for th, mk in ((1, "o"), (6, "s")):
        for fam, c in (("lts", "#dd6b20"), ("mlts", "#b7202e")):
            R = [r for r in O if r["th"] == th and ((r["s"] == "lts") == (fam == "lts"))]
            y = np.array([r["t"] / r["kappa"] for r in R]); W = np.array([r["W"] for r in R]); C = np.array([r["C"] for r in R])
            A = np.column_stack([W / y, C / y]); sol, *_ = np.linalg.lstsq(A, np.ones_like(y), rcond=None)
            pred = sol[0] * W + sol[1] * C
            G_meas = np.array([4 * r["N"] * r["ratio"] * r["kappa"] / r["t"] for r in R])
            G_mod = np.array([4 * r["N"] * r["ratio"] / pv for r, pv in zip(R, pred)])
            ax.plot(G_mod, G_meas, mk, color=c, ms=2.6, mfc=c if th == 1 else "white", mew=0.6,
                    label=f"{fam.upper()}, {th} thr.")
    g = np.array([0.1, 40]); ax.plot(g, g, "-", color="#6a6a66", lw=0.6)
    style(ax, r"model $\kappa(aW+n_0C)$", "measured gain")
    ax.set_xlim(0.1, 40); ax.set_ylim(0.1, 40)
    ax.legend(loc="upper left", frameon=False, fontsize=5.5)
    fig.tight_layout(pad=0.2, w_pad=0.8)
    save(fig, "multilevel_measured.pdf")


if __name__ == "__main__":
    for f in (stability, two_levels, li_vs_lts, multilevel, multilevel_measured):
        try:
            f()
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f.__name__, "failed:", repr(e))
