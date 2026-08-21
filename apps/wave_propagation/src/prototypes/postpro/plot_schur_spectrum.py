#!/usr/bin/env python3
"""
Visualisation du spectre du complement de Schur S_SCHUR = Kcc - Kcf*Sff_inv*Kfc
(LHS_spectrum_acou.hpp, fichiers "level_*.txt") dans le plan
complexe, pour une serie de maillages grossiers.

S_SCHUR est de taille n_c x n_c (faces eliminees par condensation statique) :
c'est l'operateur cellule-cellule effectivement vu par les schemas explicites
de ce code, par opposition au bloc brut LHS.topLeftCorner(n_c,n_c) qui ignore
la contribution des faces.

Usage:
    python3 plot_schur_spectrum.py [dir]

    dir : repertoire "spectrum/condensed" contenant les fichiers level_*.txt
          (par defaut : repertoire courant)

Une image "level_<lvl>.png" est enregistree par maillage, et une figure
recapitulative "all_levels.png" affiche tous les niveaux cote a cote.
"""

import os
import re
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


def parse_spectrum_file(path):
    meta = {}
    reals, imags = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                for key in ("level", "nx", "ny", "k", "h_max", "n_c"):
                    m = re.search(rf"{key}=(-?[0-9.eE+\-]+)", line)
                    if m:
                        meta[key] = m.group(1)
                continue
            parts = line.split()
            if len(parts) >= 2:
                reals.append(float(parts[0]))
                imags.append(float(parts[1]))

    eigs = np.array(reals) + 1j * np.array(imags)
    return {
        "file": path,
        "level": int(meta.get("level", -1)),
        "nx": meta.get("nx", "?"),
        "ny": meta.get("ny", "?"),
        "k": meta.get("k", "?"),
        "h_max": meta.get("h_max", "?"),
        "n_c": meta.get("n_c", "?"),
        "eigs": eigs,
    }


def load_spectra(directory):
    files = sorted(glob.glob(os.path.join(directory, "level_*.txt")))
    spectra = [parse_spectrum_file(f) for f in files]
    spectra.sort(key=lambda s: s["level"])
    return spectra


def plot_spectrum(ax, spectrum, color="tab:purple"):
    ax.clear()
    eigs = spectrum["eigs"]
    rho = float(np.max(np.abs(eigs))) if eigs.size else float("nan")

    ax.scatter(eigs.real, eigs.imag, s=30, color=color, edgecolors="k",
               linewidths=0.3, zorder=3, label=rf"vap($S_{{SCHUR}}$), $\rho$={rho:.4g}")

    h_max = spectrum["h_max"]
    h_str = f"{float(h_max):.4g}" if h_max != "?" else "?"
    ax.set_title(
        f"Niveau {spectrum['level']}  (nx={spectrum['nx']}, ny={spectrum['ny']}, k={spectrum['k']})\n"
        f"n_c={spectrum['n_c']}  h_max={h_str}",
        fontsize=9)
    ax.set_xlabel(r"Re($\lambda$)")
    ax.set_ylabel(r"Im($\lambda$)")
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linewidth=0.3)
    ax.legend(fontsize=8, loc="best")
    return rho


def save_each_level(spectra, out_dir="."):
    for s in spectra:
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_spectrum(ax, s)
        out_name = os.path.join(out_dir, f"level_{s['level']}.png")
        fig.savefig(out_name, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  niveau {s['level']}   ->   {out_name}")


def save_all_levels(spectra, out_name="all_levels.png"):
    n = len(spectra)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 6), squeeze=False)
    axes = axes[0]

    for ax, s in zip(axes, spectra):
        plot_spectrum(ax, s)

    fig.suptitle(r"Spectre du complement de Schur $S_{SCHUR} = K_{cc} - K_{cf} S_{ff}^{-1} K_{fc}$",
                 fontsize=11)
    fig.savefig(out_name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure recapitulative -> {out_name}")


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    spectra = load_spectra(directory)

    if not spectra:
        print(f"Aucun fichier level_*.txt trouve dans '{directory}'.")
        sys.exit(1)

    print(f"{len(spectra)} maillage(s) charge(s) depuis '{directory}':")
    for s in spectra:
        rho = float(np.max(np.abs(s["eigs"]))) if s["eigs"].size else float("nan")
        print(f"  niveau {s['level']}   n_c={s['n_c']}  h_max={s['h_max']}  rho(S_SCHUR)={rho:.6g}")

    print("\nExport d'une image par maillage :")
    save_each_level(spectra, out_dir=directory)

    save_all_levels(spectra, out_name=os.path.join(directory, "all_levels.png"))


if __name__ == "__main__":
    main()
