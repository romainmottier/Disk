#!/usr/bin/env python3
"""
Visualisation du spectre complet de la matrice de rigidite (LHS_spectrum_acou.hpp)
dans le plan complexe, pour une serie de maillages grossiers.

LHS est de taille n_dof = n_c + n_f, ordonnee [dofs cellule | dofs face] :
le bloc cellule-cellule (V_T, P_T) en haut a gauche a la taille n_c x n_c.

Pour chaque valeur propre, LHS_spectrum_acou.hpp calcule frac_cell in [0,1] :
la fraction de l'energie du vecteur propre portee par les dofs cellule
(frac_cell=1 -> vecteur propre entierement dans le bloc cellule-cellule ;
 frac_cell=0 -> entierement dans le bloc face).

Coloration discrete a partir de frac_cell (memes seuils que rho_cell/rho_face
dans lhs_scaling_summary.txt) :
  frac_cell > 0.9  -> "cellule"  (rouge)
  frac_cell < 0.1  -> "face"     (bleu)
  sinon            -> "mixte"    (gris)

Usage:
    python3 plot_lhs_spectrum.py [dir]

    dir : repertoire "spectrum/dense" contenant les fichiers level_*.txt
          (par defaut : repertoire courant)

Une image "level_<lvl>.png" est enregistree par maillage, et une figure
recapitulative "all_levels.png" affiche tous les niveaux cote a cote
(meme legende partagee).
"""

import os
import re
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


def parse_spectrum_file(path):
    """Parse un fichier level_*.txt -> dict avec meta + tableaux."""
    meta = {}
    reals, imags, fracs = [], [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                for key in ("level", "nx", "ny", "k", "n_dof", "n_c", "n_f"):
                    m = re.search(rf"{key}=(-?[0-9.eE+\-]+)", line)
                    if m:
                        meta[key] = m.group(1)
                continue
            parts = line.split()
            if len(parts) >= 3:
                reals.append(float(parts[0]))
                imags.append(float(parts[1]))
                fracs.append(float(parts[2]))

    eigs = np.array(reals) + 1j * np.array(imags)
    frac_cell = np.array(fracs)
    return {
        "file": path,
        "level": int(meta.get("level", -1)),
        "nx": meta.get("nx", "?"),
        "ny": meta.get("ny", "?"),
        "k": meta.get("k", "?"),
        "n_dof": meta.get("n_dof", "?"),
        "n_c": meta.get("n_c", "?"),
        "n_f": meta.get("n_f", "?"),
        "eigs": eigs,
        "frac_cell": frac_cell,
    }


def load_spectra(directory):
    files = sorted(glob.glob(os.path.join(directory, "level_*.txt")))
    spectra = [parse_spectrum_file(f) for f in files]
    spectra.sort(key=lambda s: s["level"])
    return spectra


CELL_THRESH = 0.9
FACE_THRESH = 0.1

CATEGORIES = [
    ("cellule (frac_cell > 0.9)", "tab:red"),
    ("mixte", "0.7"),
    ("face (frac_cell < 0.1)", "tab:blue"),
]


def classify(frac_cell):
    """Retourne 0 (cellule), 1 (mixte) ou 2 (face) par valeur propre."""
    cat = np.full(frac_cell.shape, 1, dtype=int)
    cat[frac_cell > CELL_THRESH] = 0
    cat[frac_cell < FACE_THRESH] = 2
    return cat


def plot_spectrum(ax, spectrum):
    ax.clear()
    eigs = spectrum["eigs"]
    frac_cell = spectrum["frac_cell"]
    cat = classify(frac_cell)

    counts = [int(np.sum(cat == k)) for k in range(3)]
    for k, (label, color) in enumerate(CATEGORIES):
        mask = cat == k
        if not np.any(mask):
            continue
        ax.scatter(eigs.real[mask], eigs.imag[mask], s=30, color=color,
                   edgecolors="k", linewidths=0.3, zorder=3 if k != 1 else 2,
                   label=f"{label}  (n={counts[k]})")

    ax.set_title(
        f"Niveau {spectrum['level']}  (nx={spectrum['nx']}, ny={spectrum['ny']}, k={spectrum['k']})\n"
        f"n_dof={spectrum['n_dof']}  n_c={spectrum['n_c']}  n_f={spectrum['n_f']}",
        fontsize=9)
    ax.set_xlabel(r"Re($\lambda$)")
    ax.set_ylabel(r"Im($\lambda$)")
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linewidth=0.3)
    ax.legend(fontsize=7, loc="best")


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

    fig.suptitle("Spectre de LHS — cellule (rouge) vs mixte (gris) vs face (bleu)",
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
        print(f"  niveau {s['level']}   n_dof={s['n_dof']}  n_c={s['n_c']}  n_f={s['n_f']}  rho(LHS)={rho:.6g}")

    print("\nExport d'une image par maillage :")
    save_each_level(spectra, out_dir=directory)

    save_all_levels(spectra, out_name=os.path.join(directory, "all_levels.png"))


if __name__ == "__main__":
    main()
