import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

# --- LaTeX rendering ---
matplotlib.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 13,
    'legend.fontsize': 14,
})

# --- Répertoire de base ---
BASE_DIR = "../build/apps/wave_propagation/conv_tests"

# --- Paramètres ---
LEVELS  = [3, 4, 5, 6]
P_RANGE = range(1, 6)

results = {p: {'h': [], 'l2': [], 'dg': []} for p in P_RANGE}

# --- Lecture des fichiers ---
for l in LEVELS:
    l_dir = os.path.join(BASE_DIR, f"l{l}")
    if not os.path.isdir(l_dir):
        print(f"⚠️  Dossier manquant : {l_dir}")
        continue

    h_level = None

    for p in P_RANGE:
        p_dir = os.path.join(l_dir, f"p{p}")
        if not os.path.isdir(p_dir):
            print(f"⚠️  Dossier manquant : {p_dir}")
            continue

        txt_files = [f for f in os.listdir(p_dir)
                     if f.startswith("explicit_l_") and f.endswith(".txt")]
        if not txt_files:
            print(f"⚠️  Aucun fichier txt dans {p_dir}")
            continue

        filepath = os.path.join(p_dir, txt_files[0])
        h_local = None
        l2_acou = l2_elas = dg_acou = dg_elas = None

        with open(filepath, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i]
            if "Characteristic h size" in line:
                m = re.search(r"Characteristic h size\s*=\s*([0-9.eE+\-]+)", line)
                if m:
                    h_local = float(m.group(1))
            if line.strip() == "Acoustic region :" and l2_acou is None:
                if i + 1 < len(lines):
                    m = re.search(r"L2-norm error\s*=\s*([0-9.eE+\-]+)", lines[i+1])
                    if m:
                        l2_acou = float(m.group(1))
            if line.strip() == "Elastic region :" and l2_elas is None:
                if i + 1 < len(lines):
                    m = re.search(r"L2-norm error\s*=\s*([0-9.eE+\-]+)", lines[i+1])
                    if m:
                        l2_elas = float(m.group(1))
            if line.strip().startswith("L2 errors of dG unknowns:"):
                for j in range(i+1, min(i+6, len(lines))):
                    if lines[j].strip().startswith("Acoustic region :"):
                        try: dg_acou = float(lines[j].split(":")[1].strip())
                        except: pass
                    if lines[j].strip().startswith("Elastic  region :"):
                        try: dg_elas = float(lines[j].split(":")[1].strip())
                        except: pass
            i += 1

        if h_level is None and h_local is not None:
            h_level = h_local

        if (h_level is not None and
                l2_acou is not None and l2_elas is not None and
                dg_acou is not None and dg_elas is not None):
            results[p]['h'].append(h_level)
            results[p]['l2'].append(l2_acou + l2_elas)
            results[p]['dg'].append(dg_acou + dg_elas)
        else:
            print(f"⚠️  Données manquantes dans {filepath}")

# --- Couleurs (comme sur la figure de référence) ---
colors = {1: '#FFD700', 2: '#006400', 3: '#FF8C00', 4: '#8B0000', 5: '#000080'}

fig, ax = plt.subplots(figsize=(7, 7))

h_all = []
for p in P_RANGE:
    d = results[p]
    if not d['h']:
        continue
    sorted_data = sorted(zip(d['h'], d['l2'], d['dg']))
    h_s, l2_s, dg_s = zip(*sorted_data)
    h_arr = np.array(h_s)
    h_all.extend(h_s)

    # dG : trait plein + rond
    ax.loglog(h_arr, dg_s,
              marker='o', linestyle='-', color=colors[p],
              markersize=5, linewidth=1.6)
    # HHO (L2) : tirets + carré
    ax.loglog(h_arr, l2_s,
              marker='s', linestyle='--', color=colors[p],
              markersize=5, linewidth=1.6)

# --- Triangle de pente ---
# Positionné en bas à gauche, comme sur la figure de référence
h_min, h_max = min(h_all), max(h_all)

# On dessine les triangles pour p=1..5 empilés
# Base commune en x, hauteurs différentes selon la pente
x0 = h_min * 1.25
x1 = h_min * 1.8

# Récupérer les valeurs y minimales pour ancrer le triangle
y_mins = []
for p in P_RANGE:
    d = results[p]
    if d['dg']:
        y_mins.append(min(d['dg']))
y_anchor = min(y_mins) * 0.15  # un peu en dessous du minimum

for p in P_RANGE:
    order = 4
    y0 = y_anchor
    y1 = y0 * (x1 / x0) ** order
    # Triangle rempli blanc avec contour noir
    tri = Polygon([[x0, y0], [x1, y0], [x1, y1]], closed=True, facecolor='white', edgecolor='black', linewidth=0.8, zorder=5)
    ax.add_patch(tri)
    # Annotation de la pente
    ax.text(x1 * 1.04, np.sqrt(y0 * y1), str(order), fontsize=8, va='center', ha='left', zorder=6)

# --- Légende ---
# Bloc 1 : type de norme
legend_norms = [
    Line2D([0], [0], color='black', linestyle='-',  marker='o', markersize=5, linewidth=1.6, label=r'$\|\cdot\|_{\mathrm{dG}} \approx \mathcal{O}(h^{k+1})$'),
    Line2D([0], [0], color='black', linestyle='--', marker='s', markersize=5, linewidth=1.6, label=r'$\|\cdot\|_{\mathrm{HHO}} \approx \mathcal{O}(h^{k+1})$'),
]

# Bloc 2 : degrés
legend_degrees = [
    Line2D([0], [0], color=colors[p], linestyle='-', linewidth=2.2,
           label=rf'$p = 2^{p}$')
    for p in P_RANGE
]

leg1 = ax.legend(handles=legend_norms + legend_degrees,
                 loc='lower right',
                 framealpha=0.95,
                 edgecolor='grey',
                 handlelength=2.2,
                 labelspacing=0.4)

ax.add_artist(leg1)

# --- Axes ---
ax.set_xlabel(r'$h$')
ax.set_ylabel(r'$\textsc{L}^2\textsc{-ERROR}$', labelpad=8)
ax.xaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext())
ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext())
ax.grid(True, which='both', linestyle='--', alpha=0.35)
ax.tick_params(which='both', direction='in', top=True, right=True)

plt.tight_layout()
plt.savefig("../build/apps/wave_propagation/conv_tests/conv_space_explicit.png", dpi=200, bbox_inches='tight')
plt.show()
print("Figures sauvegardées : conv_space_explicit.pdf / .png")
