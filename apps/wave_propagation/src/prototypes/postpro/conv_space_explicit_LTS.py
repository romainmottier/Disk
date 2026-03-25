import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --- Répertoire de base contenant l3, l4, l5, l6 ---
BASE_DIR = "../build/apps/wave_propagation/conv_tests"   # <-- Modifier si nécessaire

# --- Paramètres ---
LEVELS = [3, 4, 5, 6]    # niveaux de raffinement spatial
P_RANGE = range(1, 6)    # p1 à p5

# Stockage : results[p] = {'h': [], 'l2': [], 'dg': []}
results = {p: {'h': [], 'l2': [], 'dg': []} for p in P_RANGE}

# --- Parcours des dossiers ---
for l in LEVELS:
    l_dir = os.path.join(BASE_DIR, f"l{l}")
    if not os.path.isdir(l_dir):
        print(f"⚠️  Dossier manquant : {l_dir}")
        continue

    h_level = None  # h lu une seule fois par niveau

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

            # h
            if "Characteristic h size" in line:
                m = re.search(r"Characteristic h size\s*=\s*([0-9.eE+\-]+)", line)
                if m:
                    h_local = float(m.group(1))

            # L2 errors (HHO) — deux blocs "Acoustic region :" et "Elastic region :"
            if line.strip() == "Acoustic region :" and l2_acou is None:
                if i + 1 < len(lines):
                    m = re.search(r"L2-norm error\s*=\s*([0-9.eE+\-]+)", lines[i + 1])
                    if m:
                        l2_acou = float(m.group(1))

            if line.strip() == "Elastic region :" and l2_elas is None:
                if i + 1 < len(lines):
                    m = re.search(r"L2-norm error\s*=\s*([0-9.eE+\-]+)", lines[i + 1])
                    if m:
                        l2_elas = float(m.group(1))

            # dG errors
            if line.strip().startswith("L2 errors of dG unknowns:"):
                for j in range(i + 1, min(i + 6, len(lines))):
                    if lines[j].strip().startswith("Acoustic region :"):
                        try:
                            dg_acou = float(lines[j].split(":")[1].strip())
                        except (IndexError, ValueError):
                            pass
                    if lines[j].strip().startswith("Elastic  region :"):
                        try:
                            dg_elas = float(lines[j].split(":")[1].strip())
                        except (IndexError, ValueError):
                            pass

            i += 1

        if h_level is None and h_local is not None:
            h_level = h_local

        # Somme acoustique + élastique
        if (h_level is not None and
                l2_acou is not None and l2_elas is not None and
                dg_acou is not None and dg_elas is not None):
            results[p]['h'].append(h_level)
            results[p]['l2'].append(l2_acou + l2_elas)
            results[p]['dg'].append(dg_acou + dg_elas)
        else:
            print(f"⚠️  Données manquantes dans {filepath}  "
                  f"(h={h_level}, L2_acou={l2_acou}, L2_elas={l2_elas}, "
                  f"dG_acou={dg_acou}, dG_elas={dg_elas})")

# --- Style identique au script original ---
styles = {1: {'color': 'darkgreen',  'marker': '^'},
          2: {'color': 'orange',     'marker': 's'},
          3: {'color': 'darkred',    'marker': 'o'},
          4: {'color': 'steelblue',  'marker': 'D'},
          5: {'color': 'purple',     'marker': 'v'}}

plt.figure(figsize=(8, 8))

for p in P_RANGE:
    d = results[p]
    if not d['h']:
        continue

    sorted_data = sorted(zip(d['h'], d['l2'], d['dg']))
    h_s, l2_s, dg_s = zip(*sorted_data)
    h_arr = np.array(h_s)

    plt.loglog(h_arr, l2_s,
               marker=styles[p]['marker'], linestyle='-',
               color=styles[p]['color'],
               label=f'k={p} (L2)')
    plt.loglog(h_arr, dg_s,
               marker=styles[p]['marker'], linestyle='--',
               color=styles[p]['color'], alpha=0.7,
               label=f'k={p} (dG)')

plt.xlabel("h")
plt.ylabel("L2 error / dG error")
plt.legend()
plt.grid(True, which='both')
plt.show()
