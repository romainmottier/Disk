#!/usr/bin/env python3
"""
plot_convergence.py — Figures et rapport des campagnes de reproduction des resultats de [CAMC22] et [JSC21] avec
wave_propagation (ordres observes : meme methode que plot_convergence.py du code P_Laplacien).

  [CAMC22] E. Burman, O. Duran, A. Ern, Hybrid high-order methods for the acoustic wave equation in
           the time domain, Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702
  [JSC21]  E. Burman, O. Duran, A. Ern, M. Steins, Convergence analysis of hybrid high-order methods
           for the wave equation, J. Sci. Comput. 87 (2021) 91, hal-02922720

Lit les fichiers de resultats JSON ecrits par le run_campaign.py d'un article dans son dossier de resultats
(figure_<n>/<cas>/<cas>.json : parametres, ordres attendus, erreurs a l'instant final, historique d'energie) et les
valeurs relevees sur les figures publiees (valeurs_article.json du dossier de l'article), puis genere, dans le
dossier figure_<n> de chaque figure n de l'article :
  - <ART>_figure_<n>.png : la figure reproduite, panneaux disposes comme dans l'article et dessines comme lui
    (erreurs L2 en tirets fins fonces, H1 en tirets epais clairs, memes couleurs, marqueurs et legendes ; pertes
    d'energie en echelle logarithmique) ; <ART>_fig<n>_<panneau>.png : chaque panneau seul ;
  - <ART>_figure_<n>_comparaison.png : la meme figure avec les valeurs relevees sur la figure publiee (cercles
    noirs) et le rapport nous/article de chaque panneau ;
  - variantes/ : les panneaux des cas « variante » (parametres de la legende quand ils different de ceux de la
    figure publiee, ...) ;
  - ordres_observes/ : une figure par panneau et par norme avec l'ordre observe (OCV) de chaque segment et la droite
    de l'ordre attendu ;
  - pour [CAMC22], figure_7/CAMC22_figure_7.png : les quatre premiers maillages PolyMesher de la Fig. 8 ;
et, a la racine du dossier de resultats, <article>_recapitulatif.png (toutes les figures reproduites) et rapport.md
(par figure : parametres, ecarts a la legende, rapports nous/article, ordres attendus et observes, statut).
PNG sur fond blanc par defaut ; --format pdf (ou pdf,png) pour des PDF a fond transparent (LaTeX).

Ordre observe retenu (comme le floor_idx de P_Laplacien) : moyenne des deux derniers OCV fiables. Le plancher final
(erreur en espace atteinte dans une etude en temps, erreur en temps ou arrondi dans une etude en espace) est exclu :
c'est la plus longue suite finale de segments dont l'OCV est < 0.8 x l'ordre attendu (aucune exclusion si tous les
segments le sont : convergence lente). Etude en temps : le premier segment (dt = 0.1, pre-asymptotique) est ignore.
Une fenetre explicite de segments peut etre donnee par cas. Statut OK si retenu >= attendu - tolerance pour un ordre
optimal annonce (une superconvergence n'est pas un echec), si |retenu - attendu| <= tolerance pour un ordre
sous-optimal annonce par l'article (tau = O(1)) ; tolerance 0.3, 0.5 en mode --quick de run_campaign.py. ECHEC
aussi si une valeur manque, n'est pas finie ou n'est pas positive (calcul instable, cas interrompu) ou si le tableau
n'a pas le nombre de lignes attendu.

Usage (depuis la copie du build de papers_Duran) :
    python3 plot_convergence.py acoustic_camc2022/
    python3 plot_convergence.py elastic_jsc2021/ --format png,pdf
"""
import argparse, collections, json, math, os, re, sys, time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
import numpy as np

HERE = Path(__file__).resolve().parent

_CMAP_NORM = {"L2": plt.cm.Blues, "H1": plt.cm.Reds}
_MARKERS = {0: "o", 1: "s", 2: "^", 3: "D", 4: "v", 5: "P"}

# Couleurs (Mathematica) des figures des articles : courbes L2 foncees, H1 claires ; k (ou s) = 0, 1, 2 -> bleu,
# rouge, orange ; marqueurs L2 (disque, carre, losange) et H1 (triangles, demi-disque)
_DARK = [(40, 96, 136), (112, 16, 8), (224, 128, 16)]
_LIGHT = [(144, 176, 192), (184, 136, 128), (240, 192, 136)]
_rgb = lambda c: tuple(v / 255.0 for v in c)
_L2_MARKERS = ["o", "s", "D"]
_H1_MARKERS = ["^", "v", "o"]

# Label d'axe (notation des articles) de chaque colonne d'erreur ecrite par les prototypes
_YLABEL = {
    "acoustic_phi":      r"$\Vert p-p_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
    "acoustic_grad_phi": r"$\Vert \nabla p-\nabla R_{\mathcal{T}}(\hat p_h)\Vert_{L^2(\Omega)}$",
    "acoustic_phi_t":    r"$\Vert \partial_t p-v_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
    "acoustic_p":        r"$\Vert v-v_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
    "acoustic_q":        r"$\Vert \sigma-\sigma_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
    "elastic_u":         r"$\Vert u-u_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
    "elastic_sigma":     r"$\Vert \sigma(u)-\sigma_h\Vert_{L^2(\Omega)}$",
    "elastic_v":         r"$\Vert v-v_{\mathcal{T}}\Vert_{L^2(\Omega)}$",
}
_XLABEL = {"h": r"$h$", "h_min": r"$h$", "dt": r"$\Delta t$"}


def _slug(text):
    return re.sub(r"[^a-zA-Z0-9]+", "_", str(text)).strip("_").lower()


def _xcol(d):
    """Colonne des abscisses d'un cas (h_min sur maillages polygonaux, sinon h ou dt)"""
    data = d.get("data") or {}
    return d["x"] if d["x"] in data or not data else "h"


def ocvs(x, e):
    """Ordres observes segment par segment (None si non calculable)"""
    out = []
    for i in range(1, len(e)):
        ok = (e[i] is not None and e[i-1] is not None and e[i] > 0 and e[i-1] > 0
              and math.isfinite(e[i]) and math.isfinite(e[i-1]) and x[i] != x[i-1])
        out.append(np.log(e[i-1]/e[i]) / np.log(x[i-1]/x[i]) if ok else None)
    return out


def load_results(paths):
    datasets = []
    for path in paths:
        try:
            data = json.loads(Path(path).read_text())
        except Exception as e:
            print(f"  [WARN] {path}: {e}"); continue
        if "case" not in data or "columns" not in data:
            continue
        data["path"] = str(path)
        data.setdefault("variant", "")
        data.setdefault("style", {})
        data.setdefault("labels", {})
        datasets.append(data)
    return datasets


FLOOR = 0.8


def observed_order(x, e, expected, kind, window=None):
    """(ordre retenu, premier segment, plancher) : voir l'en-tete"""
    o = ocvs(x, e)
    n = len(o)
    if window is None:
        a, b = (0 if kind in ("h", "h_min") else 1), n
    else:
        a, b = slice(*window).indices(n)[:2]
    bad = lambda v: v is None or v < FLOOR * expected
    floor = b
    while floor > a and bad(o[floor - 1]):
        floor -= 1
    if floor == a:                                  # convergence lente partout : pas de plancher
        floor = b
    w = [v for v in o[max(a, floor - 2):floor] if v is not None]
    if not w:
        return None, a, floor
    return float(np.mean(w)), max(a, floor - 2), floor


def _bad_values(values):
    """Valeurs manquantes, non finies ou non positives d'une colonne"""
    return [i for i, v in enumerate(values) if v is None or not math.isfinite(v) or v <= 0.0]


def evaluate(d):
    """[(colonne, attendu, ocvs, retenu, statut)] pour un cas"""
    res = []
    quick_tol = 0.2 if d.get("quick") else 0.0
    if d.get("energy_check") is not None:
        en = d.get("energy") or {}
        var = en.get("max_variation")
        finite = var is not None and math.isfinite(var) and all(math.isfinite(v) for v in en.get("E", []))
        ok = finite and var <= d["energy_check"]
        res.append(("energy", d["energy_check"], [], var, "OK" if ok else "ECHEC"))
    elif d.get("energy") is not None and d.get("status") == "ok":
        E = d["energy"].get("E") or []
        if not E or not all(math.isfinite(v) for v in E):
            res.append(("energy", None, [], None, "ECHEC (valeurs non finies)"))
    data = d.get("data") or {}
    xk = _xcol(d)
    for c in d["columns"]:
        col = c["name"]
        if col not in data or xk not in data:
            res.append((col, c["expected"], [], None, "ECHEC" if d.get("status") != "skipped" else "IGNORE"))
            continue
        bad = _bad_values(data[col])
        missing = d.get("points") is not None and len(data[col]) != d["points"]
        o = ocvs(data[xk], data[col])
        obs, first, last = observed_order(data[xk], data[col], c["expected"], d["x"], c.get("window"))
        tol = c.get("tol", 0.3) + quick_tol
        if c.get("mode", "min") == "min":
            ok = obs is not None and obs >= c["expected"] - tol
        else:
            ok = obs is not None and abs(obs - c["expected"]) <= tol
        status = "OK" if ok else "ECHEC"
        if bad:
            status = f"ECHEC (valeur(s) non finie(s) : point(s) {', '.join(str(i + 1) for i in bad)})"
        elif missing:
            status = f"ECHEC ({len(data[col])} points au lieu de {d['points']})"
        res.append((col, c["expected"], o, obs, status))
    return res


def panels(datasets, variants=True):
    """{(article, panneau, variante): [cas...]} dans l'ordre d'apparition"""
    out = {}
    for d in datasets:
        if not variants and d.get("variant"):
            continue
        out.setdefault((d["paper"], d["panel"], d.get("variant", "")), []).append(d)
    return out


def _wrap(text, width):
    """Coupe chaque ligne de l'encadre apres une virgule (jamais dans une formule) au-dela de width caracteres"""
    out = []
    for line in text.split("\n"):
        cur = ""
        for part in line.split(", "):
            if cur and len(cur) + 2 + len(part) > width:
                out.append(cur + ","); cur = part
            else:
                cur = part if not cur else cur + ", " + part
        out.append(cur)
    return "\n".join(out)


def _box(ax, text, fontsize=8.5):
    ax.text(0.035, 0.955, text, transform=ax.transAxes, fontsize=fontsize, va="top", ha="left", color="#2b2b33",
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.6", fc="#eef2f8", ec="#9aa5b1", lw=0.9, alpha=0.95))


def _hidden_by(ax, bb, renderer):
    """Une courbe (segments echantillonnes a l'ecran) ou une annotation de l'axe passe-t-elle dans bb (pixels) ?"""
    for ln in ax.get_lines():
        xy = np.column_stack([np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)])
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if not len(xy):
            continue
        p = ax.transData.transform(xy)
        if len(p) > 1:
            s = np.linspace(0.0, 1.0, 16)[:, None, None]
            p = (p[:-1] + s * (p[1:] - p[:-1])).reshape(-1, 2)
        if np.any((p[:, 0] >= bb.x0) & (p[:, 0] <= bb.x1) & (p[:, 1] >= bb.y0) & (p[:, 1] <= bb.y1)):
            return True
    return any(a.get_window_extent(renderer).overlaps(bb) for a in ax.texts
               if isinstance(a, matplotlib.text.Annotation))


def _unhide_boxes(fig, steps=15):
    """Apres la mise en page : etend vers le haut l'axe des ordonnees tant que l'encadre (_box) cache des donnees"""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        for box in [t for t in ax.texts if t.get_bbox_patch() is not None and t.get_gid() != "ratio"]:
            bb = box.get_bbox_patch().get_window_extent(renderer)
            for _ in range(steps):
                if not _hidden_by(ax, bb, renderer):
                    break
                lo, hi = ax.get_ylim()
                ax.set_ylim(lo, hi * (hi / lo) ** 0.1 if ax.get_yscale() == "log" else hi + 0.1 * (hi - lo))


def _draw_panel(ax, group, column, norm, show_rate=True):
    """Courbes d'ordres observes d'un panneau pour une colonne d'erreur : une par cas (k ou s)"""
    cmap = _CMAP_NORM.get(norm, plt.cm.Greys)
    xs_all, es_all, drawn = [], [], []
    for d in sorted(group, key=lambda d: d["series_value"]):
        data = d.get("data") or {}
        if column not in data:
            continue
        x, e = data[_xcol(d)], data[column]
        pts = [(a, b) for a, b in zip(x, e) if b is not None and b > 0 and math.isfinite(b)]
        if not pts:
            continue
        x, e = [p[0] for p in pts], [p[1] for p in pts]
        xs_all += x; es_all += e
        v = d["series_value"]
        col = cmap(0.38 + 0.16 * min(v, 3))
        mrk = _MARKERS.get(v, "o")
        expected = next(c["expected"] for c in d["columns"] if c["name"] == column)
        sym = d["series"]
        label = f"${sym}={v}$ (ordre {expected:g})" if show_rate else f"${sym}={v}$"
        if show_rate and len(x) > 1:
            x_line = np.array([min(x) * 0.6, x[0] * 1.02])
            ax.loglog(x_line, e[0] * (x_line / x[0])**expected, "--", color=col, alpha=0.45, lw=1.3, zorder=1)
        ax.loglog(x, e, "-" + mrk, color=col, lw=1.8, ms=5.5, label=label, zorder=3)
        drawn.append((x, e, col))
    if not xs_all:
        ax.text(0.5, 0.5, "(pas de données)", ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="grey", style="italic")
        return False
    ax.set_xlim(min(xs_all) * 0.45, max(xs_all) * 2.2)
    ax.set_ylim(min(es_all) * 0.15, max(es_all) * 6.)
    for x, e, col in drawn:
        o = ocvs(x, e)
        for i in range(len(o)):
            if o[i] is None:
                continue
            xm = 10**(0.5*(np.log10(x[i]) + np.log10(x[i+1])))
            em = 10**(0.5*(np.log10(e[i]) + np.log10(e[i+1])))
            ax.annotate(f"{o[i]:.2f}", xy=(xm, em), xytext=(-14, 5), textcoords="offset points",
                        fontsize=6.5, color=col, fontweight="bold", clip_on=False)
    return True


def _norm_of(group, column):
    for d in group:
        for c in d["columns"]:
            if c["name"] == column:
                return c["norm"]
    return "L2"


def _columns_of(group):
    cols = []
    for d in group:
        for c in d["columns"]:
            if c["name"] not in cols:
                cols.append(c["name"])
    return cols


def _save(fig, base, fmts, saved, pad=0.1, png_dpi=150):
    """Enregistre la figure dans chaque format : PDF (ou SVG) a fond transparent (style P_Laplacien, pour LaTeX),
    PNG sur fond blanc (lisible dans tout visualiseur, theme sombre compris)"""
    Path(base).parent.mkdir(parents=True, exist_ok=True)
    for fmt in fmts:
        fname = f"{base}.{fmt}"
        if fmt == "png":
            fig.savefig(fname, dpi=png_dpi, bbox_inches="tight", facecolor="white", transparent=False, pad_inches=pad)
        else:
            fig.savefig(fname, dpi=150, bbox_inches="tight", transparent=True, pad_inches=pad)
        saved.append(fname)
        print(f"  [OK] {fname}")


def _figure_dir(out_dir, panel):
    """Dossier de la figure de l'article d'un panneau, "Fig. 3 (gauche)" -> <out_dir>/figure_3"""
    n = re.search(r"Fig\. (\d+)", panel).group(1)
    d = Path(out_dir) / f"figure_{n}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def plot_orders(datasets, out_dir, fmts=("png",), show_rate=True):
    """Ordres observes : une figure par (panneau, variante, norme), dans figure_<n>/ordres_observes"""
    saved = []
    for (art, pan, var), group in panels(datasets).items():
        for column in _columns_of(group):
            norm = _norm_of(group, column)
            fig, ax = plt.subplots(figsize=(5.5, 4.4))
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.grid(True, which="both", ls=":", alpha=0.35, color="grey", lw=0.8)
            ax.set_xlabel(_XLABEL[group[0]["x"]], fontsize=10)
            ax.set_ylabel(_YLABEL.get(column, column), fontsize=10)
            if not _draw_panel(ax, group, column, norm, show_rate):
                plt.close(fig); continue
            _box(ax, f"[{art}] {pan}" + (f" — variante : {var}" if var else "") + f"\n{_wrap(group[0]['box'], 70)}")
            ax.legend(fontsize=7, loc="lower right", framealpha=0.9)
            fig.tight_layout()
            _unhide_boxes(fig)
            name = f"{_slug(art)}_{_slug(pan)}" + (f"_{_slug(var)}" if var else "") + f"_{norm}_{_slug(column)}"
            _save(fig, _figure_dir(out_dir, pan) / "ordres_observes" / name, fmts, saved)
            plt.close(fig)
    return saved


# ------------------------------------------------------------------------------------------------------------------
# Figures reproduites, dessinees comme dans les articles

_POS_ROWS = {"": 0, "haut": 0, "milieu": 1, "bas": 2}


def position(panel):
    """("haut", "gauche") d'un panneau "Fig. 8 (haut, gauche)" ; ("", "gauche") pour "Fig. 3 (gauche)\""""
    m = re.search(r"\(([^)]*)\)", panel)
    parts = [p.strip() for p in m.group(1).split(",")] if m else ["gauche"]
    return (parts[0], parts[1]) if len(parts) == 2 else ("", parts[0])


def _panel_slug(panel):
    """"Fig. 8 (haut, gauche)" -> "fig8_haut_gauche\""""
    n = re.search(r"Fig\. (\d+)", panel).group(1)
    row, col = position(panel)
    return f"fig{n}_" + (f"{row}_" if row else "") + col


_PAPER_VALUES = {}


def paper_values(d):
    """Valeurs relevees sur les figures publiees (valeurs_article.json du dossier de l'article)"""
    folder = d.get("paper_dir") or ("acoustic_camc2022" if d["paper"] == "CAMC22" else "elastic_jsc2021")
    if folder not in _PAPER_VALUES:
        path = HERE / folder / "valeurs_article.json"
        _PAPER_VALUES[folder] = json.loads(path.read_text())["panels"] if path.exists() else {}
    return _PAPER_VALUES[folder]


def energy_loss(d):
    """(t, |E^n - E^0| / E^0) pour t > 0, valeurs nulles exclues (echelle log)"""
    en = d.get("energy") or {}
    t, E = np.asarray(en.get("t") or [], float), np.asarray(en.get("E") or [], float)
    if len(E) < 2 or not np.isfinite(E[0]) or E[0] == 0.0:
        return np.array([]), np.array([])
    loss = np.abs(E - E[0]) / abs(E[0])
    keep = (t > 0) & np.isfinite(loss) & (loss > 0)
    return t[keep], loss[keep]


def comparison(d):
    """{serie de l'article: [(x, nous, article)]} pour un cas a comparer (ref), aux abscisses communes"""
    ref = d.get("ref")
    if not ref or d.get("status") != "ok":
        return {}
    pv = paper_values(d).get(ref["panel"])
    if not pv:
        return {}
    out = {}
    if "t" in pv:                                   # perte d'energie : instant releve le plus proche
        t, loss = energy_loss(d)
        if not len(t):
            return {}
        for col, name in ref["series"].items():
            vals = pv["series"].get(name) or []
            pts = []
            for tp, yp in zip(pv["t"], vals):
                if yp is None or tp < 0.1:          # debut de courbe tres raide : releve peu precis
                    continue
                i = int(np.argmin(np.abs(t - tp)))
                pts.append((tp, float(loss[i]), yp))
            out[name] = pts
        return out
    data = d.get("data") or {}
    xk = _xcol(d)
    for col, name in ref["series"].items():
        if col not in data or xk not in data:
            continue
        vals = pv["series"].get(name) or []
        pts = []
        for xp, yp in zip(pv["x"], vals):
            if yp is None:
                continue
            for xo, yo in zip(data[xk], data[col]):
                if xo and abs(xo - xp) / xp < 0.06 and yo is not None and math.isfinite(yo) and yo > 0:
                    pts.append((xp, yo, yp))
        out[name] = pts
    return out


def _ratio_summary(points):
    """(mediane, min, max, nombre) des rapports nous/article ; les valeurs de l'article au niveau de l'arrondi
    (< 1e-12 : energie conservee de Newmark) sont du bruit et ne sont pas comparees"""
    r = [o / p for _, o, p in points if p > 1e-12]
    if not r:
        return None
    return float(np.median(r)), min(r), max(r), len(r)


def _verdict(summary):
    if summary is None:
        return "-"
    med, lo, hi, n = summary
    if abs(med - 1.0) <= 0.1 and lo >= 0.75 and hi <= 1.33:
        return "concorde"
    if abs(med - 1.0) <= 0.3:
        return "proche"
    return "écart"


def _decimals(v):
    for k in range(7):
        if abs(round(v, k) - v) <= 1e-9 * max(1.0, abs(v)):
            return k
    return 6


def _log_axis(ax, axis, lo, hi):
    """Graduations comme les figures des articles : 1-5 (0.05, 0.10, 0.50, 1) en abscisse, decades en ordonnee
    (0.100, 0.010, 0.001, puis 10^-4 ...)"""
    if axis == "x":
        ticks = []
        for mant in ((1, 5), (1, 2, 5), (1, 1.5, 2, 2.5, 3, 5, 7)):
            ticks = [m * 10.0 ** e for e in range(-7, 2) for m in mant if lo / 1.04 <= m * 10.0 ** e <= hi * 1.04]
            if len(ticks) >= 3:
                break
        dec = max([_decimals(t) for t in ticks if t < 1] or [0])
        fmt = lambda v, pos: ("%d" % round(v)) if v >= 1 and abs(v - round(v)) < 1e-9 else f"{v:.{dec}f}"
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(fmt))
        ax.xaxis.set_minor_locator(FixedLocator([m * 10.0 ** e for e in range(-7, 2) for m in range(1, 10)]))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        e0, e1 = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
        step = max(1, int(math.ceil((e1 - e0 + 1) / 7.0)))
        es = [e for e in range(e0, e1 + 1) if (e - e1) % step == 0]
        ticks = [10.0 ** e for e in es if lo <= 10.0 ** e <= hi]
        def fmt(v, pos):
            e = int(round(math.log10(v)))
            if e >= 1:
                return "%d" % round(v)
            if e == 0:
                return "1"
            if e >= -3:
                return f"{v:.3f}"
            return rf"$10^{{{e}}}$"
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(fmt))
        ax.yaxis.set_minor_locator(FixedLocator([m * 10.0 ** e for e in range(e0 - 1, e1 + 1) for m in range(2, 10)]
                                                if step == 1 else [10.0 ** e for e in range(e0 - 1, e1 + 2)]))
        ax.yaxis.set_minor_formatter(NullFormatter())


def _paper_axes(ax, fs=10):
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=fs - 1)
    ax.tick_params(which="major", length=4.5)
    ax.tick_params(which="minor", length=2.2)
    for s in ax.spines.values():
        s.set_linewidth(1.0)
    ax.grid(True, which="major", color="#d9d9d9", lw=0.7)
    ax.set_axisbelow(True)


def _legend(ax, title, loc, fs=8):
    hs, ls = ax.get_legend_handles_labels()
    if not hs:
        return
    leg = ax.legend(hs, ls, title=title, loc=loc, fontsize=fs, title_fontsize=fs, frameon=True, fancybox=True,
                    edgecolor="#6e6e6e", framealpha=1.0, handlelength=2.6, borderpad=0.5, labelspacing=0.35)
    leg.get_frame().set_linewidth(0.8)
    for txt, h in zip(leg.get_texts(), getattr(leg, "legend_handles", None) or leg.legendHandles):
        c = h.get_color() if hasattr(h, "get_color") else None
        if c is not None:
            txt.set_color(c if not isinstance(c, str) or c != "none" else "black")


def draw_paper_conv(ax, group, overlay=False, fs=10):
    """Panneau d'erreurs comme dans l'article : L2 en tirets fins fonces, H1 en tirets epais clairs (paire foncee
    bleu/rouge pour les etudes en temps a une seule courbe par norme), legende de l'article"""
    xs_all, ys_all, ratios = [], [], []
    entries = []
    for d in group:
        data = d.get("data") or {}
        xk = _xcol(d)
        if xk not in data:
            continue
        st = d.get("style") or {}
        c = st.get("color", 0)
        cols = [cc["name"] for cc in d["columns"]]
        for j, col in enumerate(cols[:2]):
            if col not in data:
                continue
            pts = [(a, b) for a, b in zip(data[xk], data[col]) if b is not None and math.isfinite(b) and b > 0]
            if not pts:
                continue
            x, y = zip(*pts)
            xs_all += x; ys_all += y
            label = (d.get("labels") or {}).get(col, f"{'L$_2$' if j == 0 else 'H$^1$'} : {d['series']}={d['series_value']}")
            if st.get("pair"):          # etude en temps : L2 bleu (disques), H1 rouge fonce (carres), traits fins
                kw = dict(color=_rgb(_DARK[j]), lw=1.4, dashes=(5, 2.5), marker="os"[j] if st.get("markers") else None,
                          ms=4.5)
                order = j
            elif j == 0:
                kw = dict(color=_rgb(_DARK[c % 3]), lw=1.4, dashes=(5, 2.5),
                          marker=_L2_MARKERS[c % 3] if st.get("markers") else None, ms=4.5)
                order = c
            else:
                kw = dict(color=_rgb(_LIGHT[c % 3]), lw=3.6, dashes=(3.0, 1.4),
                          marker=_H1_MARKERS[c % 3] if st.get("markers") else None, ms=4.5)
                order = 10 + c
            entries.append((order, x, y, label, kw))
    if not xs_all:
        ax.text(0.5, 0.5, "(pas de données)", ha="center", va="center", transform=ax.transAxes, fontsize=9,
                color="grey", style="italic")
        return None
    for order, x, y, label, kw in sorted(entries, key=lambda e: e[0]):
        z = 3 if order < 10 else 2
        ax.loglog(x, y, label=label, zorder=z, **kw)
    if overlay:
        first = True
        for d in group:
            for name, pts in comparison(d).items():
                if pts:
                    ax.loglog([p[0] for p in pts], [p[2] for p in pts], "o", color="black", mfc="none", mec="black",
                              mew=1.0, ms=6.5, zorder=5, label="article (valeurs relevées)" if first else None)
                    first = False
                    ratios += pts
                    ys_all += [p[2] for p in pts]
    lo, hi = min(xs_all), max(xs_all)
    ax.set_xlim(lo / 1.25, hi * 1.25)
    ylo, yhi = min(ys_all), max(ys_all)
    ax.set_ylim(ylo / 2.5, yhi * 2.5)
    _log_axis(ax, "x", lo / 1.25, hi * 1.25)
    _log_axis(ax, "y", ylo / 2.5, yhi * 2.5)
    xk = group[0]["x"]
    ax.set_xlabel({"h": "h", "h_min": "h", "dt": r"$\Delta$t"}.get(xk, xk), fontsize=fs)
    ax.set_ylabel("Error", fontsize=fs)
    _paper_axes(ax, fs)
    _legend(ax, group[0].get("legend") or None, (group[0].get("style") or {}).get("loc", "lower right"), fs - 2)
    return _ratio_summary(ratios) if overlay else None


def draw_paper_energy(ax, group, overlay=False, fs=10):
    """Panneau de perte d'energie relative comme dans l'article (echelle log, temps lineaire)"""
    ys_all, ratios = [], []
    first = True
    for d in group:
        t, loss = energy_loss(d)
        if not len(t):
            continue
        st = d.get("style") or {}
        c = st.get("color", 0) % 3
        if st.get("solid"):
            kw = dict(color=_rgb(_DARK[c]), lw=1.6)
        elif st.get("thick"):
            kw = dict(color=_rgb(_LIGHT[c]), lw=3.6, dashes=(3.0, 1.4))
        else:
            kw = dict(color=_rgb(_DARK[c]), lw=1.4, dashes=(5, 2.5))
        label = (d.get("labels") or {}).get("energy", d["case"])
        ax.semilogy(t, loss, label=label, zorder=3 if not st.get("thick") else 2, **kw)
        ys_all += list(loss)
    if overlay:                                     # valeurs de l'article apres les courbes (derniere entree)
        for d in group:
            for name, pts in comparison(d).items():
                if pts:
                    ax.semilogy([p[0] for p in pts], [p[2] for p in pts], "o", color="black", mfc="none", mec="black",
                                mew=1.0, ms=6, zorder=5, label="article (valeurs relevées)" if first else None)
                    first = False
                    ratios += pts
                    ys_all += [p[2] for p in pts]
    if not ys_all:
        ax.text(0.5, 0.5, "(pas de données)", ha="center", va="center", transform=ax.transAxes, fontsize=9,
                color="grey", style="italic")
        return None
    ylo, yhi = max(min(ys_all), 1e-17), max(ys_all)
    ax.set_ylim(ylo / 3.0, yhi * 6.0)
    _log_axis(ax, "y", ylo / 3.0, yhi * 6.0)
    ax.set_xlim(-0.02, 1.02)
    ax.xaxis.set_major_locator(FixedLocator([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:.1f}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("time", fontsize=fs)
    ax.set_ylabel("Energy loss", fontsize=fs)
    _paper_axes(ax, fs)
    ax.grid(True, which="major", axis="x", color="#bdbdbd", lw=0.9)
    _legend(ax, group[0].get("legend") or None, (group[0].get("style") or {}).get("loc", "lower right"), fs - 2)
    return _ratio_summary(ratios) if overlay else None


def draw_paper_panel(ax, group, overlay=False, fs=10):
    if any(d.get("energy") for d in group) and not any(d["columns"] for d in group):
        return draw_paper_energy(ax, group, overlay, fs)
    return draw_paper_conv(ax, group, overlay, fs)


def _ratio_text(ax, summary):
    if summary is None:
        return
    med, lo, hi, n = summary
    ax.text(0.03, 0.97, f"nous/article : {med:.2f} ({lo:.2f}–{hi:.2f}, {n} points)", transform=ax.transAxes,
            fontsize=7.5, va="top", ha="left", gid="ratio",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff8e6", ec="#c9a227", lw=0.8))


def _layout(panel_list):
    """Disposition des panneaux d'une figure comme dans l'article : lignes haut/milieu/bas, colonnes gauche/droite"""
    rows = sorted({_POS_ROWS[position(p)[0]] for p in panel_list})
    return {p: (rows.index(_POS_ROWS[position(p)[0]]), 0 if position(p)[1] == "gauche" else 1) for p in panel_list}, len(rows)


def plot_paper_figures(datasets, out_dir, fmts=("png",)):
    """Figures reproduites : chaque panneau seul, la figure composee, la comparaison, les variantes ; renvoie
    (fichiers, {(article, panneau, variante): resume nous/article})"""
    saved, summaries = [], {}
    groups = panels(datasets)
    figures = collections.OrderedDict()
    for (art, pan, var), group in groups.items():
        figures.setdefault((art, int(re.search(r"Fig\. (\d+)", pan).group(1))), []).append((pan, var, group))
    for (art, n), items in figures.items():
        fdir = Path(out_dir) / f"figure_{n}"
        fdir.mkdir(parents=True, exist_ok=True)
        # panneaux seuls (figure publiee et variantes) et resumes nous/article
        for pan, var, group in items:
            fig, ax = plt.subplots(figsize=(4.8, 4.6))
            draw_paper_panel(ax, group)
            ax.set_title(f"[{art}] {pan}" + ("" if not var else f"\nvariante : {var}"), fontsize=9)
            fig.tight_layout()
            base = fdir / (("variantes/" if var else "") + f"{art}_{_panel_slug(pan)}" + (f"_{_slug(var)}" if var else ""))
            _save(fig, base, fmts, saved)
            plt.close(fig)
            if not var:
                fig, ax = plt.subplots(figsize=(4.8, 4.6))
                summaries[(art, pan, var)] = draw_paper_panel(ax, group, overlay=True)
                plt.close(fig)
        # figure composee et comparaison : panneaux de la figure publiee
        main = [(pan, group) for pan, var, group in items if not var]
        if not main:
            continue
        lay, nrows = _layout([p for p, _ in main])
        ncols = 2 if any(c == 1 for _, c in lay.values()) or len(main) > 1 else 1
        for overlay in (False, True):
            fig, axes = plt.subplots(nrows, ncols, figsize=(4.9 * ncols, 4.7 * nrows), squeeze=False)
            used = set()
            for pan, group in main:
                r, c = lay[pan]
                c = min(c, ncols - 1)
                used.add((r, c))
                ax = axes[r][c]
                summ = draw_paper_panel(ax, group, overlay=overlay)
                ax.set_title(f"{pan}\n{_wrap(group[0]['box'].split(chr(10))[0], 62)}", fontsize=8.5)
                if overlay:
                    _ratio_text(ax, summ)
            for r in range(nrows):
                for c in range(ncols):
                    if (r, c) not in used:
                        axes[r][c].axis("off")
            title = f"[{art}] Figure {n} — reproduction avec wave_propagation"
            if overlay:
                title += "\n(cercles : valeurs relevées sur la figure publiée)"
            fig.suptitle(title, fontsize=11, fontweight="bold")
            fig.tight_layout(rect=[0, 0, 1, 0.97 if not overlay else 0.955])
            _save(fig, fdir / (f"{art}_figure_{n}" + ("_comparaison" if overlay else "")), fmts, saved)
            plt.close(fig)
    return saved, summaries


def _read_polymesher(path):
    """Polygones d'un maillage PolyMesher (en-tete : sommets, cellules ; indices a partir de 1)"""
    import bz2
    opener = bz2.open if str(path).endswith(".bz2") else open
    with opener(path, "rt") as f:
        lines = [l.split() for l in f.read().splitlines() if l.strip()]
    nv, nc = int(lines[0][0]), int(lines[0][1])
    pts = np.array([[float(a), float(b)] for a, b in (l[:2] for l in lines[1:1 + nv])])
    cells = [[int(i) - 1 for i in l[1:1 + int(l[0])]] for l in lines[1 + nv:1 + nv + nc]]
    return pts, cells


def plot_fig7(out_dir, fmts=("png",)):
    """[CAMC22] Fig. 7 : les quatre premiers maillages PolyMesher (20 x 2^l cellules, l = 0..3)"""
    from matplotlib.collections import PolyCollection
    saved = []
    names = [f"unit_square_polymesh_nel_{20 * 2 ** l}.txt" for l in range(4)]
    dirs = [Path(out_dir) / "figure_8" / "meshes", HERE / "acoustic_camc2022" / "meshes"]
    files = []
    for nm in names:
        f = next((d / s for d in dirs for s in (nm, nm + ".bz2") if (d / s).exists()), None)
        if f is None:
            return saved
        files.append(f)
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.5))
    for ax, f, l in zip(axes, files, range(4)):
        pts, cells = _read_polymesher(f)
        pc = PolyCollection([pts[c] for c in cells], facecolors="#dbe9f9", edgecolors="#e3261f", linewidths=0.8)
        ax.add_collection(pc)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"{len(cells)} cellules", fontsize=9)
    fig.suptitle("[CAMC22] Figure 7 — maillages polygonaux (PolyMesher) : 20 x 2^l cellules, l = 0..3", fontsize=11,
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, Path(out_dir) / "figure_7" / "CAMC22_figure_7", fmts, saved)
    plt.close(fig)
    return saved


def plot_summary_grid(datasets, out_dir, fmts=("png",)):
    """Recapitulatif par article : tous les panneaux des figures reproduites (figure publiee), deux par ligne"""
    saved = []
    by_paper = collections.OrderedDict()
    for (art, pan, var), group in panels(datasets, variants=False).items():
        by_paper.setdefault(art, []).append((pan, group))
    for art, items in by_paper.items():
        nrows = (len(items) + 1) // 2
        fig, axes = plt.subplots(nrows, 2, figsize=(9.8, 4.6 * nrows), squeeze=False)
        for i, (pan, group) in enumerate(items):
            ax = axes[i // 2][i % 2]
            draw_paper_panel(ax, group, overlay=True)
            ax.set_title(f"{pan}\n{_wrap(group[0]['box'].split(chr(10))[0], 62)}", fontsize=8.5)
        if len(items) % 2:
            axes[-1][-1].axis("off")
        fig.suptitle(f"[{art}] figures reproduites avec wave_propagation (cercles : valeurs relevées sur les figures "
                     f"publiées)", fontsize=12, fontweight="bold", y=0.999)
        fig.tight_layout(rect=[0, 0, 1, 0.995])
        _save(fig, Path(out_dir) / f"{art}_recapitulatif", fmts, saved, png_dpi=90)
        plt.close(fig)
    return saved


def print_summary(datasets):
    n_ok, n_fail = 0, 0
    print("\n" + "=" * 118)
    print("  article  panneau                    cas                                       grandeur            attendu  retenu  statut")
    print("  " + "-" * 114)
    for d in datasets:
        modes = {c["name"]: c.get("mode", "min") for c in d["columns"]}
        for col, exp, o, obs, st in evaluate(d):
            n_ok += st == "OK"; n_fail += st.startswith("ECHEC")
            if col == "energy":
                exp_s = "-" if exp is None else f"<= {exp:.0e}"
                obs_s = "   -   " if obs is None else f"{obs:.1e}"
            else:
                exp_s = (">= " if modes.get(col) == "min" else "= ") + f"{exp:g}"
                obs_s = "   -   " if obs is None else f"{obs:7.2f}"
            print(f"  {d['paper']:7s}  {d['panel'][:26]:26s} {d['case'][:41]:41s} {col[:18]:18s} {exp_s:>8s} {obs_s:>7s}  {st}")
    print("=" * 118)
    print(f"  {n_ok} verification(s) OK, {n_fail} echec(s)\n")
    return n_ok, n_fail


def _img(files, out_dir, pred, width="49%"):
    png = [f for f in files if pred(f) and f.endswith(".png")]
    return " ".join(f'<img src="{os.path.relpath(f, out_dir)}" width="{width}">' for f in png)


def write_report(datasets, out_dir, figures, summaries):
    """rapport.md : par figure, les images, puis par panneau les parametres, les ecarts a la legende, les rapports
    nous/article et les ordres observes"""
    out_dir = Path(out_dir)
    papers = collections.OrderedDict()
    for d in datasets:
        papers.setdefault(d["paper"], d.get("paper_ref", d["paper"]))
    lines = ["# Reproduction des figures de " + " et ".join(f"[{p}]" for p in papers) + " avec wave_propagation\n",
             f"- date : {time.strftime('%Y-%m-%d %H:%M')}",
             "- exécutable(s) : " + ", ".join(f"`{e}` ({n} cas)" for e, n in
                                         collections.Counter(d.get("exe", "?") for d in datasets).items()),
             "- figures : `figure_<n>/<ART>_figure_<n>.png` (figure reproduite, dessinée comme dans l'article), "
             "`..._comparaison.png` (avec les valeurs relevées sur la figure publiée, cercles noirs, et le rapport "
             "nous/article de chaque panneau : médiane, min–max), `variantes/`, `ordres_observes/`",
             "- comparaison : rapport nous/article aux points relevés sur la figure publiée (`valeurs_article.json`, "
             "précision du relevé ~ 2 %) ; « concorde » si la médiane est à 10 % de 1 et tous les rapports dans "
             "[0.75, 1.33], « proche » si la médiane est à 30 % de 1",
             "- ordres : moyenne des deux derniers ordres observés (OCV) fiables, plancher final exclu ; statut OK si "
             "retenu >= attendu - tolérance (« >= ») ou |retenu - attendu| <= tolérance (« = », tau = O(1)) ; tolérance "
             "0.3 (0.5 en mode rapide) ; ECHEC aussi pour une valeur manquante ou non finie\n"]
    n_ok = n_fail = 0
    body = []
    groups = panels(datasets)
    for art, refname in papers.items():
        body.append(f"## [{art}] {refname}\n")
        figs = sorted({d["figure"] for d in datasets if d["paper"] == art})
        if art == "CAMC22" and any(f.endswith("CAMC22_figure_7.png") for f in figures):
            figs = sorted(set(figs) | {7})
        for n in figs:
            body.append(f"### Figure {n}\n")
            imgs = _img(figures, out_dir, lambda f: re.search(rf"figure_{n}/{art}_figure_{n}(_comparaison)?\.png$", f))
            if imgs:
                body += [imgs, ""]
            if n == 7 and art == "CAMC22":
                body.append("Les quatre premiers maillages PolyMesher de l'article (20, 40, 80, 160 cellules), lus "
                            "dans les fichiers de `meshes/` utilisés pour la Fig. 8.\n")
            for (a, pan, var), group in groups.items():
                if a != art or group[0]["figure"] != n:
                    continue
                title = f"#### {pan}" + (f" — variante : {var}" if var else "") + f" — {group[0]['box'].replace(chr(10), ', ')}"
                body.append(title + "\n")
                for note in sorted(set(d.get("note", "") for d in group if d.get("note"))):
                    body.append(f"*{note}*\n")
                if var:
                    imgs = _img(figures, out_dir, lambda f: "/variantes/" in f and Path(f).stem.startswith(
                        f"{art}_{_panel_slug(pan)}_{_slug(var)}"), width="40%")
                    if imgs:
                        body += [imgs, ""]
                summ = summaries.get((a, pan, var))
                if summ:
                    med, lo, hi, npts = summ
                    body.append(f"**Comparaison à la figure publiée : nous/article = {med:.2f} (médiane ; {lo:.2f}–"
                                f"{hi:.2f} sur {npts} points) — {_verdict(summ)}**\n")
                body.append("| cas | grandeur | nous/article (médiane, min–max) | attendu | ordres observés (OCV) | retenu | statut |")
                body.append("|---|---|---|---|---|---|---|")
                for d in group:
                    comp = comparison(d)
                    modes = {c["name"]: c.get("mode", "min") for c in d["columns"]}
                    res = evaluate(d)
                    if not res and d.get("energy"):
                        var_e = d["energy"].get("max_variation")
                        s = _ratio_summary(next(iter(comp.values()), [])) if comp else None
                        rs = "-" if s is None else f"{s[0]:.2f} ({s[1]:.2f}–{s[2]:.2f})"
                        body.append(f"| {d['case']} | énergie | {rs} | - | - | {'-' if var_e is None else f'{var_e:.2e}'} "
                                    f"| info" + ("" if d.get("status") == "ok" else f" ({d.get('status')})") + " |")
                    for col, exp, o, obs, st in res:
                        n_ok += st == "OK"; n_fail += st.startswith("ECHEC")
                        if col == "energy":
                            s = _ratio_summary(next(iter(comp.values()), [])) if comp else None
                            rs = "-" if s is None else f"{s[0]:.2f} ({s[1]:.2f}–{s[2]:.2f})"
                            body.append(f"| {d['case']} | énergie | {rs} | " + ("-" if exp is None else f"<= {exp:.0e}")
                                        + f" | variation relative max | {'-' if obs is None else f'{obs:.2e}'} | **{st}** |")
                            continue
                        name = (d.get("ref") or {}).get("series", {}).get(col)
                        s = _ratio_summary(comp.get(name, [])) if name else None
                        rs = "-" if s is None else f"{s[0]:.2f} ({s[1]:.2f}–{s[2]:.2f})"
                        oc = " ".join("-" if v is None else f"{v:.2f}" for v in o)
                        body.append(f"| {d['case']} | `{col}` | {rs} | {'>= ' if modes.get(col) == 'min' else '= '}{exp:g} | "
                                    f"{oc} | {'-' if obs is None else f'{obs:.2f}'} | **{st}**"
                                    + ("" if d.get("status") == "ok" else f" ({d.get('status')})") + " |")
                body.append("")
    verdicts = collections.Counter(_verdict(s) for s in summaries.values() if s)
    lines.append(f"**Bilan : {n_ok} vérification(s) d'ordre OK, {n_fail} échec(s) ; comparaison aux figures publiées : "
                 + ", ".join(f"{verdicts.get(v, 0)} panneau(x) {v}" for v in ("concorde", "proche", "écart")) + ".**\n")
    if figures:
        lines.append("Récapitulatif : " + ", ".join(f"[{Path(f).name}]({os.path.relpath(f, out_dir)})"
                                                    for f in figures if "_recapitulatif." in f) + "\n")
    path = out_dir / "rapport.md"
    path.write_text("\n".join(lines + body))
    print(f"  [OK] {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="dossier de resultats d'un article (figure_<n>/<cas>/<cas>.json), "
                                              "ou fichiers JSON")
    ap.add_argument("--output", "-o", default=None, help="dossier de sortie (defaut: le premier dossier donne)")
    ap.add_argument("--format", default="png", help="formats des figures, separes par des virgules : png, pdf, svg "
                                                  "(defaut: png)")
    ap.add_argument("--no-rate-line", action="store_true", help="sans droite de l'ordre attendu (ordres observes)")
    args = ap.parse_args()

    files = []
    for inp in args.inputs:
        p = Path(inp)
        if p.is_dir():
            files += sorted(p.glob("figure_*/*/*.json")) or sorted(p.glob("*/*.json"))
        elif p.exists():
            files.append(p)
    files = list(dict.fromkeys(str(f) for f in files))
    if not files:
        print("[ERREUR] Aucun fichier JSON."); sys.exit(1)
    print(f"Chargement de {len(files)} fichier(s) ...")
    ds = load_results(files)
    if not ds:
        print("[ERREUR] Aucun resultat valide."); sys.exit(1)
    order = {"CAMC22": 0, "JSC21": 1}
    ds.sort(key=lambda d: (order.get(d["paper"], 9), d.get("order", 0)))
    print(f"  -> {len(ds)} cas\n")

    out_dir = Path(args.inputs[0]) if Path(args.inputs[0]).is_dir() else Path(args.inputs[0]).parent
    out_dir = Path(args.output) if args.output else out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    n_ok, n_fail = print_summary(ds)
    fmts = [f.strip() for f in args.format.split(",") if f.strip()]
    if not fmts or any(f not in ("pdf", "png", "svg") for f in fmts):
        print(f"[ERREUR] --format {args.format} : pdf, png ou svg, separes par des virgules"); sys.exit(1)
    print("Generation des figures reproduites (dossier figure_<n> de chaque figure de l'article) ...")
    figs, summaries = plot_paper_figures(ds, out_dir, fmts)
    if any(d["paper"] == "CAMC22" and d["figure"] == 8 for d in ds):
        figs += plot_fig7(out_dir, fmts)
    print("Generation des figures d'ordres observes (figure_<n>/ordres_observes) ...")
    figs += plot_orders(ds, out_dir, fmts, show_rate=not args.no_rate_line)
    figs += plot_summary_grid(ds, out_dir, fmts)
    write_report(ds, out_dir, figs, summaries)
    print("\nComparaison aux figures publiees (nous/article, mediane et min-max par panneau) :")
    for (art, pan, var), s in summaries.items():
        if s:
            print(f"  [{art}] {pan:26s} {s[0]:6.2f}  ({s[1]:.2f}-{s[2]:.2f}, {s[3]} points)  {_verdict(s)}")
    print("\nTermine.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
