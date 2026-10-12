#!/usr/bin/env python3
"""
plot_mottier.py — figures, comparaisons avec les articles et rapport des campagnes papers_Mottier.

Chaque article fournit son module figures.py (dossier de l'article) : liste FIGURES de fonctions f(ctx) qui tracent
une figure ou un tableau de l'article a partir des resultats (ctx.case(nom), ctx.cases(regex)), des references
Gar6more2D (ctx.reference(nom, capteur)) et des valeurs relevees sur l'article (ctx.paper[cle]), enregistrent les PNG
(ctx.save) et ajoutent une section au rapport (ctx.section). Ce module fournit le contexte et les traces communs.

    python3 plot_mottier.py <dossier de l'article> [<dossier des resultats>]
"""
import importlib, json, math, re, sys, traceback
sys.dont_write_bytecode = True
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from campaign_mottier import PAPERS, item_dir     # noqa: E402

# couleurs des figures des articles (matplotlib nommees)
GREEN, ORANGE, DARKRED, NAVY = "#006400", "#ffa500", "#8b0000", "#000080"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.color": "#b0b0b0", "grid.linewidth": 0.6,
                     "legend.fontsize": 9, "figure.dpi": 100})


def _num_key(item):
    m = re.match(r"(Fig\.|Table)\s*(\d+)", item)
    return (0 if m.group(1) == "Fig." else 1, int(m.group(2)))


class Context:
    def __init__(self, paper_dir, out, fmts=("png",)):
        self.paper_dir, self.out, self.fmts = Path(paper_dir), Path(out), tuple(fmts)
        self.results = {}
        for f in sorted(self.out.glob("*/*/*.json")):
            try:
                r = json.loads(f.read_text())
            except Exception:
                continue
            if isinstance(r, dict) and "case" in r:
                r["_dir"] = str(f.parent)
                self.results[r["case"]] = r
        pv = self.paper_dir / "valeurs_article.json"
        self.paper = json.loads(pv.read_text()) if pv.exists() else {}
        self.sections = {}
        self.quick = any(r.get("quick") for r in self.results.values())

    # resultats
    def case(self, name):
        r = self.results.get(name)
        return r if r and r.get("status") == "ok" else None

    def cases(self, regex):
        return [r for n, r in sorted(self.results.items()) if re.search(regex, n) and r.get("status") == "ok"]

    def case_path(self, name, filename):
        r = self.results.get(name)
        return Path(r["_dir"]) / filename if r else None

    def reference(self, name, sensor):
        import gar6more
        p = self.out / "references" / name / f"{sensor}.csv"
        return gar6more.load(self.out / "references", name, sensor) if p.exists() else None

    # sorties
    def item_dir(self, item):
        d = self.out / item_dir(item)
        d.mkdir(parents=True, exist_ok=True)
        return d

    def save(self, fig, item, name, dpi=150):
        d = self.item_dir(item)
        paths = []
        for fmt in self.fmts:
            p = d / f"{name}.{fmt}"
            fig.savefig(p, dpi=dpi, bbox_inches="tight")
            paths.append(p)
        plt.close(fig)
        return str(paths[0].relative_to(self.out))

    def section(self, item, title, images=(), notes=(), table=None, status="", params=()):
        """Section du rapport : images (chemins relatifs), notes (texte, ecarts), table ((entetes, lignes)) de
        comparaison, statut ("reproduit", "ecart explique", "ecart"...), parametres"""
        s = self.sections.setdefault(item, {"item": item, "title": title, "images": [], "notes": [], "tables": [],
                                            "status": [], "params": []})
        s["images"] += list(images); s["notes"] += list(notes); s["params"] += list(params)
        if table:
            s["tables"].append(table)
        if status:
            s["status"].append(status)


# ------------------------------------------------------------------------------------------------------------------
# traces communs

def legend(ax, **kw):
    """Legende des courbes nommees du panneau (rien si le panneau n'en a pas : cas non calcules)"""
    if ax.get_legend_handles_labels()[0]:
        ax.legend(**kw)


def _minor_2_5(x, pos=None):
    """Etiquettes des graduations secondaires 2 et 5 x 10^e seulement (les autres se chevauchent)"""
    if x <= 0:
        return ""
    e = math.floor(math.log10(x) + 1e-9)
    m = x / 10.0 ** e
    return f"${round(m)}\\times10^{{{e}}}$" if abs(m - round(m)) < 1e-6 and round(m) in (2, 5) else ""


def loglog(ax, xlabel, ylabel):
    """Axes logarithmiques ; sur un axe de moins d'une decade et demie (une seule graduation principale ou presque),
    graduations secondaires etiquetees a 2 et 5 x 10^e seulement"""
    from matplotlib.ticker import FuncFormatter, NullFormatter
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
        if lim[0] > 0 and math.log10(lim[1] / lim[0]) < 1.5:
            axis.set_minor_formatter(FuncFormatter(_minor_2_5))
        else:
            axis.set_minor_formatter(NullFormatter())


def slope_triangles(ax, x0, x1, y0, orders, colors, label=True):
    """Triangles de pente (comme les figures des articles) : base [x0, x1] a la hauteur y0, un cote par ordre"""
    for p, c in zip(orders, colors):
        y1 = y0 * (x1 / x0) ** p
        ax.plot([x0, x1, x1], [y0, y1, y0], color=c, lw=1.5)
        if label:
            ax.text(x1 * 1.04, y1, str(p), color=c, fontsize=8, va="center")
    ax.plot([x0, x1], [y0, y0], color="k", lw=1.5)


def ratio_stats(ours, paper):
    """Rapports nous / article des points communs : (mediane, min, max, nombre)"""
    r = [o / p for o, p in zip(ours, paper) if o is not None and p is not None and p != 0 and np.isfinite(o)
         and o != 0]
    if not r:
        return None
    return float(np.median(r)), float(min(r)), float(max(r)), len(r)


def fmt_ratio(st):
    if st is None:
        return "-"
    med, lo, hi, n = st
    return f"{med:.3g} [{lo:.3g}, {hi:.3g}] ({n})"


def table_image(ctx, item, name, header, rows, title="", col_width=1.25, fontsize=9, highlight=None):
    """Tableau (liste de lignes de textes) enregistre en PNG ; highlight : {(ligne, colonne): couleur}"""
    nr, nc = len(rows) + 1, len(header)
    fig, ax = plt.subplots(figsize=(max(4.0, col_width * nc), 0.32 * nr + (0.5 if title else 0.15)))
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(fontsize); tb.scale(1.0, 1.25)
    for (i, j), c in (highlight or {}).items():
        tb[(i + 1, j)].set_facecolor(c)
    for j in range(nc):
        tb[(0, j)].set_facecolor("#e0e0e0")
    if title:
        ax.set_title(title, fontsize=fontsize + 1)
    return ctx.save(fig, item, name)


def markdown_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# ------------------------------------------------------------------------------------------------------------------
# maillages et instantanes

def read_polygons(path):
    import meshes
    return meshes.read_polygons(path)


def draw_mesh(ax, pts, cells, tags, colors, lw=0.6):
    for t in sorted(set(tags)):
        polys = [pts[c] for c, tt in zip(cells, tags) if tt == t]
        ax.add_collection(PolyCollection(polys, facecolors="none", edgecolors=colors.get(t, "k"), linewidths=lw))
    ax.set_xlim(pts[:, 0].min(), pts[:, 0].max()); ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
    ax.set_aspect("equal"); ax.axis("off")


def read_snapshot(cells_file, snap_file):
    """Cellules (milieu, sommets) et valeurs (barycentre puis sommets) d'un instantane"""
    cells = []
    for line in Path(cells_file).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        v = line.split()
        med, nv = int(v[0]), int(v[1])
        xy = np.array([float(a) for a in v[2:2 + 2 * nv]]).reshape(nv, 2)
        cells.append((med, xy))
    lines = Path(snap_file).read_text().splitlines()
    t = None
    m = re.search(r"t = (\S+)", lines[0])
    if m:
        t = float(m.group(1))
    vals = [np.array([float(a) for a in l.split()]) for l in lines if l.strip() and not l.startswith("#")]
    return cells, vals, t


def draw_snapshot(ax, cells, vals, cmap_a="seismic", cmap_e="viridis", clim_a=None, clim_e=None):
    """Instantane lisse : chaque cellule decoupee en triangles (barycentre, aretes), valeurs du polynome de la cellule
    au barycentre et aux sommets, interpolation lineaire (Gouraud) ; pression (cellules acoustiques) et norme de la
    vitesse (cellules elastiques) avec deux echelles de couleurs"""
    out = {}
    for med, cmap, clim in ((1, cmap_a, clim_a), (0, cmap_e, clim_e)):
        X, Y, Z, T = [], [], [], []
        for (m, xy), v in zip(cells, vals):
            if m != med:
                continue
            base = len(X)
            bar = xy.mean(axis=0)
            X.append(bar[0]); Y.append(bar[1]); Z.append(v[0])
            for j, p in enumerate(xy):
                X.append(p[0]); Y.append(p[1]); Z.append(v[1 + j] if len(v) > 1 + j else v[0])
            n = len(xy)
            for j in range(n):
                T.append([base, base + 1 + j, base + 1 + (j + 1) % n])
        if not T:
            continue
        tri = Triangulation(np.array(X), np.array(Y), np.array(T))
        Z = np.array(Z)
        if clim is None:
            a = np.abs(Z).max() or 1.0
            clim = (-a, a) if med == 1 else (0.0, a)
        out[med] = ax.tripcolor(tri, Z, shading="gouraud", cmap=cmap, vmin=clim[0], vmax=clim[1], rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    return out


# ------------------------------------------------------------------------------------------------------------------
# rapport

def write_report(ctx, figures_module):
    paper = getattr(figures_module, "PAPER", "")
    lines = [f"# Rapport de reproduction — [{paper}]", "", PAPERS.get(paper, ""), ""]
    lines += [getattr(figures_module, "INTRO", "").strip(), ""]
    if ctx.quick:
        lines += ["**Mode rapide** : niveaux et durees reduits, comparaisons indicatives seulement.", ""]
    n_ok = sum(1 for r in ctx.results.values() if r.get("status") == "ok")
    n_fail = [n for n, r in ctx.results.items() if r.get("status") != "ok"]
    lines += [f"{n_ok} cas reussis" + (f", {len(n_fail)} en echec : " + ", ".join(n_fail) if n_fail else "") + ".", ""]
    lines += ["| Figure / tableau | Statut |", "|---|---|"]
    for item in sorted(ctx.sections, key=_num_key):
        s = ctx.sections[item]
        lines.append(f"| {item} — {s['title']} | {'; '.join(s['status']) or '-'} |")
    lines.append("")
    for item in sorted(ctx.sections, key=_num_key):
        s = ctx.sections[item]
        lines += [f"## {item} — {s['title']}", ""]
        if s["params"]:
            lines += ["Parametres : " + " ; ".join(s["params"]), ""]
        for n in s["notes"]:
            lines += [n, ""]
        for (head, rows, cap) in s["tables"]:
            if cap:
                lines += [cap, ""]
            lines += [markdown_table(head, rows), ""]
        for im in s["images"]:
            lines += [f"![{item}]({im})", ""]
        if s["status"]:
            lines += ["Statut : " + "; ".join(s["status"]), ""]
    (ctx.out / "rapport.md").write_text("\n".join(lines))
    return ctx.out / "rapport.md"


def make_all(paper_dir, out, figures_module, fmts=("png",)):
    ctx = Context(paper_dir, out, fmts)
    rc = 0
    for fn in figures_module.FIGURES:
        try:
            fn(ctx)
            print(f"  {fn.__name__:32s} ok", flush=True)
        except Exception as e:
            rc = 1
            print(f"  {fn.__name__:32s} ECHEC : {e!r}", flush=True)
            traceback.print_exc()
            doc = (fn.__doc__ or fn.__name__).strip().splitlines()[0]
            m = re.match(r"((?:Fig\.|Table)\s*\d+)", doc)
            if m:
                ctx.section(m.group(1), doc, notes=[f"Erreur du trace : {e!r}"], status="echec du trace")
    rep = write_report(ctx, figures_module)
    print(f"Rapport : {rep}")
    return rc


if __name__ == "__main__":
    paper_dir = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else paper_dir
    sys.path.insert(0, str(paper_dir))
    figures = importlib.import_module("figures")
    sys.exit(make_all(paper_dir, out, figures))
