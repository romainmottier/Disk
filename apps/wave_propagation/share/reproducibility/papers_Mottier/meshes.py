#!/usr/bin/env python3
"""
meshes.py — maillages des campagnes [M2AN25] et [CMAME25] au format "polygons" de wave_propagation
(mesh.type = "polygons", lu par tagged_polygon_mesh_reader, src/mesh/domain.hpp) :

    # commentaires
    n_points n_cells
    x y                                   (n_points lignes)
    nv i_1 ... i_nv tag                   (n_cells lignes, sommets numerotes a partir de 0, sens trigonometrique)

Les fichiers de maillage des articles (PolyMesher, gmsh, bassin) ne sont pas disponibles : ce module les regenere.
  polymesher   maillage polygonal de Voronoi centroidal du rectangle, algorithme de PolyMesher (Talischi, Paulino,
               Pereira, Menezes 2012 : semis aleatoire, iterations de Lloyd avec reflexion des germes proches du bord,
               puis suppression des petites aretes, tolerance 0.1), germe aleatoire fixe
  triangles    maillage triangulaire non structure de gmsh (Frontal-Delaunay) du rectangle, taille caracteristique lc
  mirror       maillage du carre (0,1)^2 et son symetrique par rapport a x = 0 : maillage conforme de (-1,1) x (0,1)
               (partie elastique x < 0, etiquette 1 ; partie acoustique x > 0, etiquette 2), comme les maillages
               polygonaux et simpliciaux des figures des articles
  vtk          conversion d'un maillage gmsh au format VTK (bassin de [CMAME25]) : etiquettes = CellEntityIds
  refine       raffinement non conforme d'une partie des cellules (etiquettes donnees) : quadrangles et triangles
               coupes en 4, noeuds pendants inseres dans les cellules voisines (qui deviennent des polygones)
  plot         dessin d'un maillage (couleur par etiquette ou par milieu)

Utilisation (les campagnes appellent ces fonctions ; en ligne de commande pour les essais) :
    python3 meshes.py polymesher N fichier [--seed S] [--mirror]
    python3 meshes.py triangles LC fichier [--mirror]
    python3 meshes.py vtk fichier.vtk fichier
    python3 meshes.py refine fichier fichier_raffine --tags 2 3
    python3 meshes.py plot fichier image.png
"""
import argparse, bz2, sys
from pathlib import Path

import numpy as np


# ------------------------------------------------------------------------------------------------------------------
# lecture / ecriture

def _open(path, mode):
    path = str(path)
    return bz2.open(path, mode + "t") if path.endswith(".bz2") else open(path, mode)


def write_polygons(path, pts, cells, tags, comment=""):
    """Ecrit un maillage "polygons" (cellules orientees dans le sens trigonometrique ; .bz2 : comprime)"""
    with _open(path, "w") as f:
        for line in comment.strip().splitlines():
            f.write("# " + line + "\n")
        f.write(f"{len(pts)} {len(cells)}\n")
        for p in pts:
            f.write(f"{float(p[0])!r} {float(p[1])!r}\n")      # plus courte ecriture exacte (aller-retour)
        for c, t in zip(cells, tags):
            f.write(f"{len(c)} " + " ".join(str(int(i)) for i in c) + f" {int(t)}\n")


def read_polygons(path):
    with _open(path, "r") as f:
        lines = [l for l in f.read().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    n_pts, n_cells = map(int, lines[0].split())
    pts = np.array([[float(v) for v in l.split()[:2]] for l in lines[1:1 + n_pts]])
    cells, tags = [], []
    for l in lines[1 + n_pts:1 + n_pts + n_cells]:
        v = [int(x) for x in l.split()]
        cells.append(v[1:1 + v[0]])
        tags.append(v[1 + v[0]])
    return pts, cells, np.array(tags)


def signed_area(poly):
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def centroid(poly):
    x, y = poly[:, 0], poly[:, 1]
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    a = 0.5 * cr.sum()
    return np.array([((x + np.roll(x, -1)) * cr).sum(), ((y + np.roll(y, -1)) * cr).sum()]) / (6.0 * a), a


def orient_ccw(pts, cells):
    return [list(c) if signed_area(pts[c]) > 0 else list(c)[::-1] for c in cells]


def merge_nodes(pts, cells, tol):
    """Fusionne les noeuds a distance < tol (union-find), supprime les sommets repetes et les noeuds inutilises"""
    from scipy.spatial import cKDTree
    parent = np.arange(len(pts))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in cKDTree(pts).query_pairs(tol):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    rep = np.array([find(i) for i in range(len(pts))])
    return compact(pts, [[rep[i] for i in c] for c in cells])


def compact(pts, cells):
    """Supprime les sommets consecutifs repetes des cellules et renumerote les noeuds utilises"""
    out = []
    for c in cells:
        d = [v for k, v in enumerate(c) if v != c[k - 1]] if len(c) > 1 else list(c)
        out.append(d)
    used = sorted({v for c in out for v in c})
    new = {v: k for k, v in enumerate(used)}
    return pts[used].copy(), [[new[v] for v in c] for c in out]


# ------------------------------------------------------------------------------------------------------------------
# PolyMesher (Talischi et al., Struct. Multidisc. Optim. 45 (2012) 309-328) sur un rectangle

def _box_dist(P, box):
    x0, x1, y0, y1 = box
    d = np.stack([x0 - P[:, 0], P[:, 0] - x1, y0 - P[:, 1], P[:, 1] - y1], axis=1)
    return np.concatenate([d, d.max(axis=1, keepdims=True)], axis=1)


_NORMALS = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]])


def _reflect(P, box, alpha, eta=0.9):
    """PolyMshr_Rflct : reflexion des germes a distance < alpha d'un cote, gardees si assez loin a l'exterieur"""
    d = _box_dist(P, box)
    R = []
    for s in range(4):
        I = np.abs(d[:, s]) < alpha
        R.append(P[I] - 2.0 * _NORMALS[s][None, :] * d[I, s][:, None])
        dR = _box_dist(R[-1], box)[:, -1]
        R[-1] = R[-1][(np.abs(dR) >= eta * np.abs(d[I, s])) & (dR > 0)]
    R = np.concatenate(R, axis=0)
    return np.unique(R, axis=0) if len(R) else R


def _voronoi_cells(P, R):
    from scipy.spatial import Voronoi
    vor = Voronoi(np.concatenate([P, R], axis=0))
    cells = []
    for i in range(len(P)):
        reg = vor.regions[vor.point_region[i]]
        if not reg or -1 in reg:
            raise RuntimeError("PolyMesher : cellule de Voronoi non bornee (germe %d)" % i)
        cells.append(list(reg))
    return vor.vertices, cells


def _snap(pts, box, tol):
    x0, x1, y0, y1 = box
    for j, (a, b) in enumerate([(x0, x1), (y0, y1)]):
        for v in (a, b):
            pts[np.abs(pts[:, j] - v) < tol, j] = v
    return pts


def _boundary_rank(p, box, tol):
    x0, x1, y0, y1 = box
    on = (abs(p[0] - x0) < tol) + (abs(p[0] - x1) < tol) + (abs(p[1] - y0) < tol) + (abs(p[1] - y1) < tol)
    return int(on)          # 0 interieur, 1 bord, 2 coin


def _collapse_edges(pts, cells, box, tol=0.1):
    """PolyMshr_CllpsEdgs : contracte les aretes vues sous un angle < tol * 2 pi / nv depuis le centre de la cellule
    (cellules d'au moins 4 sommets) ; le noeud garde est un coin, sinon un noeud du bord, sinon le premier"""
    size_tol = 1e-9 * max(box[1] - box[0], box[3] - box[2])
    while True:
        cedges = set()
        for c in cells:
            if len(c) < 4:
                continue
            v = pts[c]
            m = v.mean(axis=0)
            beta = np.arctan2(v[:, 1] - m[1], v[:, 0] - m[0])
            dbeta = np.mod(np.roll(beta, -1) - beta, 2.0 * np.pi)
            ideal = 2.0 * np.pi / len(c)
            for k in np.nonzero(dbeta < tol * ideal)[0]:
                a, b = c[k], c[(k + 1) % len(c)]
                cedges.add((min(a, b), max(a, b)))
        if not cedges:
            return pts, cells
        rep = np.arange(len(pts))

        def find(i):
            while rep[i] != i:
                i = rep[i]
            return i
        for a, b in sorted(cedges):
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            ka, kb = _boundary_rank(pts[ra], box, size_tol), _boundary_rank(pts[rb], box, size_tol)
            keep, drop = (ra, rb) if ka >= kb else (rb, ra)
            rep[drop] = keep
        new_cells = [[find(i) for i in c] for c in cells]
        pts, cells = compact(pts, new_cells)


def polymesher(n, box=(0.0, 1.0, 0.0, 1.0), seed=0, max_iter=100, tol=5e-6, c=1.5, collapse=0.1):
    """Maillage polygonal de PolyMesher du rectangle box = (x0, x1, y0, y1) a n cellules : (points, cellules)"""
    x0, x1, y0, y1 = box
    rng = np.random.default_rng(seed)
    P = np.stack([x0 + (x1 - x0) * rng.random(n), y0 + (y1 - y0) * rng.random(n)], axis=1)
    area = (x1 - x0) * (y1 - y0)
    it, err = 0, 1.0
    while it <= max_iter and err > tol:
        alpha = c * np.sqrt(area / n)
        R = _reflect(P, box, alpha)
        V, cells = _voronoi_cells(P, R)
        cent = np.zeros_like(P)
        A = np.zeros(n)
        for i, cl in enumerate(cells):
            cent[i], A[i] = centroid(V[cl])
        area = np.abs(A).sum()
        err = np.sqrt(np.sum(A ** 2 * np.sum((cent - P) ** 2, axis=1))) * n / area ** 1.5
        P = cent
        it += 1
    R = _reflect(P, box, c * np.sqrt(area / n))
    V, cells = _voronoi_cells(P, R)
    scale = max(x1 - x0, y1 - y0)
    pts, cells = merge_nodes(V.copy(), cells, 1e-9 * scale)
    pts = _snap(pts, box, 1e-9 * scale)
    pts, cells = _collapse_edges(pts, cells, box, collapse)
    return pts, orient_ccw(pts, cells)


# ------------------------------------------------------------------------------------------------------------------
# gmsh

def triangles(lc, box=(0.0, 1.0, 0.0, 1.0), algorithm=6):
    """Maillage triangulaire de gmsh du rectangle, taille caracteristique lc (Frontal-Delaunay par defaut)"""
    import gmsh
    x0, x1, y0, y1 = box
    gmsh.initialize(["-noenv"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("rectangle")
        g = gmsh.model.geo
        p = [g.addPoint(x0, y0, 0, lc), g.addPoint(x1, y0, 0, lc), g.addPoint(x1, y1, 0, lc), g.addPoint(x0, y1, 0, lc)]
        l = [g.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
        g.addPlaneSurface([g.addCurveLoop(l)])
        g.synchronize()
        gmsh.option.setNumber("Mesh.Algorithm", algorithm)
        gmsh.option.setNumber("Mesh.RandomSeed", 1)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        idx = {int(t): k for k, t in enumerate(tags)}
        pts = coords.reshape(-1, 3)[:, :2].copy()
        cells = []
        for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2)):
            nn = gmsh.model.mesh.getElementProperties(etype)[3]
            for e in enodes.reshape(-1, nn):
                cells.append([idx[int(t)] for t in e[:3]])
    finally:
        gmsh.finalize()
    pts, cells = compact(pts, cells)
    pts = _snap(pts, box, 1e-12 * max(x1 - x0, y1 - y0))
    return pts, orient_ccw(pts, cells)


# ------------------------------------------------------------------------------------------------------------------
# maillage couple par symetrie

def mirror(pts, cells):
    """Maillage de (0,1)^2 -> maillage de (-1,1) x (0,1) : symetrique par rapport a x = 0 (x < 0 : etiquette 1,
    milieu elastique ; x > 0 : etiquette 2, milieu acoustique), noeuds de x = 0 communs"""
    n = len(pts)
    left = pts.copy()
    left[:, 0] = -left[:, 0]
    on_axis = np.nonzero(pts[:, 0] == 0.0)[0]
    new_pts = np.concatenate([left, pts[[i for i in range(n) if pts[i, 0] != 0.0]]], axis=0)
    right_index = np.empty(n, dtype=int)
    right_index[on_axis] = on_axis
    k = n
    for i in range(n):
        if pts[i, 0] != 0.0:
            right_index[i] = k
            k += 1
    left_cells = [list(c)[::-1] for c in cells]                 # la symetrie inverse l'orientation
    right_cells = [[int(right_index[v]) for v in c] for c in cells]
    tags = [1] * len(cells) + [2] * len(cells)
    return new_pts, left_cells + right_cells, tags


# ------------------------------------------------------------------------------------------------------------------
# VTK (gmsh) et raffinement non conforme

def read_vtk(path):
    """Maillage VTK legacy ASCII (UNSTRUCTURED_GRID) : points 2D, cellules (triangles, quadrangles), CellEntityIds"""
    with open(path) as f:
        tok = f.read().split()
    i = tok.index("POINTS")
    npts = int(tok[i + 1]); i += 3
    pts = np.array(tok[i:i + 3 * npts], float).reshape(npts, 3)[:, :2]
    i = tok.index("CELLS")
    ncells = int(tok[i + 1]); i += 3
    cells = []
    for _ in range(ncells):
        nv = int(tok[i])
        cells.append([int(v) for v in tok[i + 1:i + 1 + nv]])
        i += nv + 1
    i = tok.index("CELL_TYPES") + 2
    types = [int(v) for v in tok[i:i + ncells]]
    i = tok.index("CellEntityIds") + 5           # SCALARS CellEntityIds int 1 LOOKUP_TABLE default
    ids = np.array([int(v) for v in tok[i:i + ncells]])
    keep = [k for k, t in enumerate(types) if t in (5, 9, 7)]   # triangles, quadrangles, polygones
    pts, cells = compact(pts, [cells[k] for k in keep])
    return pts, orient_ccw(pts, cells), ids[keep]


def refine(pts, cells, tags, refine_tags):
    """Raffinement uniforme (une passe) des cellules d'etiquette dans refine_tags : triangle -> 4 triangles,
    quadrangle -> 4 quadrangles (milieux des aretes et centre) ; les milieux des aretes partagees avec une cellule non
    raffinee sont des noeuds pendants inseres dans celle-ci. Cellules non raffinees a plus de 4 sommets : erreur."""
    pts = [tuple(p) for p in pts]
    mid = {}

    def midpoint(a, b):
        key = (min(a, b), max(a, b))
        if key not in mid:
            pa, pb = pts[a], pts[b]
            pts.append((0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])))
            mid[key] = len(pts) - 1
        return mid[key]
    out_cells, out_tags = [], []
    refined = [t in refine_tags for t in tags]
    for c, t, r in zip(cells, tags, refined):
        if not r:
            continue
        m = [midpoint(c[k], c[(k + 1) % len(c)]) for k in range(len(c))]
        if len(c) == 3:
            out_cells += [[c[0], m[0], m[2]], [m[0], c[1], m[1]], [m[2], m[1], c[2]], [m[0], m[1], m[2]]]
        elif len(c) == 4:
            q = np.array([pts[v] for v in c]).mean(axis=0)
            pts.append((q[0], q[1]))
            ce = len(pts) - 1
            out_cells += [[c[0], m[0], ce, m[3]], [m[0], c[1], m[1], ce], [ce, m[1], c[2], m[2]], [m[3], ce, m[2], c[3]]]
        else:
            raise RuntimeError("refine : cellule a %d sommets" % len(c))
        out_tags += [t] * 4
    for c, t, r in zip(cells, tags, refined):
        if r:
            continue
        poly = []
        for k in range(len(c)):
            a, b = c[k], c[(k + 1) % len(c)]
            poly.append(a)
            key = (min(a, b), max(a, b))
            if key in mid:
                poly.append(mid[key])           # noeud pendant
        out_cells.append(poly)
        out_tags.append(t)
    return np.array(pts), out_cells, np.array(out_tags)


# ------------------------------------------------------------------------------------------------------------------
# dessin

def plot(pts, cells, tags, out, colors=None, lw=0.6, figsize=(6, 3), xlim=None, ylim=None, title=None, dpi=200):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    colors = colors or {1: "#b00000", 2: "#0000ff", 3: "#c8b400"}
    fig, ax = plt.subplots(figsize=figsize)
    for t in sorted(set(tags)):
        polys = [pts[c] for c, tt in zip(cells, tags) if tt == t]
        ax.add_collection(PolyCollection(polys, facecolors="none", edgecolors=colors.get(t, "k"), linewidths=lw))
    ax.set_aspect("equal")
    ax.set_xlim(*(xlim or (pts[:, 0].min(), pts[:, 0].max())))
    ax.set_ylim(*(ylim or (pts[:, 1].min(), pts[:, 1].max())))
    ax.axis("off")
    if title:
        ax.set_title(title)
    fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("polymesher"); a.add_argument("n", type=int); a.add_argument("out")
    a.add_argument("--seed", type=int, default=0); a.add_argument("--mirror", action="store_true")
    a = sub.add_parser("triangles"); a.add_argument("lc", type=float); a.add_argument("out")
    a.add_argument("--mirror", action="store_true")
    a = sub.add_parser("vtk"); a.add_argument("vtk"); a.add_argument("out")
    a = sub.add_parser("refine"); a.add_argument("mesh"); a.add_argument("out")
    a.add_argument("--tags", type=int, nargs="+", required=True)
    a = sub.add_parser("plot"); a.add_argument("mesh"); a.add_argument("out")
    args = ap.parse_args()
    if args.cmd in ("polymesher", "triangles"):
        if args.cmd == "polymesher":
            pts, cells = polymesher(args.n, seed=args.seed)
            what = f"PolyMesher, {args.n} cellules, germe {args.seed}"
        else:
            pts, cells = triangles(args.lc)
            what = f"gmsh, lc = {args.lc}"
        tags = [2] * len(cells)
        if args.mirror:
            pts, cells, tags = mirror(pts, cells)
        write_polygons(args.out, pts, cells, tags, what)
    elif args.cmd == "vtk":
        pts, cells, tags = read_vtk(args.vtk)
        write_polygons(args.out, pts, cells, tags, f"converti de {Path(args.vtk).name}")
    elif args.cmd == "refine":
        pts, cells, tags = read_polygons(args.mesh)
        pts, cells, tags = refine(pts, cells, tags, set(args.tags))
        write_polygons(args.out, pts, cells, tags, f"{Path(args.mesh).name} raffine dans {args.tags}")
    elif args.cmd == "plot":
        pts, cells, tags = read_polygons(args.mesh)
        plot(pts, cells, tags, args.out)
    print(f"{args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
