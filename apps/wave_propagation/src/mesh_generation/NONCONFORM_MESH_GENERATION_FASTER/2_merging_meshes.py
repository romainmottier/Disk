import numpy as np
from collections import defaultdict


# ============================================================
# Lecture VTK (ASCII, UNSTRUCTURED_GRID)
# ============================================================
def read_vtk(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    points = []
    cells = []
    cell_data = []

    i = 0
    while i < len(lines):
        l = lines[i].strip()

        if l.startswith("POINTS"):
            n = int(l.split()[1])
            i += 1
            for _ in range(n):
                points.append([float(x) for x in lines[i].split()])
                i += 1

        elif l.startswith("CELLS"):
            n = int(l.split()[1])
            i += 1
            for _ in range(n):
                vals = list(map(int, lines[i].split()))
                cells.append(vals[1:])
                i += 1

        elif l.startswith("CELL_DATA"):
            n = int(l.split()[1])
            i += 2  # skip SCALARS + LOOKUP
            for _ in range(n):
                while lines[i].strip().startswith("LOOKUP_TABLE"):
                    i += 1
                cell_data.append(int(lines[i].strip()))
                i += 1
        else:
            i += 1

    return np.array(points), cells, cell_data


# ============================================================
# Fusion des points doublons (O(N))
# ============================================================
def merge_duplicate_points(points, tol=1e-12):
    points = np.asarray(points)
    key_map = {}
    new_points = []
    mapping = np.empty(len(points), dtype=int)

    for i, p in enumerate(points):
        key = (round(p[0]/tol), round(p[1]/tol))
        if key in key_map:
            mapping[i] = key_map[key]
        else:
            idx = len(new_points)
            key_map[key] = idx
            mapping[i] = idx
            new_points.append(p)

    return np.array(new_points), mapping


def remap_cells(cells, mapping):
    return [[mapping[i] for i in c] for c in cells]


# ============================================================
# Détection des cellules adjacentes à l'interface
# ============================================================
def find_interface_cells(cells, cell_data):
    """
    Trouve les cellules qui ont au moins un voisin d'un matériau différent.
    On utilise les points partagés pour détecter la proximité.
    """
    # Construire le mapping point -> cellules
    point_to_cells = defaultdict(list)
    for cell_id, cell in enumerate(cells):
        for pt_id in cell:
            point_to_cells[pt_id].append(cell_id)
    
    # Identifier les cellules d'interface
    interface_cells = set()
    for pt_id, cell_ids in point_to_cells.items():
        if len(cell_ids) < 2:
            continue
        # Vérifier si différents matériaux se touchent à ce point
        materials = set(cell_data[c] for c in cell_ids)
        if len(materials) > 1:
            # Ce point est à l'interface, marquer toutes ces cellules
            interface_cells.update(cell_ids)
    
    return interface_cells


# ============================================================
# Détection des hanging nodes (version optimisée)
# ============================================================
def find_hanging_nodes(points, cells, interface_cells=None, tol=1e-12):
    """
    Trouve les hanging nodes sur les arêtes des cellules d'interface.
    Si interface_cells est None, cherche sur toutes les cellules.
    """
    points = np.asarray(points)
    hanging = {}
    
    # Si pas de filtre, traiter toutes les cellules
    if interface_cells is None:
        cells_to_process = range(len(cells))
    else:
        cells_to_process = interface_cells
    
    # Parcourir les cellules d'interface
    for cell_id in cells_to_process:
        cell = cells[cell_id]
        n = len(cell)
        
        for k in range(n):
            i = cell[k]
            j = cell[(k+1) % n]
            
            pi = points[i]
            pj = points[j]
            
            edge_vec = pj - pi
            edge_len2 = np.dot(edge_vec, edge_vec)
            
            if edge_len2 < 1e-16:
                continue
            
            # Chercher les points sur cette arête
            midpoints = []
            for idx, p in enumerate(points):
                if idx == i or idx == j:
                    continue
                
                vec = p - pi
                t = np.dot(vec, edge_vec) / edge_len2
                
                if tol < t < 1 - tol:
                    proj = pi + t * edge_vec
                    if np.linalg.norm(proj - p) < tol:
                        midpoints.append((t, idx))
            
            if midpoints:
                # Trier par position le long de l'arête
                midpoints.sort(key=lambda x: x[0])
                edge_key = tuple(sorted((i, j)))
                
                # Garder uniquement la liste la plus longue si l'arête existe déjà
                new_nodes = [idx for _, idx in midpoints]
                if edge_key not in hanging or len(new_nodes) > len(hanging[edge_key]):
                    hanging[edge_key] = new_nodes
    
    return hanging


# ============================================================
# Découpe des cellules
# ============================================================
def split_cells(cells, hanging, interface_cells=None):
    """
    Découpe les cellules en ajoutant les hanging nodes.
    Si interface_cells est fourni, ne traite que ces cellules.
    """
    new_cells = []
    
    for cell_id, cell in enumerate(cells):
        # Filtrer si nécessaire
        if interface_cells is not None and cell_id not in interface_cells:
            new_cells.append(cell)
            continue
        
        # Traiter cette cellule
        nc = []
        n = len(cell)
        for k in range(n):
            i = cell[k]
            j = cell[(k+1) % n]
            nc.append(i)
            edge_key = tuple(sorted((i, j)))
            if edge_key in hanging:
                nc.extend(hanging[edge_key])
        
        new_cells.append(nc)
    
    return new_cells


# ============================================================
# Écriture VTK polygonal
# ============================================================
def write_vtk_poly(points, cells, cell_data, filename):
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("nonconformal processed\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        f.write(f"POINTS {len(points)} double\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2] if len(p) > 2 else 0.0}\n")

        total = sum(len(c)+1 for c in cells)
        f.write(f"CELLS {len(cells)} {total}\n")
        for c in cells:
            f.write(f"{len(c)} {' '.join(map(str,c))}\n")

        f.write(f"CELL_TYPES {len(cells)}\n")
        f.write(("7\n")*len(cells))  # VTK_POLYGON

        f.write(f"CELL_DATA {len(cell_data)}\n")
        f.write("SCALARS CellEntityIds int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for v in cell_data:
            f.write(f"{v}\n")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":

    input_vtk  = "/home/mottie0000/Github/Diskpp/meshes/nonconform_square.vtk"
    output_vtk = "/home/mottie0000/Github/Diskpp/meshes/nonconform_square.vtk"

    print("📖 Lecture du maillage...")
    points, cells, cell_data = read_vtk(input_vtk)

    print("🔗 Fusion des points doublons...")
    points, mapping = merge_duplicate_points(points)
    cells = remap_cells(cells, mapping)

    print("🔍 Détection des cellules d'interface...")
    interface_cells = find_interface_cells(cells, cell_data)
    print(f"   → {len(interface_cells)} cellules touchent l'interface")

    print("🎯 Recherche des hanging nodes...")
    hanging = find_hanging_nodes(points, cells, interface_cells)
    print(f"   → {len(hanging)} arêtes avec hanging nodes")

    print("✂️  Découpe des cellules d'interface...")
    new_cells = split_cells(cells, hanging, interface_cells)

    print("💾 Écriture du fichier de sortie...")
    write_vtk_poly(points, new_cells, cell_data, output_vtk)

    print("\n✅ Maillage non conforme traité avec succès !")
    print(f"   Points: {len(points)}")
    print(f"   Cells: {len(new_cells)}")
    print(f"   Cellules d'interface: {len(interface_cells)}")
    print(f"   Arêtes avec hanging nodes: {len(hanging)}")
