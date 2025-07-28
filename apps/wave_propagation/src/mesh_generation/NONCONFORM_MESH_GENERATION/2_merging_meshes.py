import numpy as np

def read_vtk(filename):
    """Lit un fichier VTK ASCII unstructured_grid et retourne points, cells, cell_data"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    points = []
    cells = []
    cell_data = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("POINTS"):
            n_points = int(line.split()[1])
            i += 1
            for _ in range(n_points):
                pts = [float(x) for x in lines[i].strip().split()]
                points.append(pts)
                i += 1
        elif line.startswith("CELLS"):
            n_cells = int(line.split()[1])
            i += 1
            for _ in range(n_cells):
                vals = [int(x) for x in lines[i].strip().split()]
                n = vals[0]
                cells.append(vals[1:])
                i += 1
        elif line.startswith("CELL_DATA"):
            n_cd = int(line.split()[1])
            i += 1
            # Skip SCALARS and LOOKUP_TABLE
            i += 2
            for _ in range(n_cd):
                while lines[i].strip() == "" or lines[i].strip().startswith("LOOKUP_TABLE"):
                    i += 1
                cell_data.append(int(lines[i].strip()))
                i += 1
        else:
            i += 1

    return np.array(points), cells, cell_data


def merge_duplicate_points(points, tol=1e-12):
    """Supprime les points en double et retourne un tableau unique + mapping ancien->nouveau"""
    points = np.array(points)
    unique_pts = []
    mapping = np.full(len(points), -1, dtype=int)

    for idx, p in enumerate(points):
        found = False
        for u_idx, u in enumerate(unique_pts):
            if np.linalg.norm(p - u) < tol:
                mapping[idx] = u_idx
                found = True
                break
        if not found:
            mapping[idx] = len(unique_pts)
            unique_pts.append(p)

    return np.array(unique_pts), mapping


def remap_cells(cells, mapping):
    """Réécrit les indices des cellules selon le nouveau tableau de points unique"""
    new_cells = []
    for c in cells:
        new_cells.append([mapping[i] for i in c])
    return new_cells


def find_hanging_nodes(points, cells, tol=1e-12):
    """Trouve les noeuds situés sur une arête (hanging nodes)"""
    hanging_nodes = dict()  # clé=(i,j) indices des arêtes, valeur = [indices noeuds intermédiaires]
    points = np.array(points)

    for cell in cells:
        n = len(cell)
        for k in range(n):
            i = cell[k]
            j = cell[(k + 1) % n]
            pi = points[i]
            pj = points[j]

            edge_vec = pj - pi
            edge_len2 = np.dot(edge_vec, edge_vec)
            if edge_len2 < 1e-16:
                continue

            for idx, p in enumerate(points):
                if idx in cell:
                    continue
                vec = p - pi
                t = np.dot(vec, edge_vec) / edge_len2
                if tol < t < 1 - tol:
                    proj = pi + t * edge_vec
                    if np.linalg.norm(proj - p) < tol:
                        key = tuple(sorted((i, j)))
                        if key not in hanging_nodes:
                            hanging_nodes[key] = []
                        hanging_nodes[key].append(idx)

    # Trier les hanging nodes le long de l’arête
    for edge, nodes in hanging_nodes.items():
        i, j = edge
        vec = points[j] - points[i]
        hanging_nodes[edge] = sorted(nodes, key=lambda idx: np.dot(points[idx]-points[i], vec))

    return hanging_nodes


def split_edges(cells, hanging_nodes):
    """Découpe les arêtes avec des hanging nodes et crée de nouveaux polygones"""
    new_cells = []

    for cell in cells:
        n = len(cell)
        new_cell = []
        for k in range(n):
            i = cell[k]
            j = cell[(k + 1) % n]
            new_cell.append(i)
            key = tuple(sorted((i, j)))
            if key in hanging_nodes:
                new_cell.extend(hanging_nodes[key])
        new_cells.append(new_cell)

    return new_cells


def write_vtk_poly(points, cells, cell_data, filename):
    """Écrit un fichier VTK ASCII avec des polygones"""
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("nonconformal_square_poly, processed\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} double\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2] if len(p) > 2 else 0.0}\n")

        n_total = sum(len(c) + 1 for c in cells)
        f.write(f"CELLS {len(cells)} {n_total}\n")
        for c in cells:
            f.write(f"{len(c)} {' '.join(map(str, c))}\n")

        f.write(f"CELL_TYPES {len(cells)}\n")
        for _ in cells:
            f.write("7\n")  # VTK_POLYGON

        f.write(f"CELL_DATA {len(cells)}\n")
        f.write("SCALARS CellEntityIds int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for cd in cell_data:
            f.write(f"{cd}\n")


if __name__ == "__main__":
    input_vtk = "/home/mottie0000/Github/Diskpp/meshes/nonconform_square.vtk"
    output_vtk = "/home/mottie0000/Github/Diskpp/meshes/nonconform_square.vtk"

    # --- Lire le fichier VTK existant ---
    points, cells, cell_data = read_vtk(input_vtk)

    # --- Fusionner les points doublons pour éviter artefacts ---
    points, mapping = merge_duplicate_points(points)
    cells = remap_cells(cells, mapping)

    # --- Trouver les hanging nodes ---
    hanging_nodes = find_hanging_nodes(points, cells)

    # --- Créer de nouveaux polygones en intégrant les hanging nodes ---
    new_cells = split_edges(cells, hanging_nodes)

    # --- Écrire le nouveau fichier ---
    write_vtk_poly(points, new_cells, cell_data, output_vtk)

    print("✅ VTK polygonal généré sans artefacts.")
    print(f"Points: {points.shape}, Cells: {len(new_cells)}, Cell data: {len(cell_data)}")
