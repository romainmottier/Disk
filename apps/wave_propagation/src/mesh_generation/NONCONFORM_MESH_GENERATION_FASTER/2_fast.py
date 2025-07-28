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
# Hanging nodes — version O(N) (interface horizontale)
# ============================================================
def find_hanging_nodes(points, cells, y_mid, tol=1e-12):
    points = np.asarray(points)

    # 1. index points on interface y=y_mid
    interface_pts = [
        i for i,p in enumerate(points)
        if abs(p[1] - y_mid) < tol
    ]

    interface_pts.sort(key=lambda i: points[i][0])
    interface_x = np.array([points[i][0] for i in interface_pts])

    hanging = {}

    # 2. inspect edges lying on the interface
    for cell in cells:
        n = len(cell)
        for k in range(n):
            i = cell[k]
            j = cell[(k+1)%n]

            pi, pj = points[i], points[j]

            if abs(pi[1]-y_mid) < tol and abs(pj[1]-y_mid) < tol:
                x0, x1 = sorted([pi[0], pj[0]])
                if x1 - x0 < tol:
                    continue

                mask = (interface_x > x0+tol) & (interface_x < x1-tol)
                mids = [interface_pts[m] for m in np.where(mask)[0]]

                if mids:
                    hanging[tuple(sorted((i,j)))] = mids

    return hanging


# ============================================================
# Découpe des cellules
# ============================================================
def split_cells(cells, hanging):
    new_cells = []

    for cell in cells:
        nc = []
        n = len(cell)
        for k in range(n):
            i = cell[k]
            j = cell[(k+1)%n]
            nc.append(i)
            key = tuple(sorted((i,j)))
            if key in hanging:
                nc.extend(hanging[key])
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
            f.write(f"{p[0]} {p[1]} 0.0\n")

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
    output_vtk = "/home/mottie0000/Github/Diskpp/meshes/nonconform_square_poly.vtk"

    y_mid = 0.5   # INTERFACE (doit matcher le script gmsh)

    points, cells, cell_data = read_vtk(input_vtk)

    points, mapping = merge_duplicate_points(points)
    cells = remap_cells(cells, mapping)

    hanging = find_hanging_nodes(points, cells, y_mid)

    new_cells = split_cells(cells, hanging)

    write_vtk_poly(points, new_cells, cell_data, output_vtk)

    print("✅ Maillage non conforme traité")
    print(f"Points: {len(points)}")
    print(f"Cells: {len(new_cells)}")
    print(f"Hanging edges: {len(hanging)}")
