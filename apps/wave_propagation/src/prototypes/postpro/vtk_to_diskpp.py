"""VTK (ASCII UNSTRUCTURED_GRID, e.g. a gmsh export) -> DiSk++ poly_2d mesh.

Writes "<n_points> <n_cells> 4", the points (x y), the cells
(n v1 ... vn, 1-based), the boundary point indices of the left / right /
bottom / top sides of the bounding box (one line each) and, with
--cell-ids, the CellEntityIds of every cell (material / physical group,
one per line), as needed by the heterogeneous (basin, nonconforming
coupling) meshes.

Replaces parser_gmsh_diskpp.py / parser_gmsh_unstructured_simplices.py
(simplex_l7, l9_conv_test_0.005), parser_bassin.py (square2, --cell-ids)
and the 3_parser.py step of mesh_generation/nonconform_square (--cell-ids).

Usage: python3 vtk_to_diskpp.py input.vtk output.txt [--cell-ids]
"""

import argparse
import vtkmodules.all as vtk

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input_vtk")
parser.add_argument("output_txt")
parser.add_argument("--cell-ids", action="store_true", help="append the CellEntityIds of every cell")
args = parser.parse_args()

# Lire le fichier VTK
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.input_vtk)
reader.Update()

# Récupérer les données du maillage
unstructured_grid = reader.GetOutput()

# Lire le fichier VTK
with open(args.input_vtk, "r") as vtk_file:
    lines = vtk_file.readlines()

# Trouver la ligne qui commence par "POINTS"
points_line = None
for line in lines:
    if line.startswith("POINTS"):
        points_line = line
        break

# Si la ligne "POINTS" est trouvée, extraire les coordonnées des points
if points_line:
    num_points = int(points_line.split()[1])
    points_data = [line.split() for line in lines[lines.index(points_line) + 1:lines.index(points_line) + num_points + 1]]

# Trouver la ligne qui commence par "CELLS"
cells_line = None
for line in lines:
    if line.startswith("CELLS"):
        cells_line = line
        break

# Si la ligne "CELLS" est trouvée, extraire les informations des cellules
if cells_line:
    num_cells, _ = map(int, cells_line.split()[1:])  # Ne pas exclure le premier chiffre (nombre de sommets)
    cells_data = [line.split() for line in lines[lines.index(cells_line) + 1:lines.index(cells_line) + num_cells + 1]]

# Propriétés matérielles
if args.cell_ids:
    # Trouver la ligne qui commence par "CELL_DATA"
    cell_data_line = None
    for line in lines:
        if line.startswith("CELL_DATA"):
            cell_data_line = line
            break

    # Si la ligne "CELL_DATA" est trouvée, extraire le nombre de cellules associées aux données
    if cell_data_line:
        num_cells_data = int(cell_data_line.split()[1])

        # Chercher la ligne qui commence par "SCALARS CellEntityIds int 1"
        cell_entity_ids_line = None
        for line in lines:
            if line.startswith("LOOKUP_TABLE default"):
                cell_entity_ids_line = line
                break

        # Si la ligne est trouvée, extraire les valeurs associées à chaque cellule
        if cell_entity_ids_line:
            cell_entity_ids_data = [int(line) for line in lines[lines.index(cell_entity_ids_line) + 1:lines.index(cell_entity_ids_line) + num_cells_data + 1]]

# Ouvrir un fichier de sortie pour les coordonnées des points en écriture
with open(args.output_txt, "w") as output_file:

    # Écrire le nombre de points, le nombre de cellules, et "4" sur la première ligne
    output_file.write(f"{len(points_data)} {num_cells} 4\n")

    for point_data in points_data:
        # Les coordonnées sont dans les trois premières colonnes
        x, y, z = map(float, point_data[:3])
        output_file.write(f"{x} {y}\n")

    for cell_data in cells_data:
        adjusted_cell_data = [cell_data[0]] + [str(int(idx) + 1) for idx in cell_data[1:]]
        # Écrire toutes les données des cellules, y compris le nombre de sommets
        output_file.write(" ".join(adjusted_cell_data) + "\n")

    # Extract the points
    points = unstructured_grid.GetPoints()

    # Calculate the dimensions of the domain
    x_coords = [points.GetPoint(i)[0] for i in range(points.GetNumberOfPoints())]
    y_coords = [points.GetPoint(i)[1] for i in range(points.GetNumberOfPoints())]

    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    print(f"Min X: {min_x}, Max X: {max_x}")
    print(f"Min Y: {min_y}, Max Y: {max_y}")

    # Define lists to store point indices on each side
    indices_left = []
    indices_right = []
    indices_bottom = []
    indices_top = []

    # Define a threshold to determine which points are on the sides
    threshold = 1e-6

    for i in range(points.GetNumberOfPoints()):
        x, y, _ = points.GetPoint(i)
        if abs(x - min_x) < threshold:
            indices_left.append(i+1)
        if abs(x - max_x) < threshold:
            indices_right.append(i+1)
        if abs(y - min_y) < threshold:
            indices_bottom.append(i+1)
        if abs(y - max_y) < threshold:
            indices_top.append(i+1)

    # Numéros des points de la frontière
    output_file.write(" ".join(map(str, indices_left)) + "\n")
    output_file.write(" ".join(map(str, indices_right)) + "\n")
    output_file.write(" ".join(map(str, indices_bottom)) + "\n")
    output_file.write(" ".join(map(str, indices_top)) + "\n")

    if args.cell_ids:
        for entity_id in cell_entity_ids_data:
            output_file.write(f"{entity_id}\n")

print(f"Les coordonnées des points et les informations des cellules ont été écrites dans le fichier '{args.output_txt}'.")
