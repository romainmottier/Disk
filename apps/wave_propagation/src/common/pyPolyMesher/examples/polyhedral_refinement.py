import numpy as np
from pyPolyMesher import PolyMesher, Domain, mesh_assessment
from pyPolyMesher.dFunctions import dPolygon

# ---------------------------
# 1. Définir le contour du domaine (polygonal)
# ---------------------------
polygon_points = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [1.2, 0.5],
    [1.0, 1.0],
    [0.0, 1.0]
])

# Fermer le contour si nécessaire (converti en liste de tuples)
polygon_list = [tuple(pt) for pt in polygon_points]
if polygon_list[0] != polygon_list[-1]:
    polygon_list.append(polygon_list[0])

# ---------------------------
# 2. Définir la bounding box
# ---------------------------
x_min, y_min = np.min(polygon_points, axis=0)
x_max, y_max = np.max(polygon_points, axis=0)
BdBox = [x_min, x_max, y_min, y_max]

# ---------------------------
# 3. Créer la fonction de distance signée (SDF)
# ---------------------------
SDF = lambda P: dPolygon(P, polygon_list)

# ---------------------------
# 4. Créer l'objet Domain
# ---------------------------
PolygonDomain = Domain("Polygon Domain", BdBox, SDF)

# ---------------------------
# 5. Plot du domaine
# ---------------------------
PolygonDomain.Plot()

# ---------------------------
# 6. Générer le maillage
# ---------------------------
NumberofElements = 10   # nombre d'éléments désiré
MaxIterations = 1000     # nombre max d'itérations pour l'algorithme
Node, Element, Supp, Load, P = PolyMesher(PolygonDomain, NumberofElements, MaxIterations, anim=False)

# ---------------------------
# 7. Évaluer le maillage
# ---------------------------
area = PolygonDomain.CalculateArea()
metrics = mesh_assessment(Node, Element, area, verbose=True)

print("Nodes shape:", Node.shape)
print("Elements shape:", Element.shape)

# import numpy as np
# from pyPolyMesher import PolyMesher, Domain
# from pyPolyMesher.dFunctions import dPolygon


# def generate_mesh(polygon_points, NumberofElements, MaxIterations):
#     polygon_points = np.array(polygon_points)

#     # fermer le contour
#     polygon_list = [tuple(pt) for pt in polygon_points]
#     if polygon_list[0] != polygon_list[-1]:
#         polygon_list.append(polygon_list[0])

#     # bounding box
#     x_min, y_min = np.min(polygon_points, axis=0)
#     x_max, y_max = np.max(polygon_points, axis=0)
#     BdBox = [x_min, x_max, y_min, y_max]

#     SDF = lambda P: dPolygon(P, polygon_list)

#     domain = Domain("Polygon Domain", BdBox, SDF)

#     Node, Element, Supp, Load, P = PolyMesher(
#         domain,
#         NumberofElements,
#         MaxIterations,
#         anim=False
#     )

#     return Node, Element
