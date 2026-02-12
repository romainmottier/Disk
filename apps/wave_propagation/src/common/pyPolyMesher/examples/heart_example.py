# import numpy as np

# # 1. Define the Bounding Box
# BdBox = [-10, 10, -14, 8]
# #
# # 2. Create a signed distance function

# from pyPolyMesher.dFunctions import dLine, dIntersect, dCircle, dUnion


# def SDF2(p, r=5.5):
#     """
#     Calculate the signed distance from a point P to a heart shape defined by two circles and three line segments.

#     The heart shape consists of two circles representing the left and right sides, and three line segments representing the bottom tip and the tangent lines connecting the circles.

#     Parameters:
#         p (tuple): A tuple containing the coordinates (x, y) of the point P.
#         r (float): The radius of the heart shape. Default is 5.5.

#     Returns:
#         float: The signed distance from the point P to the heart shape. Negative values indicate that the point is inside the heart shape, zero indicates that the point is on the boundary, and positive values indicate that the point is outside the heart shape.

#     Note:
#         For more information please refer to: https://sadjadabedi.ir/post/constructing-complex-shapes-with-signed-distance-functions-the-heart-example/
#     """
#     # Some experimental ratio for the center of circles based on their radius
#     a = r * 3 / 4

#     circle1 = dCircle(p, -a, 0, r)  # Left Circle
#     circle2 = dCircle(p, a, 0, r)  # Right Circle

#     d2c = np.sqrt(
#         (a - 0) ** 2 + (0 - 2 * a - r) ** 2
#     )  # Distance from bottom tip to center of circles
#     d2t = np.sqrt(d2c**2 - r**2)  # Distance to tangent point of circle

#     tpr = find_tangent_point(
#         a, 0, r, 0, -2 * a - r, d2t
#     )  # Tangent point on right circle

#     line1 = dLine(p, -tpr[0], tpr[1], 0, -2 * a - r)
#     line2 = dLine(p, 0, -2 * a - r, tpr[0], tpr[1])
#     line3 = dLine(p, tpr[0], tpr[1], -tpr[0], tpr[1])

#     #    Create a triangle which base is the line that
#     #    connects tangent points and the vertex point is heart bottom tip

#     dl = dIntersect(line1, line2, line3)

#     d = dUnion(circle1, circle2, dl)

#     return d


# def find_tangent_point(x1, y1, d1, x2, y2, d2):
#     """
#     Find the tangent point that has radius(d1) distance from center of circle
#     and calculated distance (d2) from the bottom tip of the heart. The result is
#     limited to desired tangent point.


#     Inspired from following gist:
#     https://gist.github.com/jupdike/bfe5eb23d1c395d8a0a1a4ddd94882ac?permalink_comment_id=4545858#gistcomment-4545858
#     """

#     centerdx = x1 - x2
#     centerdy = y1 - y2
#     R = np.sqrt(centerdx**2 + centerdy**2)

#     if not (abs(d1 - d2) <= R <= d1 + d2):
#         # No intersections
#         return []

#     d1d2 = d1**2 - d2**2
#     a = d1d2 / (2 * R**2)

#     c = np.sqrt(2 * (d1**2 + d2**2) / R**2 - (d1d2**2) / R**4 - 1)

#     fx = (x1 + x2) / 2 + a * (x2 - x1)
#     gx = c * (y2 - y1) / 2
#     # ix1 = fx + gx
#     ix2 = fx - gx

#     fy = (y1 + y2) / 2 + a * (y2 - y1)
#     gy = c * (x1 - x2) / 2
#     # iy1 = fy + gy
#     iy2 = fy - gy

#     return [ix2, iy2]


# # 4. Create Domain object
# from pyPolyMesher import Domain

# NewDomain = Domain("My New Domain", BdBox, SDF2)

# # 5. Plot the domain
# NewDomain.Plot()

# # 6. Generate mesh
# from pyPolyMesher import PolyMesher

# Node, Element, Supp, Load, P = PolyMesher(NewDomain, 100, 300, anim=False)













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
NumberofElements = 5   # nombre d'éléments désiré
MaxIterations = 1000     # nombre max d'itérations pour l'algorithme
Node, Element, Supp, Load, P = PolyMesher(PolygonDomain, NumberofElements, MaxIterations, anim=False)

# ---------------------------
# 7. Évaluer le maillage
# ---------------------------
area = PolygonDomain.CalculateArea()
metrics = mesh_assessment(Node, Element, area, verbose=True)

print("Nodes shape:", Node.shape)
print("Elements shape:", Element.shape)
