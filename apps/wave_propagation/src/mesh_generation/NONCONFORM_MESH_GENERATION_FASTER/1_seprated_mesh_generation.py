import gmsh

gmsh.initialize()
gmsh.model.add("two_nonconformal_domains")

# Paramètres
p = 2              # facteur de raffinement
lc_inf = 0.015625       # taille éléments partie inférieure
lc_sup = lc_inf / p  # taille éléments partie supérieure

xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0
y_mid = abs((ymin + ymax) / 2)

# --- Partie inférieure : rectangle [-0.5,0.5] x [-0.5,0] ---
p_bl = gmsh.model.occ.addPoint(xmin, ymin, 0)
p_br = gmsh.model.occ.addPoint(xmax, ymin, 0)
p_tr = gmsh.model.occ.addPoint(xmax, y_mid, 0)
p_tl = gmsh.model.occ.addPoint(xmin, y_mid, 0)

l_bottom = gmsh.model.occ.addLine(p_bl, p_br)
l_right  = gmsh.model.occ.addLine(p_br, p_tr)
l_top    = gmsh.model.occ.addLine(p_tr, p_tl)
l_left   = gmsh.model.occ.addLine(p_tl, p_bl)

loop_inf = gmsh.model.occ.addCurveLoop([l_bottom, l_right, l_top, l_left])
surf_inf = gmsh.model.occ.addPlaneSurface([loop_inf])

# --- Partie supérieure : rectangle [-0.5,0.5] x [0,0.5] ---
p_bl2 = gmsh.model.occ.addPoint(xmin, y_mid, 0)
p_br2 = gmsh.model.occ.addPoint(xmax, y_mid, 0)
p_tr2 = gmsh.model.occ.addPoint(xmax, ymax, 0)
p_tl2 = gmsh.model.occ.addPoint(xmin, ymax, 0)

l_bottom2 = gmsh.model.occ.addLine(p_bl2, p_br2)
l_right2  = gmsh.model.occ.addLine(p_br2, p_tr2)
l_top2    = gmsh.model.occ.addLine(p_tr2, p_tl2)
l_left2   = gmsh.model.occ.addLine(p_tl2, p_bl2)

loop_sup = gmsh.model.occ.addCurveLoop([l_bottom2, l_right2, l_top2, l_left2])
surf_sup = gmsh.model.occ.addPlaneSurface([loop_sup])

gmsh.model.occ.synchronize()

# --- Maillage partie inférieure : cartésien ---
nx = int((xmax - xmin) / lc_inf)
ny = int((y_mid - ymin) / lc_inf)

gmsh.model.mesh.setTransfiniteCurve(l_bottom, nx+1)
gmsh.model.mesh.setTransfiniteCurve(l_top, nx+1)
gmsh.model.mesh.setTransfiniteCurve(l_left, ny+1)
gmsh.model.mesh.setTransfiniteCurve(l_right, ny+1)

cornerTags_inf = [p_bl, p_br, p_tr, p_tl]
gmsh.model.mesh.setTransfiniteSurface(surf_inf, cornerTags=cornerTags_inf)
gmsh.model.mesh.setRecombine(2, surf_inf)

# --- Maillage partie supérieure : quadrilatères raffinés ---
nx_sup = nx * p
ny_sup = int((ymax - y_mid) / lc_sup)

gmsh.model.mesh.setTransfiniteCurve(l_bottom2, nx_sup+1)
gmsh.model.mesh.setTransfiniteCurve(l_top2, nx_sup+1)
gmsh.model.mesh.setTransfiniteCurve(l_left2, ny_sup+1)
gmsh.model.mesh.setTransfiniteCurve(l_right2, ny_sup+1)

cornerTags_sup = [p_bl2, p_br2, p_tr2, p_tl2]
gmsh.model.mesh.setTransfiniteSurface(surf_sup, cornerTags=cornerTags_sup)
gmsh.model.mesh.setRecombine(2, surf_sup)

# --- Définir les propriétés physiques ---
pg_inf = gmsh.model.addPhysicalGroup(2, [surf_inf])
gmsh.model.setPhysicalName(2, pg_inf, "Material_Lower")

pg_sup = gmsh.model.addPhysicalGroup(2, [surf_sup])
gmsh.model.setPhysicalName(2, pg_sup, "Material_Upper")

# --- Générer le maillage ---
gmsh.model.mesh.generate(2)

# Sauvegarder
gmsh.write("/home/mottie0000/Github/Diskpp/meshes/nonconform_square.msh")
gmsh.write("/home/mottie0000/Github/Diskpp/meshes/nonconform_square.vtk")

# Afficher
# gmsh.fltk.run()
# gmsh.finalize()

print("Maillage généré avec facteur de raffinement p =", p)
