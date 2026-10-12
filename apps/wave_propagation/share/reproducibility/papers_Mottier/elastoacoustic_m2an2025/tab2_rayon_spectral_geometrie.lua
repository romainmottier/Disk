-- [M2AN25] Table 2 : rayon spectral (normalise par sqrt(#cellules)) du probleme couple selon la forme des cellules
-- (simplexes, quadrangles, polygones PolyMesher), ordres egal et mixte, k = 1..6, poids eta_* de (6.2). Maillages
-- "l1" du code de l'article : 8 x 4 carres, triangles gmsh de taille 0.35 et polygones PolyMesher de 32 cellules par
-- carre unite, symetriques par rapport a x = 0 (les fichiers de l'article ne sont pas disponibles : maillages
-- regeneres par ../meshes.py, dossier meshes/). Lance par run_campaign.py (space.k, space.order, maillage).

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.88, weight_elastic = 1.54 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 10 }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1 }
problem = { medium = "coupled", solution = "coupled_poly_in_time" }
study   = { type = "spectrum" }
