-- [CMAME25] Tables 2 et 3 : CFL* du probleme couple (Algorithme 1) pour ERK(s), s = 2..4, k = 1..3, ordre egal,
-- stabilisation O(1) de poids eta^f_* = 0.8, eta^s_* = 1.5 ; CFL* = c_p^s dt / h_min (h diametre des cellules).
-- Table 2 : 64 x 32 carres de (-1,1) x (0,1) (niveau 4 du code de l'article, qui reproduit le tableau a 1 %) ;
-- Table 3 : aussi triangles gmsh (lc = 0.035) et polygones PolyMesher (1024 germes par carre unite), symetriques par
-- rapport a x = 0 (maillages regeneres, dossier meshes/). Donnee initiale de la solution sinusoidale (5.2) (omega = 1,
-- theta = 10), Dirichlet homogene, sans source. Lance par run_campaign.py (space.k, time.stages, maillage).

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.8, weight_elastic = 1.5 }
time    = { scheme = "erk", stages = 2, final_time = 1.0, n_steps = 300 }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 4 }
problem = { medium = "coupled", solution = "coupled_sinusoidal", omega = 1, theta = 10 }
study   = { type = "cfl", cfl_threshold = 0.05, cfl_decrement = 0.01 }
