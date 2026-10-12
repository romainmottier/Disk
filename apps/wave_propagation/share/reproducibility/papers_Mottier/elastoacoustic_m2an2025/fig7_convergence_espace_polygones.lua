-- [M2AN25] Fig. 7 : comme la Fig. 6 sur maillages polygonaux PolyMesher symetriques par rapport a x = 0 (128, 512,
-- 2048 et 8192 cellules par carre unite pour l = 2..5 : les maillages poly_256 ... poly_16384 du code de l'article,
-- regeneres par ../meshes.py). Lance par run_campaign.py (space.k, space.order, space.scaling, time.scheme, maillages).

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.88, weight_elastic = 1.54 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "polygons", level = 1 }
problem = { medium = "coupled", solution = "coupled_poly_in_time" }
study   = { type = "space", last_level = 4, refine_dt = false }
output  = { errors = "energy" }
