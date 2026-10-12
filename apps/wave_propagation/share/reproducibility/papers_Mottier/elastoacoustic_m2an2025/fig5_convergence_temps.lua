-- [M2AN25] Fig. 5 : erreurs en fonction du pas de temps pour la solution polynomiale en espace (6.3a), ordre mixte,
-- k = 4, l = 1 (8 x 4 carres de (-1,1) x (0,1)), stabilisation O(1/h), Delta t = 0.1 x 2^-n. Gauche : SDIRK(s, s+1),
-- s = 1..3 ; droite : ERK(s), s = 2..4. Lance par run_campaign.py (time.scheme, time.stages, niveaux de temps).
-- Erreurs : composantes des normes de l'article (output.errors = "energy"), combinees par plot_mottier.py.

formulation = "first_order"
space   = { k = 4, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.80, weight_elastic = 1.38 }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 160 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1 }
problem = { medium = "coupled", solution = "coupled_poly_in_space" }
study   = { type = "time", time_levels = 6, reference = "exact" }
output  = { errors = "energy" }
