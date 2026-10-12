-- [M2AN25] Fig. 6 : erreurs en fonction du pas d'espace (maillages cartesiens) pour la solution polynomiale en temps
-- (6.3b), n = 8 (Delta t = 0.1 x 2^-8), l = 2..5 (niveaux du code de l'article : 4 x 2^l par 2 x 2^l carres de
-- (-1,1) x (0,1), h = sqrt(2) 2^-(l+1)), k = 1..4. Gauche : ordre egal, O(1), ERK ; centre : ordre mixte, O(1), ERK ;
-- droite : ordre mixte, O(1/h), SDIRK. Lance par run_campaign.py (space.k, space.order, space.scaling, time.scheme).
-- Erreurs : composantes des normes de l'article (output.errors = "energy"), combinees par plot_mottier.py.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.88, weight_elastic = 1.54 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 2 }
problem = { medium = "coupled", solution = "coupled_poly_in_time" }
study   = { type = "space", last_level = 5, refine_dt = false }
output  = { errors = "energy" }
