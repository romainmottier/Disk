-- [M2AN25] Table 1 : rayon spectral (normalise par sqrt(#cellules)) du probleme couple, ordres egal et mixte, k = 1..3,
-- poids eta^f = 2^w eta^f_*, eta^s = 2^w eta^s_*, w = -3..3, eta_* de (6.2) : 0.88 / 1.54 (ordre egal), 0.80 / 1.38
-- (ordre mixte). Meme maillage que la Fig. 3 (8 x 4 carres de (-1,1) x (0,1)). Lance par run_campaign.py (space.k,
-- space.order, poids).

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.88, weight_elastic = 1.54 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 10 }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1 }
problem = { medium = "coupled", solution = "coupled_poly_in_time" }
study   = { type = "spectrum" }
