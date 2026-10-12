-- [CMAME25] Mottier, Ern, Guillot, "Elasto-acoustic wave propagation in geophysical media using hybrid high-order
-- methods on general meshes", CMAME (2025), section 5.1, Fig. 6 : CFL* des cas purs (acoustique, elastique) en ordre
-- egal en fonction du poids eta de la stabilisation O(1), k = 1..3, par l'algorithme 1 : ERK(4), donnee initiale de la
-- solution sinusoidale (5.2) (omega = 1, theta = 10), Dirichlet homogene, sans source, seuil 5 % de croissance de
-- l'energie, pas de temps diminues de 1 % ; CFL* = c dt / h_min (c = 1 acoustique, sqrt(3) elastique, h diametre).
-- Maillage : 16 x 8 carres de (-1,1) x (0,1) (niveau 2 du code de l'article, qui reproduit la figure) remplis d'un seul
-- milieu. Lance par run_campaign.py (milieu, space.k, poids).

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 1.0, weight_elastic = 1.0 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 300 }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 2 }
problem = { medium = "acoustic", solution = "coupled_sinusoidal", omega = 1, theta = 10 }
study   = { type = "cfl", cfl_threshold = 0.05, cfl_decrement = 0.01 }
