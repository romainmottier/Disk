-- [M2AN25] Mottier, Ern, Khot, Guillot, "Hybrid high-order methods for elasto-acoustic wave propagation in the time
-- domain", ESAIM M2AN 59 (2025) 2685-2715, section 6.1, Fig. 3 : rayon spectral (racine de la plus grande valeur propre
-- de K_schur^T M_T^-1 K_schur X = gamma M_T X, K_schur complement de Schur du bloc faces) normalise par sqrt(#cellules),
-- en fonction du poids eta de la stabilisation O(1), ordres egal et mixte, k = 1..3 : acoustique pur, elastique pur,
-- couplage (eta^f = eta^s = eta). Lance par run_campaign.py, qui ajoute space.k, space.order, les poids et le milieu.
--
-- Maillage "l1" des fichiers de l'article (Coupling_rho_mo_eo_l1_k.png) : niveau 1 du code de l'article, 8 x 4 carres
-- de (-1,1) x (0,1) (32 cellules), le meme domaine pour les cas purs. Milieux de l'article : rho = 1, c_p^f = 1,
-- c_p^s = sqrt(3), c_s = 1.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 1.0, weight_elastic = 1.0 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 10 }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1 }
problem = { medium = "coupled", solution = "coupled_poly_in_time" }
study   = { type = "spectrum" }
