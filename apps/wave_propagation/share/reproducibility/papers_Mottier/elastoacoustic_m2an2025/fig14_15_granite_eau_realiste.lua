-- [M2AN25] Figs. 14 (droite) et 15 : granite-eau realiste (6.6 : rho^f = 1025, c^f = 1500, rho^s = 2690, c_p^s = 6000,
-- c_s = 3000), domaine (-2500, 2500) x (-2500, 1000) m, fluide y > 0, ondelette (flux, theta = 10, f_c = 10, lambda =
-- 150 m) centree en (0, 100), SDIRK(3,4), ordre mixte, O(1/h), k = 3, n = 8 (10 x 2^8 pas pour T = 0.425 s). Energie
-- (Fig. 14) et instantanes a t = 0.1275 et 0.3825 s (Fig. 15). Le texte donne y_c = 500 m, mais les figures publiees
-- sont celles d'une source a environ 100 m de l'interface : centre de l'anneau de pression a 90-110 m sur les deux
-- instantanes publies (rayon exterieur c t + 130 m), energie elastique non nulle des t ~ 0.03 et maximale a t ~ 0.16 ;
-- avec y_c = 500 m, rien n'atteint l'interface avant t ~ 0.28 (variante). L'article utilise 560 x 560 cellules (h_x = 8.93 m,
-- h_y = 6.25 m : 140 x 140 raffine deux fois) ; avec l'assembleur actuel ce maillage demande environ 85 Go de memoire
-- (matrice globale, condensation statique, factorisation), plus que les 86 Go de la station de la campagne : la campagne
-- utilise 280 x 280 cellules (raffine une fois, h_x = 17.9 m, h_y = 12.5 m, 8 cellules par longueur d'onde, environ
-- 20 Go) ; mesh.level = 2 redonne le maillage de l'article sur une machine plus grosse.

formulation = "first_order"
space   = { k = 3, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.80, weight_elastic = 1.38 }
time    = { scheme = "sdirk", stages = 3, final_time = 0.425, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-2500.0, 2500.0}, y = {-2500.0, 1000.0}, nx = 140, ny = 140, level = 1 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 1025.0, vp = 1500.0 }, elastic = { rho = 2690.0, vp = 6000.0, vs = 3000.0 },
            ricker = { field = "acoustic_flux", center = {0.0, 100.0}, frequency = 10, amplitude = 10, velocity = 1500.0 } }
study   = { type = "run" }
output  = { energy = true, cpu = true, snapshots = { 0.1275, 0.3825 } }
