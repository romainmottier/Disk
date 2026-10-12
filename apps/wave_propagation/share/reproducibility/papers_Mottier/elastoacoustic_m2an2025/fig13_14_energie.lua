-- [M2AN25] Figs. 13 et 14 (gauche) : historique de l'energie discrete (totale, acoustique, elastique) de l'ondelette
-- de Ricker, SDIRK(3,4), ordre mixte, O(1/h), n = 9. Fig. 13 gauche : milieux academiques, k = 3, l = 6, t dans [0, 10] ;
-- droite : perte relative d'energie (E0 - E(t)) / E0 pour t dans [0, 1], k = 1..3, l = 4..6. Fig. 14 gauche :
-- granite-eau (6.5a), k = 3, l = 6. Lance par run_campaign.py (materiaux, space.k, mesh.level, duree).

formulation = "first_order"
space   = { k = 3, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.80, weight_elastic = 1.38 }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 5120 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-0.5, 0.5}, y = {-0.5, 0.5}, nx = 2, ny = 2, level = 6 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 1.0, vp = 1.0 }, elastic = { rho = 1.0, vp = math.sqrt(3.0), vs = 1.0 },
            ricker = { field = "acoustic_flux", center = {0.0, 0.125}, frequency = 10, amplitude = 10, velocity = 1.0 } }
study   = { type = "run" }
output  = { energy = true, cpu = true }
