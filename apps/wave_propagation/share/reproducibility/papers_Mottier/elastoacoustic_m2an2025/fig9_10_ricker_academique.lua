-- [M2AN25] Figs. 9 et 10 : ondelette de Ricker en condition initiale (flux acoustique m0 = theta exp(-pi^2 r^2 /
-- lambda^2) (x - x_c, y - y_c), theta = 10, lambda = c_p^f / f_c, f_c = 10), milieux academiques (6.4), domaine
-- (-0.5,0.5)^2, fluide y > 0, solide y < 0, centre (0, 0.125), conditions de Dirichlet homogenes, SDIRK(3,4), ordre
-- mixte, stabilisation O(1/h), n = 9 (Delta t = 0.1 x 2^-9). Capteurs S^f = (-0.15, 0.1), S^s = (-0.15, -0.1).
-- Niveau l du code de l'article : 2 x 2^l carres par cote (l = 4 : h = 1/32). Lance par run_campaign.py (space.k,
-- mesh.level, duree, instantanes).

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.80, weight_elastic = 1.38 }
time    = { scheme = "sdirk", stages = 3, final_time = 0.25, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-0.5, 0.5}, y = {-0.5, 0.5}, nx = 2, ny = 2, level = 4 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 1.0, vp = 1.0 }, elastic = { rho = 1.0, vp = math.sqrt(3.0), vs = 1.0 },
            ricker = { field = "acoustic_flux", center = {0.0, 0.125}, frequency = 10, amplitude = 10, velocity = 1.0 } }
study   = { type = "run" }
output  = { cpu = true,
            sensors = { { name = "Sf", x = -0.15, y = 0.1, fields = {"acoustic_pressure", "acoustic_flux"} },
                        { name = "Ss", x = -0.15, y = -0.1, fields = {"elastic_velocity"} } } }
