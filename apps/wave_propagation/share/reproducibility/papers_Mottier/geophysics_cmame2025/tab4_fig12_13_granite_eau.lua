-- [CMAME25] Table 4, Figs. 12 et 13 : granite-eau realiste (5.7 : rho^f = 1025, c^f = 1500, rho^s = 2690, c_p^s = 6000,
-- c_s = 3000), domaine (-1800, 2000) x (-1600, 1000) m (L = 3.8 km, H^f = 1 km, H^s = 1.6 km), pression initiale
-- p0 = exp(-pi^2 r^2 / lambda^2) (lambda = 150 m, f_c = 10) centree en (0, 200), T = 0.5 s, carres de 33.3 m (114 x 78 :
-- les pas adimensionnes de la Table 4, c_p^s Delta t / h = 0.398 pour Delta t = 0.1 x 2^-5, imposent h = 47.1 m, diametre
-- de ces carres ; le "h = 23.6 m" du texte en est la moitie). Capteurs S^f = (-500, 200), S^s = (-1000, -500),
-- S^i = (-400, 0). SDIRK(3,4) (ordre mixte, O(1/h), 160 pas, solveur direct ou iteratif de tolerance 1e-2) et ERK(4)
-- (ordre egal, O(1), 400 pas pour k = 1, 1180 pour k = 3). Lance par run_campaign.py (schema, k, pas, solveur).

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.8, weight_elastic = 1.5 }
time    = { scheme = "sdirk", stages = 3, final_time = 0.5, n_steps = 160 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-1800.0, 2000.0}, y = {-1600.0, 1000.0}, nx = 114, ny = 78, level = 0 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 1025.0, vp = 1500.0 }, elastic = { rho = 2690.0, vp = 6000.0, vs = 3000.0 },
            ricker = { field = "acoustic_pressure", center = {0.0, 200.0}, frequency = 10, amplitude = 1, velocity = 1500.0 } }
study   = { type = "run" }
output  = { cpu = true,
            sensors = { { name = "Sf", x = -500.0, y = 200.0, fields = {"acoustic_pressure", "acoustic_flux"} },
                        { name = "Ss", x = -1000.0, y = -500.0, fields = {"elastic_velocity"} },
                        { name = "Si", x = -400.0, y = 0.0, fields = {"interface_velocity", "interface_pressure", "elastic_velocity", "acoustic_flux"} } } }
