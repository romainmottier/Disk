-- [M2AN25] Figs. 11 et 12 : ondelette de Ricker, contrastes granite-eau (6.5a : rho^f = 0.5, c^f = 0.5, rho^s = 1.3,
-- c_p^s = 2, c_s = 1) et granite-air (6.5b : rho^s = 2200 rho^f = 800, c_p^s = 17.5 c^f = 6.36, c_s = 9 c^f = 3.27),
-- lambda = c^f / f_c ; l = 5 et 6. Capteurs S^f = (-0.1, 0.05), S^s = (-0.1, -0.05) : ceux du code de l'article, qui
-- reproduisent les courbes publiees (le texte donne (-0.05, +-0.1), coordonnees inversees ; la reference Gar6more2D en
-- ces points coincide avec les courbes noires publiees). Le reste comme les Figs. 9-10. Lance par run_campaign.py
-- (materiaux, space.k, mesh.level).

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.80, weight_elastic = 1.38 }
time    = { scheme = "sdirk", stages = 3, final_time = 0.5, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-0.5, 0.5}, y = {-0.5, 0.5}, nx = 2, ny = 2, level = 5 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 0.5, vp = 0.5 }, elastic = { rho = 1.3, vp = 2.0, vs = 1.0 },
            ricker = { field = "acoustic_flux", center = {0.0, 0.125}, frequency = 10, amplitude = 10, velocity = 0.5 } }
study   = { type = "run" }
output  = { cpu = true,
            sensors = { { name = "Sf", x = -0.1, y = 0.05, fields = {"acoustic_pressure", "acoustic_flux"} },
                        { name = "Ss", x = -0.1, y = -0.05, fields = {"elastic_velocity"} } } }
