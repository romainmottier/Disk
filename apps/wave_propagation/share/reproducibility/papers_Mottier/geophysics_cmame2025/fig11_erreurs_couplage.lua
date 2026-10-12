-- [CMAME25] Fig. 11 : erreurs ponctuelles sur les conditions de couplage au capteur d'interface S^i, ondelette de
-- Ricker en condition initiale (flux m0 = theta exp(-pi^2 r^2 / lambda^2) (x - x_c, y - y_c), theta = 1, f_c = 10,
-- lambda = c_p^f / f_c), milieux academiques (5.6), domaine (-0.5,0.5)^2, fluide y > 0, centre (0, 0.125), T = 1,
-- SDIRK(3,4), Delta t = 0.1 x 2^-8, ordre mixte, O(1/h), poids 0.8 / 1.5, k = 1 et 3, l = 4..7 (h = 2^-l :
-- 2^(l-1) x 2^(l-1) fois 2 x 2 carres). Capteur S^i = (-0.15, 0) (position du croquis de la Fig. 10, qui reproduit les
-- courbes publiees ; le texte donne (-0.3, 0)) : p_F, v_F (inconnues de face), m_T, s_T. Lance par run_campaign.py
-- (space.k, mesh.level) ; Delta_kin = |(v_F - m_T).n|, Delta_dyn = |p_F n - s_T n| calcules par figures.py.

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.8, weight_elastic = 1.5 }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-0.5, 0.5}, y = {-0.5, 0.5}, nx = 2, ny = 2, level = 3 }
problem = { medium = "coupled", interface = "y", solution = "ricker",
            acoustic = { rho = 1.0, vp = 1.0 }, elastic = { rho = 1.0, vp = math.sqrt(3.0), vs = 1.0 },
            ricker = { field = "acoustic_flux", center = {0.0, 0.125}, frequency = 10, amplitude = 1, velocity = 1.0 } }
study   = { type = "run" }
output  = { sensors = { { name = "Si", x = -0.15, y = 0.0,
                          fields = {"interface_velocity", "interface_pressure", "acoustic_flux", "elastic_stress"} } } }
