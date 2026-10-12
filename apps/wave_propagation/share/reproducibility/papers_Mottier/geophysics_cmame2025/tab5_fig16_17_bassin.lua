-- [CMAME25] Table 5, Figs. 16 et 17 : bassin sedimentaire. Maillages gmsh de la these de l'auteur (VTK convertis par
-- ../meshes.py, etiquettes 1 = socle, 2 = sediments, 3 = atmosphere) : triangles, quadrangles, hybride, et non
-- conforme (quadrangles raffines une fois dans l'atmosphere et le bassin, noeuds pendants). Materiaux de la Fig. 14 :
-- atmosphere rho = 1.225, c = 343 ; sediments 1300, 1600, 900 ; socle 2570, 5350, 3009. Vitesse initiale elastique en
-- ondelette de Ricker (f_c = 4 Hz, lambda = 5350 / 4 m) centree en (563, 233) m, SDIRK(3,4), ordre mixte, O(1/h),
-- k = 3, Delta t = 1.5 / 320 = 4.69e-3 s (4.7e-3 dans l'article), T = 1.5 s. Capteurs S^f = (-650, 1300) (atmosphere)
-- et S^s = (-750, 1080) (bassin). Lance par run_campaign.py (maillage, k).

formulation = "first_order"
space   = { k = 3, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight_acoustic = 0.8, weight_elastic = 1.5 }
time    = { scheme = "sdirk", stages = 3, final_time = 1.5, n_steps = 320 }
solver  = { static_condensation = true }
mesh    = { type = "polygons", level = 1 }
problem = { medium = "coupled", interface = "tags", solution = "ricker",
            zones = { [1] = { medium = "elastic", rho = 2570.0, vp = 5350.0, vs = 3009.0 },
                      [2] = { medium = "elastic", rho = 1300.0, vp = 1600.0, vs = 900.0 },
                      [3] = { medium = "acoustic", rho = 1.225, vp = 343.0 } },
            ricker = { field = "elastic_velocity", center = {563.0, 233.0}, frequency = 4, amplitude = 1, velocity = 5350.0 } }
study   = { type = "run" }
output  = { cpu = true,
            sensors = { { name = "Sf", x = -650.0, y = 1300.0, fields = {"acoustic_pressure", "acoustic_flux"} },
                        { name = "Ss", x = -750.0, y = 1080.0, fields = {"elastic_velocity"} } } }
