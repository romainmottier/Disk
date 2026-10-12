-- [CMAME25] Figs. 7 a 9 : efficacite (erreur dG en fonction du temps de calcul) des schemas explicites et implicites,
-- solution sinusoidale (5.2) : erreur spatiale dominante (omega = 5, theta = sqrt(2), Figs. 7 et 8) ou temporelle
-- (omega = 1, theta = 10, Fig. 9). Etude en espace sur les niveaux l = 1..5 (h = 2^-l : (4 x 2) 2^(l-1) carres de
-- (-1,1) x (0,1), ou triangles / polygones regeneres), temps de calcul (assemblage et integration) de chaque niveau sur
-- un fil (output.cpu). Lance par run_campaign.py, qui fixe pour chaque serie :
--   ERK(s) (ordre egal, O(1), poids 0.8 / 1.5) : dt = 0.9 CFL*(s,k) h_min / sqrt(3) (time.dt_h ; CFL* des Tables 2 et 3,
--     juste sous la limite de stabilite : la valeur exacte du tableau rend le niveau le plus fin instable) ;
--   SDIRK(s,s+1) (ordre mixte, O(1/h), solveur iteratif BiCGSTAB + ILU) : dt_l = 2^-(l-1)(k+1)/(s+1) dt_0 et
--     Tol_l = 2^-(l-1)(k+1) Tol_0 (time.steps_growth, solver.tolerance_growth), dt_0 = 1/16 (s = 1) ou 1/8 (s = 2, 3) :
--     plus grand pas pour lequel l'erreur stagne a 1 % pres au niveau l = 1 (procedure de l'article), Tol_0 = 1e-2.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "least_squares", scaling = "1", weight_acoustic = 0.8, weight_elastic = 1.5 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 10 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 0 }
problem = { medium = "coupled", solution = "coupled_sinusoidal", omega = 5, theta = math.sqrt(2) }
study   = { type = "space", last_level = 4 }
output  = { errors = "energy", cpu = true }
