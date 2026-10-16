-- 6.3 stability of the fully discrete schemes, second-order formulation (setting (2E)): unit square, 16 x 16 cells, the
-- central box [3/8, 5/8]^2 refined mesh.refine_levels times with a graded transition (the meshes of stab.lua), acoustic,
-- k = 2 mixed order, weight 3 gamma*, random initial cell data, largest stable number of steps searched downwards by 1 %;
-- LTS and MLTS with the exact condensation or with the projected one (faces.projected)
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3 }
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", final_time = 1.0, n_steps = 170 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 16, ny = 16, level = 0, refine_levels = 0,
            refine_box = {0.375, 0.625, 0.375, 0.625} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat" }
