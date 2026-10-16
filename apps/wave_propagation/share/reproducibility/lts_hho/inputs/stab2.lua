-- stability, second-order formulation (Steins setting: mixed order, Lehrenfeld-Schoberl, 1/h, gamma = 3 gamma*)
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3 }
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", final_time = 1.0, n_steps = 200 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 18, ny = 18, level = 0, refine_levels = 0,
            refine_box = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat" }
