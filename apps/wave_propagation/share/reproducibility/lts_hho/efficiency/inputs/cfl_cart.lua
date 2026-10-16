-- efficiency study, stability: unit square 16 x 16, central box [3/8,5/8]^2 refined mesh.refine_levels times, graded
-- transition (meshes of 6.1), acoustic, k = 3, random initial data, largest stable number of steps (5 % growth, 1 %)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 160 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 16, ny = 16, level = 0, refine_levels = 0,
            refine_box = {0.375, 0.625, 0.375, 0.625} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat" }
