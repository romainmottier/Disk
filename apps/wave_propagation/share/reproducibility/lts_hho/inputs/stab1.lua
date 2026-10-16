-- stability of the fully discrete schemes (preprint, Section 6.1.2): unit square, 18 x 18 cells, the central square
-- [1/3, 2/3]^2 refined mesh.refine_levels times (graded transition), acoustic, first-order formulation, k = 3
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 200 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 18, ny = 18, level = 0, refine_levels = 0,
            refine_box = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat" }
