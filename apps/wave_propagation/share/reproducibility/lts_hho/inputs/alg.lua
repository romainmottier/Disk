-- algebraic comparison of the two formulations: unit square, Cartesian 8 x 8 .. 64 x 64, acoustic, ERK(4)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 160 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 8, ny = 8, level = 0 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "space", last_level = 3 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
