-- 6.2 smooth solution, k = 3: Cartesian meshes 8 x 8 .. 128 x 128 of the unit square (mesh.level 0..4), uniform or
-- with the central box [3/8, 5/8]^2 refined mesh.refine_levels times (graded transition, one level per ring)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 100 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 8, ny = 8, level = 0, refine_levels = 0,
            refine_box = {0.375, 0.625, 0.375, 0.625} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
