-- 5.2.1 smooth solution, second-order formulation (setting (2E), k = 2): Cartesian meshes 8 x 8 .. 128 x 128 of the unit square
-- (mesh.level 0..4), uniform or with the central box [3/8, 5/8]^2 refined mesh.refine_levels times (graded transition)
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3,
            stabilization_size = "cell" }   -- h = cell diameter h_T in the 1/h weight (paper, eq. stab-param-F)
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", final_time = 1.0, n_steps = 160 }
mesh    = { x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 8, ny = 8, level = 0, refine_levels = 0,
            refine_box = {0.375, 0.625, 0.375, 0.625} }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "run" }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
