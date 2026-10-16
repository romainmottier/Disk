-- efficiency study, cost: quadtree mesh (mesh.files), acoustic, k = 3, no source, homogeneous Dirichlet data (F = 0:
-- F(t) never assembled), fixed number of steps, work report and wall times of the integrator (output.profile)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 100 }
mesh    = { type = "polygons", level = 1, files = { "MESH" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", profile = true, cpu = true }
