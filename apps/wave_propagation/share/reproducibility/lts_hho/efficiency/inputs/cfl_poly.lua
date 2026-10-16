-- efficiency study, stability on a quadtree mesh (mesh.files, generate_refined_quadtree.py), acoustic, k = 3, random
-- initial data, largest stable number of steps (5 % growth, 1 %)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 160 }
mesh    = { type = "polygons", level = 1, files = { "MESH" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat" }
