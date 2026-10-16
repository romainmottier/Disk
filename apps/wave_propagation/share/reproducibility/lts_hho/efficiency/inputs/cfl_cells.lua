-- main.tex 6.4 (what sets the stable step of a cell): largest stable step of ERK4 (random initial data,
-- energy growth below 5 %, steps decreased by 1 %) on the meshes of generate_cfl_cells.py (mesh.files), acoustic,
-- first-order formulation, k = 3 (space.k = 1 appended for the thin / face / none meshes)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 1.0, n_steps = 160 }
mesh    = { type = "polygons", level = 1, files = { "MESH" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat" }
