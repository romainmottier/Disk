-- main.tex, what sets the CFL limit of a cell, setting (2E): largest stable step of ERK4 (random initial data, energy
-- growth below 5 %, steps decreased by 1 %) on the meshes of generate_cfl_cells.py (mesh.files), acoustic, second-order
-- formulation, k = 2 mixed order, O(1/h) stabilization with the weight 3 gamma* (gamma* computed on the mesh), as stab2.lua
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3,
            stabilization_size = "cell" }   -- h = cell diameter h_T in the 1/h weight (paper, eq. stab-param-F)
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", final_time = 1.0, n_steps = 170 }
mesh    = { type = "polygons", level = 1, files = { "MESH" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero", source = "zero" }
study   = { type = "cfl", cfl_random = true, cfl_decrement = 0.01, cfl_threshold = 0.05 }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat" }
