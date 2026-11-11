-- 6.6 check of the cost model (Fig. cost_examples): the quadtree meshes of make_cost_examples.py (mesh.level 1: layer
-- with an inclusion, h_max / h_min = 8; mesh.level 2: small inclusion, h_max / h_min = 256), setting (2E), k = 2 mixed,
-- ERK(4), MLTS, LI; one thread, two lengths per scheme (the setup cancels in the difference)
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3,
            stabilization_size = "cell" }   -- h = cell diameter h_T in the 1/h weight (paper, eq. stab-param-F)
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", final_time = 0.1, n_steps = 40 }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/meshes_cost/cost_a.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/meshes_cost/cost_b.txt" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "run" }
solver  = { fast = true, matrix_free = true, matrix_free_tolerance = 0 }
output  = { file = "c.txt", table = "c.dat" }
