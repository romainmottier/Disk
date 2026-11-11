-- 6.6 check of the cost model (Fig. cost_examples): the quadtree meshes of make_cost_examples.py (mesh.level 1: layer
-- with an inclusion, h_max / h_min = 8; mesh.level 2: small inclusion, h_max / h_min = 256), setting (1E), k = 3,
-- ERK(4), MLTS, LI; one thread, two lengths per scheme (the setup cancels in the difference)
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 0.1, n_steps = 40 }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/meshes_cost/cost_a.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lts_efficiency/meshes_cost/cost_b.txt" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
