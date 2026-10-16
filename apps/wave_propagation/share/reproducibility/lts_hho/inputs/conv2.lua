-- smooth convergence, second-order formulation (Steins setting, k = 2 mixed: cells of degree 3), haloscale meshes
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3 }
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "mlts", final_time = 1.0, n_steps = 57 }
mesh    = { type = "polymesher", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N0.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N1.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N2.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N3.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N4.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N5.txt" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", dirichlet = "zero" }
study   = { type = "space", last_level = 4 }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat" }
