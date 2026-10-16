-- smooth convergence, first-order formulation, k = 3, haloscale meshes N = 0..5 (mesh.level = N + 1, h_max/h_min = 4^(N+1))
-- p = cos(sqrt2 pi t) sin(pi x) sin(pi y) type solution without source (acoustic_non_polynomial), t_f = 1, dt = 0.05 h_max
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "mlts", final_time = 1.0, n_steps = 57 }
mesh    = { type = "polymesher", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N0.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N1.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N2.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N3.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N4.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/centersquarehaloscale_graded_k3_N5.txt" } }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial", source = "zero", dirichlet = "zero" }
study   = { type = "space", last_level = 5 }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
