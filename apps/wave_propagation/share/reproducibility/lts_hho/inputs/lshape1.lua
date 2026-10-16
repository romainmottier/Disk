-- L-shaped domain (0,1)^2 \ [0.5,1)x(0.5,1], u = cos(2t) r^(2/3) sin(2 theta/3) (source and Dirichlet data), meshes graded
-- towards the corner h_corner ~ h^((k+1)/gamma) (k = 3: h_max/h_min = 1, 32, 1024, 32768 for N = 0..3), t_f = 0.5
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "mlts", final_time = 0.5, n_steps = 29 }
mesh    = { type = "polymesher", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_graded_k3_N0.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_graded_k3_N1.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_graded_k3_N2.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_graded_k3_N3.txt" } }
problem = { medium = "acoustic", solution = "lshape_corner" }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
