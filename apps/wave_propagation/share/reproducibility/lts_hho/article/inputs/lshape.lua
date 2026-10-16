-- 6.6 L-shaped domain (-1/2,1/2)^2 \ [0,1/2)x(0,1/2], corner at the origin, u = sin(2 t) r^(4/3) sin(4 theta / 3),
-- graded quadtree meshes mu = 1/5, n = 4, 6, 8, 12, 16 (mesh.level 1..5), k = 3, t_f = 1/2, dt = 0.05 / n
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "mlts", final_time = 0.5, n_steps = 40 }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n4.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n6.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n8.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n12.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n16.txt" } }
problem = { medium = "acoustic", solution = "lshape_corner_sin", lshape_corner = {0.0, 0.0}, lshape_exponent = 4.0 / 3.0 }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
