-- 6.6 L-shaped domain (-1/2,1/2)^2 \ [0,1/2)x(0,1/2], corner at the origin, u = sin(2 t) r^nu sin(nu theta),
-- graded quadtree meshes mu = 1/6 (grading of the approximation theory for nu = 2/3), n = 4, 8, 16, 32 (mesh.level 1..4), k = 3, t_f = 1/2
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "mlts", final_time = 0.5, n_steps = 40 }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-6_n4.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-6_n8.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-6_n16.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-6_n32.txt" } }
problem = { medium = "acoustic", solution = "lshape_corner_sin", lshape_corner = {0.0, 0.0}, lshape_exponent = 2.0 / 3.0 }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
