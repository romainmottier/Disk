-- 5.2.2 L-shaped domain, second-order formulation (setting (2E), k = 2): meshes of lshape.lua (mu = 1/5, n = 4, 6, 8, 12, 16)
-- graded quadtree meshes mu = 1/5, n = 4, 6, 8, 12, 16 (mesh.level 1..5), k = 2 mixed order, t_f = 1/2, dt = 1/(25 n) (80 % of the CFL limit of the coarse cells in the setting (2E))
formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3,
            stabilization_size = "cell" }   -- h = cell diameter h_T in the 1/h weight (paper, eq. stab-param-F)
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "mlts", final_time = 0.5, n_steps = 40 }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n4.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n6.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n8.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n12.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/lshape_mu1-5_n16.txt" } }
problem = { medium = "acoustic", solution = "lshape_corner_sin", lshape_corner = {0.0, 0.0}, lshape_exponent = 4.0 / 3.0 }
study   = { type = "run" }
solver  = { fast = true }
output  = { file = "c.txt", table = "c.dat", cpu = true }
