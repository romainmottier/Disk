-- First-order formulation, two-level LTS-RK4 (explicit setting) on meshes refined around (+-0.5, 0.5):
--  round(h_max/h_min) = 2 fine substeps. From level 2 on, the fine region is a small part of the mesh.
--  Run: ./wave_propagation first_order_lts.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "first_order"
space   = { k = 2, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "lts", final_time = 1.0, n_steps = 160 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 1 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 4 }
