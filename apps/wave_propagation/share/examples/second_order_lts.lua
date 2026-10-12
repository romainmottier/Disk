-- Second-order formulation, two-level LTS-RK4 in the setting of Steins et al. on meshes refined twice around
--  (+-0.5, 0.5) (4 fine substeps): global ERK4 is unstable at this time step.
--  Run: ./wave_propagation second_order_lts.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3 }
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "lts", final_time = 1.0, n_steps = 320 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 2 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 3 }
