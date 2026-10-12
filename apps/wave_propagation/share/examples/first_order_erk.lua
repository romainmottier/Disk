-- First-order formulation, explicit ERK4: equal order, plain least squares, O(1) (MEG/MEKG papers).
--  Space-time convergence on uniform meshes.
--  Run: ./wave_propagation first_order_erk.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "first_order"
space   = { k = 2, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 80 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 0,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 0 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 4 }
