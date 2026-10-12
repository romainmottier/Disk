-- Second-order formulation, explicit ERK4 in the optimal setting of Steins et al. (M2AN 2023): mixed order,
--  Lehrenfeld-Schoberl, 1/h, gamma = 3 gamma*, face splitting.
--  Run: ./wave_propagation second_order_erk.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = "auto", weight_factor = 3 }
faces   = { solver = "splitting", tolerance = 1e-14 }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 120 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 0,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 0 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 4 }
