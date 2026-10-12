-- Second-order formulation, implicit SDIRK(3,4): mixed order, Lehrenfeld-Schoberl, 1/h, unit weight.
--  Space-time convergence (n_steps keeps the time error below the space error, see first_order_sdirk.lua).
--  Run: ./wave_propagation second_order_sdirk.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "second_order"
space   = { k = 2, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = 1 }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 80 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 0,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 0 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 4 }
