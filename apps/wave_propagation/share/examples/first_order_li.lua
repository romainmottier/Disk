-- First-order formulation, locally implicit LI-SDIRK(4,3) (explicit setting, SDIRK(3,4) in the fine region)
--  on meshes refined around (+-0.5, 0.5).
--  Run: ./wave_propagation first_order_li.lua   (all parameters: ../first_order.lua, ../second_order.lua)

formulation = "first_order"
space   = { k = 2, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "li", final_time = 1.0, n_steps = 160 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 4, ny = 2, level = 1,
            refine_points = {{-0.5, 0.5}, {0.5, 0.5}}, refine_levels = 1 }
problem = { solution = "interface_coupling" }
study   = { type = "space", last_level = 4 }
