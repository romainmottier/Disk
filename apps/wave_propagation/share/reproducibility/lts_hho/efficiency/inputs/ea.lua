-- LTS = MLTS with two levels in a coupled medium: elastic (x < 0) / acoustic (x > 0), manufactured solution with
-- nonzero interface traces, the elastic cells along the interface refined once
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "lts", final_time = 0.5, n_steps = 100 }
mesh    = { x = {-1.0, 1.0}, y = {0.0, 1.0}, nx = 8, ny = 4, level = 0, refine_levels = 1,
            refine_box = {-0.25, 0.0, 0.25, 0.75} }
problem = { solution = "interface_coupling" }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", profile = true }
