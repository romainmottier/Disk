-- Burman, Duran, Ern, Steins, J. Sci. Comput. 87 (2021) 91, Fig. 2 as published: its image is that of the acoustic
-- results of Burman, Duran, Ern, Commun. Appl. Math. Comput. 4 (2022) 597-633, Section 5.1 (Omega = (0,1)^2, T = 1,
-- rho = kappa = 1, quadrangular meshes h = 2^-l), whose configuration this file runs.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (top right) as published = [CAMC22] Fig. 4 (right): HHO-SDIRK(s,s+1), s = 1..3 (time.stages), equal order,
-- tau = O(1/h), h = 2^-3, k = 3, acoustic solution p = sin(sqrt2 pi t) x(1-x) y(1-y) ([CAMC22] (61)): errors vs
-- dt = 0.1 x 2^-i, i = 0..7, weight 1, 1/h_F.

formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 10 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 3 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_space" }
study   = { type = "time", time_levels = 8, reference = "exact" }
