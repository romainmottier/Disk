-- Burman, Duran, Ern, Steins, J. Sci. Comput. 87 (2021) 91, Fig. 2 as published: its image is that of the acoustic
-- results of Burman, Duran, Ern, Commun. Appl. Math. Comput. 4 (2022) 597-633, Section 5.1 (Omega = (0,1)^2, T = 1,
-- rho = kappa = 1, quadrangular meshes h = 2^-l), whose configuration this file runs.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (top left) as published = [CAMC22] Fig. 3 (right): HHO-SDIRK(2,3), first-order formulation, mixed order
-- (Lehrenfeld-Schoberl), tau = O(1/h), acoustic solution p = t^2 sin(pi x) sin(pi y) ([CAMC22] (60)), k = 0, 1, 2:
-- errors vs h (L2: velocity, H1: flux), weight 1, 1/h_F, dt = 0.1 x 2^-6 (as in ../acoustic_camc2022).

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 2, final_time = 1.0, n_steps = 640 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_time" }
study   = { type = "space", last_level = 6, refine_dt = false }
