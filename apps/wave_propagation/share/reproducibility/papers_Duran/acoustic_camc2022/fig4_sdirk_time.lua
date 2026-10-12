-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 4 (right): HHO-SDIRK(s,s+1), s = 1..3 (time.stages), equal order, tau = O(1/h), h = 2^-3, k = 3, solution (61):
-- errors vs dt = 0.1 x 2^-i, i = 0..7, expected O(dt^{s+1}). For k = 3 the solution (61), of degree 4 = k+1, is
-- reproduced by the reconstruction: the time error is the only error component. Stabilization of the code of the
-- paper: weight 1, 1/h_F.

formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 10 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 3 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_space" }
study   = { type = "time", time_levels = 8, reference = "exact" }
