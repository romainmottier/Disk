-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 1 (left): HHO-Newmark (beta = 1/4, gamma = 1/2), equal order, solution (60) quadratic in time: errors vs h,
-- expected O(h^{k+2}) in L2 (cell unknowns) and O(h^{k+1}) in H1 (reconstructed gradient). dt: 0.1 x 2^-7.

formulation = "second_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", newmark_beta = 0.25, newmark_gamma = 0.5, final_time = 1.0, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_time" }
study   = { type = "space", last_level = 6, refine_dt = false }
