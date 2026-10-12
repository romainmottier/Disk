-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 3: HHO-SDIRK(2,3), first-order formulation, mixed order (Lehrenfeld-Schoberl), solution (60): errors vs h.
-- tau = O(1) (space.scaling = "1", left): O(h^{k+1}) in L2 and H1; tau = O(1/h) ("1/h", right): O(h^{k+2}) in L2,
-- O(h^{k+1}) in H1. L2: scalar cell unknown (time derivative of the pressure), H1: flux (gradient of the pressure).
-- Stabilization of the code of the paper: tau = 1/(rho c) (weight 1), times 1/h_F for tau = O(1/h) (face diameter,
-- space.stabilization_size = "face": ratios 1.0-1.1 to the published values). dt not given in the paper: 0.1 x 2^-6.

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 2, final_time = 1.0, n_steps = 640 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_time" }
study   = { type = "space", last_level = 6, refine_dt = false }
