-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 6: HHO-ERK(4), first-order formulation, mixed order (Lehrenfeld-Schoberl), solution (62), dt = h/320 (h = 2^-l):
-- errors vs h on the meshes l = 1..5 of the figure. tau = O(1/h) ("1/h", left): O(h^{k+2}) in L2, O(h^{k+1}) in H1;
-- tau = O(1) ("1", right): O(h^{k+1}). Stabilization of the code of the paper: weight 1, 1/h_F.

formulation = "first_order"
space   = { k = 1, order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "erk", stages = 4, final_time = 1.0, n_steps = 640 }   -- dt = 1/640 = h/320 at h = 1/2
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial" }
study   = { type = "space", last_level = 5, refine_dt = true }
