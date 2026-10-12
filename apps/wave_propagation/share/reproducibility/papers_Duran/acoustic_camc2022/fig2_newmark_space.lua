-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (left): HHO-Newmark, equal order, solution (62) non-polynomial (f = 0): errors vs h, expected O(h^{k+2}) in
-- L2 and O(h^{k+1}) in H1, the L2 error of k = 2 being saturated by the time error on the finest meshes. The published
-- errors are those of dt = 0.1 x 2^-6, not of dt = 0.1 x 2^-7 as in the caption: ratios 0.97-1.10 to all the values
-- read on the figure, whose time floor for k = 2 (5.4e-7) is 4 times the one of dt = 0.1 x 2^-7. The cases
-- camc22_fig2_space_dt7_k<k> (variant) use the dt of the caption.

formulation = "second_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", final_time = 1.0, n_steps = 640 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial" }
study   = { type = "space", last_level = 6, refine_dt = false }
