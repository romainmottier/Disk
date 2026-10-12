-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 4 (left): HHO-SDIRK(3,4), first-order formulation, equal order (HHO stabilization), tau = O(1/h), solution (62):
-- errors vs h on the meshes l = 1..5 of the figure, dt = 0.1 x 2^-8, expected O(h^{k+2}) in L2 and O(h^{k+1}) in H1.
-- Stabilization of the code of the paper: weight 1, 1/h_F (space.stabilization_size = "face").

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 2560 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial" }
study   = { type = "space", last_level = 5, refine_dt = false }
