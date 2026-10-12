-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (right): HHO-Newmark (beta, gamma) = (0.25, 0.5), (0.275, 0.5), (0.275, 0.55) (time.newmark_beta,
-- time.newmark_gamma, set by run_campaign.py), equal order, k = 2, solution (62) (f = 0, homogeneous Dirichlet data),
-- h = 2^-6: relative variation of the discrete energy 1/2 v^T M v + 1/2 u^T K u, exactly conserved for (0.25, 0.5)
-- (Lemma 3.3). The published curves are those of dt = 0.1 x 2^-5, not of dt = 0.1 x 2^-7 as in the caption: ratios
-- 0.9-1.1, against 16 (0.275, 0.5: (beta - 1/4) dt^2 w^2 sin^2(wt)) and 4 (0.275, 0.55: (gamma - 1/2) w^2 dt t) times
-- smaller variations with dt = 0.1 x 2^-7 (variant cases camc22_fig2_energy_dt7_*).

formulation = "second_order"
space   = { k = 2, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", newmark_beta = 0.25, newmark_gamma = 0.5, final_time = 1.0, n_steps = 320 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 6 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial" }
study   = { type = "run" }
output  = { energy = true }
