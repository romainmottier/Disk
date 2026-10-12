-- Burman, Duran, Ern, Steins, "Convergence analysis of hybrid high-order methods for the wave equation",
-- J. Sci. Comput. 87 (2021) 91, hal-02922720, Section 6.1: Omega = (0,1)^2, T = 1, rho = 1, vS = 1, vP = sqrt(3),
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l, k = 1..3.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 1 (bottom right): HHO-Newmark (beta, gamma) = (0.25, 0.5), (0.275, 0.5), (0.275, 0.55) (time.newmark_beta,
-- time.newmark_gamma, set by run_campaign.py), equal order, k = 2, solution (80), h = 2^-6, dt = 0.1 x 2^-7: relative
-- variation of the discrete energy (69) (delta = 0). The published energy of (0.25, 0.5) is conserved to round-off
-- (1e-15), which time-dependent Dirichlet data exclude: the run of the paper kept the data of t = 0, zero for (80)
-- (data function of the original code evaluated at t = 0 and not updated), the initial data being those of (80):
-- problem.dirichlet = "zero" gives the published behaviour (conservation to 1e-15, plateau for (0.275, 0.5), growing
-- loss for (0.275, 0.55)) with 2 to 6 times smaller amplitudes. With the exact data of (80) (variant cases
-- jsc21_fig1_energy_exact_*), the energy (69) drifts by O(dt^2) (5e-6 for (0.25, 0.5)).

formulation = "second_order"
space   = { k = 2, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", newmark_beta = 0.25, newmark_gamma = 0.5, final_time = 1.0, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 6 }
problem = { medium = "elastic", solution = "elastic_non_polynomial", dirichlet = "zero" }
study   = { type = "run" }
output  = { energy = true }
