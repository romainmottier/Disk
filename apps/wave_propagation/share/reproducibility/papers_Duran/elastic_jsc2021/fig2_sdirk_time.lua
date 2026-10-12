-- Burman, Duran, Ern, Steins, "Convergence analysis of hybrid high-order methods for the wave equation",
-- J. Sci. Comput. 87 (2021) 91, hal-02922720, Section 6.1: Omega = (0,1)^2, T = 1, rho = 1, vS = 1, vP = sqrt(3),
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l, k = 1..3.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (top right), as described in the paper: HHO-SDIRK(s,s+1), s = 1..3 (time.stages), first-order formulation,
-- equal order, k = 2, solution (79) quadratic in space: errors vs dt = 0.1 x 2^-i, i = 0..7, expected O(dt^{s+1}).
-- Mesh not given in the paper: h = 2^-6, so that the error in space (the solution is of degree 4 > k + 1) stays below
-- the time error of SDIRK(3,4) down to dt = 0.1 x 2^-4. Weight 1. The published curves are those of [CAMC22] Fig. 4
-- (right), acoustic, k = 3, h = 2^-3 (identical pixel for pixel outside the legend): fig2_publiee_sdirk_time.lua.

formulation = "first_order"
space   = { k = 2, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 10 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 6 }
problem = { medium = "elastic", solution = "elastic_quadratic_in_space" }
study   = { type = "time", time_levels = 8, reference = "exact" }
