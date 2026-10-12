-- Burman, Duran, Ern, Steins, "Convergence analysis of hybrid high-order methods for the wave equation",
-- J. Sci. Comput. 87 (2021) 91, hal-02922720, Section 6.1: Omega = (0,1)^2, T = 1, rho = 1, vS = 1, vP = sqrt(3),
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l, k = 1..3.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (top left), as described in the paper: HHO-SDIRK(2,3), first-order formulation (velocity-stress), equal
-- order (HHO stabilization, 1/h), solution (78), k = 1..3: errors vs h, expected O(h^{k+1}) for the H1-type error
-- (Theorem 4.3; stress) and one order more, O(h^{k+2}), in L2 (velocity). dt not given in the paper: 0.1 x 2^-7.
-- Stabilization of the code of the paper: tau = rho vS / h (weight 1). The published image is that of [CAMC22]
-- Fig. 3 (right), acoustic, mixed order, k = 0, 1, 2 (identical pixel for pixel): fig2_publiee_sdirk23_space.lua.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 2, final_time = 1.0, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "elastic", solution = "elastic_quadratic_in_time" }
study   = { type = "space", last_level = 6, refine_dt = false }
