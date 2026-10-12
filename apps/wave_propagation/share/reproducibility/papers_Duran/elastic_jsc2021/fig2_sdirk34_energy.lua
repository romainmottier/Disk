-- Burman, Duran, Ern, Steins, "Convergence analysis of hybrid high-order methods for the wave equation",
-- J. Sci. Comput. 87 (2021) 91, hal-02922720, Section 6.1: Omega = (0,1)^2, T = 1, rho = 1, vS = 1, vP = sqrt(3),
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l, k = 1..3.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (bottom), as described in the paper: HHO-SDIRK(3,4), first-order formulation, solution (80), h = 2^-5,
-- dt = 0.1 x 2^-7, k = 1..3, equal order (left, HHO stabilization) and mixed order (right, Lehrenfeld-Schoberl,
-- set by run_campaign.py): relative energy loss (E^0 - E^n)/E^0 of E = 1/2 (v, v)_rho + 1/2 (s, s)_{A^-1}
-- (output.energy). Weight 1. The published curves are the tau = O(1) curves of [CAMC22] Fig. 5 (acoustic, solution
-- (62), k = 0, 1, 2 labelled k = 1, 2, 3): fig2_publiee_sdirk34_energy.lua.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 5 }
problem = { medium = "elastic", solution = "elastic_non_polynomial" }
study   = { type = "run" }
output  = { energy = true }
