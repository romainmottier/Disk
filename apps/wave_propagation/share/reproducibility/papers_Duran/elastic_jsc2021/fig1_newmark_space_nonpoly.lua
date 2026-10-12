-- Burman, Duran, Ern, Steins, "Convergence analysis of hybrid high-order methods for the wave equation",
-- J. Sci. Comput. 87 (2021) 91, hal-02922720, Section 6.1: Omega = (0,1)^2, T = 1, rho = 1, vS = 1, vP = sqrt(3),
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l, k = 1..3.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 1 (bottom left): HHO-Newmark, equal order, solution (80) non-polynomial (f = 0), dt = 0.1 x 2^-10: errors vs h,
-- expected O(h^{k+2}) in L2 and O(h^{k+1}) in H1. Final time of the code of the paper for this solution: T = 0.125
-- (IHHOSecondOrder of elastodynamics.cpp, tf = 0.125), 1280 steps of 0.1 x 2^-10: ratios 0.98-1.06 to the published
-- values, against 1.6-1.9 times larger H1 errors at T = 1.

formulation = "second_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", final_time = 0.125, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 1 }
problem = { medium = "elastic", solution = "elastic_non_polynomial" }
study   = { type = "space", last_level = 6, refine_dt = false }
