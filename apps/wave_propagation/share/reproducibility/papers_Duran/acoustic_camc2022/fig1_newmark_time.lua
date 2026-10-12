-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 1 (right): HHO-Newmark, equal order, k = 2, h = 2^-6 (mesh not given in the paper): errors vs dt = 0.1 x 2^-i,
-- i = 0..7, expected O(dt^2). The published errors are those of p = cos(sqrt2 pi t) x(1-x) y(1-y), not of (61) as
-- printed (sin): ratios 1.00-1.02 to the 16 values read on the figure, against |tan(sqrt2 pi)| = 3.6 times smaller
-- errors with sin (Newmark phase error w^3 T dt^2/12 times |cos wT| for sin, |sin wT| for cos). The case
-- camc22_fig1_time_sin (variant) runs (61) as printed.

formulation = "second_order"
space   = { k = 2, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1 }
time    = { scheme = "newmark", final_time = 1.0, n_steps = 10 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 6 }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_space_cos" }
study   = { type = "time", time_levels = 8, reference = "exact" }
