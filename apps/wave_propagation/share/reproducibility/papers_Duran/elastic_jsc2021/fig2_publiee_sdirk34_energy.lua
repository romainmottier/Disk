-- Burman, Duran, Ern, Steins, J. Sci. Comput. 87 (2021) 91, Fig. 2 as published: its image is that of the acoustic
-- results of Burman, Duran, Ern, Commun. Appl. Math. Comput. 4 (2022) 597-633, Section 5.1 (Omega = (0,1)^2, T = 1,
-- rho = kappa = 1, quadrangular meshes h = 2^-l), whose configuration this file runs.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 2 (bottom) as published = the tau = O(1) curves of [CAMC22] Fig. 5 (k = 0, 1, 2, labelled k = 1, 2, 3 in
-- [JSC21]: ratios 0.95-1.05 between the values read on both figures): HHO-SDIRK(3,4), first-order formulation, acoustic
-- solution p = sin(sqrt2 pi t) sin(pi x) sin(pi y) / (sqrt2 pi) ([CAMC22] (62)), tau = O(1), h = 2^-6 (as the
-- published [CAMC22] Fig. 5), dt = 0.1 x 2^-7, equal order (left) and mixed order (right, set by run_campaign.py):
-- relative energy loss (E^0 - E^n)/E^0 of E = 1/2 (|v_T|^2 + |sigma_T|^2). Weight 1.

formulation = "first_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1", weight = 1, stabilization_size = "face" }
time    = { scheme = "sdirk", stages = 3, final_time = 1.0, n_steps = 1280 }
solver  = { static_condensation = true }
mesh    = { type = "cartesian", x = {0.0, 1.0}, y = {0.0, 1.0}, nx = 1, ny = 1, level = 6 }
problem = { medium = "acoustic", solution = "acoustic_non_polynomial" }
study   = { type = "run" }
output  = { energy = true }
