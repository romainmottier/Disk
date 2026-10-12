-- Burman, Duran, Ern, "Hybrid high-order methods for the acoustic wave equation in the time domain",
-- Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, Section 5.1: Omega = (0,1)^2, T = 1, rho = kappa = 1,
-- quadrangular meshes h = 2^-l, l = 1..6, dt = 0.1 x 2^-l.
-- Run by run_campaign.py (this folder), which sets space.k (and the varying parameters) at the end of a complete copy
-- per case, figure_<n>/<case>/<case>.lua in the build, which is what is run; standalone:
-- ./wave_propagation <this file> "space.k = 1"
--
-- Fig. 8: polygonal meshes, solution (60), same rates as on quadrangular meshes. The PolyMesher meshes of the paper
-- (meshes/unit_square_polymesh_nel_<n>.txt.bz2, n = 20 x 2^l cells, l = 0..8; Fig. 7 draws the first four):
-- run_campaign.py decompresses them in figure_8/meshes and sets mesh.files and the scheme of each panel. h on the
-- figure: smallest cell diameter (h_min, as the code of the paper). Newmark and SDIRK(3,4): dt = 0.1 x 2^-5 on the 9
-- meshes; ERK(4): dt = h_min/320 (time.dt_h) on the first 6, h_min/640 for tau = O(1/h) and k = 2 (stability limit in
-- h^2 exceeded by h_min/320 on meshes 5 and 6; same errors on meshes 1 to 4). Stabilization of the code of the paper:
-- weight 1, 1/h_F. Defaults below: HHO-Newmark, equal order.

formulation = "second_order"
space   = { k = 1, order = "equal", stabilization = "hho", scaling = "1/h", weight = 1, stabilization_size = "face" }
time    = { scheme = "newmark", final_time = 1.0, n_steps = 320 }
solver  = { static_condensation = true }
mesh    = { type = "polymesher", level = 1, files = {} }
problem = { medium = "acoustic", solution = "acoustic_quadratic_in_time" }
study   = { type = "space", last_level = 9, refine_dt = false }
