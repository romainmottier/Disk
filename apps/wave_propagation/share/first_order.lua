-- =====================================================================================================================
--  wave_propagation input: elasto-acoustic coupling, FIRST-ORDER formulation
--  (velocity-stress in the solid, velocity-pressure in the fluid; Mottier-Ern-Guillot papers)
--
--      ./wave_propagation first_order.lua ["lua statements executed after this file" ...]
--      e.g. ./wave_propagation first_order.lua "time.scheme = 'sdirk'; space.order = 'mixed'; space.scaling = '1/h'"
--
--  Every entry is optional: a missing entry takes the default given in the comment; unknown entries are reported.
--  This file is plain Lua: variables, loops and expressions (2^l, math.sqrt, ...) can be used.
--  Annotated template, copied next to the executable by CMake: edit a copy with another name.
-- =====================================================================================================================

formulation = "first_order"         -- "first_order" | "second_order"

-- Space discretization --------------------------------------------------------------------------------------------------
--  Papers' settings: explicit integrators (erk, lts, li): order = "equal", stabilization = "least_squares", scaling = "1"
--                    implicit integrator (sdirk):        order = "mixed", stabilization = "least_squares", scaling = "1/h"
space = {
    k = 2,                          -- face polynomial degree (default 1)
    order = "equal",                -- "equal": cell degree k | "mixed": cell degree k+1 (default "equal")
    stabilization = "least_squares",-- "least_squares": plain least squares in equal order, Lehrenfeld-Schoberl in mixed
                                    --   order (block-diagonal face block, required by erk, lts, li)
                                    -- "hho": HHO stabilization with its high-order correction (sdirk only, always 1/h)
                                    -- (default "least_squares")
    scaling = "1",                  -- weight tau = eta (rho c)^(+-1) h^(-alpha): "1" (alpha = 0, keeps the CFL of the
                                    --   explicit schemes) | "1/h" (alpha = 1, implicit schemes) (default "1")
                                    --   | "graded" (O(1) times max(1, h_max / h_T)) | "hybrid" (1/h in the cells smaller
                                    --   than 0.75 h_max, the implicit region of li, O(1) elsewhere)
}

-- Time discretization ---------------------------------------------------------------------------------------------------
time = {
    scheme = "erk",                 -- "erk": ERK(s) | "sdirk": SDIRK(s, s+1) | "lts": two-level LTS-RK4 (substeps in the
                                    --   fine region) | "li": LI-SDIRK(4,3), ERK4 coarse / SDIRK(3,4) fine | "mlts":
                                    --   multilevel LTS-RK4, levels from the cell sizes (default "erk")
    stages = 4,                     -- erk: 1..4, sdirk: 1..3; 0: 4 for erk, 3 for sdirk (default 0)
    final_time = 1.0,               -- (default 1)
    n_steps = 80,                   -- time steps on the first mesh of the study, dt = final_time / n_steps;
                                    --   doubled with each level of a space study (default 10)
    substeps = 0,                   -- lts: fine substeps, 0: round(h_max / h_min) (default 0)
    fine_size = 0,                  -- lts, li: cells smaller than this (+ one layer of vertex neighbours) are fine;
                                    --   0: 0.75 h_max, everything on a uniform mesh (default 0); also the threshold
                                    --   of space.scaling = "hybrid"
    -- fine_layers = 1,             -- lts, li: layers of vertex neighbours added to the fine region (default 1)
    -- fine_measure = "diameter",   -- lts, li: size of a cell compared with fine_size, "diameter" | "width" (smallest
                                    --   altitude: the stable step of a cell follows its width) (default "diameter")
    -- level_ratio = 2,             -- mlts: ratio of the steps of two consecutive levels (default 2)
    -- level_ratios = {2, 2, 8},    -- mlts: ratios per level, coarsest first, the last one repeated (instead of
                                    --   level_ratio; e.g. a finest level covering three cell sizes)
    -- level_measure = "h",         -- mlts: levels from h_max / h_T ("h") or from max(h/c) / (h_T/c_T) ("h/c",
                                    --   heterogeneous media) (default "h")
    -- level_size = "diameter",     -- mlts: size of a cell in the level measure, "diameter" | "area" (4 |T| / |dT|) |
                                    --   "width" (smallest altitude: thin cells, e.g. pieces of cells cut along fitted
                                    --   interfaces) (default "diameter")
    -- level_tolerance = 0.0,       -- mlts: measure of a cell divided by 1 + level_tolerance before its level is chosen
                                    --   (media of nearly equal speeds on the same level) (default 0)
    -- li_coupling = "taylor",      -- li (matrix_free): coarse values in the implicit stages, "taylor" | "stage" (stage
                                    --   values of the SDIRK tableau: no order reduction) (default "taylor")
    -- li_condensation = false,     -- li (matrix_free): static condensation of the fine cell unknowns (default false)
    -- li_system = "dynamic",       -- li (matrix_free): "dynamic" | "constraint" (interior fine faces by their static
                                    --   condensation: well conditioned on very stiff fine regions) (default "dynamic")
    -- overlap = 1,                 -- mlts: rings of vertex neighbours annexed by each level (default 1)
    -- overlaps = {1, 0},           -- mlts: rings annexed by the levels 1, 2, ... (the last one repeated), instead of
                                    --   overlap
    -- mlts_reuse = false,          -- mlts (matrix_free): own coefficients of a level from the chain a_{i+1} = B a_i on
                                    --   its deep unknowns (about 4 instead of 7 products per unknown of a thick level)
    -- mlts_refresh = "window",     -- mlts: "window" (Algorithm 4) | "substep" (prototypes) (default "window")
}

solver = {
    static_condensation = false,    -- sdirk: static condensation of the cell unknowns in the stage systems (default false)
    matrix_free = false,            -- fast mode (PARAMETRES.txt, section 7), same results up to round-off: erk, lts, li,
                                    --   mlts, cfl study on the matrix-free explicit operator (shared local blocks, OpenMP); sdirk:
                                    --   parallel products, Pardiso without iterative refinement; parallel assembly, cached
                                    --   F(t) data (default false)
    matrix_free_tolerance = 1e-10,  -- local blocks shared to this relative tolerance; 0: bitwise, < 0: one per cell
    matrix_free_batch = 64,         -- cells per block product (default 64)
    -- parallel_threshold = 20000,  -- lts, mlts, li (matrix_free): products on sets of more unknowns on the threads
    matrix_free_precision = "double", -- "single": compute blocks of the cells stored in single precision (products in
                                    --   double, operator perturbed at 1e-7): faster on unstructured meshes (default
                                    --   "double")
    fast = false,                   -- fast mode of every scheme (both formulations); first order: = matrix_free
                                    --   (default false)
}

-- Mesh ------------------------------------------------------------------------------------------------------------------
--  Cartesian mesh of the rectangle x x y, elastic for x < 0, acoustic for x > 0, interface x = 0 (a mesh line).
mesh = {
    x = {-1.0, 1.0},                -- (default {-1, 1})
    y = {0.0, 1.0},                 -- (default {0, 1})
    nx = 4, ny = 2,                 -- cells of the coarsest mesh (default 4, 2)
    level = 0,                      -- uniform refinements: (nx 2^level) x (ny 2^level) cells (default 0)
    refine_points = {{-0.5, 0.5}, {0.5, 0.5}},  -- the cells within one cell size of these points ...
    refine_levels = 0,              -- ... are refined this many times (hanging nodes): fine region of lts, li (default 0)
}

-- Problem ---------------------------------------------------------------------------------------------------------------
problem = {
    solution = "interface_coupling",-- manufactured solution, Dirichlet data on the whole boundary:
                                    --   "interface_coupling": nonzero interface traces (tests the coupling terms)
                                    --   "non_polynomial": vanishes with its gradient on the interface
                                    -- materials: elastic rho = 1, vp = sqrt(3), vs = 1; acoustic rho = 1, vp = 1
                                    -- (default "interface_coupling")
}

-- Study -----------------------------------------------------------------------------------------------------------------
study = {
    type = "space",                 -- "run": one computation on mesh.level, errors at final_time
                                    -- "space": levels mesh.level .. last_level, h and dt halved together
                                    -- "time": mesh.level, n_steps 2^i (i < time_levels) against a reference solution
                                    --   with reference_factor times the finest n_steps (default "run")
                                    -- "cfl": largest stable number of steps; "rho": spectral radius of one step from
                                    --   the energy growth over n_steps steps (random data)
    last_level = 4,                 -- (default mesh.level)
    time_levels = 5,                -- (default 5)
    reference_factor = 4,           -- (default 4)
}

output = {
    file = "",                      -- convergence table; "": name built from the parameters (default "")
    -- profile = false,             -- matrix_free lts / mlts / li: work per macro step and wall times of the parts
    -- snapshot_samples = 0,        -- > 0: snapshots also sampled on an m x m grid of each cell
    -- snapshot_elastic = "velocity", -- samples in the elastic cells: "velocity" (|v|) | "stress" (-tr(sigma)/2, and the
                                    --   shear stress in <table>_snapshot_<i>_shear.dat)
    -- peak_samples = 0,            -- > 0: max over time of |p| / |tr(sigma)|/2 on an m x m grid of each cell
}
