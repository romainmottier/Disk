-- =====================================================================================================================
--  wave_propagation input: elasto-acoustic coupling, SECOND-ORDER formulation
--  (displacement in the solid, velocity potential in the fluid; M u'' + C u' + K u = F, solved as a first-order
--  system with the auxiliary cell velocity v = u_T' and the dynamic interface face unknowns)
--
--      ./wave_propagation second_order.lua ["lua statements executed after this file" ...]
--      e.g. ./wave_propagation second_order.lua "time.scheme = 'sdirk'; space.weight = 1; time.n_steps = 20"
--
--  Every entry is optional: a missing entry takes the default given in the comment; unknown entries are reported.
--  This file is plain Lua: variables, loops and expressions can be used.
--  Annotated template, copied next to the executable by CMake: edit a copy with another name.
-- =====================================================================================================================

formulation = "second_order"        -- "first_order" | "second_order"

-- Space discretization --------------------------------------------------------------------------------------------------
--  Explicit integrators (erk, lts, li): optimal setting of Steins et al. (M2AN 57, 2023): order = "mixed",
--      stabilization = "least_squares" (Lehrenfeld-Schoberl), scaling = "1/h", weight = "auto" (weight_factor in [3, 5]),
--      faces.solver = "splitting"
--  Implicit integrator (sdirk): order = "mixed", stabilization = "least_squares", scaling = "1/h", weight = 1
space = {
    k = 2,                          -- face polynomial degree (default 1)
    order = "mixed",                -- "equal": cell degree k | "mixed": cell degree k+1 (default "equal")
    stabilization = "least_squares",-- "least_squares": plain least squares in equal order, Lehrenfeld-Schoberl in mixed
                                    --   order | "hho": HHO stabilization with its high-order correction (always 1/h)
                                    -- (default "least_squares")
    scaling = "1/h",                -- "1/h" | "1": O(1) weight (default "1")
    weight = "auto",                -- stabilization weight gamma: a number, or "auto" = weight_factor gamma*, gamma* the
                                    --   critical weight of the face splitting computed cell by cell (Steins et al.,
                                    --   Lemma 4.3, Table 1; independent of h) (default 1)
    weight_factor = 3,              -- (default 3)
}

-- Algebraic face unknowns of the explicit (erk, lts) and locally implicit (li) integrators ------------------------------
faces = {
    solver = "splitting",           -- "splitting": D a^{m+1} = r - (K_aa - D) a^m, D = gamma sigma* (block-diagonal),
                                    --   fully explicit, contraction gamma*/gamma (needs gamma > gamma*; impossible in
                                    --   equal order with the "hho" stabilization for k = 1, 3)
                                    -- "direct": sparse factorization of K_aa (default "splitting")
    tolerance = 1e-14,              -- relative increment stopping the splitting (default 1e-14)
    max_iterations = 20000,         -- (default 20000)
}

-- Time discretization ---------------------------------------------------------------------------------------------------
time = {
    scheme = "erk",                 -- "erk": ERK(s) | "sdirk": SDIRK(s, s+1) | "lts": two-level LTS-RK4 | "li":
                                    --   LI-SDIRK(4,3), ERK4 coarse / SDIRK(3,4) fine | "mlts": multilevel LTS-RK4
                                    --   (keys level_ratio, level_ratios, overlap, mlts_refresh of first_order.lua)
                                    --   (default "erk")
    stages = 4,                     -- erk: 1..4, sdirk: 1..3; 0: 4 for erk, 3 for sdirk (default 0)
    final_time = 1.0,               -- (default 1)
    n_steps = 120,                  -- time steps on the first mesh of the study, doubled with each level of a space
                                    --   study; the explicit integrators print omega_max and warn above their stability
                                    --   limit (default 10)
    substeps = 0,                   -- lts: fine substeps, 0: round(h_max / h_min) (default 0)
    fine_size = 0,                  -- lts, li: cells smaller than this (+ one layer of vertex neighbours) are fine;
                                    --   0: 0.75 h_max, everything on a uniform mesh (default 0)
}

-- Linear solvers --------------------------------------------------------------------------------------------------------
solver = {
    static_condensation = false,    -- newmark: static condensation of the cell unknowns (default false)
    fast = false,                   -- fast mode (PARAMETRES.txt, section 7), same results up to round-off (algebraic
                                    --   faces: up to faces.tolerance): parallel assembly by classes of cells, cached F(t)
                                    --   data, parallel products, faces by conjugate gradients (splitting) or Pardiso
                                    --   (direct), stage systems condensed and solved by Pardiso, shared LTS / LI with
                                    --   fewer evaluations (default false)
}

-- Mesh ------------------------------------------------------------------------------------------------------------------
--  Cartesian mesh of the rectangle x x y, elastic for x < 0, acoustic for x > 0, interface x = 0 (a mesh line).
--  In this formulation one local refinement level (h/2) hardly raises the largest frequency: lts and li need
--  refine_levels >= 2 for the fine region to limit the CFL.
mesh = {
    x = {-1.0, 1.0},                -- (default {-1, 1})
    y = {0.0, 1.0},                 -- (default {0, 1})
    nx = 4, ny = 2,                 -- cells of the coarsest mesh (default 4, 2)
    level = 0,                      -- uniform refinements (default 0)
    refine_points = {{-0.5, 0.5}, {0.5, 0.5}},  -- the cells within one cell size of these points ...
    refine_levels = 0,              -- ... are refined this many times (hanging nodes): fine region of lts, li (default 0)
}

-- Problem ---------------------------------------------------------------------------------------------------------------
problem = {
    solution = "interface_coupling",-- "interface_coupling" (nonzero interface traces) | "non_polynomial" (vanishes on
                                    --   the interface); elastic rho = 1, vp = sqrt(3), vs = 1; acoustic rho = 1, vp = 1
                                    -- (default "interface_coupling")
}

-- Study -----------------------------------------------------------------------------------------------------------------
study = {
    type = "space",                 -- "run" | "space" (levels mesh.level .. last_level) | "time" (n_steps 2^i,
                                    --   i < time_levels, against a reference solution) (default "run")
    last_level = 4,                 -- (default mesh.level)
    time_levels = 5,                -- (default 5)
    reference_factor = 4,           -- (default 4)
}

output = {
    file = "",                      -- convergence table; "": name built from the parameters (default "")
}
