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
}

-- Time discretization ---------------------------------------------------------------------------------------------------
time = {
    scheme = "erk",                 -- "erk": ERK(s) | "sdirk": SDIRK(s, s+1) | "lts": two-level LTS-RK4 (substeps in the
                                    --   fine region) | "li": LI-SDIRK(4,3), ERK4 coarse / SDIRK(3,4) fine (default "erk")
    stages = 4,                     -- erk: 1..4, sdirk: 1..3; 0: 4 for erk, 3 for sdirk (default 0)
    final_time = 1.0,               -- (default 1)
    n_steps = 80,                   -- time steps on the first mesh of the study, dt = final_time / n_steps;
                                    --   doubled with each level of a space study (default 10)
    substeps = 0,                   -- lts: fine substeps, 0: round(h_max / h_min) (default 0)
    fine_size = 0,                  -- lts, li: cells smaller than this (+ one layer of vertex neighbours) are fine;
                                    --   0: 0.75 h_max, everything on a uniform mesh (default 0)
}

solver = {
    static_condensation = false,    -- sdirk: static condensation of the cell unknowns in the stage systems (default false)
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
    last_level = 4,                 -- (default mesh.level)
    time_levels = 5,                -- (default 5)
    reference_factor = 4,           -- (default 4)
}

output = {
    file = "",                      -- convergence table; "": name built from the parameters (default "")
}
