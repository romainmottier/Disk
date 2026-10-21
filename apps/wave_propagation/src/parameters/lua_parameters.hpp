//
//  lua_parameters.hpp
//
//  Parameters of the two elasto-acoustic prototypes (HHOFirstOrder, HHOSecondOrder), read from a Lua input file with
//  sol2, as in the other diskpp applications. The script sets the string formulation and fills the tables space,
//  faces, solver, time, mesh, problem, study and output; missing entries keep their default (the values below) and unknown
//  entries are reported. The arguments following the file on the command line are executed as Lua statements after
//  it, e.g.
//      ./wave_propagation first_order.lua "time.scheme = 'sdirk'; time.n_steps = 20"
//  Annotated templates: share/first_order.lua, share/second_order.lua.
//

#pragma once
#ifndef lua_parameters_hpp
#define lua_parameters_hpp

#include "sol/sol.hpp"
#include <array>
#include <cmath>
#include <map>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct wave_parameters {

    std::string formulation = "first_order";            // "first_order" | "second_order"

    struct space_data {
        size_t      k = 1;                              // face degree
        std::string order = "equal";                    // "equal": cells of degree k | "mixed": cells of degree k+1
        std::string stabilization = "least_squares";    // "least_squares" (plain LS / Lehrenfeld-Schoberl) | "hho"
        std::string scaling = "1";                      // "1": O(1) weight | "1/h": weight in 1/h | "graded" (first order):
                                                        //   O(1) weight times max(1, h_max / h_T) | "hybrid" (first
                                                        //   order): 1/h on the cells smaller than 0.75 h_max, O(1) on
                                                        //   the others
        std::string elastic_formulation = "velocity_stress";  // first order, elastic medium: "velocity_stress" (coupled
                                                        //   assembler) | "three_fields" (Burman-Duran-Ern-Steins 2021)
        std::string stabilization_size = "diskpp";      // h of the 1/h weight: "diskpp" (as diskpp: cell diameter for
                                                        //   the scalar HHO stabilization, see the assemblers), "cell"
                                                        //   (cell diameter h_T) | "face" (face diameter h_F)
        double      weight = 1.0;                       // second order: stabilization weight gamma; first order: weight
                                                        //   eta of the stabilization if weight_given (default 0.8
                                                        //   acoustic, 1.5 elastic, as in MEKG)
        bool        weight_given = false;               // space.weight set in the input
        double      weight_acoustic = -1.0;             // first order: weights eta^f, eta^s of the acoustic and elastic
        double      weight_elastic = -1.0;              //   stabilizations (space.weight_acoustic / weight_elastic, they
                                                        //   override space.weight); < 0: not given
        bool        weight_auto = false;                // second order: gamma = weight_factor gamma* (Steins et al.)
        double      weight_factor = 3.0;
        std::string interface_quadrature = "exact";     // elasto-acoustic coupling: "exact" (degree 2k rule) | "legacy"
                                                        //   (degree k rule of the papers' code, rank-deficient coupling:
                                                        //   only to reproduce published results)
    } space;

    struct faces_data {                                 // second order, explicit and locally implicit integrators
        std::string solver = "splitting";               // "splitting" (Steins et al.) | "direct"
        double      tolerance = 1.0e-14;                // relative increment of the splitting
        size_t      max_iterations = 20000;
        bool        localize = false;                   // lts, mlts with F = 0: local evaluations (faces from zero on
                                                        //   the faces reached; B y, B^2 y, B^3 y of the coarse
                                                        //   coefficients truncated to R rings, R from the contraction)
        size_t      localize_margin = 2;                //   rings added to R
        bool        projected = false;                  // lts: splitting B = A + (B - A), A = B P with the faces solved on
                                                        //   the fine region and projected_rings rings only (local
                                                        //   substeps), B - A in the coarse Taylor polynomial
        size_t      projected_rings = 2;
    } faces;

    struct solver_data {
        bool static_condensation = false;               // sdirk (first order), newmark: static condensation of the cells
        bool iterative = false;                         // sdirk (first order): BiCGSTAB with incomplete LU instead of
                                                        //   the direct solver (Pardiso)
        double tolerance = 1.0e-6;                      // relative residual of the iterative solver
        double tolerance_growth = 1.0;                  // space study: tolerance multiplied by this factor from one
                                                        //   level to the next (Mottier-Ern-Guillot: 2^-(k+1))
        bool matrix_free = false;                       // first order: erk, lts, li, cfl study on the matrix-free
                                                        //   explicit operator (time_schemes/mf_erk_operator.hpp, shared
                                                        //   blocks); sdirk: fast mode (parallel products, Pardiso without
                                                        //   iterative refinement); parallel assembly, cached F(t) data
        double matrix_free_tolerance = 1.0e-10;         // local blocks shared when equal to this relative tolerance (0:
                                                        //   equal bit for bit, < 0: one block per cell)
        size_t matrix_free_batch = 64;                  // cells per block product
        size_t parallel_threshold = 20000;              // matrix-free lts / mlts / li: products on sets of more unknowns
                                                        //   than this on the threads, the others on one thread
        std::string matrix_free_precision = "double";   // "single": compute blocks of the cells stored in single
                                                        //   precision (products in double; operator perturbed at 1e-7)
        bool fast = false;                              // fast mode of every scheme (results equal up to round-off):
                                                        //   first order = matrix_free; second order: parallel assembly,
                                                        //   cached F(t) data, parallel products, algebraic faces by
                                                        //   conjugate gradients (splitting) or Pardiso (direct), stage
                                                        //   systems condensed and solved by Pardiso
    } solver;

    struct time_data {
        std::string scheme = "erk";                     // "erk" | "sdirk" | "lts" | "li" | "mlts" (first order) |
                                                        //   "newmark" (second order)
        size_t      stages = 0;                         // erk 1..4, sdirk 1..3; 0: 4 for erk, 3 for sdirk
        double      final_time = 1.0;
        size_t      n_steps = 10;                       // time steps on the first mesh of the study
        double      dt_h = 0.0;                         // > 0: dt = dt_h h_min on each mesh (h_min: smallest cell
                                                        //   diameter), instead of n_steps
        double      steps_growth = 2.0;                 // space study with study.refine_dt: number of steps multiplied
                                                        //   by this factor from one level to the next (rounded; 2: dt
                                                        //   halved with h ; Mottier-Ern-Guillot: 2^((k+1)/(s+1)))
        size_t      substeps = 0;                       // lts: fine substeps, 0: round(h_max / h_min)
        double      fine_size = 0.0;                    // lts, li: cells smaller than this are fine, 0: automatic
        std::string fine_measure = "diameter";          // lts, li: size of a cell compared with fine_size: its diameter or its
                                                        //   smallest width ("width": thin cut cells, slivers)
        size_t      fine_layers = 1;                    // lts, li: rings of vertex neighbours added to the fine region
        double      newmark_beta = 0.25;                // newmark: parameters beta, gamma
        double      newmark_gamma = 0.5;
        size_t      level_ratio = 2;                    // mlts: ratio r of the time steps of two consecutive levels;
                                                        //   level of a cell ceil(log_r(h_max / h_T))
        std::vector<size_t> level_ratios;               // mlts: ratios p_0, p_1, ... of the levels (coarsest first, the
                                                        //   last one repeated) instead of level_ratio; level of a cell:
                                                        //   first l with p_0 ... p_{l-1} >= h_max / h_T
        std::string level_measure = "h";                // mlts: level of a cell from h_max / h_T ("h") or from
                                                        //   max(h / c) / (h_T / c_T) ("h/c": heterogeneous media)
        std::string level_size = "diameter";            // mlts: size of a cell in its level measure: its diameter,
                                                        //   4 |T| / |dT| ("area") or its smallest width ("width": the
                                                        //   measure that sets the stable step of thin or cut cells)
        double      level_tolerance = 0.0;              // mlts: measure of a cell divided by 1 + level_tolerance before
                                                        //   its level is chosen (media of nearly equal speeds, e.g. soft
                                                        //   tissues 1.50 / 1.55, stay on the same level)
        size_t      overlap = 1;                        // mlts: rings of vertex neighbours annexed by each level from
                                                        //   the coarser ones (overlap region of the preprint)
        std::vector<size_t> overlaps;                   // mlts: rings annexed by the levels 1, 2, ... (finest last, the
                                                        //   last one repeated) instead of overlap
        bool        mlts_reuse = false;                 // mlts (matrix-free): own coefficients of a level from the chain
                                                        //   a_{i+1} = B a_i on its deep unknowns, masked products on the
                                                        //   boundary layer only (same scheme, about 4 instead of 7
                                                        //   products per unknown of a thick level)
        std::string li_coupling = "taylor";             // li (matrix-free): coarse values in the implicit stages, the
                                                        //   Taylor polynomial at t_n + c_i dt ("taylor") or the stage
                                                        //   values of the SDIRK tableau, y_n + dt sum_l a_il Y'(c_l dt)
                                                        //   ("stage": consistent with the fine stages, no order
                                                        //   reduction from the coarse / fine interface)
        bool        li_condensation = false;            // li (matrix-free): static condensation of the fine cell
                                                        //   unknowns, the stage systems solved on the fine faces only
        std::string li_system = "dynamic";              // li (matrix-free): fine face unknowns inside the fine region
                                                        //   by their evolution equation, k_f = G k_c ("dynamic"), or by
                                                        //   their static constraint u_f = G u_c, with the stage slopes
                                                        //   (U_i - R_i) / (gamma dt) ("constraint": the same scheme, well
                                                        //   conditioned when the fine region is very stiff)
        std::string mlts_refresh = "window";            // mlts: own Taylor coefficients of a level computed once per
                                                        //   window of the level ("window", Algorithm 4) or at each
                                                        //   of its sub-windows ("substep", prototypes)
    } time;

    struct mesh_data {
        std::string type = "cartesian";                 // "cartesian" | "fvca5" (files of the FVCA5 benchmark format) |
                                                        //   "polymesher" (PolyMesher meshes, polygon_2d_mesh_reader) |
                                                        //   "polygons" (polygonal meshes with one integer tag per cell,
                                                        //   see mesh/domain.hpp)
        std::array<double, 2> x = {-1.0, 1.0};          // cartesian: rectangle [x0,x1]x[y0,y1]
        std::array<double, 2> y = {0.0, 1.0};
        size_t nx = 4, ny = 2;                          // cartesian: cells of the coarsest mesh
        size_t level = 0;                               // cartesian: uniform refinements; fvca5, polymesher: index in
                                                        //   files (>= 1)
        std::vector<std::array<double, 2>> refine_points = {{-0.5, 0.5}, {0.5, 0.5}};
        size_t refine_levels = 0;                       // cartesian: local refinements around refine_points
        std::vector<double> refine_box;                 // cartesian: {x0, x1, y0, y1}: cells with barycenter in the box
                                                        //   refined refine_levels times, graded transition (one level
                                                        //   less per ring of neighbours), instead of refine_points
        std::vector<std::string> files;                 // fvca5, polymesher: files[l] is the mesh of level l (Lua
                                                        //   indexing from 1)
    } mesh;

    struct problem_data {
        std::string medium = "coupled";                 // "coupled": elastic x < 0, acoustic x > 0 | "acoustic" | "elastic"
        std::string solution = "interface_coupling";    // manufactured solution, see solution_names()
        std::string dirichlet = "exact";                // Dirichlet data: "exact" (the solution) | "zero" (homogeneous
                                                        //   data, the initial data remaining those of the solution)
        std::string source = "solution";                // first order: "solution" (source of the manufactured solution)
                                                        //   | "zero"
        std::string interface = "x";                    // coupled: "x" (elastic x < 0, acoustic x > 0), "y" (elastic
                                                        //   y < 0, acoustic y > 0) | "tags" (medium and material of each
                                                        //   cell from its mesh tag, problem.zones)
        double rho_acoustic = 1.0, vp_acoustic = 1.0;   // problem.acoustic = {rho = .., vp = ..}
        double rho_elastic = 1.0, vp_elastic = std::sqrt(3.0), vs_elastic = 1.0;   // problem.elastic = {rho, vp, vs}
        struct zone { std::string medium; double rho = 1.0, vp = 1.0, vs = 0.0; };
        std::map<int, zone> zones;                      // interface = "tags": problem.zones = {[tag] = {medium =
                                                        //   "acoustic" | "elastic", rho = .., vp = .., vs = ..}, ...}
        double omega = 1.0, theta = 10.0;               // coupled_sinusoidal: spatial and temporal frequencies
        std::array<double, 2> lshape_corner = {0.5, 0.5};   // lshape_corner(_sin): reentrant corner of the L-shaped domain
        double lshape_exponent = 2.0 / 3.0;             // lshape_corner(_sin): exponent nu of phi = r^nu sin(nu theta)
                                                        //   (2/3: first singular function, 4/3: second one)
        struct ricker_data {                            // solution = "ricker": initial wavelet theta exp(-pi^2 r^2 /
            std::string field = "acoustic_flux";        //   lambda^2) (x - xc, y - yc) of the acoustic flux or of the
            std::array<double, 2> center = {0.0, 0.125};//   elastic velocity, or theta exp(-pi^2 r^2 / lambda^2) of the
            double frequency = 10.0;                    //   acoustic pressure ("acoustic_pressure"), lambda = velocity /
            double amplitude = 10.0;                    //   frequency
            double velocity = 1.0;
            std::string shape = "gaussian";             // "focused": converging tone burst of a focused transducer (arc of
            double radius = 0.0;                        //   radius `radius` around the focus `center`, axis `axis` degrees
            double axis = 0.0;                          //   from the focus to the transducer, half aperture `aperture`
            double aperture = 30.0;                     //   degrees, `cycles` cycles at `frequency`, Gaussian envelope),
            double cycles = 3.0;                        //   pressure and acoustic flux of a one-way wave towards the focus
            double density = 1.0;                       //   (impedance density x velocity of the medium of the transducer)
        } ricker;
    } problem;

    struct study_data {
        std::string type = "run";                       // "run" | "space" | "time"
        size_t last_level = 0;                          // space: levels mesh.level .. last_level
        bool   refine_dt = true;                        // space: dt halved with each level (true) or fixed (false)
        size_t time_levels = 5;                         // time: n_steps 2^i, i < time_levels
        std::string reference = "numerical";            // time: errors against a reference solution of the same
                                                        //   integrator ("numerical") or against the exact solution only
        size_t reference_factor = 4;                    // time: reference with this many times the finest n_steps
        double cfl_threshold = 0.05;                    // cfl: unstable when the energy grows by more than this fraction
                                                        //   of the initial or of the previous energy (MEG, Algorithm 1)
        double cfl_decrement = 0.01;                    // cfl: number of steps decreased by this fraction (+ 1) after
                                                        //   each stable run
        size_t cfl_runs = 2000;                         // cfl: at most this many runs
        bool cfl_random = false;                        // cfl: random initial cell data (all the modes) instead of the
                                                        //   initial data of the solution
        bool cfl_absolute = false;                      // cfl: relative variations |E^n - E^0| / E^0 and
                                                        //   |E^n - E^n-1| / E^n-1 (text of Algorithm 1) instead of the
                                                        //   signed growth (code of the paper)
    } study;

    struct output_data {
        std::string file = "";                          // convergence table, "" = name built from the parameters
        std::string table = "";                         // the same errors without orders, one row per level, columns
                                                        //   named in a "#" header line; "" = file with extension .dat
        bool energy = false;                            // discrete energy history (file _energy.dat): newmark, first
                                                        //   order (1/2 x^T M x on the cell unknowns, total, acoustic
                                                        //   and elastic parts)
        struct sensor { std::string name; std::array<double, 2> point = {0.0, 0.0}; std::vector<std::string> fields; };
        std::vector<sensor> sensors;                    // first order: output.sensors = {{name = "S", x = .., y = ..,
                                                        //   fields = {..}}, ...}, fields among sensor_fields(); one file
                                                        //   <table>_<name>.csv per sensor
        size_t sensor_every = 1;                        // sensors recorded every sensor_every time steps
        size_t snapshot_samples = 0;                    // first order: > 0: snapshots also sampled on an m x m grid of
                                                        //   the bounding box of each cell (<table>_snapshot_<i>_samples.dat)
        std::string snapshot_elastic = "velocity";      // samples in the elastic cells: |velocity| ("velocity") or the
                                                        //   mean stress -tr(sigma)/2 ("stress", the pressure of the
                                                        //   solid), with the shear stress in <table>_snapshot_<i>_shear.dat
        size_t peak_samples = 0;                        // first order: > 0: max over time (every sensor_every steps) of
                                                        //   |p| (acoustic) and |tr(sigma)|/2 (elastic) on an m x m grid
                                                        //   of each cell (<table>_peak_samples.dat)
        std::vector<double> snapshots;                  // first order: times of the snapshots (cell values at the
                                                        //   nearest time step, files <table>_snapshot_<i>.dat)
        bool profile = false;                           // matrix-free lts / mlts / li: work per macro step (unknowns x
                                                        //   products) and wall times of the parts of the integrator
        bool cpu = false;                               // space study: columns cpu (seconds of assembly and time
                                                        //   integration of each level), dofs and face_dofs (unknowns
                                                        //   before and after the static condensation) in the tables
        std::string errors = "fields";                  // first order: "fields" (L2 errors of each unknown) or
                                                        //   "energy" (components of the norms of Mottier-Ern-Khot-
                                                        //   Guillot, weighted by 1/kappa, rho, rho, C^-1: cell unknowns
                                                        //   p, v against the projections, dual unknowns m, s, against
                                                        //   the projections, plain L2 norms of s; HHOFirstOrder.hpp)
    } output;

    // Fields of a sensor: cell values at the point (pressure, flux, velocity, stress) and, at an interface point, the
    // face unknowns of the acoustic pressure and of the elastic velocity
    static const std::set<std::string>& sensor_fields() {
        static const std::set<std::string> f = {"acoustic_pressure", "acoustic_flux", "elastic_velocity", "elastic_stress",
                                                "interface_pressure", "interface_velocity"};
        return f;
    }

    // Manufactured solutions: name, medium ("" = any), description
    static const std::vector<std::array<std::string, 3>>& solution_names() {
        static const std::vector<std::array<std::string, 3>> names = {
            {"interface_coupling", "", "EFunctionInterfaceCoupling: coupled solution with nonzero interface traces"},
            {"non_polynomial", "", "EFunctionNonPolynomial: vanishes with its gradient on x = 0"},
            {"acoustic_quadratic_in_time", "acoustic", "p = t^2 sin(pi x) sin(pi y), Burman-Duran-Ern 2022 eq. (60)"},
            {"acoustic_quadratic_in_space", "acoustic", "p = sin(sqrt2 pi t) x(1-x) y(1-y), eq. (61)"},
            {"acoustic_non_polynomial", "acoustic", "p = sin(sqrt2 pi t) sin(pi x) sin(pi y) / (sqrt2 pi), eq. (62)"},
            {"acoustic_quadratic_in_space_cos", "acoustic", "p = cos(sqrt2 pi t) x(1-x) y(1-y): eq. (61) with cos (published Fig. 1 right)"},
            {"acoustic_non_polynomial_cos", "acoustic", "p = cos(sqrt2 pi t) sin(pi x) sin(pi y) / (sqrt2 pi): eq. (62) with cos"},
            {"lshape_corner", "acoustic", "u = cos(2 t) r^(2/3) sin(2 theta / 3) on the L-shaped domain (0,1)^2 \\ [0.5,1)x(0.5,1] (corner (0.5,0.5))"},
            {"lshape_corner_sin", "acoustic", "u = sin(2 t) r^(2/3) sin(2 theta / 3) on the L-shaped domain: v(0) = 0, compatible initial data"},
            {"elastic_quadratic_in_time", "elastic", "u = t^2 (-sin(pi x) cos(pi y), cos(pi x) sin(pi y)), Burman-Duran-Ern-Steins 2021 eq. (78)"},
            {"elastic_quadratic_in_space", "elastic", "u = sin(sqrt2 pi t) x(1-x) y(1-y) (1,1), eq. (79)"},
            {"elastic_non_polynomial", "elastic", "u = sin(sqrt2 pi t) (-sin(pi x) cos(pi y), cos(pi x) sin(pi y)), eq. (80)"},
            {"coupled_sinusoidal", "", "Mottier-Ern-Guillot: u = x^2 sin(omega pi x) sin(omega pi y) sin(theta pi t), u_x = u_y = x^2 cos(omega pi x / 2) sin(omega pi y) cos(theta pi t) (problem.omega, problem.theta)"},
            {"coupled_poly_in_space", "", "Mottier-Ern-Khot-Guillot: u = (1-x) x^2 (1-y) y sin(sqrt2 pi t), u_x = u_y = (1+x) x^2 (1-y) y sin(sqrt2 pi t)"},
            {"coupled_poly_in_time", "", "Mottier-Ern-Khot-Guillot: u = u_x = u_y = x sin(pi x) sin(pi y) t^2"},
            {"ricker", "", "no source, homogeneous Dirichlet data, initial Ricker wavelet (problem.ricker)"}};
        return names;
    }

    bool explicit_scheme() const {
        return time.scheme == "erk" || time.scheme == "lts" || time.scheme == "li" || time.scheme == "mlts";
    }
    bool mixed_order() const { return space.order == "mixed"; }
    bool scaled_stabilization() const { return space.scaling == "1/h"; }
    int  stabilization_size() const { return space.stabilization_size == "cell" ? 1 : space.stabilization_size == "face" ? 2 : 0; }
    bool least_squares() const { return space.stabilization == "least_squares"; }
    // First order: stabilization weights eta^f (acoustic) and eta^s (elastic): space.weight_acoustic / weight_elastic,
    // else space.weight, else the MEKG values 0.8 / 1.5
    double acoustic_weight() const { return space.weight_acoustic > 0.0 ? space.weight_acoustic : space.weight_given ? space.weight : 0.8; }
    double elastic_weight() const { return space.weight_elastic > 0.0 ? space.weight_elastic : space.weight_given ? space.weight : 1.5; }
    bool has_elastic() const { return problem.medium != "acoustic"; }
    bool has_acoustic() const { return problem.medium != "elastic"; }

    size_t stages() const {
        if (time.stages > 0) return time.stages;
        return time.scheme == "sdirk" ? 3 : 4;
    }

    std::string scheme_label() const {
        if (time.scheme == "erk")   return "ERK(" + std::to_string(stages()) + ")";
        if (time.scheme == "sdirk") return "SDIRK(" + std::to_string(stages()) + "," + std::to_string(stages() + 1) + ")";
        if (time.scheme == "lts")   return "LTS-RK4";
        if (time.scheme == "mlts")  return "MLTS-RK4";
        if (time.scheme == "newmark") {
            std::ostringstream os;
            os << "Newmark(beta = " << time.newmark_beta << ", gamma = " << time.newmark_gamma << ")";
            return os.str();
        }
        return "LI-SDIRK(4,3)";
    }

    std::string space_label() const {
        std::string s = (mixed_order() ? "mixed order (k+1,k)" : "equal order (k,k)");
        s += least_squares() ? (mixed_order() ? ", Lehrenfeld-Schoberl" : ", plain least squares") : ", HHO stabilization";
        s += scaled_stabilization() ? ", 1/h" : ", O(1)";
        if (scaled_stabilization() && space.stabilization_size != "diskpp") s += " (h = " + space.stabilization_size + " diameter)";
        if (space.elastic_formulation == "three_fields") s += ", three-field elastic formulation";
        return s;
    }

    // Default name of the convergence table
    std::string output_file() const {
        if (!output.file.empty()) return output.file;
        std::ostringstream os;
        os << formulation << "_" << time.scheme;
        if (time.scheme == "erk" || time.scheme == "sdirk") os << stages();
        os << "_k" << space.k << "_" << space.order << "_" << study.type << "_l" << mesh.level << "_p" << mesh.refine_levels << ".txt";
        return os.str();
    }

    // Name of the data table: output.table, or output_file() with the extension .dat
    std::string table_file() const {
        if (!output.table.empty()) return output.table;
        std::string f = output_file();
        size_t dot = f.find_last_of('.');
        size_t slash = f.find_last_of('/');
        if (dot != std::string::npos && (slash == std::string::npos || dot > slash))
            f = f.substr(0, dot);
        return f + ".dat";
    }

    std::string energy_file() const {
        std::string f = table_file();
        return f.substr(0, f.size() - 4) + "_energy.dat";
    }

    void print(std::ostream& os) const {
        os << "   formulation : " << formulation << std::endl;
        os << "   space       : k = " << space.k << ", " << space_label();
        if (formulation == "second_order") {
            if (space.weight_auto) os << ", weight = " << space.weight_factor << " gamma*";
            else                   os << ", weight = " << space.weight;
        }
        else if (space.weight_acoustic > 0.0 || space.weight_elastic > 0.0)
            os << ", weights " << acoustic_weight() << " (acoustic) / " << elastic_weight() << " (elastic)";
        else
            os << ", weight " << (space.weight_given ? std::to_string(space.weight) : std::string("0.8 / 1.5 (MEKG)"));
        if (space.interface_quadrature == "legacy")
            os << ", legacy interface quadrature (degree k)";
        os << std::endl;
        if (formulation == "second_order" && explicit_scheme())
            os << "   faces       : " << faces.solver << (faces.solver == "splitting" ? ", tolerance " : "")
               << (faces.solver == "splitting" ? std::to_string(faces.tolerance) : "")
               << (faces.localize ? ", local evaluations (lts, mlts), margin " + std::to_string(faces.localize_margin) : "")
               << (faces.projected ? ", projected condensation (lts), " + std::to_string(faces.projected_rings) + " ring(s)" : "")
               << std::endl;
        if ((formulation == "first_order" && time.scheme == "sdirk") || time.scheme == "newmark")
            os << "   solver      : static condensation " << (solver.static_condensation ? "yes" : "no")
               << (solver.iterative ? ", BiCGSTAB-ILU, tolerance " + std::to_string(solver.tolerance) : std::string(""))
               << std::endl;
        if (solver.fast && formulation == "second_order")
            os << "   solver      : fast mode (parallel assembly, cached F(t), "
               << (solver.matrix_free ? "matrix-free products" : "parallel products") << ", faces: "
               << (faces.solver == "splitting" ? "conjugate gradients preconditioned by the splitting" : "Pardiso")
               << ", stages: condensed, Pardiso)" << std::endl;
        else if (solver.matrix_free && time.scheme == "sdirk")
            os << "   solver      : fast mode (parallel products, Pardiso without iterative refinement, parallel assembly)"
               << std::endl;
        else if (solver.matrix_free) {
            std::ostringstream tol;
            tol << solver.matrix_free_tolerance;
            os << "   solver      : matrix-free explicit operator, blocks shared "
               << (solver.matrix_free_tolerance < 0.0 ? std::string("never") : "to " + tol.str()) << ", " << solver.matrix_free_batch
               << " cells per product" << (solver.matrix_free_precision == "single" ? ", single precision blocks" : "") << std::endl;
        }
        os << "   time        : " << scheme_label() << ", final time " << time.final_time << ", ";
        if (time.dt_h > 0.0) os << "dt = " << time.dt_h << " h_min";
        else                 os << time.n_steps << " steps";
        if (time.scheme == "lts" || time.scheme == "li") {
            os << ", substeps " << (time.substeps > 0 ? std::to_string(time.substeps) : "auto")
               << ", fine size " << (time.fine_size > 0 ? std::to_string(time.fine_size) : "auto");
        }
        if (time.scheme == "mlts") {
            if (time.level_ratios.empty()) os << ", level ratio " << time.level_ratio;
            else {
                os << ", level ratios";
                for (size_t r : time.level_ratios) os << " " << r;
            }
            if (time.level_measure != "h") os << " (levels from " << time.level_measure << ")";
            if (time.level_size == "area") os << " (cell size 4|T|/|dT|)";
            if (time.level_size == "width") os << " (cell size: smallest width)";
            if (time.level_tolerance > 0) os << " (tolerance " << time.level_tolerance << ")";
            if (time.overlaps.empty()) os << ", overlap " << time.overlap << " ring(s)";
            else {
                os << ", overlaps";
                for (size_t o : time.overlaps) os << " " << o;
            }
            os << ", coefficients per " << (time.mlts_refresh == "window" ? "window" : "sub-window");
            if (time.mlts_reuse) os << " (chain reused)";
        }
        os << std::endl;
        if (mesh.type != "cartesian") {
            os << "   mesh        : " << mesh.type << ", " << mesh.files.size() << " files, level " << mesh.level << std::endl;
        }
        else {
            os << "   mesh        : [" << mesh.x[0] << "," << mesh.x[1] << "]x[" << mesh.y[0] << "," << mesh.y[1] << "], "
               << mesh.nx << "x" << mesh.ny << " cells, level " << mesh.level << ", " << mesh.refine_levels;
            if (mesh.refine_box.empty()) os << " local refinement levels around " << mesh.refine_points.size() << " points" << std::endl;
            else os << " graded refinement levels in the box [" << mesh.refine_box[0] << "," << mesh.refine_box[1] << "]x["
                    << mesh.refine_box[2] << "," << mesh.refine_box[3] << "]" << std::endl;
        }
        os << "   problem     : " << problem.medium << ", " << problem.solution
           << (problem.dirichlet == "zero" ? ", homogeneous Dirichlet data" : "")
           << (problem.source == "zero" ? ", no source" : "") << std::endl;
        if (problem.solution == "coupled_sinusoidal")
            os << "                 omega = " << problem.omega << ", theta = " << problem.theta << std::endl;
        if (problem.medium == "coupled" || problem.interface == "tags") {
            if (problem.interface == "tags") {
                os << "   materials   : zones of the mesh tags";
                for (auto& z : problem.zones)
                    os << ", [" << z.first << "] " << z.second.medium << " (rho " << z.second.rho << ", vp " << z.second.vp
                       << (z.second.medium == "elastic" ? ", vs " + std::to_string(z.second.vs) : std::string("")) << ")";
                os << std::endl;
            }
            else
                os << "   materials   : interface " << problem.interface << " = 0, elastic " << problem.interface
                   << " < 0 (rho " << problem.rho_elastic << ", vp " << problem.vp_elastic << ", vs " << problem.vs_elastic
                   << "), acoustic " << problem.interface << " > 0 (rho " << problem.rho_acoustic << ", vp " << problem.vp_acoustic
                   << ")" << std::endl;
        }
        if (problem.solution == "ricker" && problem.ricker.shape == "focused")
            os << "   ricker      : focused tone burst, focus (" << problem.ricker.center[0] << ", " << problem.ricker.center[1]
               << "), radius " << problem.ricker.radius << ", axis " << problem.ricker.axis << " deg, half aperture "
               << problem.ricker.aperture << " deg, " << problem.ricker.cycles << " cycles at " << problem.ricker.frequency
               << ", amplitude " << problem.ricker.amplitude << ", velocity " << problem.ricker.velocity << ", density "
               << problem.ricker.density << std::endl;
        else if (problem.solution == "ricker")
            os << "   ricker      : " << problem.ricker.field << ", center (" << problem.ricker.center[0] << ", "
               << problem.ricker.center[1] << "), frequency " << problem.ricker.frequency << ", amplitude "
               << problem.ricker.amplitude << ", velocity " << problem.ricker.velocity << std::endl;
        os << "   study       : " << study.type;
        if (study.type == "space") os << ", levels " << mesh.level << ".." << study.last_level << (study.refine_dt ? ", dt halved with h" : ", fixed dt");
        if (study.type == "spectrum") os << ", levels " << mesh.level << ".." << study.last_level;
        if (study.type == "time")  {
            os << ", " << study.time_levels << " time levels, errors against the " << (study.reference == "exact" ? "exact solution" : "reference solution");
            if (study.reference == "numerical") os << " (x" << study.reference_factor << " steps)";
        }
        os << std::endl;
        if (study.type == "cfl")
            os << "   cfl         : energy " << (study.cfl_absolute ? "variation" : "growth") << " > " << study.cfl_threshold << ", steps decreased by "
               << study.cfl_decrement << " + 1, at most " << study.cfl_runs << " runs" << std::endl;
        os << "   output      : " << output_file() << ", " << table_file() << (output.energy ? ", " + energy_file() : "")
           << (output.errors == "energy" ? ", components of the errors in the energy norms" : "") << std::endl;
        for (auto& sn : output.sensors) {
            os << "   sensor      : " << sn.name << " (" << sn.point[0] << ", " << sn.point[1] << "):";
            for (auto& f : sn.fields) os << " " << f;
            os << std::endl;
        }
        if (!output.snapshots.empty()) {
            os << "   snapshots   :";
            for (auto t : output.snapshots) os << " " << t;
            os << std::endl;
        }
    }
};

namespace lua_parameters_detail {

    inline void check_keys(const sol::table& t, const std::string& name, const std::set<std::string>& allowed) {
        for (const auto& kv : t) {
            if (kv.first.get_type() != sol::type::string)
                continue;
            std::string key = kv.first.as<std::string>();
            if (!allowed.count(key))
                std::cout << "   WARNING: unknown parameter " << name << "." << key << " (ignored)" << std::endl;
        }
    }

    inline double get_number(const sol::table& t, const std::string& name, const std::string& key, double def) {
        sol::object o = t[key];
        if (!o.valid()) return def;
        if (o.get_type() != sol::type::number)
            throw std::runtime_error(name + "." + key + " must be a number");
        return o.as<double>();
    }

    inline size_t get_count(const sol::table& t, const std::string& name, const std::string& key, size_t def) {
        double v = get_number(t, name, key, double(def));
        if (v < 0 || std::floor(v) != v)
            throw std::runtime_error(name + "." + key + " must be a nonnegative integer");
        return size_t(v);
    }

    inline bool get_bool(const sol::table& t, const std::string& name, const std::string& key, bool def) {
        sol::object o = t[key];
        if (!o.valid()) return def;
        if (o.get_type() != sol::type::boolean)
            throw std::runtime_error(name + "." + key + " must be true or false");
        return o.as<bool>();
    }

    inline std::string get_string(const sol::table& t, const std::string& name, const std::string& key,
                                  const std::string& def, const std::set<std::string>& values) {
        sol::object o = t[key];
        if (!o.valid()) return def;
        if (o.get_type() != sol::type::string)
            throw std::runtime_error(name + "." + key + " must be a string");
        std::string s = o.as<std::string>();
        if (!values.empty() && !values.count(s)) {
            std::string list;
            for (auto& v : values) list += " \"" + v + "\"";
            throw std::runtime_error(name + "." + key + " = \"" + s + "\" is not one of" + list);
        }
        return s;
    }

    inline std::array<double, 2> get_pair(const sol::table& t, const std::string& name, const std::string& key,
                                          const std::array<double, 2>& def) {
        sol::object o = t[key];
        if (!o.valid()) return def;
        if (o.get_type() != sol::type::table)
            throw std::runtime_error(name + "." + key + " must be a table {a, b}");
        sol::table p = o.as<sol::table>();
        sol::object a = p[1], b = p[2];
        if (!a.valid() || !b.valid() || a.get_type() != sol::type::number || b.get_type() != sol::type::number)
            throw std::runtime_error(name + "." + key + " must be a table of two numbers {a, b}");
        return {a.as<double>(), b.as<double>()};
    }

    inline sol::table get_table(sol::state& lua, const std::string& name) {
        sol::object o = lua[name];
        if (!o.valid())
            return lua.create_table();
        if (o.get_type() != sol::type::table)
            throw std::runtime_error(name + " must be a table");
        return o.as<sol::table>();
    }
}

// Reads the input file (then the extra Lua statements) and checks the parameters
inline wave_parameters load_wave_parameters(const std::string& file, const std::vector<std::string>& statements = {}) {

    using namespace lua_parameters_detail;
    sol::state lua;
    lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::string, sol::lib::table);
    for (auto name : {"space", "faces", "solver", "time", "mesh", "problem", "study", "output"})
        lua[name] = lua.create_table();
    lua.script_file(file);
    for (auto& s : statements)
        lua.script(s);

    wave_parameters p;
    sol::object f = lua["formulation"];
    if (!f.valid() || f.get_type() != sol::type::string)
        throw std::runtime_error("formulation must be set to \"first_order\" or \"second_order\"");
    p.formulation = f.as<std::string>();
    if (p.formulation != "first_order" && p.formulation != "second_order")
        throw std::runtime_error("formulation = \"" + p.formulation + "\" is not \"first_order\" or \"second_order\"");

    sol::table space = get_table(lua, "space");
    check_keys(space, "space", {"k", "order", "stabilization", "scaling", "stabilization_size", "elastic_formulation", "weight",
                                "weight_factor", "weight_acoustic", "weight_elastic", "interface_quadrature"});
    p.space.k             = get_count(space, "space", "k", p.space.k);
    p.space.order         = get_string(space, "space", "order", p.space.order, {"equal", "mixed"});
    p.space.stabilization = get_string(space, "space", "stabilization", p.space.stabilization, {"least_squares", "hho"});
    p.space.scaling       = get_string(space, "space", "scaling", p.space.scaling, {"1", "1/h", "graded", "hybrid"});
    p.space.stabilization_size = get_string(space, "space", "stabilization_size", p.space.stabilization_size, {"diskpp", "cell", "face"});
    p.space.elastic_formulation = get_string(space, "space", "elastic_formulation", p.space.elastic_formulation, {"velocity_stress", "three_fields"});
    sol::object w = space["weight"];
    if (w.valid()) {
        p.space.weight_given = true;
        if (w.get_type() == sol::type::string && w.as<std::string>() == "auto")
            p.space.weight_auto = true;
        else if (w.get_type() == sol::type::number)
            p.space.weight = w.as<double>();
        else
            throw std::runtime_error("space.weight must be a number or \"auto\"");
    }
    p.space.weight_factor = get_number(space, "space", "weight_factor", p.space.weight_factor);
    p.space.weight_acoustic = get_number(space, "space", "weight_acoustic", p.space.weight_acoustic);
    p.space.weight_elastic  = get_number(space, "space", "weight_elastic", p.space.weight_elastic);
    p.space.interface_quadrature = get_string(space, "space", "interface_quadrature", p.space.interface_quadrature,
                                              {"exact", "legacy"});

    sol::table faces = get_table(lua, "faces");
    check_keys(faces, "faces", {"solver", "tolerance", "max_iterations", "localize", "localize_margin", "projected", "projected_rings"});
    p.faces.solver         = get_string(faces, "faces", "solver", p.faces.solver, {"splitting", "direct"});
    p.faces.tolerance      = get_number(faces, "faces", "tolerance", p.faces.tolerance);
    p.faces.max_iterations = get_count(faces, "faces", "max_iterations", p.faces.max_iterations);
    p.faces.localize        = get_bool(faces, "faces", "localize", p.faces.localize);
    p.faces.localize_margin = get_count(faces, "faces", "localize_margin", p.faces.localize_margin);
    p.faces.projected       = get_bool(faces, "faces", "projected", p.faces.projected);
    p.faces.projected_rings = get_count(faces, "faces", "projected_rings", p.faces.projected_rings);

    sol::table solver = get_table(lua, "solver");
    check_keys(solver, "solver", {"static_condensation", "iterative", "tolerance", "tolerance_growth", "matrix_free",
                                  "matrix_free_tolerance", "matrix_free_batch", "matrix_free_precision", "fast",
                                  "parallel_threshold"});
    p.solver.static_condensation = get_bool(solver, "solver", "static_condensation", p.solver.static_condensation);
    p.solver.iterative = get_bool(solver, "solver", "iterative", p.solver.iterative);
    p.solver.tolerance = get_number(solver, "solver", "tolerance", p.solver.tolerance);
    p.solver.tolerance_growth = get_number(solver, "solver", "tolerance_growth", p.solver.tolerance_growth);
    p.solver.matrix_free = get_bool(solver, "solver", "matrix_free", p.solver.matrix_free);
    p.solver.matrix_free_tolerance = get_number(solver, "solver", "matrix_free_tolerance", p.solver.matrix_free_tolerance);
    p.solver.matrix_free_batch = get_count(solver, "solver", "matrix_free_batch", p.solver.matrix_free_batch);
    p.solver.matrix_free_precision = get_string(solver, "solver", "matrix_free_precision", p.solver.matrix_free_precision,
                                                {"double", "single"});
    p.solver.fast = get_bool(solver, "solver", "fast", p.solver.fast);
    p.solver.parallel_threshold = get_count(solver, "solver", "parallel_threshold", p.solver.parallel_threshold);
    if (p.solver.fast && p.formulation == "first_order")
        p.solver.matrix_free = true;
    if (p.solver.matrix_free && p.formulation == "second_order")     // matrix-free products of the fast mode
        p.solver.fast = true;

    sol::table time = get_table(lua, "time");
    check_keys(time, "time", {"scheme", "stages", "final_time", "n_steps", "dt_h", "substeps", "fine_size", "fine_layers", "newmark_beta",
                              "newmark_gamma", "steps_growth", "level_ratio", "level_ratios", "level_measure", "level_tolerance", "overlap",
                              "mlts_refresh", "li_system", "li_coupling", "li_condensation", "overlaps", "mlts_reuse",
                              "level_size", "fine_measure"});
    p.time.scheme        = get_string(time, "time", "scheme", p.time.scheme, {"erk", "sdirk", "lts", "li", "mlts", "newmark"});
    p.time.stages        = get_count(time, "time", "stages", p.time.stages);
    p.time.final_time    = get_number(time, "time", "final_time", p.time.final_time);
    p.time.n_steps       = get_count(time, "time", "n_steps", p.time.n_steps);
    p.time.dt_h          = get_number(time, "time", "dt_h", p.time.dt_h);
    p.time.steps_growth  = get_number(time, "time", "steps_growth", p.time.steps_growth);
    p.time.substeps      = get_count(time, "time", "substeps", p.time.substeps);
    p.time.fine_size     = get_number(time, "time", "fine_size", p.time.fine_size);
    p.time.fine_layers   = get_count(time, "time", "fine_layers", p.time.fine_layers);
    p.time.newmark_beta  = get_number(time, "time", "newmark_beta", p.time.newmark_beta);
    p.time.newmark_gamma = get_number(time, "time", "newmark_gamma", p.time.newmark_gamma);
    p.time.level_ratio   = get_count(time, "time", "level_ratio", p.time.level_ratio);
    p.time.overlap       = get_count(time, "time", "overlap", p.time.overlap);
    p.time.level_measure = get_string(time, "time", "level_measure", p.time.level_measure, {"h", "h/c"});
    p.time.level_tolerance = get_number(time, "time", "level_tolerance", p.time.level_tolerance);
    p.time.level_size = get_string(time, "time", "level_size", p.time.level_size, {"diameter", "area", "width"});
    p.time.fine_measure = get_string(time, "time", "fine_measure", p.time.fine_measure, {"diameter", "width"});
    if (p.time.level_tolerance < 0.0) throw std::runtime_error("time.level_tolerance must be nonnegative");
    sol::object ratios = time["level_ratios"];
    if (ratios.valid()) {
        if (ratios.get_type() != sol::type::table)
            throw std::runtime_error("time.level_ratios must be a table of integers {p_0, p_1, ...}");
        p.time.level_ratios.clear();
        sol::table t = ratios.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            if (o.get_type() != sol::type::number || o.as<double>() != std::floor(o.as<double>()) || o.as<double>() < 2.0)
                throw std::runtime_error("time.level_ratios[" + std::to_string(i) + "] must be an integer >= 2");
            p.time.level_ratios.push_back(size_t(o.as<double>()));
        }
        if (p.time.level_ratios.empty())
            throw std::runtime_error("time.level_ratios must not be empty");
    }
    sol::object overlaps = time["overlaps"];
    if (overlaps.valid()) {
        if (overlaps.get_type() != sol::type::table)
            throw std::runtime_error("time.overlaps must be a table of integers {o_1, o_2, ...}");
        p.time.overlaps.clear();
        sol::table t = overlaps.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            if (o.get_type() != sol::type::number || o.as<double>() != std::floor(o.as<double>()) || o.as<double>() < 0.0)
                throw std::runtime_error("time.overlaps[" + std::to_string(i) + "] must be an integer >= 0");
            p.time.overlaps.push_back(size_t(o.as<double>()));
        }
        if (p.time.overlaps.empty())
            throw std::runtime_error("time.overlaps must not be empty");
    }
    p.time.mlts_reuse = get_bool(time, "time", "mlts_reuse", p.time.mlts_reuse);
    p.time.mlts_refresh  = get_string(time, "time", "mlts_refresh", p.time.mlts_refresh, {"window", "substep"});
    p.time.li_system = get_string(time, "time", "li_system", p.time.li_system, {"dynamic", "constraint"});
    p.time.li_coupling = get_string(time, "time", "li_coupling", p.time.li_coupling, {"taylor", "stage"});
    p.time.li_condensation = get_bool(time, "time", "li_condensation", p.time.li_condensation);

    sol::table mesh = get_table(lua, "mesh");
    check_keys(mesh, "mesh", {"type", "x", "y", "nx", "ny", "level", "refine_points", "refine_levels", "refine_box", "files"});
    p.mesh.type          = get_string(mesh, "mesh", "type", p.mesh.type, {"cartesian", "fvca5", "polymesher", "polygons"});
    p.mesh.x             = get_pair(mesh, "mesh", "x", p.mesh.x);
    p.mesh.y             = get_pair(mesh, "mesh", "y", p.mesh.y);
    p.mesh.nx            = get_count(mesh, "mesh", "nx", p.mesh.nx);
    p.mesh.ny            = get_count(mesh, "mesh", "ny", p.mesh.ny);
    p.mesh.level         = get_count(mesh, "mesh", "level", p.mesh.level);
    p.mesh.refine_levels = get_count(mesh, "mesh", "refine_levels", p.mesh.refine_levels);
    sol::object box = mesh["refine_box"];
    if (box.valid()) {
        sol::table t = box.get_type() == sol::type::table ? box.as<sol::table>() : lua.create_table();
        p.mesh.refine_box.clear();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            if (o.get_type() != sol::type::number)
                throw std::runtime_error("mesh.refine_box must be {x0, x1, y0, y1}");
            p.mesh.refine_box.push_back(o.as<double>());
        }
        if (p.mesh.refine_box.size() != 4)
            throw std::runtime_error("mesh.refine_box must be {x0, x1, y0, y1}");
    }
    sol::object pts = mesh["refine_points"];
    if (pts.valid()) {
        if (pts.get_type() != sol::type::table)
            throw std::runtime_error("mesh.refine_points must be a table of points {{x, y}, ...}");
        p.mesh.refine_points.clear();
        sol::table t = pts.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            sol::table pt = o.get_type() == sol::type::table ? o.as<sol::table>() : lua.create_table();
            sol::object a = pt[1], b = pt[2];
            if (!a.valid() || !b.valid() || a.get_type() != sol::type::number || b.get_type() != sol::type::number)
                throw std::runtime_error("mesh.refine_points[" + std::to_string(i) + "] must be a point {x, y}");
            p.mesh.refine_points.push_back({a.as<double>(), b.as<double>()});
        }
    }
    sol::object files = mesh["files"];
    if (files.valid()) {
        if (files.get_type() != sol::type::table)
            throw std::runtime_error("mesh.files must be a table of file names {\"a\", \"b\", ...}");
        sol::table t = files.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            if (o.get_type() != sol::type::string)
                throw std::runtime_error("mesh.files[" + std::to_string(i) + "] must be a file name");
            p.mesh.files.push_back(o.as<std::string>());
        }
    }

    sol::table problem = get_table(lua, "problem");
    check_keys(problem, "problem", {"medium", "solution", "dirichlet", "source", "interface", "acoustic", "elastic", "zones",
                                    "omega", "theta", "ricker", "lshape_corner", "lshape_exponent"});
    p.problem.medium = get_string(problem, "problem", "medium", p.problem.medium, {"coupled", "acoustic", "elastic"});
    std::set<std::string> solutions;
    for (auto& s : wave_parameters::solution_names())
        solutions.insert(s[0]);
    p.problem.solution = get_string(problem, "problem", "solution", p.problem.solution, solutions);
    p.problem.dirichlet = get_string(problem, "problem", "dirichlet", p.problem.dirichlet, {"exact", "zero"});
    p.problem.source    = get_string(problem, "problem", "source", p.problem.source, {"solution", "zero"});
    p.problem.interface = get_string(problem, "problem", "interface", p.problem.interface, {"x", "y", "tags"});
    p.problem.omega     = get_number(problem, "problem", "omega", p.problem.omega);
    p.problem.theta     = get_number(problem, "problem", "theta", p.problem.theta);
    auto material_table = [&](const char* key) -> sol::table {
        sol::object o = problem[key];
        if (!o.valid()) return lua.create_table();
        if (o.get_type() != sol::type::table)
            throw std::runtime_error(std::string("problem.") + key + " must be a table {rho = .., vp = .., ...}");
        return o.as<sol::table>();
    };
    sol::table ac = material_table("acoustic"), el = material_table("elastic");
    check_keys(ac, "problem.acoustic", {"rho", "vp"});
    check_keys(el, "problem.elastic", {"rho", "vp", "vs"});
    p.problem.rho_acoustic = get_number(ac, "problem.acoustic", "rho", p.problem.rho_acoustic);
    p.problem.vp_acoustic  = get_number(ac, "problem.acoustic", "vp", p.problem.vp_acoustic);
    p.problem.rho_elastic  = get_number(el, "problem.elastic", "rho", p.problem.rho_elastic);
    p.problem.vp_elastic   = get_number(el, "problem.elastic", "vp", p.problem.vp_elastic);
    p.problem.vs_elastic   = get_number(el, "problem.elastic", "vs", p.problem.vs_elastic);
    sol::object zones = problem["zones"];
    if (zones.valid()) {
        if (zones.get_type() != sol::type::table)
            throw std::runtime_error("problem.zones must be a table {[tag] = {medium = .., rho = .., vp = .., vs = ..}, ...}");
        for (const auto& kv : zones.as<sol::table>()) {
            if (kv.first.get_type() != sol::type::number || kv.second.get_type() != sol::type::table)
                throw std::runtime_error("problem.zones: entries [tag] = {medium = .., rho = .., vp = .., vs = ..}");
            int tag = kv.first.as<int>();
            sol::table z = kv.second.as<sol::table>();
            std::string zn = "problem.zones[" + std::to_string(tag) + "]";
            check_keys(z, zn, {"medium", "rho", "vp", "vs"});
            wave_parameters::problem_data::zone zone;
            zone.medium = get_string(z, zn, "medium", "", {"acoustic", "elastic"});
            if (zone.medium.empty())
                throw std::runtime_error(zn + ".medium must be \"acoustic\" or \"elastic\"");
            zone.rho = get_number(z, zn, "rho", zone.rho);
            zone.vp  = get_number(z, zn, "vp", zone.vp);
            zone.vs  = get_number(z, zn, "vs", zone.vs);
            p.problem.zones[tag] = zone;
        }
    }
    p.problem.lshape_corner = get_pair(problem, "problem", "lshape_corner", p.problem.lshape_corner);
    p.problem.lshape_exponent = get_number(problem, "problem", "lshape_exponent", p.problem.lshape_exponent);
    {
        const double j = 1.5 * p.problem.lshape_exponent;
        if (p.problem.lshape_exponent <= 0.0 || std::abs(j - std::round(j)) > 1e-12)
            throw std::runtime_error("problem.lshape_exponent must be 2 j / 3, j = 1, 2, ...");
    }
    sol::object rk = problem["ricker"];
    if (rk.valid()) {
        if (rk.get_type() != sol::type::table)
            throw std::runtime_error("problem.ricker must be a table {field = .., center = {x, y}, frequency = .., ...}");
        sol::table r = rk.as<sol::table>();
        check_keys(r, "problem.ricker", {"field", "center", "frequency", "amplitude", "velocity", "shape", "radius", "axis",
                                         "aperture", "cycles", "density"});
        p.problem.ricker.field = get_string(r, "problem.ricker", "field", p.problem.ricker.field,
                                            {"acoustic_flux", "acoustic_pressure", "elastic_velocity"});
        p.problem.ricker.center    = get_pair(r, "problem.ricker", "center", p.problem.ricker.center);
        p.problem.ricker.frequency = get_number(r, "problem.ricker", "frequency", p.problem.ricker.frequency);
        p.problem.ricker.amplitude = get_number(r, "problem.ricker", "amplitude", p.problem.ricker.amplitude);
        p.problem.ricker.velocity  = get_number(r, "problem.ricker", "velocity", p.problem.ricker.velocity);
        p.problem.ricker.shape     = get_string(r, "problem.ricker", "shape", p.problem.ricker.shape, {"gaussian", "focused"});
        p.problem.ricker.radius    = get_number(r, "problem.ricker", "radius", p.problem.ricker.radius);
        p.problem.ricker.axis      = get_number(r, "problem.ricker", "axis", p.problem.ricker.axis);
        p.problem.ricker.aperture  = get_number(r, "problem.ricker", "aperture", p.problem.ricker.aperture);
        p.problem.ricker.cycles    = get_number(r, "problem.ricker", "cycles", p.problem.ricker.cycles);
        p.problem.ricker.density   = get_number(r, "problem.ricker", "density", p.problem.ricker.density);
        if (p.problem.ricker.shape == "focused" && p.problem.ricker.radius <= 0.0)
            throw std::runtime_error("problem.ricker.radius must be positive with shape = \"focused\"");
    }

    sol::table study = get_table(lua, "study");
    check_keys(study, "study", {"type", "last_level", "refine_dt", "time_levels", "reference", "reference_factor",
                                "cfl_threshold", "cfl_decrement", "cfl_runs", "cfl_absolute", "cfl_random"});
    p.study.type             = get_string(study, "study", "type", p.study.type, {"run", "space", "time", "cfl", "spectrum", "rho"});
    p.study.cfl_threshold    = get_number(study, "study", "cfl_threshold", p.study.cfl_threshold);
    p.study.cfl_decrement    = get_number(study, "study", "cfl_decrement", p.study.cfl_decrement);
    p.study.cfl_runs         = get_count(study, "study", "cfl_runs", p.study.cfl_runs);
    p.study.cfl_random       = get_bool(study, "study", "cfl_random", p.study.cfl_random);
    p.study.cfl_absolute     = get_bool(study, "study", "cfl_absolute", p.study.cfl_absolute);
    p.study.last_level       = get_count(study, "study", "last_level", p.mesh.level);
    p.study.refine_dt        = get_bool(study, "study", "refine_dt", p.study.refine_dt);
    p.study.time_levels      = get_count(study, "study", "time_levels", p.study.time_levels);
    p.study.reference        = get_string(study, "study", "reference", p.study.reference, {"numerical", "exact"});
    p.study.reference_factor = get_count(study, "study", "reference_factor", p.study.reference_factor);

    sol::table output = get_table(lua, "output");
    check_keys(output, "output", {"file", "table", "energy", "sensors", "sensor_every", "snapshots", "snapshot_samples", "cpu", "errors",
                                  "profile", "snapshot_elastic", "peak_samples"});
    p.output.file   = get_string(output, "output", "file", p.output.file, {});
    p.output.table  = get_string(output, "output", "table", p.output.table, {});
    p.output.energy = get_bool(output, "output", "energy", p.output.energy);
    p.output.sensor_every = get_count(output, "output", "sensor_every", p.output.sensor_every);
    p.output.cpu = get_bool(output, "output", "cpu", p.output.cpu);
    p.output.profile = get_bool(output, "output", "profile", p.output.profile);
    p.output.errors = get_string(output, "output", "errors", p.output.errors, {"fields", "energy"});
    sol::object sens = output["sensors"];
    if (sens.valid()) {
        if (sens.get_type() != sol::type::table)
            throw std::runtime_error("output.sensors must be a table {{name = .., x = .., y = .., fields = {..}}, ...}");
        sol::table t = sens.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            std::string sn = "output.sensors[" + std::to_string(i) + "]";
            if (o.get_type() != sol::type::table)
                throw std::runtime_error(sn + " must be a table {name = .., x = .., y = .., fields = {..}}");
            sol::table st = o.as<sol::table>();
            check_keys(st, sn, {"name", "x", "y", "fields"});
            wave_parameters::output_data::sensor sensor;
            sensor.name = get_string(st, sn, "name", "S" + std::to_string(i), {});
            sensor.point = {get_number(st, sn, "x", 0.0), get_number(st, sn, "y", 0.0)};
            sol::object fl = st["fields"];
            if (!fl.valid() || fl.get_type() != sol::type::table)
                throw std::runtime_error(sn + ".fields must be a table of field names");
            sol::table ft = fl.as<sol::table>();
            for (size_t j = 1; j <= ft.size(); j++) {
                sol::object fo = ft[j];
                if (fo.get_type() != sol::type::string || !wave_parameters::sensor_fields().count(fo.as<std::string>())) {
                    std::string list;
                    for (auto& v : wave_parameters::sensor_fields()) list += " \"" + v + "\"";
                    throw std::runtime_error(sn + ".fields[" + std::to_string(j) + "] must be one of" + list);
                }
                sensor.fields.push_back(fo.as<std::string>());
            }
            p.output.sensors.push_back(sensor);
        }
    }
    p.output.snapshot_samples = get_count(output, "output", "snapshot_samples", p.output.snapshot_samples);
    p.output.snapshot_elastic = get_string(output, "output", "snapshot_elastic", p.output.snapshot_elastic, {"velocity", "stress"});
    p.output.peak_samples = get_count(output, "output", "peak_samples", p.output.peak_samples);
    sol::object snaps = output["snapshots"];
    if (snaps.valid()) {
        if (snaps.get_type() != sol::type::table)
            throw std::runtime_error("output.snapshots must be a table of times {t1, t2, ...}");
        sol::table t = snaps.as<sol::table>();
        for (size_t i = 1; i <= t.size(); i++) {
            sol::object o = t[i];
            if (o.get_type() != sol::type::number)
                throw std::runtime_error("output.snapshots[" + std::to_string(i) + "] must be a time");
            p.output.snapshots.push_back(o.as<double>());
        }
    }

    // consistency
    if (p.time.scheme == "erk" && (p.stages() < 1 || p.stages() > 4))
        throw std::runtime_error("time.stages must be 1..4 for erk");
    if (p.time.scheme == "sdirk" && (p.stages() < 1 || p.stages() > 3))
        throw std::runtime_error("time.stages must be 1..3 for sdirk");
    if (p.time.level_ratio < 2)
        throw std::runtime_error("time.level_ratio must be at least 2");
    if (p.time.scheme == "newmark" && p.formulation != "second_order")
        throw std::runtime_error("time.scheme = \"newmark\" needs formulation = \"second_order\"");
    if (p.time.scheme == "newmark" && !(p.time.newmark_beta > 0.0 && p.time.newmark_gamma >= 0.5))
        throw std::runtime_error("newmark: newmark_beta > 0 (implicit scheme) and newmark_gamma >= 1/2 are required");
    if (p.time.final_time <= 0.0 || p.time.n_steps == 0 || p.time.dt_h < 0.0)
        throw std::runtime_error("time.final_time and time.n_steps must be positive, time.dt_h nonnegative");
    if (p.time.dt_h > 0.0 && p.study.type == "time")
        throw std::runtime_error("time.dt_h: the time study sets dt from time.n_steps");
    if (p.formulation == "second_order" && !p.least_squares() && !p.scaled_stabilization())
        throw std::runtime_error("space.stabilization = \"hho\" is always scaled by 1/h in diskpp: set space.scaling = \"1/h\"");
    if (p.formulation == "first_order" && p.space.weight_given && (p.space.weight_auto || p.space.weight <= 0.0))
        throw std::runtime_error("first order: space.weight must be a positive number");
    if (p.space.elastic_formulation == "three_fields" && (p.formulation != "first_order" || p.problem.medium != "elastic" ||
                                                         p.time.scheme != "sdirk" || !p.scaled_stabilization()))
        throw std::runtime_error("space.elastic_formulation = \"three_fields\": first order, problem.medium = \"elastic\", "
                                 "time.scheme = \"sdirk\" and space.scaling = \"1/h\" (assembler of Burman-Duran-Ern-Steins 2021)");
    if (p.problem.interface == "tags" && p.mesh.type != "polygons")
        throw std::runtime_error("problem.interface = \"tags\" needs mesh.type = \"polygons\" (meshes with cell tags)");
    if (p.problem.interface == "tags" && p.problem.zones.empty())
        throw std::runtime_error("problem.interface = \"tags\": problem.zones must give the medium and material of each tag");
    if (p.problem.interface != "x" && p.formulation != "first_order")
        throw std::runtime_error("problem.interface = \"" + p.problem.interface + "\": first order only");
    if (p.study.type == "spectrum" && p.formulation != "first_order")
        throw std::runtime_error("study.type = \"" + p.study.type + "\": first order only");
    if (p.study.type == "cfl" && !p.explicit_scheme())
        throw std::runtime_error("study.type = \"cfl\": an explicit time.scheme (erk, lts, li, mlts)");
    if (p.study.type == "cfl" && (p.study.cfl_threshold <= 0.0 || p.study.cfl_decrement <= 0.0 || p.study.cfl_decrement >= 1.0))
        throw std::runtime_error("study.cfl_threshold > 0 and 0 < study.cfl_decrement < 1 are required");
    if (p.solver.iterative && !(p.formulation == "first_order" && p.time.scheme == "sdirk"))
        throw std::runtime_error("solver.iterative: first order with time.scheme = \"sdirk\" only");
    if (p.solver.matrix_free && !(p.formulation == "first_order" && (p.time.scheme == "erk" || p.time.scheme == "lts"
                                                                     || p.time.scheme == "li" || p.time.scheme == "sdirk"
                                                                     || p.time.scheme == "mlts"))
                             && !(p.formulation == "second_order" && (p.time.scheme == "erk" || p.time.scheme == "lts"
                                                                      || p.time.scheme == "mlts")))
        throw std::runtime_error("solver.matrix_free: first order with time.scheme = \"erk\", \"lts\", \"li\", \"mlts\" or "
                                 "\"sdirk\", second order with \"erk\", \"lts\" or \"mlts\" only");
    if (p.solver.matrix_free_batch < 1)
        throw std::runtime_error("solver.matrix_free_batch must be positive");
    if (p.solver.tolerance <= 0.0 || p.solver.tolerance_growth <= 0.0)
        throw std::runtime_error("solver.tolerance and solver.tolerance_growth must be positive");
    if (p.time.steps_growth <= 0.0)
        throw std::runtime_error("time.steps_growth must be positive");
    if ((!p.output.sensors.empty() || !p.output.snapshots.empty()) && p.formulation != "first_order")
        throw std::runtime_error("output.sensors, output.snapshots: first order only");
    if (p.output.sensor_every == 0)
        throw std::runtime_error("output.sensor_every must be positive");
    if (p.output.errors == "energy" && (p.formulation != "first_order" || p.space.elastic_formulation == "three_fields"))
        throw std::runtime_error("output.errors = \"energy\": first order, velocity-stress formulation only");
    if (p.formulation == "first_order" && ((p.space.weight_acoustic >= 0.0 && p.space.weight_acoustic == 0.0) ||
                                           (p.space.weight_elastic >= 0.0 && p.space.weight_elastic == 0.0)))
        throw std::runtime_error("space.weight_acoustic, space.weight_elastic must be positive");
    if (p.output.energy && p.formulation == "second_order" && p.time.scheme != "newmark")
        throw std::runtime_error("output.energy: second order with time.scheme = \"newmark\" only");
    if (p.formulation == "first_order" && p.explicit_scheme() && !p.least_squares())
        throw std::runtime_error("the first-order explicit and locally implicit integrators eliminate the faces with a "
                                 "block-diagonal face block: space.stabilization must be \"least_squares\"");
    if (p.formulation == "second_order" && !p.space.weight_auto && p.space.weight <= 0.0)
        throw std::runtime_error("space.weight must be positive (or \"auto\")");
    if (p.faces.tolerance <= 0.0 || p.faces.max_iterations == 0)
        throw std::runtime_error("faces.tolerance and faces.max_iterations must be positive");
    for (auto& s : wave_parameters::solution_names())
        if (s[0] == p.problem.solution && !s[1].empty() && s[1] != p.problem.medium)
            throw std::runtime_error("problem.solution = \"" + s[0] + "\" needs problem.medium = \"" + s[1] + "\"");
    if (p.mesh.type == "cartesian") {
        if (!(p.mesh.x[0] < p.mesh.x[1]) || !(p.mesh.y[0] < p.mesh.y[1]) || p.mesh.nx == 0 || p.mesh.ny == 0)
            throw std::runtime_error("mesh: x[1] < x[2], y[1] < y[2] and nx, ny > 0 are required");
        if (p.problem.medium == "coupled" && p.problem.interface == "x") {
            if (!(p.mesh.x[0] < 0.0 && 0.0 < p.mesh.x[1]))
                throw std::runtime_error("mesh: the rectangle must contain the interface x = 0 (x[1] < 0 < x[2])");
            double ix = -p.mesh.x[0] / (p.mesh.x[1] - p.mesh.x[0]) * double(p.mesh.nx);
            if (std::abs(ix - std::round(ix)) > 1.0e-10)
                throw std::runtime_error("mesh: the interface x = 0 must be a line of the coarsest mesh (nx cells on [x[1], x[2]])");
        }
        if (p.problem.medium == "coupled" && p.problem.interface == "y") {
            if (!(p.mesh.y[0] < 0.0 && 0.0 < p.mesh.y[1]))
                throw std::runtime_error("mesh: the rectangle must contain the interface y = 0 (y[1] < 0 < y[2])");
            double iy = -p.mesh.y[0] / (p.mesh.y[1] - p.mesh.y[0]) * double(p.mesh.ny);
            if (std::abs(iy - std::round(iy)) > 1.0e-10)
                throw std::runtime_error("mesh: the interface y = 0 must be a line of the coarsest mesh (ny cells on [y[1], y[2]])");
        }
    }
    else {
        size_t last = p.study.type == "space" || p.study.type == "spectrum" ? std::max(p.study.last_level, p.mesh.level) : p.mesh.level;
        if (p.mesh.level < 1 || last > p.mesh.files.size())
            throw std::runtime_error("mesh.type = \"" + p.mesh.type + "\": mesh.files[l] is read for the levels l = mesh.level .. (last level), "
                                     "which must be in 1.." + std::to_string(p.mesh.files.size()));
        if (p.mesh.refine_levels > 0)
            throw std::runtime_error("mesh.refine_levels: local refinement is only available on cartesian meshes");
    }
    if ((p.study.type == "space" || p.study.type == "spectrum") && p.study.last_level < p.mesh.level)
        throw std::runtime_error("study.last_level must be >= mesh.level");
    if (p.study.type == "time" && (p.study.time_levels < 2 || p.study.reference_factor < 2))
        throw std::runtime_error("study.time_levels and study.reference_factor must be >= 2");
    return p;
}

#endif /* lua_parameters_hpp */
