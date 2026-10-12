//  Created by Romain Mottier
//
// GENUINE classical two-level LTS -- NOT the general multilevel Pell-sparse
// recursion used by ERK4_MLTS_CenterSquare_TwoLevel_conv_test (which forces
// L=2 out of the recursive band machinery shared with the deep-L drivers).
// two_level_square_sweep() uses
// elastoacoustic_four_fields_assembler::assemble_P(msh, h_c, nb_layer) to
// build the diagonal Pcoarse/Pfine projectors (threshold cell diameter vs
// h_c, with a vertex-sharing "layer" dilation around the fine region), then
// per macro step ONE erk_weight_LTS_coarse call (frozen coarse Taylor
// polynomial) followed by a fine-region integrator:
//
//   - classical_fine (this file): `p` explicit RK4 sub-steps
//     erk_weight_LTS_fine, no band-index bookkeeping, no recursion -- the
//     methodology of ERK4_LTS_conv_test.hpp / ERK4_LTS_Lshape_conv_test.hpp
//     (the latter confirms this machinery works unmodified on a
//     pure-acoustic domain, empty e_material, nb_layer=3 recommended for a
//     first-order mixed system as p grows);
//   - lisdirk_fine (LTS_ERK_SDIRK/LI_SDIRK_CenterSquare_conv_tests.hpp):
//     one implicit SDIRK(3,4) step over the whole macro step.
//
// All sweeps run on center_square_fixed meshes; the functions at the end of
// the file are the individual sweeps (mesh family, exact solution and
// output file of each one are unchanged).
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0
#pragma once

#include "../../lts_common.hpp"

struct two_level_square_config {
    std::string title;          // sweep banner
    std::string mesh_prefix;    // center_square_fixed/meshes/<prefix>_k<k>_N<N>.txt
    std::string level_tag;      // appended to the "MESH LEVEL N=" banner
    std::string conv_stem;      // lshape/results/<stem>_k_<k>.txt
    std::string p_global_note;  // appended to the n_cells line
    bool L_column      = true;  // L_levels column (always 2)
    bool cpu_time      = true;  // cpu_time_s column (getrusage)
    bool nt_dt         = false; // tf overridable by MLTS_TF, "nt dt" columns
    bool coarse_always = false; // coarse precompute even when p_global == 1
    // MLTS_DEBUG_IC / MLTS_DEBUG_OP / MLTS_DEBUG_PTS diagnostics.
    bool debug_ic = false, debug_op = false, debug_pts = false;
};

// make_fine(pb, eval_F, dt_macro) is called once per mesh level and returns
// the fine-region step fine(x_dof_n, w, tn) applied after the coarse
// precompute of each macro step.
template<typename Solution, typename MakeFine>
void two_level_square_sweep(int argc, char **argv, const two_level_square_config& cfg, const Solution& sol0, MakeFine&& make_fine)
{
    using namespace lts;
    std::cout << std::endl << bold << red << cfg.title << std::endl << std::endl;

    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);
    const size_t mesh_k = sim_data.m_k_degree;
    std::vector<int> level_N = levels(0, 4);

    // nb_layer: vertex-sharing dilation passes around the fine region in
    // assemble_P -- ERK4_LTS_Lshape_conv_test.hpp's comment explains a
    // first-order mixed acoustic system needs more than the default 1 as p
    // grows; the haloscale family reaches p_global=1024, so default to 3
    // for safety (overridable).
    const int nb_layer = env_int("MLTS_NB_LAYER", 3);

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);
    std::ostringstream conv_fname;
    conv_fname << out_dir << "/" << cfg.conv_stem << "_k_" << mesh_k << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# N  h_max  h_min  p_global" << (cfg.L_column ? "  L_levels" : "")
             << "  L2_error_pressure  wall_time_s  L2_error_velocity"
             << (cfg.cpu_time ? "  cpu_time_s" : "") << (cfg.nt_dt ? "  nt  dt" : "") << "\n";

    for (int N : level_N) {

        std::cout << bold << red << "\n   ================ MESH LEVEL N=" << N << cfg.level_tag << " ================" << reset << std::endl;

        std::ostringstream mesh_fname;
        mesh_fname << mesh_root << "center_square_fixed/meshes/" << cfg.mesh_prefix << "_k" << mesh_k << "_N" << N << ".txt";
        acoustic_problem<Solution> pb(mesh_fname.str(), sol0);
        mesh_type& msh = pb.msh;
        const RealType h_max = pb.h_max, h_min = pb.h_min;
        const int p_global = pb.p_global;
        RealType h_c = (p_global == 1) ? 1.25 * h_max : 0.75 * h_max;

        std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p_global=" << p_global << cfg.p_global_note << reset << std::endl;

        pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
        assembler_type& assembler = pb.assembler();
        erk_type& erk_an = pb.erk();
        vector_type& x_dof = pb.x_dof;

        if (cfg.debug_ic && std::getenv("MLTS_DEBUG_IC")) {
            // L2 error of the CELL pressure dof against the exact field at
            // t=0 immediately after IC projection + refresh_faces_unknowns,
            // with ZERO time steps taken -- isolates whether the flat,
            // h-independent residual seen at small tf already exists before
            // any time integration at all.
            std::cout << bold << green << "      [DEBUG IC] pressure L2 error vs s_v_fun0 at t=0, BEFORE any time step = "
                      << std::setprecision(10) << pressure_l2_error(msh, hho_di, x_dof, pb.sol.p0()) << reset << std::endl;
        }

        if constexpr (std::is_same_v<Solution, analytic_solution>) {
            if (cfg.debug_op && std::getenv("MLTS_DEBUG_OP")) {
                // Apply the ASSEMBLED discrete operator once (k =
                // erk_weight(x_dof), i.e. Mc_inv*(Fc - Kcc*y_c - Kcf*y_f)
                // then face back-substitution -- the exact SAME B operator
                // used throughout erk_weight_LTS_coarse/fine and
                // dirk_weight_LTS_fine) to the EXACT-projected IC at t=0, and
                // compare the CELL PRESSURE part of k (the discrete dp_h/dt
                // implied by the assembled system) against the TRUE dp/dt
                // (central finite difference in TIME, independent of the
                // mesh/scheme entirely) at a few cell barycenters. Isolates
                // whether the ASSEMBLED PDE OPERATOR itself is consistent
                // with the solution's own time derivative, with NO BC, no
                // time-stepping and no P-splitting involved at all.
                erk_an.ZeroFc();
                vector_type k_op;
                erk_an.erk_weight(x_dof, k_op);

                // Evaluate_s_v(double&t) returns a lambda that captures t BY
                // REFERENCE -- each snapshot needs its OWN persistent time
                // variable, evaluated immediately, not a shared/reused one.
                RealType hfd = 1e-6;
                RealType t_plus  = hfd;
                RealType t_minus = -hfd;
                auto s_v_fun_plus  = pb.sol.functions().Evaluate_s_v(t_plus);
                auto s_v_fun_minus = pb.sol.functions().Evaluate_s_v(t_minus);

                std::cout << bold << green << "      [DEBUG OP] discrete k=B*y0 (cell pressure part) vs TRUE dp/dt (central FD, h=1e-6), first 5 cells:" << reset << std::endl;
                size_t ci = 0;
                for (auto & cell : msh) {
                    if (ci >= 5) break;
                    auto bar = barycenter(msh, cell);
                    vector_type kp_dof = k_op.block(ci * pb.cell_dof + pb.n_vec_dof, 0, pb.n_scal_dof, 1);
                    auto cell_basis_dbg = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                    auto phi = cell_basis_dbg.eval_functions(bar);
                    RealType kp_h = 0.0;
                    for (int i = 0; i < phi.size(); i++) kp_h += kp_dof(i) * phi(i);
                    RealType dpdt_true = (s_v_fun_plus(bar) - s_v_fun_minus(bar)) / (2.0*hfd);
                    std::cout << bold << green << "         cell " << ci << " bar=(" << bar.x() << "," << bar.y()
                              << ")  k_h(dp/dt)=" << kp_h << "  dp/dt_true=" << dpdt_true
                              << "  ratio=" << (kp_h/dpdt_true) << reset << std::endl;
                    ci++;
                }
            }
        }

        // Classical two-level P projectors: threshold cell diameter vs h_c,
        // NOT the multilevel Pell/band index bookkeeping.
        assembler.assemble_P(msh, h_c, nb_layer);

        // ---------------- Time discretization ----------------
        // MLTS_TF (nt_dt drivers): the default tf=0.1 already gives clean
        // errors converging at the expected rate on the sinusoidal-in-time
        // solutions (not initial-projection noise); left overridable to
        // re-check at a larger fraction of the period.
        RealType cfl_factor = env_real("MLTS_CFL_FACTOR", 0.05);
        RealType ti = 0.0, tf = 0.1;
        if (cfg.nt_dt) tf = env_real("MLTS_TF", tf);
        RealType dt_macro = cfl_factor * h_max;
        size_t nt = macro_steps(tf, dt_macro);
        if (cfg.nt_dt)
            std::cout << bold << cyan << "      nt=" << nt << "  dt=" << dt_macro << "  tf=" << tf << reset << std::endl;
        else
            std::cout << bold << cyan << "      nt=" << nt << reset << std::endl;

        auto eval_F = [&](RealType t_abs) -> vector_type { return pb.rhs(t_abs); };
        auto fine = make_fine(pb, eval_F, dt_macro);
        const size_t n_dof = x_dof.rows();

        // ---------------- Time marching ----------------
        timecounter tc_run;
        tc_run.tic();
        RealType cpu_start = cpu_seconds();
        for (size_t it = 1; it <= nt; it++) {
            RealType tn = ti + dt_macro * (it - 1);
            auto x_dof_n = x_dof;

            std::vector<vector_type> w(4);
            for (int j = 0; j < 4; ++j) { w[j].resize(n_dof); w[j].setZero(); }

            erk_an.ZeroFc();
            if (cfg.coarse_always || p_global != 1) {
                RealType tn12 = tn + 0.5*dt_macro;
                RealType tn1  = tn + dt_macro;
                vector_type Fn   = eval_F(tn);
                vector_type Fn12 = eval_F(tn12);
                vector_type Fn1  = eval_F(tn1);
                erk_an.ZeroFc();
                erk_an.erk_weight_LTS_coarse(x_dof_n, assembler.Pcoarse, w, Fn, Fn12, Fn1, dt_macro);
            }

            fine(x_dof_n, w, tn);

            x_dof = x_dof_n;
            if (it % std::max<size_t>(1, nt/10) == 0)
                std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
        }
        tc_run.toc();
        RealType cpu_time = cpu_seconds() - cpu_start;

        // ---------------- Error computation (pressure + velocity) ----------------
        if (cfg.debug_pts && std::getenv("MLTS_DEBUG_PTS")) {
            auto s_v_fun_T = pb.sol.p_exact(tf);
            RealType bx_min=1e10,bx_max=-1e10,by_min=1e10,by_max=-1e10;
            for (auto & cell : msh) {
                auto bar = barycenter(msh, cell);
                bx_min = std::min(bx_min, (RealType)bar.x()); bx_max = std::max(bx_max, (RealType)bar.x());
                by_min = std::min(by_min, (RealType)bar.y()); by_max = std::max(by_max, (RealType)bar.y());
            }
            std::cout << bold << green << "      [DEBUG PTS] domain barycenter bbox: x in [" << bx_min << "," << bx_max
                      << "]  y in [" << by_min << "," << by_max << "]" << reset << std::endl;
            size_t dbg_i = 0;
            for (auto & cell : msh) {
                if (dbg_i >= 5) break;
                auto bar = barycenter(msh, cell);
                vector_type p_dof_dbg = x_dof.block(dbg_i * pb.cell_dof + pb.n_vec_dof, 0, pb.n_scal_dof, 1);
                auto cell_basis_dbg = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                auto phi = cell_basis_dbg.eval_functions(bar);
                RealType p_h = 0.0;
                for (int i = 0; i < phi.size(); i++) p_h += p_dof_dbg(i) * phi(i);
                RealType p_exact = s_v_fun_T(bar);
                std::cout << bold << green << "         cell " << dbg_i << " bar=(" << bar.x() << "," << bar.y()
                          << ")  p_h=" << p_h << "  p_exact=" << p_exact
                          << "  diff=" << (p_h - p_exact) << reset << std::endl;
                dbg_i++;
            }
        }

        RealType l2_error   = pb.pressure_error(x_dof, tf);
        RealType v_l2_error = pb.velocity_error(x_dof, tf);

        std::cout << bold << red << "      wall time: " << tc_run << " s   ";
        if (cfg.cpu_time) std::cout << "cpu time: " << cpu_time << " s   ";
        std::cout << "L2_error(pressure) = " << std::setprecision(10) << l2_error
                  << "   L2_error(velocity) = " << v_l2_error << reset << std::endl;

        // Columns: N h_max h_min p_global [L_levels] L2_error_pressure
        // wall_time_s L2_error_velocity [cpu_time_s] [nt dt] -- cpu_time_s is
        // the exact process CPU time (user+system, via getrusage), immune to
        // wall-clock contention from unrelated processes.
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << p_global;
        if (cfg.L_column) conv_log << " " << 2;
        conv_log << " " << l2_error << " " << tc_run.elapsed() << " " << v_l2_error;
        if (cfg.cpu_time) conv_log << " " << cpu_time;
        if (cfg.nt_dt) conv_log << " " << nt << " " << dt_macro;
        conv_log << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}

// Fine region: p_global explicit RK4 sub-steps of size dt_macro/p_global.
inline auto classical_fine = [](auto& pb, auto& eval_F, lts::RealType dt_macro) {
    lts::RealType dtau = dt_macro / pb.p_global;
    return [&pb, &eval_F, dtau](lts::vector_type& x_dof_n, std::vector<lts::vector_type>& w, lts::RealType tn) {
        for (int m = 0; m < pb.p_global; m++) {
            lts::RealType tm  =  m      * dtau;
            lts::RealType tmh = (m+0.5) * dtau;
            lts::RealType tm1 = (m+1.0) * dtau;
            lts::vector_type Fm  = eval_F(tn + tm);
            lts::vector_type Fmh = eval_F(tn + tmh);
            lts::vector_type Fm1 = eval_F(tn + tm1);
            pb.erk().erk_weight_LTS_fine(x_dof_n, pb.assembler().Pfine, w, Fm, Fmh, Fm1, tm, dtau);
        }
    };
};

static two_level_square_config classical_haloscale_config(const std::string& title, const std::string& tag, const std::string& stem) {
    two_level_square_config cfg;
    cfg.title       = title;
    cfg.mesh_prefix = "centersquarehaloscale_graded";
    cfg.level_tag   = tag;
    cfg.conv_stem   = stem;
    return cfg;
}

// Fixed-p mesh families centersq_pfix{pfam} (MLTS_PFAM, default 1; p_global
// = 2^pfam, constant across the whole N sweep), same square/sinusoidal
// problem and tf/cfl_factor as ERK4_MLTS_CenterSquare_TwoLevel_conv_test, so
// the two drivers are directly comparable -- one "classical", one "MLTS
// machinery at L=2".
void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_conv_test(int argc, char **argv) {
    const int pfam = lts::env_int("MLTS_PFAM", 1);
    two_level_square_config cfg;
    cfg.title       = "   CLASSICAL TWO-LEVEL LTS-RK4 SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP";
    cfg.mesh_prefix = "centersq_pfix" + std::to_string(pfam);
    cfg.level_tag   = " (classical p=2^" + std::to_string(pfam) + ")";
    cfg.conv_stem   = "centersquare_classicaltwolevel_p" + std::to_string(1 << pfam) + "_convergence";
    cfg.cpu_time    = false;
    two_level_square_sweep(argc, argv, cfg,
        lts::analytic_solution(scal_vec_analytic_functions::EFunctionType::reproduction_acoustic), classical_fine);
}

// Classical-vs-multilevel EXACT comparison: reads the SAME mesh family
// (centersquarehaloscale_graded, generate_haloscale_family.py) as the
// multilevel ERK4_MLTS_CenterSquare_Adaptive_conv_test -- both drivers run
// on IDENTICAL meshes (same h_max, same p_global per N), making wall/CPU
// time and accuracy directly comparable, not just order-of-magnitude.
void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_conv_test(int argc, char **argv) {
    two_level_square_config cfg = classical_haloscale_config(
        "   CLASSICAL TWO-LEVEL LTS-RK4 SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP",
        " (classical two-level, haloscale mesh)", "centersquare_classicaltwolevel_haloscale_convergence");
    two_level_square_sweep(argc, argv, cfg,
        lts::analytic_solution(scal_vec_analytic_functions::EFunctionType::reproduction_acoustic), classical_fine);
}

// Same as the HaloScale sweep with ONLY the analytical solution changed:
// reproduction_acoustic (t^2*sin(pi x)*sin(pi y), sinusoidal in SPACE but
// only QUADRATIC in TIME) is replaced by a solution sinusoidal in x, y AND t
// (EFunctionNonPolynomial; first written for EFunctionPlaneWaveAcoustic,
// theta=x+y-sqrt(2)*t, hence the name). Regenerates the classical two-level
// LTS-RK4 curve for a fair comparison against LI-SDIRK(4,3) on the
// identical solution (LI_SDIRK_CenterSquare_HaloScale_PlaneWave_conv_test).
void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_PlaneWave_conv_test(int argc, char **argv) {
    two_level_square_config cfg = classical_haloscale_config(
        "   CLASSICAL TWO-LEVEL LTS-RK4 SQUARE PLANE-WAVE (sin in x,y,t) CONVERGENCE SWEEP",
        " (classical two-level, haloscale mesh, plane wave)", "centersquare_classicaltwolevel_haloscale_planewave_convergence");
    cfg.nt_dt = cfg.debug_ic = cfg.debug_op = cfg.debug_pts = true;
    two_level_square_sweep(argc, argv, cfg,
        lts::analytic_solution(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial), classical_fine);
}

// Same with the slow sine (lts::slow_sine_solution, omega = pi/2, ~4x slower
// in time than EFunctionNonPolynomial's sqrt(2)*pi).
void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_SlowSine_conv_test(int argc, char **argv) {
    two_level_square_config cfg = classical_haloscale_config(
        "   CLASSICAL TWO-LEVEL LTS-RK4 SQUARE SLOW-SINE (sin in x,y,t, omega=pi/2) CONVERGENCE SWEEP",
        " (classical two-level, haloscale mesh, slow sine)", "centersquare_classicaltwolevel_haloscale_slowsine_convergence");
    cfg.nt_dt = cfg.debug_ic = cfg.debug_pts = true;
    two_level_square_sweep(argc, argv, cfg, lts::slow_sine_solution(), classical_fine);
}
