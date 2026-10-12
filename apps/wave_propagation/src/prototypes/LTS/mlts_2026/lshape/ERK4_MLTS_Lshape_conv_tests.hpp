//  Created by Romain Mottier
//
// Multi-level LTS-RK4 convergence sweeps for the L-shape re-entrant corner
// benchmark (lts::lshape_corner_solution) on the graded meshes of
// mesh_generation/lshape_graded (lshape_graded_k<k>_N<N>.txt, N=0..4 by
// default, MLTS_N_MAX overrides), one function per integrator of
// ../mlts_steppers.hpp.
//
// NOTE on cost: multi-level reduces the PER-STEP cost (each step only
// touches a band's own+halo dofs, not the whole mesh), but NOT the total
// NUMBER of finest-scale steps, which is fixed by p_global=h_max/h_min --
// a property of the mesh, not of how many levels we use. For this k=3
// grading, p_global grows ~32x per level N (1, 32, 1024, 32768, ~1e6 for
// N=0..4), so total cost (~nt*p_global) still grows steeply: N=3 takes
// hours, N=4 is likely impractical.
//
// CFL: the default 0.05 was recalibrated from 0.002 by a binary search on
// the true RK4 stability limit (see report, sec:cflfix), found near
// 0.08-0.085; 0.05 keeps a safety margin while being ~25x larger than the
// old default, with L2 error identical to 6 significant digits (spatial
// error dominates at this order). HybridTerminal still uses 0.002.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0
#pragma once

#include "../mlts_steppers.hpp"

static lts::mlts_sweep_config mlts_lshape_config(const std::string& stem, bool read_force_l) {
    lts::mlts_sweep_config cfg;
    cfg.title       = "   MULTI-LEVEL LTS L-SHAPE CONVERGENCE SWEEP (N=0..4)";
    cfg.mesh_family = "lshape_graded";
    cfg.mesh_prefix = "lshape_graded";
    cfg.levels      = lts::level_labels(lts::levels(0, 4));
    cfg.conv_stem   = stem;
    cfg.tf          = 0.05;
    if (read_force_l)
        cfg.thresholds = [](lts::RealType h_max, lts::RealType h_min, int) { return lts::env_thresholds(h_max, h_min); };
    else
        cfg.thresholds = [](lts::RealType h_max, lts::RealType h_min, int) { return lts::natural_thresholds(h_max, h_min); };
    return cfg;
}

static auto mlts_global_exact = [](auto& pb, const lts::bands& bd) { return mlts::make_global_exact(pb, bd); };

// Proven-exact reference: GLOBAL-matrix Diaz-Grote recursion
// (_coarse_global), O(n_dof) per call, no truncation -- validated to
// machine precision against the trusted 2-level reference at L=2..5.
void ERK4_MLTS_Lshape_GlobalExact_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence", true),
        lts::lshape_corner_solution(), mlts_global_exact);
}

// Restricted recursion on own + 8-ring halo sub-blocks (build_LTS_subblock):
// ~5-9x residual vs GlobalExact, traced to the fixed, guessed halo width
// (8 rings; 20 rings tested WORSE) not always being a large-enough
// neighbourhood for the B-chain (up to B^3) to reproduce the global
// operator at a band's own dofs. (Validation: error ratio vs the trusted
// unrestricted reference 1.22x at L=3/N=1, 9.1x at L=5/N=2, see
// ERK4_MLTS_timing_test / ERK4_MLTS_L3_v2design_validation_test.)
// Writes lshape_mlts_convergence_restrictedexact_k_<k>.txt (it used to
// share lshape_mlts_convergence_k_<k>.txt with GlobalExact, HybridTerminal
// and DiagN2, each run overwriting the others' results).
void ERK4_MLTS_Lshape_RestrictedExact_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence_restrictedexact", false),
        lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) { return mlts::make_restricted_exact(pb, bd, 8); });
}

// EXPERIMENTAL "P_ell-faithful masking" hypothesis (refuted: no better).
// _coarse_v2 computes each band's own-dynamics B-chain from an UNMASKED
// halo (finer AND coarser neighbours) -- unlike Almquist-Mehlin eq.
// (20)-(22), where B_L = B*P_L masks OUT every coarser (ancestor) level at
// EVERY application of B. The unmasked chain re-reads its parent's CURRENT
// value, double-counting a contribution already carried by the
// Taylor-shifted term -- a plausible cause of the ~5.9x residual at L=5
// (N=2), of the NaN at L=8 (N=3), and of widening the halo (8->20 rings)
// making things WORSE. build_LTS_subblock_leveled + _coarse_v3 split the
// halo into finer/coarser and mask the B-chain to "own+finer" (P_ell) at
// every step, keeping the one-time final output leak into the coarser halo
// (the correct (P_ell-P_{ell+1}) term). Masked and unmasked chains turned
// out to give essentially IDENTICAL error once overlap is added to both.
// Run with halo_rings = 8 (the 4 rings first planned were not used).
void ERK4_MLTS_Lshape_RestrictedLeveled_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence_leveled", false),
        lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) { return mlts::make_restricted_leveled(pb, bd, 8); });
}

// LITERAL reproduction of Almquist-Mehlin's Algorithm 3 (eq. 20-22),
// _coarse_mehlin instead of _coarse_global: masks the B-chain to P_ell =
// "level ell AND every finer level" at EVERY application of B (not just
// once at the end, the way the reformulated _coarse_global does); the
// ancestor cross-term sum of eq. (22) is exactly reproduced by the
// incremental Taylor-shift accumulation. PURPOSE: pure fidelity
// validation, not performance -- if GlobalExact is a valid algebraic
// simplification of the paper's literal construction, this reproduces the
// EXACT SAME numbers (it matches GlobalExact to <0.5%). GLOBAL matrices
// throughout: P_ell masking to "all finer levels" cannot be represented by
// any FIXED-width local halo.
void ERK4_MLTS_Lshape_MehlinExact_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence_mehlin", true),
        lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) { return mlts::make_mehlin_exact(pb, bd); });
}

// RestrictedExact's recursion with the per-band sub-matrices built from
// the exact P_ell domain + overlap (build_LTS_subblock_Pell) instead of a
// guessed ring count -- the fix suggested by GlobalExact/MehlinExact
// matching the reference with FULL matrices (only the sub-matrix REACH was
// wrong). Still ~5x residual: the truncated sub-matrix itself drops
// genuine Kcc/Kcf/Kfc coupling.
void ERK4_MLTS_Lshape_MehlinRestricted_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence_mehlinrestricted", true),
        lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) { return mlts::make_mehlin_restricted(pb, bd); });
}

// MehlinExact's math with a sparse-matvec implementation
// (_coarse_Pell_sparse): exact and faster.
void ERK4_MLTS_Lshape_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_lshape_config("lshape_mlts_convergence_pellsparse", true),
        lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) { return mlts::make_pell_sparse(pb, bd); });
}

// Global non-terminal levels + terminal restricted to a 3-ring sub-block,
// ancestor own dofs outside its coverage advanced in closed form (see
// make_hybrid_terminal). Writes lshape_mlts_convergence_hybridterminal_k_<k>.txt
// (formerly shared with GlobalExact, see RestrictedExact).
void ERK4_MLTS_Lshape_HybridTerminal_conv_test(int argc, char **argv) {
    lts::mlts_sweep_config cfg = mlts_lshape_config("lshape_mlts_convergence_hybridterminal", false);
    cfg.cfl_factor = 0.002;
    lts::mlts_sweep(argc, argv, cfg, lts::lshape_corner_solution(),
        [](auto& pb, const lts::bands& bd) { return mlts::make_hybrid_terminal(pb, bd, 3, true); });
}

// Superseded 5-point mesh family attempt: 5 points LOG-EQUIDISTANT in
// h_max (geometric spacing), spanning N=0..N=2's h-range only (h_max in
// [0.088, 0.354]) -- NOT out to N=3/N=4, whose p_global explodes
// exponentially (32^N) as an intrinsic consequence of the graded-mesh
// theory (matching the corner singularity's regularity for k=3); reaching
// a log-equidistant point near N=3's h would cost about as much as N=3
// itself.
//   N=0  (n0=4, ell=0): h_max=0.353553  p_global=1
//   0b   (n0=6, ell=0): h_max=0.235702  p_global~16
//   N=1  (n0=4, ell=1): h_max=0.176777  p_global=32
//   1b   (n0=6, ell=1): h_max=0.117851  p_global~512
//   N=2  (n0=4, ell=2): h_max=0.088388  p_global=1024
// Successive ratios 1.50 / 1.333 / 1.50 / 1.333 (geometric mean per step =
// sqrt(2) exactly, the closest achievable to perfect log-equidistance
// given n0 must be an even integer). The "0b"/"1b" meshes use a
// theory-corrected corner depth, total_depth =
// round(1.5*(k+1)*(ell + log2(n0/4))), anchored to the validated n0=4
// reference so the SAME order-(k+1) grading applies at every point (see
// generate_lshape_graded_mesh.py's docstring). All points use the GLOBAL
// EXACT algorithm, so the observed order is the TRUE discretization order.
// MLTS_LEVELS (space-separated labels) overrides the point list.
void ERK4_MLTS_Lshape_FivePoint_conv_test(int argc, char **argv) {
    lts::mlts_sweep_config cfg = mlts_lshape_config("lshape_mlts_convergence_5pt", true);
    cfg.levels = {"0", "0b", "1", "1b", "2"};
    if (const char* env_n = std::getenv("MLTS_LEVELS")) {
        cfg.levels.clear();
        std::istringstream iss(env_n);
        std::string tok;
        while (iss >> tok) cfg.levels.push_back(tok);
    }
    cfg.level_column = "label";
    lts::mlts_sweep(argc, argv, cfg, lts::lshape_corner_solution(), mlts_global_exact);
}

// Diagnostic on the real N=2 mesh: one restricted macro step (8-ring
// sub-blocks, _coarse_v2 + _advance_uncovered_full, MLTS_DIAG_NT steps,
// default 1) against the GLOBAL exact recursion run on a copy of the SAME
// initial state over the SAME step(s) -- so any divergence is due ONLY to
// the restriction itself, not to mesh/threshold/CFL differences. Prints
// band0's restricted vs exact own_w_i (cells, faces, halo faces) at the
// first sub-step and the per-band max|restricted - global| afterwards.
// Writes lshape_mlts_diagN2_k_<k>.txt (formerly shared with GlobalExact).
void ERK4_MLTS_Lshape_DiagN2_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   MULTI-LEVEL LTS L-SHAPE CONVERGENCE SWEEP (N=0..4)" << std::endl << std::endl;

    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);
    const size_t mesh_k = sim_data.m_k_degree;

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);
    std::ostringstream conv_fname;
    conv_fname << out_dir << "/lshape_mlts_diagN2_k_" << mesh_k << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# N  h_max  h_min  p_global  L_levels  L2_error_pressure  wall_time_s\n";

    for (int N : {2}) {

        std::cout << bold << red << "\n   ================ MESH LEVEL N=" << N << " ================" << reset << std::endl;

        std::ostringstream mesh_fname;
        mesh_fname << mesh_root << "lshape_graded/meshes/lshape_graded_k" << mesh_k << "_N" << N << ".txt";
        acoustic_problem<lshape_corner_solution> pb(mesh_fname.str());
        const RealType h_max = pb.h_max, h_min = pb.h_min;
        const int p_global = pb.p_global;

        std::vector<RealType> thresholds = natural_thresholds(h_max, h_min);
        const int L = (int)thresholds.size() + 1;
        std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p_global=" << p_global << "  L=" << L << " levels" << reset << std::endl;

        pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
        erk_type& erk_an = pb.erk();
        vector_type& x_dof = pb.x_dof;
        bands bd = make_bands(pb, thresholds);

        std::vector<erk_type::LTS_subblock_set> blocks(L);
        for (int b = 0; b < L; ++b)
            blocks[b] = erk_an.build_LTS_subblock(bd.own_c[b], bd.own_f[b], 8);
        print_bands(bd);

        std::vector<RealType> dt_level;
        std::vector<int> p_level;
        level_steps(env_real("MLTS_CFL_FACTOR", 0.002), h_max, h_min, thresholds, dt_level, p_level);
        RealType ti = 0.0, tf = 0.05;
        RealType dt_macro = dt_level[0];
        size_t nt = (size_t)env_int("MLTS_DIAG_NT", 1);
        macro_steps(tf, dt_macro);
        std::cout << bold << cyan << "      nt=" << nt << " (DIAG single-step) dt_macro=" << dt_macro << reset << std::endl;

        const size_t n_dof = x_dof.rows();
        SparseMatrix<RealType> P_term = band_projector(bd.own_c[L-1], bd.own_f[L-1], erk_an.n_c_dof(), n_dof);
        std::vector<std::vector<int>> active(L);
        for (int b = 0; b < L; ++b) active[b] = merged_dofs(bd.own_c[b], bd.own_f[b], erk_an.n_c_dof());

        auto rhs3 = [&](RealType t0, RealType dtau) { return mlts::rk4_coarse_rhs(pb, t0, dtau); };

        // Band0's restricted own_w_i vs the exact (_coarse_global) one at
        // the first sub-step of every macro step, cell+face and halo faces
        // (the leak into an adjacent band's own interface face).
        auto diag_band0 = [&](const vector_type& xn, const std::vector<vector_type>& own_w_i, const std::array<vector_type, 3>& F, RealType dtau_i) {
            std::vector<vector_type> own_w_i_exact = zero_weights(n_dof);
            erk_an.erk_weight_LTS_coarse_global(xn, active[0], own_w_i_exact, F[0], F[1], F[2], dtau_i);
            std::cout << bold << yellow << "      [DIAG] band0 own_w_i restricted vs exact, cell+face:" << reset << std::endl;
            for (int j = 0; j < 4; ++j) {
                RealType max_c = 0.0, max_f = 0.0;
                for (int g : bd.own_c[0]) max_c = std::max(max_c, std::abs(own_w_i[j](g) - own_w_i_exact[j](g)));
                for (int g : bd.own_f[0]) { int gg=(int)erk_an.n_c_dof()+g; max_f = std::max(max_f, std::abs(own_w_i[j](gg) - own_w_i_exact[j](gg))); }
                RealType mag_c = 0.0, mag_f = 0.0;
                for (int g : bd.own_c[0]) mag_c = std::max(mag_c, std::abs(own_w_i_exact[j](g)));
                for (int g : bd.own_f[0]) { int gg=(int)erk_an.n_c_dof()+g; mag_f = std::max(mag_f, std::abs(own_w_i_exact[j](gg))); }
                std::cout << bold << cyan << "         w[" << j << "]: max|diff|_c=" << max_c << " (mag " << mag_c
                          << ")  max|diff|_f=" << max_f << " (mag " << mag_f << ")" << reset << std::endl;
            }
            RealType max_halo_f = 0.0, mag_halo_f = 0.0;
            for (size_t hi = blocks[0].n_own_f; hi < blocks[0].active_f.size(); ++hi) {
                int gg = (int)erk_an.n_c_dof() + blocks[0].active_f[hi];
                for (int j = 0; j < 4; ++j) {
                    max_halo_f = std::max(max_halo_f, std::abs(own_w_i[j](gg) - own_w_i_exact[j](gg)));
                    mag_halo_f = std::max(mag_halo_f, std::abs(own_w_i_exact[j](gg)));
                }
            }
            std::cout << bold << cyan << "         halo faces (n=" << (blocks[0].active_f.size()-blocks[0].n_own_f)
                      << "): max|diff|=" << max_halo_f << " (mag " << mag_halo_f << ")" << reset << std::endl;
        };

        // GLOBAL exact recursion on a copy of the same initial state.
        vector_type x_dof_global = x_dof;
        {
            auto coarse_g = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
                auto F = rhs3(t0, dtau_i);
                erk_an.erk_weight_LTS_coarse_global(xn, active[i], own_w_i, F[0], F[1], F[2], dtau_i);
            };
            auto fine_g = [&](vector_type& xn, const std::vector<vector_type>& w, RealType t_start, RealType len) {
                vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
                erk_an.erk_weight_LTS_fine(xn, P_term, w, Fm, Fmh, Fm1, 0.0, len);
            };
            for (size_t it = 1; it <= nt; it++)
                multilevel_step(L, p_level, n_dof, true, coarse_g, fine_g, x_dof_global, ti + dt_macro * (it - 1), dt_macro);
        }

        // Restricted recursion (ERK4_MLTS_Lshape_RestrictedExact_conv_test's).
        size_t level0_calls = 0;
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>& shifted) {
            auto F = rhs3(t0, dtau_i);
            erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w_i, F[0], F[1], F[2], dtau_i, &shifted);
            if (i == 0 && level0_calls++ == 0) diag_band0(xn, own_w_i, F, dtau_i);
        };
        auto fine = [&](vector_type& xn, const std::vector<vector_type>& w, RealType t_start, RealType len) {
            vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
            erk_an.erk_weight_LTS_fine_v2(xn, blocks[L-1], w, Fm, Fmh, Fm1, 0.0, len);
        };
        auto advance = [&](int i, vector_type& xn, const std::vector<vector_type>& own_w_i, RealType dtau_i) {
            erk_an.erk_weight_LTS_coarse_advance_uncovered_full(xn, blocks[i], own_w_i, blocks[L-1].active_c, blocks[L-1].active_f, dtau_i);
        };

        timecounter tc_run;
        tc_run.tic();
        for (size_t it = 1; it <= nt; it++) {
            RealType tn = ti + dt_macro * (it - 1);
            level0_calls = 0;
            multilevel_step(L, p_level, n_dof, false, coarse, fine, x_dof, tn, dt_macro, no_hook(), advance);
            if (it % std::max<size_t>(1, nt/10) == 0)
                std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
        }
        tc_run.toc();

        std::cout << bold << yellow << "      [DIAG] per-band max|restricted - global| after 1 macro-step:" << reset << std::endl;
        for (int b = 0; b < L; ++b) {
            RealType max_c = 0.0, max_f = 0.0;
            int argmax_c = -1, argmax_f = -1;
            for (int g : bd.own_c[b]) {
                RealType d = std::abs(x_dof(g) - x_dof_global(g));
                if (d > max_c) { max_c = d; argmax_c = g; }
            }
            for (int g : bd.own_f[b]) {
                int gg = (int)erk_an.n_c_dof() + g;
                RealType d = std::abs(x_dof(gg) - x_dof_global(gg));
                if (d > max_f) { max_f = d; argmax_f = gg; }
            }
            std::cout << bold << cyan << "         band" << b << ": max|diff|_c=" << max_c
                      << " (dof " << argmax_c << ")  max|diff|_f=" << max_f
                      << " (dof " << argmax_f << ")" << reset << std::endl;
        }

        RealType l2_error = pb.pressure_error(x_dof, tf);
        std::cout << bold << red << "      wall time: " << tc_run << " s   L2_error(pressure) = "
                  << std::setprecision(10) << l2_error << reset << std::endl;
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << p_global << " " << L
                  << " " << l2_error << " " << tc_run.elapsed() << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}
