//  Created by Romain Mottier
//
// Multi-level LTS-RK4 convergence sweeps on square domains with a SMOOTH
// manufactured solution (reproduction_acoustic: u = t^2 sin(pi x)
// sin(pi y), C-infinity everywhere, homogeneous Dirichlet on [0,1]^2 by
// construction), on graded quadtree meshes from mesh_generation/
// (square_graded, corner_square_graded, corner_square_fixed,
// center_square_fixed -- same quadtree/hanging-node machinery and poly_2d
// format as the L-shape generator, without the notch).
//
// PURPOSE of the square problems: isolate whether the restricted drivers'
// residual (~5x on the L-shape corner problem) comes from the algorithm
// itself or from the corner singularity's steep gradients interacting
// badly with a TRUNCATED local reach. Same band/threshold/CFL/recursion
// machinery and same p_global sequence as the L-shape drivers; only the
// mesh's geometry (graded but non-singular) and the solution (smooth)
// change.
//
// The integrators are the steppers of ../mlts_steppers.hpp; most sweeps
// use PellSparse (Diaz-Grote (2015)-style recursion, Algorithm 4/6: nested
// P_l projectors, own-tier term recomputed fresh at every visit,
// translated to our RK4/Taylor "coarse role" building block in place of
// their leap-frog A*z step): erk_weight_LTS_coarse_Pell_sparse (P_l-masked
// but against the FULL GLOBAL Kcc/Kcf/Kfc -- exact by construction) for
// every non-terminal band, the GLOBAL erk_weight_LTS_fine for the
// terminal band. It verifies that this EXACT (not truncated-submatrix)
// construction recovers order ~k+1=4 on the smooth square problem.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0
#pragma once

#include "../mlts_steppers.hpp"

static lts::mlts_sweep_config mlts_square_config(const std::string& family, const std::string& prefix,
                                                 const std::string& stem, int n_max) {
    lts::mlts_sweep_config cfg;
    cfg.title       = "   MULTI-LEVEL LTS SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP (PellSparse)";
    cfg.mesh_family = family;
    cfg.mesh_prefix = prefix;
    cfg.levels      = lts::level_labels(lts::levels(0, n_max));
    cfg.conv_stem   = stem;
    cfg.thresholds  = [](lts::RealType h_max, lts::RealType h_min, int) { return lts::env_thresholds(h_max, h_min); };
    return cfg;
}

static const auto mlts_square_solution = lts::analytic_solution(scal_vec_analytic_functions::EFunctionType::reproduction_acoustic);

static auto mlts_pell_sparse = [](auto& pb, const lts::bands& bd) { return mlts::make_pell_sparse(pb, bd); };

// Hybrid architecture (global coarse levels, terminal restricted to a
// MLTS_TERMINAL_HALO-ring sub-block, default 8, own dofs advanced in
// closed form after the recursion): broken at deep L.
void ERK4_MLTS_Square_conv_test(int argc, char **argv) {
    lts::mlts_sweep_config cfg = mlts_square_config("square_graded", "square_graded", "square_mlts_convergence", 3);
    cfg.title = "   MULTI-LEVEL LTS SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP";
    lts::mlts_sweep(argc, argv, cfg, mlts_square_solution, [](auto& pb, const lts::bands& bd) {
        return mlts::make_hybrid_terminal(pb, bd, lts::env_int("MLTS_TERMINAL_HALO", 8), false);
    });
}

// The MehlinRestricted machinery verbatim (build_LTS_subblock_Pell +
// erk_weight_LTS_coarse_v2/_fine_v2 + _advance_uncovered_full): truncated
// sub-matrices, broken the same way as the L-shape RestrictedLeveled.
void ERK4_MLTS_Square_Pell_conv_test(int argc, char **argv) {
    lts::mlts_sweep_config cfg = mlts_square_config("square_graded", "square_graded", "square_pell_mlts_convergence", 3);
    cfg.title = "   MULTI-LEVEL LTS SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP";
    lts::mlts_sweep(argc, argv, cfg, mlts_square_solution, [](auto& pb, const lts::bands& bd) {
        return mlts::make_mehlin_restricted(pb, bd);
    });
}

// CENTER-graded square (square_graded, abrupt grading): exact sparse P_ell,
// anomalous N=1->N=2 error increase.
void ERK4_MLTS_Square_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("square_graded", "square_graded",
        "square_pellsparse_mlts_convergence", 3), mlts_square_solution, mlts_pell_sparse);
}

// CORNER-graded square (corner_square_graded, abrupt corner_boost grading):
// the graded region touches 2 boundary sides, like the L-shape, to test
// whether the CENTER-graded anomaly is a topology artifact.
void ERK4_MLTS_CornerSquare_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("corner_square_graded", "cornersquare_graded",
        "cornersquare_pellsparse_mlts_convergence", 3), mlts_square_solution, mlts_pell_sparse);
}

// CORNER-graded, fixed local depth (p_global constant).
void ERK4_MLTS_CornerSquareFixed_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("corner_square_fixed", "cornersquarefixed_graded",
        "cornersquarefixed_pellsparse_mlts_convergence", 4), mlts_square_solution, mlts_pell_sparse);
}

// CORNER-graded, local depth ramping +1 per N.
void ERK4_MLTS_CornerSquareRamp_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("corner_square_fixed", "cornersquareramp_graded",
        "cornersquareramp_pellsparse_mlts_convergence", 4), mlts_square_solution, mlts_pell_sparse);
}

// CORNER-graded, coarse start (n0=2) + progressive depth: clean order 4.
void ERK4_MLTS_CornerSquareCoarse_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("corner_square_fixed", "cornersquarecoarse_graded",
        "cornersquarecoarse_pellsparse_mlts_convergence", 4), mlts_square_solution, mlts_pell_sparse);
}

// CENTER-graded, same progressive construction: also clean order 4; the
// source of the "p=2^N growing" curve of the PLOT ARTICLES section.
void ERK4_MLTS_CenterSquareCoarse_PellSparse_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_square_config("center_square_fixed", "centersquarecoarse_graded",
        "centersquarecoarse_pellsparse_mlts_convergence", 4), mlts_square_solution, mlts_pell_sparse);
}

// Two-level (L=2) p-sweep: `pfam` (MLTS_PFAM, default 1) selects a FIXED-p
// mesh family (p_global = 2^pfam, constant across the whole N sweep -- see
// mesh_generation/center_square_fixed's centersq_pfix{pfam}_k3_N{N}.txt),
// and the coarse/fine threshold is hardcoded to a single split
// (background vs. everything locally refined, regardless of how many
// internal quadtree sub-sizes the refined region itself contains).
//
// NOTE: this still runs through the SAME general multilevel recursive
// machinery as the rest of this investigation (PellSparse projectors,
// overlap rings, the generic recursion) -- it is NOT a separate,
// purpose-built two-level-only implementation like
// ERK4_MLTS_CenterSquare_ClassicalTwoLevel_conv_test. What makes this
// "two-level" is only that `thresholds` holds a single split, forcing L=2
// regardless of `pfam`. Do not confuse `p` (the mesh ratio h_max/h_min,
// i.e. how many RK4 substeps happen inside the fine band per coarse step)
// with `L` (the number of levels/bands) -- L is held fixed at 2 for every
// curve in this sweep; only p varies.
void ERK4_MLTS_CenterSquare_TwoLevel_conv_test(int argc, char **argv) {
    const int pfam = lts::env_int("MLTS_PFAM", 1);
    lts::mlts_sweep_config cfg = mlts_square_config("center_square_fixed", "centersq_pfix" + std::to_string(pfam),
        "centersquare_twolevel_p" + std::to_string(1 << pfam) + "_mlts_convergence", 4);
    cfg.level_tag      = " (p=2^" + std::to_string(pfam) + ")";
    cfg.thresholds     = [](lts::RealType h_max, lts::RealType, int) { return std::vector<lts::RealType>{h_max * 0.999}; };
    cfg.velocity_error = true;
    lts::mlts_sweep(argc, argv, cfg, mlts_square_solution, mlts_pell_sparse);
}

// FIXED-level comparison curve: `pfam` (MLTS_PFAM, default 2) selects a
// mesh family whose p_global = 2^pfam is CONSTANT across the whole N sweep
// (same centersq_pfix{pfam}_k3_N{N}.txt meshes as the two-level p-sweep),
// run at its NATURAL (unforced) multilevel L -- for a FIXED extra depth
// this stays roughly constant across N too, so this is the "fixed levels
// tested" curve requested alongside the "variable level" (Coarse/Ramp,
// extra growing with N) curve. MLTS_FORCE_L > 0 forces L instead and is
// folded into the output filename, so several MLTS_FORCE_L runs at the
// same MLTS_PFAM (e.g. the Version B L-sweep, all at pfam=10) do not
// collide if run concurrently.
static lts::mlts_sweep_config mlts_fixedlevel_config(const std::string& rhs_mode) {
    const int pfam = lts::env_int("MLTS_PFAM", 2);
    const int force_l = lts::env_int("MLTS_FORCE_L", 0);
    lts::mlts_sweep_config cfg = mlts_square_config("center_square_fixed", "centersq_pfix" + std::to_string(pfam),
        "centersquare_fixedlevel_p" + std::to_string(1 << pfam)
        + (force_l > 0 ? "_Lforce" + std::to_string(force_l) : "")
        + (rhs_mode.empty() ? "" : "_rhsopt_" + rhs_mode) + "_mlts_convergence", 4);
    cfg.level_tag      = " (fixed p=2^" + std::to_string(pfam) + ")";
    cfg.thresholds     = [force_l](lts::RealType h_max, lts::RealType h_min, int) {
        return force_l > 0 ? lts::forced_thresholds(h_max, force_l) : lts::natural_thresholds(h_max, h_min);
    };
    cfg.velocity_error = true;
    return cfg;
}

void ERK4_MLTS_CenterSquare_FixedLevel_conv_test(int argc, char **argv) {
    lts::mlts_sweep(argc, argv, mlts_fixedlevel_config(""), mlts_square_solution, mlts_pell_sparse);
}

// EXPERIMENTAL (not used by any campaign): isolates two candidate
// RHS-assembly optimizations for non-terminal bands (whose own_w_i is
// masked to Pell_c[i]/Pell_f[i] downstream, so assembling the RHS on the
// OTHER cells is pure waste): (1) restrict assemble_rhs to Pell_c[i]'s
// cells, (2) run assemble_rhs's acoustic cell loop under OpenMP. Selected
// via MLTS_RHS_MODE = "baseline" (default, identical to FixedLevel) |
// "restrict" | "parallel" | "both". Terminal band and ancestor-forcing
// eval_F calls are UNCHANGED in every mode (out of scope here, not
// audited for correctness under restriction/apply_bc interaction -- see
// assemble_rhs_restricted's comment).
// Concluded: NEITHER helps -- restrict gave +11.8% slower at N=3/L=4
// (isolated, no contention), parallel gave x3.1 slower; band0
// (unmaskable by construction) remains the real bottleneck.
void ERK4_MLTS_CenterSquare_FixedLevel_RhsOpt_conv_test(int argc, char **argv) {
    const char* env_mode = std::getenv("MLTS_RHS_MODE");
    const std::string rhs_mode = env_mode ? std::string(env_mode) : "baseline";
    lts::mlts_sweep_config cfg = mlts_fixedlevel_config(rhs_mode);
    cfg.title      = "   MULTI-LEVEL LTS SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP (PellSparse, RHS-OPT EXPERIMENT)";
    cfg.title_note = "      MLTS_RHS_MODE=" + rhs_mode;
    lts::mlts_sweep(argc, argv, cfg, mlts_square_solution, [rhs_mode](auto& pb, const lts::bands& bd) {
        return mlts::make_pell_sparse(pb, bd, [&pb, rhs_mode, L = bd.L](const std::vector<std::vector<int>>& Pell_c) {
            auto active_a_cells = std::make_shared<std::vector<std::vector<char>>>(L, std::vector<char>(pb.msh.cells_size(), 0));
            for (int b = 0; b < L - 1; ++b)
                for (int dof : Pell_c[b])
                    (*active_a_cells)[b][(size_t)dof / pb.cell_dof] = 1;
            return [&pb, rhs_mode, active_a_cells](lts::RealType t_abs, int band_i) -> lts::vector_type {
                auto& assembler = pb.assembler();
                auto s_f_fun_t = pb.sol.set_time(assembler, t_abs);
                if (rhs_mode == "restrict")
                    assembler.assemble_rhs_restricted(pb.msh, lts::null_fun, s_f_fun_t, true, (*active_a_cells)[band_i]);
                else if (rhs_mode == "parallel")
                    assembler.assemble_rhs_parallel(pb.msh, lts::null_fun, s_f_fun_t, true);
                else if (rhs_mode == "both")
                    assembler.assemble_rhs_restricted_parallel(pb.msh, lts::null_fun, s_f_fun_t, true, (*active_a_cells)[band_i]);
                else
                    assembler.assemble_rhs(pb.msh, lts::null_fun, s_f_fun_t, true);
                return assembler.RHS;
            };
        });
    });
}

// Variable p and L: the centersquarehaloscale_graded family
// (generate_haloscale_family.py: ell = N, background halves every N, small
// fixed halo 0.01 for every local pass), L derived from p_global ("p and L
// chosen intelligently, with a relation between the two"): band0
// (background) is ALWAYS exactly ONE raw quadtree level --
// `_uniform_refine` makes every background cell the SAME size regardless
// of how many ell passes happened -- so the very first threshold gap is
// 2x, not 4x. Every band AFTER that spans TWO raw local-refinement levels
// (a 4x gap), matching the 2-cell-thick radial shell established as
// necessary (a 1-level/1-cell-thick shell gets fully swallowed by
// MLTS_OVERLAP_RINGS's dilation, confirmed empirically). p_global =
// 2^extra exactly (ell cancels out of h_max/h_min), so extra =
// log2(p_global) and L = 1 (background) + extra/2 (local pairs),
// overridable via MLTS_FORCE_L. Strict '>' band test: with power-of-2
// thresholds a cell's diameter can land EXACTLY on a threshold value (it
// is the actual size of the next raw quadtree level); '>=' would put it in
// the coarser band. cpu_time_s (getrusage) gives an exact
// classical-vs-multilevel CPU comparison with
// ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_conv_test, which
// runs on the SAME meshes. MLTS_N_MIN/MLTS_N_MAX narrow the sweep
// (N=0..6 by default); a truncated sweep gets a distinct "_Nmin" filename.
// The output keeps the historical "centersquarefixedbg" stem: this driver
// started on the fixed-background family (generate_fixedbg_family.py)
// before moving to the haloscale meshes.
void ERK4_MLTS_CenterSquare_Adaptive_conv_test(int argc, char **argv) {
    lts::mlts_sweep_config cfg = mlts_square_config("center_square_fixed", "centersquarehaloscale_graded",
        "centersquarefixedbg_pellsparse_mlts_convergence", 6);
    cfg.levels = lts::level_labels(lts::levels(0, 6, true));
    if (const char* env_nmin = std::getenv("MLTS_N_MIN")) cfg.conv_suffix = std::string("_Nmin") + env_nmin;
    cfg.thresholds = [](lts::RealType h_max, lts::RealType, int p_global) {
        int L_target = std::getenv("MLTS_FORCE_L")
            ? std::atoi(std::getenv("MLTS_FORCE_L"))
            : 1 + (int)std::round(std::log2((lts::RealType)p_global) / 2.0);
        std::vector<lts::RealType> thresholds;
        if (L_target > 1) thresholds.push_back(h_max / 2.0);
        for (int i = 1; i < L_target - 1; ++i)
            thresholds.push_back(thresholds.back() / 4.0);
        return thresholds;
    };
    cfg.strict_bands   = true;
    cfg.velocity_error = true;
    cfg.cpu_time       = true;
    lts::mlts_sweep(argc, argv, cfg, mlts_square_solution, mlts_pell_sparse);
}
