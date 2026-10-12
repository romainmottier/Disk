//  Created by Claude (Sonnet 5), for Romain Mottier
//
// LI-SDIRK(4,3) analogue of the classical two-level square sweeps
// (ERK4_MLTS_ClassicalTwoLevel_square_conv_tests.hpp): same mesh family
// (centersquarehaloscale_graded), same solutions, same CFL-based
// dt_macro = cfl_factor*h_max rule, same nb_layer=3, same error routines --
// ONLY the fine-region time integrator differs:
//
//   - erk_weight_LTS_coarse : UNCHANGED, called exactly as in the classical
//     driver (one fixed Pcoarse for this level's whole run), here at every
//     macro step even when p_global == 1.
//   - fine region: the classical driver's `for(m<p_global) erk_weight_LTS_
//     fine(...)` explicit RK4 sub-stepping loop is replaced by ONE call to
//     dirk_weight_LTS_fine (dirk_coupling_hho_scheme.hpp) -- a single
//     implicit SDIRK(3,4) stage recursion over the FULL macro step
//     dt_macro, no sub-stepping (p_global is informational only).
//   - dirk_op.compose(erk_an, Pfine, gamma, dt_macro) once per mesh level
//     (gamma, dt_macro constant within that level's run).
//
// The implicit solve goes through one of four operator classes, all with
// the same compose()/solve() interface: dirk_coupling_hho_scheme (direct
// full (n_c+n_f) factorization), _condensed (static condensation of the
// fine CELL unknowns, paper Sec. 5.2), _condensed_local (halo-restricted
// Schur complement) and _condensed_local_v2 (faster local build). The
// condensed variants must match the direct one to numerical precision;
// their sweeps exist to validate that claim and to measure the speedup.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0
#pragma once

#include "../LTS/mlts_2026/square/ERK4_MLTS_ClassicalTwoLevel_square_conv_tests.hpp"
#include "dirk_coupling_hho_scheme.hpp"
#include "dirk_coupling_hho_scheme_condensed.hpp"
#include "dirk_coupling_hho_scheme_condensed_local.hpp"
#include "dirk_coupling_hho_scheme_condensed_local_v2.hpp"

// Fine region: one SDIRK(3,4) step of size dt_macro through DirkOp, composed
// once per mesh level. compose_report: 0 = silent, 1 = compose() time,
// 2 = compose() time + active cell/face dof counts (halo-restricted ops).
template<template<typename> class DirkOp, int compose_report>
auto lisdirk_fine() {
    return [](auto& pb, auto& eval_F, lts::RealType dt_macro) {
        using lts::RealType;
        auto tab = std::make_shared<std::tuple<Matrix<double, Dynamic, Dynamic>, Matrix<double, Dynamic, 1>, Matrix<double, Dynamic, 1>>>();
        auto& [a_tab, b_tab, c_tab] = *tab;
        dirk_butcher_tableau::sdirk_tables(3, a_tab, b_tab, c_tab);
        RealType gamma = a_tab(0,0);

        timecounter tc_compose;
        tc_compose.tic();
        auto dirk_op = std::make_shared<DirkOp<RealType>>();
        dirk_op->compose(pb.erk(), pb.assembler().Pfine, gamma, dt_macro);
        tc_compose.toc();
        if constexpr (compose_report == 1)
            std::cout << bold << cyan << "      compose() time: " << tc_compose << " s" << reset << std::endl;
        if constexpr (compose_report == 2)
            std::cout << bold << cyan << "      compose() time: " << tc_compose << " s" << "  active_c=" << dirk_op->n_active_c() << "/" << pb.erk().n_c_dof() << "  active_f=" << dirk_op->n_active_f() << "/" << pb.erk().n_f_dof() << reset << std::endl;

        std::function<lts::vector_type(RealType)> eval_F_std = eval_F;
        return [&pb, tab, dirk_op, eval_F_std, gamma, dt_macro](lts::vector_type& x_dof_n, std::vector<lts::vector_type>& w, RealType tn) {
            auto& [a_tab, b_tab, c_tab] = *tab;
            dirk_weight_LTS_fine<RealType>(pb.erk(), *dirk_op, x_dof_n, pb.assembler().Pfine, w, tn, dt_macro, eval_F_std, a_tab, b_tab, c_tab, gamma);
        };
    };
}

static two_level_square_config lisdirk_config(const std::string& title, const std::string& tag, const std::string& stem) {
    two_level_square_config cfg;
    cfg.title         = title;
    cfg.mesh_prefix   = "centersquarehaloscale_graded";
    cfg.level_tag     = tag;
    cfg.conv_stem     = stem;
    cfg.p_global_note = " (informational only -- LI-SDIRK does no sub-stepping)";
    cfg.L_column      = false;
    cfg.nt_dt         = true;
    cfg.coarse_always = true;
    return cfg;
}

static const auto lisdirk_sine = lts::analytic_solution(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);

void LI_SDIRK_CenterSquare_HaloScale_PlaneWave_conv_test(int argc, char **argv) {
    two_level_square_sweep(argc, argv, lisdirk_config(
        "   LI-SDIRK(4,3) SQUARE PLANE-WAVE (sin in x,y,t) CONVERGENCE SWEEP",
        " (LI-SDIRK(4,3), haloscale mesh, plane wave)", "centersquare_lisdirk_haloscale_planewave_convergence"),
        lisdirk_sine, lisdirk_fine<dirk_coupling_hho_scheme, 0>());
}

void LI_SDIRK_CenterSquare_HaloScale_PlaneWave_Condensed_conv_test(int argc, char **argv) {
    two_level_square_sweep(argc, argv, lisdirk_config(
        "   LI-SDIRK(4,3) CONDENSED SQUARE (sin in x,y,t) CONVERGENCE SWEEP",
        " (LI-SDIRK(4,3) CONDENSED, haloscale mesh)", "centersquare_lisdirk_condensed_haloscale_convergence"),
        lisdirk_sine, lisdirk_fine<dirk_coupling_hho_scheme_condensed, 1>());
}

void LI_SDIRK_CenterSquare_HaloScale_PlaneWave_CondensedLocal_conv_test(int argc, char **argv) {
    two_level_square_sweep(argc, argv, lisdirk_config(
        "   LI-SDIRK(4,3) CONDENSED-LOCAL (halo-restricted) SQUARE (sin in x,y,t) CONVERGENCE SWEEP",
        " (LI-SDIRK(4,3) CONDENSED-LOCAL, haloscale mesh)", "centersquare_lisdirk_condensedlocal_haloscale_convergence"),
        lisdirk_sine, lisdirk_fine<dirk_coupling_hho_scheme_condensed_local, 2>());
}

void LI_SDIRK_CenterSquare_HaloScale_PlaneWave_CondensedLocalV2_conv_test(int argc, char **argv) {
    two_level_square_sweep(argc, argv, lisdirk_config(
        "   LI-SDIRK(4,3) CONDENSED-LOCAL-V2 (fast, halo-restricted) SQUARE (sin in x,y,t) CONVERGENCE SWEEP",
        " (LI-SDIRK(4,3) CONDENSED-LOCAL-V2, haloscale mesh)", "centersquare_lisdirk_condensedlocalv2_haloscale_convergence"),
        lisdirk_sine, lisdirk_fine<dirk_coupling_hho_scheme_condensed_local_v2, 2>());
}

void LI_SDIRK_CenterSquare_HaloScale_SlowSine_conv_test(int argc, char **argv) {
    two_level_square_sweep(argc, argv, lisdirk_config(
        "   LI-SDIRK(4,3) SQUARE SLOW-SINE (sin in x,y,t, omega=pi/2) CONVERGENCE SWEEP",
        " (LI-SDIRK(4,3), haloscale mesh, slow sine)", "centersquare_lisdirk_haloscale_slowsine_convergence"),
        lts::slow_sine_solution(), lisdirk_fine<dirk_coupling_hho_scheme, 0>());
}
