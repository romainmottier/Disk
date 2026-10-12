//  Created by Romain Mottier
//
// Multi-level LTS-RK4 integrators of the 2026 investigation, as steppers
// for lts::mlts_sweep: each make_*(pb, bands) does the integrator's own
// preparation (projectors, Pell sets, sub-blocks) and band report, and
// returns step(x, tn, levels) advancing x over one macro step.
//
// Two recursion families (lts::multilevel_step):
//   - GLOBAL-matrix Diaz-Grote recursion (accumulate = true): every level
//     i < L-1 recomputes its own contribution fresh from the current
//     global state and SUMS it with the Taylor-shifted inherited one;
//     only the terminal level gets genuine RK4 dynamics, every other
//     level's own dofs being updated through the unrestricted leak of the
//     global Kcc/Kcf/Kfc operators inside erk_weight_LTS_fine's Taylor_w
//     term. Validated to machine precision against the trusted 2-level
//     reference at L=2..5, no dependence on L.
//     Coarse kernels: _global (GlobalExact), _mehlin (MehlinExact),
//     _Pell_sparse (PellSparse).
//   - RESTRICTED recursion (accumulate = false): per-band own+halo
//     sub-blocks, coarse kernel _v2/_v3 folding the inherited weights in
//     through ancestor_w (adding them again would double-count their
//     order-0 term), terminal erk_weight_LTS_fine_v2, and every
//     non-terminal level's own dofs advanced in closed form over its FULL
//     active set by erk_weight_LTS_coarse_advance_uncovered_full. Root
//     cause of the old ~20x residual: _coarse_v2's own_w_i has a genuinely
//     nonzero FACE part at HALO faces (a band's leak into an ADJACENT
//     band's OWN interface face), which the OLD _advance_uncovered
//     silently dropped; the _full version closes that gap (ratio
//     ~20-23x -> 0.99-1.01x at L=2..5 on a small controlled mesh). Every
//     computation stays O(band size) instead of O(n_dof); on the real
//     L-shape meshes a ~5-9x residual remains, traced to genuine
//     Kcc/Kcf/Kfc coupling dropped beyond a fixed-width halo
//     (Mc_inv/Sff_inv are block-diagonal so extracting them is exact; the
//     stiffness blocks are not).
#pragma once

#include "../lts_common.hpp"

namespace mlts {

using lts::RealType;
using lts::vector_type;
using lts::bands;
using lts::mlts_levels;

template<typename Problem>
auto rk4_coarse_rhs(Problem& pb, RealType t0, RealType dtau) {
    return std::array<vector_type, 3>{ pb.rhs(t0), pb.rhs(t0 + 0.5*dtau), pb.rhs(t0 + dtau) };
}

// Terminal level: genuine RK4 step against a GLOBAL diagonal projector
// (erk_weight_LTS_fine), no truncation.
template<typename Problem>
auto global_fine(Problem& pb, std::shared_ptr<SparseMatrix<RealType>> P_term) {
    return [&pb, P_term](vector_type& xn, const std::vector<vector_type>& w, RealType t_start, RealType len) {
        vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
        pb.erk().erk_weight_LTS_fine(xn, *P_term, w, Fm, Fmh, Fm1, 0.0, len);
    };
}

// Terminal level restricted to its own+halo sub-block (erk_weight_LTS_fine_v2).
template<typename Problem, typename Block>
auto restricted_fine(Problem& pb, std::shared_ptr<Block> term) {
    return [&pb, term](vector_type& xn, const std::vector<vector_type>& w, RealType t_start, RealType len) {
        vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
        pb.erk().erk_weight_LTS_fine_v2(xn, *term, w, Fm, Fmh, Fm1, 0.0, len);
    };
}

template<typename Problem>
std::shared_ptr<SparseMatrix<RealType>> terminal_projector(Problem& pb, const bands& bd) {
    return std::make_shared<SparseMatrix<RealType>>(lts::band_projector(bd.own_c[bd.L-1], bd.own_f[bd.L-1], pb.erk().n_c_dof(), pb.x_dof.rows()));
}

// Active-dofs list per band (own_c ++ own_f shifted by n_c_dof, sorted),
// passed directly to erk_weight_LTS_coarse_global's vector-taking overload
// so it never has to re-scan a whole diagonal (O(n_dof)) to rediscover the
// SAME list at every one of the many thousands of calls per macro-step --
// this scan was the single largest cost of the per-step computation at
// N>=2 mesh scale.
template<typename Problem>
std::shared_ptr<std::vector<std::vector<int>>> band_active_dofs(Problem& pb, const bands& bd) {
    auto active = std::make_shared<std::vector<std::vector<int>>>(bd.L);
    for (int b = 0; b < bd.L; ++b)
        (*active)[b] = lts::merged_dofs(bd.own_c[b], bd.own_f[b], pb.erk().n_c_dof());
    return active;
}

// GlobalExact: _coarse_global on each band's own dofs, full Kcc/Kcf/Kfc
// (the proven-exact reference; O(n_dof) per call, no truncation).
template<typename Problem>
auto make_global_exact(Problem& pb, const bands& bd) {
    auto P_term = terminal_projector(pb, bd);
    auto active = band_active_dofs(pb, bd);
    lts::print_bands(bd);
    auto fine = global_fine(pb, P_term);
    return [&pb, active, fine](vector_type& x, RealType tn, const mlts_levels& lv) {
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
            auto F = rk4_coarse_rhs(pb, t0, dtau_i);
            pb.erk().erk_weight_LTS_coarse_global(xn, (*active)[i], own_w_i, F[0], F[1], F[2], dtau_i);
        };
        lts::multilevel_step(lv.L, lv.p_level, x.rows(), true, coarse, fine, x, tn, lv.dt_macro);
    };
}

// MehlinExact: the literal Almquist-Mehlin P_ell masking + overlap
// (_coarse_mehlin on the nested Pell dofs, own-band output), still against
// the GLOBAL matrices.
template<typename Problem>
auto make_mehlin_exact(Problem& pb, const bands& bd) {
    auto P_term = terminal_projector(pb, bd);
    auto active = band_active_dofs(pb, bd);
    lts::print_bands(bd);
    std::vector<std::vector<int>> Pell_c, Pell_f;
    lts::pell_sets(pb.erk(), bd, lts::env_int("MLTS_OVERLAP_RINGS", 3), Pell_c, Pell_f);
    auto Pell_dofs = std::make_shared<std::vector<std::vector<int>>>(bd.L);
    for (int b = 0; b < bd.L; ++b)
        (*Pell_dofs)[b] = lts::merged_dofs(Pell_c[b], Pell_f[b], pb.erk().n_c_dof());
    auto fine = global_fine(pb, P_term);
    return [&pb, active, Pell_dofs, fine](vector_type& x, RealType tn, const mlts_levels& lv) {
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
            auto F = rk4_coarse_rhs(pb, t0, dtau_i);
            pb.erk().erk_weight_LTS_coarse_mehlin(xn, (*Pell_dofs)[i], (*active)[i], own_w_i, F[0], F[1], F[2], dtau_i);
        };
        lts::multilevel_step(lv.L, lv.p_level, x.rows(), true, coarse, fine, x, tn, lv.dt_macro);
    };
}

// PellSparse: MehlinExact's math with a sparse-matvec implementation
// (_coarse_Pell_sparse: P_l-masked, against the FULL GLOBAL Kcc/Kcf/Kfc,
// exact by construction and validated bit-identical to the dense
// reference). rhs_band(i, t) gives the forcing of non-terminal band i
// (plain RHS by default; the RhsOpt experiment restricts/parallelizes it).
template<typename Problem, typename RhsBand>
auto make_pell_sparse(Problem& pb, const bands& bd, RhsBand rhs_band) {
    auto P_term = terminal_projector(pb, bd);
    lts::print_bands(bd);
    auto Pell_c = std::make_shared<std::vector<std::vector<int>>>();
    auto Pell_f = std::make_shared<std::vector<std::vector<int>>>();
    lts::pell_sets(pb.erk(), bd, lts::env_int("MLTS_OVERLAP_RINGS", 3), *Pell_c, *Pell_f);
    auto own_c = std::make_shared<std::vector<std::vector<int>>>(bd.own_c);
    auto own_f = std::make_shared<std::vector<std::vector<int>>>(bd.own_f);
    auto fine = global_fine(pb, P_term);
    auto rhs = rhs_band(*Pell_c);
    return [&pb, Pell_c, Pell_f, own_c, own_f, fine, rhs](vector_type& x, RealType tn, const mlts_levels& lv) {
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
            vector_type Fn = rhs(t0, i), Fn12 = rhs(t0 + 0.5*dtau_i, i), Fn1 = rhs(t0 + dtau_i, i);
            pb.erk().erk_weight_LTS_coarse_Pell_sparse(xn, (*Pell_c)[i], (*Pell_f)[i], (*own_c)[i], (*own_f)[i], own_w_i, Fn, Fn12, Fn1, dtau_i);
        };
        lts::multilevel_step(lv.L, lv.p_level, x.rows(), true, coarse, fine, x, tn, lv.dt_macro);
    };
}

template<typename Problem>
auto make_pell_sparse(Problem& pb, const bands& bd) {
    return make_pell_sparse(pb, bd, [&pb](const std::vector<std::vector<int>>&) {
        return [&pb](RealType t, int) { return pb.rhs(t); };
    });
}

// Global coarse levels + a terminal restricted to its own+halo sub-block
// (halo_rings): ancestor levels keep the GLOBAL exact computation
// (_coarse_global, visited relatively few times, O(n_dof) affordable),
// only the TERMINAL level (visited nt*p_global times -- the dominant cost
// driver) is restricted, via erk_weight_LTS_fine_v2. Since the restricted
// terminal's reach never extends to a far-away ancestor band (confirmed:
// zero overlap even at wide halos), each ancestor band's own dofs outside
// the terminal's coverage are advanced directly by a closed-form integral
// of its own_w_i -- before the recursion (HybridTerminal) or after it
// (the square hybrid driver).
template<typename Problem>
auto make_hybrid_terminal(Problem& pb, const bands& bd, int halo_rings, bool advance_before) {
    auto active = band_active_dofs(pb, bd);
    auto term = std::make_shared<lts::erk_type::LTS_subblock_set>(pb.erk().build_LTS_subblock(bd.own_c[bd.L-1], bd.own_f[bd.L-1], halo_rings));
    auto cov_c = std::make_shared<std::set<int>>(term->active_c.begin(), term->active_c.end());
    auto cov_f = std::make_shared<std::set<int>>(term->active_f.begin(), term->active_f.end());
    auto own_c = std::make_shared<std::vector<std::vector<int>>>(bd.own_c);
    auto own_f = std::make_shared<std::vector<std::vector<int>>>(bd.own_f);
    lts::print_bands(bd);
    auto fine = restricted_fine(pb, term);
    return [&pb, active, cov_c, cov_f, own_c, own_f, fine, advance_before](vector_type& x, RealType tn, const mlts_levels& lv) {
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
            auto F = rk4_coarse_rhs(pb, t0, dtau_i);
            pb.erk().erk_weight_LTS_coarse_global(xn, (*active)[i], own_w_i, F[0], F[1], F[2], dtau_i);
        };
        auto advance = [&](int i, vector_type& xn, const std::vector<vector_type>& own_w_i, RealType dtau_i) {
            lts::advance_own(xn, own_w_i, dtau_i, (*own_c)[i], (*own_f)[i], pb.erk().n_c_dof(), *cov_c, *cov_f);
        };
        if (advance_before)
            lts::multilevel_step(lv.L, lv.p_level, x.rows(), true, coarse, fine, x, tn, lv.dt_macro, advance, lts::no_hook());
        else
            lts::multilevel_step(lv.L, lv.p_level, x.rows(), true, coarse, fine, x, tn, lv.dt_macro, lts::no_hook(), advance);
    };
}

// RESTRICTED recursion on per-band sub-blocks (see the file header); the
// coarse kernel is _coarse_v2 for LTS_subblock_set, _coarse_v3 for the
// P_ell-masked LTS_subblock_set_leveled.
template<typename Problem, typename Block>
auto restricted_stepper(Problem& pb, std::shared_ptr<std::vector<Block>> blocks) {
    auto term = std::make_shared<Block>(blocks->back());
    auto fine = restricted_fine(pb, term);
    return [&pb, blocks, fine](vector_type& x, RealType tn, const mlts_levels& lv) {
        const int L = lv.L;
        auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>& shifted) {
            auto F = rk4_coarse_rhs(pb, t0, dtau_i);
            if constexpr (std::is_same_v<Block, lts::erk_type::LTS_subblock_set_leveled>)
                pb.erk().erk_weight_LTS_coarse_v3(xn, (*blocks)[i], own_w_i, F[0], F[1], F[2], dtau_i, &shifted);
            else
                pb.erk().erk_weight_LTS_coarse_v2(xn, (*blocks)[i], own_w_i, F[0], F[1], F[2], dtau_i, &shifted);
        };
        auto advance = [&](int i, vector_type& xn, const std::vector<vector_type>& own_w_i, RealType dtau_i) {
            pb.erk().erk_weight_LTS_coarse_advance_uncovered_full(xn, (*blocks)[i], own_w_i, (*blocks)[L-1].active_c, (*blocks)[L-1].active_f, dtau_i);
        };
        lts::multilevel_step(L, lv.p_level, x.rows(), false, coarse, fine, x, tn, lv.dt_macro, lts::no_hook(), advance);
    };
}

// RestrictedExact: own + fixed-width BFS halo (halo_rings) sub-blocks.
template<typename Problem>
auto make_restricted_exact(Problem& pb, const bands& bd, int halo_rings) {
    auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set>>(bd.L);
    for (int b = 0; b < bd.L; ++b)
        (*blocks)[b] = pb.erk().build_LTS_subblock(bd.own_c[b], bd.own_f[b], halo_rings);
    lts::print_bands(bd);
    return restricted_stepper(pb, blocks);
}

// MehlinRestricted: sub-blocks built from the EXACT Pell_c/Pell_f domain
// (own + every finer band + overlap, the terminal band included so it gets
// the interface halo its own leak-into-neighbour write needs) instead of a
// guessed fixed-width halo (build_LTS_subblock_Pell). Reach shrinks with
// the band index, so this stays O(band size) for i >= 1.
template<typename Problem>
auto make_mehlin_restricted(Problem& pb, const bands& bd) {
    std::vector<std::vector<int>> Pell_c, Pell_f;
    lts::pell_sets(pb.erk(), bd, lts::env_int("MLTS_OVERLAP_RINGS", 3), Pell_c, Pell_f, true);
    auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set>>(bd.L);
    for (int b = 0; b < bd.L; ++b) {
        (*blocks)[b] = pb.erk().build_LTS_subblock_Pell(bd.own_c[b], bd.own_f[b], Pell_c[b], Pell_f[b]);
        std::cout << bold << cyan << "      band" << b << ": own_c=" << bd.own_c[b].size()
                  << " own_f=" << bd.own_f[b].size() << "  active_c=" << (*blocks)[b].active_c.size()
                  << " active_f=" << (*blocks)[b].active_f.size() << reset << std::endl;
    }
    return restricted_stepper(pb, blocks);
}

// RestrictedLeveled: the halo split into finer/coarser parts
// (build_LTS_subblock_leveled), the B-chain masked to "own+finer" (P_ell)
// at every step by _coarse_v3, the one-time final output leak into the
// coarser halo kept.
template<typename Problem>
auto make_restricted_leveled(Problem& pb, const bands& bd, int halo_rings) {
    lts::erk_type& erk_an = pb.erk();
    std::vector<int> dof_band_c(erk_an.n_c_dof());
    for (int b = 0; b < bd.L; ++b) for (int g : bd.own_c[b]) dof_band_c[g] = b;
    std::vector<int> dof_band_f(erk_an.n_f_dof());
    for (int b = 0; b < bd.L; ++b) for (int g : bd.own_f[b]) dof_band_f[g] = b;

    const int overlap_rings = lts::env_int("MLTS_OVERLAP_RINGS", 3);
    auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set_leveled>>(bd.L);
    for (int b = 0; b < bd.L; ++b) {
        (*blocks)[b] = erk_an.build_LTS_subblock_leveled(bd.own_c[b], bd.own_f[b], halo_rings,
                                                         dof_band_c, dof_band_f, b, overlap_rings);
        std::cout << bold << cyan << "      band" << b << ": own_c=" << bd.own_c[b].size()
                  << " own_f=" << bd.own_f[b].size()
                  << "  (+finer_c=" << (*blocks)[b].n_finer_c << " +coarser_c="
                  << ((*blocks)[b].active_c.size() - (*blocks)[b].n_own_c - (*blocks)[b].n_finer_c) << ")"
                  << reset << std::endl;
    }
    return restricted_stepper(pb, blocks);
}

} // namespace mlts
