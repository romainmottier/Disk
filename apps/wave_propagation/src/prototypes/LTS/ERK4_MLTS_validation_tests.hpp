//  Created by Romain Mottier
//
// Development history of the multi-level LTS-RK4 driver on the L-shape
// corner benchmark (lts::lshape_corner_solution, lshape_graded meshes),
// one test per validation step:
//
//   1. ERK4_LTS_v2_L2_validation_test / ERK4_LTS_v2_timing_test: the
//      restricted building blocks (LTS_subblock_set, build_LTS_subblock,
//      erk_weight_LTS_coarse_v2 / _fine_v2 / _advance_uncovered) at L=2.
//   2. ERK4_MLTS_L2/L3/L5small_validation_test: the first generic
//      recursive design (mlts_v1_step) against the unrestricted 2-level
//      reference.
//   3. ERK4_MLTS_timing_test, ERK4_MLTS_L3/L5small_v2design_validation_test,
//      ERK4_MLTS_Lshape_conv_test: the revised design (mlts_v2design_step,
//      every level >= 1 gets a real RK4 step).
//   4. ERK4_MLTS_GlobalGD/RestrictedGD/Lsweep_small/DiagCompare_test: error
//      ratio vs L at fixed p_global (mlts_L_sweep), leading to the
//      Grote-Diaz-faithful recursions of mlts_2026/mlts_steppers.hpp.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0
#pragma once

#include "lts_common.hpp"
#include "mlts_2026/mlts_steppers.hpp"

using mlts_lshape_problem = lts::acoustic_problem<lts::lshape_corner_solution>;

// L-shape problem on lshape_graded_k<k>_N<N>.txt (mesh read, not set up).
static std::unique_ptr<mlts_lshape_problem> mlts_lshape_mesh(const simulation_data& sim_data, int N) {
    std::ostringstream mesh_fname;
    mesh_fname << lts::mesh_root << "lshape_graded/meshes/lshape_graded_k" << sim_data.m_k_degree << "_N" << N << ".txt";
    return std::make_unique<mlts_lshape_problem>(mesh_fname.str());
}

static void mlts_print_level_steps(const std::vector<lts::RealType>& thresholds, lts::RealType h_max, lts::RealType h_min,
                                   const std::vector<lts::RealType>& dt_level, const std::vector<int>& p_level) {
    const int L = (int)thresholds.size() + 1;
    std::vector<lts::RealType> scale(L);
    scale[0] = h_max;
    for (int i = 1; i < L - 1; ++i) scale[i] = thresholds[i-1];
    scale[L-1] = h_min;
    for (int i = 0; i < L - 1; ++i)
        std::cout << bold << cyan << "      level " << i << ": scale=" << scale[i] << " dt=" << dt_level[i] << " p=" << p_level[i] << reset << std::endl;
}

static std::vector<lts::erk_type::LTS_subblock_set> mlts_blocks(lts::erk_type& erk_an, const lts::bands& bd, int halo_rings, bool print_halo) {
    std::vector<lts::erk_type::LTS_subblock_set> blocks(bd.L);
    for (int i = 0; i < bd.L; ++i) {
        blocks[i] = erk_an.build_LTS_subblock(bd.own_c[i], bd.own_f[i], halo_rings);
        std::cout << bold << cyan << "      band" << i << ": own_c=" << bd.own_c[i].size() << " own_f=" << bd.own_f[i].size();
        if (print_halo) std::cout << " halo_c=" << (blocks[i].active_c.size() - blocks[i].n_own_c);
        std::cout << reset << std::endl;
    }
    return blocks;
}

// ---------------------------------------------------------------------------
// Recursive designs
// ---------------------------------------------------------------------------

// First generic design (ERK4_MLTS_timing_test's original recurse()): at
// level i < L-1, add this level's OWN "coarse role" cubic-Taylor fit
// (_coarse_v2, level i's own dt) to the ancestor contribution (both in the
// same local-time coordinate, tau=0 <-> t_start), advance band i's own dofs
// NOT covered by the TERMINAL band's halo in closed form
// (_advance_uncovered) -- coverage must be checked against the terminal
// band: only it actually scatters via _fine_v2, an intermediate band's halo
// is read-only context for its own B-chain -- then loop p_i times,
// Taylor-shifting the combined polynomial to each sub-interval and either
// calling _fine_v2 (base case) or recursing.
template<typename Problem>
static void mlts_v1_step(Problem& pb, const std::vector<lts::erk_type::LTS_subblock_set>& blocks, const std::vector<int>& p_level,
                         lts::vector_type& x, lts::RealType tn, lts::RealType dt_macro) {
    using namespace lts;
    erk_type& erk_an = pb.erk();
    const int L = (int)blocks.size();
    const size_t n_dof = x.rows();
    std::function<void(int, RealType, RealType, vector_type&, const std::vector<vector_type>&)> recurse;
    recurse = [&](int i, RealType t_start, RealType len, vector_type& xn, const std::vector<vector_type>& ancestor_w) {
        std::vector<vector_type> combined = ancestor_w;
        if (i < L - 1) {
            std::vector<vector_type> own_w = zero_weights(n_dof);
            vector_type Fn = pb.rhs(t_start), Fn12 = pb.rhs(t_start + 0.5*len), Fn1 = pb.rhs(t_start + len);
            erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w, Fn, Fn12, Fn1, len);
            for (int j = 0; j < 4; ++j) combined[j] = combined[j] + own_w[j];
            erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, blocks[i], combined, blocks[L-1].active_c, blocks[L-1].active_f, len);
        }
        if (i == L - 1) return;
        int p_i = p_level[i];
        RealType dtau_i = len / p_i;
        for (int m = 0; m < p_i; ++m) {
            RealType tm = m * dtau_i;
            auto shifted = taylor_shift(combined, tm);
            if (i + 1 == L - 1) {
                vector_type Fm = pb.rhs(t_start + tm), Fmh = pb.rhs(t_start + tm + 0.5*dtau_i), Fm1 = pb.rhs(t_start + tm + dtau_i);
                erk_an.erk_weight_LTS_fine_v2(xn, blocks[i+1], shifted, Fm, Fmh, Fm1, 0.0, dtau_i);
            } else {
                recurse(i + 1, t_start + tm, dtau_i, xn, shifted);
            }
        }
    };
    recurse(0, tn, dt_macro, x, zero_weights(n_dof));
}

// REVISED design (validated in ERK4_MLTS_L3_v2design_validation_test: error
// ratio 8.7x -> 1.22x on the controlled L=3 case): level 0 (coarsest) is
// NEVER genuinely time-stepped -- only its own cubic-Taylor fit,
// closed-form-integrated for whatever isn't covered by level 1's halo.
// EVERY level i >= 1 gets a REAL RK4 step (_fine_v2) of size `len` (its own
// dt_level[i] by construction of p_level), forced by the frozen w handed
// down from level i-1; levels i < L-1 ALSO compute their own w_i (BEFORE
// their own _fine_v2 update mutates the state) to hand down to level i+1.
// pass_ancestor: band i's own fit also reflects the ancestor chain's
// ongoing forcing (ancestor_w passed to _coarse_v2, unit-consistent; an
// earlier attempt folding it into the RAW force Fn blew up -- ancestor_w is
// already Mc_inv-scaled). L == 1 (p_global == 1, no threshold): band0 is
// the only band, plain single RK4 step with zero ancestor forcing, as the
// original scheme treats p == 1.
template<typename Problem>
static void mlts_v2design_step(Problem& pb, const std::vector<lts::erk_type::LTS_subblock_set>& blocks, const std::vector<int>& p_level,
                               bool pass_ancestor, lts::vector_type& x, lts::RealType tn, lts::RealType dt_macro) {
    using namespace lts;
    erk_type& erk_an = pb.erk();
    const int L = (int)blocks.size();
    const size_t n_dof = x.rows();
    std::function<void(int, RealType, RealType, vector_type&, const std::vector<vector_type>&)> recurse;
    recurse = [&](int i, RealType t_start, RealType len, vector_type& xn, const std::vector<vector_type>& ancestor_w) {
        if (i == 0 && L == 1) {
            vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
            erk_an.erk_weight_LTS_fine_v2(xn, blocks[0], ancestor_w, Fm, Fmh, Fm1, 0.0, len);
            return;
        }
        if (i == 0) {
            std::vector<vector_type> own_w = zero_weights(n_dof);
            vector_type Fn = pb.rhs(t_start), Fn12 = pb.rhs(t_start + 0.5*len), Fn1 = pb.rhs(t_start + len);
            erk_an.erk_weight_LTS_coarse_v2(xn, blocks[0], own_w, Fn, Fn12, Fn1, len);
            erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, blocks[0], own_w, blocks[1].active_c, blocks[1].active_f, len);
            int p_0 = p_level[0];
            RealType dtau_0 = len / p_0;
            for (int m = 0; m < p_0; ++m) {
                RealType tm = m * dtau_0;
                auto shifted = taylor_shift(own_w, tm);
                recurse(1, t_start + tm, dtau_0, xn, shifted);
            }
            return;
        }
        std::vector<vector_type> own_w_i;
        bool has_own_w = (i < L - 1);
        if (has_own_w) {
            own_w_i = zero_weights(n_dof);
            vector_type Fn = pb.rhs(t_start), Fn12 = pb.rhs(t_start + 0.5*len), Fn1 = pb.rhs(t_start + len);
            if (pass_ancestor)
                erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w_i, Fn, Fn12, Fn1, len, &ancestor_w);
            else
                erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w_i, Fn, Fn12, Fn1, len);
        }
        vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
        erk_an.erk_weight_LTS_fine_v2(xn, blocks[i], ancestor_w, Fm, Fmh, Fm1, 0.0, len);
        if (!has_own_w) return;
        int p_i = p_level[i];
        RealType dtau_i = len / p_i;
        for (int m = 0; m < p_i; ++m) {
            RealType tm = m * dtau_i;
            auto shifted = taylor_shift(own_w_i, tm);
            recurse(i + 1, t_start + tm, dtau_i, xn, shifted);
        }
    };
    recurse(0, tn, dt_macro, x, zero_weights(n_dof));
}

// ---------------------------------------------------------------------------
// 1. Restricted building blocks at L=2
// ---------------------------------------------------------------------------

// Validation gate for the restricted building blocks (LTS_subblock_set,
// build_LTS_subblock, erk_weight_LTS_coarse_v2, erk_weight_LTS_fine_v2): at
// L=2 (one "coarse" band = assembler.Pcoarse, one "fine" band =
// assembler.Pfine, exactly the existing partition), the v2 functions must
// reproduce the trusted ERK4_LTS_Lshape_MMS_conv_test result bit-for-bit
// (up to roundoff) -- N=1, k=3: L2 error = 2.76864549890122e-4. A mismatch
// means a real bug in the v2 restriction/masking.
//
// nb_layer=0 (no buffer dilation), unlike the production nb_layer=3: on
// this small N=1 mesh the dilation eats the ENTIRE coarse region (36/228
// cells are geometrically >= h_c, but 3 rings of vertex-adjacency dilation
// from the other 192 fine cells consume all of them), so the "trusted"
// reference never exercised a nonempty coarse band. nb_layer=0 keeps a
// genuine ~36-cell coarse region so this validation stresses the coarse
// role / halo / masking logic.
void ERK4_LTS_v2_L2_validation_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   ERK4-LTS v2 BUILDING-BLOCK VALIDATION (L=2, N=1)" << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 1);
    mlts_lshape_problem& pb = *pb_ptr;
    mesh_type& msh = pb.msh;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p = pb.p_global;
    RealType h_c = (p == 1) ? 1.25 * h_max : 0.75 * h_max;
    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p=" << p << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    assembler_type& assembler = pb.assembler();
    erk_type& erk_an = pb.erk();
    vector_type& x_dof = pb.x_dof;
    assembler.assemble_P(msh, h_c, 0);

    // ---------------- Build the L=2 LTS_subblock_set pair from Pcoarse/Pfine ----------------
    auto extract_active = [&](const SparseMatrix<RealType>& P, size_t n_c_dof) {
        std::vector<int> c_idx, f_idx;
        for (int k = 0; k < P.outerSize(); ++k)
            for (typename SparseMatrix<RealType>::InnerIterator it(P, k); it; ++it) {
                int i = (int)it.row();
                if (i < (int)n_c_dof) c_idx.push_back(i);
                else f_idx.push_back(i - (int)n_c_dof);
            }
        return std::make_pair(c_idx, f_idx);
    };
    auto [coarse_c, coarse_f] = extract_active(assembler.Pcoarse, erk_an.n_c_dof());
    auto [fine_c,   fine_f]   = extract_active(assembler.Pfine,   erk_an.n_c_dof());
    std::cout << bold << cyan << "      band0 (coarse): " << coarse_c.size() << " cells, " << coarse_f.size() << " faces" << reset << std::endl;
    std::cout << bold << cyan << "      band1 (fine):   " << fine_c.size()   << " cells, " << fine_f.size()   << " faces" << reset << std::endl;
    auto blocks_coarse = erk_an.build_LTS_subblock(coarse_c, coarse_f, 6);
    auto blocks_fine   = erk_an.build_LTS_subblock(fine_c,   fine_f,   6);

    // ---------------- Time discretization (identical to the trusted MMS test) ----------------
    RealType ti = 0.0, tf = 0.05;
    RealType dt = 0.002 * h_max;
    size_t nt = macro_steps(tf, dt);
    RealType dtau = dt / p;
    auto eval_F = [&](RealType t_abs) -> vector_type { return pb.rhs(t_abs); };

    // ---------------- DIAGNOSTIC: direct wk[] comparison, coarse role only ----------------
    {
        vector_type Fn   = eval_F(ti);
        vector_type Fn12 = eval_F(ti + 0.5 * dt);
        vector_type Fn1  = eval_F(ti + dt);
        std::vector<vector_type> wk_orig = zero_weights(x_dof.rows()), wk_v2 = zero_weights(x_dof.rows());
        erk_an.ZeroFc();
        erk_an.erk_weight_LTS_coarse(x_dof, assembler.Pcoarse, wk_orig, Fn, Fn12, Fn1, dt);
        erk_an.erk_weight_LTS_coarse_v2(x_dof, blocks_coarse, wk_v2, Fn, Fn12, Fn1, dt);
        std::cout << bold << red << "\n   [DIAGNOSTIC] wk[] comparison (coarse role only, first macro-step):" << reset << std::endl;
        for (int j = 0; j < 4; ++j) {
            vector_type diff = wk_orig[j] - wk_v2[j];
            RealType max_diff = diff.cwiseAbs().maxCoeff();
            RealType max_orig = wk_orig[j].cwiseAbs().maxCoeff();
            // Where is the worst mismatch, cell or face?
            int worst_idx = 0; RealType worst_val = 0;
            for (int i = 0; i < diff.rows(); ++i) if (std::abs(diff(i)) > worst_val) { worst_val = std::abs(diff(i)); worst_idx = i; }
            bool worst_is_face = (size_t)worst_idx >= erk_an.n_c_dof();
            bool worst_in_coarse_own = false, worst_in_coarse_halo = false;
            if (!worst_is_face) {
                for (size_t i = 0; i < blocks_coarse.n_own_c; ++i) if (blocks_coarse.active_c[i] == worst_idx) worst_in_coarse_own = true;
                for (size_t i = blocks_coarse.n_own_c; i < blocks_coarse.active_c.size(); ++i) if (blocks_coarse.active_c[i] == worst_idx) worst_in_coarse_halo = true;
            } else {
                int local_f = worst_idx - (int)erk_an.n_c_dof();
                for (size_t i = 0; i < blocks_coarse.n_own_f; ++i) if (blocks_coarse.active_f[i] == local_f) worst_in_coarse_own = true;
                for (size_t i = blocks_coarse.n_own_f; i < blocks_coarse.active_f.size(); ++i) if (blocks_coarse.active_f[i] == local_f) worst_in_coarse_halo = true;
            }
            std::cout << bold << cyan << "      w[" << j << "]: max|diff|=" << std::setprecision(6) << max_diff
                      << "  max|orig|=" << max_orig
                      << "  worst_idx=" << worst_idx << (worst_is_face ? " (face)" : " (cell)")
                      << "  in_coarse_own=" << worst_in_coarse_own << "  in_coarse_halo=" << worst_in_coarse_halo
                      << reset << std::endl;
        }
    }

    // ---------------- Time marching: ORIGINAL (trusted) vs NEW v2, side by side ----------------
    vector_type x_dof_orig = x_dof;
    vector_type x_dof_v2   = x_dof;
    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt * (it - 1);
        two_level_reference_step(pb, assembler.Pcoarse, assembler.Pfine, p, x_dof_orig, tn, dt, true);
        // ---- new v2 (restricted, multi-level-ready) algorithm ----
        {
            auto x_dof_n = x_dof_v2;
            std::vector<vector_type> wk = zero_weights(x_dof_v2.rows());
            if (p != 1) {
                vector_type Fn   = eval_F(tn);
                vector_type Fn12 = eval_F(tn + 0.5 * dt);
                vector_type Fn1  = eval_F(tn + dt);
                erk_an.erk_weight_LTS_coarse_v2(x_dof_n, blocks_coarse, wk, Fn, Fn12, Fn1, dt);
                erk_an.erk_weight_LTS_coarse_advance_uncovered(x_dof_n, blocks_coarse, wk,
                    blocks_fine.active_c, blocks_fine.active_f, dt);
            }
            for (int m = 0; m < p; m++) {
                RealType tm  =  m        * dtau;
                RealType tmh = (m + 0.5) * dtau;
                RealType tm1 = (m + 1.0) * dtau;
                vector_type Fm  = eval_F(tn + tm);
                vector_type Fmh = eval_F(tn + tmh);
                vector_type Fm1 = eval_F(tn + tm1);
                erk_an.erk_weight_LTS_fine_v2(x_dof_n, blocks_fine, wk, Fm, Fmh, Fm1, tm, dtau);
            }
            x_dof_v2 = x_dof_n;
        }
        if (it == 1) {
            vector_type diff = (x_dof_orig - x_dof_v2).cwiseAbs();
            RealType d = diff.maxCoeff();
            std::cout << bold << cyan << "      [step " << it << "] max|orig-v2| = " << std::setprecision(6) << d << reset << std::endl;
            std::set<int> band0_own(coarse_c.begin(), coarse_c.end());
            std::set<int> band0_own_f(coarse_f.begin(), coarse_f.end());
            std::set<int> band1_own(fine_c.begin(), fine_c.end());
            std::set<int> band1_own_f(fine_f.begin(), fine_f.end());
            std::set<int> band1_halo_c(blocks_fine.active_c.begin() + blocks_fine.n_own_c, blocks_fine.active_c.end());
            std::set<int> band1_halo_f(blocks_fine.active_f.begin() + blocks_fine.n_own_f, blocks_fine.active_f.end());
            std::cout << bold << yellow << "      band1 halo size: " << band1_halo_c.size() << " cells, " << band1_halo_f.size() << " faces" << reset << std::endl;
            // Categorize every dof by its |diff| magnitude and band membership
            std::map<std::string, RealType> max_by_cat;
            std::map<std::string, int> count_by_cat;
            for (int i = 0; i < diff.rows(); ++i) {
                if (diff(i) < 1e-12) continue;
                bool is_face = (size_t)i >= erk_an.n_c_dof();
                int local = is_face ? (i - (int)erk_an.n_c_dof()) : i;
                std::string cat;
                if (!is_face) {
                    if (band0_own.count(local)) cat = "band0_own_cell";
                    else if (band1_own.count(local)) cat = "band1_own_cell";
                    else cat = "UNCLASSIFIED_cell";
                } else {
                    if (band0_own_f.count(local)) cat = "band0_own_face";
                    else if (band1_own_f.count(local)) cat = "band1_own_face";
                    else cat = "UNCLASSIFIED_face";
                }
                if (diff(i) > max_by_cat[cat]) max_by_cat[cat] = diff(i);
                count_by_cat[cat]++;
            }
            for (auto& [cat, mx] : max_by_cat) {
                std::cout << bold << cyan << "         " << cat << ": count=" << count_by_cat[cat] << "  max|diff|=" << std::setprecision(6) << mx << reset << std::endl;
            }
        }
    }

    // ---------------- Direct DOF-vector comparison (strongest possible check) ----------------
    RealType max_abs_diff = (x_dof_orig - x_dof_v2).cwiseAbs().maxCoeff();
    RealType max_abs_val  = x_dof_orig.cwiseAbs().maxCoeff();
    std::cout << bold << red << "\n   max|x_dof_orig - x_dof_v2|                 = "
              << std::setprecision(6) << max_abs_diff << reset << std::endl;
    std::cout << bold << red << "   max|x_dof_orig| (scale reference)          = "
              << std::setprecision(6) << max_abs_val << reset << std::endl;

    // ---------------- Error at t=tf vs closed-form solution, both runs ----------------
    RealType l2_orig = pb.pressure_error(x_dof_orig, tf);
    RealType l2_v2   = pb.pressure_error(x_dof_v2, tf);
    std::cout << bold << red   << "\n   L2_error(pressure, t=" << tf << ") [original, unrestricted] = " << std::setprecision(15) << l2_orig << reset << std::endl;
    std::cout << bold << yellow << "   L2_error(pressure, t=" << tf << ") [new v2, restricted]     = " << std::setprecision(15) << l2_v2   << reset << std::endl;
    bool gate_pass = (max_abs_diff < 1e-9 * std::max<RealType>(1.0, max_abs_val));
    std::cout << bold << (gate_pass ? green : red)
              << "   " << (gate_pass ? "MATCH (gate passes -- v2 building blocks are exact)" : "MISMATCH (gate fails -- do not proceed to L>2)")
              << reset << std::endl << std::endl;
}

// Timing experiment for the validated, cost-restricted 2-level building
// blocks, on N=2 (p=1024) with a DELIBERATELY TIGHT h_c = 16 h_min so the
// "fine" band covers only the last couple of octaves (near h_min), not
// almost the whole mesh like the production h_c=0.75*h_max. The number of
// fine RK4 sub-steps (p) does NOT change with this choice, but the PER-CALL
// cost (restricted to the small fine band + halo) should drop sharply.
// (Outcome: a single coarse/fine split can only span ~1-2 octaves before the
// cubic-in-time Taylor polynomial stops being valid -- NaN even at
// h_c=16*h_min -- which motivated the multi-level driver.)
//
// The bands are classified directly (unrounded diameter vs h_c, a face is
// fine if any adjacent cell is fine) rather than with assemble_P, which
// rounds both the cell diameter and h_c to 3 decimal places
// (elastoacoustic_four_fields_assembler.hpp, "round(diameter*1000)/1000")
// and so silently zeroes out any h_c below ~5e-4 -- every graded mesh here
// has h_min well under that (8.6e-5 at N=2, 1.3e-6 at N=3), band1 came back
// with 0 dofs.
void ERK4_LTS_v2_timing_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   ERK4-LTS v2 TIMING EXPERIMENT (N=2, tight h_c)" << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 2);
    mlts_lshape_problem& pb = *pb_ptr;
    mesh_type& msh = pb.msh;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p = pb.p_global;
    RealType h_c = 16.0 * h_min;
    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p=" << p << "  h_c=" << h_c << " (4*h_min)" << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    erk_type& erk_an = pb.erk();
    vector_type& x_dof = pb.x_dof;

    bands bd = make_bands(pb, {h_c});
    const std::vector<int>& coarse_c = bd.own_c[0];
    const std::vector<int>& fine_c   = bd.own_c[1];
    std::cout << bold << cyan << "      band0 (coarse): " << coarse_c.size() << " cell-dofs" << reset << std::endl;
    std::cout << bold << cyan << "      band1 (fine):   " << fine_c.size()   << " cell-dofs (vs "
              << msh.cells_size() << " mesh cells total)" << reset << std::endl;

    timecounter tc;
    tc.tic();
    auto blocks_coarse = erk_an.build_LTS_subblock(bd.own_c[0], bd.own_f[0], 6);
    auto blocks_fine   = erk_an.build_LTS_subblock(bd.own_c[1], bd.own_f[1], 6);
    tc.toc();
    std::cout << bold << cyan << "      band1 halo reach: " << (blocks_fine.active_c.size() - blocks_fine.n_own_c)
              << " cell-dofs beyond own  (block-build time: " << tc << " s)" << reset << std::endl;

    RealType ti = 0.0, tf = 0.05;
    RealType dt = 0.002 * h_max;
    size_t nt = macro_steps(tf, dt);
    RealType dtau = dt / p;
    std::cout << bold << cyan << "      nt=" << nt << "  p=" << p << "  (total fine sub-steps = " << nt*p << ")" << reset << std::endl;
    auto eval_F = [&](RealType t_abs) -> vector_type { return pb.rhs(t_abs); };

    timecounter tc_run;
    tc_run.tic();
    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt * (it - 1);
        auto x_dof_n = x_dof;
        std::vector<vector_type> wk = zero_weights(x_dof.rows());
        if (p != 1) {
            vector_type Fn   = eval_F(tn);
            vector_type Fn12 = eval_F(tn + 0.5 * dt);
            vector_type Fn1  = eval_F(tn + dt);
            erk_an.erk_weight_LTS_coarse_v2(x_dof_n, blocks_coarse, wk, Fn, Fn12, Fn1, dt);
            erk_an.erk_weight_LTS_coarse_advance_uncovered(x_dof_n, blocks_coarse, wk,
                blocks_fine.active_c, blocks_fine.active_f, dt);
        }
        for (int m = 0; m < p; m++) {
            RealType tm  =  m        * dtau;
            RealType tmh = (m + 0.5) * dtau;
            RealType tm1 = (m + 1.0) * dtau;
            vector_type Fm  = eval_F(tn + tm);
            vector_type Fmh = eval_F(tn + tmh);
            vector_type Fm1 = eval_F(tn + tm1);
            erk_an.erk_weight_LTS_fine_v2(x_dof_n, blocks_fine, wk, Fm, Fmh, Fm1, tm, dtau);
        }
        x_dof = x_dof_n;
        if (it % std::max<size_t>(1, nt/10) == 0)
            std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
    }
    tc_run.toc();

    RealType l2_error = pb.pressure_error(x_dof, tf);
    std::cout << bold << red << "\n   TOTAL TIME-MARCHING WALL TIME: " << tc_run << " s" << reset << std::endl;
    std::cout << bold << red << "   L2_error(pressure, t=" << tf << ") = " << std::setprecision(10) << l2_error << reset << std::endl << std::endl;
}

// ---------------------------------------------------------------------------
// 2./3. Generic recursive designs vs the unrestricted 2-level reference
// ---------------------------------------------------------------------------

// Regression gate for the first generic recursive driver specialized to
// L=2 (single threshold, exactly ERK4_LTS_v2_L2_validation_test's
// bit-exact setup: N=1, h_c=0.75*h_max, nb_layer=0): it must STILL match the
// original, unrestricted 2-level scheme bit-exactly. If the general L-level
// code correctly reduces to the validated L=2 case, the recursion adds no
// new bug for L>2 -- the accuracy gap at L=5 in ERK4_MLTS_timing_test would
// then be a tuning question (octaves per level, cfl_factor).
void ERK4_MLTS_L2_validation_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   MULTI-LEVEL DRIVER REGRESSION GATE (L=2, N=1)" << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 1);
    mlts_lshape_problem& pb = *pb_ptr;
    mesh_type& msh = pb.msh;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p_global = pb.p_global;
    RealType h_c = (p_global == 1) ? 1.25 * h_max : 0.75 * h_max; // matches ERK4_LTS_v2_L2_validation_test exactly
    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << "  h_c=" << h_c << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    erk_type& erk_an = pb.erk();
    const size_t n_dof = pb.x_dof.rows();
    bands bd = make_bands(pb, {h_c});
    auto blocks = mlts_blocks(erk_an, bd, 6, false);
    std::vector<int> p_level = {p_global};

    RealType ti = 0.0, tf = 0.05;
    RealType dt_macro = 0.002 * h_max;
    size_t nt = macro_steps(tf, dt_macro);

    // ---------------- ORIGINAL (trusted) vs GENERIC-RECURSIVE-AT-L=2, side by side ----------------
    SparseMatrix<RealType> Pcoarse = bands_projector(bd, 0, 1, erk_an.n_c_dof(), n_dof);
    SparseMatrix<RealType> Pfine   = bands_projector(bd, 1, 2, erk_an.n_c_dof(), n_dof);
    vector_type x_dof_orig = pb.x_dof;
    vector_type x_dof_mlts = pb.x_dof;
    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt_macro * (it - 1);
        two_level_reference_step(pb, Pcoarse, Pfine, p_global, x_dof_orig, tn, dt_macro);
        mlts_v1_step(pb, blocks, p_level, x_dof_mlts, tn, dt_macro);
    }
    RealType max_abs_diff = (x_dof_orig - x_dof_mlts).cwiseAbs().maxCoeff();
    RealType max_abs_val  = x_dof_orig.cwiseAbs().maxCoeff();
    std::cout << bold << red << "\n   max|x_dof_orig - x_dof_mlts_at_L2| = " << std::setprecision(6) << max_abs_diff << reset << std::endl;
    std::cout << bold << red << "   max|x_dof_orig| (scale reference)  = " << std::setprecision(6) << max_abs_val << reset << std::endl;
    RealType l2_orig = pb.pressure_error(x_dof_orig, tf);
    RealType l2_mlts = pb.pressure_error(x_dof_mlts, tf);
    std::cout << bold << red    << "\n   L2_error [original]              = " << std::setprecision(15) << l2_orig << reset << std::endl;
    std::cout << bold << yellow << "   L2_error [generic recursive, L=2] = " << std::setprecision(15) << l2_mlts << reset << std::endl;
    bool gate_pass = (max_abs_diff < 1e-9 * std::max<RealType>(1.0, max_abs_val));
    std::cout << bold << (gate_pass ? green : red) << "   " << (gate_pass ? "MATCH (generic recursive driver is exact at L=2)" : "MISMATCH (bug in the generic driver)") << reset << std::endl << std::endl;
}

// N=1 (p_global=32), fixed thresholds, generic recursion step(x, tn, dt,
// blocks, p_level) against the ORIGINAL, unrestricted 2-level scheme with
// Pcoarse=band0, Pfine={every finer band combined}. That reference is NOT
// the same algorithm (it gives every fine band real RK4 dynamics at the tiny
// global dtau), so an exact match is not expected; with this deliberately
// CHEAP, SMALL p_global, "expected approximation accumulation" predicts
// only a mild discrepancy, making a real bug (order-of-magnitude or
// blow-up) easy to tell apart.
template<typename Step>
static void mlts_vs_two_level(int argc, char **argv, const std::string& title, std::vector<lts::RealType> (*split)(lts::RealType, lts::RealType),
                              int halo_rings, Step&& step)
{
    using namespace lts;
    std::cout << std::endl << bold << red << title << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 1);
    mlts_lshape_problem& pb = *pb_ptr;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p_global = pb.p_global;
    std::vector<RealType> thresholds = split(h_max, h_min);
    const int L = (int)thresholds.size() + 1;
    std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    erk_type& erk_an = pb.erk();
    const size_t n_dof = pb.x_dof.rows();
    bands bd = make_bands(pb, thresholds);
    auto blocks = mlts_blocks(erk_an, bd, halo_rings, true);

    std::vector<RealType> dt_level;
    std::vector<int> p_level;
    level_steps(0.002, h_max, h_min, thresholds, dt_level, p_level);
    mlts_print_level_steps(thresholds, h_max, h_min, dt_level, p_level);
    RealType ti = 0.0, tf = 0.05;
    RealType dt_macro = dt_level[0];
    size_t nt = macro_steps(tf, dt_macro);
    dt_level[0] = dt_macro;
    std::cout << bold << cyan << "      nt=" << nt << "  (p0*p1=" << p_level[0]*p_level[1] << " vs p_global=" << p_global << ")" << reset << std::endl;

    // ---------------- reference: ORIGINAL 2-level scheme, Pcoarse=band0, Pfine={band1+...} ----------------
    vector_type x_dof_ref  = pb.x_dof;
    vector_type x_dof_mlts = pb.x_dof;
    SparseMatrix<RealType> Pcoarse_ref = bands_projector(bd, 0, 1, erk_an.n_c_dof(), n_dof);
    SparseMatrix<RealType> Pfine_ref   = bands_projector(bd, 1, L, erk_an.n_c_dof(), n_dof);
    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt_macro * (it - 1);
        two_level_reference_step(pb, Pcoarse_ref, Pfine_ref, p_global, x_dof_ref, tn, dt_macro);
        step(pb, blocks, p_level, x_dof_mlts, tn, dt_macro);
    }

    RealType max_abs_diff = (x_dof_ref - x_dof_mlts).cwiseAbs().maxCoeff();
    RealType max_abs_val  = x_dof_ref.cwiseAbs().maxCoeff();
    print_band_max_diff((x_dof_ref - x_dof_mlts).cwiseAbs(), bd, erk_an.n_c_dof());

    RealType l2_ref  = pb.pressure_error(x_dof_ref, tf);
    RealType l2_mlts = pb.pressure_error(x_dof_mlts, tf);
    std::cout << bold << red << "\n   max|x_dof_ref - x_dof_mlts_L3| = " << std::setprecision(6) << max_abs_diff << reset << std::endl;
    std::cout << bold << red << "   max|x_dof_ref| (scale)         = " << std::setprecision(6) << max_abs_val << reset << std::endl;
    std::cout << bold << red   << "\n   L2_error [2-level reference, Pfine=band1+band2] = " << std::setprecision(10) << l2_ref  << reset << std::endl;
    std::cout << bold << yellow << "   L2_error [generic recursive, L=3]                = " << std::setprecision(10) << l2_mlts << reset << std::endl;
    RealType ratio = l2_mlts / std::max(l2_ref, (RealType)1e-300);
    std::cout << bold << (ratio < 10.0 ? green : red)
              << "   ratio = " << ratio << "  (" << (ratio < 10.0 ? "plausible -- no gross bug" : "SUSPICIOUS -- likely a real bug in L>=3 accumulation")
              << reset << std::endl << std::endl;
}

// L=3: two thresholds splitting the p_global=32 (5 octave) range into 3
// roughly-equal parts (~1.7 octaves each). threshold[0]=h_max/4 gives the
// EXACT integer p0*p1=p_global (p0=4, p1=8, scale[1]=h_max/4=8*h_min), so
// the comparison isn't confounded by a quantization mismatch with the
// reference's total fine-scale resolution; threshold[1] (the geometric
// band1/band2 boundary) does not feed into p0/p1, just above h_min.
static std::vector<lts::RealType> mlts_split_L3(lts::RealType h_max, lts::RealType h_min) { return { h_max / 4.0, 2.0 * h_min }; }
// L=5 at ~1 octave/band.
static std::vector<lts::RealType> mlts_split_L5small(lts::RealType h_max, lts::RealType h_min) { return { h_max / 2.0, h_max / 4.0, h_max / 8.0, 2.0 * h_min }; }

static auto mlts_v1 = [](auto& pb, const auto& blocks, const std::vector<int>& p_level, lts::vector_type& x, lts::RealType tn, lts::RealType dt) {
    mlts_v1_step(pb, blocks, p_level, x, tn, dt);
};

// First case that actually exercises multi-level ancestor accumulation
// (with L=2 there is only ever ONE "coarse role" contributor, so
// ERK4_MLTS_L2_validation_test, bit-exact, structurally cannot catch an
// L>=3 accumulation bug). Ratio ~8.7x.
void ERK4_MLTS_L3_validation_test(int argc, char **argv) {
    mlts_vs_two_level(argc, argv, "   MULTI-LEVEL DRIVER VALIDATION (L=3, N=1)", mlts_split_L3, 6, mlts_v1);
}

void ERK4_MLTS_L5small_validation_test(int argc, char **argv) {
    mlts_vs_two_level(argc, argv, "   MULTI-LEVEL DRIVER VALIDATION (L=5, N=1, small p_global)", mlts_split_L5small, 6, mlts_v1);
}

// Revised design (mlts_v2design_step): ratio 1.22x at L=3.
void ERK4_MLTS_L3_v2design_validation_test(int argc, char **argv) {
    mlts_vs_two_level(argc, argv, "   MULTI-LEVEL DRIVER VALIDATION -- REVISED DESIGN (every level>=1 gets real RK4) (L=3, N=1)",
        mlts_split_L3, 3, [](auto& pb, const auto& blocks, const std::vector<int>& p_level, lts::vector_type& x, lts::RealType tn, lts::RealType dt) {
            mlts_v2design_step(pb, blocks, p_level, false, x, tn, dt);
        });
}

// Revised design at L=5 with the ancestor forcing folded into the
// intermediate levels' own fit: ratio ~9.1x (~9.08x -> ~9.02x thanks to the
// ancestor term, NOT the dominant contributor).
void ERK4_MLTS_L5small_v2design_validation_test(int argc, char **argv) {
    mlts_vs_two_level(argc, argv, "   MULTI-LEVEL DRIVER VALIDATION (L=5, N=1, small p_global)",
        mlts_split_L5small, 3, [](auto& pb, const auto& blocks, const std::vector<int>& p_level, lts::vector_type& x, lts::RealType tn, lts::RealType dt) {
            mlts_v2design_step(pb, blocks, p_level, true, x, tn, dt);
        });
}

// Genuine multi-level (L>2) LTS-RK4 on N=2, revised design. KEY FIX vs the
// failed "just shrink h_c" 2-level experiment (ERK4_LTS_v2_timing_test):
// EVERY level i gets its OWN dt_i, scaled to ITS OWN characteristic cell
// size (thresholds spanning ~2 octaves each), so no level's cubic-Taylor
// fit is ever asked to represent more than a couple of octaves' worth of
// dynamics. (The original unrestricted 2-level scheme gives
// 1.14498093168055e-05 at N=2.)
void ERK4_MLTS_timing_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   MULTI-LEVEL LTS TIMING/STABILITY TEST (N=2)" << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 2);
    mlts_lshape_problem& pb = *pb_ptr;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p_global = pb.p_global;
    std::vector<RealType> thresholds = natural_thresholds(h_max, h_min);
    const int L = (int)thresholds.size() + 1;
    std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << "  L=" << L << " levels" << reset << std::endl;
    for (size_t i = 0; i < thresholds.size(); ++i)
        std::cout << bold << cyan << "        threshold[" << i << "] = " << thresholds[i] << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    bands bd = make_bands(pb, thresholds);
    auto blocks = mlts_blocks(pb.erk(), bd, 3, true);

    std::vector<RealType> dt_level;
    std::vector<int> p_level;
    level_steps(0.002, h_max, h_min, thresholds, dt_level, p_level);
    mlts_print_level_steps(thresholds, h_max, h_min, dt_level, p_level);
    RealType ti = 0.0, tf = 0.05;
    RealType dt_macro = dt_level[0];
    size_t nt = macro_steps(tf, dt_macro);
    dt_level[0] = dt_macro;
    std::cout << bold << cyan << "      nt=" << nt << reset << std::endl;

    timecounter tc_run;
    tc_run.tic();
    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt_macro * (it - 1);
        mlts_v2design_step(pb, blocks, p_level, false, pb.x_dof, tn, dt_macro);
        if (it % std::max<size_t>(1, nt/10) == 0)
            std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
    }
    tc_run.toc();

    RealType l2_error = pb.pressure_error(pb.x_dof, tf);
    std::cout << bold << red << "\n   TOTAL TIME-MARCHING WALL TIME: " << tc_run << " s" << reset << std::endl;
    std::cout << bold << red << "   L2_error(pressure, t=" << tf << ") = " << std::setprecision(10) << l2_error << reset << std::endl;
    std::cout << bold << yellow << "   (for reference, the original unrestricted 2-level scheme gave 1.14498093168055e-05 at N=2)" << reset << std::endl << std::endl;
}

// Convergence sweep N=0..4 of the revised design (error ratio vs the
// trusted unrestricted reference: 1.22x at L=3/N=1, 9.1x at L=5/N=2).
// Writes lshape_mlts_convergence_v2design_k_<k>.txt (it used to share
// lshape_mlts_convergence_k_<k>.txt with mlts_2026's GlobalExact).
void ERK4_MLTS_Lshape_conv_test(int argc, char **argv)
{
    lts::mlts_sweep_config cfg;
    cfg.title       = "   MULTI-LEVEL LTS L-SHAPE CONVERGENCE SWEEP (N=0..4)";
    cfg.mesh_family = "lshape_graded";
    cfg.mesh_prefix = "lshape_graded";
    cfg.levels      = lts::level_labels({0, 1, 2, 3, 4});
    cfg.conv_stem   = "lshape_mlts_convergence_v2design";
    cfg.thresholds  = [](lts::RealType h_max, lts::RealType h_min, int) { return lts::natural_thresholds(h_max, h_min); };
    cfg.cfl_factor  = 0.002;
    cfg.tf          = 0.05;
    lts::mlts_sweep(argc, argv, cfg, lts::lshape_corner_solution(), [](auto& pb, const lts::bands& bd) {
        auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set>>(mlts_blocks(pb.erk(), bd, 3, true));
        return [&pb, blocks](lts::vector_type& x, lts::RealType tn, const lts::mlts_levels& lv) {
            mlts_v2design_step(pb, *blocks, lv.p_level, false, x, tn, lv.dt_macro);
        };
    });
}

// ---------------------------------------------------------------------------
// 4. Error ratio vs L at fixed p_global
// ---------------------------------------------------------------------------

// How the multi-level driver's error ratio (vs the trusted 2-level
// reference, Pcoarse=band0/Pfine=rest) grows with L ALONE, at a FIXED, tiny
// p_global=32 (N=1 mesh): the [h_max/2, h_min] range is split into L-1
// "coarse role" sub-levels, band0 (h >= h_max/2) is IDENTICAL for every L,
// so the 2-level reference is computed once. Motivation: two independent
// single-L tests (L=3 at ~1.7 octaves/band: ratio ~1.22x; L=5 at ~1
// octave/band: ratio ~9.1x, matching the separate L=5/large-p_global
// result) showed the growth tracks L, not p_global or band width; folding
// the ancestor forcing into the intermediate fits moved L=5 only from
// ~9.08x to ~9.02x. Does the ratio grow smoothly with L, or jump at a
// specific L? MLTS_NT1 runs a single macro step.
//
// make_step(pb, bands, p_level, dt_macro) is the variant's per-L
// preparation (and diagnostics), returning step(x, tn).
template<typename MakeStep>
static void mlts_L_sweep(int argc, char **argv, bool band0_change, MakeStep&& make_step)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   MULTI-LEVEL DRIVER RATIO-vs-L SWEEP (N=1, fixed p_global=32)" << std::endl << std::endl;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    auto pb_ptr = mlts_lshape_mesh(sim_data, 1);
    mlts_lshape_problem& pb = *pb_ptr;
    const RealType h_max = pb.h_max, h_min = pb.h_min;
    const int p_global = pb.p_global;
    std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << reset << std::endl;

    pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
    erk_type& erk_an = pb.erk();
    const vector_type x_dof0 = pb.x_dof;
    const size_t n_dof = x_dof0.rows();

    // dt_macro/nt only depend on scale[0]=h_max (same cfl_factor for every
    // L), and the reference's split uses ONLY the h_max/2 threshold: both are
    // identical across the sweep, computed once.
    const RealType cfl_factor = 0.002;
    RealType dt_macro = cfl_factor * h_max;
    RealType ti = 0.0, tf = 0.05;
    size_t nt = macro_steps(tf, dt_macro);
    if (std::getenv("MLTS_NT1")) { nt = 1; tf = dt_macro; }

    bands bd_ref = make_bands(pb, {h_max / 2.0});
    SparseMatrix<RealType> Pcoarse_ref = bands_projector(bd_ref, 0, 1, erk_an.n_c_dof(), n_dof);
    SparseMatrix<RealType> Pfine_ref   = bands_projector(bd_ref, 1, 2, erk_an.n_c_dof(), n_dof);
    vector_type x_dof_ref = x_dof0;
    for (size_t it = 1; it <= nt; it++)
        two_level_reference_step(pb, Pcoarse_ref, Pfine_ref, p_global, x_dof_ref, ti + dt_macro * (it - 1), dt_macro);
    RealType l2_ref = pb.pressure_error(x_dof_ref, tf);
    std::cout << bold << red << "\n   L2_error [2-level reference, band0 vs rest] = " << std::setprecision(10) << l2_ref << reset << std::endl << std::endl;

    // ---------------- sweep over L, same band0 threshold every time ----------------
    std::vector<std::vector<RealType>> threshold_sets = {
        { h_max / 2.0 },                                    // L=2
        { h_max / 2.0, h_max / 8.0 },                        // L=3
        { h_max / 2.0, h_max / 4.0, h_max / 8.0 },            // L=4
        { h_max / 2.0, h_max / 4.0, h_max / 8.0, 2.0*h_min }, // L=5
    };
    for (auto & thresholds : threshold_sets) {
        const int L = (int)thresholds.size() + 1;
        std::cout << bold << magenta << "   ---- L=" << L << " ----" << reset << std::endl;
        bands bd = make_bands(pb, thresholds);

        std::vector<RealType> dt_level;
        std::vector<int> p_level;
        level_steps(cfl_factor, h_max, h_min, thresholds, dt_level, p_level);
        dt_level[0] = dt_macro; // exact match to the reference's macro-dt/nt
        p_level[0] = std::max(1, (int)std::round(dt_level[0] / dt_level[1]));

        auto step = make_step(pb, bd, p_level, dt_macro);
        vector_type x_dof_mlts = x_dof0;
        for (size_t it = 1; it <= nt; it++)
            step(x_dof_mlts, ti + dt_macro * (it - 1));

        if (band0_change) {
            RealType max_change = 0.0;
            for (int g : bd.own_c[0]) max_change = std::max(max_change, std::abs(x_dof_mlts(g) - x_dof0(g)));
            for (int g : bd.own_f[0]) { int gg = (int)erk_an.n_c_dof() + g; max_change = std::max(max_change, std::abs(x_dof_mlts(gg) - x_dof0(gg))); }
            std::cout << bold << yellow << "      [DIAG] band0 max|x_final - x_initial| = " << max_change << reset << std::endl;
        }
        print_band_max_diff((x_dof_ref - x_dof_mlts).cwiseAbs(), bd, erk_an.n_c_dof());
        RealType l2_mlts = pb.pressure_error(x_dof_mlts, tf);
        RealType ratio = l2_mlts / std::max(l2_ref, (RealType)1e-300);
        std::cout << bold << yellow << "      L2_error [L=" << L << "] = " << std::setprecision(10) << l2_mlts << reset << std::endl;
        std::cout << bold << (ratio < 3.0 ? green : red) << "      ratio  [L=" << L << "] = " << ratio << reset << std::endl << std::endl;
    }
}

// Grote-Diaz-faithful, GLOBAL-matrix recursion (Algorithm 4, RK4 variant):
// every level i < L-1 recomputes its OWN contribution fresh, from the
// CURRENT global state, at every short sub-step, and ACCUMULATES it with
// what it inherited before Taylor-shifting and recursing; band b's own
// indicator P_band[b] is the direct analog of their (P_l - P_{l+1}) (bands
// are DISJOINT, unlike Grote-Diaz's NESTED P_l). No restricted sub-blocks,
// no "advance uncovered": every level's own dofs are updated purely through
// the leak of the GLOBAL Kcc/Kcf/Kfc operators, exactly as in the paper,
// where only the base-case LTS2/LTS4 call updates the state vector.
void ERK4_MLTS_GlobalGD_test(int argc, char **argv) {
    mlts_L_sweep(argc, argv, true, [](auto& pb, const lts::bands& bd, const std::vector<int>& p_level, lts::RealType dt_macro) {
        auto P_band = std::make_shared<std::vector<SparseMatrix<lts::RealType>>>(bd.L);
        for (int b = 0; b < bd.L; ++b)
            (*P_band)[b] = lts::band_projector(bd.own_c[b], bd.own_f[b], pb.erk().n_c_dof(), pb.x_dof.rows());
        return [&pb, P_band, p_level, dt_macro, L = bd.L](lts::vector_type& x, lts::RealType tn) {
            auto coarse = [&](int i, lts::vector_type& xn, std::vector<lts::vector_type>& own_w_i, lts::RealType t0, lts::RealType dtau_i, const std::vector<lts::vector_type>&) {
                auto F = mlts::rk4_coarse_rhs(pb, t0, dtau_i);
                pb.erk().erk_weight_LTS_coarse_global(xn, (*P_band)[i], own_w_i, F[0], F[1], F[2], dtau_i);
            };
            auto fine = [&](lts::vector_type& xn, const std::vector<lts::vector_type>& w, lts::RealType t_start, lts::RealType len) {
                lts::vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
                pb.erk().erk_weight_LTS_fine(xn, (*P_band)[L-1], w, Fm, Fmh, Fm1, 0.0, len);
            };
            lts::multilevel_step(L, p_level, x.rows(), true, coarse, fine, x, tn, dt_macro);
        };
    });
}

// The same algorithm on RESTRICTED own+3-ring sub-blocks (the recursion of
// mlts::restricted_stepper: _coarse_v2 with ancestor_w,
// _advance_uncovered_full). At L=3, reports how much of band0 the terminal
// band's halo reaches (confirmed: zero overlap even at wide halos).
void ERK4_MLTS_RestrictedGD_test(int argc, char **argv) {
    mlts_L_sweep(argc, argv, true, [](auto& pb, const lts::bands& bd, const std::vector<int>& p_level, lts::RealType dt_macro) {
        const int L = bd.L;
        auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set>>(L);
        for (int b = 0; b < L; ++b) (*blocks)[b] = pb.erk().build_LTS_subblock(bd.own_c[b], bd.own_f[b], 3);
        if (L == 3) {
            const auto& term = (*blocks)[L-1];
            std::set<int> own_c0(bd.own_c[0].begin(), bd.own_c[0].end());
            std::set<int> own_f0(bd.own_f[0].begin(), bd.own_f[0].end());
            int overlap_c = 0, overlap_f = 0;
            for (int g : term.active_c) if (own_c0.count(g)) overlap_c++;
            for (int g : term.active_f) if (own_f0.count(g)) overlap_f++;
            std::cout << bold << yellow << "      [DIAG] blocks[L-1].active_c size=" << term.active_c.size()
                      << " overlap with own_c[0]=" << overlap_c << reset << std::endl;
            std::cout << bold << yellow << "      [DIAG] blocks[L-1].active_f size=" << term.active_f.size()
                      << " overlap with own_f[0]=" << overlap_f << reset << std::endl;
        }
        auto step = mlts::restricted_stepper(pb, blocks);
        return [step, p_level, dt_macro, L](lts::vector_type& x, lts::RealType tn) mutable {
            lts::mlts_levels lv;
            lv.L = L; lv.p_level = p_level; lv.dt_macro = dt_macro;
            step(x, tn, lv);
        };
    });
}

// Grote-Diaz-faithful recursion on restricted sub-blocks: the root level is
// recomputed fresh every short sub-step too (not one fit extrapolated
// across the whole macro step -- Grote-Diaz's base-case LTS2, chained p
// times via its own leap-frog recurrence, achieves p-hop reach for its OWN
// coarse rows, a single frozen fit does not); intermediate levels get NO
// dedicated _fine_v2 step (only the terminal level integrates), their own
// dofs evolving via the leak from deeper levels plus their own closed-form
// "uncovered" update, with coverage checked against the TERMINAL level's
// halo -- intermediate levels no longer leak, so checking against
// blocks[i+1] would wrongly exclude cells nothing then touches (the bug that
// blew up band0/band1 in the first Grote-Diaz-faithful attempt). Note: with
// every level >= 1 leaking into its own halo (the revised design), band0's
// "uncovered" check would need the UNION of every descendant halo
// (blocks[1..L-1]) -- checking blocks[1] alone applies the closed-form
// integral (exact ONLY for cells with zero fine coupling) to interface cells
// that DO have such coupling, the band0-own-dofs error jumping from ~2e-5
// (L=3) to ~8e-3 (L=4/5).
void ERK4_MLTS_Lsweep_small_test(int argc, char **argv) {
    mlts_L_sweep(argc, argv, false, [](auto& pb, const lts::bands& bd, const std::vector<int>& p_level, lts::RealType dt_macro) {
        const int L = bd.L;
        auto blocks = std::make_shared<std::vector<lts::erk_type::LTS_subblock_set>>(L);
        for (int i = 0; i < L; ++i) (*blocks)[i] = pb.erk().build_LTS_subblock(bd.own_c[i], bd.own_f[i], 3);
        return [&pb, blocks, p_level, dt_macro, L](lts::vector_type& x, lts::RealType tn) {
            using namespace lts;
            erk_type& erk_an = pb.erk();
            const auto& bl = *blocks;
            const size_t n_dof = x.rows();
            std::function<void(int, RealType, RealType, vector_type&, const std::vector<vector_type>&)> recurse;
            recurse = [&](int i, RealType t_start, RealType len, vector_type& xn, const std::vector<vector_type>& ancestor_w) {
                if (i == 0) {
                    int p_0 = p_level[0];
                    RealType dtau_0 = len / p_0;
                    for (int m = 0; m < p_0; ++m) {
                        RealType tm = m * dtau_0;
                        std::vector<vector_type> own_w = zero_weights(n_dof);
                        vector_type Fn = pb.rhs(t_start + tm), Fn12 = pb.rhs(t_start + tm + 0.5*dtau_0), Fn1 = pb.rhs(t_start + tm + dtau_0);
                        erk_an.erk_weight_LTS_coarse_v2(xn, bl[0], own_w, Fn, Fn12, Fn1, dtau_0);
                        if (L > 1) erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, bl[0], own_w, bl[L-1].active_c, bl[L-1].active_f, dtau_0);
                        recurse(1, t_start + tm, dtau_0, xn, own_w);
                    }
                    return;
                }
                if (i == L - 1) {
                    vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
                    erk_an.erk_weight_LTS_fine_v2(xn, bl[i], ancestor_w, Fm, Fmh, Fm1, 0.0, len);
                    return;
                }
                int p_i = p_level[i];
                RealType dtau_i = len / p_i;
                for (int m = 0; m < p_i; ++m) {
                    RealType tm = m * dtau_i;
                    auto shifted_ancestor = taylor_shift(ancestor_w, tm);
                    std::vector<vector_type> own_w_i = zero_weights(n_dof);
                    vector_type Fn = pb.rhs(t_start + tm), Fn12 = pb.rhs(t_start + tm + 0.5*dtau_i), Fn1 = pb.rhs(t_start + tm + dtau_i);
                    erk_an.erk_weight_LTS_coarse_v2(xn, bl[i], own_w_i, Fn, Fn12, Fn1, dtau_i, &shifted_ancestor);
                    erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, bl[i], own_w_i, bl[L-1].active_c, bl[L-1].active_f, dtau_i);
                    recurse(i + 1, t_start + tm, dtau_i, xn, own_w_i);
                }
            };
            recurse(0, tn, dt_macro, x, zero_weights(n_dof));
        };
    });
}

// GlobalGD's accumulating recursion with the restricted kernels
// (_coarse_v2 without ancestor, _fine_v2). At L=3, compares band0's first
// restricted fit with the global one (_coarse_global).
void ERK4_MLTS_DiagCompare_test(int argc, char **argv) {
    mlts_L_sweep(argc, argv, false, [](auto& pb, const lts::bands& bd, const std::vector<int>& p_level, lts::RealType dt_macro) {
        using namespace lts;
        const int L = bd.L;
        erk_type& erk_an = pb.erk();
        const size_t n_dof = pb.x_dof.rows();
        auto blocks = std::make_shared<std::vector<erk_type::LTS_subblock_set>>(L);
        for (int b = 0; b < L; ++b) (*blocks)[b] = erk_an.build_LTS_subblock(bd.own_c[b], bd.own_f[b], 3);
        if (L == 3) {
            SparseMatrix<RealType> P0 = band_projector(bd.own_c[0], bd.own_f[0], erk_an.n_c_dof(), n_dof);
            RealType ti = 0.0;
            RealType dtau0 = dt_macro / p_level[0];
            vector_type Fn = pb.rhs(ti), Fn12 = pb.rhs(ti+0.5*dtau0), Fn1 = pb.rhs(ti+dtau0);
            std::vector<vector_type> w_glob = zero_weights(n_dof), w_restr = zero_weights(n_dof);
            vector_type xn0 = pb.x_dof;
            erk_an.erk_weight_LTS_coarse_global(xn0, P0, w_glob, Fn, Fn12, Fn1, dtau0);
            erk_an.erk_weight_LTS_coarse_v2(xn0, (*blocks)[0], w_restr, Fn, Fn12, Fn1, dtau0);
            for (int j=0;j<4;++j) {
                RealType d = (w_glob[j]-w_restr[j]).cwiseAbs().maxCoeff();
                RealType mg = w_glob[j].cwiseAbs().maxCoeff();
                RealType mr = w_restr[j].cwiseAbs().maxCoeff();
                std::cout << bold << yellow << "      [DIAG] j=" << j << " diff=" << d << " max|global|=" << mg << " max|restr|=" << mr << reset << std::endl;
            }
        }
        return [&pb, blocks, p_level, dt_macro, L](vector_type& x, RealType tn) {
            const auto& bl = *blocks;
            auto coarse = [&](int i, vector_type& xn, std::vector<vector_type>& own_w_i, RealType t0, RealType dtau_i, const std::vector<vector_type>&) {
                auto F = mlts::rk4_coarse_rhs(pb, t0, dtau_i);
                pb.erk().erk_weight_LTS_coarse_v2(xn, bl[i], own_w_i, F[0], F[1], F[2], dtau_i);
            };
            auto fine = [&](vector_type& xn, const std::vector<vector_type>& w, RealType t_start, RealType len) {
                vector_type Fm = pb.rhs(t_start), Fmh = pb.rhs(t_start + 0.5*len), Fm1 = pb.rhs(t_start + len);
                pb.erk().erk_weight_LTS_fine_v2(xn, bl[L-1], w, Fm, Fmh, Fm1, 0.0, len);
            };
            multilevel_step(L, p_level, x.rows(), true, coarse, fine, x, tn, dt_macro);
        };
    });
}
