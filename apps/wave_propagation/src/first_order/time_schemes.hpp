//
//  time_schemes.hpp
//
//  First-order formulation (HHOFirstOrder.hpp): the time integrators, ERK(s) and SDIRK(s, s+1) on the shared drivers
//  (time_schemes/runge_kutta.hpp), two-level LTS-RK4 and LI-SDIRK(4,3) on the shared ones (time_schemes/lts.hpp) with
//  the operator of erk_coupling_hho_scheme, their matrix-free versions (solver.matrix_free,
//  time_schemes/mf_erk_operator.hpp), and the choice of the integrator from the parameters (integrate).
//

#pragma once
#ifndef first_order_time_schemes_hpp
#define first_order_time_schemes_hpp

#include "observers.hpp"
#include "../time_schemes/runge_kutta.hpp"
#include "../time_schemes/lts.hpp"

namespace hho_first_order {

    // Called after each time step with the time and the state (energy history)
    typedef std::function<void(RealType, const vector_type&)> observer_type;

    // ERK(s)
    inline vector_type integrate_erk(problem& pb, vector_type x, size_t nt, RealType tf, size_t s, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        runge_kutta::erk(x, nt, tf, a, b, c, [&](RealType t, vector_type& y, vector_type& k) {
            vector_type F = pb.eval_F(t);
            pb.erk().SetFg(F);
            pb.erk().erk_weight(y, k);
        }, [&](RealType t, const vector_type& y) { if (observe) observe(t, y); return true; });
        return x;
    }

    // SDIRK(s, s+1): (M + a_ii dt K) k_i = F(t_i) - K y_i, one factorization since a_ii is constant
    // iterative_tol > 0: BiCGSTAB with incomplete LU (relative residual iterative_tol) instead of the direct solver
    template<typename Problem>
    inline vector_type integrate_sdirk(Problem& pb, vector_type x, size_t nt, RealType tf, size_t s, bool sc_Q,
                                       const observer_type& observe = {}, RealType iterative_tol = -1.0, bool fast = false) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(int(s), a, b, c);
        RealType dt = tf / nt;
        auto& as = pb.assembler();
        // reference mode: no copies of K and M (memory of the large meshes), same arithmetic
        dirk_hho_scheme<RealType> dirk_an(as.LHS, as.RHS, as.MASS, typename dirk_hho_scheme<RealType>::by_reference_t{});
        if (sc_Q)
            dirk_an.set_static_condensation_data(pb.cell_blocks(), pb.n_face_dof());
        if (fast)
            dirk_an.set_fast(true);
        RealType scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        dirk_an.ComposeMatrix();
        if (iterative_tol > 0.0)
            dirk_an.setIterativeSolver(iterative_tol);
        dirk_an.DecomposeMatrix();
        runge_kutta::sdirk(x, nt, tf, a, b, c, [](RealType, const vector_type&) {},
            [&](size_t i, RealType t, vector_type& y, vector_type& k) {
                vector_type F = pb.eval_F(t);
                dirk_an.SetFg(F);
                dirk_an.irk_weight(y, k, dt, a(i,i), true);
            }, [&](RealType t, const vector_type& y) { if (observe) observe(t, y); return true; });
        return x;
    }

    // Two-level LTS-RK4 and LI-SDIRK(3,4) (time_schemes/lts.hpp) on the ERK operator: force masked before the
    // operator, LI slopes through the operator (erk_coupling_hho_scheme::lts_operator)
    inline std::vector<int> coarse_dofs(problem& pb) {
        const SparseMatrix<RealType>& Pcoarse = pb.assembler().Pcoarse;
        std::vector<int> dofs;
        for (int i = 0; i < Pcoarse.rows(); ++i)
            if (Pcoarse.coeff(i, i) > 0.5)
                dofs.push_back(i);
        return dofs;
    }

    // Two-level LTS-RK4: coarse Taylor coefficients w once per step, p fine RK4 substeps
    inline vector_type integrate_lts(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        std::vector<int> coarse = coarse_dofs(pb);
        erk_type::lts_operator op(pb.erk(), &coarse, &pb.assembler().Pfine);
        pb.erk().ZeroFc();
        return local_time_stepping::integrate_lts(op, x, nt, tf, pb.p, [&](RealType t) { return pb.eval_F(t); }, observe);
    }

    // LI-SDIRK(3,4): coarse ERK4 Taylor coefficients, one implicit SDIRK(3,4) step in the fine region
    inline vector_type integrate_li(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(3, a, b, c);
        std::vector<int> coarse = coarse_dofs(pb);
        erk_type::lts_operator op(pb.erk(), &coarse, &pb.assembler().Pfine);
        dirk_coupling_hho_scheme<RealType> dirk_op;
        dirk_op.compose(pb.erk(), pb.assembler().Pfine, a(0,0), tf / nt);
        pb.erk().ZeroFc();
        return local_time_stepping::integrate_li(op, x, nt, tf, a, b, c, [&](RealType t) { return pb.eval_F(t); },
                                                 [&](const vector_type& r) { return dirk_op.solve(r); }, observe);
    }

    // Multilevel LTS-RK4 (time.scheme = "mlts", time_schemes/lts.hpp): levels from the cell sizes (problem::mlts_cell_levels,
    // unknowns of a face at the level of its finest cell), ratio time.level_ratio (or time.level_ratios) between the steps
    // of two consecutive levels; a single level is ERK(4)
    inline std::vector<int> mlts_dof_levels(problem& pb, const wave_parameters& prm, size_t& L) {
        std::vector<int> cell_level;
        L = pb.mlts_cell_levels(prm, cell_level);
        std::vector<size_t> cells(L, 0);
        for (int l : cell_level) cells[size_t(l)]++;
        std::cout << "      MLTS: " << L << " levels, cells per level (coarsest first):";
        for (size_t c : cells) std::cout << " " << c;
        if (prm.time.level_ratios.empty()) std::cout << ", ratio " << prm.time.level_ratio << std::endl;
        else {
            std::cout << ", ratios";
            for (size_t r : mlts_ratios(prm, L)) std::cout << " " << r;
            std::cout << std::endl;
        }
        return pb.assembler().dof_levels(pb.msh, cell_level);
    }

    inline vector_type integrate_mlts(problem& pb, const wave_parameters& prm, vector_type x, size_t nt, RealType tf,
                                      const observer_type& observe = {}) {
        size_t L = 0;
        const std::vector<int> dof_level = mlts_dof_levels(pb, prm, L);
        if (L == 1) return integrate_erk(pb, x, nt, tf, 4, observe);
        std::vector<std::vector<int>> band(L);
        for (size_t i = 0; i < dof_level.size(); i++) band[size_t(dof_level[i])].push_back(int(i));
        SparseMatrix<RealType> P(dof_level.size(), dof_level.size());
        std::vector<Triplet<RealType>> t;
        for (int i : band[L - 1]) t.emplace_back(i, i, 1.0);
        P.setFromTriplets(t.begin(), t.end());
        std::vector<erk_type::lts_operator> ops;
        ops.reserve(L - 1);
        for (size_t l = 0; l + 1 < L; l++) ops.emplace_back(pb.erk(), &band[l], nullptr);
        std::vector<erk_type::lts_operator*> coarse;
        for (auto& op : ops) coarse.push_back(&op);
        erk_type::lts_operator fine(pb.erk(), nullptr, &P);
        pb.erk().ZeroFc();
        return local_time_stepping::integrate_mlts(coarse, fine, x, nt, tf, mlts_ratios(prm, L),
                                                   [&](RealType t) { return pb.eval_F(t); }, observe,
                                                   prm.time.mlts_refresh == "substep");
    }

    // Matrix-free versions (solver.matrix_free, time_schemes/mf_erk_operator.hpp): same schemes, results equal up to round-off.
    // F = 0 (no source, homogeneous Dirichlet data) is detected once and F is then never assembled.
    inline bool zero_force(problem& pb) {
        return pb.constant_F && pb.eval_F(0.0).isZero(0.0);
    }

    inline vector_type integrate_erk_mf(problem& pb, vector_type x, size_t nt, RealType tf, size_t s, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        mf_type& op = pb.mf();
        op.erk(x, nt, tf, a, b, c, [&](RealType t) { return pb.eval_F(t); }, zero_force(pb),
               [&](RealType t, const vector_type& y) { if (observe) observe(t, y); return true; });
        return x;
    }

    inline vector_type integrate_lts_mf(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        mf_type& op = pb.mf();
        op.lts(x, nt, tf, pb.p, pb.fine_unknowns(), [&](RealType t) { return pb.eval_F(t); }, zero_force(pb), observe);
        return x;
    }

    inline vector_type integrate_li_mf(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(3, a, b, c);
        mf_type& op = pb.mf();
        op.li(x, nt, tf, a, b, c, pb.fine_unknowns(), [&](RealType t) { return pb.eval_F(t); }, zero_force(pb), observe);
        return x;
    }

    inline vector_type integrate_mlts_mf(problem& pb, const wave_parameters& prm, vector_type x, size_t nt, RealType tf,
                                         const observer_type& observe = {}) {
        size_t L = 0;
        const std::vector<int> dof_level = mlts_dof_levels(pb, prm, L);
        if (L == 1) return integrate_erk_mf(pb, x, nt, tf, 4, observe);
        mf_type::local_force_type local_F;
        if (pb.separable_F)
            local_F = [&pb](RealType t, const std::vector<int>& dofs, RealType* out) { pb.eval_F_local(t, dofs, out); };
        pb.mf().mlts(x, nt, tf, dof_level, int(L), mlts_ratios(prm, L), [&](RealType t) { return pb.eval_F(t); }, zero_force(pb),
                     observe, prm.time.mlts_refresh == "substep", local_F);
        return x;
    }

    // iterative solver tolerance of the mesh level (space study: multiplied by solver.tolerance_growth per level)
    inline RealType level_tolerance(const wave_parameters& prm, size_t level) {
        if (prm.study.type != "space" || prm.solver.tolerance_growth == 1.0) return prm.solver.tolerance;
        return prm.solver.tolerance * std::pow(prm.solver.tolerance_growth, RealType(level - prm.mesh.level));
    }

    template<typename Problem>
    inline vector_type integrate(Problem& pb, const wave_parameters& prm, vector_type x, size_t nt, const observer_type& observe = {},
                                 size_t level = 0) {
        RealType tf = prm.time.final_time;
        if (prm.time.scheme == "sdirk")
            return integrate_sdirk(pb, x, nt, tf, prm.stages(), prm.solver.static_condensation, observe,
                                   prm.solver.iterative ? level_tolerance(prm, std::max(level, prm.mesh.level)) : -1.0,
                                   prm.solver.matrix_free);
        if constexpr (std::is_same_v<Problem, problem>) {
            if (prm.solver.matrix_free) {
                if (prm.time.scheme == "erk") return integrate_erk_mf(pb, x, nt, tf, prm.stages(), observe);
                if (prm.time.scheme == "lts") return integrate_lts_mf(pb, x, nt, tf, observe);
                if (prm.time.scheme == "mlts") return integrate_mlts_mf(pb, prm, x, nt, tf, observe);
                return integrate_li_mf(pb, x, nt, tf, observe);
            }
            if (prm.time.scheme == "erk")   return integrate_erk(pb, x, nt, tf, prm.stages(), observe);
            if (prm.time.scheme == "lts")   return integrate_lts(pb, x, nt, tf, observe);
            if (prm.time.scheme == "mlts")  return integrate_mlts(pb, prm, x, nt, tf, observe);
            return integrate_li(pb, x, nt, tf, observe);
        }
        throw std::runtime_error("space.elastic_formulation = \"three_fields\": time.scheme = \"sdirk\" only");
    }

}

#endif /* first_order_time_schemes_hpp */
