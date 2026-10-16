//
//  time_schemes.hpp
//
//  Second-order formulation (HHOSecondOrder.hpp): the time integrators on the first-order system, ERK(s), SDIRK(s,
//  s+1) on the shared drivers (time_schemes/runge_kutta.hpp), two-level LTS-RK4 and LI-SDIRK(4,3) on the shared ones
//  (time_schemes/lts.hpp: force masked after the operator, LI slopes from the stage solves), Newmark on all the
//  unknowns, and the choice of the integrator from the parameters (integrate). faces.localize (F = 0): LTS and MLTS
//  with local evaluations of B (local_operator, second_order_system::eval_local).
//

#pragma once
#ifndef second_order_time_schemes_hpp
#define second_order_time_schemes_hpp

#include "problem.hpp"
#include "../time_schemes/runge_kutta.hpp"
#include "../time_schemes/lts.hpp"
#include <array>

namespace hho_second_order {

    // ERK(s)
    inline vector_type integrate_erk(problem& pb, vector_type y, size_t nt, RealType tf, size_t s) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        second_order_system& sys = pb.system();
        vector_type a_ws;
        runge_kutta::erk(y, nt, tf, a, b, c, [&](RealType t, vector_type& Y, vector_type& k) {
            k = sys.eval(pb.eval_F(t), Y, a_ws);
        }, [](RealType, const vector_type&) { return true; });
        return y;
    }

    // SDIRK(s, s+1): stage Y_i = R_i + a_ii dt f(t_i, Y_i), one factorization since a_ii is constant
    inline vector_type integrate_sdirk(problem& pb, vector_type y, size_t nt, RealType tf, size_t s) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(int(s), a, b, c);
        second_order_system& sys = pb.system();
        RealType dt = tf / nt;
        RealType h = a(0,0) * dt;
        sys.factorize_stage(h, vector_type::Ones(sys.state_size()));
        runge_kutta::sdirk(y, nt, tf, a, b, c, [](RealType, const vector_type&) {}, [&](size_t, RealType t, vector_type& R, vector_type& k) {
            vector_type Y = sys.stage_solve(pb.eval_F(t), R);
            k = (Y - R) / h;
        }, [](RealType, const vector_type&) { return true; });
        return y;
    }

    // Operator of the LTS / LI integrators (time_schemes/lts.hpp): f(F, y) = second_order_system::eval, force masked
    // after the operator (P g(F), g(F) = f(F, 0)), LI slopes from the stage solves; one warm start of the algebraic
    // faces per evaluation site
    struct lts_operator {
        static constexpr local_time_stepping::force_mask  mask  = local_time_stepping::force_mask::after_operator;
        static constexpr local_time_stepping::stage_slope slope = local_time_stepping::stage_slope::from_solve;

        second_order_system& sys;
        const vector_type& Py;
        vector_type Qp, zero, ys;
        std::vector<vector_type> ws;
        bool fine_only;

        lts_operator(problem& pb) : sys(pb.system()), Py(pb.Py), Qp(vector_type::Ones(sys.state_size()) - pb.Py),
            zero(vector_type::Zero(sys.N)), ys(vector_type::Zero(sys.state_size())), ws(local_time_stepping::n_sites),
            fine_only(pb.Py.minCoeff() > 0.5) {}

        vector_type f(const vector_type& F, const vector_type& y, size_t site) { return sys.eval(F, y, ws[site]); }
        vector_type B(const vector_type& y, size_t site) { return sys.eval(zero, y, ws[site]); }
        vector_type g(const vector_type& F, size_t site) { return sys.eval(F, ys, ws[site]); }
        vector_type fine(const vector_type& v) const { return Py.cwiseProduct(v); }
        vector_type coarse(const vector_type& v) const { return Qp.cwiseProduct(v); }
        bool all_fine() const { return fine_only; }
        vector_type zero_state() const { return ys; }
    };

    // Two-level LTS-RK4: p fine RK4 substeps of z' = B P z + P g(t) + w(tau)
    inline vector_type integrate_lts(problem& pb, vector_type y, size_t nt, RealType tf) {
        lts_operator op(pb);
        return local_time_stepping::integrate_lts(op, y, nt, tf, pb.p, [&](RealType t) { return pb.eval_F(t); }, {},
                                                  pb.system().fast());
    }

    // LI-SDIRK(4,3): coarse ERK4 Taylor coefficients, one SDIRK(3,4) step of z' = B P z + P g(t) + w(tau) whose
    // stages (I - a_ii dt B P) Z = R + a_ii dt (P g + w) are solved with the projected stage matrix
    inline vector_type integrate_li(problem& pb, vector_type y, size_t nt, RealType tf) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(3, a, b, c);
        second_order_system& sys = pb.system();
        sys.factorize_stage(a(0,0) * (tf / nt), pb.Py);
        lts_operator op(pb);
        return local_time_stepping::integrate_li(op, y, nt, tf, a, b, c, [&](RealType t) { return pb.eval_F(t); },
                                                 [&](const vector_type& r) { return sys.stage_solve(op.zero, r); }, {},
                                                 sys.fast());
    }

    // Operator of one level of the multilevel LTS-RK4 (local_time_stepping::integrate_mlts): lts_operator with the
    // projectors of the level, Q on the unknowns of a level l < L-1 (its coarse operator) or P on those of the finest
    // level (the fine operator), and the complementary projector
    struct level_operator {
        static constexpr local_time_stepping::force_mask  mask  = local_time_stepping::force_mask::after_operator;
        static constexpr local_time_stepping::stage_slope slope = local_time_stepping::stage_slope::from_solve;

        second_order_system& sys;
        vector_type P, Q, zero, ys;
        std::vector<vector_type> ws;

        level_operator(problem& pb, vector_type P_in, vector_type Q_in) : sys(pb.system()), P(std::move(P_in)),
            Q(std::move(Q_in)), zero(vector_type::Zero(sys.N)), ys(vector_type::Zero(sys.state_size())),
            ws(local_time_stepping::n_sites) {}

        vector_type f(const vector_type& F, const vector_type& y, size_t site) { return sys.eval(F, y, ws[site]); }
        vector_type B(const vector_type& y, size_t site) { return sys.eval(zero, y, ws[site]); }
        vector_type g(const vector_type& F, size_t site) { return sys.eval(F, ys, ws[site]); }
        vector_type fine(const vector_type& v) const { return P.cwiseProduct(v); }
        vector_type coarse(const vector_type& v) const { return Q.cwiseProduct(v); }
        bool all_fine() const { return Q.maxCoeff() < 0.5; }
        vector_type zero_state() const { return ys; }
    };

    // Multilevel LTS-RK4 (time.scheme = "mlts", Algorithm 4 of the preprint LTS_HHO, time_schemes/lts.hpp): levels from
    // the cell sizes (problem::mlts_state_levels), ratios coupling_common::mlts_ratios; a single level is ERK(4)
    inline vector_type integrate_mlts(problem& pb, const wave_parameters& prm, vector_type y, size_t nt, RealType tf) {
        std::vector<int> level;
        const size_t L = pb.mlts_state_levels(prm, level);
        const std::vector<size_t> ratios = mlts_ratios(prm, L);
        std::cout << "      MLTS: " << L << " levels, ratios";
        for (size_t r : ratios) std::cout << " " << r;
        std::cout << std::endl;
        if (L == 1) return integrate_erk(pb, y, nt, tf, 4);
        const size_t n = level.size();
        auto mask = [&](auto&& in) {
            vector_type m = vector_type::Zero(n);
            for (size_t i = 0; i < n; i++) if (in(level[i])) m(i) = 1.0;
            return m;
        };
        const vector_type ones = vector_type::Ones(n);
        std::vector<level_operator> ops;
        ops.reserve(L - 1);
        for (size_t l = 0; l + 1 < L; l++) {
            vector_type Q = mask([l](int v) { return v == int(l); });
            ops.emplace_back(pb, ones - Q, Q);
        }
        vector_type P = mask([L](int v) { return v == int(L) - 1; });
        level_operator fine(pb, P, ones - P);
        std::vector<level_operator*> coarse;
        for (auto& op : ops) coarse.push_back(&op);
        return local_time_stepping::integrate_mlts(coarse, fine, y, nt, tf, ratios, [&](RealType t) { return pb.eval_F(t); },
                                                   {}, prm.time.mlts_refresh == "substep");
    }

    // --------------------------------------------------------------------------------------- faces.localize
    // F = 0 (no source, homogeneous Dirichlet data): checked at three times, F is then never assembled
    inline bool zero_force(problem& pb, RealType tf) {
        for (RealType t : {RealType(0), RealType(0.37) * tf, tf})
            if (!pb.eval_F(t).isZero(0.0)) return false;
        return true;
    }

    // Truncation radius R (rings of cells) of the local evaluations of B y, B^2 y, B^3 y: rho^R / (1 - rho) <= tol,
    // rho the spectral radius of the splitting (contraction of the faces per ring), plus faces.localize_margin
    inline size_t local_radius(problem& pb, const wave_parameters& prm) {
        second_order_system& sys = pb.system();
        const RealType rho = sys.last_splitting_radius > 0.0 ? sys.last_splitting_radius : sys.splitting_spectral_radius();
        if (!(rho > 0.0 && rho < 1.0))
            throw std::runtime_error("faces.localize: the splitting does not contract (rho = " + std::to_string(rho) + ")");
        const RealType r = std::log(sys.tol * (1.0 - rho)) / std::log(rho);
        const size_t R = size_t(std::ceil(r)) + prm.faces.localize_margin;
        std::cout << "      local evaluations: splitting rho = " << rho << ", truncation R = " << R << " rings" << std::endl;
        return R;
    }

    // Operator of the LTS / MLTS integrators with local evaluations (faces.localize, F = 0): P, Q the projectors of the
    // level (lts_operator / level_operator). B P z (fine slopes) and B Q a_i (Taylor coefficients) read only the
    // unknowns of P, Q, the faces solved from zero on the faces reached (exact); B y, B^2 y, B^3 y of the Taylor
    // coefficients are only needed on Q and are computed from the cells of Q and 3R, 2R, R rings around it, faces
    // capped to these regions (truncation error of the order of the tolerance); the force terms vanish.
    struct local_operator {
        static constexpr local_time_stepping::force_mask  mask  = local_time_stepping::force_mask::after_operator;
        static constexpr local_time_stepping::stage_slope slope = local_time_stepping::stage_slope::from_solve;

        second_order_system& sys;
        vector_type P, Q, ys;
        std::vector<int> pc, pd, qc, qd;                 // cells and interface blocks of P, Q
        std::array<std::vector<int>, 3> yc, yd;          // regions read by B y, B^2 y, B^3 y: Q + 3R, 2R, R rings
        std::array<std::vector<char>, 3> cap;            // their algebraic faces
        bool fine_only;

        local_operator(second_order_system& s, vector_type P_in, vector_type Q_in, size_t R)
            : sys(s), P(std::move(P_in)), Q(std::move(Q_in)), ys(vector_type::Zero(s.state_size())), fine_only(Q.maxCoeff() < 0.5) {
            pc = second_order_system::list_of(sys.mask_cells(P)); pd = sys.mask_dblocks(P);
            const std::vector<char> qin = sys.mask_cells(Q);
            qc = second_order_system::list_of(qin); qd = sys.mask_dblocks(Q);
            for (size_t i = 0; i < 3; i++) {
                const std::vector<char> region = sys.grow_cells(qin, (3 - i) * R);
                yc[i] = second_order_system::list_of(region);
                yd[i] = sys.cells_dblocks(region);
                cap[i] = sys.cells_cap(region);
            }
        }

        vector_type f(const vector_type&, const vector_type& y, size_t site) {
            if (site != local_time_stepping::fine_slope_f)
                throw std::logic_error("local_operator::f: unexpected evaluation site");
            return sys.eval_local(y, pc, pd, nullptr);
        }
        vector_type B(const vector_type& y, size_t site) {
            using namespace local_time_stepping;
            switch (site) {
                case B_y:  return sys.eval_local(y, yc[0], yd[0], &cap[0]);
                case B2_y: return sys.eval_local(y, yc[1], yd[1], &cap[1]);
                case B3_y: return sys.eval_local(y, yc[2], yd[2], &cap[2]);
                case w_0: case w_1: case w_2: case w_3: return sys.eval_local(y, qc, qd, nullptr);
                default:   return ys;      // B g_i: g_i = 0 without force
            }
        }
        vector_type g(const vector_type&, size_t) { return ys; }
        vector_type fine(const vector_type& v) const { return P.cwiseProduct(v); }
        vector_type coarse(const vector_type& v) const { return Q.cwiseProduct(v); }
        bool all_fine() const { return fine_only; }
        vector_type zero_state() const { return ys; }
        size_t region_cells(size_t i) const { return yc[i].size(); }

        // Check (environment variable WAVE_LOCALIZE_CHECK): local evaluations against eval from zero faces on the whole
        // mesh, for the state y: B P y and B Q y (exact: round-off expected), the chain B y, B^2 y, B^3 y on Q + 2R, R,
        // 0 rings (truncation: tolerance expected). Prints the largest differences relative to the largest entry.
        void check(const vector_type& y, size_t R, const std::string& name) {
            const vector_type F0 = vector_type::Zero(sys.N);
            auto glob = [&](const vector_type& v) { vector_type a; return sys.eval(F0, v, a); };
            auto rel = [](const vector_type& a, const vector_type& b, const vector_type* rows) {
                RealType d = 0.0, m = 0.0;
                for (Eigen::Index i = 0; i < a.size(); i++) {
                    if (rows && (*rows)(i) < 0.5) continue;
                    d = std::max(d, std::abs(a(i) - b(i))); m = std::max(m, std::abs(b(i)));
                }
                return m > 0.0 ? d / m : d;
            };
            const vector_type Py = fine(y), Qy = coarse(y);
            std::cout << "      check " << name << ": B P y " << rel(f(F0, Py, local_time_stepping::fine_slope_f), glob(Py), nullptr);
            if (!all_fine()) {
                std::cout << ", B Q y " << rel(B(Qy, local_time_stepping::w_0), glob(Qy), nullptr);
                const std::vector<char> qin = sys.mask_cells(Q);
                vector_type v = y, vg = y;
                const size_t sites[3] = {local_time_stepping::B_y, local_time_stepping::B2_y, local_time_stepping::B3_y};
                for (size_t i = 0; i < 3; i++) {
                    v = B(v, sites[i]);
                    vg = glob(vg);
                    // rows of the cells of Q + (2 - i) R rings (u_T, v_T) and of their interface blocks
                    const std::vector<char> out = sys.grow_cells(qin, (2 - i) * R);
                    vector_type rows = vector_type::Zero(sys.state_size());
                    const auto& g = sys.graph();
                    for (size_t c = 0; c < out.size(); c++) if (out[c])
                        for (size_t k = 0; k < g.c_size[c]; k++) { rows(Eigen::Index(g.c_ofs[c] + k)) = 1.0; rows(Eigen::Index(sys.n_c + g.c_ofs[c] + k)) = 1.0; }
                    for (int b : sys.cells_dblocks(out))
                        for (size_t k = 0; k < g.d_size[b]; k++) rows(Eigen::Index(2 * sys.n_c + g.d_ofs[b] + k)) = 1.0;
                    std::cout << ", B^" << i + 1 << " y " << rel(v, vg, &rows);
                }
            }
            std::cout << std::endl;
        }
    };

    inline bool localize_check() { return std::getenv("WAVE_LOCALIZE_CHECK") != nullptr; }

    inline void print_local_stats(second_order_system& sys) {
        if (sys.local_evaluations == 0) return;
        std::cout << "      local evaluations: " << sys.local_evaluations << ", algebraic faces reached "
                  << RealType(sys.local_faces) / RealType(sys.local_evaluations) << " / " << sys.n_a << ", rows "
                  << RealType(sys.local_rows) / RealType(sys.local_evaluations) << " / " << sys.n_c + sys.n_d
                  << " on average" << std::endl;
        sys.local_evaluations = sys.local_faces = sys.local_rows = 0;
    }

    // Two-level LTS-RK4 with local evaluations (coarse coefficients of local_time_stepping::coarse_coefficients)
    inline vector_type integrate_lts_local(problem& pb, const wave_parameters& prm, vector_type y, size_t nt, RealType tf) {
        second_order_system& sys = pb.system();
        sys.build_local();
        const size_t R = local_radius(pb, prm);
        local_operator op(sys, pb.Py, vector_type::Ones(sys.state_size()) - pb.Py, R);
        if (localize_check()) op.check(pb.random_state(), R, "LTS");
        std::cout << "      local LTS: fine cells " << op.pc.size() << ", regions of B y, B^2 y, B^3 y: " << op.region_cells(0)
                  << ", " << op.region_cells(1) << ", " << op.region_cells(2) << " / " << sys.graph().c_ofs.size() << " cells" << std::endl;
        const vector_type F0 = vector_type::Zero(sys.N);
        y = local_time_stepping::integrate_lts(op, y, nt, tf, pb.p, [&](RealType) { return F0; }, {}, false);
        print_local_stats(sys);
        return y;
    }

    // Multilevel LTS-RK4 with local evaluations
    inline vector_type integrate_mlts_local(problem& pb, const wave_parameters& prm, vector_type y, size_t nt, RealType tf) {
        std::vector<int> level;
        const size_t L = pb.mlts_state_levels(prm, level);
        const std::vector<size_t> ratios = mlts_ratios(prm, L);
        std::cout << "      MLTS: " << L << " levels, ratios";
        for (size_t r : ratios) std::cout << " " << r;
        std::cout << std::endl;
        if (L == 1) return integrate_erk(pb, y, nt, tf, 4);
        second_order_system& sys = pb.system();
        sys.build_local();
        const size_t R = local_radius(pb, prm);
        const size_t n = level.size();
        auto mask = [&](auto&& in) {
            vector_type m = vector_type::Zero(n);
            for (size_t i = 0; i < n; i++) if (in(level[i])) m(i) = 1.0;
            return m;
        };
        const vector_type ones = vector_type::Ones(n);
        std::vector<local_operator> ops;
        ops.reserve(L - 1);
        for (size_t l = 0; l + 1 < L; l++) {
            vector_type Q = mask([l](int v) { return v == int(l); });
            ops.emplace_back(sys, ones - Q, Q, R);
            std::cout << "      level " << l << ": cells " << ops.back().qc.size() << ", regions of B y, B^2 y, B^3 y: "
                      << ops.back().region_cells(0) << ", " << ops.back().region_cells(1) << ", " << ops.back().region_cells(2)
                      << " / " << sys.graph().c_ofs.size() << " cells" << std::endl;
        }
        vector_type P = mask([L](int v) { return v == int(L) - 1; });
        local_operator fine(sys, P, ones - P, R);
        std::cout << "      level " << L - 1 << ": cells " << fine.pc.size() << " (fine substeps)" << std::endl;
        if (localize_check()) {
            const vector_type yr = pb.random_state();
            for (size_t l = 0; l < ops.size(); l++) ops[l].check(yr, R, "level " + std::to_string(l));
            fine.check(yr, R, "finest level");
        }
        std::vector<local_operator*> coarse;
        for (auto& op : ops) coarse.push_back(&op);
        const vector_type F0 = vector_type::Zero(sys.N);
        y = local_time_stepping::integrate_mlts(coarse, fine, y, nt, tf, ratios, [&](RealType) { return F0; },
                                                {}, prm.time.mlts_refresh == "substep");
        print_local_stats(sys);
        return y;
    }

    // faces.projected: two-level LTS-RK4 of the splitting B = A + (B - A), A = B~ P the fine operator whose algebraic faces
    // are solved on the faces of the fine region and of faces.projected_rings rings of cells around it only (zero beyond:
    // a local problem), B - A in the coarse Taylor polynomial through local_time_stepping::coarse_coefficients_fast,
    // w_i = a_{i+1} - A a_i - P g_i with a_{i+1} = B a_i + g_i (global, warm-started faces): the splitting is exact, the
    // fine substeps local. With faces.projected_rings large enough to hold all the faces reached, A = B P and the scheme
    // is the LTS-RK4 of integrate_lts.
    struct projected_operator {
        static constexpr local_time_stepping::force_mask  mask  = local_time_stepping::force_mask::after_operator;
        static constexpr local_time_stepping::stage_slope slope = local_time_stepping::stage_slope::from_solve;

        second_order_system& sys;
        const vector_type& Py;
        vector_type Qp, zero, ys;
        std::vector<vector_type> ws;
        bool fine_only, zero_F;
        std::vector<int> pc, pd;                     // cells and interface blocks of the fine region
        std::vector<char> cap;                       // algebraic faces solved by A
        size_t cap_faces = 0;
        std::unique_ptr<second_order_system::local_block> blk;   // A: faces of cap factorized once

        projected_operator(problem& pb, size_t rings, bool zF) : sys(pb.system()), Py(pb.Py),
            Qp(vector_type::Ones(sys.state_size()) - pb.Py), zero(vector_type::Zero(sys.N)), ys(vector_type::Zero(sys.state_size())),
            ws(local_time_stepping::n_sites), fine_only(pb.Py.minCoeff() > 0.5), zero_F(zF) {
            const std::vector<char> fin = sys.mask_cells(Py);
            pc = second_order_system::list_of(fin);
            pd = sys.mask_dblocks(Py);
            cap = sys.cells_cap(sys.grow_cells(fin, rings));
            for (char c : cap) cap_faces += c ? 1 : 0;
            blk = sys.make_local_block(pc, pd, cap);
        }

        vector_type A(const vector_type& y) { return sys.apply_local_block(*blk, y); }
        vector_type f(const vector_type& F, const vector_type& y, size_t site) {
            if (site == local_time_stepping::fine_slope_f) {
                vector_type k = A(y);
                if (!zero_F) k += sys.eval(F, ys, ws[site]);
                return k;
            }
            return sys.eval(F, y, ws[site]);
        }
        vector_type B(const vector_type& y, size_t site) {
            if (site >= local_time_stepping::w_0 && site <= local_time_stepping::w_3) return A(y);
            return sys.eval(zero, y, ws[site]);
        }
        vector_type g(const vector_type& F, size_t site) { return zero_F ? ys : sys.eval(F, ys, ws[site]); }
        vector_type fine(const vector_type& v) const { return Py.cwiseProduct(v); }
        vector_type coarse(const vector_type& v) const { return Qp.cwiseProduct(v); }
        bool all_fine() const { return fine_only; }
        vector_type zero_state() const { return ys; }
    };

    inline vector_type integrate_lts_projected(problem& pb, const wave_parameters& prm, vector_type y, size_t nt, RealType tf) {
        second_order_system& sys = pb.system();
        sys.build_local();
        const bool zF = zero_force(pb, tf);
        projected_operator op(pb, prm.faces.projected_rings, zF);
        static bool told = false;
        if (!told) {
            std::cout << "      projected LTS: fine cells " << op.pc.size() << ", faces of A " << op.cap_faces << " / "
                      << sys.graph().a_ofs.size() << " algebraic faces (" << prm.faces.projected_rings << " ring(s))"
                      << (zF ? ", F = 0" : "") << std::endl;
            told = true;
        }
        const vector_type F0 = vector_type::Zero(sys.N);
        if (zF)
            y = local_time_stepping::integrate_lts(op, y, nt, tf, pb.p, [&](RealType) { return F0; }, {}, true);
        else
            y = local_time_stepping::integrate_lts(op, y, nt, tf, pb.p, [&](RealType t) { return pb.eval_F(t); }, {}, true);
        return y;
    }

    // faces.projected, multilevel: nested splitting B = A_{L-1} + sum_{l < L-1} (A_l - A_{l+1}), A_0 = B and, for l >= 1,
    // A_l the operator of the unknowns of the levels >= l whose faces are solved on their cells and faces.projected_rings
    // rings around them only. Level l < L-1 freezes (A_l - A_{l+1}) a_i over its window, a_i the derivatives of its own
    // problem z' = A_l z + W(tau) (W: Taylor polynomial of the coarser levels, a_0 = z, a_{i+1} = A_l a_i + W^(i)), adds
    // them to W and hands the sum, re-expanded at the start of each sub-window, to level l+1; the finest level takes RK4
    // steps of z' = A_{L-1} z + W(tau). With two levels, the scheme of integrate_lts_projected. F = 0 only.
    inline vector_type integrate_mlts_projected(problem& pb, const wave_parameters& prm, vector_type y, size_t nt, RealType tf) {
        std::vector<int> level;
        const size_t L = pb.mlts_state_levels(prm, level);
        const std::vector<size_t> ratios = mlts_ratios(prm, L);
        std::cout << "      MLTS (projected condensation): " << L << " levels, ratios";
        for (size_t r : ratios) std::cout << " " << r;
        std::cout << std::endl;
        if (L == 1) return integrate_erk(pb, y, nt, tf, 4);
        if (!zero_force(pb, tf))
            throw std::runtime_error("faces.projected with time.scheme = \"mlts\": F = 0 only (no source, homogeneous Dirichlet data)");
        second_order_system& sys = pb.system();
        sys.build_local();
        struct level_set { std::vector<int> c, d; std::vector<char> cap; size_t faces = 0; std::unique_ptr<second_order_system::local_block> blk; };
        std::vector<level_set> A(L);
        for (size_t l = 1; l < L; l++) {
            vector_type m = vector_type::Zero(Eigen::Index(level.size()));
            for (size_t i = 0; i < level.size(); i++) if (level[i] >= int(l)) m(Eigen::Index(i)) = 1.0;
            const std::vector<char> in = sys.mask_cells(m);
            A[l].c = second_order_system::list_of(in);
            A[l].d = sys.mask_dblocks(m);
            A[l].cap = sys.cells_cap(sys.grow_cells(in, prm.faces.projected_rings));
            for (char c : A[l].cap) A[l].faces += c ? 1 : 0;
            A[l].blk = sys.make_local_block(A[l].c, A[l].d, A[l].cap);
            std::cout << "      levels >= " << l << ": cells " << A[l].c.size() << ", faces of A_" << l << " " << A[l].faces
                      << " / " << sys.graph().a_ofs.size() << std::endl;
        }
        const vector_type F0 = vector_type::Zero(sys.N);
        std::vector<vector_type> ws(4);
        auto apply = [&](size_t l, const vector_type& v, size_t i) -> vector_type {
            return l == 0 ? sys.eval(F0, v, ws[i]) : sys.apply_local_block(*A[l].blk, v);
        };
        const double dt = tf / double(nt);
        std::function<void(size_t, double, const std::vector<vector_type>&)> step = [&](size_t l, double h, const std::vector<vector_type>& w_in) {
            if (l == L - 1) {
                y += runge_kutta::rk4_increment(y, h, [&](const vector_type& z, int st) {
                    vector_type k = apply(L - 1, z, 0);
                    k += runge_kutta::taylor_polynomial(w_in, st == 0 ? 0.0 : st < 3 ? 0.5 * h : h);
                    return k;
                });
                return;
            }
            std::vector<vector_type> w = w_in;
            vector_type a = y;
            for (size_t i = 0; i < 4; i++) {
                vector_type Aa = apply(l, a, i);
                w[i] += Aa - apply(l + 1, a, i);
                if (i < 3) a = Aa + w_in[i];
            }
            const double hs = h / double(ratios[l]);
            for (size_t m = 0; m < ratios[l]; m++)
                step(l + 1, hs, local_time_stepping::taylor_shift(w, double(m) * hs));
        };
        for (size_t it = 0; it < nt; it++)
            step(0, dt, std::vector<vector_type>(4, vector_type::Zero(sys.state_size())));
        print_local_stats(sys);
        return y;
    }

    // Newmark(beta, gamma) on M u'' + C u' + K u = F (cells and faces, M zero on the faces): predictor u*, v*, solve
    // (M + gamma dt C + beta dt^2 K) a = F(t_{n+1}) - K u* - C v*, corrector. Discrete energy (69) of [JSC21]
    // (problem::energy, Dirichlet faces = data) written to energy_os at each step if given (exactly conserved for
    // beta = 1/4, gamma = 1/2 without sources). Returns the state y = (u_T, v_T, d) at the final time.
    inline vector_type integrate_newmark(problem& pb, const wave_parameters& prm, size_t nt, std::ostream* energy_os = nullptr) {
        assembler_type& as = pb.assembler();
        const RealType tf = prm.time.final_time, dt = tf / nt;
        const RealType beta = prm.time.newmark_beta, gam = prm.time.newmark_gamma;
        const sparse_type& K = pb.K;
        const sparse_type& M = as.MASS;
        const sparse_type& C = as.COUPLING;
        sparse_type A = M + gam * dt * C + beta * dt * dt * K;
        linear_solver<RealType> analysis;
        analysis.set_fast(prm.solver.fast);     // fast mode: parallel condensation, Pardiso settings of the fast mode
        if (prm.solver.static_condensation) {
            analysis.set_Kg(A, as.get_n_face_dof());
            std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(2);
            vec_cell_basis_data[0] = std::make_pair(as.get_e_material_data().size(), as.get_e_cell_basis_data());
            vec_cell_basis_data[1] = std::make_pair(as.get_a_material_data().size(), as.get_a_cell_basis_data());
            analysis.condense_equations(vec_cell_basis_data);
        }
        else
            analysis.set_Kg(A);
        analysis.factorize();

        auto x0 = pb.newmark_initial_state();
        vector_type u = x0[0], v = x0[1], a = x0[2];
        if (energy_os)
            *energy_os << std::setprecision(16) << std::scientific << 0.0 << " " << pb.energy(u, v, 0.0) << std::endl;
        for (size_t it = 0; it < nt; it++) {
            vector_type F = pb.eval_F_unknowns(dt * (it + 1));
            u += dt * v + 0.5 * dt * dt * (1.0 - 2.0 * beta) * a;
            v += dt * (1.0 - gam) * a;
            vector_type rhs = F - K * u - C * v;
            a = analysis.solve(rhs);
            u += beta * dt * dt * a;
            v += gam * dt * a;
            if (energy_os)
                *energy_os << dt * (it + 1) << " " << pb.energy(u, v, dt * (it + 1)) << std::endl;
        }
        return pb.state(u, v);
    }

    inline vector_type integrate(problem& pb, const wave_parameters& prm, vector_type y, size_t nt) {
        RealType tf = prm.time.final_time;
        if (prm.time.scheme == "newmark") return integrate_newmark(pb, prm, nt);
        if (prm.time.scheme == "erk")   return integrate_erk(pb, y, nt, tf, prm.stages());
        if (prm.time.scheme == "sdirk") return integrate_sdirk(pb, y, nt, tf, prm.stages());
        if (prm.faces.projected && prm.time.scheme == "lts")  return integrate_lts_projected(pb, prm, y, nt, tf);
        if (prm.faces.projected && prm.time.scheme == "mlts") return integrate_mlts_projected(pb, prm, y, nt, tf);
        if (prm.faces.localize && (prm.time.scheme == "lts" || prm.time.scheme == "mlts")) {
            if (zero_force(pb, tf) || localize_check())
                return prm.time.scheme == "lts" ? integrate_lts_local(pb, prm, y, nt, tf) : integrate_mlts_local(pb, prm, y, nt, tf);
            std::cout << "      faces.localize: F is not zero, global evaluations" << std::endl;
        }
        if (prm.time.scheme == "lts")   return integrate_lts(pb, y, nt, tf);
        if (prm.time.scheme == "mlts")  return integrate_mlts(pb, prm, y, nt, tf);
        return integrate_li(pb, y, nt, tf);
    }

}

#endif /* second_order_time_schemes_hpp */
