//
//  lts.hpp
//
//  Two-level local time stepping shared by the two formulations (first_order/time_schemes.hpp,
//  second_order/time_schemes.hpp) and by the LTS routines of erk_coupling_hho_scheme / dirk_coupling_hho_scheme used
//  by the prototypes. For y' = B y + g(F(t)), P the 0/1 projector on the fine region and Q = I - P:
//    coarse  w_i = B Q a_i + m_i, i = 0..3, a = (y, B y + g0, B^2 y + B g0 + g1, B^3 y + B^2 g0 + B g1 + g2),
//            g0 + tau g1 + tau^2/2 g2 = g(F(t_n + tau)) (quadratic in time), m_i: forcing of the coarse region
//            (m_3 = 0); frozen over the step as the Taylor polynomial w(tau) = w0 + tau w1 + tau^2/2 w2 + tau^3/6 w3
//    LTS     p RK4 substeps of z' = B P z + (forcing of the fine region) + w(tau)
//    LI      one SDIRK step of the same equation, stages (I - h B P) Z = R + h (forcing of the fine region + w(tau))
//  Two conventions, given by the operator (Op::mask, Op::slope):
//    force_mask::before_operator    m_i = g(Q F_i), fine forcing g(P F): the faces respond to the masked force
//                                   (first-order formulation)
//    force_mask::after_operator     m_i = Q g_i, fine forcing P g(F) (second-order formulation)
//    stage_slope::through_operator  LI slope k = B P Z + (fine forcing) + w(tau) (first-order formulation)
//    stage_slope::from_solve        LI slope k = (Z - R) / h (second-order formulation)
//  The two masks agree in the interior of each region and differ near the fine/coarse boundary when there is a force;
//  the two LI slopes are algebraically equal (Z = R + h k) and differ by round-off only.
//
//  Operator interface (V: state and force vectors):
//    V f(const V& F, const V& y, size_t site)    f(F, y) = B y + g(F)
//    V B(const V& y, size_t site)                B y
//    V g(const V& F, size_t site)                g(F) = f(F, 0)
//    V fine(const V& v), V coarse(const V& v)    P v, Q v on the state
//    V fine_force(const V& F), V coarse_force(const V& F)     P, Q on the force (force_mask::before_operator)
//    bool all_fine()                             Q = 0 (LTS: no coarse coefficients)
//    V zero_state()
//  site: the evaluation (one warm start per site for the iterative operators, the algebraic faces of the second order)
//

#pragma once
#ifndef lts_hpp
#define lts_hpp

#include <Eigen/Dense>
#include <deque>
#include <functional>
#include <type_traits>
#include <vector>

#include "runge_kutta.hpp"

namespace local_time_stepping {

    enum class force_mask { before_operator, after_operator };
    enum class stage_slope { through_operator, from_solve };

    // Evaluation sites (a_4, m_0..m_2: fast mode)
    enum site : size_t { g_n, g_h, g_e, B_y, B2_y, B3_y, B_g0, B2_g0, B_g1, w_0, w_1, w_2, w_3,
                         fine_slope_f, fine_forcing, a_4, m_0, m_1, m_2, n_sites };

    // Coarse Taylor coefficients w (4 vectors) at y, force F on [t_n, t_n + dt] from its values Fn, Fh, Fe at t_n,
    // t_n + dt/2, t_n + dt. The masked quantity is quadratic in time: F (before_operator, w_i = f(Q F_i, Q arg_i)) or
    // g (after_operator, w_i = B Q arg_i + Q g_i)
    template<typename Op, typename V>
    void coarse_coefficients(Op& op, const V& y, const V& Fn, const V& Fh, const V& Fe, double dt, std::vector<V>& w) {
        V F1, F2, g0, g1, g2;
        if constexpr (Op::mask == force_mask::before_operator) {
            runge_kutta::quadratic_in_time(Fn, Fh, Fe, dt, F1, F2);
            g0 = op.g(Fn, g_n);
            g1 = op.g(F1, g_h);
            g2 = op.g(F2, g_e);
        }
        else {
            g0 = op.g(Fn, g_n);
            V gh = op.g(Fh, g_h), ge = op.g(Fe, g_e);
            runge_kutta::quadratic_in_time(g0, gh, ge, dt, g1, g2);
        }
        V By = op.B(y, B_y), B2y = op.B(By, B2_y), B3y = op.B(B2y, B3_y);
        V Bg0 = op.B(g0, B_g0), B2g0 = op.B(Bg0, B2_g0), Bg1 = op.B(g1, B_g1);
        const V arg[4] = {y, By + g0, B2y + Bg0 + g1, B3y + B2g0 + Bg1 + g2};
        [[maybe_unused]] const V* F[3] = {&Fn, &F1, &F2};
        [[maybe_unused]] const V* g[3] = {&g0, &g1, &g2};
        w.resize(4);
        for (size_t i = 0; i < 4; i++) {
            V Qa = op.coarse(arg[i]);
            if (i == 3)
                w[i] = op.B(Qa, w_0 + i);
            else if constexpr (Op::mask == force_mask::before_operator)
                w[i] = op.f(op.coarse_force(*F[i]), Qa, w_0 + i);
            else
                w[i] = op.B(Qa, w_0 + i) + op.coarse(*g[i]);
        }
    }

    // Fast mode (solver.fast) of coarse_coefficients, through a_0 = y, a_{i+1} = B a_i + g_i (g_3 = 0):
    // w_i = B Q a_i + m_i = a_{i+1} - B P a_i - (forcing of the fine region at F_i), 11 evaluations instead of 13 (equal
    // up to round-off)
    template<typename Op, typename V>
    void coarse_coefficients_fast(Op& op, const V& y, const V& Fn, const V& Fh, const V& Fe, double dt, std::vector<V>& w) {
        V F1, F2, g0, g1, g2;
        if constexpr (Op::mask == force_mask::before_operator) {
            runge_kutta::quadratic_in_time(Fn, Fh, Fe, dt, F1, F2);
            g0 = op.g(Fn, g_n);
            g1 = op.g(F1, g_h);
            g2 = op.g(F2, g_e);
        }
        else {
            g0 = op.g(Fn, g_n);
            V gh = op.g(Fh, g_h), ge = op.g(Fe, g_e);
            runge_kutta::quadratic_in_time(g0, gh, ge, dt, g1, g2);
        }
        V a[5];
        a[0] = y;
        a[1] = op.B(y, B_y) + g0;
        a[2] = op.B(a[1], B2_y) + g1;
        a[3] = op.B(a[2], B3_y) + g2;
        a[4] = op.B(a[3], a_4);
        [[maybe_unused]] const V* F[3] = {&Fn, &F1, &F2};
        [[maybe_unused]] const V* g[3] = {&g0, &g1, &g2};
        w.resize(4);
        for (size_t i = 0; i < 4; i++) {
            w[i] = a[i + 1] - op.B(op.fine(a[i]), w_0 + i);
            if (i == 3) continue;
            if constexpr (Op::mask == force_mask::before_operator)
                w[i] -= op.g(op.fine_force(*F[i]), m_0 + i);
            else
                w[i] -= op.fine(*g[i]);
        }
    }

    // B P z + forcing of the fine region at the force F
    template<typename Op, typename V>
    V fine_slope(Op& op, const V& F, const V& z) {
        if constexpr (Op::mask == force_mask::before_operator)
            return op.f(op.fine_force(F), op.fine(z), fine_slope_f);
        else {
            // f(F, P z) - Q g(F) = B P z + P g(F)
            V k = op.f(F, op.fine(z), fine_slope_f);
            if (!op.all_fine())
                k -= op.coarse(op.g(F, fine_forcing));
            return k;
        }
    }

    // Forcing of the fine region at the force F
    template<typename Op, typename V>
    V fine_forcing_term(Op& op, const V& F) {
        if constexpr (Op::mask == force_mask::before_operator)
            return op.g(op.fine_force(F), fine_forcing);
        else
            return op.fine(op.g(F, fine_forcing));
    }

    // One RK4 fine substep [t_n + tm, t_n + tm + dtau] of x, force Fm, Fh, Fe at its start, middle and end
    template<typename Op, typename V>
    void fine_substep(Op& op, V& x, const std::vector<V>& w, const V& Fm, const V& Fh, const V& Fe, double tm, double dtau) {
        double tmh = tm + 0.5 * dtau, tme = tm + dtau;
        x += runge_kutta::rk4_increment(x, dtau, [&](const V& z, int st) {
            V k = fine_slope(op, st == 0 ? Fm : st < 3 ? Fh : Fe, z);
            k += runge_kutta::taylor_polynomial(w, st == 0 ? tm : st < 3 ? tmh : tme);
            return k;
        });
    }

    // One LI step [t_n, t_n + dt] of x: SDIRK (a, b, c) whose stages (I - h B P) Z = rhs are solved by solve(rhs),
    // h = a_ii dt; force(t) = F(t)
    template<typename Op, typename V, typename Force, typename Solve>
    void li_step(Op& op, V& x, const std::vector<V>& w, double tn, double dt, const runge_kutta::tableau_matrix& a,
                 const runge_kutta::tableau_vector& b, const runge_kutta::tableau_vector& c, double h, Force&& force,
                 Solve&& solve) {
        runge_kutta::sdirk(x, 1, dt, a, b, c, [](double, const V&) {}, [&](size_t i, double t, V& R, V& k) {
            V F = force(tn + t);
            V wt = runge_kutta::taylor_polynomial(w, c(i) * dt);
            V rhs = R + h * (fine_forcing_term(op, F) + wt);
            V Z = solve(rhs);
            if constexpr (Op::slope == stage_slope::from_solve)
                k = (Z - R) / h;
            else {
                k = fine_slope(op, F, Z);
                k += wt;
            }
        }, [](double, const V&) { return true; });
    }

    // Fast mode: forces F(t) of the last distinct times and their corrections of the fine slopes (after_operator: Q g(F),
    // computed once per force), reused when a time comes back bit for bit
    template<typename Op, typename V, typename Force>
    struct force_memo {
        struct entry { double t; V F, q; bool q_ready = false; };
        Op& op;
        Force& force;
        std::deque<entry> entries;
        static constexpr size_t capacity = 8;

        force_memo(Op& op, Force& force) : op(op), force(force) {}

        entry& at(double t) {
            for (auto& e : entries)
                if (e.t == t) return e;
            entries.push_back({t, force(t), V(), false});
            if (entries.size() > capacity) entries.pop_front();
            return entries.back();
        }
        const V& F(double t) { return at(t).F; }
        // Q g(F(t)) (after_operator, not all fine), empty otherwise
        const V& correction(double t) {
            entry& e = at(t);
            if constexpr (Op::mask == force_mask::after_operator)
                if (!e.q_ready && !op.all_fine()) e.q = op.coarse(op.g(e.F, fine_forcing));
            e.q_ready = true;
            return e.q;
        }
    };

    // fine_slope with the correction q = Q g(F) given (after_operator)
    template<typename Op, typename V>
    V fine_slope_fast(Op& op, const V& F, const V& q, const V& z) {
        if constexpr (Op::mask == force_mask::before_operator)
            return op.f(op.fine_force(F), op.fine(z), fine_slope_f);
        else {
            V k = op.f(F, op.fine(z), fine_slope_f);
            if (q.size()) k -= q;
            return k;
        }
    }

    // Two-level LTS-RK4 from x, nt steps up to tf, p fine substeps; observe(t, x) after each step
    // Fast mode: coarse_coefficients_fast; forces at the half substeps t_n + j dtau / 2 (t_n + dt / 2 and t_n + dt
    // exactly at j = p and 2p: the forces of the coarse coefficients), each force and its correction computed once
    template<typename Op, typename V, typename Force>
    V integrate_lts(Op& op, V x, size_t nt, double tf, size_t p, Force&& force,
                    const std::function<void(double, const std::type_identity_t<V>&)>& observe = {}, bool fast = false) {
        double dt = tf / nt;
        double dtau = dt / p;
        std::vector<V> w(4);
        if (fast) {
            force_memo<Op, V, std::remove_reference_t<Force>> memo(op, force);
            for (size_t it = 0; it < nt; it++) {
                double tn = dt * it;
                auto tau = [&](size_t j) { return j == 2 * p ? dt : j == p ? 0.5 * dt : j * (0.5 * dtau); };
                if (op.all_fine())
                    for (auto& wi : w) wi = op.zero_state();
                else {
                    V Fn = memo.F(tn), Fh = memo.F(tn + 0.5 * dt), Fe = memo.F(tn + dt);
                    coarse_coefficients_fast(op, x, Fn, Fh, Fe, dt, w);
                }
                for (size_t m = 0; m < p; m++) {
                    const double ts[3] = {tau(2 * m), tau(2 * m + 1), tau(2 * m + 2)};
                    V F[3], q[3];
                    for (size_t j = 0; j < 3; j++) {
                        F[j] = memo.F(tn + ts[j]);
                        q[j] = memo.correction(tn + ts[j]);
                    }
                    x += runge_kutta::rk4_increment(x, dtau, [&](const V& z, int st) {
                        size_t j = st == 0 ? 0 : st < 3 ? 1 : 2;
                        V k = fine_slope_fast(op, F[j], q[j], z);
                        k += runge_kutta::taylor_polynomial(w, ts[j]);
                        return k;
                    });
                }
                if (observe) observe(dt * (it + 1), x);
            }
            return x;
        }
        for (size_t it = 0; it < nt; it++) {
            double tn = dt * it;
            if (op.all_fine())
                for (auto& wi : w) wi = op.zero_state();
            else
                coarse_coefficients(op, x, force(tn), force(tn + 0.5 * dt), force(tn + dt), dt, w);
            for (size_t m = 0; m < p; m++) {
                double tm = m * dtau, tmh = tm + 0.5 * dtau, tme = tm + dtau;
                fine_substep(op, x, w, force(tn + tm), force(tn + tmh), force(tn + tme), tm, dtau);
            }
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

    // Taylor coefficients of w(t) = sum_j t^j / j! w^(j) (degree 3) re-expanded around a:
    // w^(j)(a) = sum_{q >= j} a^(q - j) / (q - j)! w^(q)
    template<typename V>
    std::vector<V> taylor_shift(const std::vector<V>& w, double a) {
        const double a2 = a * a, a3 = a2 * a;
        std::vector<V> out(4);
        out[0] = w[0] + a * w[1] + (a2 / 2.0) * w[2] + (a3 / 6.0) * w[3];
        out[1] = w[1] + a * w[2] + (a2 / 2.0) * w[3];
        out[2] = w[2] + a * w[3];
        out[3] = w[3];
        return out;
    }

    // Multilevel LTS-RK4 (Algorithm 4 of the preprint LTS_HHO, recursion of Diaz-Grote 2015 for Runge-Kutta): levels
    // l = 0 (coarsest) .. L-1 (finest), level l advanced with dt_l = dt / (p_0 ... p_{l-1}). coarse[l], l < L-1: operator
    // whose coarse projector selects the unknowns of level l; fine: operator whose fine projector selects those of level
    // L-1. Level l < L-1 computes its own Taylor coefficients (coarse_coefficients) from the current state over its window
    // of length dt_l (refresh_each_substep: over each of its p_l sub-windows instead, the variant of the prototypes), adds
    // them to the inherited ones and hands the sum, re-expanded at the start of each sub-window, to level l+1; level L-1
    // takes one RK4 step of z' = B P z + (forcing of level L-1) + w(tau). With L = 2 it is the two-level LTS-RK4.
    template<typename Op, typename V, typename Force>
    V integrate_mlts(const std::vector<Op*>& coarse, Op& fine, V x, size_t nt, double tf, const std::vector<size_t>& p,
                     Force&& force, const std::function<void(double, const std::type_identity_t<V>&)>& observe = {},
                     bool refresh_each_substep = false) {
        const size_t L = coarse.size() + 1;
        const double dt = tf / nt;
        std::vector<V> own;
        std::function<void(size_t, double, double, const std::vector<V>&)> step =
            [&](size_t l, double t0, double h, const std::vector<V>& w_in) {
            if (l == L - 1) {
                fine_substep(fine, x, w_in, force(t0), force(t0 + 0.5 * h), force(t0 + h), 0.0, h);
                return;
            }
            const double hs = h / double(p[l]);
            auto add_own = [&](double ts, double len, std::vector<V>& acc) {
                coarse_coefficients(*coarse[l], x, force(ts), force(ts + 0.5 * len), force(ts + len), len, own);
                for (size_t j = 0; j < 4; j++) acc[j] += own[j];
            };
            std::vector<V> w = w_in;
            if (!refresh_each_substep) add_own(t0, h, w);
            for (size_t m = 0; m < p[l]; m++) {
                std::vector<V> ws = taylor_shift(w, double(m) * hs);
                if (refresh_each_substep) add_own(t0 + double(m) * hs, hs, ws);
                step(l + 1, t0 + double(m) * hs, hs, ws);
            }
        };
        for (size_t it = 0; it < nt; it++) {
            const double tn = dt * double(it);
            step(0, tn, dt, std::vector<V>(4, fine.zero_state()));
            if (observe) observe(dt * double(it + 1), x);
        }
        return x;
    }

    // LI-SDIRK from x, nt steps up to tf (stage solve: see li_step); observe(t, x) after each step
    // Fast mode: coarse_coefficients_fast, forces of the coarse coefficients reused by the stages at the same times
    template<typename Op, typename V, typename Force, typename Solve>
    V integrate_li(Op& op, V x, size_t nt, double tf, const runge_kutta::tableau_matrix& a,
                   const runge_kutta::tableau_vector& b, const runge_kutta::tableau_vector& c, Force&& force, Solve&& solve,
                   const std::function<void(double, const std::type_identity_t<V>&)>& observe = {}, bool fast = false) {
        double dt = tf / nt;
        double h = a(0,0) * dt;
        std::vector<V> w(4);
        force_memo<Op, V, std::remove_reference_t<Force>> memo(op, force);
        auto memo_force = [&](double t) { return V(memo.F(t)); };
        for (size_t it = 0; it < nt; it++) {
            double tn = dt * it;
            if (fast) {
                V Fn = memo.F(tn), Fh = memo.F(tn + 0.5 * dt), Fe = memo.F(tn + dt);
                coarse_coefficients_fast(op, x, Fn, Fh, Fe, dt, w);
                li_step(op, x, w, tn, dt, a, b, c, h, memo_force, solve);
                if (observe) observe(dt * (it + 1), x);
                continue;
            }
            coarse_coefficients(op, x, force(tn), force(tn + 0.5 * dt), force(tn + dt), dt, w);
            li_step(op, x, w, tn, dt, a, b, c, h, force, solve);
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

}

#endif /* lts_hpp */
