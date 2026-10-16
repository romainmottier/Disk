//
//  runge_kutta.hpp
//
//  Runge-Kutta drivers shared by the two formulations (first_order/time_schemes.hpp, second_order/time_schemes.hpp)
//  and the building blocks of their LTS / LI integrators. A formulation gives its slope (ERK) or its stage solve
//  (SDIRK); the drivers make the stage combinations and the updates, the same operations for both:
//    erk    x_{n+1} = x_n + dt sum_i b_i k_i,  k_i = f(t_n + c_i dt, x_n + dt sum_{j<i} a_ij k_j)
//    sdirk  k_i = f(t_n + c_i dt, r_i + a_ii dt k_i),  r_i = x_n + dt sum_{j<i} a_ij k_j
//  LTS / LI: quadratic interpolation of the force on [t_n, t_n + dt], Taylor polynomial of the coarse coefficients,
//  classical RK4 increment of the fine substeps.
//

#pragma once
#ifndef runge_kutta_hpp
#define runge_kutta_hpp

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace runge_kutta {

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tableau_matrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1>              tableau_vector;

    // ERK (a, b, c) from x, nt steps up to tf: slope(t, y, k) computes k = f(t, y); after_step(t, x) after each step,
    // the integration stops when it returns false. Returns the number of steps done.
    template<typename Vector, typename Slope, typename AfterStep>
    size_t erk(Vector& x, size_t nt, double tf, const tableau_matrix& a, const tableau_vector& b, const tableau_vector& c,
               Slope&& slope, AfterStep&& after_step) {
        const size_t s = b.size();
        double dt = tf / nt;
        std::vector<Vector> k(s);
        for (size_t it = 0; it < nt; it++) {
            double tn = dt * it;
            Vector x_n = x;
            for (size_t i = 0; i < s; i++) {
                Vector y = x_n;
                for (size_t j = 0; j < i; j++)
                    y += a(i,j) * dt * k[j];
                slope(tn + c(i) * dt, y, k[i]);
                x += dt * b(i) * k[i];
            }
            if (!after_step(dt * (it + 1), x)) return it + 1;
        }
        return nt;
    }

    // SDIRK (a, b, c) from x, nt steps up to tf: stage(i, t_i, r_i, k_i) computes k_i = f(t_i, r_i + a_ii dt k_i);
    // before_step(t_n, x) before each step; after_step(t, x) after each step, the integration stops when it returns
    // false. Returns the number of steps done.
    template<typename Vector, typename BeforeStep, typename Stage, typename AfterStep>
    size_t sdirk(Vector& x, size_t nt, double tf, const tableau_matrix& a, const tableau_vector& b, const tableau_vector& c,
                 BeforeStep&& before_step, Stage&& stage, AfterStep&& after_step) {
        const size_t s = b.size();
        double dt = tf / nt;
        std::vector<Vector> k(s);
        for (size_t it = 0; it < nt; it++) {
            double tn = dt * it;
            before_step(tn, x);
            Vector x_n = x;
            for (size_t i = 0; i < s; i++) {
                Vector r = x_n;
                for (size_t j = 0; j < i; j++)
                    r += a(i,j) * dt * k[j];
                stage(i, tn + c(i) * dt, r, k[i]);
                x += dt * b(i) * k[i];
            }
            if (!after_step(dt * (it + 1), x)) return it + 1;
        }
        return nt;
    }

    // Force on [t_n, t_n + dt] from its values at t_n, t_n + dt/2 and t_n + dt: f(t_n + tau) = f0 + tau f1 + tau^2/2 f2,
    // f1 = (-3 f(t_n) + 4 f(t_n + dt/2) - f(t_n + dt)) / dt, f2 = 4 (f(t_n) - 2 f(t_n + dt/2) + f(t_n + dt)) / dt^2
    template<typename Vector>
    void quadratic_in_time(const Vector& fn, const Vector& fh, const Vector& fe, double dt, Vector& f1, Vector& f2) {
        f1 = (-3*fn + 4*fh - fe) / dt;
        f2 = ( 4*fn - 8*fh + 4*fe) / (dt*dt);
    }

    // Weights 1, tau, tau^2/2, tau^3/6 of the Taylor polynomial w0 + tau w1 + tau^2/2 w2 + tau^3/6 w3
    inline std::array<double, 4> taylor_weights(double tau) {
        double tau2 = tau*tau, tau3 = tau*tau2;
        return {1.0, tau, tau2/2.0, tau3/6.0};
    }

    // Taylor polynomial of the coarse coefficients w (4 vectors) at the local time tau
    template<typename W>
    typename W::value_type taylor_polynomial(const W& w, double tau) {
        const std::array<double, 4> c = taylor_weights(tau);
        return w[0] + c[1]*w[1] + c[2]*w[2] + c[3]*w[3];
    }

    // Classical RK4 increment h (k1 + 2 k2 + 2 k3 + k4) / 6 of z' = f from z: k_{s+1} = slope(z_s, s) at the times
    // t, t + h/2, t + h/2, t + h for s = 0, 1, 2, 3
    template<typename Vector, typename Slope>
    Vector rk4_increment(const Vector& z, double h, Slope&& slope) {
        Vector k1 = slope(z, 0);
        Vector k2 = slope(z + 0.5*h*k1, 1);
        Vector k3 = slope(z + 0.5*h*k2, 2);
        Vector k4 = slope(z +     h*k3, 3);
        return h * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    }

}

#endif /* runge_kutta_hpp */
