/*
 *       /\        Guillaume Delay 2018,2019
 *      /__\       guillaume.delay@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is ProtoN, a library for fast Prototyping of
 *  /_\/_\/_\/_\   Numerical methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <cmath>
#include <complex>
#include <Eigen/Dense>

#include "diskpp/mesh/point.hpp"

/* Level set function of the interface, in dimension DIM: the negative side
 * is Omega_1, the positive side Omega_2. */
template<typename T, size_t DIM = 2>
struct level_set {
    static constexpr size_t dimension = DIM;

    virtual ~level_set() = default;

    virtual T operator()(const disk::point<T,DIM>& pt) const = 0;

    virtual Eigen::Matrix<T,DIM,1> gradient(const disk::point<T,DIM>& pt) const = 0;

    Eigen::Matrix<T,DIM,1> normal(const disk::point<T,DIM>& pt) const {
        Eigen::Matrix<T,DIM,1> ret;
        ret = gradient(pt);
        return ret/ret.norm();
    }
};

/* |x - c|^2 - R^2: circle in 2D, sphere in 3D */
template<typename T, size_t DIM = 2>
struct circle_level_set: public level_set<T, DIM>
{
    T                radius;
    disk::point<T,DIM> center;

    circle_level_set(T r, const disk::point<T,DIM>& c) : radius(r), center(c) {
    }

    /* circle of center (a, b) */
    circle_level_set(T r, T a, T b) requires (DIM == 2) : radius(r), center(a, b) {
    }

    T operator()(const disk::point<T,DIM>& pt) const {
        T r2 = 0.0;
        for (size_t i = 0; i < DIM; i++)
            r2 += (pt[i] - center[i]) * (pt[i] - center[i]);
        return r2 - radius*radius;
    }

    Eigen::Matrix<T,DIM,1> gradient(const disk::point<T,DIM>& pt) const {
        Eigen::Matrix<T,DIM,1> ret;
        for (size_t i = 0; i < DIM; i++)
            ret(i) = 2*(pt[i] - center[i]);
        return ret;
    }
};

template<typename T>
struct line_level_set: public level_set<T>
{
    T cut_y;

    line_level_set(T cy)
        : cut_y(cy)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto y = pt.y();

        return y - cut_y;
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = 0;
        ret(1) = 1;
        return ret;
    }
};



template<typename T>
struct square_level_set: public level_set<T>
{
    public:
    T y_top, y_bot, x_left, x_right;

    square_level_set(T yt, T yb, T xl, T xr)
        : y_top(yt), y_bot(yb), x_left(xl), x_right(xr)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        T in = 1;
        if(x > x_left && x < x_right && y > y_bot && y < y_top)
            in = 1;
        else
            in = -1;

        T dist_x = std::min( abs(x-x_left), abs(x-x_right));
        T dist_y = std::min( abs(y-y_bot), abs(y-y_top));


        return - in * std::min(dist_x , dist_y);
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;


        auto x = pt.x();
        auto y = pt.y();

        T dist = abs(x - x_left);
        ret(0) = -1;
        ret(1) = 0;

        if(abs(x - x_right) < dist )
        {
            dist = abs(x - x_right);
            ret(0) = 1;
            ret(1) = 0;
        }
        if(abs(y - y_bot) < dist )
        {
            dist = abs(y - y_bot);
            ret(0) = 0;
            ret(1) = -1;
        }
        if(abs(y - y_top) < dist)
        {
            ret(0) = 0;
            ret(1) = 1;
        }

        return ret;
    }
};



template<typename T>
struct flower_level_set: public level_set<T>
{
    T radius, alpha, beta, a;
    size_t N;

    flower_level_set(T r, T al, T b, size_t N_, T a_)
        : radius(r), alpha(al), beta(b), N(N_), a(a_)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        T theta;
        if(x == alpha && y < beta)
            theta = - M_PI / 2.0;
        else if(x == alpha && y >= beta)
            theta = M_PI / 2.0;
        else
            theta = atan((y-beta)/(x-alpha));

        if(x < alpha)
            theta = theta + M_PI;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - radius*radius
            - a * std::cos(N*theta);
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        auto X = pt.x() - alpha;
        auto Y = pt.y() - beta;

        T theta;
        if(X == 0 && Y < 0)
            theta = - M_PI / 2.0;
        else if(X == 0 && Y >= 0)
            theta = M_PI / 2.0;
        else
            theta = atan( Y / X );

        if(pt.x() < alpha)
            theta = theta + M_PI;

        ret(0) = 2*X - a * N * std::sin(N * theta) * Y / (X*X + Y*Y);
        ret(1) = 2*Y + a * N * std::sin(N * theta) * X / (X*X + Y*Y);
        return ret;
    }
};

/* Coordinates y = Q^T (x - c) of a point in a frame of center c and of
 * orthonormal axes the columns of Q */
template<typename T, size_t DIM, typename Rotation>
Eigen::Matrix<T,DIM,1>
frame_coordinates(const disk::point<T,DIM>& pt, const disk::point<T,DIM>& c, const Rotation& Q)
{
    Eigen::Matrix<T,DIM,1> x;
    for (size_t i = 0; i < DIM; i++)
        x(i) = pt[i] - c[i];
    return Q.transpose() * x;
}

/* Rotation of angle "a" in 2D; rotations of angles a, b, c around the axes
 * z, y, x (Rz(a) Ry(b) Rx(c)) in 3D */
template<typename T, size_t DIM>
Eigen::Matrix<T,DIM,DIM>
rotation_matrix(T a, T b = 0, T c = 0)
{
    Eigen::Matrix<T,DIM,DIM> Q;
    if constexpr (DIM == 2) {
        Q << std::cos(a), -std::sin(a),
             std::sin(a),  std::cos(a);
    }
    else {
        Eigen::Matrix<T,3,3> Rz, Ry, Rx;
        Rz << std::cos(a), -std::sin(a), 0, std::sin(a), std::cos(a), 0, 0, 0, 1;
        Ry << std::cos(b), 0, std::sin(b), 0, 1, 0, -std::sin(b), 0, std::cos(b);
        Rx << 1, 0, 0, 0, std::cos(c), -std::sin(c), 0, std::sin(c), std::cos(c);
        Q = Rz * Ry * Rx;
    }
    return Q;
}

/* Superellipse (superellipsoid in 3D) of semi-axes a_d along the columns of
 * Q, centered at c: s (sum_d (y_d/a_d)^p - 1), y = Q^T (x - c), p even. p = 2
 * gives an ellipse (ellipsoid), p = 4 a rounded rectangle (box). The scale s
 * = 2 mean(a)^2/p makes |grad phi| of the order of that of the circle. */
template<typename T, size_t DIM = 2>
struct superellipse_level_set: public level_set<T, DIM>
{
    disk::point<T,DIM>          center;
    Eigen::Matrix<T,DIM,1>      axes;
    Eigen::Matrix<T,DIM,DIM>    Q;
    size_t                      p;
    T                           scale;

    superellipse_level_set(const disk::point<T,DIM>& c, const Eigen::Matrix<T,DIM,1>& a,
                           const Eigen::Matrix<T,DIM,DIM>& rot, size_t exponent)
        : center(c), axes(a), Q(rot), p(exponent)
    {
        T am = a.mean();
        scale = 2 * am * am / p;
    }

    T operator()(const disk::point<T,DIM>& pt) const {
        auto y = frame_coordinates(pt, center, Q);
        T s = 0;
        for (size_t d = 0; d < DIM; d++)
            s += std::pow(y(d) / axes(d), T(p));
        return scale * (s - 1);
    }

    Eigen::Matrix<T,DIM,1> gradient(const disk::point<T,DIM>& pt) const {
        auto y = frame_coordinates(pt, center, Q);
        Eigen::Matrix<T,DIM,1> g;
        for (size_t d = 0; d < DIM; d++)
            g(d) = scale * p * std::pow(y(d) / axes(d), T(p - 1)) / axes(d);
        return Q * g;
    }
};

/* (sqrt(y_1^2 + y_2^2) - R)^2 + sum_{d >= 3} y_d^2 - r^2, y = Q^T (x - c): the
 * annulus R - r < |x - c| < R + r in 2D (interface made of two circles), a
 * torus of axis Q e_3 in 3D */
template<typename T, size_t DIM = 2>
struct torus_level_set: public level_set<T, DIM>
{
    disk::point<T,DIM>          center;
    T                           R, r;
    Eigen::Matrix<T,DIM,DIM>    Q;

    torus_level_set(const disk::point<T,DIM>& c, T R_, T r_, const Eigen::Matrix<T,DIM,DIM>& rot)
        : center(c), R(R_), r(r_), Q(rot)
    {}

    T operator()(const disk::point<T,DIM>& pt) const {
        auto y = frame_coordinates(pt, center, Q);
        T rho = std::sqrt(y(0)*y(0) + y(1)*y(1));
        T s = (rho - R) * (rho - R) - r * r;
        for (size_t d = 2; d < DIM; d++)
            s += y(d) * y(d);
        return s;
    }

    Eigen::Matrix<T,DIM,1> gradient(const disk::point<T,DIM>& pt) const {
        auto y = frame_coordinates(pt, center, Q);
        T rho = std::sqrt(y(0)*y(0) + y(1)*y(1));
        Eigen::Matrix<T,DIM,1> g = Eigen::Matrix<T,DIM,1>::Zero();
        if (rho > 0) {
            g(0) = 2 * (rho - R) * y(0) / rho;
            g(1) = 2 * (rho - R) * y(1) / rho;
        }
        for (size_t d = 2; d < DIM; d++)
            g(d) = 2 * y(d);
        return Q * g;
    }
};

/* |x - c|^2 - R^2 + a Re(w^n), w = ((x_1 - c_1) + i (x_2 - c_2)) / |x - c|:
 * the flower (55b) in 2D (Re(w^n) = cos(n theta)), in 3D a sphere with n
 * petals around the axis x_3 (Re(w^n) = sin(psi)^n cos(n theta)) */
template<typename T, size_t DIM = 2>
struct star_level_set: public level_set<T, DIM>
{
    disk::point<T,DIM>  center;
    T                   radius, a;
    size_t              n;

    star_level_set(const disk::point<T,DIM>& c, T R, T a_, size_t n_)
        : center(c), radius(R), a(a_), n(n_)
    {}

    T operator()(const disk::point<T,DIM>& pt) const {
        Eigen::Matrix<T,DIM,1> y;
        for (size_t i = 0; i < DIM; i++)
            y(i) = pt[i] - center[i];
        T rho = y.norm();
        T ret = rho * rho - radius * radius;
        if (rho > 0)
            ret += a * std::real(std::pow(std::complex<T>(y(0), y(1)) / rho, int(n)));
        return ret;
    }

    Eigen::Matrix<T,DIM,1> gradient(const disk::point<T,DIM>& pt) const {
        Eigen::Matrix<T,DIM,1> y;
        for (size_t i = 0; i < DIM; i++)
            y(i) = pt[i] - center[i];
        T rho = y.norm();
        Eigen::Matrix<T,DIM,1> g = 2 * y;
        if (rho > 0) {
            /* P = Re(z^n), z = y_1 + i y_2: grad P = n (Re z^(n-1), -Im z^(n-1)),
             * grad (P/rho^n) = grad P / rho^n - n P y / rho^(n+2) */
            std::complex<T> z(y(0), y(1));
            std::complex<T> zn1 = std::pow(z, int(n) - 1);
            T P = std::real(zn1 * z);
            T rn = std::pow(rho, T(n));
            Eigen::Matrix<T,DIM,1> gP = Eigen::Matrix<T,DIM,1>::Zero();
            gP(0) = n * std::real(zn1);
            gP(1) = -T(n) * std::imag(zn1);
            g += a * (gP / rn - n * P * y / (rn * rho * rho));
        }
        return g;
    }
};
