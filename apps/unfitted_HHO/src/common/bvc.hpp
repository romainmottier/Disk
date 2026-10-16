/*
 * Unfitted HHO method with boundary value correction (prototype -p2) for the
 * elliptic interface problem:
 *
 *   E. Burman, "An unfitted hybrid high-order method with boundary value
 *   correction for elliptic interface problems".
 *
 * Port to diskpp, in 2D and 3D, of the ProtoN implementation cuthho_bvc.cpp
 * (branch megamerge), whose notations are kept.
 *
 * The interface is replaced by a polygonal one, made of flat pieces whose
 * vertices lie on the interface (disk::cut_geometry_options::flat: one chord
 * per cut cell in 2D; in 3D the chords of the cut faces and a triangulation
 * of their loop, see disk::cutgeo::build_flat_cut_polyhedron), and every
 * integral is computed on straight subcells, faces and pieces of the
 * interface. The geometric error is compensated at the quadrature nodes x of
 * the interface pieces, of normal n_h (from Omega_1 to Omega_2):
 *
 *   rho_h(x): phi(x + rho_h n_h) = 0 (root of smallest modulus),
 *   p_h(x) = x + rho_h(x) n_h on Gamma, n*(x) = grad phi(p_h)/|grad phi(p_h)|,
 *   T_m w(x) = sum_{j <= m} rho_h^j/j! d^j_{n_h} w(x)   (Taylor operator),
 *
 * applied to the trial function only (non-symmetric method). Mixed order:
 * cells P^{k+1}, faces P^k; small cuts are controlled by cell agglomeration
 * (threshold -a). With [v] = v_1 - v_2 and
 *
 *   (G_i v, q)_{T^i} = (grad v_{T^i}, q) + sum_F (v_F - v_{T^i}, q.n_T)_F,
 *   (L g, q)_{T^1}   = (g, q.n_h)_{Gamma_h},
 *   Ghat_1(v) = G_1(v) - L([T_m v_T]),  Gring_1(w) = G_1(w) - L([w_T]),
 *
 * the local forms of a cut cell are
 *
 *   a_T(v,w) = k1 (Ghat_1 v, Gring_1 w)_{T^1} + k2 (G_2 v, G_2 w)_{T^2}
 *            + sum_i k_i/h_T (Pi_F v_{T^i} - v_F, Pi_F w_{T^i} - w_F)_{dT^i}
 *            + eta k1/h_T ([T_m v_T], [w_T])_{Gamma_h} + gamma_T(v,w),
 *   gamma_T(v,w) = (T_m(Q v).n* - (Q v).n_h, w_{T^2})_{Gamma_h},
 *   Q v = k1 G_1 v - k2 G_2 v,
 *   l_T(w) = sum_i (f_i, w_{T^i}) + k1 (gD(p_h), eta h_T^{-1} [w_T] - Gring_1(w).n_h)_{Gamma_h}
 *          + (gN(p_h), w_{T^2})_{Gamma_h},  gN(p_h) = [k grad u](p_h).n*.
 *
 * Variants (-v): full (m = m(k), gamma_T), taylor (m = m(k), no gamma_T),
 * gam (m = 0, gamma_T), naive (m = 0, no gamma_T: the method of [BCDE21] on
 * the polygonal interface), ext (naive with the data of the extended solution
 * at x with n_h: consistency benchmark). m(k) = 0, 1, 1, 2, 2 for
 * k = 0, ..., 4 (-T overrides it).
 *
 * Test cases (-t): 0 contrast u_i = rho^6/kappa_i (+ const) around the circle
 * or the sphere (-i 0), 1 discontinuous solution u_1 = prod cos(pi x_d),
 * u_2 = prod sin(pi x_d) (any interface; -i 7: flower of the paper),
 * 2 patch test: plane interface (-i 8) and a polynomial of degree k+1 on each
 * side, reproduced up to round-off.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#ifdef HAVE_INTEL_MKL
#include <Eigen/PardisoSupport>
#include <unsupported/Eigen/IterativeSolvers>
#endif

#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/bases/cut_bases.hpp"

#ifdef HAVE_INTEL_MKL
/* Pardiso LU for a structurally symmetric matrix (mtype 1): symmetric
 * ordering, no matching nor scaling, cheaper than the general unsymmetric
 * factorization (mtype 11) of Eigen::PardisoLU. The matrices of the method
 * have a symmetric pattern (dense local Schur complements). */
template<typename MatrixType> class pardiso_structsym_lu;
namespace Eigen { namespace internal {
template<typename MT>
struct pardiso_traits<pardiso_structsym_lu<MT>> {
    typedef MT MatrixType;
    typedef typename MT::Scalar Scalar;
    typedef typename MT::RealScalar RealScalar;
    typedef typename MT::StorageIndex StorageIndex;
};
} }
template<typename MatrixType>
class pardiso_structsym_lu : public Eigen::PardisoImpl<pardiso_structsym_lu<MatrixType>> {
protected:
    using Base = Eigen::PardisoImpl<pardiso_structsym_lu<MatrixType>>;
    using Base::pardisoInit;
    using Base::m_matrix;
    friend class Eigen::PardisoImpl<pardiso_structsym_lu<MatrixType>>;
public:
    using Base::compute;
    using Base::solve;
    pardiso_structsym_lu() : Base() { pardisoInit(1); }
protected:
    void getMatrix(const MatrixType& matrix) {
        m_matrix = matrix;
        m_matrix.makeCompressed();
    }
};
#endif

#include "preprocessor.hpp"
#include "agglo_operators.hpp"
#include "operators.hpp"

namespace bvc {

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T, size_t DIM>
using vec = Eigen::Matrix<T, DIM, 1>;

inline disk::location side_location(size_t s) {
    return s == 0 ? disk::location::IN_NEGATIVE_SIDE : disk::location::IN_POSITIVE_SIDE;
}

/* dimension of P^deg in DIM variables */
inline size_t dim_pk(size_t deg, size_t DIM) {
    size_t n = 1;
    for (size_t i = 1; i <= DIM; i++)
        n = n * (deg + i) / i;
    return n;
}

/* exponents of the monomials of degree <= deg, ordered by total degree */
template<size_t DIM>
std::vector<std::array<int,DIM>>
priv_exponents(size_t deg)
{
    std::vector<std::array<int,DIM>> ret;
    for (int d = 0; d <= int(deg); d++) {
        if constexpr (DIM == 2) {
            for (int a = d; a >= 0; a--)
                ret.push_back({a, d - a});
        }
        else {
            for (int a = d; a >= 0; a--)
                for (int b = d - a; b >= 0; b--)
                    ret.push_back({a, b, d - a - b});
        }
    }
    return ret;
}

/***************************************************************************/
/*                               TEST PROBLEMS                              */
/***************************************************************************/

/* The solution, its gradient and the load are given on both sides (0 =
 * Omega_1 = {phi < 0}, 1 = Omega_2 = {phi > 0}) by expressions valid around
 * the interface: the method evaluates them on the polygonal subcells, which
 * stick out of Omega_i (Assumption 2.1 of the paper). */
template<typename T, size_t DIM>
struct problem {
    using point_type = disk::point<T,DIM>;

    std::string name;
    T kappa[2] = {1.0, 1.0};
    std::shared_ptr<level_set<T,DIM>> ls;
    std::function<T(const point_type&, size_t)>             u;
    std::function<vec<T,DIM>(const point_type&, size_t)>    grad_u;
    std::function<T(const point_type&, size_t)>             f;

    T gD(const point_type& x) const { return u(x, 0) - u(x, 1); }

    T gN(const point_type& x, const vec<T,DIM>& n) const {
        return (kappa[0] * grad_u(x, 0) - kappa[1] * grad_u(x, 1)).dot(n);
    }
};

/* u_i = rho^6/kappa_i, + R^6 (1/kappa_1 - 1/kappa_2) on side 2: no jumps,
 * f = -Delta rho^6 = -6 (DIM+4) rho^4 ([BCDE21, (36)]) */
template<typename T, size_t DIM>
problem<T,DIM>
make_contrast_problem(std::shared_ptr<level_set<T,DIM>> ls, T kap1, T kap2)
{
    auto circ = std::dynamic_pointer_cast<circle_level_set<T,DIM>>(ls);
    if (!circ)
        throw std::invalid_argument("-p2 -t0 needs the circle or the sphere (-i 0)");
    T R = circ->radius;
    auto c = circ->center;
    T off = std::pow(R, 6) * (1.0/kap1 - 1.0/kap2);

    problem<T,DIM> p;
    std::ostringstream ss;
    ss << "contrast rho^6 (kappa2/kappa1 = " << kap2/kap1 << ")";
    p.name = ss.str();
    p.kappa[0] = kap1;
    p.kappa[1] = kap2;
    p.ls = ls;
    auto r2 = [=](const disk::point<T,DIM>& x) {
        T s = 0;
        for (size_t d = 0; d < DIM; d++)
            s += (x[d] - c[d]) * (x[d] - c[d]);
        return s;
    };
    p.u = [=](const disk::point<T,DIM>& x, size_t s) {
        T r = r2(x);
        return r*r*r / (s == 0 ? kap1 : kap2) + (s == 1 ? off : 0.0);
    };
    p.grad_u = [=](const disk::point<T,DIM>& x, size_t s) {
        T r = r2(x);
        vec<T,DIM> g;
        for (size_t d = 0; d < DIM; d++)
            g(d) = 6.0 * r * r * (x[d] - c[d]) / (s == 0 ? kap1 : kap2);
        return g;
    };
    p.f = [=](const disk::point<T,DIM>& x, size_t) {
        T r = r2(x);
        return -6.0 * (DIM + 4.0) * r * r;
    };
    return p;
}

/* u_1 = prod cos(pi x_d), u_2 = prod sin(pi x_d): discontinuous solution
 * and flux, f_i = kappa_i DIM pi^2 u_i ([BCDE21, (37)]) */
template<typename T, size_t DIM>
problem<T,DIM>
make_discontinuous_problem(std::shared_ptr<level_set<T,DIM>> ls, T kap1, T kap2)
{
    problem<T,DIM> p;
    std::ostringstream ss;
    ss << "discontinuous prod cos / prod sin (kappa2/kappa1 = " << kap2/kap1 << ")";
    p.name = ss.str();
    p.kappa[0] = kap1;
    p.kappa[1] = kap2;
    p.ls = ls;
    auto fun = [](const disk::point<T,DIM>& x, size_t s) {
        T v = 1;
        for (size_t d = 0; d < DIM; d++)
            v *= (s == 0) ? std::cos(M_PI * x[d]) : std::sin(M_PI * x[d]);
        return v;
    };
    p.u = fun;
    p.grad_u = [](const disk::point<T,DIM>& x, size_t s) {
        vec<T,DIM> g;
        for (size_t d = 0; d < DIM; d++) {
            T v = M_PI * ((s == 0) ? -std::sin(M_PI * x[d]) : std::cos(M_PI * x[d]));
            for (size_t e = 0; e < DIM; e++)
                if (e != d)
                    v *= (s == 0) ? std::cos(M_PI * x[e]) : std::sin(M_PI * x[e]);
            g(d) = v;
        }
        return g;
    };
    p.f = [=](const disk::point<T,DIM>& x, size_t s) {
        return (s == 0 ? kap1 : kap2) * DIM * M_PI * M_PI * fun(x, s);
    };
    return p;
}

/* Eq. (59) of the polynomial extension paper (-p0, -p1 -t4): u_1 = e^x cos y
 * (harmonic), u_2 = prod sin(pi x_d), kappa = 1: discontinuous solution and
 * flux on any interface */
template<typename T, size_t DIM>
problem<T,DIM>
make_jumps59_problem(std::shared_ptr<level_set<T,DIM>> ls)
{
    problem<T,DIM> p;
    p.name = "non-polynomial jumps (59), u_1 = e^x cos y, u_2 = prod sin(pi x_d)";
    p.kappa[0] = 1.0;
    p.kappa[1] = 1.0;
    p.ls = ls;
    auto sp = [](const disk::point<T,DIM>& x) {
        T v = 1;
        for (size_t d = 0; d < DIM; d++)
            v *= std::sin(M_PI * x[d]);
        return v;
    };
    p.u = [=](const disk::point<T,DIM>& x, size_t s) {
        return (s == 0) ? std::exp(x[0]) * std::cos(x[1]) : sp(x);
    };
    p.grad_u = [](const disk::point<T,DIM>& x, size_t s) {
        vec<T,DIM> g = vec<T,DIM>::Zero();
        if (s == 0) {
            g(0) = std::exp(x[0]) * std::cos(x[1]);
            g(1) = -std::exp(x[0]) * std::sin(x[1]);
            return g;
        }
        for (size_t d = 0; d < DIM; d++) {
            T v = M_PI * std::cos(M_PI * x[d]);
            for (size_t e = 0; e < DIM; e++)
                if (e != d)
                    v *= std::sin(M_PI * x[e]);
            g(d) = v;
        }
        return g;
    };
    p.f = [=](const disk::point<T,DIM>& x, size_t s) {
        return (s == 0) ? T(0) : DIM * M_PI * M_PI * sp(x);
    };
    return p;
}

/* Patch test: plane interface and a polynomial of degree deg on each side
 * (deterministic coefficients), kappa = (1, 7) by default: for deg <= k+1 the
 * method reproduces it up to round-off */
template<typename T, size_t DIM>
problem<T,DIM>
make_patch_problem(std::shared_ptr<level_set<T,DIM>> ls, size_t deg, T kap1, T kap2)
{
    if (!std::dynamic_pointer_cast<plane_level_set<T,DIM>>(ls))
        throw std::invalid_argument("-p2 -t2 needs the plane interface (-i 8)");

    problem<T,DIM> p;
    std::ostringstream ss;
    ss << "patch test, polynomials of degree " << deg << " (kappa2/kappa1 = " << kap2/kap1 << ")";
    p.name = ss.str();
    p.kappa[0] = kap1;
    p.kappa[1] = kap2;
    p.ls = ls;

    std::vector<std::array<int,DIM>> ex;
    for (auto& e : priv_exponents<DIM>(deg))
        ex.push_back(e);
    std::vector<T> c0(ex.size()), c1(ex.size());
    for (size_t i = 0; i < ex.size(); i++) {
        c0[i] = std::sin(1.3 + 2.1*i);
        c1[i] = std::cos(0.7 + 1.7*i);
    }
    auto ipow = [](T x, int e) { T r = 1; for (int i = 0; i < e; i++) r *= x; return r; };
    p.u = [=](const disk::point<T,DIM>& x, size_t s) {
        const auto& c = (s == 0) ? c0 : c1;
        T r = 0;
        for (size_t i = 0; i < ex.size(); i++) {
            T m = c[i];
            for (size_t d = 0; d < DIM; d++)
                m *= ipow(x[d], ex[i][d]);
            r += m;
        }
        return r;
    };
    p.grad_u = [=](const disk::point<T,DIM>& x, size_t s) {
        const auto& c = (s == 0) ? c0 : c1;
        vec<T,DIM> g = vec<T,DIM>::Zero();
        for (size_t i = 0; i < ex.size(); i++)
            for (size_t d = 0; d < DIM; d++) {
                if (ex[i][d] == 0)
                    continue;
                T m = c[i] * ex[i][d];
                for (size_t e = 0; e < DIM; e++)
                    m *= ipow(x[e], ex[i][e] - (e == d ? 1 : 0));
                g(d) += m;
            }
        return g;
    };
    p.f = [=](const disk::point<T,DIM>& x, size_t s) {
        const auto& c = (s == 0) ? c0 : c1;
        T lap = 0;
        for (size_t i = 0; i < ex.size(); i++)
            for (size_t d = 0; d < DIM; d++) {
                if (ex[i][d] < 2)
                    continue;
                T m = c[i] * ex[i][d] * (ex[i][d] - 1);
                for (size_t e = 0; e < DIM; e++)
                    m *= ipow(x[e], ex[i][e] - (e == d ? 2 : 0));
                lap += m;
            }
        return -(s == 0 ? kap1 : kap2) * lap;
    };
    return p;
}

/***************************************************************************/
/*                                  BASES                                   */
/***************************************************************************/

/* L2(T^i)-orthonormal hierarchical basis of P^deg(T^i): the scaled
 * monomials ((x - xc)/h)^alpha ordered by total degree, orthonormalized by
 * two Cholesky factorizations of their mass matrix on the region (the first
 * dim P^j functions span P^j). The Taylor operator is exact:
 * T_m psi(x) = sum_{|beta| <= m} (rho nu)^beta / beta! d^beta psi(x). */
template<typename T, size_t DIM>
class side_basis {

    using point_type = disk::point<T,DIM>;

    point_type                      xc;
    T                               hs = 1;
    std::vector<std::array<int,DIM>> ex;
    Mat<T>                          C;      // psi = mono * C (upper triangular)

    static T ipow(T x, int e) { T r = 1; for (int i = 0; i < e; i++) r *= x; return r; }
    static T falling(int a, int p) { T r = 1; for (int t = 0; t < p; t++) r *= (a - t); return r; }

    Mat<T> mono(const std::vector<point_type>& pts) const {
        Mat<T> V(pts.size(), ex.size());
        for (size_t q = 0; q < pts.size(); q++) {
            std::array<T,DIM> xi;
            for (size_t d = 0; d < DIM; d++)
                xi[d] = (pts[q][d] - xc[d]) / hs;
            for (size_t i = 0; i < ex.size(); i++) {
                T v = 1;
                for (size_t d = 0; d < DIM; d++)
                    v *= ipow(xi[d], ex[i][d]);
                V(q, i) = v;
            }
        }
        return V;
    }

public:
    side_basis() = default;

    /* the scaled monomials themselves (no orthonormalization) */
    side_basis(const point_type& center, T h, size_t degree)
        : xc(center), hs(h), ex(priv_exponents<DIM>(degree)), C(Mat<T>::Identity(ex.size(), ex.size()))
    {}

    side_basis(const std::vector<disk::quadrature_point<T,DIM>>& qps, const point_type& center, T h,
               size_t degree)
        : xc(center), hs(h), ex(priv_exponents<DIM>(degree))
    {
        size_t n = ex.size();
        std::vector<point_type> pts;
        Vec<T> w(qps.size());
        for (size_t q = 0; q < qps.size(); q++) {
            pts.push_back(qps[q].point());
            w(q) = qps[q].weight();
        }
        Mat<T> V = mono(pts);
        C = Mat<T>::Identity(n, n);
        for (size_t pass = 0; pass < 2; pass++) {
            Mat<T> P = V * C;
            Mat<T> M = P.transpose() * w.asDiagonal() * P;
            Eigen::LLT<Mat<T>> llt(M);
            if (llt.info() != Eigen::Success)
                throw std::runtime_error("bvc::side_basis: mass matrix not SPD (degenerate subcell)");
            Mat<T> Linv = llt.matrixL().solve(Mat<T>::Identity(n, n));
            C = C * Linv.transpose();
        }
    }

    size_t size() const { return ex.size(); }

    /* values at the points: one row per point */
    Mat<T> values(const std::vector<point_type>& pts) const { return mono(pts) * C; }

    /* derivatives: one matrix per direction */
    std::array<Mat<T>, DIM> gradients(const std::vector<point_type>& pts) const {
        std::array<Mat<T>, DIM> G;
        for (size_t d = 0; d < DIM; d++) {
            Mat<T> D(pts.size(), ex.size());
            for (size_t q = 0; q < pts.size(); q++) {
                std::array<T,DIM> xi;
                for (size_t e = 0; e < DIM; e++)
                    xi[e] = (pts[q][e] - xc[e]) / hs;
                for (size_t i = 0; i < ex.size(); i++) {
                    if (ex[i][d] == 0) {
                        D(q, i) = 0;
                        continue;
                    }
                    T v = ex[i][d] / hs;
                    for (size_t e = 0; e < DIM; e++)
                        v *= ipow(xi[e], ex[i][e] - (e == d ? 1 : 0));
                    D(q, i) = v;
                }
            }
            G[d] = D * C;
        }
        return G;
    }

    /* T_m psi at the points x_q, for the displacements rho_q nu_q */
    Mat<T> taylor(const std::vector<point_type>& pts, const std::vector<vec<T,DIM>>& nus,
                  const std::vector<T>& rhos, int m) const
    {
        if (m <= 0)
            return values(pts);
        auto betas = priv_exponents<DIM>(m);
        std::vector<T> inv_fact(betas.size());
        for (size_t b = 0; b < betas.size(); b++) {
            T f = 1;
            for (size_t d = 0; d < DIM; d++)
                for (int t = 2; t <= betas[b][d]; t++)
                    f *= t;
            inv_fact[b] = 1.0 / f;
        }
        Mat<T> D = Mat<T>::Zero(pts.size(), ex.size());
        for (size_t q = 0; q < pts.size(); q++) {
            std::array<T,DIM> xi, dl;
            for (size_t d = 0; d < DIM; d++) {
                xi[d] = (pts[q][d] - xc[d]) / hs;
                dl[d] = rhos[q] * nus[q](d) / hs;
            }
            for (size_t b = 0; b < betas.size(); b++) {
                T co = inv_fact[b];
                for (size_t d = 0; d < DIM; d++)
                    co *= ipow(dl[d], betas[b][d]);
                if (co == 0)
                    continue;
                for (size_t i = 0; i < ex.size(); i++) {
                    T v = co;
                    for (size_t d = 0; d < DIM && v != 0; d++) {
                        if (ex[i][d] < betas[b][d]) {
                            v = 0;
                            break;
                        }
                        v *= falling(ex[i][d], betas[b][d]) * ipow(xi[d], ex[i][d] - betas[b][d]);
                    }
                    D(q, i) += v;
                }
            }
        }
        return D * C;
    }
};

/***************************************************************************/
/*                          GEOMETRY OF A MACRO CELL                        */
/***************************************************************************/

/* A quadrature node of the polygonal interface */
template<typename T, size_t DIM>
struct interface_node {
    disk::point<T,DIM>  x;
    T                   w;
    vec<T,DIM>          nh;      // normal of the piece, from Omega_1 to Omega_2
    T                   rho;     // offset to Gamma along nh
    disk::point<T,DIM>  ph;      // x + rho nh, on Gamma
    vec<T,DIM>          nstar;   // normal of Gamma at ph
};

/* Root of smallest modulus of phi(x + r nu) = 0: Newton from r = 0 with
 * bounded steps, bisection as fallback (ProtoN's compute_rho) */
template<typename T, size_t DIM>
T
compute_rho(const level_set<T,DIM>& ls, const disk::point<T,DIM>& x, const vec<T,DIM>& nu, T hmax)
{
    auto at = [&](T r) { return x + disk::point<T,DIM>(vec<T,DIM>(nu * r)); };
    T r = 0;
    for (size_t it = 0; it < 60; it++) {
        auto y = at(r);
        T fv = ls(y);
        T g = ls.gradient(y).dot(nu);
        if (std::abs(g) < 1e-14)
            g = 1e-14;
        T step = std::max(-0.5*hmax, std::min(0.5*hmax, -fv/g));
        r = std::max(-2.0*hmax, std::min(2.0*hmax, r + step));
        if (std::abs(step) < 1e-13 * std::max(hmax, T(1)))
            break;
    }
    T scale = std::max(T(ls.gradient(x).norm()), T(1e-14));
    if (std::abs(ls(at(r))) / scale <= 1e-9 * std::max(hmax, T(1)))
        return r;
    T f0 = ls(x);
    for (T sgn : {1.0, -1.0}) {
        T hi = sgn * hmax;
        if ((ls(at(hi)) < 0) != (f0 < 0)) {
            T lo = 0;
            for (size_t it = 0; it < 80; it++) {
                T mid = 0.5 * (lo + hi);
                if ((ls(at(mid)) < 0) == (f0 < 0))
                    lo = mid;
                else
                    hi = mid;
            }
            return 0.5 * (lo + hi);
        }
    }
    return 0;
}

/* A piece of the boundary of a side of the macro cell on a face of the mesh */
template<typename Mesh>
struct face_piece {
    using T = typename Mesh::coordinate_type;
    typename Mesh::face_type    fc;
    size_t                      side;
    bool                        cut;    // the face is cut by Gamma_h
    vec<T, Mesh::dimension>     n;      // outward normal of the macro cell
    std::vector<disk::quadrature_point<T, Mesh::dimension>> qps;
};

template<typename Mesh>
struct macro_geometry {
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using qp_type = disk::quadrature_point<T, DIM>;

    bool                                cut = false;
    std::vector<size_t>                 sides;
    T                                   h = 0;
    std::array<std::vector<qp_type>, 2> vol;
    std::array<disk::point<T,DIM>, 2>   centroid;
    std::vector<face_piece<Mesh>>       pieces;
    std::vector<interface_node<T,DIM>>  inodes;
};

struct quadrature_degrees {
    size_t vol, face, itf;
};

template<typename Mesh>
macro_geometry<Mesh>
make_macro_geometry(const Mesh& msh, const macro_cell<Mesh>& Q, const level_set<typename Mesh::coordinate_type,
                    Mesh::dimension>& ls, const quadrature_degrees& qd)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    macro_geometry<Mesh> g;
    g.cut = Q.cut;
    g.h = Q.h;
    if (g.cut)
        g.sides = {0, 1};
    else
        g.sides = {Q.loc == disk::location::IN_NEGATIVE_SIDE ? size_t(0) : size_t(1)};

    for (auto s : g.sides) {
        for (auto m : Q.members) {
            auto mc = msh[m];
            if (disk::is_cut(msh, mc)) {
                auto q = integrate(msh, mc, qd.vol, side_location(s));
                g.vol[s].insert(g.vol[s].end(), q.begin(), q.end());
            }
            else if (disk::locate(msh, mc) == side_location(s)) {
                auto q = integrate(msh, mc, qd.vol);
                g.vol[s].insert(g.vol[s].end(), q.begin(), q.end());
            }
        }
        T a = 0;
        disk::point<T,DIM> c;
        for (auto& qp : g.vol[s]) {
            a += qp.weight();
            c = c + qp.point() * qp.weight();
        }
        g.centroid[s] = c / a;
    }

    for (auto& fm : Q.faces) {
        const auto& fc = fm.first;
        auto n = normal(msh, msh[fm.second], fc);
        auto floc = disk::locate(msh, fc);
        for (auto s : g.sides) {
            face_piece<Mesh> fp;
            fp.fc = fc;
            fp.side = s;
            fp.n = n;
            if (floc == disk::location::ON_INTERFACE) {
                fp.cut = true;
                fp.qps = integrate(msh, fc, qd.face, side_location(s));
            }
            else if (floc == side_location(s)) {
                fp.cut = false;
                fp.qps = integrate(msh, fc, qd.face);
            }
            else
                continue;
            if (!fp.qps.empty())
                g.pieces.push_back(fp);
        }
    }

    if (g.cut) {
        for (auto m : Q.members) {
            auto mc = msh[m];
            if (!disk::is_cut(msh, mc))
                continue;
            auto iq = disk::integrate_interface(msh, mc, qd.itf, disk::location::IN_NEGATIVE_SIDE);
            auto in = disk::interface_normals(msh, mc, qd.itf);
            for (size_t q = 0; q < iq.size(); q++) {
                interface_node<T,DIM> nd;
                nd.x = iq[q].point();
                nd.w = iq[q].weight();
                for (size_t d = 0; d < DIM; d++)
                    nd.nh(d) = in[q](d);
                nd.nh /= nd.nh.norm();
                nd.rho = compute_rho<T,DIM>(ls, nd.x, nd.nh, g.h);
                nd.ph = nd.x + disk::point<T,DIM>(vec<T,DIM>(nd.nh * nd.rho));
                vec<T,DIM> gp = ls.gradient(nd.ph);
                nd.nstar = gp / gp.norm();
                g.inodes.push_back(nd);
            }
        }
    }
    return g;
}

/***************************************************************************/
/*                               LOCAL SYSTEM                               */
/***************************************************************************/

enum class data_mode { exact, extended };

struct method_options {
    std::string name = "full";
    int         m = 0;          // order of the Taylor correction
    bool        gamma = true;   // correction gamma_T of the normal
    data_mode   data = data_mode::exact;
    double      eta = 1.0;      // interface penalty eta kappa_1 / h_T
    double      stab = 1.0;     // weight of s°: stab kappa_i / h_T
};

/* smallest m with h^{-1/2} delta_h^{m+1} <= h^{k+1} for delta_h = O(h^2) */
inline int taylor_order(size_t k) {
    return int(std::ceil(std::max(0.0, (double(k) - 0.5) / 2.0)));
}

inline method_options
make_variant(const std::string& name, size_t k)
{
    method_options o;
    o.name = name;
    if (name == "full")        { o.m = taylor_order(k); o.gamma = true; }
    else if (name == "taylor") { o.m = taylor_order(k); o.gamma = false; }
    else if (name == "gam")    { o.m = 0; o.gamma = true; }
    else if (name == "naive")  { o.m = 0; o.gamma = false; }
    else if (name == "ext")    { o.m = 0; o.gamma = false; o.data = data_mode::extended; }
    else
        throw std::invalid_argument("unknown variant " + name + " (full, taylor, gam, naive, ext)");
    return o;
}

/* Local unknowns: [ v_{T^i} for the sides of the cell | v_F for the face
 * pieces ], matrices (test rows, trial columns) */
template<typename Mesh>
struct local_system {
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

    Mat<T>                          A;
    Vec<T>                          b;
    std::array<side_basis<T,DIM>, 2> basis;
    std::array<size_t, 2>           coff = {0, 0};
    size_t                          cbs = 0, fbs = 0, nc = 0, ntot = 0;
    macro_geometry<Mesh>            geo;

    size_t poff(size_t piece) const { return nc + piece * fbs; }
};

template<typename Mesh>
std::vector<typename Mesh::point_type>
qp_points(const std::vector<disk::quadrature_point<typename Mesh::coordinate_type, Mesh::dimension>>& qps)
{
    std::vector<typename Mesh::point_type> ret;
    ret.reserve(qps.size());
    for (auto& qp : qps)
        ret.push_back(qp.point());
    return ret;
}

template<typename T, typename QPS>
Vec<T> qp_weights(const QPS& qps)
{
    Vec<T> w(qps.size());
    for (size_t q = 0; q < qps.size(); q++)
        w(q) = qps[q].weight();
    return w;
}

/* Face basis of a piece: L2-orthonormal on the face (uncut) or on its part */
template<typename Mesh>
disk::cut_face_basis<Mesh>
piece_basis(const Mesh& msh, const face_piece<Mesh>& fp, size_t k,
            const disk::cut_basis_options<typename Mesh::coordinate_type, Mesh::dimension>& bopt)
{
    auto where = fp.cut ? side_location(fp.side) : disk::locate(msh, fp.fc);
    return disk::cut_face_basis<Mesh>(msh, fp.fc, k, where, bopt);
}

template<typename Mesh>
local_system<Mesh>
make_local_system(const Mesh& msh, const macro_cell<Mesh>& Q, const problem<typename Mesh::coordinate_type,
                  Mesh::dimension>& prob, size_t k, const method_options& opt, const quadrature_degrees& qd,
                  const disk::cut_basis_options<typename Mesh::coordinate_type, Mesh::dimension>& bopt,
                  bool with_matrix = true)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

    local_system<Mesh> ls;
    ls.geo = make_macro_geometry(msh, Q, *prob.ls, qd);
    const auto& g = ls.geo;
    ls.cbs = dim_pk(k+1, DIM);
    ls.fbs = dim_pk(k, DIM-1);
    const size_t nsk = dim_pk(k, DIM);  // gradients in P^k(T^i)^DIM
    const size_t cbs = ls.cbs, fbs = ls.fbs;
    const T hT = g.h;
    const T* kap = prob.kappa;

    ls.nc = cbs * g.sides.size();
    ls.ntot = ls.nc + g.pieces.size() * fbs;
    for (size_t i = 0; i < g.sides.size(); i++)
        ls.coff[g.sides[i]] = i * cbs;
    for (auto s : g.sides)
        ls.basis[s] = side_basis<T,DIM>(g.vol[s], g.centroid[s], hT, k+1);

    const size_t ntot = ls.ntot;
    ls.A = Mat<T>::Zero(ntot, ntot);
    ls.b = Vec<T>::Zero(ntot);

    /* volume load */
    for (auto s : g.sides) {
        auto pts = qp_points<Mesh>(g.vol[s]);
        Vec<T> wf(pts.size());
        for (size_t q = 0; q < pts.size(); q++)
            wf(q) = g.vol[s][q].weight() * prob.f(pts[q], s);
        ls.b.segment(ls.coff[s], cbs) += ls.basis[s].values(pts).transpose() * wf;
    }
    if (!with_matrix)
        return ls;

    /* one-sided gradient reconstructions: Grec[s] (DIM*nsk x ntot), row
     * d*nsk + a: coefficient of psi_a e_d (orthonormal basis) */
    std::array<Mat<T>, 2> Grec;
    for (auto s : g.sides) {
        Grec[s] = Mat<T>::Zero(DIM*nsk, ntot);
        auto pts = qp_points<Mesh>(g.vol[s]);
        Vec<T> w = qp_weights<T>(g.vol[s]);
        Mat<T> Pk = ls.basis[s].values(pts).leftCols(nsk);
        auto dP = ls.basis[s].gradients(pts);
        for (size_t d = 0; d < DIM; d++)
            Grec[s].block(d*nsk, ls.coff[s], nsk, cbs) += Pk.transpose() * w.asDiagonal() * dP[d];
    }
    std::vector<Mat<T>> MF(g.pieces.size()), TF(g.pieces.size());
    for (size_t i = 0; i < g.pieces.size(); i++) {
        const auto& fp = g.pieces[i];
        size_t s = fp.side;
        auto pts = qp_points<Mesh>(fp.qps);
        Vec<T> w = qp_weights<T>(fp.qps);
        auto fb = piece_basis(msh, fp, k, bopt);
        Mat<T> PF = fb.eval_basis_batch(pts);
        Mat<T> PT = ls.basis[s].values(pts);
        Mat<T> PTk = PT.leftCols(nsk);
        for (size_t d = 0; d < DIM; d++) {
            Vec<T> wn = w * fp.n(d);
            Grec[s].block(d*nsk, ls.poff(i), nsk, fbs) += PTk.transpose() * wn.asDiagonal() * PF;
            Grec[s].block(d*nsk, ls.coff[s], nsk, cbs) -= PTk.transpose() * wn.asDiagonal() * PT;
        }
        MF[i] = PF.transpose() * w.asDiagonal() * PF;
        TF[i] = PF.transpose() * w.asDiagonal() * PT;
    }

    if (g.cut) {
        const size_t nq = g.inodes.size();
        std::vector<typename Mesh::point_type> xs;
        std::vector<vec<T,DIM>> nus;
        std::vector<T> rhos;
        Vec<T> W(nq);
        for (auto& in : g.inodes) {
            xs.push_back(in.x);
            nus.push_back(in.nh);
            rhos.push_back(in.rho);
        }
        for (size_t q = 0; q < nq; q++)
            W(q) = g.inodes[q].w;
        Mat<T> P0 = ls.basis[0].values(xs), P1 = ls.basis[1].values(xs);
        Mat<T> P0T = (opt.m > 0) ? ls.basis[0].taylor(xs, nus, rhos, opt.m) : P0;
        Mat<T> P1T = (opt.m > 0) ? ls.basis[1].taylor(xs, nus, rhos, opt.m) : P1;

        Mat<T> Jp = Mat<T>::Zero(nq, ntot), Jt = Mat<T>::Zero(nq, ntot);
        Jp.block(0, ls.coff[0], nq, cbs) = P0;
        Jp.block(0, ls.coff[1], nq, cbs) = -P1;
        Jt.block(0, ls.coff[0], nq, cbs) = P0T;
        Jt.block(0, ls.coff[1], nq, cbs) = -P1T;

        Mat<T> P0k = P0.leftCols(nsk), P1k = P1.leftCols(nsk);
        Mat<T> P0Tk = P0T.leftCols(nsk), P1Tk = P1T.leftCols(nsk);
        std::array<Vec<T>, DIM> wn, nhd, nsd;
        for (size_t d = 0; d < DIM; d++) {
            wn[d].resize(nq);
            nhd[d].resize(nq);
            nsd[d].resize(nq);
            for (size_t q = 0; q < nq; q++) {
                wn[d](q) = W(q) * g.inodes[q].nh(d);
                nhd[d](q) = g.inodes[q].nh(d);
                nsd[d](q) = g.inodes[q].nstar(d);
            }
        }
        /* lifting L_{T^1}(R) = sum_q w psi_a n_h,d R(q,:) */
        auto lift = [&](const Mat<T>& R) {
            Mat<T> ret(DIM*nsk, ntot);
            for (size_t d = 0; d < DIM; d++)
                ret.block(d*nsk, 0, nsk, ntot) = P0k.transpose() * wn[d].asDiagonal() * R;
            return ret;
        };
        Mat<T> Ghat0 = Grec[0] - lift(Jt);
        Mat<T> Gring0 = Grec[0] - lift(Jp);

        ls.A += kap[0] * Gring0.transpose() * Ghat0;
        ls.A += kap[1] * Grec[1].transpose() * Grec[1];
        ls.A += (opt.eta * kap[0] / hT) * Jp.transpose() * W.asDiagonal() * Jt;

        if (opt.gamma) {
            /* gamma_T(v,w) = (T_m(Q v).n* - (Q v).n_h, w_{T^2})_{Gamma_h} */
            Mat<T> Gam = Mat<T>::Zero(nq, ntot);
            for (size_t d = 0; d < DIM; d++) {
                Mat<T> G0d = Grec[0].block(d*nsk, 0, nsk, ntot);
                Mat<T> G1d = Grec[1].block(d*nsk, 0, nsk, ntot);
                Mat<T> Qd  = kap[0] * P0k * G0d - kap[1] * P1k * G1d;
                Mat<T> QTd = kap[0] * P0Tk * G0d - kap[1] * P1Tk * G1d;
                Gam += nsd[d].asDiagonal() * QTd - nhd[d].asDiagonal() * Qd;
            }
            ls.A.block(ls.coff[1], 0, cbs, ntot) += P1.transpose() * W.asDiagonal() * Gam;
        }

        /* loads */
        Vec<T> gD(nq), gN(nq);
        for (size_t q = 0; q < nq; q++) {
            const auto& in = g.inodes[q];
            if (opt.data == data_mode::exact) {
                gD(q) = prob.gD(in.ph);
                gN(q) = prob.gN(in.ph, in.nstar);
            }
            else {
                gD(q) = prob.gD(in.x);
                gN(q) = prob.gN(in.x, in.nh);
            }
        }
        Vec<T> WgD = W.cwiseProduct(gD);
        ls.b += (opt.eta * kap[0] / hT) * Jp.transpose() * WgD;
        Mat<T> Gn = Mat<T>::Zero(nq, ntot);
        for (size_t d = 0; d < DIM; d++)
            Gn += nhd[d].asDiagonal() * P0k * Gring0.block(d*nsk, 0, nsk, ntot);
        ls.b -= kap[0] * Gn.transpose() * WgD;
        ls.b.segment(ls.coff[1], cbs) += P1.transpose() * W.cwiseProduct(gN);
    }
    else {
        size_t s = g.sides[0];
        ls.A += kap[s] * Grec[s].transpose() * Grec[s];
    }

    /* stabilization: k_i/h_T (Pi_F v_{T^i} - v_F, Pi_F w_{T^i} - w_F)_F */
    for (size_t i = 0; i < g.pieces.size(); i++) {
        size_t s = g.pieces[i].side;
        Mat<T> R = Mat<T>::Zero(fbs, ntot);
        R.block(0, ls.coff[s], fbs, cbs) = MF[i].ldlt().solve(TF[i]);
        R.block(0, ls.poff(i), fbs, fbs) -= Mat<T>::Identity(fbs, fbs);
        ls.A += (opt.stab * kap[s] / hT) * R.transpose() * MF[i] * R;
    }
    return ls;
}

/* Interpolant of the exact solution: L2 projections on the cell unknowns and
 * on the face pieces */
template<typename Mesh>
Vec<typename Mesh::coordinate_type>
local_interpolant(const Mesh& msh, const problem<typename Mesh::coordinate_type, Mesh::dimension>& prob,
                  size_t k, const local_system<Mesh>& ls,
                  const disk::cut_basis_options<typename Mesh::coordinate_type, Mesh::dimension>& bopt)
{
    using T = typename Mesh::coordinate_type;
    Vec<T> I = Vec<T>::Zero(ls.ntot);
    for (auto s : ls.geo.sides) {
        auto pts = qp_points<Mesh>(ls.geo.vol[s]);
        Vec<T> wu(pts.size());
        for (size_t q = 0; q < pts.size(); q++)
            wu(q) = ls.geo.vol[s][q].weight() * prob.u(pts[q], s);
        I.segment(ls.coff[s], ls.cbs) = ls.basis[s].values(pts).transpose() * wu;
    }
    for (size_t i = 0; i < ls.geo.pieces.size(); i++) {
        const auto& fp = ls.geo.pieces[i];
        auto pts = qp_points<Mesh>(fp.qps);
        auto fb = piece_basis(msh, fp, k, bopt);
        Mat<T> PF = fb.eval_basis_batch(pts);
        Vec<T> w = qp_weights<T>(fp.qps), wu(pts.size());
        for (size_t q = 0; q < pts.size(); q++)
            wu(q) = w(q) * prob.u(pts[q], fp.side);
        Mat<T> MF = PF.transpose() * w.asDiagonal() * PF;
        I.segment(ls.poff(i), ls.fbs) = MF.ldlt().solve(PF.transpose() * wu);
    }
    return I;
}

/***************************************************************************/
/*                                ASSEMBLY                                  */
/***************************************************************************/

/* Global unknowns: P^k on each face piece, i.e. one set per uncut face and
 * one per side on a cut face, for the faces of the macro cells; the faces
 * inside a macro cell carry no unknown. Dirichlet data are imposed strongly
 * (L2 projection of the solution of the side of the piece), also on the cut
 * boundary faces. The cell unknowns are eliminated with an LU factorization:
 * the method is not symmetric. */
template<typename Mesh>
class assembler {
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

    size_t                              fbs, k;
    std::vector<std::array<long,2>>     fdof;   // global offset per (face, side), -1 if none
    std::vector<std::array<Vec<T>,2>>   fdir;   // Dirichlet values
    std::vector<bool>                   is_dir;
    std::vector<Eigen::Triplet<T>>      triplets;
    size_t                              ndof = 0;

public:
    Eigen::SparseMatrix<T>  LHS;
    Vec<T>                  RHS;

    assembler(const Mesh& msh, const std::vector<macro_cell<Mesh>>& macros,
              const problem<T,DIM>& prob, size_t degree,
              const disk::cut_basis_options<T,DIM>& bopt)
        : fbs(dim_pk(degree, DIM-1)), k(degree)
    {
        size_t nf = msh.faces_size();
        std::vector<std::array<bool,2>> used(nf, {false, false});
        for (auto& Q : macros) {
            std::vector<size_t> sides;
            if (Q.cut)
                sides = {0, 1};
            else
                sides = {Q.loc == disk::location::IN_NEGATIVE_SIDE ? size_t(0) : size_t(1)};
            for (auto& fm : Q.faces) {
                auto floc = disk::locate(msh, fm.first);
                auto f = offset(msh, fm.first);
                for (auto s : sides)
                    if (floc == disk::location::ON_INTERFACE || floc == side_location(s))
                        used[f][s] = true;
            }
        }
        fdof.assign(nf, {-1, -1});
        fdir.resize(nf);
        is_dir.assign(nf, false);
        size_t f = 0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, f++) {
            auto fc = *itor;
            bool cut = disk::locate(msh, fc) == disk::location::ON_INTERFACE;
            if (msh.is_boundary(fc)) {
                is_dir[f] = true;
                for (size_t s = 0; s < 2; s++) {
                    if (!used[f][s])
                        continue;
                    face_piece<Mesh> fp;
                    fp.fc = fc;
                    fp.side = s;
                    fp.cut = cut;
                    fp.qps = cut ? integrate(msh, fc, 2*k + 6, side_location(s)) : integrate(msh, fc, 2*k + 6);
                    auto fb = piece_basis(msh, fp, k, bopt);
                    auto pts = qp_points<Mesh>(fp.qps);
                    Mat<T> PF = fb.eval_basis_batch(pts);
                    Vec<T> w = qp_weights<T>(fp.qps), wu(pts.size());
                    for (size_t q = 0; q < pts.size(); q++)
                        wu(q) = w(q) * prob.u(pts[q], s);
                    Mat<T> MF = PF.transpose() * w.asDiagonal() * PF;
                    fdir[f][s] = MF.ldlt().solve(PF.transpose() * wu);
                }
                continue;
            }
            for (size_t s = 0; s < 2; s++)
                if (used[f][s]) {
                    fdof[f][s] = ndof;
                    ndof += fbs;
                }
        }
        RHS = Vec<T>::Zero(ndof);
        LHS = Eigen::SparseMatrix<T>(ndof, ndof);
    }

    size_t num_dofs() const { return ndof; }

    /* local face unknown -> global index (>= 0) or Dirichlet (-1) */
    void local_map(const Mesh& msh, const local_system<Mesh>& ls, std::vector<long>& map, Vec<T>& dval) const
    {
        size_t nf = ls.ntot - ls.nc;
        map.assign(nf, -1);
        dval = Vec<T>::Zero(nf);
        for (size_t i = 0; i < ls.geo.pieces.size(); i++) {
            const auto& fp = ls.geo.pieces[i];
            auto f = offset(msh, fp.fc);
            for (size_t j = 0; j < fbs; j++) {
                size_t lo = i*fbs + j;
                if (is_dir[f])
                    dval(lo) = fdir[f][fp.side](j);
                else
                    map[lo] = fdof[f][fp.side] + j;
            }
        }
    }

    /* static condensation of the local system and assembly */
    void assemble(const Mesh& msh, const local_system<Mesh>& ls, const Mat<T>& A, const Vec<T>& b,
                  const Eigen::PartialPivLU<Mat<T>>& lu)
    {
        size_t nc = ls.nc, nf = ls.ntot - ls.nc;
        Mat<T> AL = lu.solve(A.topRightCorner(nc, nf));
        Vec<T> bL = lu.solve(b.head(nc));
        Mat<T> S = A.bottomRightCorner(nf, nf) - A.bottomLeftCorner(nf, nc) * AL;
        Vec<T> r = b.tail(nf) - A.bottomLeftCorner(nf, nc) * bL;
        std::vector<long> map;
        Vec<T> dval;
        local_map(msh, ls, map, dval);
        for (size_t i = 0; i < nf; i++) {
            if (map[i] < 0)
                continue;
            for (size_t j = 0; j < nf; j++) {
                if (map[j] >= 0)
                    triplets.push_back(Eigen::Triplet<T>(map[i], map[j], S(i, j)));
                else
                    RHS(map[i]) -= S(i, j) * dval(j);
            }
            RHS(map[i]) += r(i);
        }
    }

    void finalize() {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        triplets.shrink_to_fit();
    }

    /* local unknowns (cells and faces) from the global solution */
    Vec<T> local_solution(const Mesh& msh, const local_system<Mesh>& ls, const Mat<T>& A, const Vec<T>& b,
                          const Eigen::PartialPivLU<Mat<T>>& lu, const Vec<T>& sol) const
    {
        size_t nc = ls.nc, nf = ls.ntot - ls.nc;
        std::vector<long> map;
        Vec<T> dval;
        local_map(msh, ls, map, dval);
        Vec<T> solF(nf);
        for (size_t i = 0; i < nf; i++)
            solF(i) = (map[i] >= 0) ? sol(map[i]) : dval(i);
        Vec<T> ret(ls.ntot);
        ret.tail(nf) = solF;
        ret.head(nc) = lu.solve(b.head(nc) - A.topRightCorner(nc, nf) * solF);
        return ret;
    }
};

/***************************************************************************/
/*                                  DRIVER                                  */
/***************************************************************************/

/* Local matrices of the uncut cells, computed once per shape (they do not
 * depend on the position of the cell; the load does) */
template<typename T>
struct uncut_entry {
    Mat<T>                      A;
    Eigen::PartialPivLU<Mat<T>> lu;
};

template<typename T>
struct run_result {
    size_t  cells = 0, cut = 0, macros = 0, agglomerated = 0, largest = 0, ndof = 0;
    T       h = 0, H1 = 0, L2 = 0, EN = 0, rho_max = 0, nonsym = 0;
    T       en2 = 0;    // a_h(e, e), negative if a_h is not coercive (coarse meshes)
    T       H1_ext = 0, L2_ext = 0;  // -p3: extension of the paired polynomial on the ill-cut subcells
    double  t_mesh = 0, t_asm = 0, t_solve = 0, t_err = 0;
};

/* One level with the stabilization by polynomial extension (-p3,
 * bvc_extension.hpp) */
template<typename Mesh, typename Timer>
run_result<typename Mesh::coordinate_type>
run_level_ext(const simulation_options& opts, std::shared_ptr<level_set<typename Mesh::coordinate_type,
              Mesh::dimension>> ls, size_t k, size_t l, const method_options& mopt);

#ifdef HAVE_INTEL_MKL
/* Preconditioner of the iterative solvers of Eigen: Cholesky factorization
 * (Pardiso) of the symmetric part of the matrix */
template<typename T>
class symmetric_part_preconditioner {
    std::function<Vec<T>(const Vec<T>&)> fact;
public:
    using StorageIndex = typename Eigen::SparseMatrix<T>::StorageIndex;
    enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic };
    symmetric_part_preconditioner() {}
    template<typename M> explicit symmetric_part_preconditioner(const M&) {}
    template<typename M> symmetric_part_preconditioner& analyzePattern(const M&) { return *this; }
    template<typename M> symmetric_part_preconditioner& factorize(const M&) { return *this; }
    template<typename M> symmetric_part_preconditioner& compute(const M&) { return *this; }
    void set(std::function<Vec<T>(const Vec<T>&)> f) { fact = f; }
    template<typename Rhs> Vec<T> solve(const Rhs& b) const { return fact(Vec<T>(b)); }
    Eigen::ComputationInfo info() const { return Eigen::Success; }
};
#endif

/* Solve of the unsymmetric condensed system (-S): 0 Pardiso LU (general
 * unsymmetric), 1 Pardiso LU for a structurally symmetric matrix, 2 GMRES
 * preconditioned by the Cholesky factorization of the symmetric part (the
 * unsymmetric part, due to the corrections, is small), LU if the symmetric
 * part is not positive definite or GMRES does not converge. SparseLU without
 * Intel MKL. */
template<typename T>
Vec<T>
solve_unsymmetric(const Eigen::SparseMatrix<T>& A, const Vec<T>& b, size_t kind = 0)
{
#ifdef HAVE_INTEL_MKL
    auto refine = [&](auto& solver, Vec<T> x) {
        /* Pardiso perturbs the small pivots: iterative refinement while the
         * relative residual exceeds 1e-12 */
        T nb = b.norm(), res0 = (b - A * x).norm() / nb, res = res0;
        size_t it = 0;
        for (; it < 10 && res > 1e-12; it++) {
            x += solver.solve(Vec<T>(b - A * x));
            res = (b - A * x).norm() / nb;
        }
        if (it > 0)
            std::cout << "         Pardiso: relative residual " << res0 << ", " << res << " after " << it
                      << " steps of iterative refinement" << std::endl;
        return x;
    };
    if (kind == 2) {
        Eigen::SparseMatrix<T> S = T(0.5) * (A + Eigen::SparseMatrix<T>(A.transpose()));
        /* Cholesky factorization of the symmetric part. If it is not
         * positive definite (k = 3 on the flowers, where the corrections make
         * the form non coercive on a few cells), the unsymmetric part is not
         * small: GMRES preconditioned by the LDL^T factorization of the
         * symmetric part needs 70 to 200 iterations, 5 to 7 times the time
         * of LU, which is used instead. */
        std::function<Vec<T>(const Vec<T>&)> fact;
        std::string name = "Cholesky";
        auto llt = std::make_shared<Eigen::PardisoLLT<Eigen::SparseMatrix<T>>>();
        llt->compute(S);
        if (llt->info() == Eigen::Success)
            fact = [llt](const Vec<T>& r) { return Vec<T>(llt->solve(r)); };
        if (fact) {
            /* tolerance on the true relative residual (Eigen's criterion is
             * on the preconditioned residual: 1e-10 there gave 10 times too
             * large errors for k = 3 on the finest 2D meshes): 1e-12 gives
             * the errors of LU to all digits in 2D and 3D, 1e-13 doubles or
             * triples the iterations in 3D for nothing */
            const T tol = std::getenv("BVC_GMRES_TOL") ? T(std::atof(std::getenv("BVC_GMRES_TOL"))) : T(1e-12);
            Eigen::GMRES<Eigen::SparseMatrix<T>, symmetric_part_preconditioner<T>> gmres;
            gmres.setTolerance(tol);
            gmres.setMaxIterations(100);
            gmres.set_restart(50);
            gmres.preconditioner().set(fact);
            gmres.compute(A);
            Vec<T> x = gmres.solve(b);
            size_t its = gmres.iterations();
            T nb = b.norm(), res = (b - A * x).norm() / nb;
            for (size_t pass = 0; pass < 3 && res > tol; pass++) {
                x = gmres.solveWithGuess(b, x);
                its += gmres.iterations();
                res = (b - A * x).norm() / nb;
            }
            if (res <= 10 * tol) {
                std::cout << "         GMRES (" << name << " of the symmetric part): " << its
                          << " iterations, relative residual " << res << std::endl;
                return x;
            }
            std::cout << "         GMRES (" << name << " of the symmetric part) did not converge (relative residual "
                      << res << "), LU used" << std::endl;
        }
        else
            std::cout << "         symmetric part not positive definite, LU used" << std::endl;
        kind = 1;
    }
    if (kind == 1) {
        pardiso_structsym_lu<Eigen::SparseMatrix<T>> solver;
        solver.compute(A);
        if (solver.info() == Eigen::Success) {
            Vec<T> x = solver.solve(b);
            if (solver.info() == Eigen::Success)
                return refine(solver, x);
        }
        std::cout << "         Pardiso (structurally symmetric) failed, unsymmetric LU used" << std::endl;
    }
    Eigen::PardisoLU<Eigen::SparseMatrix<T>> solver;
    solver.compute(A);
    if (solver.info() == Eigen::Success) {
        Vec<T> x = solver.solve(b);
        if (solver.info() == Eigen::Success)
            return refine(solver, x);
    }
    std::cout << "         PardisoLU failed, SparseLU used" << std::endl;
#endif
    Eigen::SparseLU<Eigen::SparseMatrix<T>> lu;
    lu.analyzePattern(A);
    lu.factorize(A);
    if (lu.info() != Eigen::Success)
        throw std::runtime_error("SparseLU: factorization failed");
    return lu.solve(b);
}

template<typename Mesh>
problem<typename Mesh::coordinate_type, Mesh::dimension>
make_problem(const simulation_options& opts, std::shared_ptr<level_set<typename Mesh::coordinate_type,
             Mesh::dimension>> ls, size_t k)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    T kap1 = 1.0, kap2 = std::pow(10.0, T(opts.contrast));
    switch (opts.test_case) {
        case 0: return make_contrast_problem<T,DIM>(ls, kap1, kap2);
        case 1: return make_discontinuous_problem<T,DIM>(ls, kap1, kap2);
        case 2: return make_patch_problem<T,DIM>(ls, k+1, kap1, opts.contrast == 0 ? T(7) : kap2);
        case 4: return make_jumps59_problem<T,DIM>(ls);
        default:
            throw std::invalid_argument("-p2: unknown test case " + std::to_string(opts.test_case) +
                                        " (0 contrast rho^6, 1 discontinuous, 2 patch test, 3 geometry, 4 jumps (59))");
    }
}

/* -f: geometry of the cut mesh of the level l for plotting: cells (location,
 * cut type, stabilizing cell, agglomeration group) and flat interface (2D:
 * chords in cells_l*.txt; 3D: triangles in flat_interface_l*.txt, one
 * triangle per line) */
template<typename Pre, typename Mesh>
void
write_geometry(const Mesh& msh, size_t l, const disk::cut_agglomeration& agglo = disk::cut_agglomeration())
{
    static const size_t DIM = Mesh::dimension;
    Pre::write_cells_info(msh, "cells_l" + std::to_string(l) + ".txt", agglo);
    if constexpr (DIM == 3) {
        std::ofstream ofs("flat_interface_l" + std::to_string(l) + ".txt");
        ofs << std::setprecision(10);
        for (auto& cl : msh) {
            if (!disk::is_cut(msh, cl))
                continue;
            for (auto& tri : disk::cut_info(msh, cl).poly.flat_interface) {
                for (auto& p : tri)
                    for (size_t d = 0; d < DIM; d++)
                        ofs << p[d] << " ";
                ofs << offset(msh, cl) << "\n";
            }
        }
    }
}

template<typename Mesh, typename Timer>
run_result<typename Mesh::coordinate_type>
run_level(const simulation_options& opts, std::shared_ptr<level_set<typename Mesh::coordinate_type,
          Mesh::dimension>> ls, size_t k, size_t l, const method_options& mopt)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using pre = preprocessor<Mesh>;
    using cut_mesh_type = typename pre::cut_mesh_type;

    run_result<T> res;
    Timer tc;

    /* mesh, flat interface, agglomeration */
    tc.tic();
    auto msh = pre::msh_generation(ls, l, opts, geometry_options(opts, k), false);
    auto agglo = disk::make_cell_agglomeration(msh, opts.debug != 0, opts.pairing != 0);
    auto macros = make_macro_cells(msh, agglo);
    res.t_mesh = tc.toc();
    if (opts.debug)
        write_geometry<pre>(msh, l, agglo);

    auto prob = make_problem<cut_mesh_type>(opts, ls, k);
    quadrature_degrees qd;
    if constexpr (DIM == 2)
        qd = {std::min<size_t>(2*(k+1) + 4, 16), 2*(k+1) + 2, 2*(k+1) + 6};
    else
        qd = {2*(k+1) + 2, 2*(k+1) + 2, 2*(k+1) + 2};
    qd.vol += opts.bvc_quad_extra;
    qd.face += opts.bvc_quad_extra;
    qd.itf += opts.bvc_quad_extra;

    disk::cut_basis_options<T,DIM> bopt;
    bopt.cache = std::make_shared<disk::cut_basis_cache<T,DIM>>();
    bopt.cache->degree = k+1;
    bopt.cache->face_degree = k;

    /* h: the largest diameter of the cells of the background mesh */
    res.cells = msh.cells_size();
    for (auto& cl : msh)
        res.h = std::max(res.h, diameter(msh, cl));
    res.macros = macros.size();
    for (auto& cl : msh)
        if (disk::is_cut(msh, cl))
            res.cut++;
    for (auto& Q : macros)
        if (Q.members.size() > 1) {
            res.agglomerated++;
            res.largest = std::max(res.largest, Q.members.size());
        }

    /* assembly */
    tc.tic();
    assembler<cut_mesh_type> as(msh, macros, prob, k, bopt);
    std::map<std::vector<long long>, uncut_entry<T>> cache;
    auto shape_key = [&](const macro_cell<cut_mesh_type>& Q, std::vector<long long>& key) {
        if (Q.cut || Q.members.size() != 1 || opts.cache_uncut == 0)
            return false;
        size_t s = (Q.loc == disk::location::IN_NEGATIVE_SIDE) ? 0 : 1;
        return uncut_cell_key(msh, msh[Q.members[0]], prob.kappa[s], key);
    };
    for (auto& Q : macros) {
        std::vector<long long> key;
        bool cacheable = shape_key(Q, key);
        if (cacheable) {
            auto it = cache.find(key);
            if (it != cache.end()) {
                auto ls_ = make_local_system(msh, Q, prob, k, mopt, qd, bopt, false);
                as.assemble(msh, ls_, it->second.A, ls_.b, it->second.lu);
                continue;
            }
        }
        auto ls_ = make_local_system(msh, Q, prob, k, mopt, qd, bopt);
        Eigen::PartialPivLU<Mat<T>> lu(ls_.A.topLeftCorner(ls_.nc, ls_.nc));
        as.assemble(msh, ls_, ls_.A, ls_.b, lu);
        if (Q.cut) {
            for (auto& in : ls_.geo.inodes)
                res.rho_max = std::max(res.rho_max, std::abs(in.rho));
            T nA = ls_.A.cwiseAbs().maxCoeff();
            if (nA > 0)
                res.nonsym = std::max(res.nonsym, (ls_.A - ls_.A.transpose()).cwiseAbs().maxCoeff() / nA);
        }
        if (cacheable)
            cache.emplace(key, uncut_entry<T>{ls_.A, lu});
    }
    as.finalize();
    res.t_asm = tc.toc();
    res.ndof = as.num_dofs();

    /* solve */
    tc.tic();
    Vec<T> sol = solve_unsymmetric(as.LHS, as.RHS, opts.bvc_solver);
    res.t_solve = tc.toc();

    /* errors on the polygonal subcells, energy error a_h(e, e) with
     * e = uhat_h - Ihat_h u */
    tc.tic();
    T e0 = 0, e1 = 0, en = 0;
    for (auto& Q : macros) {
        std::vector<long long> key;
        bool cacheable = shape_key(Q, key);
        local_system<cut_mesh_type> ls_;
        const Mat<T>* A = nullptr;
        const Eigen::PartialPivLU<Mat<T>>* lu = nullptr;
        Eigen::PartialPivLU<Mat<T>> lu_own;
        if (cacheable && cache.count(key)) {
            ls_ = make_local_system(msh, Q, prob, k, mopt, qd, bopt, false);
            A = &cache.at(key).A;
            lu = &cache.at(key).lu;
        }
        else {
            ls_ = make_local_system(msh, Q, prob, k, mopt, qd, bopt);
            lu_own.compute(ls_.A.topLeftCorner(ls_.nc, ls_.nc));
            A = &ls_.A;
            lu = &lu_own;
        }
        Vec<T> loc = as.local_solution(msh, ls_, *A, ls_.b, *lu, sol);
        for (auto s : ls_.geo.sides) {
            Vec<T> v = loc.segment(ls_.coff[s], ls_.cbs);
            auto pts = qp_points<cut_mesh_type>(ls_.geo.vol[s]);
            Vec<T> uh = ls_.basis[s].values(pts) * v;
            auto dP = ls_.basis[s].gradients(pts);
            for (size_t q = 0; q < pts.size(); q++) {
                T w = ls_.geo.vol[s][q].weight();
                T du = uh(q) - prob.u(pts[q], s);
                auto ge = prob.grad_u(pts[q], s);
                T dg = 0;
                for (size_t d = 0; d < DIM; d++) {
                    T gh = dP[d].row(q).dot(v);
                    dg += (gh - ge(d)) * (gh - ge(d));
                }
                e0 += w * du * du;
                e1 += w * dg;
            }
        }
        Vec<T> err = loc - local_interpolant(msh, prob, k, ls_, bopt);
        en += err.dot((*A) * err);
    }
    res.t_err = tc.toc();
    res.H1 = std::sqrt(e1);
    res.L2 = std::sqrt(e0);
    res.en2 = en;
    res.EN = std::sqrt(std::max(en, T(0)));
    return res;
}

/* Volume of {phi < 0} if it is star-shaped about the center of the domain:
 * the integral over the directions omega of r(omega)^DIM / DIM, r(omega)
 * the root of phi(c + r omega) (bisection); -1 if not available */
template<typename T, size_t DIM>
T
star_volume(const level_set<T,DIM>& ls)
{
    disk::point<T,DIM> c;
    for (size_t d = 0; d < DIM; d++)
        c[d] = 0.5;
    if (ls(c) >= 0)
        return -1;
    bool ok = true;
    auto radius = [&](const vec<T,DIM>& w) {
        auto at = [&](T r) { return c + disk::point<T,DIM>(vec<T,DIM>(w * r)); };
        T lo = 0, hi = 0.5 * std::sqrt(T(DIM));
        if (ls(at(hi)) < 0)
            ok = false;
        for (size_t it = 0; it < 100; it++) {
            T mid = 0.5 * (lo + hi);
            if (ls(at(mid)) < 0)
                lo = mid;
            else
                hi = mid;
        }
        return 0.5 * (lo + hi);
    };
    const size_t nt = 2048;
    T vol = 0;
    if constexpr (DIM == 2) {
        for (size_t j = 0; j < nt; j++) {
            T th = 2 * M_PI * j / nt;
            T r = radius(vec<T,2>(std::cos(th), std::sin(th)));
            vol += (2 * M_PI / nt) * r * r / 2;
        }
    }
    else {
        for (auto& tw : disk::cutgeo::gauss_legendre_01<T>(256)) {
            T ct = 2 * tw.first - 1, st = std::sqrt(std::max(T(0), 1 - ct * ct));
            for (size_t j = 0; j < nt / 4; j++) {
                T th = 2 * M_PI * j / (nt / 4);
                T r = radius(vec<T,3>(st * std::cos(th), st * std::sin(th), ct));
                vol += 2 * tw.second * (2 * M_PI / (nt / 4)) * r * r * r / 3;
            }
        }
    }
    return ok ? vol : T(-1);
}

/* Geometric verification (-t3, Table 1 and checks V5, V6 of the paper): for
 * each mesh, the cut cells, the small cuts (threshold -a) and the
 * agglomeration, the smallest subcell measure relative to the measure of a
 * cell, before and after agglomeration, delta_h = max |rho_h| and
 * theta_h = max |n_h - n*| at the quadrature nodes of Gamma_h and their root
 * mean squares on Gamma_h, the volume error |Omega_1h| - |Omega_1| and
 * max |T_m w - w(p_h)| on Gamma_h for a polynomial w of degree 4, which
 * converges as delta_h^{m+1} (Lemma 3.2) */
template<typename Mesh>
void
geometry_check(const simulation_options& opts)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using pre = preprocessor<Mesh>;
    using cut_mesh_type = typename pre::cut_mesh_type;

    auto ls = pre::make_level_set_function(opts);
    T vol_exact = star_volume<T,DIM>(*ls);
    if (vol_exact > 0)
        std::cout << "         |Omega_1| = " << std::setprecision(15) << vol_exact << std::setprecision(6) << std::endl;

    /* w = sum c_alpha x^alpha, degree 4 */
    side_basis<T,DIM> w4(disk::point<T,DIM>(), 1.0, 4);
    Vec<T> cw(w4.size());
    for (size_t i = 0; i < w4.size(); i++)
        cw(i) = std::sin(1.3 + 2.1 * i);

    std::ostringstream table;
    table << "#  l          h  cells    cut  small  aggl  max   min|T^i|/|T|   after      delta  delta/h^2"
          << "      theta  theta/h  rms(rho)/h^2  rms(theta)/h  |O1h|-|O1|     QG   max|T0w-w(ph)|   max|T1w-w(ph)|"
          << "   max|T2w-w(ph)|\n";
    size_t l_min = std::min(opts.l_min, opts.l_divs);
    for (size_t l = l_min; l <= opts.l_divs; l++) {
        auto msh = pre::msh_generation(ls, l, opts, geometry_options(opts, 1), false);
        T h = 0;
        for (auto& cl : msh)
            h = std::max(h, diameter(msh, cl));

        size_t ncut = 0, nsmall = 0;
        T min_before = 1, vol1 = 0;
        for (auto& cl : msh) {
            auto c = offset(msh, cl);
            if (!disk::is_cut(msh, cl)) {
                if (disk::locate(msh, cl) == disk::location::IN_NEGATIVE_SIDE)
                    vol1 += measure(msh, cl);
                continue;
            }
            ncut++;
            if (msh.cut_cell_data[c].cut != disk::cut_type::T_OK)
                nsmall++;
            T mn = measure(msh, cl, disk::location::IN_NEGATIVE_SIDE);
            T mp = measure(msh, cl, disk::location::IN_POSITIVE_SIDE);
            vol1 += mn;
            min_before = std::min(min_before, std::min(mn, mp) / measure(msh, cl));
        }

        auto agglo = disk::make_cell_agglomeration(msh, opts.debug != 0, opts.pairing != 0);
        auto macros = make_macro_cells(msh, agglo);
        size_t naggl = 0, maxm = 1;
        T min_after = 1;
        for (auto& Q : macros) {
            if (Q.members.size() > 1) {
                naggl++;
                maxm = std::max(maxm, Q.members.size());
            }
            if (!Q.cut)
                continue;
            T a[2] = {0, 0}, tot = 0;
            for (auto m : Q.members) {
                auto mc = msh[m];
                tot += measure(msh, mc);
                for (size_t s = 0; s < 2; s++)
                    a[s] += disk::is_cut(msh, mc) ? measure(msh, mc, side_location(s))
                          : (disk::locate(msh, mc) == side_location(s) ? measure(msh, mc) : T(0));
            }
            /* relative to the measure of a cell of the mesh (|T| = h^DIM) */
            min_after = std::min(min_after, std::min(a[0], a[1]) / (tot / Q.members.size()));
        }

        /* largest values and root mean squares on Gamma_h */
        T delta = 0, theta = 0, rho2 = 0, theta2 = 0, area = 0;
        std::array<T,3> terr = {0, 0, 0};
        size_t nq = 0;
        size_t qitf = (DIM == 2) ? 2*2 + 6 : 2*2 + 2;   // k = 1
        for (auto& cl : msh) {
            if (!disk::is_cut(msh, cl))
                continue;
            auto iq = disk::integrate_interface(msh, cl, qitf, disk::location::IN_NEGATIVE_SIDE);
            auto in = disk::interface_normals(msh, cl, qitf);
            T hT = diameter(msh, cl);
            std::vector<disk::point<T,DIM>> xs, ps;
            std::vector<vec<T,DIM>> nus;
            std::vector<T> rhos;
            for (size_t q = 0; q < iq.size(); q++) {
                vec<T,DIM> nh;
                for (size_t d = 0; d < DIM; d++)
                    nh(d) = in[q](d);
                nh /= nh.norm();
                T rho = compute_rho<T,DIM>(*ls, iq[q].point(), nh, hT);
                auto ph = iq[q].point() + disk::point<T,DIM>(vec<T,DIM>(nh * rho));
                vec<T,DIM> gp = ls->gradient(ph);
                T th = (nh - gp / gp.norm()).norm();
                delta = std::max(delta, std::abs(rho));
                theta = std::max(theta, th);
                T w = std::abs(iq[q].weight());
                rho2 += w * rho * rho;
                theta2 += w * th * th;
                area += w;
                xs.push_back(iq[q].point());
                ps.push_back(ph);
                nus.push_back(nh);
                rhos.push_back(rho);
            }
            nq += xs.size();
            Vec<T> wp = w4.values(ps) * cw;
            for (int m = 0; m <= 2; m++) {
                Vec<T> tw = w4.taylor(xs, nus, rhos, m) * cw;
                terr[m] = std::max(terr[m], (tw - wp).cwiseAbs().maxCoeff());
            }
        }

        table << std::setw(4) << l << std::scientific << std::setprecision(3) << std::setw(11) << h
              << std::setw(7) << msh.cells_size() << std::setw(7) << ncut << std::setw(7) << nsmall
              << std::setw(6) << naggl << std::setw(5) << maxm
              << std::setw(15) << min_before << std::setw(8) << std::fixed << std::setprecision(3) << min_after
              << std::scientific << std::setw(11) << delta << std::fixed << std::setw(11) << delta / (h*h)
              << std::scientific << std::setw(11) << theta << std::fixed << std::setw(9) << theta / h
              << std::setw(14) << std::sqrt(rho2 / area) / (h*h) << std::setw(14) << std::sqrt(theta2 / area) / h
              << std::scientific << std::setw(12) << (vol_exact > 0 ? vol1 - vol_exact : T(0))
              << std::setw(7) << nq << std::setw(17) << terr[0] << std::setw(17) << terr[1]
              << std::setw(17) << terr[2] << "\n";
        std::cout << "         Level " << l << ": " << ncut << " cut cells, " << nsmall << " small cuts, "
                  << naggl << " agglomerated cells, delta_h = " << delta << ", theta_h = " << theta << std::endl;
    }
    std::cout << std::endl << table.str();
    std::ofstream("bvc_geometry.txt") << table.str();
}

template<typename Mesh, typename Timer>
void
run(const simulation_options& opts)
{
    using T = typename Mesh::coordinate_type;
    using pre = preprocessor<Mesh>;

    if (opts.test_case == 3) {
        geometry_check<Mesh>(opts);
        return;
    }
    auto ls = pre::make_level_set_function(opts);

    std::string error_file_txt = "solution_error_file.txt";
    std::ofstream error_file(error_file_txt);
    error_file << std::setprecision(14);

    size_t k_min = std::min(opts.k_min, opts.degree);
    size_t l_min = std::min(opts.l_min, opts.l_divs);
    std::ostringstream table;

    for (size_t k = k_min; k <= opts.degree; k++) {
        auto mopt = make_variant(opts.bvc_variant, k);
        if (opts.bvc_taylor >= 0)
            mopt.m = opts.bvc_taylor;
        mopt.eta = opts.bvc_eta;
        mopt.stab = (opts.bvc_stab >= 0) ? opts.bvc_stab : std::pow(double(k+1), -opts.bvc_stab);
        auto prob0 = make_problem<typename pre::cut_mesh_type>(opts, ls, k);

        std::cout << std::endl << bold << red << "   Polynomial degree k : " << k << reset
                  << " (cells P^" << k+1 << ", faces P^" << k << "), variant " << mopt.name << " (m = " << mopt.m
                  << ", gamma_T " << (mopt.gamma ? "on" : "off") << ", data "
                  << (mopt.data == data_mode::exact ? "at p_h with n*" : "extended, at x with n_h") << ")"
                  << std::endl;
        std::cout << "         Test case: " << prob0.name << ", " << (size_t(1) << opts.bvc_chords)
                  << " chord(s) per cut cell (2D)" << std::endl;
        error_file << std::endl << "Polynomial degree k : " << k << std::endl;

        table << "# " << prob0.name << ", k = " << k << ", variant " << mopt.name << " (m = " << mopt.m
              << ", gamma " << (mopt.gamma ? 1 : 0) << "), "
              << (opts.prototype == 3 ? "polynomial extension" : "cell agglomeration") << "\n";
        table << "#    N  l          h     ndof          H1   rate          L2   rate          EN   rate\n";
        run_result<T> prev;
        for (size_t l = l_min; l <= opts.l_divs; l++) {
            std::cout << bold << cyan << "      Space refinment level -l : " << l << reset << std::endl;
            auto r = (opts.prototype == 3) ? run_level_ext<Mesh, Timer>(opts, ls, k, l, mopt)
                                           : run_level<Mesh, Timer>(opts, ls, k, l, mopt);
            std::cout << "         Mesh: " << r.cells << " cells, " << r.cut << " cut, " << r.agglomerated
                      << (opts.prototype == 3 ? " clusters of cells stabilized by extension" : " agglomerated cells")
                      << " (largest: " << r.largest << " cells), " << r.ndof << " face unknowns"
                      << ", max |rho_h|/h^2 = " << r.rho_max / (r.h*r.h) << ", non-symmetry " << r.nonsym << std::endl;
            std::cout << bold << yellow << "         H1-Error: " << r.H1 << "   L2-Error: " << r.L2
                      << "   Energy error: " << r.EN << reset << std::endl;
            if (opts.prototype == 3)
                std::cout << "         With the extension of the paired polynomial on the ill-cut subcells: H1-Error: "
                          << r.H1_ext << "   L2-Error: " << r.L2_ext << std::endl;
            if (r.en2 < 0)
                std::cout << "         Warning: a_h(e, e) = " << r.en2 << " < 0 (a_h is not coercive on this mesh):"
                          << " the energy error is reported as 0" << std::endl;
            /* order between the last two meshes, "-" if an error vanishes */
            auto rate = [&](T a, T b) -> std::string {
                if (!(prev.h > 0 && a > 0 && b > 0))
                    return "-";
                std::ostringstream os;
                os << std::fixed << std::setprecision(2) << std::log(a/b) / std::log(prev.h/r.h);
                return os.str();
            };
            if (prev.h > 0)
                std::cout << bold << yellow << "         order H1: " << rate(prev.H1, r.H1) << "   order L2: "
                          << rate(prev.L2, r.L2) << "   order energy: " << rate(prev.EN, r.EN) << reset << std::endl;
            std::cout << "         Times (s): mesh " << r.t_mesh << ", assembly " << r.t_asm << ", solve "
                      << r.t_solve << ", errors " << r.t_err << std::endl;

            error_file << "Space refinment level -l : " << l << std::endl;
            error_file << "Characteristic h size = " << r.h << std::endl;
            error_file << "L2-norm error = " << r.L2 << std::endl;
            error_file << "H1-norm error = " << r.H1 << std::endl;
            error_file << "Energy error = " << r.EN << std::endl;
            error_file << "a_h(e, e) = " << r.en2 << std::endl;
            if (opts.prototype == 3) {
                error_file << "H1-norm error with extension on ill-cut subcells = " << r.H1_ext << std::endl;
                error_file << "L2-norm error with extension on ill-cut subcells = " << r.L2_ext << std::endl;
            }
            error_file << "Face unknowns = " << r.ndof << std::endl;
            if (prev.h > 0) {
                error_file << "order L2 = " << rate(prev.L2, r.L2) << std::endl;
                error_file << "order H1 = " << rate(prev.H1, r.H1) << std::endl;
                error_file << "order energy = " << rate(prev.EN, r.EN) << std::endl;
            }
            error_file << std::endl;

            /* cells per direction of the structured meshes (0: hexagons, prisms) */
            size_t N = 0;
            if (opts.mesh != size_t(mesh_kind::polygons) && opts.mesh != size_t(mesh_kind::polyhedra))
                N = (opts.base_cells > 0 ? opts.base_cells : (Mesh::dimension == 2 ? 10 : 5)) * (size_t(1) << l);
            table << std::setw(6) << N << std::setw(3) << l << std::scientific << std::setprecision(3)
                  << std::setw(11) << r.h
                  << std::setw(9) << r.ndof;
            for (auto [a, b] : {std::make_pair(prev.H1, r.H1), std::make_pair(prev.L2, r.L2),
                                std::make_pair(prev.EN, r.EN)}) {
                table << std::scientific << std::setprecision(4) << std::setw(12) << b;
                table << std::setw(7) << rate(a, b);
            }
            table << "\n";
            prev = r;
        }
        table << "\n";
    }
    std::cout << std::endl << table.str();
    std::ofstream("bvc_convergence.txt") << table.str();
}

} // namespace bvc
