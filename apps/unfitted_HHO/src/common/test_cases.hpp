/*
 * Test cases of Section 5 of the paper (ported from ProtoN), in dimension 2
 * or 3 (rho is the distance to the center of the circle or of the sphere).
 * Omega_1 is the negative side of the level set, Omega_2 the positive one.
 */

#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include "diskpp/mesh/level_set.hpp"

template<typename T>
struct params {
    T kappa_1, kappa_2;
    params() : kappa_1(1.0), kappa_2(1.0) {}
    params(T kap1, T kap2) : kappa_1(kap1), kappa_2(kap2) {}
};

template<typename Mesh>
struct test_case_laplacian {

    using T          = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;
    using scalar_fun = std::function<T(const point_type&)>;
    using grad_fun   = std::function<Eigen::Matrix<T, 1, DIM>(const point_type&)>;

    std::string                         name;
    std::shared_ptr<level_set<T, DIM>>  level_set_;
    params<T>                           parms;
    scalar_fun                          sol_fun;
    scalar_fun                          rhs_fun;
    scalar_fun                          bcs_fun;
    grad_fun                            sol_grad;
    scalar_fun                          dirichlet_jump;   // g_D = u_1 - u_2
    scalar_fun                          neumann_jump;     // g_N = (kappa_1 grad u_1 - kappa_2 grad u_2).n_Gamma
};

namespace priv_tc {

/* prod_d sin(pi x_d) and its gradient */
template<typename T, size_t DIM>
T sin_product(const disk::point<T,DIM>& pt) {
    T v = 1;
    for (size_t d = 0; d < DIM; d++)
        v *= std::sin(M_PI*pt[d]);
    return v;
}

template<typename T, size_t DIM>
Eigen::Matrix<T, 1, DIM> sin_product_grad(const disk::point<T,DIM>& pt) {
    Eigen::Matrix<T, 1, DIM> ret;
    for (size_t d = 0; d < DIM; d++) {
        T v = M_PI * std::cos(M_PI*pt[d]);
        for (size_t e = 0; e < DIM; e++)
            if (e != d)
                v *= std::sin(M_PI*pt[e]);
        ret(d) = v;
    }
    return ret;
}

/* |x - c|^2 */
template<typename T, size_t DIM>
T rho2(const disk::point<T,DIM>& pt, const disk::point<T,DIM>& c) {
    T r2 = 0;
    for (size_t d = 0; d < DIM; d++)
        r2 += (pt[d] - c[d]) * (pt[d] - c[d]);
    return r2;
}

} // namespace priv_tc

/* Eq. (56): u = prod_d sin(pi x_d) in the whole domain, kappa_1 = kappa_2 = 1 */
template<typename Mesh>
test_case_laplacian<Mesh>
make_test_case_laplacian_sin_sin(std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    test_case_laplacian<Mesh> tc;
    tc.name = (DIM == 2) ? "sin(pi x) sin(pi y)" : "sin(pi x) sin(pi y) sin(pi z)";
    tc.level_set_ = ls;
    tc.parms = params<T>();
    if constexpr (DIM == 2) {
        tc.sol_fun = [](const point_type& pt) -> T {
            return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        };
        tc.rhs_fun = [](const point_type& pt) -> T {
            return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        };
        tc.sol_grad = [](const point_type& pt) -> Eigen::Matrix<T, 1, 2> {
            Eigen::Matrix<T, 1, 2> ret;
            ret(0) = M_PI * std::cos(M_PI*pt.x()) * std::sin(M_PI*pt.y());
            ret(1) = M_PI * std::sin(M_PI*pt.x()) * std::cos(M_PI*pt.y());
            return ret;
        };
    }
    else {
        tc.sol_fun  = [](const point_type& pt) -> T { return priv_tc::sin_product(pt); };
        tc.rhs_fun  = [](const point_type& pt) -> T { return DIM * M_PI * M_PI * priv_tc::sin_product(pt); };
        tc.sol_grad = [](const point_type& pt) { return priv_tc::sin_product_grad(pt); };
    }
    tc.bcs_fun = tc.sol_fun;
    tc.dirichlet_jump = [](const point_type&) -> T { return 0.0; };
    tc.neumann_jump   = [](const point_type&) -> T { return 0.0; };

    return tc;
}

template<typename T, size_t DIM>
const circle_level_set<T, DIM>&
as_circle(const std::shared_ptr<level_set<T, DIM>>& ls) {
    auto c = std::dynamic_pointer_cast<circle_level_set<T, DIM>>(ls);
    if (!c)
        throw std::invalid_argument("This test case requires a circular (spherical) interface (-i 0)");
    return *c;
}

/* Eq. (57): u_1 = rho^6/kappa_1, u_2 = rho^6/kappa_2 + R^6 (1/kappa_1 - 1/kappa_2), no jumps */
template<typename Mesh>
test_case_laplacian<Mesh>
make_test_case_laplacian_contrast_6(std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls,
                                    params<typename Mesh::coordinate_type> parms)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    const auto& circle = as_circle(ls);
    T R = circle.radius;
    point_type c = circle.center;
    T lap = 6.0 * (4.0 + DIM); /* Delta rho^6 = 6 (6 + DIM - 2) rho^4 */

    test_case_laplacian<Mesh> tc;
    tc.name = "contrast, rho^6 (57)";
    tc.level_set_ = ls;
    tc.parms = parms;
    tc.sol_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        if (r2 < R*R)
            return r2*r2*r2 / parms.kappa_1;
        return r2*r2*r2 / parms.kappa_2 + (R*R*R*R*R*R) * (1.0/parms.kappa_1 - 1.0/parms.kappa_2);
    };
    tc.rhs_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        return -lap * r2 * r2;
    };
    tc.bcs_fun = tc.sol_fun;
    tc.sol_grad = [=](const point_type& pt) -> Eigen::Matrix<T, 1, DIM> {
        Eigen::Matrix<T, 1, DIM> ret;
        T r2 = priv_tc::rho2(pt, c);
        T kappa = (r2 < R*R) ? parms.kappa_1 : parms.kappa_2;
        for (size_t d = 0; d < DIM; d++)
            ret(d) = 6 * r2 * r2 * (pt[d] - c[d]) / kappa;
        return ret;
    };
    tc.dirichlet_jump = [](const point_type&) -> T { return 0.0; };
    tc.neumann_jump   = [](const point_type&) -> T { return 0.0; };

    return tc;
}

/* Eq. (58a): u_1 = rho^6/kappa_1, u_2 = (rho^8 - R^8)/kappa_2 + R^6/kappa_1,
 * g_D = 0, g_N = 2R^5 (3 - 4R^2) */
template<typename Mesh>
test_case_laplacian<Mesh>
make_test_case_laplacian_contrast_jump_gN(std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls,
                                          params<typename Mesh::coordinate_type> parms)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    const auto& circle = as_circle(ls);
    T R = circle.radius;
    point_type c = circle.center;
    T R2 = R*R, R5 = R2*R2*R, R6 = R5*R, R7 = R6*R, R8 = R7*R;
    T lap6 = 6.0 * (4.0 + DIM), lap8 = 8.0 * (6.0 + DIM);

    test_case_laplacian<Mesh> tc;
    tc.name = "contrast, Neumann jump (58a)";
    tc.level_set_ = ls;
    tc.parms = parms;
    tc.sol_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        if (r2 < R2)
            return r2*r2*r2 / parms.kappa_1;
        return (r2*r2*r2*r2 - R8)/parms.kappa_2 + R6/parms.kappa_1;
    };
    tc.rhs_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        if (r2 < R2)
            return -lap6*r2*r2;
        return -lap8*r2*r2*r2;
    };
    tc.bcs_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        return (r2*r2*r2*r2 - R8)/parms.kappa_2 + R6/parms.kappa_1;
    };
    tc.sol_grad = [=](const point_type& pt) -> Eigen::Matrix<T, 1, DIM> {
        Eigen::Matrix<T, 1, DIM> ret;
        T r2 = priv_tc::rho2(pt, c);
        for (size_t d = 0; d < DIM; d++) {
            if (r2 < R2)
                ret(d) = 6 * r2 * r2 * (pt[d] - c[d]) / parms.kappa_1;
            else
                ret(d) = 8 * r2 * r2 * r2 * (pt[d] - c[d]) / parms.kappa_2;
        }
        return ret;
    };
    tc.dirichlet_jump = [](const point_type&) -> T { return 0.0; };
    tc.neumann_jump   = [=](const point_type&) -> T { return 6*R5 - 8*R7; };

    return tc;
}

/* Eq. (58b): u_1 = rho^6/kappa_1, u_2 = rho^6/kappa_2,
 * g_D = R^6 (1/kappa_1 - 1/kappa_2), g_N = 0 */
template<typename Mesh>
test_case_laplacian<Mesh>
make_test_case_laplacian_contrast_jump_gD(std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls,
                                          params<typename Mesh::coordinate_type> parms)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    const auto& circle = as_circle(ls);
    T R = circle.radius;
    point_type c = circle.center;
    T R6 = R*R*R*R*R*R;
    T lap = 6.0 * (4.0 + DIM);

    test_case_laplacian<Mesh> tc;
    tc.name = "contrast, Dirichlet jump (58b)";
    tc.level_set_ = ls;
    tc.parms = parms;
    tc.sol_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        T kappa = (r2 < R*R) ? parms.kappa_1 : parms.kappa_2;
        return r2*r2*r2/kappa;
    };
    tc.rhs_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        return -lap*r2*r2;
    };
    tc.bcs_fun = [=](const point_type& pt) -> T {
        T r2 = priv_tc::rho2(pt, c);
        return r2*r2*r2/parms.kappa_2;
    };
    tc.sol_grad = [=](const point_type& pt) -> Eigen::Matrix<T, 1, DIM> {
        Eigen::Matrix<T, 1, DIM> ret;
        T r2 = priv_tc::rho2(pt, c);
        T kappa = (r2 < R*R) ? parms.kappa_1 : parms.kappa_2;
        for (size_t d = 0; d < DIM; d++)
            ret(d) = 6 * r2 * r2 * (pt[d] - c[d]) / kappa;
        return ret;
    };
    tc.dirichlet_jump = [=](const point_type&) -> T { return R6*(1/parms.kappa_1 - 1/parms.kappa_2); };
    tc.neumann_jump   = [](const point_type&) -> T { return 0.0; };

    return tc;
}

/* Eq. (59): u_1 = cos(y) e^x, u_2 = prod_d sin(pi x_d), kappa_1 = kappa_2 = 1 */
template<typename Mesh>
test_case_laplacian<Mesh>
make_test_case_laplacian_jumps_2(std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    test_case_laplacian<Mesh> tc;
    tc.name = "non-polynomial jumps (59)";
    tc.level_set_ = ls;
    tc.parms = params<T>();
    if constexpr (DIM == 2) {
        tc.sol_fun = [ls](const point_type& pt) -> T {
            if ((*ls)(pt) > 0)
                return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
            return std::exp(pt.x()) * std::cos(pt.y());
        };
        tc.rhs_fun = [ls](const point_type& pt) -> T {
            if ((*ls)(pt) > 0)
                return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
            return 0.0;
        };
        tc.sol_grad = [ls](const point_type& pt) -> Eigen::Matrix<T, 1, 2> {
            Eigen::Matrix<T, 1, 2> ret;
            if ((*ls)(pt) > 0) {
                ret(0) = M_PI * std::cos(M_PI*pt.x()) * std::sin(M_PI*pt.y());
                ret(1) = M_PI * std::sin(M_PI*pt.x()) * std::cos(M_PI*pt.y());
            }
            else {
                ret(0) =   std::exp(pt.x()) * std::cos(pt.y());
                ret(1) = - std::exp(pt.x()) * std::sin(pt.y());
            }
            return ret;
        };
        tc.dirichlet_jump = [](const point_type& pt) -> T {
            return - std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y()) + std::exp(pt.x()) * std::cos(pt.y());
        };
        tc.neumann_jump = [ls](const point_type& pt) -> T {
            Eigen::Matrix<T, 2, 1> normal = ls->normal(pt);
            return -(M_PI * std::cos(M_PI*pt.x()) * std::sin(M_PI*pt.y()) - std::exp(pt.x()) * std::cos(pt.y())) * normal(0)
                   -(M_PI * std::sin(M_PI*pt.x()) * std::cos(M_PI*pt.y()) + std::exp(pt.x()) * std::sin(pt.y())) * normal(1);
        };
    }
    else {
        auto u1_grad = [](const point_type& pt) {
            Eigen::Matrix<T, 1, DIM> ret = Eigen::Matrix<T, 1, DIM>::Zero();
            ret(0) =   std::exp(pt.x()) * std::cos(pt.y());
            ret(1) = - std::exp(pt.x()) * std::sin(pt.y());
            return ret;
        };
        tc.sol_fun = [ls](const point_type& pt) -> T {
            if ((*ls)(pt) > 0)
                return priv_tc::sin_product(pt);
            return std::exp(pt.x()) * std::cos(pt.y());
        };
        tc.rhs_fun = [ls](const point_type& pt) -> T {
            if ((*ls)(pt) > 0)
                return DIM * M_PI * M_PI * priv_tc::sin_product(pt);
            return 0.0;
        };
        tc.sol_grad = [ls, u1_grad](const point_type& pt) -> Eigen::Matrix<T, 1, DIM> {
            if ((*ls)(pt) > 0)
                return priv_tc::sin_product_grad(pt);
            return u1_grad(pt);
        };
        tc.dirichlet_jump = [](const point_type& pt) -> T {
            return std::exp(pt.x()) * std::cos(pt.y()) - priv_tc::sin_product(pt);
        };
        tc.neumann_jump = [ls, u1_grad](const point_type& pt) -> T {
            Eigen::Matrix<T, DIM, 1> normal = ls->normal(pt);
            return (u1_grad(pt) - priv_tc::sin_product_grad(pt)).dot(normal.transpose());
        };
    }
    tc.bcs_fun = tc.sol_fun;

    return tc;
}

template<typename Mesh>
test_case_laplacian<Mesh>
select_test_case(size_t id, std::shared_ptr<level_set<typename Mesh::coordinate_type, Mesh::dimension>> ls,
                 params<typename Mesh::coordinate_type> parms)
{
    switch (id) {
        case 0: return make_test_case_laplacian_sin_sin<Mesh>(ls);
        case 1: return make_test_case_laplacian_contrast_6<Mesh>(ls, parms);
        case 2: return make_test_case_laplacian_contrast_jump_gN<Mesh>(ls, parms);
        case 3: return make_test_case_laplacian_contrast_jump_gD<Mesh>(ls, parms);
        case 4: return make_test_case_laplacian_jumps_2<Mesh>(ls);
        default:
            throw std::invalid_argument("Unknown test case id: " + std::to_string(id));
    }
}
