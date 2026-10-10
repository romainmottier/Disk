/*
 *       /\        Matteo Cicuttin (C) 2017,2018; Guillaume Delay 2018,2019
 *      /__\       matteo.cicuttin@enpc.fr        guillaume.delay@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Bases for unfitted HHO methods, ported from ProtoN, in dimension 2 or 3.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "diskpp/bases/bases_all.hpp"
#include "diskpp/mesh/cut_mesh.hpp"

namespace disk {

/* Cell bases of the unfitted HHO methods. Each basis is attached to a region:
 * the cell if it is not cut, else its subcell, an ill-cut subcell being merged
 * with the same side of its stabilizing cell (Section 5.1 of the paper).
 *  - cell_monomial: scaled monomials centered at the barycenter of the whole
 *    cell (dashed lines of Fig. 7);
 *  - subcell_monomial: scaled monomials centered at the barycenter of the
 *    region and scaled by its diameter (ProtoN, figures of the paper);
 *  - orthonormal: tensor Legendre polynomials in the bounding box of the
 *    region in its principal axes of inertia, orthonormalized for the mean L2
 *    product of the region (1/|R|)(u,v)_R. The basis is well conditioned
 *    whatever the shape of the region (thin subcells) and the degree. */
enum class cut_basis_kind { cell_monomial, subcell_monomial, orthonormal };

/* Local frame of an orthonormal basis: x -> xi = axes (x - center) maps the
 * region into [-1,1]^DIM, and the basis is transform * (Legendre products). */
template<typename T, size_t DIM = 2>
struct cut_basis_frame {
    point<T,DIM>                                     center;
    Eigen::Matrix<T, DIM, DIM>                       axes;
    size_t                                           degree = 0;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> transform; // lower triangular, empty: identity
};

/* Frame of an orthonormal basis of a face of a 3D mesh: x -> u = U^T (x -
 * origin) in the coordinates of its plane, then the frame of these
 * coordinates */
template<typename T>
struct cut_face_frame {
    point<T,3>                  origin;
    Eigen::Matrix<T, 3, 2>      U;
    cut_basis_frame<T, 2>       frame;
};

/* Bounded map keeping the most recently used values */
template<typename Key, typename Value>
class lru_map {
    using list_type = std::list<Key>;
    std::map<Key, std::pair<Value, typename list_type::iterator>> m_map;
    list_type   m_order;
    size_t      m_capacity;

public:
    explicit lru_map(size_t capacity = 4096) : m_capacity(capacity) {}

    const Value*
    find(const Key& k) {
        auto it = m_map.find(k);
        if (it == m_map.end())
            return nullptr;
        m_order.splice(m_order.begin(), m_order, it->second.second);
        return &it->second.first;
    }

    void
    insert(const Key& k, const Value& v) {
        if (m_map.find(k) != m_map.end())
            return;
        while (m_map.size() >= m_capacity && !m_order.empty()) {
            m_map.erase(m_order.back());
            m_order.pop_back();
        }
        m_order.push_front(k);
        m_map[k] = std::make_pair(v, m_order.begin());
    }
};

/* Frames of the regions of a mesh, computed once (keys: 2 cell + side for the
 * cut cells, other keys for the agglomerated cells; the frames of the uncut
 * cells and of the faces of 3D meshes are kept in bounded maps). Must be
 * reset when the mesh changes. */
template<typename T, size_t DIM = 2>
struct cut_basis_cache {
    size_t                                                              degree = 0;      // degree of the cell frames
    size_t                                                              face_degree = 0; // degree of the face frames
    std::map<size_t, std::shared_ptr<const cut_basis_frame<T,DIM>>>     frames;
    lru_map<size_t, std::shared_ptr<const cut_basis_frame<T,DIM>>>      uncut_frames;
    std::map<size_t, std::shared_ptr<const cut_face_frame<T>>>          face_frames;
    lru_map<size_t, std::shared_ptr<const cut_face_frame<T>>>           uncut_face_frames;
};

template<typename T, size_t DIM = 2>
struct cut_basis_options {
    cut_basis_kind                          kind = cut_basis_kind::orthonormal;
    std::shared_ptr<cut_basis_cache<T,DIM>> cache;
};

/* dimension of P^degree in DIM variables */
inline size_t
cut_basis_frame_size(size_t degree, size_t dim = 2) {
    size_t n = 1;
    for (size_t i = 1; i <= dim; i++)
        n = n * (degree + i) / i;
    return n;
}

namespace priv {

/* Exponents of the monomials of total degree <= degree, ordered by total
 * degree as in ProtoN: x^(k-i) y^i in 2D, x^(k-i) y^(i-j) z^j in 3D */
template<size_t DIM>
const std::vector<std::array<size_t, DIM>>&
cut_monomial_exponents(size_t degree)
{
    static std::map<size_t, std::vector<std::array<size_t, DIM>>> cache;
    auto it = cache.find(degree);
    if (it != cache.end())
        return it->second;
    std::vector<std::array<size_t, DIM>> ret;
    for (size_t k = 0; k <= degree; k++) {
        if constexpr (DIM == 2) {
            for (size_t i = 0; i <= k; i++)
                ret.push_back({k-i, i});
        }
        else {
            for (size_t i = 0; i <= k; i++)
                for (size_t j = 0; j <= i; j++)
                    ret.push_back({k-i, i-j, j});
        }
    }
    cache[degree] = ret;
    return cache[degree];
}

/* sqrt(2n+1) P_n(x) (mean square 1 on [-1,1]) and derivatives, n = 0..degree */
template<typename T>
void
normalized_legendre(T x, size_t degree, std::vector<T>& val, std::vector<T>& der)
{
    val.resize(degree+1);
    der.resize(degree+1);
    val[0] = 1.0;
    der[0] = 0.0;
    if (degree >= 1) {
        val[1] = x;
        der[1] = 1.0;
    }
    for (size_t n = 1; n < degree; n++) {
        val[n+1] = ((2.0*n + 1.0) * x * val[n] - n * val[n-1]) / (n + 1.0);
        der[n+1] = der[n-1] + (2.0*n + 1.0) * val[n];
    }
    for (size_t n = 0; n <= degree; n++) {
        auto s = std::sqrt(2.0*n + 1.0);
        val[n] *= s;
        der[n] *= s;
    }
}

/* Legendre products psi_a = prod_d L_(a_d)(xi_d), |a| <= degree, in the order
 * of cut_monomial_exponents */
template<typename T, size_t DIM, typename Grad>
void
legendre_products(const cut_basis_frame<T,DIM>& fr, const point<T,DIM>& pt, size_t degree,
                  Eigen::Matrix<T, Eigen::Dynamic, 1>* val, Grad* grad)
{
    Eigen::Matrix<T, DIM, 1> xi = fr.axes * (pt - fr.center).to_vector();

    thread_local std::vector<T> lv[DIM], dlv[DIM];
    for (size_t d = 0; d < DIM; d++)
        normalized_legendre(xi(d), degree, lv[d], dlv[d]);

    const auto& expo = cut_monomial_exponents<DIM>(degree);
    for (size_t pos = 0; pos < expo.size(); pos++) {
        const auto& a = expo[pos];
        if (val) {
            T v = 1;
            for (size_t d = 0; d < DIM; d++)
                v *= lv[d][a[d]];
            (*val)(pos) = v;
        }
        if (grad) {
            /* d psi / d xi_d, then chain rule */
            Eigen::Matrix<T, DIM, 1> dxi;
            for (size_t d = 0; d < DIM; d++) {
                T v = dlv[d][a[d]];
                for (size_t e = 0; e < DIM; e++)
                    if (e != d)
                        v *= lv[e][a[e]];
                dxi(d) = v;
            }
            if constexpr (DIM == 2) {
                (*grad)(pos,0) = dxi(0)*fr.axes(0,0) + dxi(1)*fr.axes(1,0);
                (*grad)(pos,1) = dxi(0)*fr.axes(0,1) + dxi(1)*fr.axes(1,1);
            }
            else
                grad->row(pos) = (fr.axes.transpose() * dxi).transpose();
        }
    }
}

/* Legendre products at the points pts (rows), and their derivatives with
 * respect to x_d (if grads is not null) */
template<typename T, size_t DIM, typename Points>
void
legendre_products_batch(const cut_basis_frame<T,DIM>& fr, const Points& pts, size_t degree,
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* vals,
                        std::array<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, DIM>* grads)
{
    const auto& expo = cut_monomial_exponents<DIM>(degree);
    size_t n = expo.size(), np = pts.size();
    if (vals)
        vals->resize(np, n);
    if (grads)
        for (size_t d = 0; d < DIM; d++)
            (*grads)[d].resize(np, n);

    std::vector<T> lv[DIM], dlv[DIM];
    for (size_t q = 0; q < np; q++) {
        Eigen::Matrix<T, DIM, 1> xi = fr.axes * (pts[q] - fr.center).to_vector();
        for (size_t d = 0; d < DIM; d++)
            normalized_legendre(xi(d), degree, lv[d], dlv[d]);
        for (size_t pos = 0; pos < n; pos++) {
            const auto& a = expo[pos];
            if (vals) {
                T v = 1;
                for (size_t d = 0; d < DIM; d++)
                    v *= lv[d][a[d]];
                (*vals)(q, pos) = v;
            }
            if (grads) {
                Eigen::Matrix<T, DIM, 1> dxi;
                for (size_t d = 0; d < DIM; d++) {
                    T v = dlv[d][a[d]];
                    for (size_t e = 0; e < DIM; e++)
                        if (e != d)
                            v *= lv[e][a[e]];
                    dxi(d) = v;
                }
                Eigen::Matrix<T, DIM, 1> g = fr.axes.transpose() * dxi;
                for (size_t d = 0; d < DIM; d++)
                    (*grads)[d](q, pos) = g(d);
            }
        }
    }
}

/* Inverse of the Cholesky factor of a SPD matrix (empty matrix on failure) */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
inverse_cholesky_factor(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& G)
{
    Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> llt(G);
    if (llt.info() != Eigen::Success)
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Id
        = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(G.rows(), G.cols());
    return llt.matrixL().solve(Id);
}

/* Frame of the orthonormal basis of degree "degree" of a region of R^DIM given
 * by a quadrature and points of its boundary (bounding box) */
template<typename T, size_t DIM>
std::shared_ptr<cut_basis_frame<T,DIM>>
make_frame(const std::vector<quadrature_point<T, DIM>>& qps, const std::vector<point<T,DIM>>& vertices, size_t degree)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto fr = std::make_shared<cut_basis_frame<T,DIM>>();
    fr->degree = degree;

    /* Principal axes of inertia */
    T area = 0.0;
    point<T,DIM> c;
    for (auto& qp : qps) {
        area += qp.weight();
        c = c + qp.point() * qp.weight();
    }
    c = c / area;
    Eigen::Matrix<T, DIM, DIM> C = Eigen::Matrix<T, DIM, DIM>::Zero();
    for (auto& qp : qps) {
        Eigen::Matrix<T, DIM, 1> r = (qp.point() - c).to_vector();
        C += qp.weight() * r * r.transpose();
    }
    Eigen::Matrix<T, DIM, DIM> U = Eigen::Matrix<T, DIM, DIM>::Identity(); // columns: axes
    T offdiag = 0;
    for (size_t i = 0; i < DIM; i++)
        for (size_t j = i+1; j < DIM; j++)
            offdiag = std::max(offdiag, std::abs(C(i,j)));
    if (offdiag > 1e-12 * C.trace()) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, DIM, DIM>> es(C);
        U = es.eigenvectors();
        /* canonical sign of the axes (first significant component positive):
         * the translates of a region get the same basis (cache of the local
         * operators of the uncut cells) */
        for (size_t j = 0; j < DIM; j++) {
            size_t i = 0;
            while (i + 1 < DIM && std::abs(U(i,j)) < 1e-6)
                i++;
            if (U(i,j) < 0)
                U.col(j) = -U.col(j);
        }
    }

    /* Bounding box of the region in these axes */
    Eigen::Matrix<T, DIM, 1> lo, hi;
    lo.setConstant(std::numeric_limits<T>::max());
    hi.setConstant(-std::numeric_limits<T>::max());
    for (auto& v : vertices) {
        Eigen::Matrix<T, DIM, 1> p = U.transpose() * (v - c).to_vector();
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    Eigen::Matrix<T, DIM, 1> mid  = 0.5 * (lo + hi);
    Eigen::Matrix<T, DIM, 1> half = 0.5 * (hi - lo);
    Eigen::Matrix<T, DIM, 1> cm = U * mid;
    fr->center = c + point<T,DIM>(cm);
    for (size_t i = 0; i < DIM; i++)
        fr->axes.row(i) = U.col(i).transpose() / half(i);

    /* Orthonormalization for (1/|R|)(.,.)_R: Cholesky QR, twice */
    auto n = cut_basis_frame_size(degree, DIM);
    matrix_type G = matrix_type::Zero(n, n);
    vector_type psi(n);
    for (auto& qp : qps) {
        legendre_products(*fr, qp.point(), degree, &psi, (Eigen::Matrix<T, Eigen::Dynamic, DIM>*)nullptr);
        G += qp.weight() * psi * psi.transpose();
    }
    G /= area;

    matrix_type A = inverse_cholesky_factor(G);
    if (A.size() == 0) {
        std::cout << "Warning: orthonormalization of a cut basis failed, diagonal scaling used" << std::endl;
        fr->transform = G.diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
        return fr;
    }
    matrix_type G2 = A * G * A.transpose();
    matrix_type A2 = inverse_cholesky_factor(G2);
    if (A2.size() != 0)
        A = A2 * A;
    fr->transform = A.template triangularView<Eigen::Lower>();

    return fr;
}

/* Frame of an axis-aligned box given by its vertices (nullptr if they are not
 * the vertices of such a box): the Legendre products are orthonormal. */
template<typename T, size_t DIM>
std::shared_ptr<cut_basis_frame<T,DIM>>
box_frame(const std::vector<point<T,DIM>>& pts, size_t degree)
{
    if (pts.size() != (size_t(1) << DIM))
        return nullptr;
    Eigen::Matrix<T, DIM, 1> lo = pts[0].to_vector(), hi = pts[0].to_vector();
    for (auto& p : pts) {
        lo = lo.cwiseMin(p.to_vector());
        hi = hi.cwiseMax(p.to_vector());
    }
    T tol = 1e-12 * (hi - lo).maxCoeff();
    for (auto& p : pts)
        for (size_t d = 0; d < DIM; d++)
            if (std::abs(p[d] - lo(d)) > tol && std::abs(p[d] - hi(d)) > tol)
                return nullptr;

    auto fr = std::make_shared<cut_basis_frame<T,DIM>>();
    fr->degree = degree;
    fr->center = point<T,DIM>(Eigen::Matrix<T, DIM, 1>(0.5*(lo + hi)));
    fr->axes = Eigen::Matrix<T, DIM, DIM>::Zero();
    for (size_t d = 0; d < DIM; d++)
        fr->axes(d,d) = 2.0/(hi(d) - lo(d));
    return fr;
}

} // namespace priv

/* Frame of the orthonormal basis of degree "degree" on the union of the
 * pieces (cell, side) (the whole cell if it is not cut). */
template<CutMesh Mesh>
std::shared_ptr<const cut_basis_frame<typename Mesh::coordinate_type, Mesh::dimension>>
make_cut_basis_frame(const Mesh& msh,
                     const std::vector<std::pair<typename Mesh::cell_type, location>>& pieces,
                     size_t degree)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

    std::vector<quadrature_point<T, DIM>> qps;
    std::vector<point<T,DIM>> vertices;
    for (auto& pc : pieces) {
        const auto& cl = pc.first;
        auto where = is_cut(msh, cl) ? pc.second : locate(msh, cl);
        auto q = integrate(msh, cl, std::max<size_t>(2, 2*degree), where);
        qps.insert(qps.end(), q.begin(), q.end());
        auto rp = region_points(msh, cl, where);
        vertices.insert(vertices.end(), rp.begin(), rp.end());
    }
    return priv::make_frame(qps, vertices, degree);
}

/* Region of the basis of the side "where" of the cell cl: the cell if it is
 * not cut, else its subcell, merged with its stabilizing cell if it is ill-cut */
template<CutMesh Mesh>
std::vector<std::pair<typename Mesh::cell_type, location>>
cut_basis_region(const Mesh& msh, const typename Mesh::cell_type& cl, location where)
{
    std::vector<std::pair<typename Mesh::cell_type, location>> ret;
    if (!is_cut(msh, cl)) {
        ret.push_back( std::make_pair(cl, locate(msh, cl)) );
        return ret;
    }

    ret.push_back( std::make_pair(cl, where) );
    const auto& ccd = cut_info(msh, cl);
    bool ill_cut_side = (ccd.cut == cut_type::T_KO_NEG and where == location::IN_NEGATIVE_SIDE)
                     or (ccd.cut == cut_type::T_KO_POS and where == location::IN_POSITIVE_SIDE);
    if (ill_cut_side and ccd.stabilizing_cell.has_value())
        ret.push_back( std::make_pair(msh[ccd.stabilizing_cell.value()], where) );

    return ret;
}

/* Frame of an uncut cell which is an axis-aligned box (nullptr if it is
 * not): the Legendre products of its bounding box are already orthonormal. */
template<CutMesh Mesh>
std::shared_ptr<const cut_basis_frame<typename Mesh::coordinate_type, Mesh::dimension>>
rectangle_basis_frame(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    auto pts = points(msh, cl);
    return priv::box_frame(std::vector<typename Mesh::point_type>(pts.begin(), pts.end()), degree);
}

/* Frame of the orthonormal basis of the side "where" of cl. The frames of the
 * cut cells are kept in the cache, if any; those of the uncut cells are cheap
 * (axis-aligned boxes) or kept in a bounded map. */
template<CutMesh Mesh>
std::shared_ptr<const cut_basis_frame<typename Mesh::coordinate_type, Mesh::dimension>>
cut_cell_basis_frame(const Mesh& msh, const typename Mesh::cell_type& cl, location where,
                     size_t degree, cut_basis_cache<typename Mesh::coordinate_type, Mesh::dimension>* cache)
{
    using frame_ptr = std::shared_ptr<const cut_basis_frame<typename Mesh::coordinate_type, Mesh::dimension>>;
    bool cut = is_cut(msh, cl);
    if (!cut) {
        auto fr = rectangle_basis_frame(msh, cl, degree);
        if (fr)
            return fr;
    }
    if (!cache)
        return make_cut_basis_frame(msh, cut_basis_region(msh, cl, where), degree);

    size_t cl_id = offset(msh, cl);
    if (!cut) {
        if (auto f = cache->uncut_frames.find(cl_id))
            if ((*f)->degree >= degree)
                return *f;
        frame_ptr fr = make_cut_basis_frame(msh, cut_basis_region(msh, cl, where), std::max(degree, cache->degree));
        cache->uncut_frames.insert(cl_id, fr);
        return fr;
    }

    size_t key = 2*cl_id + (where == location::IN_POSITIVE_SIDE ? 1 : 0);
    auto it = cache->frames.find(key);
    if (it != cache->frames.end() && it->second->degree >= degree)
        return it->second;

    frame_ptr fr = make_cut_basis_frame(msh, cut_basis_region(msh, cl, where), std::max(degree, cache->degree));
    cache->frames[key] = fr;
    return fr;
}

/* Scalar basis of P^degree on a region (see cut_basis_kind) */
template<CutMesh Mesh, typename VT = typename Mesh::coordinate_type>
class cut_scalar_basis {

    using coordinate_type = typename Mesh::coordinate_type;
    using point_type      = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    /* scaled monomials */
    point_type      cell_bar;
    coordinate_type cell_h = 0.0;
    /* orthonormal basis */
    std::shared_ptr<const cut_basis_frame<coordinate_type, DIM>> frame;

    size_t          basis_degree, basis_size;

    void
    init_monomial(const Mesh& msh, const typename Mesh::cell_type& cl, location where, bool subcell_centering)
    {
        if (!is_cut(msh, cl) or !subcell_centering) {
            cell_bar = barycenter(msh, cl);
            cell_h   = diameter(msh, cl);
            return;
        }

        const auto& ccd = cut_info(msh, cl);
        bool ill_cut_side = (ccd.cut == cut_type::T_KO_NEG and where == location::IN_NEGATIVE_SIDE)
                         or (ccd.cut == cut_type::T_KO_POS and where == location::IN_POSITIVE_SIDE);

        if (ill_cut_side and ccd.stabilizing_cell.has_value()) {
            auto cl_stab = msh[ccd.stabilizing_cell.value()];
            cell_bar = 0.5 * (barycenter(msh, cl, where) + barycenter(msh, cl_stab, where));
            cell_h   = diameter(msh, cl, where) + diameter(msh, cl_stab, where);
        }
        else {
            cell_bar = barycenter(msh, cl, where);
            cell_h   = diameter(msh, cl, where);
        }
    }

public:
    cut_scalar_basis(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree,
                     location where, const cut_basis_options<coordinate_type, DIM>& opt)
        : basis_degree(degree), basis_size(size(degree))
    {
        if (opt.kind == cut_basis_kind::orthonormal)
            frame = cut_cell_basis_frame(msh, cl, where, degree, opt.cache.get());
        else
            init_monomial(msh, cl, where, opt.kind == cut_basis_kind::subcell_monomial);
    }

    /* Scaled monomials (ProtoN) */
    cut_scalar_basis(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree,
                     location where, bool subcell_centering = true)
        : basis_degree(degree), basis_size(size(degree))
    {
        init_monomial(msh, cl, where, subcell_centering);
    }

    /* Scaled monomials centered at bar and scaled by h/2 (e.g. agglomerated cells) */
    cut_scalar_basis(const point_type& bar, coordinate_type h, size_t degree)
        : cell_bar(bar), cell_h(h), basis_degree(degree), basis_size(size(degree))
    {}

    /* Orthonormal basis of a given frame (e.g. agglomerated cells) */
    cut_scalar_basis(std::shared_ptr<const cut_basis_frame<coordinate_type, DIM>> fr, size_t degree)
        : frame(fr), basis_degree(degree), basis_size(size(degree))
    {
        assert(fr->degree >= degree);
    }

    Eigen::Matrix<VT, Eigen::Dynamic, 1>
    eval_basis(const point_type& pt) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, 1> ret(basis_size);

        if (frame) {
            priv::legendre_products(*frame, pt, basis_degree, &ret, (Eigen::Matrix<VT, Eigen::Dynamic, DIM>*)nullptr);
            if (frame->transform.size() != 0)
                ret = frame->transform.topLeftCorner(basis_size, basis_size).template triangularView<Eigen::Lower>() * ret;
            return ret;
        }

        std::array<coordinate_type, DIM> b;
        for (size_t d = 0; d < DIM; d++)
            b[d] = (pt[d] - cell_bar[d]) / (0.5*cell_h);

        const auto& expo = priv::cut_monomial_exponents<DIM>(basis_degree);
        for (size_t pos = 0; pos < basis_size; pos++) {
            VT v = 1;
            for (size_t d = 0; d < DIM; d++)
                v *= iexp_pow(b[d], expo[pos][d]);
            ret(pos) = v;
        }

        return ret;
    }

    Eigen::Matrix<VT, Eigen::Dynamic, DIM>
    eval_gradients(const point_type& pt) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, DIM> ret(basis_size, DIM);

        if (frame) {
            priv::legendre_products(*frame, pt, basis_degree, (Eigen::Matrix<VT, Eigen::Dynamic, 1>*)nullptr, &ret);
            if (frame->transform.size() != 0)
                ret = frame->transform.topLeftCorner(basis_size, basis_size).template triangularView<Eigen::Lower>() * ret;
            return ret;
        }

        std::array<coordinate_type, DIM> b;
        for (size_t d = 0; d < DIM; d++)
            b[d] = (pt[d] - cell_bar[d]) / (0.5*cell_h);
        auto ih = 2.0/cell_h;

        const auto& expo = priv::cut_monomial_exponents<DIM>(basis_degree);
        for (size_t pos = 0; pos < basis_size; pos++) {
            for (size_t d = 0; d < DIM; d++) {
                auto p = expo[pos][d];
                VT v = (p == 0) ? 0 : p*ih*iexp_pow(b[d], p-1);
                for (size_t e = 0; e < DIM; e++)
                    if (e != d)
                        v *= iexp_pow(b[e], expo[pos][e]);
                ret(pos, d) = v;
            }
        }

        return ret;
    }

    /* Values at the points (one row per point) */
    template<typename Points>
    Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>
    eval_basis_batch(const Points& pts) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic> ret(pts.size(), basis_size);
        if (frame) {
            Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic> P;
            priv::legendre_products_batch(*frame, pts, basis_degree, &P,
                (std::array<Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>, DIM>*)nullptr);
            if (frame->transform.size() == 0)
                return P;
            ret.noalias() = P * frame->transform.topLeftCorner(basis_size, basis_size)
                                     .template triangularView<Eigen::Lower>().transpose();
            return ret;
        }
        for (size_t q = 0; q < pts.size(); q++)
            ret.row(q) = eval_basis(pts[q]).transpose();
        return ret;
    }

    /* Derivatives at the points: one matrix per direction, one row per point */
    template<typename Points>
    std::array<Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>, DIM>
    eval_gradients_batch(const Points& pts) const
    {
        std::array<Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>, DIM> ret;
        if (frame) {
            std::array<Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>, DIM> D;
            priv::legendre_products_batch(*frame, pts, basis_degree,
                (Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>*)nullptr, &D);
            if (frame->transform.size() == 0)
                return D;
            for (size_t d = 0; d < DIM; d++)
                ret[d].noalias() = D[d] * frame->transform.topLeftCorner(basis_size, basis_size)
                                              .template triangularView<Eigen::Lower>().transpose();
            return ret;
        }
        for (size_t d = 0; d < DIM; d++)
            ret[d].resize(pts.size(), basis_size);
        for (size_t q = 0; q < pts.size(); q++) {
            auto g = eval_gradients(pts[q]);
            for (size_t d = 0; d < DIM; d++)
                ret[d].row(q) = g.col(d).transpose();
        }
        return ret;
    }

    size_t size() const { return basis_size; }
    size_t degree() const { return basis_degree; }

    static size_t size(size_t degree) {
        return cut_basis_frame_size(degree, DIM);
    }
};

/* Vector-valued version of cut_scalar_basis, the components of each scalar
 * function being interleaved. */
template<CutMesh Mesh, typename VT = typename Mesh::coordinate_type>
class cut_vector_basis {

    using point_type = typename Mesh::point_type;
    using coordinate_type = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

    cut_scalar_basis<Mesh, VT> scalar_basis;

public:
    cut_vector_basis(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree,
                     location where, const cut_basis_options<coordinate_type, DIM>& opt)
        : scalar_basis(msh, cl, degree, where, opt)
    {}

    cut_vector_basis(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree,
                     location where, bool subcell_centering = true)
        : scalar_basis(msh, cl, degree, where, subcell_centering)
    {}

    cut_vector_basis(const point_type& bar, coordinate_type h, size_t degree)
        : scalar_basis(bar, h, degree)
    {}

    cut_vector_basis(std::shared_ptr<const cut_basis_frame<coordinate_type, DIM>> fr, size_t degree)
        : scalar_basis(fr, degree)
    {}

    Eigen::Matrix<VT, Eigen::Dynamic, DIM>
    eval_basis(const point_type& pt) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, DIM> ret = Eigen::Matrix<VT, Eigen::Dynamic, DIM>::Zero(size(), DIM);

        const auto phi = scalar_basis.eval_basis(pt);
        for (size_t i = 0; i < scalar_basis.size(); i++)
            for (size_t d = 0; d < DIM; d++)
                ret(DIM*i + d, d) = phi(i);

        return ret;
    }

    size_t size() const { return DIM*scalar_basis.size(); }
    size_t degree() const { return scalar_basis.degree(); }

    /* the scalar basis whose functions are the components (DIM i + d, d) */
    const cut_scalar_basis<Mesh, VT>& scalar_part() const { return scalar_basis; }

    static size_t size(size_t degree) {
        return DIM*cut_scalar_basis<Mesh, VT>::size(degree);
    }
};

/* Frame of the orthonormal basis of the part "where" of a face of a 3D mesh
 * (the whole face if it is not cut), in the coordinates of its plane */
template<CutMesh Mesh>
std::shared_ptr<const cut_face_frame<typename Mesh::coordinate_type>>
make_cut_face_frame(const Mesh& msh, const typename Mesh::face_type& fc, location where, size_t degree)
    requires (Mesh::dimension == 3)
{
    using T = typename Mesh::coordinate_type;
    auto ff = std::make_shared<cut_face_frame<T>>();
    auto loc = locate(msh, fc);
    bool part = (loc == location::ON_INTERFACE);
    if (part) {
        const auto& fd = cut_info(msh, fc);
        ff->origin = fd.origin;
        ff->U = fd.U;
    }
    else
        cutgeo::polygon_frame(points(msh, fc), ff->origin, ff->U);
    cutgeo::planar_level_set<T,3> plane(*msh.level_set_function, ff->origin, ff->U);

    auto qps3 = part ? integrate(msh, fc, std::max<size_t>(2, 2*degree), where)
                     : integrate(msh, fc, std::max<size_t>(2, 2*degree));
    std::vector<point<T,3>> pts3;
    if (part)
        pts3 = face_region_points(msh, fc, where);
    else {
        auto fp = points(msh, fc);
        pts3.assign(fp.begin(), fp.end());
    }
    std::vector<quadrature_point<T,2>> qps;
    std::vector<point<T,2>> pts;
    for (auto& qp : qps3)
        qps.push_back(make_qp(plane.to_local(qp.point()), qp.weight()));
    for (auto& p : pts3)
        pts.push_back(plane.to_local(p));

    auto box = part ? nullptr : priv::box_frame(pts, degree);
    ff->frame = box ? *box : *priv::make_frame(qps, pts, degree);
    return ff;
}

/* L2-orthonormal basis on the part of a face lying in the subdomain "where"
 * (on the whole face if the face is not cut): Legendre polynomials on the
 * segment (2D), orthonormalized Legendre products in the plane of the face
 * (3D). */
template<CutMesh Mesh, typename VT = typename Mesh::coordinate_type>
class cut_face_basis {

    using coordinate_type = typename Mesh::coordinate_type;
    using point_type      = typename Mesh::point_type;
    static const size_t DIM = Mesh::dimension;

    /* 2D */
    point_type      face_bar;
    point_type      base;
    coordinate_type face_h = 0;
    /* 3D */
    std::shared_ptr<const cut_face_frame<coordinate_type>> fframe;
    coordinate_type scale = 1; // 1/sqrt(|F|)

    size_t          basis_degree, basis_size;

    void
    init(const Mesh& msh, const typename Mesh::face_type& fc, location where,
         cut_basis_cache<coordinate_type, DIM>* cache)
    {
        auto loc = locate(msh, fc);
        bool part = (loc == where || loc == location::ON_INTERFACE);
        if constexpr (DIM == 2) {
            if (!part) {
                face_bar = barycenter(msh, fc);
                face_h   = diameter(msh, fc);
                base     = face_bar - points(msh, fc)[0];
            }
            else {
                face_bar = barycenter(msh, fc, where);
                face_h   = diameter(msh, fc, where);
                base     = face_bar - points(msh, fc, where)[0];
            }
        }
        else {
            if (!part)
                where = loc;
            bool cut = (loc == location::ON_INTERFACE);
            size_t fc_id = offset(msh, fc);
            size_t key = 2*fc_id + (where == location::IN_POSITIVE_SIDE ? 1 : 0);
            size_t deg = cache ? std::max(basis_degree, cache->face_degree) : basis_degree;
            if (cache && cut) {
                auto it = cache->face_frames.find(key);
                if (it != cache->face_frames.end() && it->second->frame.degree >= basis_degree)
                    fframe = it->second;
                else
                    fframe = cache->face_frames[key] = make_cut_face_frame(msh, fc, where, deg);
            }
            else if (cache) {
                auto f = cache->uncut_face_frames.find(fc_id);
                if (f && (*f)->frame.degree >= basis_degree)
                    fframe = *f;
                else {
                    fframe = make_cut_face_frame(msh, fc, where, deg);
                    cache->uncut_face_frames.insert(fc_id, fframe);
                }
            }
            else
                fframe = make_cut_face_frame(msh, fc, where, basis_degree);

            coordinate_type meas = 0;
            for (auto& qp : integrate(msh, fc, 0, where))
                meas += qp.weight();
            if (!cut)
                meas = measure(msh, fc);
            scale = 1.0 / std::sqrt(meas);
        }
    }

public:
    cut_face_basis(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree, location where)
        : basis_degree(degree), basis_size(size(degree))
    {
        init(msh, fc, where, nullptr);
    }

    cut_face_basis(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree, location where,
                   const cut_basis_options<coordinate_type, DIM>& opt)
        : basis_degree(degree), basis_size(size(degree))
    {
        init(msh, fc, where, opt.cache.get());
    }

    Eigen::Matrix<VT, Eigen::Dynamic, 1>
    eval_basis(const point_type& pt) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, 1> ret(basis_size);

        if constexpr (DIM == 2) {
            auto v = base.to_vector();
            auto t = (pt - face_bar).to_vector();
            auto ep = 4.0*v.dot(t)/(face_h*face_h);
            auto coeff = std::sqrt(face_h / 2.0);

            /* Legendre polynomials, three-term recurrence */
            coordinate_type p_prev = 0.0, p_curr = 1.0;
            for (size_t n = 0; n <= basis_degree; n++) {
                ret(n) = p_curr * std::sqrt((2.0*n + 1.0) / 2.0) / coeff;
                auto p_next = ((2.0*n + 1.0) * ep * p_curr - n * p_prev) / (n + 1.0);
                p_prev = p_curr;
                p_curr = p_next;
            }
        }
        else {
            Eigen::Matrix<coordinate_type, 2, 1> u = fframe->U.transpose() * (pt - fframe->origin).to_vector();
            const auto& fr = fframe->frame;
            priv::legendre_products(fr, point<coordinate_type,2>(u(0), u(1)), basis_degree, &ret,
                                    (Eigen::Matrix<VT, Eigen::Dynamic, 2>*)nullptr);
            if (fr.transform.size() != 0)
                ret = fr.transform.topLeftCorner(basis_size, basis_size).template triangularView<Eigen::Lower>() * ret;
            ret *= scale;
        }

        return ret;
    }

    /* Values at the points (one row per point) */
    template<typename Points>
    Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>
    eval_basis_batch(const Points& pts) const
    {
        Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic> ret(pts.size(), basis_size);
        if constexpr (DIM == 3) {
            std::vector<point<coordinate_type,2>> u(pts.size());
            for (size_t q = 0; q < pts.size(); q++) {
                Eigen::Matrix<coordinate_type, 2, 1> uq = fframe->U.transpose() * (pts[q] - fframe->origin).to_vector();
                u[q] = point<coordinate_type,2>(uq(0), uq(1));
            }
            const auto& fr = fframe->frame;
            Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic> P;
            priv::legendre_products_batch(fr, u, basis_degree, &P,
                (std::array<Eigen::Matrix<VT, Eigen::Dynamic, Eigen::Dynamic>, 2>*)nullptr);
            if (fr.transform.size() == 0)
                ret = P * scale;
            else
                ret.noalias() = scale * P * fr.transform.topLeftCorner(basis_size, basis_size)
                                                .template triangularView<Eigen::Lower>().transpose();
            return ret;
        }
        for (size_t q = 0; q < pts.size(); q++)
            ret.row(q) = eval_basis(pts[q]).transpose();
        return ret;
    }

    size_t size() const { return basis_size; }
    size_t degree() const { return basis_degree; }

    static size_t size(size_t degree) {
        return cut_basis_frame_size(degree, DIM-1);
    }
};

} // namespace disk
