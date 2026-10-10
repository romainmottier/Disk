/*
 * Local operators of the unfitted HHO method stabilized by polynomial
 * extension (Burman, Ern, Mottier), ported from ProtoN.
 *
 * Local unknowns of a subcell P = (T, i, N_i^{-1}(T)), see cut_subcell:
 *  - T uncut:  [ u_T | u_dT ]
 *  - T cut:    [ u_T1 | u_T2 | u_(dT)1 | u_(dT)2 ]
 *  - followed, for each dependent ill-cut cell S, by
 *              [ u_S1 | u_S2 | u_(dS)1 | u_(dS)2 ]
 * Faces are numbered as in faces(msh, cl). Negative side = Omega_1.
 */

#pragma once

#include <array>
#include <map>
#include <memory>
#include <Eigen/Dense>

#include "diskpp/methods/hho"
#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/bases/cut_bases.hpp"
#include "quadrature_matrices.hpp"

template<typename T, size_t DIM = 2>
struct unfitted_hho_parameters {
    disk::hho_degree_info hdi;
    T    eta_gamma;          // penalty of s^Gamma and of the Dirichlet jump
    T    eta_stab;           // scaling of s°
    T    eta_N;              // scaling of s^N
    disk::cut_basis_options<T, DIM> basis; // cell and face bases (orthonormal by default) and their cache

    /* Local matrix and gradient operator of the uncut cells without dependent
     * ill-cut cells: they only depend on the shape of the cell (up to a
     * translation), on the order of its faces and on kappa, hence are
     * computed once per shape on a structured mesh (shared by the copies of
     * the parameters; null: no caching). */
    using uncut_key = std::vector<long long>;
    using local_operators = std::pair<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>,
                                      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
    std::shared_ptr<std::map<uncut_key, local_operators>> uncut_cache;

    /* s° is scaled by the cell degree, as in ProtoN. s^N is scaled by
     * 20 (k+1)^2 (ProtoN: (k+1)^2): this law keeps the errors of Theorem 4.6
     * and the condition number in the worst cut configuration within a factor
     * 2 of the best eta_N for k = 0..3 (see the README). */
    explicit unfitted_hho_parameters(size_t k)
        : hdi(k+1, k), eta_gamma(1.0), eta_stab(k+1), eta_N(20.0*(k+1)*(k+1))
    {}
};

template<typename Mesh>
using mesh_parameters = unfitted_hho_parameters<typename Mesh::coordinate_type, Mesh::dimension>;

/* Key of the local operators of an uncut cell: the positions of its vertices
 * and of the vertices of its faces (in their local order) relative to its
 * first vertex, at the precision 1e-9 h, and kappa. Two cells with the same
 * key are translates of each other with the same numbering, hence have the
 * same local operators (the bases only depend on the shape of the cell). */
template<disk::CutMesh Mesh>
bool
uncut_cell_key(const Mesh& msh, const typename Mesh::cell_type& cl, typename Mesh::coordinate_type kappa,
               std::vector<long long>& key)
{
    auto pts = points(msh, cl);
    auto h = diameter(msh, cl);
    key.clear();
    auto add = [&](const typename Mesh::point_type& p) {
        for (size_t d = 0; d < Mesh::dimension; d++)
            key.push_back(std::llround((p[d] - pts[0][d]) / h * 1e9));
    };
    key.push_back(std::llround(h * 1e12));
    for (auto& p : pts)
        add(p);
    for (auto& fc : faces(msh, cl)) {
        auto fp = points(msh, fc);
        key.push_back(-1);
        for (auto& p : fp)
            add(p);
    }
    key.push_back(std::llround(kappa * 1e6));
    return true;
}

template<disk::CutMesh Mesh>
class subcell_layout {

    using T = typename Mesh::coordinate_type;

public:
    size_t cbs, fbs, gbs, num_faces;
    bool   cut;
    size_t current_dofs, local_dofs;
    std::vector<size_t> dependent_faces;   // number of faces of each dependent cell
    std::vector<size_t> dependent_offsets; // offset of the block of each dependent cell

    subcell_layout(const Mesh& msh, const disk::cut_subcell& P, const disk::hho_degree_info& hdi)
    {
        auto cl   = msh[P.cell];
        cbs       = disk::cut_scalar_basis<Mesh>::size(hdi.cell_degree());
        fbs       = disk::cut_face_basis<Mesh>::size(hdi.face_degree());
        gbs       = disk::cut_vector_basis<Mesh>::size(hdi.grad_degree());
        num_faces = faces(msh, cl).size();
        cut       = disk::is_cut(msh, cl);

        current_dofs = cut ? 2*(cbs + num_faces*fbs) : cbs + num_faces*fbs;
        local_dofs   = current_dofs;
        for (auto dp : P.dependent_cells) {
            dependent_faces.push_back(faces(msh, msh[dp]).size());
            dependent_offsets.push_back(local_dofs);
            local_dofs += extended_dofs(dependent_faces.size() - 1);
        }
    }

    /* Size of the block of the dependent cell j (a cut cell) */
    size_t extended_dofs(size_t j) const {
        return 2*(cbs + dependent_faces[j]*fbs);
    }

    /* Offset of the cell unknowns of side "side" in a block of a cut cell */
    size_t cut_cell_offset(disk::location side) const {
        return (side == disk::location::IN_NEGATIVE_SIDE) ? 0 : cbs;
    }

    /* Offset of the face unknowns of side "side" in a block of a cut cell
     * with nf faces */
    size_t cut_face_offset(disk::location side, size_t nf) const {
        return (side == disk::location::IN_NEGATIVE_SIDE) ? 2*cbs : 2*cbs + nf*fbs;
    }

    size_t cut_face_offset(disk::location side) const {
        return cut_face_offset(side, num_faces);
    }

    /* Offsets in the current block */
    size_t cell_offset(disk::location side) const {
        return cut ? cut_cell_offset(side) : 0;
    }

    size_t face_offset(disk::location side) const {
        return cut ? cut_face_offset(side) : cbs;
    }

    size_t dependent_offset(size_t j) const {
        return dependent_offsets[j];
    }
};

/* Scatter a block of a one-side operator (columns [cell | faces]) into a
 * block of a cut cell with nf faces (columns [cell_1 | cell_2 | faces_1 | faces_2]). */
template<typename Mat, typename Layout>
void
place_one_side(Mat& dst, const Mat& src, const Layout& lay, disk::location side, size_t nf)
{
    auto rows = src.rows();
    dst.block(0, lay.cut_cell_offset(side), rows, lay.cbs) += src.block(0, 0, rows, lay.cbs);
    dst.block(0, lay.cut_face_offset(side, nf), rows, nf*lay.fbs)
        += src.block(0, lay.cbs, rows, nf*lay.fbs);
}

template<typename Mat, typename Layout>
void
place_one_side(Mat& dst, const Mat& src, const Layout& lay, disk::location side)
{
    place_one_side(dst, src, lay, side, lay.num_faces);
}

template<typename T>
T side_coeff(disk::location side) {
    /* The interface terms are only present on the side of the smallest
     * diffusivity, i.e. Omega_1 (kappa_1 <= kappa_2) */
    return (side == disk::location::IN_NEGATIVE_SIDE) ? 1.0 : 0.0;
}

/***************************************************************************/
/*                         GRADIENT RECONSTRUCTION                          */
/***************************************************************************/

/* Contribution of the unknowns of the cell T to the rhs of the reconstruction
 * of the gradient on (T, side), eq. (24): one block of current_dofs columns. */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
gradrec_local_contribution(const Mesh& msh, const typename Mesh::cell_type& cl,
                           disk::location side, const subcell_layout<Mesh>& lay,
                           const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    static const size_t DIM = Mesh::dimension;

    const auto celdeg  = prm.hdi.cell_degree();
    const auto facdeg  = prm.hdi.face_degree();
    const auto graddeg = prm.hdi.grad_degree();
    const auto cbs = lay.cbs, fbs = lay.fbs, gbs = lay.gbs;
    const auto num_faces = lay.num_faces;

    disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, side, prm.basis);
    disk::cut_vector_basis<Mesh> gb(msh, cl, graddeg, side, prm.basis);
    const auto& sgb = gb.scalar_part();

    matrix_type rhs_tmp = matrix_type::Zero(gbs, cbs + num_faces*fbs);

    // (grad u_T, q)_T
    {
        const auto qps = integrate(msh, cl, celdeg - 1 + facdeg, side);
        auto w  = quadrature_weights(qps);
        auto G  = basis_values(sgb, qps);
        auto dC = basis_gradients<DIM>(cb, qps);
        for (size_t d = 0; d < DIM; d++)
            add_component_rows<DIM>(rhs_tmp, 0, weighted_product(G, w, dC[d]), d);
    }

    // (u_dT - u_T, q.n_T)_dT
    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++) {
        const auto fc = fcs[i];
        const Eigen::Matrix<T, Mesh::dimension, 1> n = normal(msh, cl, fc);
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, side, prm.basis);
        const auto qps_f = integrate(msh, fc, facdeg + std::max(facdeg, celdeg), side);
        if (qps_f.empty())
            continue;
        auto w  = quadrature_weights(qps_f);
        auto G  = basis_values(sgb, qps_f);
        matrix_type GF = weighted_product(G, w, basis_values(fb, qps_f));
        matrix_type GC = weighted_product(G, w, basis_values(cb, qps_f));
        for (size_t d = 0; d < DIM; d++) {
            add_component_rows<DIM>(rhs_tmp, cbs + i*fbs, GF, d, T(n(d)));
            add_component_rows<DIM>(rhs_tmp, 0, GC, d, T(-n(d)));
        }
    }

    if (!disk::is_cut(msh, cl))
        return rhs_tmp;

    matrix_type gr_rhs = matrix_type::Zero(gbs, lay.current_dofs);
    place_one_side(gr_rhs, rhs_tmp, lay, side);

    // - delta_i1 ([u_T]_Gamma, q.n_Gamma)_TGamma
    if (side_coeff<T>(side) != 0) {
        disk::cut_scalar_basis<Mesh> cb_neg(msh, cl, celdeg, disk::location::IN_NEGATIVE_SIDE, prm.basis);
        disk::cut_scalar_basis<Mesh> cb_pos(msh, cl, celdeg, disk::location::IN_POSITIVE_SIDE, prm.basis);
        matrix_type interface_term = matrix_type::Zero(gbs, 2*cbs);
        auto dn   = disk::interface_normals(msh, cl, celdeg + graddeg);
        auto iqps = disk::integrate_interface(msh, cl, celdeg + graddeg, disk::location::IN_NEGATIVE_SIDE);
        auto w  = quadrature_weights(iqps);
        auto G  = basis_values(sgb, iqps);
        auto Cn = basis_values(cb_neg, iqps);
        auto Cp = basis_values(cb_pos, iqps);
        for (size_t d = 0; d < DIM; d++) {
            qm_vector<T> wd(w.size());
            for (Eigen::Index q = 0; q < w.size(); q++)
                wd(q) = w(q) * dn[q](d);
            add_component_rows<DIM>(interface_term, 0, weighted_product(G, wd, Cn), d, T(-1));
            add_component_rows<DIM>(interface_term, cbs, weighted_product(G, wd, Cp), d, T(1));
        }
        gr_rhs.block(0, 0, gbs, 2*cbs) += side_coeff<T>(side) * interface_term;
    }

    return gr_rhs;
}

/* Contribution of the unknowns of a dependent ill-cut cell S to the rhs of the
 * reconstruction of the gradient on (T, side): the test functions of T^i are
 * extended to S, eq. (24). One block of extended_dofs columns. */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
gradrec_extended_contribution(const Mesh& msh, const typename Mesh::cell_type& cl,
                              const typename Mesh::cell_type& dp_cl,
                              disk::location side, const subcell_layout<Mesh>& lay,
                              const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    static const size_t DIM = Mesh::dimension;

    const auto celdeg  = prm.hdi.cell_degree();
    const auto facdeg  = prm.hdi.face_degree();
    const auto graddeg = prm.hdi.grad_degree();
    const auto cbs = lay.cbs, fbs = lay.fbs, gbs = lay.gbs;
    const auto fcs = faces(msh, dp_cl);
    const size_t num_faces = fcs.size();

    disk::cut_scalar_basis<Mesh> cb_dp(msh, dp_cl, celdeg, side, prm.basis);
    disk::cut_vector_basis<Mesh> gb(msh, cl, graddeg, side, prm.basis);
    const auto& sgb = gb.scalar_part();

    matrix_type rhs_tmp = matrix_type::Zero(gbs, cbs + num_faces*fbs);
    matrix_type gr_rhs  = matrix_type::Zero(gbs, 2*(cbs + num_faces*fbs));

    // (u_(dS)i - u_Si, q+.n_S)_(dS)i
    for (size_t i = 0; i < fcs.size(); i++) {
        const auto fc = fcs[i];
        const Eigen::Matrix<T, Mesh::dimension, 1> n = normal(msh, dp_cl, fc);
        disk::cut_face_basis<Mesh> fb_dp(msh, fc, facdeg, side, prm.basis);
        const auto qps_f = integrate(msh, fc, facdeg + std::max(facdeg, celdeg), side);
        if (qps_f.empty())
            continue;
        auto w  = quadrature_weights(qps_f);
        auto G  = basis_values(sgb, qps_f);
        matrix_type GF = weighted_product(G, w, basis_values(fb_dp, qps_f));
        matrix_type GC = weighted_product(G, w, basis_values(cb_dp, qps_f));
        for (size_t d = 0; d < DIM; d++) {
            add_component_rows<DIM>(rhs_tmp, 0, GC, d, T(-n(d)));
            add_component_rows<DIM>(rhs_tmp, cbs + i*fbs, GF, d, T(n(d)));
        }
    }
    place_one_side(gr_rhs, rhs_tmp, lay, side, num_faces);

    // - delta_i1 ([u_S]_Gamma, q+.n_Gamma)_SGamma
    if (side_coeff<T>(side) != 0) {
        disk::cut_scalar_basis<Mesh> cb_dp_neg(msh, dp_cl, celdeg, disk::location::IN_NEGATIVE_SIDE, prm.basis);
        disk::cut_scalar_basis<Mesh> cb_dp_pos(msh, dp_cl, celdeg, disk::location::IN_POSITIVE_SIDE, prm.basis);
        matrix_type interface_term = matrix_type::Zero(gbs, 2*cbs);
        auto dn   = disk::interface_normals(msh, dp_cl, celdeg + graddeg);
        auto iqps = disk::integrate_interface(msh, dp_cl, celdeg + graddeg, disk::location::IN_NEGATIVE_SIDE);
        auto w  = quadrature_weights(iqps);
        auto G  = basis_values(sgb, iqps);
        auto Cn = basis_values(cb_dp_neg, iqps);
        auto Cp = basis_values(cb_dp_pos, iqps);
        for (size_t d = 0; d < DIM; d++) {
            qm_vector<T> wd(w.size());
            for (Eigen::Index q = 0; q < w.size(); q++)
                wd(q) = w(q) * dn[q](d);
            add_component_rows<DIM>(interface_term, 0, weighted_product(G, wd, Cn), d, T(-1));
            add_component_rows<DIM>(interface_term, cbs, weighted_product(G, wd, Cp), d, T(1));
        }
        gr_rhs.block(0, 0, gbs, 2*cbs) += side_coeff<T>(side) * interface_term;
    }

    return gr_rhs;
}

/* Gradient reconstruction on a subcell of P^OK, eq. (24).
 * Returns (operator, G^T M G). */
template<disk::CutMesh Mesh>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>>
make_gradrec_POK(const Mesh& msh, const disk::cut_subcell& P,
                 const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto cl  = msh[P.cell];
    auto loc = P.side;
    subcell_layout<Mesh> lay(msh, P, prm.hdi);

    const auto celdeg  = prm.hdi.cell_degree();
    const auto graddeg = prm.hdi.grad_degree();

    disk::cut_vector_basis<Mesh> gb(msh, cl, graddeg, loc, prm.basis);

    /* mass matrix of the scalar basis of the gradient (the components are
     * independent) */
    const auto qps = integrate(msh, cl, celdeg - 1 + graddeg, loc);
    auto G = basis_values(gb.scalar_part(), qps);
    matrix_type M = weighted_product(G, quadrature_weights(qps), G);

    matrix_type gr_rhs = matrix_type::Zero(lay.gbs, lay.local_dofs);
    gr_rhs.block(0, 0, lay.gbs, lay.current_dofs) += gradrec_local_contribution(msh, cl, loc, lay, prm);
    for (size_t j = 0; j < P.dependent_cells.size(); j++) {
        auto dp_cl = msh[P.dependent_cells[j]];
        gr_rhs.block(0, lay.dependent_offset(j), lay.gbs, lay.extended_dofs(j))
            += gradrec_extended_contribution(msh, cl, dp_cl, loc, lay, prm);
    }

    matrix_type oper = solve_vector_mass<Mesh::dimension>(M, gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

/* Gradient reconstruction on a subcell of P^KO, eq. (25): G = grad u_Ti.
 * With the orthonormal bases, (G, G)_Ti is computed directly as the stiffness
 * matrix of u_Ti: the basis of an ill-cut side is orthonormal on the subcell
 * merged with its stabilizing cell, and its mass matrix on the (possibly tiny)
 * subcell alone would be ill-conditioned. The operator itself is only used by
 * the right-hand side of P^OK (lifting of the Dirichlet jump). */
template<disk::CutMesh Mesh>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>>
make_gradrec_PKO(const Mesh& msh, const disk::cut_subcell& P,
                 const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto cl  = msh[P.cell];
    auto loc = P.side;
    subcell_layout<Mesh> lay(msh, P, prm.hdi);
    assert(lay.cut);

    const auto celdeg  = prm.hdi.cell_degree();
    const auto graddeg = prm.hdi.grad_degree();

    disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, loc, prm.basis);

    static const size_t DIM = Mesh::dimension;

    if (prm.basis.kind == disk::cut_basis_kind::orthonormal) {
        matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);
        const auto qps = integrate(msh, cl, 2*celdeg - 2, loc);
        auto w  = quadrature_weights(qps);
        auto dC = basis_gradients<DIM>(cb, qps);
        auto co = lay.cell_offset(loc);
        for (size_t d = 0; d < DIM; d++)
            data.block(co, co, lay.cbs, lay.cbs) += weighted_product(dC[d], w, dC[d]);
        return std::make_pair(matrix_type(), data);
    }

    disk::cut_vector_basis<Mesh> gb(msh, cl, graddeg, loc, prm.basis);

    matrix_type gr_rhs = matrix_type::Zero(lay.gbs, lay.local_dofs);

    const auto qps = integrate(msh, cl, celdeg - 1 + graddeg, loc);
    auto w  = quadrature_weights(qps);
    auto G  = basis_values(gb.scalar_part(), qps);
    auto dC = basis_gradients<DIM>(cb, qps);
    matrix_type M = weighted_product(G, w, G);
    for (size_t d = 0; d < DIM; d++)
        add_component_rows<DIM>(gr_rhs, lay.cell_offset(loc), weighted_product(G, w, dC[d]), d);

    matrix_type oper = solve_vector_mass<DIM>(M, gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

/***************************************************************************/
/*                              STABILIZATIONS                              */
/***************************************************************************/

/* Lehrenfeld-Schöberl stabilization on one side of a cell, with columns
 * [cell | faces] (eq. (28), weight 1/h_T) */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_one_side_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, disk::location where,
                            const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto celdeg = prm.hdi.cell_degree();
    auto facdeg = prm.hdi.face_degree();
    auto cbs = disk::cut_scalar_basis<Mesh>::size(celdeg);
    auto fbs = disk::cut_face_basis<Mesh>::size(facdeg);
    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();
    matrix_type data = matrix_type::Zero(cbs + num_faces*fbs, cbs + num_faces*fbs);
    matrix_type If   = matrix_type::Identity(fbs, fbs);

    disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, where, prm.basis);

    auto hT = diameter(msh, cl);
    for (size_t i = 0; i < num_faces; i++) {
        auto fc = fcs[i];
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, where, prm.basis);

        matrix_type oper  = matrix_type::Zero(fbs, cbs + num_faces*fbs);
        matrix_type mass  = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, cbs + i*fbs, fbs, fbs) = -If;
        auto qps = integrate(msh, fc, facdeg + celdeg, where);
        if (qps.size() == 0) /* Avoid to invert a zero matrix */
            continue;

        auto w = quadrature_weights(qps);
        auto F = basis_values(fb, qps);
        mass  = weighted_product(F, w, F);
        trace = weighted_product(F, w, basis_values(cb, qps));

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * mass * oper * (1./hT);
    }

    return data;
}

/* s° restricted to the subcell P, in the local numbering of P */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_stabilization(const Mesh& msh, const disk::cut_subcell& P,
                   const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto cl  = msh[P.cell];
    auto loc = P.side;
    subcell_layout<Mesh> lay(msh, P, prm.hdi);
    auto cbs = lay.cbs;
    auto nfb = lay.num_faces*lay.fbs;

    matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);
    const auto stab = make_one_side_stabilization(msh, cl, loc, prm);

    auto co = lay.cell_offset(loc);
    auto fo = lay.face_offset(loc);
    data.block(co, co, cbs, cbs) += stab.block(  0,   0, cbs, cbs);
    data.block(co, fo, cbs, nfb) += stab.block(  0, cbs, cbs, nfb);
    data.block(fo, co, nfb, cbs) += stab.block(cbs,   0, nfb, cbs);
    data.block(fo, fo, nfb, nfb) += stab.block(cbs, cbs, nfb, nfb);

    return data;
}

/* s^Gamma, eq. (29): assembled once per cut cell, from its negative side */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_interface_penalty(const Mesh& msh, const disk::cut_subcell& P,
                       const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto cl = msh[P.cell];
    subcell_layout<Mesh> lay(msh, P, prm.hdi);
    matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);

    if (!lay.cut || P.side != disk::location::IN_NEGATIVE_SIDE)
        return data;

    auto celdeg = prm.hdi.cell_degree();
    auto cbs = lay.cbs;
    disk::cut_scalar_basis<Mesh> cb_n(msh, cl, celdeg, disk::location::IN_NEGATIVE_SIDE, prm.basis);
    disk::cut_scalar_basis<Mesh> cb_p(msh, cl, celdeg, disk::location::IN_POSITIVE_SIDE, prm.basis);

    auto hT = diameter(msh, cl);
    auto iqps = disk::integrate_interface(msh, cl, 2*celdeg, disk::location::IN_NEGATIVE_SIDE);
    qm_vector<T> w = quadrature_weights(iqps) * (prm.eta_gamma / hT);
    auto Cn = basis_values(cb_n, iqps);
    auto Cp = basis_values(cb_p, iqps);
    matrix_type Mnp = weighted_product(Cn, w, Cp);
    data.block(  0,   0, cbs, cbs) += weighted_product(Cn, w, Cn);
    data.block(cbs, cbs, cbs, cbs) += weighted_product(Cp, w, Cp);
    data.block(  0, cbs, cbs, cbs) -= Mnp;
    data.block(cbs,   0, cbs, cbs) -= Mnp.transpose();

    return data;
}

/* s^N, eq. (30), unscaled (multiplied by eta_N by the caller) */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_extension_stabilization(const Mesh& msh, const disk::cut_subcell& P,
                             const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto cl  = msh[P.cell];
    auto loc = P.side;
    subcell_layout<Mesh> lay(msh, P, prm.hdi);
    auto cbs = lay.cbs;
    auto celdeg = prm.hdi.cell_degree();

    matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);

    if (P.dependent_cells.empty())
        return data;

    disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, loc, prm.basis);
    auto h = diameter(msh, cl);
    auto offset_cl = lay.cell_offset(loc);
    auto qps = integrate(msh, cl, 2*celdeg, loc);
    qm_vector<T> w = quadrature_weights(qps) / (h*h);
    auto C = basis_values(cb, qps);
    matrix_type MCC = weighted_product(C, w, C);

    for (size_t j = 0; j < P.dependent_cells.size(); j++) {
        auto dp_cl = msh[P.dependent_cells[j]];
        auto offset_dp = lay.dependent_offset(j) + lay.cut_cell_offset(loc);
        disk::cut_scalar_basis<Mesh> cb_dp(msh, dp_cl, celdeg, loc, prm.basis);

        /* (v_Si - v_Ti^+, w_Si - w_Ti^+)_Ti */
        auto D = basis_values(cb_dp, qps);
        matrix_type MCD = weighted_product(C, w, D);
        data.block(offset_cl, offset_cl, cbs, cbs) += MCC;
        data.block(offset_cl, offset_dp, cbs, cbs) -= MCD;
        data.block(offset_dp, offset_cl, cbs, cbs) -= MCD.transpose();
        data.block(offset_dp, offset_dp, cbs, cbs) += weighted_product(D, w, D);
    }

    return data;
}

/***************************************************************************/
/*                               RIGHT-HAND SIDE                            */
/***************************************************************************/

/* (f, phi)_{T^where} */
template<disk::CutMesh Mesh, typename Function>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>
make_cut_rhs(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, const Function& f,
             disk::location where, const disk::cut_basis_options<typename Mesh::coordinate_type, Mesh::dimension>& basis)
{
    using T = typename Mesh::coordinate_type;

    disk::cut_scalar_basis<Mesh> cb(msh, cl, degree, where, basis);
    auto qps = integrate(msh, cl, 2*degree, where);
    if (qps.empty())
        return Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(cb.size());
    qm_vector<T> wf(qps.size());
    for (size_t q = 0; q < qps.size(); q++)
        wf(q) = qps[q].weight() * f(qps[q].point());

    return basis_values(cb, qps).transpose() * wf;
}

/* Right-hand side of the subcell P, eq. (35). oper_gr is the gradient
 * reconstruction operator of P (used for the lifting of the Dirichlet jump). */
template<disk::CutMesh Mesh, typename TestCase>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>
make_rhs_jumps(const Mesh& msh, const disk::cut_subcell& P, bool POK,
               const Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>& oper_gr,
               const TestCase& test_case,
               const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto cl  = msh[P.cell];
    auto loc = P.side;
    subcell_layout<Mesh> lay(msh, P, prm.hdi);
    auto cbs = lay.cbs;
    auto celdeg  = prm.hdi.cell_degree();
    auto graddeg = prm.hdi.grad_degree();
    auto kappa_1 = test_case.parms.kappa_1;
    auto hT = diameter(msh, cl);

    vector_type f = vector_type::Zero(lay.local_dofs);

    // (f, w_Ti)_Ti
    f.segment(lay.cell_offset(loc), cbs) += make_cut_rhs(msh, cl, celdeg, test_case.rhs_fun, loc, prm.basis);

    if (lay.cut) {
        disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, loc, prm.basis);
        auto iqps = disk::integrate_interface(msh, cl, 2*celdeg, disk::location::IN_NEGATIVE_SIDE);
        for (auto& qp : iqps) {
            auto phi = cb.eval_basis(qp.point());
            auto gD  = test_case.dirichlet_jump(qp.point());
            // kappa_1 h^-1 (g_D, [w_T]_Gamma)_TGamma
            if (loc == disk::location::IN_NEGATIVE_SIDE)
                f.segment(0, cbs) += kappa_1 * qp.weight() * gD * phi * prm.eta_gamma / hT;
            else {
                f.segment(cbs, cbs) -= kappa_1 * qp.weight() * gD * phi * prm.eta_gamma / hT;
                // (g_N, w_T2)_TGamma
                f.segment(cbs, cbs) += qp.weight() * test_case.neumann_jump(qp.point()) * phi;
            }
        }
    }

    // - kappa_1 (L(g_D), G(w))_T1 : lifting of the Dirichlet jump, eq. (26)
    if (POK && loc == disk::location::IN_NEGATIVE_SIDE) {
        disk::cut_vector_basis<Mesh> gb(msh, cl, graddeg, loc, prm.basis);
        vector_type F_bis = vector_type::Zero(lay.gbs);
        auto add_interface = [&](const typename Mesh::cell_type& c) {
            auto dn   = disk::interface_normals(msh, c, 2*graddeg);
            auto iqps = disk::integrate_interface(msh, c, 2*graddeg, disk::location::IN_NEGATIVE_SIDE);
            for (size_t q = 0; q < iqps.size(); q++) {
                const auto& qp = iqps[q];
                F_bis += qp.weight() * test_case.dirichlet_jump(qp.point()) * gb.eval_basis(qp.point()) * dn[q];
            }
        };
        if (lay.cut)
            add_interface(cl);
        for (auto dp : P.dependent_cells)
            add_interface(msh[dp]);
        f -= kappa_1 * oper_gr.transpose() * F_bis;
    }

    return f;
}

/***************************************************************************/
/*                           LOCAL CONTRIBUTIONS                            */
/***************************************************************************/

template<typename T, typename TestCase>
T subcell_kappa(const TestCase& test_case, disk::location side) {
    return (side == disk::location::IN_NEGATIVE_SIDE) ? test_case.parms.kappa_1 : test_case.parms.kappa_2;
}

/* Local matrix and rhs of a subcell of P^OK */
template<disk::CutMesh Mesh, typename TestCase>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>>
make_contrib_POK(const Mesh& msh, const disk::cut_subcell& P, const TestCase& test_case,
                 const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto kappa = subcell_kappa<T>(test_case, P.side);

    /* uncut cell without dependent cells: local operators from the cache */
    std::vector<long long> key;
    bool cacheable = prm.uncut_cache && !disk::is_cut(msh, msh[P.cell]) && P.dependent_cells.empty()
                     && uncut_cell_key(msh, msh[P.cell], kappa, key);
    if (cacheable) {
        auto it = prm.uncut_cache->find(key);
        if (it != prm.uncut_cache->end())
            return std::make_pair(it->second.first, make_rhs_jumps(msh, P, true, it->second.second, test_case, prm));
    }

    auto gr = make_gradrec_POK(msh, P, prm);                        // G
    matrix_type stab_usual = make_stabilization(msh, P, prm);           // s°
    matrix_type stab_cut   = make_interface_penalty(msh, P, prm);       // s^Gamma
    matrix_type stab_ext   = make_extension_stabilization(msh, P, prm); // s^N

    matrix_type lc = kappa * (gr.second + prm.eta_stab*stab_usual + prm.eta_N*stab_ext + stab_cut);
    auto f = make_rhs_jumps(msh, P, true, gr.first, test_case, prm);
    if (cacheable)
        (*prm.uncut_cache)[key] = std::make_pair(lc, gr.first);

    return std::make_pair(lc, f);
}

/* Local matrix and rhs of a subcell of P^KO */
template<disk::CutMesh Mesh, typename TestCase>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>>
make_contrib_PKO(const Mesh& msh, const disk::cut_subcell& P, const TestCase& test_case,
                 const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto kappa = subcell_kappa<T>(test_case, P.side);

    auto gr = make_gradrec_PKO(msh, P, prm);                     // G
    matrix_type stab_usual = make_stabilization(msh, P, prm);     // s°
    matrix_type stab_cut   = make_interface_penalty(msh, P, prm); // s^Gamma

    matrix_type lc = kappa * (gr.second + prm.eta_stab*stab_usual + stab_cut);
    auto f = make_rhs_jumps(msh, P, false, gr.first, test_case, prm);

    return std::make_pair(lc, f);
}
