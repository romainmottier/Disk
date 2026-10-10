/*
 * Local operators of the unfitted HHO method stabilized by cell agglomeration
 * ([12] Burman, Cicuttin, Delay, Ern), ported from ProtoN.
 *
 * The mesh is not modified: an agglomerated cell ("macro cell") is the union
 * of a few cells of the mesh, whose interior faces carry no unknown. Its
 * local unknowns are
 *  - uncut: [ u_Q | u_dQ ]
 *  - cut:   [ u_Q1 | u_Q2 | u_(dQ)1 | u_(dQ)2 ]
 * where the faces of dQ are the boundary faces of the union.
 */

#pragma once

#include <map>
#include <vector>

#include <Eigen/Dense>

#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/bases/cut_bases.hpp"
#include "operators.hpp"
#include "quadrature_matrices.hpp"

template<disk::CutMesh Mesh>
struct macro_cell {
    using T         = typename Mesh::coordinate_type;
    using face_type = typename Mesh::face_type;

    std::vector<size_t>                     members;  // cells of the mesh
    std::vector<std::pair<face_type, size_t>> faces;  // boundary faces and the member owning them
    bool                                    cut;
    disk::location                          loc;      // location if not cut
    T                                       h;        // diameter of the union
};

/* Agglomerated cells (groups of the agglomeration) and the other cells of the
 * mesh, ordered by their first cell. */
template<disk::CutMesh Mesh>
std::vector<macro_cell<Mesh>>
make_macro_cells(const Mesh& msh, const disk::cut_agglomeration& agglo)
{
    using T = typename Mesh::coordinate_type;

    std::vector<macro_cell<Mesh>> ret;
    std::vector<bool> done(msh.cells_size(), false);

    for (auto& cl : msh) {
        auto c = offset(msh, cl);
        if (done[c])
            continue;

        macro_cell<Mesh> Q;
        if (agglo.group_of_cell.size() == 0 || agglo.group_of_cell[c] == -1)
            Q.members = {c};
        else {
            Q.members = agglo.groups[agglo.group_of_cell[c]];
            std::sort(Q.members.begin(), Q.members.end());
        }

        std::map<size_t, int> count;
        std::vector<std::pair<typename Mesh::face_type, size_t>> all;
        std::vector<typename Mesh::point_type> pts;
        Q.cut = false;
        for (auto m : Q.members) {
            done[m] = true;
            auto mc = msh[m];
            for (auto& fc : faces(msh, mc)) {
                count[offset(msh, fc)]++;
                all.push_back( std::make_pair(fc, m) );
            }
            for (auto& p : points(msh, mc))
                pts.push_back(p);
            Q.cut = Q.cut || disk::is_cut(msh, mc);
        }
        for (auto& fm : all)
            if (count[offset(msh, fm.first)] == 1)
                Q.faces.push_back(fm);

        Q.loc = Q.cut ? disk::location::ON_INTERFACE : disk::locate(msh, msh[Q.members[0]]);

        Q.h = 0.0;
        for (size_t i = 0; i < pts.size(); i++)
            for (size_t j = i+1; j < pts.size(); j++)
                Q.h = std::max(Q.h, (pts[i]-pts[j]).to_vector().norm());

        ret.push_back(Q);
    }

    return ret;
}

/* Centroid and diameter of the side "side" of a macro cell, used to center
 * and scale its polynomial bases (whole macro cell if !subcell_centering). */
template<disk::CutMesh Mesh>
std::pair<typename Mesh::point_type, typename Mesh::coordinate_type>
macro_center_diameter(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side, bool subcell_centering)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    T area = 0.0;
    point_type bar;
    std::vector<point_type> pts;

    for (auto m : Q.members) {
        auto mc = msh[m];
        if (!subcell_centering || !disk::is_cut(msh, mc)) {
            if (subcell_centering && disk::locate(msh, mc) != side)
                continue;
            auto A = measure(msh, mc);
            area += A;
            bar = bar + barycenter(msh, mc) * A;
            for (auto& p : points(msh, mc))
                pts.push_back(p);
        }
        else {
            auto A = measure(msh, mc, side);
            area += A;
            bar = bar + barycenter(msh, mc, side) * A;
            auto tp = disk::region_points(msh, mc, side);
            pts.insert(pts.end(), tp.begin(), tp.end());
        }
    }

    T h = 0.0;
    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            h = std::max(h, (pts[i]-pts[j]).to_vector().norm());

    return std::make_pair(bar/area, h);
}

/* Frame of the orthonormal bases of the side "side" of a macro cell: the
 * union of the parts of its cells on that side (cached with keys beyond those
 * of the cells). A macro cell made of one cell has the basis of that cell. */
template<disk::CutMesh Mesh>
std::shared_ptr<const disk::cut_basis_frame<typename Mesh::coordinate_type, Mesh::dimension>>
macro_basis_frame(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
                  const mesh_parameters<Mesh>& prm)
{
    auto deg   = prm.hdi.cell_degree();
    auto cache = prm.basis.cache.get();
    if (Q.members.size() == 1)
        return disk::cut_cell_basis_frame(msh, msh[Q.members[0]], side, deg, cache);

    size_t key = 2*(msh.cells_size() + Q.members[0]) + (side == disk::location::IN_POSITIVE_SIDE ? 1 : 0);
    if (cache) {
        auto it = cache->frames.find(key);
        if (it != cache->frames.end() && it->second->degree >= deg)
            return it->second;
        deg = std::max(deg, cache->degree);
    }

    std::vector<std::pair<typename Mesh::cell_type, disk::location>> pieces;
    for (auto m : Q.members) {
        auto mc = msh[m];
        if (disk::is_cut(msh, mc) || disk::locate(msh, mc) == side)
            pieces.push_back( std::make_pair(mc, side) );
    }
    auto fr = disk::make_cut_basis_frame(msh, pieces, deg);
    if (cache)
        cache->frames[key] = fr;
    return fr;
}

template<disk::CutMesh Mesh>
class macro_layout {
public:
    size_t cbs, fbs, gbs, num_faces;
    bool   cut;
    size_t local_dofs;

    macro_layout(const macro_cell<Mesh>& Q, const disk::hho_degree_info& hdi)
    {
        cbs       = disk::cut_scalar_basis<Mesh>::size(hdi.cell_degree());
        fbs       = disk::cut_face_basis<Mesh>::size(hdi.face_degree());
        gbs       = disk::cut_vector_basis<Mesh>::size(hdi.grad_degree());
        num_faces = Q.faces.size();
        cut       = Q.cut;
        local_dofs = cut ? 2*(cbs + num_faces*fbs) : cbs + num_faces*fbs;
    }

    size_t cut_cell_offset(disk::location side) const {
        return (side == disk::location::IN_NEGATIVE_SIDE) ? 0 : cbs;
    }
    size_t cut_face_offset(disk::location side, size_t nf) const {
        return (side == disk::location::IN_NEGATIVE_SIDE) ? 2*cbs : 2*cbs + nf*fbs;
    }
    size_t cut_face_offset(disk::location side) const { return cut_face_offset(side, num_faces); }
    size_t cell_offset(disk::location side) const { return cut ? cut_cell_offset(side) : 0; }
    size_t face_offset(disk::location side) const { return cut ? cut_face_offset(side) : cbs; }
};

template<disk::CutMesh Mesh>
struct macro_bases {
    disk::cut_scalar_basis<Mesh> cb;
    disk::cut_vector_basis<Mesh> gb;

    macro_bases(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
                const mesh_parameters<Mesh>& prm)
        : cb(make(msh, Q, side, prm, prm.hdi.cell_degree())),
          gb(makev(msh, Q, side, prm, prm.hdi.grad_degree()))
    {}

private:
    static disk::cut_scalar_basis<Mesh>
    make(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
         const mesh_parameters<Mesh>& prm, size_t deg) {
        if (prm.basis.kind == disk::cut_basis_kind::orthonormal)
            return disk::cut_scalar_basis<Mesh>(macro_basis_frame(msh, Q, side, prm), deg);
        auto bh = macro_center_diameter(msh, Q, side, prm.basis.kind == disk::cut_basis_kind::subcell_monomial);
        return disk::cut_scalar_basis<Mesh>(bh.first, bh.second, deg);
    }
    static disk::cut_vector_basis<Mesh>
    makev(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
          const mesh_parameters<Mesh>& prm, size_t deg) {
        if (prm.basis.kind == disk::cut_basis_kind::orthonormal)
            return disk::cut_vector_basis<Mesh>(macro_basis_frame(msh, Q, side, prm), deg);
        auto bh = macro_center_diameter(msh, Q, side, prm.basis.kind == disk::cut_basis_kind::subcell_monomial);
        return disk::cut_vector_basis<Mesh>(bh.first, bh.second, deg);
    }
};

/* Gradient reconstruction on the side "side" of a macro cell (eq. (17)).
 * Returns (operator, G^T M G). */
template<disk::CutMesh Mesh>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>>
make_agglo_gradrec(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
                   const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto celdeg  = prm.hdi.cell_degree();
    const auto facdeg  = prm.hdi.face_degree();
    const auto graddeg = prm.hdi.grad_degree();
    macro_layout<Mesh> lay(Q, prm.hdi);
    const auto cbs = lay.cbs, fbs = lay.fbs, gbs = lay.gbs;

    static const size_t DIM = Mesh::dimension;
    macro_bases<Mesh> mb(msh, Q, side, prm);
    const auto& sgb = mb.gb.scalar_part();

    matrix_type M = matrix_type::Zero(gbs/DIM, gbs/DIM);
    matrix_type rhs_tmp = matrix_type::Zero(gbs, cbs + lay.num_faces*fbs);

    for (auto m : Q.members) {
        auto qps = integrate(msh, msh[m], celdeg - 1 + graddeg, side);
        if (!qps.empty()) {
            auto G = basis_values(sgb, qps);
            M += weighted_product(G, quadrature_weights(qps), G);
        }
        auto qpr = integrate(msh, msh[m], celdeg - 1 + facdeg, side);
        if (!qpr.empty()) {
            auto w  = quadrature_weights(qpr);
            auto G  = basis_values(sgb, qpr);
            auto dC = basis_gradients<DIM>(mb.cb, qpr);
            for (size_t d = 0; d < DIM; d++)
                add_component_rows<DIM>(rhs_tmp, 0, weighted_product(G, w, dC[d]), d);
        }
    }

    for (size_t i = 0; i < Q.faces.size(); i++) {
        const auto& fc = Q.faces[i].first;
        const Eigen::Matrix<T, Mesh::dimension, 1> n = normal(msh, msh[Q.faces[i].second], fc);
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, side, prm.basis);
        auto qps = integrate(msh, fc, facdeg + std::max(facdeg, celdeg), side);
        if (qps.empty())
            continue;
        auto w  = quadrature_weights(qps);
        auto G  = basis_values(sgb, qps);
        matrix_type GF = weighted_product(G, w, basis_values(fb, qps));
        matrix_type GC = weighted_product(G, w, basis_values(mb.cb, qps));
        for (size_t d = 0; d < DIM; d++) {
            add_component_rows<DIM>(rhs_tmp, cbs + i*fbs, GF, d, T(n(d)));
            add_component_rows<DIM>(rhs_tmp, 0, GC, d, T(-n(d)));
        }
    }

    matrix_type gr_rhs;
    if (!lay.cut)
        gr_rhs = rhs_tmp;
    else {
        gr_rhs = matrix_type::Zero(gbs, lay.local_dofs);
        place_one_side(gr_rhs, rhs_tmp, lay, side);

        macro_bases<Mesh> mb_neg(msh, Q, disk::location::IN_NEGATIVE_SIDE, prm);
        macro_bases<Mesh> mb_pos(msh, Q, disk::location::IN_POSITIVE_SIDE, prm);
        matrix_type interface_term = matrix_type::Zero(gbs, 2*cbs);
        for (auto m : Q.members) {
            auto mc = msh[m];
            if (!disk::is_cut(msh, mc))
                continue;
            auto dn   = disk::interface_normals(msh, mc, celdeg + graddeg);
            auto iqps = disk::integrate_interface(msh, mc, celdeg + graddeg, disk::location::IN_NEGATIVE_SIDE);
            auto w  = quadrature_weights(iqps);
            auto G  = basis_values(sgb, iqps);
            auto Cn = basis_values(mb_neg.cb, iqps);
            auto Cp = basis_values(mb_pos.cb, iqps);
            for (size_t d = 0; d < DIM; d++) {
                qm_vector<T> wd(w.size());
                for (Eigen::Index q = 0; q < w.size(); q++)
                    wd(q) = w(q) * dn[q](d);
                add_component_rows<DIM>(interface_term, 0, weighted_product(G, wd, Cn), d, T(-1));
                add_component_rows<DIM>(interface_term, cbs, weighted_product(G, wd, Cp), d, T(1));
            }
        }
        gr_rhs.block(0, 0, gbs, 2*cbs) += side_coeff<T>(side) * interface_term;
    }

    matrix_type oper = solve_vector_mass<DIM>(M, gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

/* s° on the side "side" of a macro cell, in its local numbering */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_agglo_stabilization(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
                         const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto celdeg = prm.hdi.cell_degree();
    auto facdeg = prm.hdi.face_degree();
    macro_layout<Mesh> lay(Q, prm.hdi);
    auto cbs = lay.cbs, fbs = lay.fbs, nf = lay.num_faces;

    macro_bases<Mesh> mb(msh, Q, side, prm);

    matrix_type stab = matrix_type::Zero(cbs + nf*fbs, cbs + nf*fbs);
    matrix_type If   = matrix_type::Identity(fbs, fbs);
    for (size_t i = 0; i < nf; i++) {
        const auto& fc = Q.faces[i].first;
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, side, prm.basis);

        matrix_type oper  = matrix_type::Zero(fbs, cbs + nf*fbs);
        matrix_type mass  = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, cbs + i*fbs, fbs, fbs) = -If;
        auto qps = integrate(msh, fc, facdeg + celdeg, side);
        if (qps.size() == 0)
            continue;
        auto w = quadrature_weights(qps);
        auto F = basis_values(fb, qps);
        mass  = weighted_product(F, w, F);
        trace = weighted_product(F, w, basis_values(mb.cb, qps));
        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        stab += oper.transpose() * mass * oper * (1./Q.h);
    }

    matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);
    auto co = lay.cell_offset(side);
    auto fo = lay.face_offset(side);
    auto nfb = nf*fbs;
    data.block(co, co, cbs, cbs) += stab.block(  0,   0, cbs, cbs);
    data.block(co, fo, cbs, nfb) += stab.block(  0, cbs, cbs, nfb);
    data.block(fo, co, nfb, cbs) += stab.block(cbs,   0, nfb, cbs);
    data.block(fo, fo, nfb, nfb) += stab.block(cbs, cbs, nfb, nfb);

    return data;
}

/* s^Gamma of a cut macro cell, assembled from its negative side */
template<disk::CutMesh Mesh>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>
make_agglo_interface_penalty(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
                             const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    macro_layout<Mesh> lay(Q, prm.hdi);
    matrix_type data = matrix_type::Zero(lay.local_dofs, lay.local_dofs);
    if (!lay.cut || side != disk::location::IN_NEGATIVE_SIDE)
        return data;

    auto celdeg = prm.hdi.cell_degree();
    auto cbs = lay.cbs;
    macro_bases<Mesh> mb_neg(msh, Q, disk::location::IN_NEGATIVE_SIDE, prm);
    macro_bases<Mesh> mb_pos(msh, Q, disk::location::IN_POSITIVE_SIDE, prm);

    for (auto m : Q.members) {
        auto mc = msh[m];
        if (!disk::is_cut(msh, mc))
            continue;
        auto iqps = disk::integrate_interface(msh, mc, 2*celdeg, disk::location::IN_NEGATIVE_SIDE);
        qm_vector<T> w = quadrature_weights(iqps) * (prm.eta_gamma / Q.h);
        auto Cn = basis_values(mb_neg.cb, iqps);
        auto Cp = basis_values(mb_pos.cb, iqps);
        matrix_type Mnp = weighted_product(Cn, w, Cp);
        data.block(  0,   0, cbs, cbs) += weighted_product(Cn, w, Cn);
        data.block(cbs, cbs, cbs, cbs) += weighted_product(Cp, w, Cp);
        data.block(  0, cbs, cbs, cbs) -= Mnp;
        data.block(cbs,   0, cbs, cbs) -= Mnp.transpose();
    }

    return data;
}

/* Right-hand side of the side "side" of a macro cell */
template<disk::CutMesh Mesh, typename TestCase>
Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>
make_agglo_rhs(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side,
               const Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>& oper_gr,
               const TestCase& test_case,
               const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    auto celdeg  = prm.hdi.cell_degree();
    auto graddeg = prm.hdi.grad_degree();
    macro_layout<Mesh> lay(Q, prm.hdi);
    auto cbs = lay.cbs;
    auto kappa_1 = test_case.parms.kappa_1;

    macro_bases<Mesh> mb(msh, Q, side, prm);
    vector_type f = vector_type::Zero(lay.local_dofs);

    // (f, w_Qi)_Qi
    for (auto m : Q.members)
        for (auto& qp : integrate(msh, msh[m], 2*celdeg, side))
            f.segment(lay.cell_offset(side), cbs) += qp.weight() * mb.cb.eval_basis(qp.point()) * test_case.rhs_fun(qp.point());

    if (!lay.cut)
        return f;

    vector_type F_bis = vector_type::Zero(lay.gbs);
    for (auto m : Q.members) {
        auto mc = msh[m];
        if (!disk::is_cut(msh, mc))
            continue;
        // kappa_1 h^-1 (g_D, [w_Q]_Gamma) and (g_N, w_Q2)
        for (auto& qp : disk::integrate_interface(msh, mc, 2*celdeg, disk::location::IN_NEGATIVE_SIDE)) {
            auto phi = mb.cb.eval_basis(qp.point());
            auto gD  = test_case.dirichlet_jump(qp.point());
            if (side == disk::location::IN_NEGATIVE_SIDE)
                f.segment(0, cbs) += kappa_1 * qp.weight() * gD * phi * prm.eta_gamma / Q.h;
            else {
                f.segment(cbs, cbs) -= kappa_1 * qp.weight() * gD * phi * prm.eta_gamma / Q.h;
                f.segment(cbs, cbs) += qp.weight() * test_case.neumann_jump(qp.point()) * phi;
            }
        }
        // lifting of the Dirichlet jump
        if (side == disk::location::IN_NEGATIVE_SIDE) {
            auto dn   = disk::interface_normals(msh, mc, 2*graddeg);
            auto iqps = disk::integrate_interface(msh, mc, 2*graddeg, disk::location::IN_NEGATIVE_SIDE);
            for (size_t q = 0; q < iqps.size(); q++) {
                const auto& qp = iqps[q];
                F_bis += qp.weight() * test_case.dirichlet_jump(qp.point()) * mb.gb.eval_basis(qp.point()) * dn[q];
            }
        }
    }
    if (side == disk::location::IN_NEGATIVE_SIDE)
        f -= kappa_1 * oper_gr.transpose() * F_bis;

    return f;
}

/* Local matrix and rhs of the side "side" of a macro cell */
template<disk::CutMesh Mesh, typename TestCase>
std::pair<Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<typename Mesh::coordinate_type, Eigen::Dynamic, 1>>
make_contrib_agglo(const Mesh& msh, const macro_cell<Mesh>& Q, disk::location side, const TestCase& test_case,
                   const mesh_parameters<Mesh>& prm)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto kappa = subcell_kappa<T>(test_case, side);

    /* agglomerated cell made of one uncut cell: local operators from the cache */
    std::vector<long long> key;
    bool cacheable = prm.uncut_cache && Q.members.size() == 1 && !Q.cut
                     && uncut_cell_key(msh, msh[Q.members[0]], kappa, key);
    if (cacheable) {
        auto it = prm.uncut_cache->find(key);
        if (it != prm.uncut_cache->end())
            return std::make_pair(it->second.first, make_agglo_rhs(msh, Q, side, it->second.second, test_case, prm));
    }

    auto gr = make_agglo_gradrec(msh, Q, side, prm);
    matrix_type stab    = make_agglo_stabilization(msh, Q, side, prm);
    matrix_type penalty = make_agglo_interface_penalty(msh, Q, side, prm);

    matrix_type lc = kappa * (gr.second + prm.eta_stab*stab + penalty);
    auto f = make_agglo_rhs(msh, Q, side, gr.first, test_case, prm);
    if (cacheable)
        (*prm.uncut_cache)[key] = std::make_pair(lc, gr.first);

    return std::make_pair(lc, f);
}

/* Sides of the macro cell carrying unknowns */
template<disk::CutMesh Mesh>
std::vector<disk::location>
macro_sides(const macro_cell<Mesh>& Q) {
    if (Q.cut)
        return {disk::location::IN_NEGATIVE_SIDE, disk::location::IN_POSITIVE_SIDE};
    return {Q.loc};
}
