/*
 * Boundary value correction with the stabilization by polynomial extension
 * (prototype -p3): the method of bvc.hpp (flat interface, Taylor operator,
 * correction gamma_T, data at p_h) in which the small cuts are stabilized by
 * the polynomial extension of
 *
 *   E. Burman, A. Ern, R. Mottier, "Unfitted HHO methods stabilized by
 *   polynomial extension for elliptic interface problems"
 *
 * instead of cell agglomeration. Each subcell (T,i) keeps its own unknowns.
 * An ill-cut subcell (S,i) (P^KO) is paired with the subcell (T,i) = N_i(S)
 * of P^OK, whose gradient reconstruction takes the faces and the interface
 * of S^i into account with the extensions q^+ of its test functions:
 *
 *   (G_{T^i} v, q)_{T^i} = (grad v_{T^i}, q)_{T^i} + (v_dT - v_{T^i}, q.n_T)_{dT^i}
 *        + sum_{S in N_i^-1(T)} (v_dS - v_{S^i}, q^+.n_S)_{dS^i},
 *
 * the interface terms being those of the boundary value correction, on T^Gamma_h
 * and on the S^Gamma_h: Ghat_1 = G_1 - L([T_m v]), Gring_1 = G_1 - L([w]),
 * (L g, q)_{T^1} = (g, q.n_h)_{T^Gamma_h} + sum_S (g, q^+.n_h)_{S^Gamma_h}.
 * The forms are
 *
 *   a_h(v,w) = sum_{(T,1) in P^OK} k1 (Ghat_1 v, Gring_1 w)_{T^1}
 *            + sum_{(T,2) in P^OK} k2 (G_2 v, G_2 w)_{T^2}
 *            + sum_{(S,i) in P^KO} k_i (grad v_{S^i}, grad w_{S^i})_{S^i}
 *            + sum_{(T,i)} k_i/h_T (Pi_F v_{T^i} - v_F, Pi_F w_{T^i} - w_F)_{dT^i}       (s°)
 *            + sum_{(T,i) in P^OK} sum_S eta_N k_i/h_T^2 (v_{S^i} - v_{T^i}^+, ...)_{T^i}  (s^N)
 *            + sum_{T cut} eta k1/h_T ([T_m v_T], [w_T])_{T^Gamma_h} + gamma_T(v, w),
 *
 * gamma_T(v,w) = (T_m(Q v).n* - (Q v).n_h, w_{T^2})_{T^Gamma_h}, Q v = k1 G_1 v - k2 G_2 v,
 * the reconstruction G_i on the interface of a cut cell being the one of the
 * subcell that carries its flux: (T,i) itself if it is in P^OK, else the
 * extension of G_{N_i(T)^i}. The right-hand side is as in bvc.hpp, the
 * lifting of g_D collecting the S^Gamma_h. Mixed order, s° and the interface
 * penalty as in -p2, s^N weighted by eta_N (option -e, as -p0).
 *
 * The cell unknowns of a subcell and of the subcells it stabilizes are
 * eliminated together (LU: the method is not symmetric). The bases are the
 * cut bases of diskpp (-b), whose frame on an ill-cut subcell is built on its
 * union with its stabilizing subcell; the Taylor operator is evaluated from
 * point values (exact interpolation along n_h), for any basis.
 */

#pragma once

#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include "bvc.hpp"

namespace bvc {

/* T_m b(x_q) for a basis of polynomials of degree p evaluated at points
 * (eval_basis_batch): the Taylor polynomial of order m at t = 0 of
 * t -> b(x_q + t nu_q), evaluated at t = rho_q, computed exactly from its
 * values at p+1 Chebyshev points of [-h, h] */
template<typename T, size_t DIM, typename Basis>
Mat<T>
taylor_values(const Basis& b, size_t p, const std::vector<disk::point<T,DIM>>& xs,
              const std::vector<vec<T,DIM>>& nus, const std::vector<T>& rhos, int m, T h)
{
    if (m <= 0 || p == 0)
        return b.eval_basis_batch(xs);
    const size_t np = p + 1;
    Vec<T> c(np);
    for (size_t j = 0; j < np; j++)
        c(j) = std::cos(M_PI * (2.0*j + 1.0) / (2.0*np));
    Mat<T> V(np, np);
    for (size_t j = 0; j < np; j++)
        for (size_t n = 0; n < np; n++)
            V(j, n) = std::pow(c(j), T(n));
    Eigen::PartialPivLU<Mat<T>> Vlu(V);

    std::vector<disk::point<T,DIM>> pts;
    pts.reserve(xs.size() * np);
    for (size_t q = 0; q < xs.size(); q++)
        for (size_t j = 0; j < np; j++)
            pts.push_back(xs[q] + disk::point<T,DIM>(vec<T,DIM>(nus[q] * (h * c(j)))));
    Mat<T> vals = b.eval_basis_batch(pts);

    Mat<T> ret = Mat<T>::Zero(xs.size(), vals.cols());
    for (size_t q = 0; q < xs.size(); q++) {
        Mat<T> a = Vlu.solve(vals.block(q*np, 0, np, vals.cols()));   // g(t) = sum a_n (t/h)^n
        T s = rhos[q] / h, sn = 1;
        for (size_t n = 0; n <= std::min<size_t>(m, p); n++, sn *= s)
            ret.row(q) += sn * a.row(n);
    }
    return ret;
}

/* A subcell (cell, side): quadrature, faces and bases */
template<typename Mesh>
struct ext_subcell {
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using basis_type = disk::cut_scalar_basis<Mesh>;

    size_t                                          cell = 0, side = 0;
    T                                               h = 0;      // diameter of the cell
    std::vector<disk::quadrature_point<T,DIM>>      vol;
    std::vector<face_piece<Mesh>>                   pieces;
    std::shared_ptr<basis_type>                     cb, gb;     // P^{k+1}, P^k

    std::vector<typename Mesh::point_type> vol_points() const { return qp_points<Mesh>(vol); }
};

/* A local contribution: dense matrix and right-hand side on blocks of cell
 * unknowns (cell, side) followed by blocks of face unknowns (face, side) */
template<typename T>
struct ext_element {
    std::vector<std::array<size_t,2>> cells, faces;
    size_t  cbs = 0, fbs = 0;
    Mat<T>  A;
    Vec<T>  b;

    size_t size() const { return cells.size()*cbs + faces.size()*fbs; }
    size_t coff(size_t c, size_t s) const {
        for (size_t i = 0; i < cells.size(); i++)
            if (cells[i][0] == c && cells[i][1] == s)
                return i*cbs;
        throw std::logic_error("bvc::ext_element: missing cell block");
    }
    size_t foff(size_t f, size_t s) const {
        for (size_t i = 0; i < faces.size(); i++)
            if (faces[i][0] == f && faces[i][1] == s)
                return cells.size()*cbs + i*fbs;
        throw std::logic_error("bvc::ext_element: missing face block");
    }
    void add_cell(size_t c, size_t s) {
        for (auto& x : cells)
            if (x[0] == c && x[1] == s)
                return;
        cells.push_back({c, s});
    }
    void add_face(size_t f, size_t s) {
        for (auto& x : faces)
            if (x[0] == f && x[1] == s)
                return;
        faces.push_back({f, s});
    }
    void allocate() {
        A = Mat<T>::Zero(size(), size());
        b = Vec<T>::Zero(size());
    }
};

/* Construction of the local contributions of a cluster of cells */
template<typename Mesh>
class ext_builder {
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using subcell = ext_subcell<Mesh>;
    using element = ext_element<T>;

    const Mesh&                                 msh;
    const problem<T,DIM>&                       prob;
    size_t                                      k;
    method_options                              mopt;
    quadrature_degrees                          qd;
    disk::cut_basis_options<T,DIM>              bopt;
    T                                           eta_N;

    std::map<std::array<size_t,2>, std::shared_ptr<subcell>>    subcells;  // of the current cluster
    std::map<size_t, std::vector<interface_node<T,DIM>>>       nodes;

public:
    size_t cbs, fbs, nsk;

    ext_builder(const Mesh& m, const problem<T,DIM>& p, size_t degree, const method_options& mo,
                const quadrature_degrees& q, const disk::cut_basis_options<T,DIM>& bo, T etaN)
        : msh(m), prob(p), k(degree), mopt(mo), qd(q), bopt(bo), eta_N(etaN),
          cbs(dim_pk(degree+1, DIM)), fbs(dim_pk(degree, DIM-1)), nsk(dim_pk(degree, DIM))
    {}

    void clear() { subcells.clear(); nodes.clear(); }

    const subcell& sub(size_t c, size_t s) {
        auto key = std::array<size_t,2>{c, s};
        auto it = subcells.find(key);
        if (it != subcells.end())
            return *it->second;
        auto sc = std::make_shared<subcell>();
        auto cl = msh[c];
        auto where = side_location(s);
        bool cut = disk::is_cut(msh, cl);
        sc->cell = c;
        sc->side = s;
        sc->h = diameter(msh, cl);
        sc->vol = cut ? integrate(msh, cl, qd.vol, where) : integrate(msh, cl, qd.vol);
        sc->cb = std::make_shared<disk::cut_scalar_basis<Mesh>>(msh, cl, k+1, where, bopt);
        sc->gb = std::make_shared<disk::cut_scalar_basis<Mesh>>(msh, cl, k, where, bopt);
        for (auto& fc : faces(msh, cl)) {
            face_piece<Mesh> fp;
            fp.fc = fc;
            fp.side = s;
            fp.n = normal(msh, cl, fc);
            auto floc = disk::locate(msh, fc);
            if (floc == disk::location::ON_INTERFACE) {
                fp.cut = true;
                fp.qps = integrate(msh, fc, qd.face, where);
            }
            else if (floc == where) {
                fp.cut = false;
                fp.qps = integrate(msh, fc, qd.face);
            }
            else
                continue;
            if (!fp.qps.empty())
                sc->pieces.push_back(fp);
        }
        subcells[key] = sc;
        return *sc;
    }

    const std::vector<interface_node<T,DIM>>& interface(size_t c) {
        auto it = nodes.find(c);
        if (it != nodes.end())
            return it->second;
        auto cl = msh[c];
        T h = diameter(msh, cl);
        std::vector<interface_node<T,DIM>> ret;
        auto iq = disk::integrate_interface(msh, cl, qd.itf, disk::location::IN_NEGATIVE_SIDE);
        auto in = disk::interface_normals(msh, cl, qd.itf);
        for (size_t q = 0; q < iq.size(); q++) {
            interface_node<T,DIM> nd;
            nd.x = iq[q].point();
            nd.w = iq[q].weight();
            for (size_t d = 0; d < DIM; d++)
                nd.nh(d) = in[q](d);
            nd.nh /= nd.nh.norm();
            nd.rho = compute_rho<T,DIM>(*prob.ls, nd.x, nd.nh, h);
            nd.ph = nd.x + disk::point<T,DIM>(vec<T,DIM>(nd.nh * nd.rho));
            vec<T,DIM> gp = prob.ls->gradient(nd.ph);
            nd.nstar = gp / gp.norm();
            ret.push_back(nd);
        }
        return nodes[c] = ret;
    }

private:
    struct itf_data {
        std::vector<disk::point<T,DIM>> xs;
        std::vector<vec<T,DIM>>         nus;
        std::vector<T>                  rhos;
        Vec<T>                          W;
        T                               h;
    };

    itf_data interface_data(size_t c) {
        itf_data d;
        const auto& nds = interface(c);
        d.W.resize(nds.size());
        for (size_t q = 0; q < nds.size(); q++) {
            d.xs.push_back(nds[q].x);
            d.nus.push_back(nds[q].nh);
            d.rhos.push_back(nds[q].rho);
            d.W(q) = nds[q].w;
        }
        d.h = diameter(msh, msh[c]);
        return d;
    }

    /* data of the interface conditions at the nodes of the cut cell c */
    void interface_loads(size_t c, Vec<T>& gD, Vec<T>& gN) {
        const auto& nds = interface(c);
        gD.resize(nds.size());
        gN.resize(nds.size());
        for (size_t q = 0; q < nds.size(); q++) {
            const auto& in = nds[q];
            if (mopt.data == data_mode::exact) {
                gD(q) = prob.gD(in.ph);
                gN(q) = prob.gN(in.ph, in.nstar);
            }
            else {
                gD(q) = prob.gD(in.x);
                gN(q) = prob.gN(in.x, in.nh);
            }
        }
    }

    /* s° of the subcell sc and (f, w) on it */
    void add_subcell_terms(element& E, const subcell& sc, bool with_matrix) {
        auto pts = sc.vol_points();
        Vec<T> wf(pts.size());
        for (size_t q = 0; q < pts.size(); q++)
            wf(q) = sc.vol[q].weight() * prob.f(pts[q], sc.side);
        size_t co = E.coff(sc.cell, sc.side);
        E.b.segment(co, cbs) += sc.cb->eval_basis_batch(pts).transpose() * wf;
        if (!with_matrix)
            return;
        T kap = prob.kappa[sc.side];
        for (auto& fp : sc.pieces) {
            auto fpts = qp_points<Mesh>(fp.qps);
            Vec<T> w = qp_weights<T>(fp.qps);
            auto fb = piece_basis(msh, fp, k, bopt);
            Mat<T> PF = fb.eval_basis_batch(fpts);
            Mat<T> PT = sc.cb->eval_basis_batch(fpts);
            Mat<T> MF = PF.transpose() * w.asDiagonal() * PF;
            Mat<T> R = Mat<T>::Zero(fbs, E.size());
            R.block(0, co, fbs, cbs) = MF.ldlt().solve(PF.transpose() * w.asDiagonal() * PT);
            R.block(0, E.foff(offset(msh, fp.fc), sc.side), fbs, fbs) -= Mat<T>::Identity(fbs, fbs);
            E.A += (mopt.stab * kap / sc.h) * R.transpose() * MF * R;
        }
    }

    /* jumps [v] (columns of E) at the interface nodes of the cut cell c,
     * plain and with the Taylor operator */
    void jumps(const element& E, size_t c, const itf_data& id, Mat<T>& Jp, Mat<T>& Jt) {
        const auto& s0 = sub(c, 0);
        const auto& s1 = sub(c, 1);
        size_t nq = id.xs.size();
        Jp = Mat<T>::Zero(nq, E.size());
        Jt = Mat<T>::Zero(nq, E.size());
        Mat<T> P0 = s0.cb->eval_basis_batch(id.xs), P1 = s1.cb->eval_basis_batch(id.xs);
        Jp.block(0, E.coff(c, 0), nq, cbs) = P0;
        Jp.block(0, E.coff(c, 1), nq, cbs) = -P1;
        if (mopt.m > 0) {
            Jt.block(0, E.coff(c, 0), nq, cbs) = taylor_values<T,DIM>(*s0.cb, k+1, id.xs, id.nus, id.rhos, mopt.m, id.h);
            Jt.block(0, E.coff(c, 1), nq, cbs) = -taylor_values<T,DIM>(*s1.cb, k+1, id.xs, id.nus, id.rhos, mopt.m, id.h);
        }
        else
            Jt = Jp;
    }

public:
    /* element of a subcell (T,i) of P^OK and of the ill-cut subcells it stabilizes */
    element ok_element(const disk::cut_subcell& P, bool with_matrix = true) {
        const size_t c = P.cell, i = (P.side == disk::location::IN_NEGATIVE_SIDE) ? 0 : 1;
        const bool tcut = disk::is_cut(msh, msh[c]);
        const auto& st = sub(c, i);
        element E;
        E.cbs = cbs;
        E.fbs = fbs;
        E.add_cell(c, i);
        if (i == 0 && tcut)
            E.add_cell(c, 1);
        for (auto s : P.dependent_cells) {
            E.add_cell(s, i);
            if (i == 0)
                E.add_cell(s, 1);
        }
        for (auto& fp : st.pieces)
            E.add_face(offset(msh, fp.fc), i);
        for (auto s : P.dependent_cells)
            for (auto& fp : sub(s, i).pieces)
                E.add_face(offset(msh, fp.fc), i);
        E.allocate();
        add_subcell_terms(E, st, with_matrix);
        if (!with_matrix)
            return E;

        const size_t n = E.size();
        const T kap = prob.kappa[i];
        const auto& gb = *st.gb;

        /* gradient reconstruction without the interface terms: R0 */
        Mat<T> R0 = Mat<T>::Zero(DIM*nsk, n);
        {
            auto pts = st.vol_points();
            Vec<T> w = qp_weights<T>(st.vol);
            Mat<T> Pk = gb.eval_basis_batch(pts);
            auto dC = st.cb->eval_gradients_batch(pts);
            for (size_t d = 0; d < DIM; d++)
                R0.block(d*nsk, E.coff(c, i), nsk, cbs) += Pk.transpose() * w.asDiagonal() * dC[d];
        }
        auto add_faces = [&](const subcell& sc) {
            for (auto& fp : sc.pieces) {
                auto fpts = qp_points<Mesh>(fp.qps);
                Vec<T> w = qp_weights<T>(fp.qps);
                auto fb = piece_basis(msh, fp, k, bopt);
                Mat<T> PF = fb.eval_basis_batch(fpts);
                Mat<T> PC = sc.cb->eval_basis_batch(fpts);
                Mat<T> Pk = gb.eval_basis_batch(fpts);
                size_t fo = E.foff(offset(msh, fp.fc), i), co = E.coff(sc.cell, i);
                for (size_t d = 0; d < DIM; d++) {
                    Vec<T> wn = w * fp.n(d);
                    R0.block(d*nsk, fo, nsk, fbs) += Pk.transpose() * wn.asDiagonal() * PF;
                    R0.block(d*nsk, co, nsk, cbs) -= Pk.transpose() * wn.asDiagonal() * PC;
                }
            }
        };
        add_faces(st);
        for (auto s : P.dependent_cells)
            add_faces(sub(s, i));

        Mat<T> M;
        {
            auto pts = st.vol_points();
            Vec<T> w = qp_weights<T>(st.vol);
            Mat<T> Pk = gb.eval_basis_batch(pts);
            M = Pk.transpose() * w.asDiagonal() * Pk;
        }
        Eigen::LDLT<Mat<T>> Mldlt(M);
        auto Minv = [&](const Mat<T>& R) {
            Mat<T> ret(R.rows(), R.cols());
            for (size_t d = 0; d < DIM; d++)
                ret.block(d*nsk, 0, nsk, R.cols()) = Mldlt.solve(R.block(d*nsk, 0, nsk, R.cols()));
            return ret;
        };

        /* the cut cells whose interface is seen from (T,i) */
        std::vector<size_t> owned;
        if (tcut)
            owned.push_back(c);
        for (auto s : P.dependent_cells)
            owned.push_back(s);

        if (i == 0) {
            /* Ghat_1, Gring_1 and the lifting of g_D */
            Mat<T> Rhat = R0, Rring = R0;
            Vec<T> lvec = Vec<T>::Zero(DIM*nsk);
            for (auto cc : owned) {
                auto id = interface_data(cc);
                Mat<T> Jp, Jt;
                jumps(E, cc, id, Jp, Jt);
                Mat<T> Pk = gb.eval_basis_batch(id.xs);
                Vec<T> gD, gN;
                interface_loads(cc, gD, gN);
                const auto& nds = interface(cc);
                for (size_t d = 0; d < DIM; d++) {
                    Vec<T> wn(id.W.size());
                    for (size_t q = 0; q < nds.size(); q++)
                        wn(q) = id.W(q) * nds[q].nh(d);
                    Rhat.block(d*nsk, 0, nsk, n) -= Pk.transpose() * wn.asDiagonal() * Jt;
                    Rring.block(d*nsk, 0, nsk, n) -= Pk.transpose() * wn.asDiagonal() * Jp;
                    lvec.segment(d*nsk, nsk) += Pk.transpose() * wn.cwiseProduct(gD);
                }
            }
            E.A += kap * Rring.transpose() * Minv(Rhat);
            E.b -= kap * Rring.transpose() * Minv(lvec);
        }
        else
            E.A += kap * R0.transpose() * Minv(R0);

        /* gamma_T on the interface of the owned cut cells, flux of side i */
        if (mopt.gamma) {
            Mat<T> G0 = Minv(R0);
            T sgn = (i == 0) ? prob.kappa[0] : -prob.kappa[1];
            for (auto cc : owned) {
                auto id = interface_data(cc);
                const auto& nds = interface(cc);
                size_t nq = id.xs.size();
                Mat<T> Pk = gb.eval_basis_batch(id.xs);
                Mat<T> PkT = (mopt.m > 0) ? taylor_values<T,DIM>(gb, k, id.xs, id.nus, id.rhos, mopt.m, id.h) : Pk;
                Mat<T> Gam = Mat<T>::Zero(nq, n);
                for (size_t d = 0; d < DIM; d++) {
                    Mat<T> Gd = Pk * G0.block(d*nsk, 0, nsk, n);
                    Mat<T> GTd = PkT * G0.block(d*nsk, 0, nsk, n);
                    for (size_t q = 0; q < nq; q++)
                        Gam.row(q) += nds[q].nstar(d) * GTd.row(q) - nds[q].nh(d) * Gd.row(q);
                }
                Mat<T> P1 = sub(cc, 1).cb->eval_basis_batch(id.xs);
                E.A.block(E.coff(cc, 1), 0, cbs, n) += sgn * P1.transpose() * id.W.asDiagonal() * Gam;
            }
        }

        /* s^N: eta_N k_i/h_T^2 (v_{S^i} - v_{T^i}^+, w_{S^i} - w_{T^i}^+)_{T^i} */
        if (!P.dependent_cells.empty()) {
            auto pts = st.vol_points();
            Vec<T> w = qp_weights<T>(st.vol) * (eta_N * kap / (st.h * st.h));
            Mat<T> C = st.cb->eval_basis_batch(pts);
            for (auto s : P.dependent_cells) {
                Mat<T> D = sub(s, i).cb->eval_basis_batch(pts);
                size_t oc = E.coff(c, i), od = E.coff(s, i);
                Mat<T> MCD = C.transpose() * w.asDiagonal() * D;
                E.A.block(oc, oc, cbs, cbs) += C.transpose() * w.asDiagonal() * C;
                E.A.block(oc, od, cbs, cbs) -= MCD;
                E.A.block(od, oc, cbs, cbs) -= MCD.transpose();
                E.A.block(od, od, cbs, cbs) += D.transpose() * w.asDiagonal() * D;
            }
        }
        return E;
    }

    /* element of an ill-cut subcell (S,i): k_i (grad v, grad w)_{S^i}, s°, f */
    element ko_element(const disk::cut_subcell& P) {
        const size_t c = P.cell, i = (P.side == disk::location::IN_NEGATIVE_SIDE) ? 0 : 1;
        const auto& sc = sub(c, i);
        element E;
        E.cbs = cbs;
        E.fbs = fbs;
        E.add_cell(c, i);
        for (auto& fp : sc.pieces)
            E.add_face(offset(msh, fp.fc), i);
        E.allocate();
        add_subcell_terms(E, sc, true);
        auto pts = sc.vol_points();
        Vec<T> w = qp_weights<T>(sc.vol);
        auto dC = sc.cb->eval_gradients_batch(pts);
        for (size_t d = 0; d < DIM; d++)
            E.A.block(0, 0, cbs, cbs) += prob.kappa[i] * dC[d].transpose() * w.asDiagonal() * dC[d];
        return E;
    }

    /* element of a cut cell: interface penalty and interface data */
    element interface_element(size_t c) {
        element E;
        E.cbs = cbs;
        E.fbs = fbs;
        E.add_cell(c, 0);
        E.add_cell(c, 1);
        E.allocate();
        auto id = interface_data(c);
        Mat<T> Jp, Jt;
        jumps(E, c, id, Jp, Jt);
        Vec<T> gD, gN;
        interface_loads(c, gD, gN);
        T pen = mopt.eta * prob.kappa[0] / id.h;
        E.A += pen * Jp.transpose() * id.W.asDiagonal() * Jt;
        E.b += pen * Jp.transpose() * id.W.cwiseProduct(gD);
        Mat<T> P1 = sub(c, 1).cb->eval_basis_batch(id.xs);
        E.b.segment(E.coff(c, 1), cbs) += P1.transpose() * id.W.cwiseProduct(gN);
        return E;
    }

    /* coefficients in the basis of (c, s) of the polynomial p of the basis of
     * (t, s), by least squares on the points of the subcell (t, s) */
    Vec<T> convert(size_t c, size_t t, size_t s, const Vec<T>& p) {
        const auto& st = sub(t, s);
        auto pts = st.vol_points();
        Vec<T> w = qp_weights<T>(st.vol);
        Mat<T> B = sub(c, s).cb->eval_basis_batch(pts);
        Mat<T> BT = st.cb->eval_basis_batch(pts);
        Mat<T> M = B.transpose() * w.asDiagonal() * B;
        return M.ldlt().solve(B.transpose() * w.asDiagonal() * (BT * p));
    }

    /* L2 projection of u_s on the subcell (c, s) and on the subcells (d, s)
     * of the cells d it stabilizes, in the basis of (c, s) */
    Vec<T> project_cell(size_t c, size_t s, const std::vector<size_t>& dependents = {}) {
        const auto& sc = sub(c, s);
        Mat<T> M = Mat<T>::Zero(cbs, cbs);
        Vec<T> r = Vec<T>::Zero(cbs);
        auto add = [&](const subcell& sd) {
            auto pts = sd.vol_points();
            Vec<T> w = qp_weights<T>(sd.vol), wu(pts.size());
            for (size_t q = 0; q < pts.size(); q++)
                wu(q) = w(q) * prob.u(pts[q], s);
            Mat<T> B = sc.cb->eval_basis_batch(pts);
            M += B.transpose() * w.asDiagonal() * B;
            r += B.transpose() * wu;
        };
        add(sc);
        for (auto d : dependents)
            add(sub(d, s));
        return M.ldlt().solve(r);
    }

    /* L2 projection of u_s on the face piece (f, s) of the cell c */
    Vec<T> project_face(size_t c, size_t f, size_t s) {
        for (auto& fp : sub(c, s).pieces) {
            if (offset(msh, fp.fc) != f)
                continue;
            auto pts = qp_points<Mesh>(fp.qps);
            Vec<T> w = qp_weights<T>(fp.qps), wu(pts.size());
            for (size_t q = 0; q < pts.size(); q++)
                wu(q) = w(q) * prob.u(pts[q], s);
            auto fb = piece_basis(msh, fp, k, bopt);
            Mat<T> PF = fb.eval_basis_batch(pts);
            Mat<T> MF = PF.transpose() * w.asDiagonal() * PF;
            return MF.ldlt().solve(PF.transpose() * wu);
        }
        throw std::logic_error("bvc::project_face: missing face");
    }

    /* squared errors H1 and L2 of the polynomial p of the basis of (b, s) on the subcell (c, s) */
    std::pair<T,T> errors(size_t c, size_t b, size_t s, const Vec<T>& p) {
        const auto& sc = sub(c, s);
        const auto& sb = sub(b, s);
        auto pts = sc.vol_points();
        Vec<T> uh = sb.cb->eval_basis_batch(pts) * p;
        auto dP = sb.cb->eval_gradients_batch(pts);
        T e1 = 0, e0 = 0;
        for (size_t q = 0; q < pts.size(); q++) {
            T w = sc.vol[q].weight();
            T du = uh(q) - prob.u(pts[q], s);
            auto ge = prob.grad_u(pts[q], s);
            for (size_t d = 0; d < DIM; d++) {
                T gh = dP[d].row(q).dot(p);
                e1 += w * (gh - ge(d)) * (gh - ge(d));
            }
            e0 += w * du * du;
        }
        return {e1, e0};
    }
};

/* One level of -p3: mesh with the pairing of the polynomial extension, local
 * contributions and static condensation by clusters of cells, solve, errors */
template<typename Mesh, typename Timer>
run_result<typename Mesh::coordinate_type>
run_level_ext(const simulation_options& opts, std::shared_ptr<level_set<typename Mesh::coordinate_type,
              Mesh::dimension>> ls, size_t k, size_t l, const method_options& mopt)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using pre = preprocessor<Mesh>;
    using cut_mesh_type = typename pre::cut_mesh_type;
    using element = ext_element<T>;

    run_result<T> res;
    Timer tc;

    tc.tic();
    auto msh = pre::msh_generation(ls, l, opts, geometry_options(opts, k), true);
    auto partition = disk::make_subcell_partition(msh);
    if (opts.debug)
        write_geometry<pre>(msh, l);
    const auto& P_OK = partition.first;
    const auto& P_KO = partition.second;
    res.t_mesh = tc.toc();

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
    bopt.kind = (opts.basis == 0) ? disk::cut_basis_kind::cell_monomial
              : (opts.basis == 1) ? disk::cut_basis_kind::subcell_monomial
                                  : disk::cut_basis_kind::orthonormal;
    bopt.cache = std::make_shared<disk::cut_basis_cache<T,DIM>>();
    bopt.cache->degree = k+1;
    bopt.cache->face_degree = k;
    T eta_N = (opts.eta_N >= 0) ? T(opts.eta_N) : T(-opts.eta_N * (k+1) * (k+1));

    const size_t nc = msh.cells_size(), nf = msh.faces_size();
    res.cells = nc;
    for (auto& cl : msh) {
        res.h = std::max(res.h, diameter(msh, cl));
        if (disk::is_cut(msh, cl))
            res.cut++;
    }

    /* clusters: a cell and the ill-cut cells it stabilizes (union-find) */
    std::vector<size_t> parent(nc);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<size_t(size_t)> find = [&](size_t x) { return parent[x] == x ? x : parent[x] = find(parent[x]); };
    for (auto& P : P_OK)
        for (auto s : P.dependent_cells) {
            auto a = find(P.cell), b = find(s);
            if (a != b)
                parent[std::max(a, b)] = std::min(a, b);
        }
    std::map<size_t, std::vector<size_t>> cluster_map;
    for (size_t c = 0; c < nc; c++)
        cluster_map[find(c)].push_back(c);
    std::vector<std::vector<size_t>> clusters;
    for (auto& cm : cluster_map) {
        clusters.push_back(cm.second);
        if (cm.second.size() > 1) {
            res.agglomerated++;   // clusters of several cells
            res.largest = std::max(res.largest, cm.second.size());
        }
    }
    res.macros = clusters.size();

    /* subcells of P^OK and P^KO by (cell, side) */
    std::vector<std::array<long,2>> ok_of(nc, {-1, -1}), ko_of(nc, {-1, -1});
    for (size_t j = 0; j < P_OK.size(); j++)
        ok_of[P_OK[j].cell][P_OK[j].side == disk::location::IN_NEGATIVE_SIDE ? 0 : 1] = j;
    for (size_t j = 0; j < P_KO.size(); j++)
        ko_of[P_KO[j].cell][P_KO[j].side == disk::location::IN_NEGATIVE_SIDE ? 0 : 1] = j;

    /* face unknowns: (face, side) for the faces of the subcells, Dirichlet
     * values on the boundary */
    ext_builder<cut_mesh_type> B(msh, prob, k, mopt, qd, bopt, eta_N);
    const size_t cbs = B.cbs, fbs = B.fbs;
    std::vector<std::array<bool,2>> used(nf, {false, false});
    for (auto& cl : msh) {
        bool cut = disk::is_cut(msh, cl);
        for (size_t s = 0; s < 2; s++) {
            if (!cut && disk::locate(msh, cl) != side_location(s))
                continue;
            for (auto& fc : faces(msh, cl)) {
                auto floc = disk::locate(msh, fc);
                if (floc == disk::location::ON_INTERFACE || floc == side_location(s))
                    used[offset(msh, fc)][s] = true;
            }
        }
    }
    std::vector<std::array<long,2>> fdof(nf, {-1, -1});
    std::vector<bool> is_dir(nf, false);
    size_t ndof = 0;
    {
        size_t f = 0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++, f++) {
            if (msh.is_boundary(*itor)) {
                is_dir[f] = true;
                continue;
            }
            for (size_t s = 0; s < 2; s++)
                if (used[f][s]) {
                    fdof[f][s] = ndof;
                    ndof += fbs;
                }
        }
    }
    res.ndof = ndof;

    /* the local contributions of a cluster */
    std::map<std::vector<long long>, Mat<T>> cache;
    auto cluster_elements = [&](const std::vector<size_t>& cells) {
        std::vector<element> E;
        for (auto c : cells) {
            for (size_t s = 0; s < 2; s++) {
                if (ok_of[c][s] >= 0) {
                    const auto& P = P_OK[ok_of[c][s]];
                    std::vector<long long> key;
                    bool cacheable = opts.cache_uncut && !disk::is_cut(msh, msh[c]) && P.dependent_cells.empty()
                                     && uncut_cell_key(msh, msh[c], prob.kappa[s], key);
                    if (cacheable && cache.count(key)) {
                        auto e = B.ok_element(P, false);
                        e.A = cache.at(key);
                        E.push_back(std::move(e));
                        continue;
                    }
                    E.push_back(B.ok_element(P));
                    if (cacheable)
                        cache[key] = E.back().A;
                }
                if (ko_of[c][s] >= 0)
                    E.push_back(B.ko_element(P_KO[ko_of[c][s]]));
            }
            if (disk::is_cut(msh, msh[c]))
                E.push_back(B.interface_element(c));
        }
        return E;
    };

    /* condensation of a cluster: layout of its cell and face unknowns */
    struct cluster_system {
        std::vector<std::array<size_t,2>> cells, faces;
        Mat<T> K;
        Vec<T> r;
        std::vector<std::vector<size_t>> maps;   // element -> cluster positions
    };
    auto build_system = [&](const std::vector<element>& E) {
        cluster_system cs;
        auto add = [](std::vector<std::array<size_t,2>>& v, const std::array<size_t,2>& x) {
            for (size_t i = 0; i < v.size(); i++)
                if (v[i] == x)
                    return i;
            v.push_back(x);
            return v.size() - 1;
        };
        for (auto& e : E) {
            for (auto& x : e.cells) add(cs.cells, x);
            for (auto& x : e.faces) add(cs.faces, x);
        }
        size_t ncd = cs.cells.size()*cbs, n = ncd + cs.faces.size()*fbs;
        cs.K = Mat<T>::Zero(n, n);
        cs.r = Vec<T>::Zero(n);
        for (auto& e : E) {
            std::vector<size_t> map(e.size());
            for (size_t i = 0; i < e.cells.size(); i++) {
                size_t p = add(cs.cells, e.cells[i]);
                for (size_t j = 0; j < cbs; j++)
                    map[i*cbs + j] = p*cbs + j;
            }
            for (size_t i = 0; i < e.faces.size(); i++) {
                size_t p = add(cs.faces, e.faces[i]);
                for (size_t j = 0; j < fbs; j++)
                    map[e.cells.size()*cbs + i*fbs + j] = ncd + p*fbs + j;
            }
            for (size_t a = 0; a < e.size(); a++) {
                cs.r(map[a]) += e.b(a);
                for (size_t b = 0; b < e.size(); b++)
                    cs.K(map[a], map[b]) += e.A(a, b);
            }
            cs.maps.push_back(map);
        }
        return cs;
    };
    /* values of the face unknowns of a cluster: from the global solution, or
     * Dirichlet values */
    auto face_values = [&](const cluster_system& cs, const Vec<T>* sol, std::vector<long>& gidx, Vec<T>& vals) {
        size_t nfd = cs.faces.size()*fbs;
        gidx.assign(nfd, -1);
        vals = Vec<T>::Zero(nfd);
        for (size_t i = 0; i < cs.faces.size(); i++) {
            size_t f = cs.faces[i][0], s = cs.faces[i][1];
            if (is_dir[f]) {
                /* a cell of the cluster owning this face piece */
                size_t owner = 0;
                bool found = false;
                for (auto& c : cs.cells) {
                    if (c[1] != s)
                        continue;
                    for (auto& fc : faces(msh, msh[c[0]]))
                        if (offset(msh, fc) == f) {
                            owner = c[0];
                            found = true;
                        }
                    if (found)
                        break;
                }
                vals.segment(i*fbs, fbs) = B.project_face(owner, f, s);
            }
            else {
                for (size_t j = 0; j < fbs; j++)
                    gidx[i*fbs + j] = fdof[f][s] + j;
                if (sol)
                    vals.segment(i*fbs, fbs) = sol->segment(fdof[f][s], fbs);
            }
        }
    };

    /* assembly of the condensed system */
    tc.tic();
    std::vector<Eigen::Triplet<T>> trip;
    Vec<T> RHS = Vec<T>::Zero(ndof);
    for (auto& cells : clusters) {
        B.clear();
        auto E = cluster_elements(cells);
        auto cs = build_system(E);
        size_t ncd = cs.cells.size()*cbs, nfd = cs.faces.size()*fbs;
        std::vector<long> gidx;
        Vec<T> dval;
        face_values(cs, nullptr, gidx, dval);
        Eigen::PartialPivLU<Mat<T>> lu(cs.K.topLeftCorner(ncd, ncd));
        Mat<T> X = lu.solve(cs.K.topRightCorner(ncd, nfd));
        Vec<T> y = lu.solve(cs.r.head(ncd));
        Mat<T> S = cs.K.bottomRightCorner(nfd, nfd) - cs.K.bottomLeftCorner(nfd, ncd) * X;
        Vec<T> g = cs.r.tail(nfd) - cs.K.bottomLeftCorner(nfd, ncd) * y;
        for (size_t a = 0; a < nfd; a++) {
            if (gidx[a] < 0)
                continue;
            T ra = g(a);
            for (size_t b = 0; b < nfd; b++) {
                if (gidx[b] >= 0)
                    trip.emplace_back(gidx[a], gidx[b], S(a, b));
                else
                    ra -= S(a, b) * dval(b);
            }
            RHS(gidx[a]) += ra;
        }
        for (auto& e : E)
            if (e.cells.size() > 1 || !e.faces.empty()) {
                T nA = e.A.cwiseAbs().maxCoeff();
                if (nA > 0 && e.cells.size() > 1)
                    res.nonsym = std::max(res.nonsym, (e.A - e.A.transpose()).cwiseAbs().maxCoeff() / nA);
            }
        for (auto c : cells)
            if (disk::is_cut(msh, msh[c]))
                for (auto& nd : B.interface(c))
                    res.rho_max = std::max(res.rho_max, std::abs(nd.rho));
    }
    Eigen::SparseMatrix<T> LHS(ndof, ndof);
    LHS.setFromTriplets(trip.begin(), trip.end());
    trip.clear();
    trip.shrink_to_fit();
    res.t_asm = tc.toc();

    tc.tic();
    Vec<T> sol = solve_unsymmetric(LHS, RHS, opts.bvc_solver);
    res.t_solve = tc.toc();

    /* recovery of the cell unknowns, errors with the own unknowns of the
     * subcells (and with the extension of the paired polynomial on the
     * ill-cut subcells), energy error a_h(e, e), e = uhat_h - Ihat_h u */
    tc.tic();
    T e0 = 0, e1 = 0, x0 = 0, x1 = 0, en = 0;
    for (auto& cells : clusters) {
        B.clear();
        auto E = cluster_elements(cells);
        auto cs = build_system(E);
        size_t ncd = cs.cells.size()*cbs, nfd = cs.faces.size()*fbs;
        std::vector<long> gidx;
        Vec<T> fval;
        face_values(cs, &sol, gidx, fval);
        Eigen::PartialPivLU<Mat<T>> lu(cs.K.topLeftCorner(ncd, ncd));
        Vec<T> cval = lu.solve(cs.r.head(ncd) - cs.K.topRightCorner(ncd, nfd) * fval);

        /* interpolant: L2 projection on a subcell of P^OK and on the ill-cut subcells it stabilizes, extended to them */
        Vec<T> I(ncd + nfd);
        std::map<std::array<size_t,2>, Vec<T>> proj;
        for (size_t i = 0; i < cs.cells.size(); i++) {
            size_t c = cs.cells[i][0], s = cs.cells[i][1];
            if (ok_of[c][s] >= 0)
                proj[cs.cells[i]] = B.project_cell(c, s, P_OK[ok_of[c][s]].dependent_cells);
        }
        for (size_t i = 0; i < cs.cells.size(); i++) {
            size_t c = cs.cells[i][0], s = cs.cells[i][1];
            if (ok_of[c][s] >= 0)
                I.segment(i*cbs, cbs) = proj.at(cs.cells[i]);
            else {
                size_t t = msh.cut_cell_data[c].stabilizing_cell.value();
                I.segment(i*cbs, cbs) = B.convert(c, t, s, proj.at({t, s}));
            }
        }
        for (size_t i = 0; i < cs.faces.size(); i++) {
            size_t f = cs.faces[i][0], s = cs.faces[i][1];
            for (auto& c : cs.cells) {
                if (c[1] != s)
                    continue;
                bool has = false;
                for (auto& fp : B.sub(c[0], s).pieces)
                    has = has || (offset(msh, fp.fc) == f);
                if (has) {
                    I.segment(ncd + i*fbs, fbs) = B.project_face(c[0], f, s);
                    break;
                }
            }
        }
        Vec<T> x(ncd + nfd);
        x.head(ncd) = cval;
        x.tail(nfd) = fval;
        Vec<T> err = x - I;
        for (size_t j = 0; j < E.size(); j++) {
            Vec<T> el(E[j].size());
            for (size_t a = 0; a < E[j].size(); a++)
                el(a) = err(cs.maps[j][a]);
            en += el.dot(E[j].A * el);
        }

        /* errors */
        for (size_t i = 0; i < cs.cells.size(); i++) {
            size_t c = cs.cells[i][0], s = cs.cells[i][1];
            Vec<T> v = cval.segment(i*cbs, cbs);
            auto ow = B.errors(c, c, s, v);
            e1 += ow.first;
            e0 += ow.second;
            if (ok_of[c][s] >= 0) {
                x1 += ow.first;
                x0 += ow.second;
            }
            else {
                size_t t = msh.cut_cell_data[c].stabilizing_cell.value();
                size_t it = 0;
                while (cs.cells[it] != std::array<size_t,2>{t, s})
                    it++;
                auto ex = B.errors(c, t, s, Vec<T>(cval.segment(it*cbs, cbs)));
                x1 += ex.first;
                x0 += ex.second;
            }
        }
    }
    res.t_err = tc.toc();
    res.H1 = std::sqrt(e1);
    res.L2 = std::sqrt(e0);
    res.H1_ext = std::sqrt(x1);
    res.L2_ext = std::sqrt(x0);
    res.en2 = en;
    res.EN = std::sqrt(std::max(en, T(0)));
    return res;
}

} // namespace bvc
