//
//  dirk_coupling_hho_scheme_condensed_local_v2.hpp
//
//  Fourth, SEPARATE prototype: fixes dirk_coupling_hho_scheme_condensed_
//  local.hpp's remaining inefficiency (that file is untouched -- new,
//  additive twin, same public interface). There, M1..M4 were computed as
//  FULL n_c/n_f-sized sparse products (Mcinv*(Kcc*Pc) etc.) and only
//  reduced to the small active submatrix AFTERWARD via extract() -- so
//  compose() still paid for multiplying against Sffinv/Mcinv at full
//  mesh size even though the final result only needed a tiny block.
//
//  Fix: extract the SMALL submatrices of Kcc, Kcf, Kfc, Mcinv, Sffinv
//  FIRST (restricted to active_c/active_f from the very start), THEN do
//  the small matrix products. A full-mesh-sized matrix or matrix product
//  is never touched.
//
//  Why this is EXACTLY equivalent (not an approximation) -- the same
//  block-locality argument used in _condensed_local.hpp, one step
//  further:
//    - Mcinv, Kcc are block-diagonal PER CELL: for i in active_c, M1(i,j)
//      = sum_k Mcinv(i,k)*Kcc(k,j)*Pc(j,j) has k forced into i's own
//      cell (Mcinv locality) and j forced into k's cell too (Kcc
//      locality) -- i.e. j is automatically in the SAME active cell,
//      hence Pc(j,j)=1 always for surviving terms. So M1_active =
//      Mcinv_active*Kcc_active EXACTLY (the Pc factor is redundant once
//      restricted to active_c on both sides -- never needs multiplying).
//    - Similarly for M2=Mcinv*(Kcf*Pf): a fine cell's faces are ALL
//      marked Pf=1 by assemble_P's own rule ("fine if any neighbour cell
//      is fine"), EXCEPT Dirichlet boundary faces (assemble_P explicitly
//      skips marking those) -- but such a face is then simply NOT in
//      active_f (built from Pf's support) and thus never appears as a
//      column of Kcf_active either, reproducing the same masking as
//      multiplying by Pf, without ever forming Pf explicitly.
//    - Sffinv is block-diagonal PER FACE (per interface-pair in the
//      general elastoacoustic case, but plain per-face here since this
//      problem is acoustic-only) -- so restricting it to (active_f,
//      active_f) exactly captures the needed diagonal blocks, and M3/M4
//      built from Sffinv_active*(Kfc_active*M1_active /
//      Kfc_active*M2_active) are exact by the same row/column-locality
//      argument as above.
//
//  So this file computes IDENTICAL numbers to _condensed_local.hpp (up
//  to solver round-off), just without ever forming or multiplying a
//  full-mesh-sized matrix in compose().
//

#ifndef dirk_coupling_hho_scheme_condensed_local_v2_hpp
#define dirk_coupling_hho_scheme_condensed_local_v2_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <functional>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "dirk_coupling_hho_scheme.hpp" // erk_coupling_hho_scheme + dirk_weight_LTS_fine (any DirkOp)

template<typename T>
class dirk_coupling_hho_scheme_condensed_local_v2 {

    private:

    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_Acc;
    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_S;
    Eigen::SparseMatrix<T> m_Acf_active, m_Afc_active;
    std::vector<int> m_active_c;
    std::vector<int> m_active_f;
    size_t m_n_c_dof;
    size_t m_n_f_dof;
    bool m_composed_Q;

    static Eigen::SparseMatrix<T> extract(const Eigen::SparseMatrix<T> & M,
                                           const std::vector<int> & rows,
                                           const std::vector<int> & cols) {
        std::unordered_map<int,int> row_map; row_map.reserve(rows.size());
        for (int i = 0; i < (int)rows.size(); ++i) row_map[rows[i]] = i;
        std::unordered_map<int,int> col_map; col_map.reserve(cols.size());
        for (int j = 0; j < (int)cols.size(); ++j) col_map[cols[j]] = j;

        std::vector<Eigen::Triplet<T>> trips;
        for (int k = 0; k < M.outerSize(); ++k)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, k); it; ++it) {
                auto ri = row_map.find((int)it.row());
                if (ri == row_map.end()) continue;
                auto ci = col_map.find((int)it.col());
                if (ci == col_map.end()) continue;
                trips.emplace_back(ri->second, ci->second, it.value());
            }
        Eigen::SparseMatrix<T> out((int)rows.size(), (int)cols.size());
        out.setFromTriplets(trips.begin(), trips.end());
        return out;
    }

    public:

    dirk_coupling_hho_scheme_condensed_local_v2() : m_n_c_dof(0), m_n_f_dof(0), m_composed_Q(false) {}

    void compose(erk_coupling_hho_scheme<T> & erk_an,
                 const Eigen::SparseMatrix<T> & Pfine,
                 T gamma,
                 T dt) {

        m_n_c_dof = erk_an.n_c_dof();
        m_n_f_dof = erk_an.n_f_dof();
        size_t nc = m_n_c_dof;
        size_t nf = m_n_f_dof;

        Eigen::SparseMatrix<T> Pc = Pfine.block(0, 0, nc, nc);
        Eigen::SparseMatrix<T> Pf = Pfine.block(nc, nc, nf, nf);

        m_active_c.clear();
        for (int i = 0; i < Pc.outerSize(); ++i)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(Pc, i); it; ++it)
                if (it.row() == it.col() && it.value() > 0.5) m_active_c.push_back((int)it.row());
        m_active_f.clear();
        for (int i = 0; i < Pf.outerSize(); ++i)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(Pf, i); it; ++it)
                if (it.row() == it.col() && it.value() > 0.5) m_active_f.push_back((int)it.row());

        size_t nac = m_active_c.size();
        size_t naf = m_active_f.size();

        // BUGFIX (found via MLTS_DEBUG_M below): M2 = Mcinv*(Kcf*Pf) is
        // masked on the FACE side only (Pf), NOT on the cell side -- so a
        // COARSE cell adjacent to a fine-marked boundary face (assemble_P
        // marks a face fine if ANY neighbour cell is fine, even if the
        // OTHER neighbour is coarse) can have a nonzero M2 ROW there, even
        // though that cell itself is not in active_c. M1 has no such
        // issue (M1=Mcinv*Kcc*Pc IS masked on the cell/row side, via Kcc's
        // pure cell-diagonality, so coarse-cell rows of M1 are exactly
        // zero) -- confirmed empirically: MLTS_DEBUG_M showed M1/M2/M3
        // matching the reference exactly, M4 (the only one built from
        // Kfc*M2, which needs M2's cell-row sum to be complete) off by a
        // relative ~0.3% at N=2, enough to corrupt the whole solve.
        // Fix: widen the CELL side to "halo_c" = active_c UNION coarse
        // cells touching an active face (found via one pass over Kfc,
        // which is cheap/sparse -- NOT a full Mcinv/Kcc/Sffinv touch)
        // ONLY for the M4 computation; A_cf (from M2 restricted to
        // active_c rows only) is correct as-is since it only ever feeds
        // the active (fine) cells' own equations.
        std::unordered_set<int> active_f_set(m_active_f.begin(), m_active_f.end());
        std::vector<int> halo_c;
        {
            std::unordered_set<int> halo_c_set(m_active_c.begin(), m_active_c.end());
            const Eigen::SparseMatrix<T> & Kfc_full = erk_an.Kfc();
            for (int k = 0; k < Kfc_full.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(Kfc_full, k); it; ++it)
                    if (active_f_set.count((int)it.row())) { halo_c_set.insert(k); break; }
            halo_c.assign(halo_c_set.begin(), halo_c_set.end());
            std::sort(halo_c.begin(), halo_c.end());
        }

        // Extract the SMALL submatrices FIRST -- never touch a full
        // n_c/n_f-sized matrix or matrix product from here on.
        Eigen::SparseMatrix<T> Kcc_active   = extract(erk_an.Kcc(),   m_active_c, m_active_c);
        Eigen::SparseMatrix<T> Kcf_active   = extract(erk_an.Kcf(),   m_active_c, m_active_f);
        Eigen::SparseMatrix<T> Kfc_active   = extract(erk_an.Kfc(),   m_active_f, m_active_c);
        Eigen::SparseMatrix<T> Mcinv_active = extract(erk_an.invMc(), m_active_c, m_active_c);
        Eigen::SparseMatrix<T> Sffinv_active= extract(erk_an.SffInv(),m_active_f, m_active_f);

        // Wide (halo) versions -- needed ONLY for M4's cell-sum.
        Eigen::SparseMatrix<T> Kcf_halo   = extract(erk_an.Kcf(),   halo_c, m_active_f); // |halo_c| x naf
        Eigen::SparseMatrix<T> Mcinv_halo = extract(erk_an.invMc(), halo_c, halo_c);     // |halo_c| x |halo_c|, block-diag
        Eigen::SparseMatrix<T> Kfc_halo   = extract(erk_an.Kfc(),   m_active_f, halo_c); // naf x |halo_c|

        // Small matrix products only (sizes nac/naf/|halo_c|, NOT nc/nf).
        Eigen::SparseMatrix<T> M1_active = Mcinv_active * Kcc_active;                    // nac x nac
        Eigen::SparseMatrix<T> M2_active = Mcinv_active * Kcf_active;                    // nac x naf  (for A_cf)
        Eigen::SparseMatrix<T> M2_halo   = Mcinv_halo * Kcf_halo;                        // |halo_c| x naf (for M4 only)
        Eigen::SparseMatrix<T> M3_active = Sffinv_active * (Kfc_active * M1_active);     // naf x nac
        Eigen::SparseMatrix<T> M4_active = Sffinv_active * (Kfc_halo * M2_halo);         // naf x naf  (FIXED)

        if (std::getenv("MLTS_DEBUG_M")) {
            Eigen::SparseMatrix<T> M1_full = erk_an.invMc() * (erk_an.Kcc() * Pc);
            Eigen::SparseMatrix<T> M2_full = erk_an.invMc() * (erk_an.Kcf() * Pf);
            Eigen::SparseMatrix<T> M3_full = erk_an.SffInv() * (erk_an.Kfc() * M1_full);
            Eigen::SparseMatrix<T> M4_full = erk_an.SffInv() * (erk_an.Kfc() * M2_full);
            Eigen::SparseMatrix<T> M1_ref = extract(M1_full, m_active_c, m_active_c);
            Eigen::SparseMatrix<T> M2_ref = extract(M2_full, m_active_c, m_active_f);
            Eigen::SparseMatrix<T> M3_ref = extract(M3_full, m_active_f, m_active_c);
            Eigen::SparseMatrix<T> M4_ref = extract(M4_full, m_active_f, m_active_f);
            auto fro = [](const Eigen::SparseMatrix<T> & A) -> T {
                T s = 0; for (int k=0;k<A.outerSize();++k) for (typename Eigen::SparseMatrix<T>::InnerIterator it(A,k); it; ++it) s += it.value()*it.value();
                return std::sqrt(s);
            };
            Eigen::SparseMatrix<T> d1 = M1_active - M1_ref;
            Eigen::SparseMatrix<T> d2 = M2_active - M2_ref;
            Eigen::SparseMatrix<T> d3 = M3_active - M3_ref;
            Eigen::SparseMatrix<T> d4 = M4_active - M4_ref;
            std::cout << "[DEBUG_M] |M1_active-ref|=" << fro(d1) << " / |M1_ref|=" << fro(M1_ref)
                      << "  |M2diff|=" << fro(d2) << " / " << fro(M2_ref)
                      << "  |M3diff|=" << fro(d3) << " / " << fro(M3_ref)
                      << "  |M4diff|=" << fro(d4) << " / " << fro(M4_ref) << std::endl;
        }

        Eigen::SparseMatrix<T> A_cc_active(nac, nac);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M1_active.nonZeros() + nac);
            for (int k = 0; k < M1_active.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M1_active, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), gamma*dt*it.value());
            for (size_t i = 0; i < nac; i++) trips.emplace_back((int)i, (int)i, T(1));
            A_cc_active.setFromTriplets(trips.begin(), trips.end());
            A_cc_active.makeCompressed();
        }

        m_Acf_active = M2_active * (gamma*dt);
        m_Acf_active.makeCompressed();

        m_Afc_active = M3_active * (-gamma*dt);
        m_Afc_active.makeCompressed();

        Eigen::SparseMatrix<T> A_ff_active(naf, naf);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M4_active.nonZeros() + naf);
            for (int k = 0; k < M4_active.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M4_active, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), -gamma*dt*it.value());
            for (size_t i = 0; i < naf; i++) trips.emplace_back((int)i, (int)i, T(1));
            A_ff_active.setFromTriplets(trips.begin(), trips.end());
            A_ff_active.makeCompressed();
        }

        m_solver_Acc.analyzePattern(A_cc_active);
        m_solver_Acc.factorize(A_cc_active);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Acf_dense = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(m_Acf_active);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> X_active = m_solver_Acc.solve(Acf_dense);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AfcX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(m_Afc_active) * X_active;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S_dense = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(A_ff_active) - AfcX;

        Eigen::SparseMatrix<T> S_active = S_dense.sparseView();
        S_active.makeCompressed();

        m_solver_S.analyzePattern(S_active);
        m_solver_S.factorize(S_active);

        m_composed_Q = true;
    }

    Matrix<T, Dynamic, 1> solve(const Matrix<T, Dynamic, 1> & RHS) {

        Matrix<T, Dynamic, 1> result = RHS;

        size_t nac = m_active_c.size();
        size_t naf = m_active_f.size();
        size_t nc  = m_n_c_dof;

        Matrix<T, Dynamic, 1> rhs_c_active(nac);
        for (size_t i = 0; i < nac; i++) rhs_c_active(i) = RHS(m_active_c[i]);

        Matrix<T, Dynamic, 1> rhs_f_active(naf);
        for (size_t i = 0; i < naf; i++) rhs_f_active(i) = RHS(nc + m_active_f[i]);

        Matrix<T, Dynamic, 1> Acc_inv_rhsc = m_solver_Acc.solve(rhs_c_active);
        Matrix<T, Dynamic, 1> rhs_f_reduced = rhs_f_active - m_Afc_active * Acc_inv_rhsc;
        Matrix<T, Dynamic, 1> u_f_active = m_solver_S.solve(rhs_f_reduced);
        Matrix<T, Dynamic, 1> u_c_active = m_solver_Acc.solve(rhs_c_active - m_Acf_active * u_f_active);

        for (size_t i = 0; i < nac; i++) result(m_active_c[i]) = u_c_active(i);
        for (size_t i = 0; i < naf; i++) result(nc + m_active_f[i]) = u_f_active(i);

        return result;
    }

    bool composed() const { return m_composed_Q; }
    size_t n_active_c() const { return m_active_c.size(); }
    size_t n_active_f() const { return m_active_f.size(); }
};

#endif /* dirk_coupling_hho_scheme_condensed_local_v2_hpp */
