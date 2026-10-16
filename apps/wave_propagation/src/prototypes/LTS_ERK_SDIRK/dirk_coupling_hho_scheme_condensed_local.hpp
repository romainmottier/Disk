//
//  dirk_coupling_hho_scheme_condensed_local.hpp
//
//  Third, SEPARATE prototype: HALO-RESTRICTED statically-condensed LI-SDIRK
//  operator. Neither dirk_coupling_hho_scheme.hpp (full global factorization)
//  nor dirk_coupling_hho_scheme_condensed.hpp (plain condensation, full n_f
//  Schur complement) is modified -- this is a new, additive twin, same
//  public interface (compose(erk_an,Pfine,gamma,dt) + solve(RHS) returning
//  the FULL (n_c+n_f) vector).
//
//  Why this is needed: dirk_coupling_hho_scheme_condensed already reduces
//  A_cc away (cheap, block-diagonal), but then assembles/factorizes the
//  FULL n_f x n_f Schur complement S = A_ff - A_fc*A_cc^{-1}*A_cf. Outside
//  the "active" faces -- those touching at least one fine cell -- S is
//  PROVABLY just the identity:
//    - A_ff = I_f - gamma*dt*Sffinv*Kfc*Mcinv*Kcf*Pf is I wherever Pf=0
//      (Pf zeroes the whole M2=Mcinv*Kcf*Pf term there).
//    - A_fc*A_cc^{-1}*A_cf's row f is nonzero only if face f touches a
//      cell c with Pc_c=1 (a fine cell) -- i.e. only for faces touching a
//      fine cell, which is EXACTLY assemble_P's own definition of Pf=1
//      (assemble_P's own "mark faces fine if any neighbour cell is fine"
//      step, elastoacoustic_four_fields_assembler.hpp).
//  So "active faces" = support of Pf, and factorizing the full n_f system
//  wastes almost all of the work factorizing a matrix that's provably
//  just I outside that support. Same story for cells: A_cc is I wherever
//  Pc=0 (coarse cells), so only active (fine) cells need a real solve.
//
//  compose() therefore:
//    1. identifies active_c = support of Pc (fine cell dofs), active_f =
//       support of Pf (fine face dofs) -- mirrors erk_weight_LTS_coarse's
//       own m_coarse_active_dofs masking pattern, just for P instead of
//       (I-P);
//    2. extracts the SMALL active-only submatrices A_cc, A_cf, A_fc, A_ff
//       (sizes = #active_c, #active_f -- the fine region + its halo, NOT
//       n_c/n_f);
//    3. factorizes those small submatrices only, and builds/factorizes
//       the correspondingly SMALL Schur complement S_active
//       (#active_f x #active_f).
//    4. Batches X_active = A_cc_active^{-1}*A_cf_active as ONE dense
//       multi-RHS solve (feasible now that both dimensions are small,
//       instead of the plain condensed version's column-by-column loop).
//
//  solve(RHS) starts from `result = RHS` (handles the "identity
//  everywhere outside the active set" case for both cells AND faces with
//  no extra bookkeeping), then overwrites only the active entries via the
//  small Schur back-substitution.
//

#ifndef dirk_coupling_hho_scheme_condensed_local_hpp
#define dirk_coupling_hho_scheme_condensed_local_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <functional>
#include <vector>
#include <unordered_map>

#include "../../time_schemes/dirk_coupling_hho_scheme.hpp" // erk_coupling_hho_scheme + dirk_weight_LTS_fine (any DirkOp)

template<typename T>
class dirk_coupling_hho_scheme_condensed_local {

    private:

    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_Acc;
    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_S;
    Eigen::SparseMatrix<T> m_Acf_active, m_Afc_active;
    std::vector<int> m_active_c;  // global cell-dof indices (into [0,n_c))
    std::vector<int> m_active_f;  // global face-dof indices (into [0,n_f)), offset by n_c to become global
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

    dirk_coupling_hho_scheme_condensed_local() : m_n_c_dof(0), m_n_f_dof(0), m_composed_Q(false) {}

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

        // Active sets = support of Pc / Pf (fine cells / faces touching
        // at least one fine cell -- assemble_P's own definition).
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

        Eigen::SparseMatrix<T> & Kcc = erk_an.Kcc();
        Eigen::SparseMatrix<T> & Kcf = erk_an.Kcf();
        Eigen::SparseMatrix<T> & Kfc = erk_an.Kfc();
        Eigen::SparseMatrix<T> & Mcinv = erk_an.invMc();
        Eigen::SparseMatrix<T> & Sffinv = erk_an.SffInv();

        Eigen::SparseMatrix<T> M1 = Mcinv * (Kcc * Pc);   // nc x nc
        Eigen::SparseMatrix<T> M2 = Mcinv * (Kcf * Pf);   // nc x nf
        Eigen::SparseMatrix<T> M3 = Sffinv * (Kfc * M1);  // nf x nc
        Eigen::SparseMatrix<T> M4 = Sffinv * (Kfc * M2);  // nf x nf

        // Extract the SMALL active-only blocks and add the identity
        // contribution (only meaningful on the active diagonal -- outside
        // the active set A_cc/A_ff are identically I, handled implicitly
        // in solve() by starting from RHS itself).
        Eigen::SparseMatrix<T> M1_active = extract(M1, m_active_c, m_active_c);
        Eigen::SparseMatrix<T> A_cf_active_raw = extract(M2, m_active_c, m_active_f);
        Eigen::SparseMatrix<T> M3_active = extract(M3, m_active_f, m_active_c);
        Eigen::SparseMatrix<T> M4_active = extract(M4, m_active_f, m_active_f);

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

        m_Acf_active = A_cf_active_raw * gamma * dt;
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

        // X_active = A_cc_active^{-1} * A_cf_active -- ONE batched dense
        // multi-RHS solve (both dims now small: nac x naf).
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

        // Outside the active set, A is exactly I -- so start from a copy
        // of RHS (handles both coarse cells and non-fine-adjacent faces
        // with zero extra work), then overwrite only the active entries.
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

#endif /* dirk_coupling_hho_scheme_condensed_local_hpp */
