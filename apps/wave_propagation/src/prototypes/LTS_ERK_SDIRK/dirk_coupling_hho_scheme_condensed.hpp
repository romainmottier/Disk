//
//  dirk_coupling_hho_scheme_condensed.hpp
//
//  Second, SEPARATE prototype: statically-condensed version of
//  dirk_coupling_hho_scheme (that file is NOT modified -- this is a new,
//  additive twin). Same public interface (compose(erk_an,Pfine,gamma,dt)
//  + solve(RHS) returning the FULL (n_c+n_f) vector), so it is a drop-in
//  replacement for dirk_coupling_hho_scheme from the caller's point of
//  view, but internally eliminates the fine CELL unknowns first (paper
//  Sec. 5.2), instead of factorizing the full (n_c+n_f) system directly.
//
//  Why the cell block condenses cleanly (see dirk_coupling_hho_scheme.hpp
//  for the A_cc/A_cf/A_fc/A_ff block definitions, reused verbatim here):
//
//      A_cc = I_c + gamma*dt*Mc_inv*Kcc*Pc
//
//  Mc_inv and Kcc are BOTH exactly block-diagonal per cell (standard HHO
//  structure -- no direct cell-cell coupling, all coupling goes through
//  faces via Kcf/Kfc), and Pc is diagonal (0/1 per cell dof), so A_cc is
//  itself EXACTLY block-diagonal, one small dense block per cell (trivial
//  I block for coarse cells where Pc_i=0). Consequences exploited below:
//
//    - A_cc is cheap to factorize (SparseLU on a block-diagonal matrix
//      has zero fill-in beyond the existing blocks -- algorithmically
//      identical cost to inverting each cell block separately).
//    - Because A_cc is block-diagonal, solving A_cc*x=b for a RHS b
//      supported on a SINGLE cell's dof range gives an x supported on
//      that SAME cell's dof range ONLY (no spread to other cells). Since
//      A_cf's columns (one per face dof) are nonzero only on the 1-2
//      cells touching that face, X := A_cc^{-1}*A_cf keeps EXACTLY the
//      same sparsity pattern as A_cf -- just rescaled. Built here via a
//      column-by-column sparse solve (looping only over A_cf's nonzero
//      columns -- i.e. only fine-region-adjacent faces, a small subset
//      of all faces by construction of the coarse/fine LTS split).
//
//  Standard Schur complement, using the SAME A_cc/A_cf/A_fc/A_ff block
//  definitions as dirk_coupling_hho_scheme::compose (just not assembled
//  into one big matrix):
//
//      u_c = A_cc^{-1}*(rhs_c - A_cf*u_f)                         [cheap]
//      (A_ff - A_fc*A_cc^{-1}*A_cf) u_f = rhs_f - A_fc*A_cc^{-1}*rhs_c
//                                                        [reduced n_f x n_f]
//
//  S := A_ff - A_fc*A_cc^{-1}*A_cf is built ONCE (gamma,dt constant) and
//  factorized ONCE, exactly mirroring dirk_coupling_hho_scheme's
//  "compose once per run" pattern.
//
//  Plain (non-halo-restricted) version: the Schur solve is over ALL n_f
//  faces. A halo-restricted reduction (mirroring erk_weight_LTS_coarse's
//  m_coarse_active_dofs masking) is a possible further optimization, not
//  implemented here.
//

#ifndef dirk_coupling_hho_scheme_condensed_hpp
#define dirk_coupling_hho_scheme_condensed_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <functional>
#include <vector>

#include "../../common/erk_coupling_hho_scheme.hpp"

template<typename T>
class dirk_coupling_hho_scheme_condensed {

    private:

    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_Acc;
    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver_S;
    Eigen::SparseMatrix<T> m_Acf, m_Afc; // kept for use in solve()
    size_t m_n_c_dof;
    size_t m_n_f_dof;
    bool m_composed_Q;

    public:

    dirk_coupling_hho_scheme_condensed() : m_n_c_dof(0), m_n_f_dof(0), m_composed_Q(false) {}

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

        Eigen::SparseMatrix<T> & Kcc = erk_an.Kcc();
        Eigen::SparseMatrix<T> & Kcf = erk_an.Kcf();
        Eigen::SparseMatrix<T> & Kfc = erk_an.Kfc();
        Eigen::SparseMatrix<T> & Mcinv = erk_an.invMc();
        Eigen::SparseMatrix<T> & Sffinv = erk_an.SffInv();

        // Same M1..M4 building blocks as dirk_coupling_hho_scheme::compose.
        Eigen::SparseMatrix<T> M1 = Mcinv * (Kcc * Pc);   // nc x nc
        Eigen::SparseMatrix<T> M2 = Mcinv * (Kcf * Pf);   // nc x nf
        Eigen::SparseMatrix<T> M3 = Sffinv * (Kfc * M1);  // nf x nc
        Eigen::SparseMatrix<T> M4 = Sffinv * (Kfc * M2);  // nf x nf

        // A_cc = I_c + gamma*dt*M1   (EXACTLY block-diagonal per cell)
        Eigen::SparseMatrix<T> A_cc(nc, nc);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M1.nonZeros() + nc);
            for (int k = 0; k < M1.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M1, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), gamma*dt*it.value());
            for (size_t i = 0; i < nc; i++)
                trips.emplace_back((int)i, (int)i, T(1));
            A_cc.setFromTriplets(trips.begin(), trips.end());
            A_cc.makeCompressed();
        }

        // A_cf = gamma*dt*M2  (nc x nf)
        Eigen::SparseMatrix<T> A_cf(nc, nf);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M2.nonZeros());
            for (int k = 0; k < M2.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M2, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), gamma*dt*it.value());
            A_cf.setFromTriplets(trips.begin(), trips.end());
            A_cf.makeCompressed();
        }

        // A_fc = -gamma*dt*M3  (nf x nc)
        Eigen::SparseMatrix<T> A_fc(nf, nc);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M3.nonZeros());
            for (int k = 0; k < M3.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M3, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), -gamma*dt*it.value());
            A_fc.setFromTriplets(trips.begin(), trips.end());
            A_fc.makeCompressed();
        }

        // A_ff = I_f - gamma*dt*M4  (nf x nf)
        Eigen::SparseMatrix<T> A_ff(nf, nf);
        {
            std::vector<Eigen::Triplet<T>> trips;
            trips.reserve(M4.nonZeros() + nf);
            for (int k = 0; k < M4.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(M4, k); it; ++it)
                    trips.emplace_back((int)it.row(), (int)it.col(), -gamma*dt*it.value());
            for (size_t i = 0; i < nf; i++)
                trips.emplace_back((int)i, (int)i, T(1));
            A_ff.setFromTriplets(trips.begin(), trips.end());
            A_ff.makeCompressed();
        }

        // Factorize A_cc ONCE -- block-diagonal, so this is exactly as
        // cheap as factorizing each cell's small dense block separately.
        m_solver_Acc.analyzePattern(A_cc);
        m_solver_Acc.factorize(A_cc);

        // X = A_cc^{-1} * A_cf, column-by-column over A_cf's NONZERO
        // columns only (i.e. only fine-region-adjacent face dofs -- a
        // small subset of all nf faces). Because A_cc is block-diagonal,
        // X keeps EXACTLY A_cf's sparsity pattern (no fill-in): a RHS
        // supported on one cell's block stays supported on that same
        // block after solving. Values below a tight relative threshold
        // are dropped to reconstruct that sparsity exactly (up to
        // round-off) instead of storing a dense nc-length result per
        // column.
        std::vector<Eigen::Triplet<T>> X_trips;
        {
            Eigen::Matrix<T, Eigen::Dynamic, 1> col_dense(nc);
            for (int j = 0; j < A_cf.outerSize(); ++j) {
                typename Eigen::SparseMatrix<T>::InnerIterator it0(A_cf, j);
                if (!it0) continue; // empty column -- skip entirely
                col_dense.setZero();
                T col_max = 0;
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A_cf, j); it; ++it) {
                    col_dense(it.row()) = it.value();
                    col_max = std::max(col_max, std::abs(it.value()));
                }
                Eigen::Matrix<T, Eigen::Dynamic, 1> x = m_solver_Acc.solve(col_dense);
                T xmax = x.cwiseAbs().maxCoeff();
                T thresh = std::max(xmax, col_max) * T(1e-13);
                for (int i = 0; i < x.rows(); i++)
                    if (std::abs(x(i)) > thresh)
                        X_trips.emplace_back(i, j, x(i));
            }
        }
        Eigen::SparseMatrix<T> X(nc, nf);
        X.setFromTriplets(X_trips.begin(), X_trips.end());
        X.makeCompressed();

        // S = A_ff - A_fc*X  (nf x nf), factorized ONCE.
        Eigen::SparseMatrix<T> AfcX = A_fc * X;
        Eigen::SparseMatrix<T> S = A_ff - AfcX;
        S.makeCompressed();

        m_solver_S.analyzePattern(S);
        m_solver_S.factorize(S);

        m_Acf = A_cf;
        m_Afc = A_fc;
        m_composed_Q = true;
    }

    Matrix<T, Dynamic, 1> solve(const Matrix<T, Dynamic, 1> & RHS) {

        size_t nc = m_n_c_dof;
        size_t nf = m_n_f_dof;

        Matrix<T, Dynamic, 1> rhs_c = RHS.head(nc);
        Matrix<T, Dynamic, 1> rhs_f = RHS.tail(nf);

        // u_c_partial = A_cc^{-1}*rhs_c
        Matrix<T, Dynamic, 1> Acc_inv_rhsc = m_solver_Acc.solve(rhs_c);

        // reduced face RHS = rhs_f - A_fc*(A_cc^{-1}*rhs_c)
        Matrix<T, Dynamic, 1> rhs_f_reduced = rhs_f - m_Afc * Acc_inv_rhsc;

        // Solve the condensed face-only system
        Matrix<T, Dynamic, 1> u_f = m_solver_S.solve(rhs_f_reduced);

        // Back-substitute for the cell unknowns
        Matrix<T, Dynamic, 1> u_c = m_solver_Acc.solve(rhs_c - m_Acf * u_f);

        Matrix<T, Dynamic, 1> result(nc + nf);
        result.head(nc) = u_c;
        result.tail(nf) = u_f;
        return result;
    }

    bool composed() const { return m_composed_Q; }
};


// =====================================================================
// dirk_weight_LTS_fine_condensed : EXACT twin of dirk_weight_LTS_fine
// (dirk_coupling_hho_scheme.hpp, untouched), just templated on the
// operator type so it can take a dirk_coupling_hho_scheme_condensed<T>
// (or, generically, anything exposing the same solve(RHS) interface)
// instead of dirk_coupling_hho_scheme<T>. Algorithm is byte-for-byte
// identical -- same coarse Taylor term usage, same SDIRK(3,4) stage
// recursion, same final combination.
// =====================================================================
template<typename T, typename DirkOp>
void dirk_weight_LTS_fine_condensed(erk_coupling_hho_scheme<T> & erk_an,
                                     DirkOp & dirk_op,
                                     Matrix<T, Dynamic, 1> & x_dof_n,
                                     const Eigen::SparseMatrix<T> & Pfine,
                                     const std::vector<Matrix<T, Dynamic, 1>> & w,
                                     const T tn,
                                     const T dt,
                                     const std::function<Matrix<T, Dynamic, 1>(T)> & eval_F,
                                     const Matrix<double, Dynamic, Dynamic> & a_tab,
                                     const Matrix<double, Dynamic, 1> & b_tab,
                                     const Matrix<double, Dynamic, 1> & c_tab,
                                     const T gamma) {

    const size_t s = 3; // SDIRK(3,4)
    const size_t ndof = x_dof_n.rows();

    auto Taylor_w = [&](T tau) -> Matrix<T, Dynamic, 1> {
        T tau2 = tau*tau, tau3 = tau*tau2;
        return w[0] + tau*w[1] + (tau2/2.0)*w[2] + (tau3/6.0)*w[3];
    };

    Matrix<T, Dynamic, 1> Un = x_dof_n;
    std::vector<Matrix<T, Dynamic, 1>> k(s);

    for (size_t i = 0; i < s; i++) {

        T tau_i = (T)c_tab(i,0) * dt;

        Matrix<T, Dynamic, 1> Ffull = eval_F(tn + tau_i);
        Matrix<T, Dynamic, 1> PF    = Pfine * Ffull;

        erk_an.SetFg(PF);
        Matrix<T, Dynamic, 1> zero = Matrix<T, Dynamic, 1>::Zero(ndof);
        Matrix<T, Dynamic, 1> MinvPF;
        erk_an.erk_weight(zero, MinvPF);
        erk_an.ZeroFc();

        Matrix<T, Dynamic, 1> w_tau = Taylor_w(tau_i);

        Matrix<T, Dynamic, 1> RHS = Un;
        for (size_t l = 0; l < i; l++)
            RHS += dt * (T)a_tab(i,l) * k[l];
        RHS += gamma*dt*MinvPF + gamma*dt*w_tau;

        Matrix<T, Dynamic, 1> Util = dirk_op.solve(RHS);

        Matrix<T, Dynamic, 1> Py = Pfine * Util;
        erk_an.SetFg(PF);
        Matrix<T, Dynamic, 1> ki;
        erk_an.erk_weight(Py, ki);
        erk_an.ZeroFc();
        ki += w_tau;

        k[i] = ki;
    }

    Matrix<T, Dynamic, 1> Unp1 = Un;
    for (size_t i = 0; i < s; i++)
        Unp1 += dt * (T)b_tab(i,0) * k[i];

    x_dof_n = Unp1;
}

#endif /* dirk_coupling_hho_scheme_condensed_hpp */
