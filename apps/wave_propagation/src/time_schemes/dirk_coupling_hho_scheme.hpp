//
//  dirk_coupling_hho_scheme.hpp
//
//  LI-SDIRK(4,3) prototype: implicit SDIRK(3,4) fine-region time-stepping,
//  coupled with the EXISTING explicit LTS-ERK coarse-region "frozen Taylor
//  polynomial" precompute (erk_weight_LTS_coarse, unchanged, see
//  erk_coupling_hho_scheme.hpp).
//
//  This file does NOT modify erk_coupling_hho_scheme.hpp -- it only
//  #includes it and adds new, additive functionality:
//
//    1. dirk_coupling_hho_scheme<T> : builds and factorizes, ONCE per run
//       (gamma and dt are constant across the whole simulation, exactly
//       like dirk_hho_scheme::ComposeMatrix/DecomposeMatrix), the global
//       linear operator
//
//           A := I - gamma*dt*B*P
//
//       where B is the SAME homogeneous linear operator implicitly defined
//       by erk_coupling_hho_scheme::erk_weight (k_c = Mc_inv*(-Kcc*y_c -
//       Kcf*y_f), k_f = -Sff_inv*Kfc*k_c), and P = Pfine is the boolean
//       diagonal "fine region" projector built by
//       elastoacoustic_four_fields_assembler::assemble_P.
//
//       B*P block form (Pc,Pf = diagonal blocks of P restricted to cell/
//       face dofs):
//
//           (B*P)_cc = -Mc_inv*Kcc*Pc          (B*P)_cf = -Mc_inv*Kcf*Pf
//           (B*P)_fc =  Sff_inv*Kfc*Mc_inv*Kcc*Pc
//           (B*P)_ff =  Sff_inv*Kfc*Mc_inv*Kcf*Pf
//
//       so A is assembled directly as a sparse matrix and factorized once
//       with a general (non-symmetric) sparse LU. This is "option (b)" of
//       the task brief: the FULL global system, using the same B/erk_weight
//       machinery already validated for the explicit LTS scheme, rather
//       than a halo-restricted local solve (a later optimization -- see
//       TODO below).
//
//    2. dirk_weight_LTS_fine(...) : the new fine-region routine, replacing
//       erk_weight_LTS_fine's p-substep explicit RK4 loop with a SINGLE
//       implicit SDIRK(3,4) stage recursion of size dt (i=1,2,3):
//
//           tau_i = c_i*dt
//           (I - gamma*dt*B*P) Util^i = U^n + dt*sum_{l<i} a_il*k^l
//                                          + gamma*dt*Minv*P*F(tn+tau_i)
//                                          + gamma*dt*w(tau_i)
//           k^i = B*(P*Util^i) + Minv*P*F(tn+tau_i) + w(tau_i)
//
//       final update:  U^{n+1} = U^n + dt*sum_i b_i*k^i
//
//       w(tau) = w[0] + tau*w[1] + tau^2/2*w[2] + tau^3/6*w[3] is the SAME
//       frozen coarse Taylor polynomial produced by erk_weight_LTS_coarse
//       (called once per macro step, BEFORE this routine, exactly as in
//       ERK4_LTS.hpp -- reused verbatim, not touched here).
//

#ifndef dirk_coupling_hho_scheme_hpp
#define dirk_coupling_hho_scheme_hpp

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <functional>
#include <vector>

#include "erk_coupling_hho_scheme.hpp"

template<typename T>
class dirk_coupling_hho_scheme {

    private:

    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> m_solver;
    size_t m_n_c_dof;
    size_t m_n_f_dof;
    T m_gamma;
    T m_dt;
    bool m_composed_Q;

    public:

    dirk_coupling_hho_scheme() : m_n_c_dof(0), m_n_f_dof(0), m_gamma(0.0), m_dt(0.0), m_composed_Q(false) {}

    // Builds A = I - gamma*dt*B*P and factorizes it ONCE. erk_an must
    // already have Mcc_inverse()/Sff_inverse() computed (as in every
    // existing LTS driver). Pfine is the boolean diagonal fine-region
    // projector (assembler.Pfine, size n_c_dof+n_f_dof).
    //
    // TODO (optimization, not correctness): this builds/factorizes the
    // FULL global (n_c_dof+n_f_dof) system. A fully localized version
    // restricted to the fine region + its immediate halo (mirroring
    // erk_weight_LTS_coarse's m_coarse_active_dofs masking, done for P
    // instead of (I-P)) would be cheaper for a genuinely local fine
    // region, but correctness is established first here; the full
    // global operator gives the exact same k^i as any exact local
    // reduction would, since B itself is never truncated.
    void compose(erk_coupling_hho_scheme<T> & erk_an,
                 const Eigen::SparseMatrix<T> & Pfine,
                 T gamma,
                 T dt) {

        m_gamma = gamma;
        m_dt    = dt;
        m_n_c_dof = erk_an.n_c_dof();
        m_n_f_dof = erk_an.n_f_dof();
        size_t nc = m_n_c_dof;
        size_t nf = m_n_f_dof;
        size_t n  = nc + nf;

        Eigen::SparseMatrix<T> Pc = Pfine.block(0, 0, nc, nc);
        Eigen::SparseMatrix<T> Pf = Pfine.block(nc, nc, nf, nf);

        Eigen::SparseMatrix<T> & Kcc = erk_an.Kcc();
        Eigen::SparseMatrix<T> & Kcf = erk_an.Kcf();
        Eigen::SparseMatrix<T> & Kfc = erk_an.Kfc();
        Eigen::SparseMatrix<T> & Mcinv = erk_an.invMc();
        Eigen::SparseMatrix<T> & Sffinv = erk_an.SffInv();

        // M1 = Mc_inv*Kcc*Pc  (nc x nc)
        // M2 = Mc_inv*Kcf*Pf  (nc x nf)
        // M3 = Sff_inv*Kfc*M1 (nf x nc)
        // M4 = Sff_inv*Kfc*M2 (nf x nf)
        Eigen::SparseMatrix<T> M1 = Mcinv * (Kcc * Pc);
        Eigen::SparseMatrix<T> M2 = Mcinv * (Kcf * Pf);
        Eigen::SparseMatrix<T> M3 = Sffinv * (Kfc * M1);
        Eigen::SparseMatrix<T> M4 = Sffinv * (Kfc * M2);

        // A = I - gamma*dt*(B*P), block form:
        //   A_cc =  I_c + gamma*dt*M1     A_cf =  gamma*dt*M2
        //   A_fc = -gamma*dt*M3           A_ff =  I_f - gamma*dt*M4
        std::vector<Eigen::Triplet<T>> trips;
        trips.reserve(M1.nonZeros() + M2.nonZeros() + M3.nonZeros() + M4.nonZeros() + n);

        for (int k = 0; k < M1.outerSize(); ++k)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(M1, k); it; ++it)
                trips.emplace_back((int)it.row(), (int)it.col(), gamma*dt*it.value());

        for (int k = 0; k < M2.outerSize(); ++k)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(M2, k); it; ++it)
                trips.emplace_back((int)it.row(), (int)(nc+it.col()), gamma*dt*it.value());

        for (int k = 0; k < M3.outerSize(); ++k)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(M3, k); it; ++it)
                trips.emplace_back((int)(nc+it.row()), (int)it.col(), -gamma*dt*it.value());

        for (int k = 0; k < M4.outerSize(); ++k)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(M4, k); it; ++it)
                trips.emplace_back((int)(nc+it.row()), (int)(nc+it.col()), -gamma*dt*it.value());

        for (size_t i = 0; i < n; i++)
            trips.emplace_back((int)i, (int)i, T(1));

        Eigen::SparseMatrix<T> A(n, n);
        A.setFromTriplets(trips.begin(), trips.end());
        A.makeCompressed();

        m_solver.analyzePattern(A);
        m_solver.factorize(A);
        m_composed_Q = true;
    }

    Matrix<T, Dynamic, 1> solve(const Matrix<T, Dynamic, 1> & RHS) {
        return m_solver.solve(RHS);
    }

    bool composed() const { return m_composed_Q; }
};


// =====================================================================
// dirk_weight_LTS_fine : the new implicit fine-region routine.
//
// Signature deliberately mirrors erk_weight_LTS_fine's role (mutate
// x_dof_n IN PLACE from U^n to U^{n+1}), but takes:
//   - erk_an       : the SAME erk_coupling_hho_scheme used for the coarse
//                     precompute (reused for B = erk_weight)
//   - dirk_op      : the pre-factorized (I - gamma*dt*B*P) operator: any
//                     class exposing solve(RHS) -> full (n_c+n_f) vector
//                     (dirk_coupling_hho_scheme above, or its condensed /
//                     halo-restricted twins dirk_coupling_hho_scheme_
//                     condensed, _condensed_local, _condensed_local_v2)
//   - Pfine        : boolean diagonal fine-region projector
//   - w            : frozen coarse Taylor coefficients from
//                     erk_weight_LTS_coarse (called ONCE per step, BEFORE
//                     this routine, unchanged)
//   - tn, dt       : macro step start time / size (single global step,
//                     NO p/dtau substepping -- that is the whole point of
//                     LI-SDIRK)
//   - eval_F       : evaluates the FULL global RHS force vector at an
//                     absolute time (same lambda pattern as in
//                     ERK4_LTS.hpp's eval_F)
//   - a_tab,b_tab,c_tab : SDIRK(3,4) Butcher tableau, e.g. from
//                     dirk_butcher_tableau::sdirk_tables(3, a, b, c)
//   - gamma        : the SDIRK diagonal entry (a_tab(i,i) for all i)
// =====================================================================
template<typename T, typename DirkOp>
void dirk_weight_LTS_fine(erk_coupling_hho_scheme<T> & erk_an,
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

    // Stages (time_schemes/lts.hpp, LI slopes through the operator):
    //   RHS  = U^n + dt*sum_{l<i} a_il*k^l + gamma*dt*(Minv*P*F(tn+tau_i) + w(tau_i))
    //   Util = (I - gamma*dt*B*P)^{-1} RHS
    //   k^i  = B*(P*Util) + Minv*P*F(tn+tau_i) + w(tau_i)
    // U^{n+1} = U^n + dt*sum_i b_i*k^i
    typename erk_coupling_hho_scheme<T>::lts_operator op(erk_an, nullptr, &Pfine);
    local_time_stepping::li_step(op, x_dof_n, w, tn, dt, a_tab, b_tab, c_tab, gamma*dt, eval_F,
                                 [&](const Matrix<T, Dynamic, 1> & RHS) { return Matrix<T, Dynamic, 1>(dirk_op.solve(RHS)); });
}

#endif /* dirk_coupling_hho_scheme_hpp */
