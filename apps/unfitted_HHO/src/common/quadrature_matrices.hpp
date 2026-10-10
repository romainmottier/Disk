/*
 * Matrices of the values of the bases at the points of a quadrature, so that
 * the integrals of the local operators are computed by matrix products
 * (B1^T diag(w) B2) instead of rank-one updates point by point.
 */

#pragma once

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>

template<typename T>
using qm_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using qm_vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/* scalar type and point type of a quadrature point */
template<typename QP>
using qp_scalar = std::decay_t<decltype(std::declval<QP>().weight())>;

template<typename QP>
using qp_point = std::decay_t<decltype(std::declval<QP>().point())>;

/* Weights of a quadrature */
template<typename QP>
qm_vector<qp_scalar<QP>>
quadrature_weights(const std::vector<QP>& qps)
{
    qm_vector<qp_scalar<QP>> w(qps.size());
    for (size_t q = 0; q < qps.size(); q++)
        w(q) = qps[q].weight();
    return w;
}

/* Values of a scalar basis: one row per quadrature point */
template<typename Basis, typename QP>
qm_matrix<qp_scalar<QP>>
basis_values(const Basis& b, const std::vector<QP>& qps)
{
    if constexpr (requires(const std::vector<qp_point<QP>>& p) { b.eval_basis_batch(p); }) {
        std::vector<qp_point<QP>> pts(qps.size());
        for (size_t q = 0; q < qps.size(); q++)
            pts[q] = qps[q].point();
        return b.eval_basis_batch(pts);
    }
    else {
        qm_matrix<qp_scalar<QP>> V(qps.size(), b.size());
        for (size_t q = 0; q < qps.size(); q++)
            V.row(q) = b.eval_basis(qps[q].point()).transpose();
        return V;
    }
}

/* Derivatives of a scalar basis: one matrix per direction, one row per point */
template<size_t DIM, typename Basis, typename QP>
std::array<qm_matrix<qp_scalar<QP>>, DIM>
basis_gradients(const Basis& b, const std::vector<QP>& qps)
{
    if constexpr (requires(const std::vector<qp_point<QP>>& p) { b.eval_gradients_batch(p); }) {
        std::vector<qp_point<QP>> pts(qps.size());
        for (size_t q = 0; q < qps.size(); q++)
            pts[q] = qps[q].point();
        return b.eval_gradients_batch(pts);
    }
    else {
        std::array<qm_matrix<qp_scalar<QP>>, DIM> G;
        for (size_t d = 0; d < DIM; d++)
            G[d].resize(qps.size(), b.size());
        for (size_t q = 0; q < qps.size(); q++) {
            auto g = b.eval_gradients(qps[q].point());
            for (size_t d = 0; d < DIM; d++)
                G[d].row(q) = g.col(d).transpose();
        }
        return G;
    }
}

/* B1^T diag(w) B2 */
template<typename T>
qm_matrix<T>
weighted_product(const qm_matrix<T>& B1, const qm_vector<T>& w, const qm_matrix<T>& B2)
{
    return B1.transpose() * (w.asDiagonal() * B2);
}

/* dst(DIM a + d, col0 + j) += scale S(a, j): the rows of the vector basis
 * whose functions are those of the scalar basis in the component d */
template<size_t DIM, typename T>
void
add_component_rows(qm_matrix<T>& dst, size_t col0, const qm_matrix<T>& S, size_t d, T scale = 1)
{
    for (Eigen::Index a = 0; a < S.rows(); a++)
        dst.block(DIM*a + d, col0, 1, S.cols()) += scale * S.row(a);
}

/* Mass matrix of the vector basis from the mass matrix M of its scalar basis */
template<size_t DIM, typename T>
qm_matrix<T>
vector_mass(const qm_matrix<T>& M)
{
    qm_matrix<T> ret = qm_matrix<T>::Zero(DIM*M.rows(), DIM*M.cols());
    for (Eigen::Index a = 0; a < M.rows(); a++)
        for (Eigen::Index b = 0; b < M.cols(); b++)
            for (size_t d = 0; d < DIM; d++)
                ret(DIM*a + d, DIM*b + d) = M(a, b);
    return ret;
}

/* Solve (vector_mass(M)) X = R, the components being independent */
template<size_t DIM, typename T>
qm_matrix<T>
solve_vector_mass(const qm_matrix<T>& M, const qm_matrix<T>& R)
{
    Eigen::LDLT<qm_matrix<T>> ldlt(M);
    qm_matrix<T> X(R.rows(), R.cols());
    qm_matrix<T> Rd(M.rows(), R.cols());
    for (size_t d = 0; d < DIM; d++) {
        for (Eigen::Index a = 0; a < M.rows(); a++)
            Rd.row(a) = R.row(DIM*a + d);
        qm_matrix<T> Xd = ldlt.solve(Rd);
        for (Eigen::Index a = 0; a < M.rows(); a++)
            X.row(DIM*a + d) = Xd.row(a);
    }
    return X;
}
