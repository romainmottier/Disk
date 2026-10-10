/*
 * Static condensation of the cell unknowns of the unfitted HHO methods
 * (derivation in doc/static_condensation.tex).
 *
 * local_condensation: the cell unknowns are eliminated patch by patch before
 * any global assembly. A patch is a set of mesh cells whose cell unknowns
 * are coupled: one cell (or agglomerated cell), both sides of a cut cell, and,
 * with the polynomial extension, a cell T together with the ill-cut cells
 * S in N^-1(T) (connected components of the relation S ~ N(S)). The local
 * contributions of the elements of a patch are summed in a dense patch
 * matrix, its cell unknowns are eliminated, and only the Schur complement on
 * the face unknowns is assembled.
 *
 * algebraic_condensation: the same elimination performed on the assembled
 * global system. The global system is assembled with the cell unknowns first (n_c rows),
 * then the face unknowns. The cell unknowns are coupled to each other only
 * in small clusters: one cell (or agglomerated cell), the two sides of a cut
 * cell (interface terms), a cell and the ill-cut cells it stabilizes (s^N and
 * extended gradient of the polynomial extension). Each cluster K is
 * eliminated with a dense factorization of A_KK, which leaves the Schur
 * complement on the face unknowns
 *      S = A_FF - sum_K A_FK A_KK^-1 A_KF,   g = b_F - sum_K A_FK A_KK^-1 b_K.
 * The clusters are the connected components of the graph of A_CC, so that
 * the procedure applies to both stabilizations without any knowledge of the
 * pairing or of the agglomeration.
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "assembly_index.hpp"

template<typename T>
class algebraic_condensation {

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    size_t                           n_c, n_f;
    std::vector<std::vector<size_t>> cluster_dofs;  // cell unknowns of each cluster
    std::vector<std::vector<size_t>> cluster_faces; // face unknowns coupled to it (numbered from 0)
    std::vector<matrix_type>         X;             // A_KK^-1 A_KF
    std::vector<vector_type>         y;             // A_KK^-1 b_K

    static size_t find(std::vector<size_t>& parent, size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }

public:
    Eigen::SparseMatrix<T> S;   // Schur complement on the face unknowns
    vector_type            g;   // condensed right-hand side
    size_t                 max_cluster = 0;

    algebraic_condensation(const Eigen::SparseMatrix<T>& A, const vector_type& b, size_t num_cell_dofs)
        : n_c(num_cell_dofs), n_f(A.rows() - num_cell_dofs)
    {
        /* clusters: connected components of the graph of A_CC (union-find) */
        std::vector<size_t> parent(n_c);
        std::iota(parent.begin(), parent.end(), 0);
        for (Eigen::Index j = 0; j < Eigen::Index(n_c); j++)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it)
                if (size_t(it.row()) < n_c) {
                    auto a = find(parent, it.row()), c = find(parent, j);
                    if (a != c)
                        parent[std::max(a, c)] = std::min(a, c);
                }
        std::vector<int> cluster_of(n_c, -1);
        for (size_t i = 0; i < n_c; i++) {
            auto r = find(parent, i);
            if (cluster_of[r] < 0) {
                cluster_of[r] = cluster_dofs.size();
                cluster_dofs.push_back({});
            }
            cluster_of[i] = cluster_of[r];
            cluster_dofs[cluster_of[i]].push_back(i);
        }

        /* face unknowns coupled to each cluster */
        cluster_faces.resize(cluster_dofs.size());
        for (Eigen::Index j = 0; j < Eigen::Index(n_c); j++)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it)
                if (size_t(it.row()) >= n_c)
                    cluster_faces[cluster_of[j]].push_back(it.row() - n_c);
        for (auto& f : cluster_faces) {
            std::sort(f.begin(), f.end());
            f.erase(std::unique(f.begin(), f.end()), f.end());
        }

        /* elimination of the clusters */
        std::vector<Eigen::Triplet<T>> trip;
        for (Eigen::Index j = 0; j < Eigen::Index(n_f); j++)
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, n_c + j); it; ++it)
                if (size_t(it.row()) >= n_c)
                    trip.emplace_back(it.row() - n_c, j, it.value());
        g = b.tail(n_f);

        X.resize(cluster_dofs.size());
        y.resize(cluster_dofs.size());
        for (size_t K = 0; K < cluster_dofs.size(); K++) {
            const auto& cd = cluster_dofs[K];
            const auto& cf = cluster_faces[K];
            max_cluster = std::max(max_cluster, cd.size());
            matrix_type AKK(cd.size(), cd.size()), AKF(cd.size(), cf.size());
            vector_type bK(cd.size());
            for (size_t r = 0; r < cd.size(); r++) {
                bK(r) = b(cd[r]);
                for (size_t c = 0; c < cd.size(); c++)
                    AKK(r, c) = A.coeff(cd[r], cd[c]);
                for (size_t c = 0; c < cf.size(); c++)
                    AKF(r, c) = A.coeff(cd[r], n_c + cf[c]);
            }
            Eigen::LDLT<matrix_type> ldlt(AKK);
            X[K] = ldlt.solve(AKF);
            y[K] = ldlt.solve(bK);
            matrix_type SK = AKF.transpose() * X[K];
            vector_type gK = AKF.transpose() * y[K];
            for (size_t r = 0; r < cf.size(); r++) {
                g(cf[r]) -= gK(r);
                for (size_t c = 0; c < cf.size(); c++)
                    trip.emplace_back(cf[r], cf[c], -SK(r, c));
            }
        }

        S = Eigen::SparseMatrix<T>(n_f, n_f);
        S.setFromTriplets(trip.begin(), trip.end());
    }

    size_t num_clusters() const { return cluster_dofs.size(); }

    /* Full solution [cells; faces] from the face unknowns */
    vector_type
    recover(const vector_type& xF) const
    {
        vector_type x(n_c + n_f);
        x.tail(n_f) = xF;
        for (size_t K = 0; K < cluster_dofs.size(); K++) {
            const auto& cd = cluster_dofs[K];
            const auto& cf = cluster_faces[K];
            vector_type xFK(cf.size());
            for (size_t c = 0; c < cf.size(); c++)
                xFK(c) = xF(cf[c]);
            vector_type xK = y[K] - X[K] * xFK;
            for (size_t r = 0; r < cd.size(); r++)
                x(cd[r]) = xK(r);
        }
        return x;
    }
};


template<typename T>
class local_condensation {

    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    struct patch_data {
        std::vector<size_t> cell_dofs;  // global indices of the cell unknowns
        std::vector<size_t> face_dofs;  // face unknowns, numbered from 0
        matrix_type         X;          // A_KK^-1 A_KF
        vector_type         y;          // A_KK^-1 b_K
    };

    size_t                         n_c, n_f;
    std::vector<patch_data>        patches;
    std::vector<Eigen::Triplet<T>> triplets;

public:
    /* local contribution of an element: global numbering (assembled or
     * Dirichlet unknowns), Dirichlet values, local matrix and right-hand side */
    struct element {
        std::vector<assembly_index> asm_map;
        std::vector<T>              dirichlet;
        matrix_type                 lhs;
        vector_type                 rhs;
    };

    Eigen::SparseMatrix<T> S;   // Schur complement on the face unknowns
    vector_type            g;   // condensed right-hand side
    size_t                 max_patch = 0;

    local_condensation(size_t num_cell_dofs, size_t system_size)
        : n_c(num_cell_dofs), n_f(system_size - num_cell_dofs), g(vector_type::Zero(system_size - num_cell_dofs))
    {}

    /* Eliminate the cell unknowns of one patch, given the contributions of all
     * the elements involving them */
    void
    condense_patch(const std::vector<element>& elems)
    {
        /* local numbering: cell unknowns, then face unknowns */
        std::vector<size_t> kg, fg;
        for (auto& e : elems)
            for (auto& ai : e.asm_map)
                if (ai.assemble())
                    (ai.vidx() < n_c ? kg : fg).push_back(ai.vidx());
        for (auto* v : {&kg, &fg}) {
            std::sort(v->begin(), v->end());
            v->erase(std::unique(v->begin(), v->end()), v->end());
        }
        auto nk = kg.size(), nf = fg.size();
        auto pos = [&](size_t gi) -> size_t {
            if (gi < n_c)
                return std::lower_bound(kg.begin(), kg.end(), gi) - kg.begin();
            return nk + (std::lower_bound(fg.begin(), fg.end(), gi) - fg.begin());
        };

        matrix_type A = matrix_type::Zero(nk + nf, nk + nf);
        vector_type b = vector_type::Zero(nk + nf);
        for (auto& e : elems) {
            std::vector<size_t> lp(e.asm_map.size());
            for (size_t i = 0; i < e.asm_map.size(); i++)
                lp[i] = e.asm_map[i].assemble() ? pos(e.asm_map[i].vidx()) : 0;
            for (size_t i = 0; i < e.asm_map.size(); i++) {
                if (!e.asm_map[i].assemble())
                    continue;
                for (size_t j = 0; j < e.asm_map.size(); j++) {
                    if (e.asm_map[j].assemble())
                        A(lp[i], lp[j]) += e.lhs(i, j);
                    else
                        b(lp[i]) -= e.lhs(i, j) * e.dirichlet[j];
                }
                b(lp[i]) += e.rhs(i);
            }
        }

        /* elimination of the cell unknowns of the patch */
        patch_data pd;
        pd.cell_dofs = kg;
        for (auto f : fg)
            pd.face_dofs.push_back(f - n_c);
        Eigen::LDLT<matrix_type> ldlt(A.topLeftCorner(nk, nk));
        pd.X = ldlt.solve(A.topRightCorner(nk, nf));
        pd.y = ldlt.solve(b.head(nk));
        matrix_type SF = A.bottomRightCorner(nf, nf) - A.bottomLeftCorner(nf, nk) * pd.X;
        vector_type gF = b.tail(nf) - A.bottomLeftCorner(nf, nk) * pd.y;
        for (size_t r = 0; r < nf; r++) {
            g(pd.face_dofs[r]) += gF(r);
            for (size_t c = 0; c < nf; c++)
                triplets.emplace_back(pd.face_dofs[r], pd.face_dofs[c], SF(r, c));
        }
        max_patch = std::max(max_patch, nk);
        patches.push_back(std::move(pd));
    }

    void
    finalize()
    {
        S = Eigen::SparseMatrix<T>(n_f, n_f);
        S.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        triplets.shrink_to_fit();
    }

    size_t num_patches() const { return patches.size(); }

    /* Full solution [cells; faces] from the face unknowns */
    vector_type
    recover(const vector_type& xF) const
    {
        vector_type x(n_c + n_f);
        x.tail(n_f) = xF;
        for (auto& pd : patches) {
            vector_type xFK(pd.face_dofs.size());
            for (size_t c = 0; c < pd.face_dofs.size(); c++)
                xFK(c) = xF(pd.face_dofs[c]);
            vector_type xK = pd.y - pd.X * xFK;
            for (size_t r = 0; r < pd.cell_dofs.size(); r++)
                x(pd.cell_dofs[r]) = xK(r);
        }
        return x;
    }
};
