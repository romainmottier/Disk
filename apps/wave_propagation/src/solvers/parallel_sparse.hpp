//
//  parallel_sparse.hpp
//
//  Sparse matrix tools of the fast mode (solver.fast, solver.matrix_free), built by the OpenMP threads.
//

#pragma once
#ifndef parallel_sparse_hpp
#define parallel_sparse_hpp

#include <Eigen/Sparse>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace parallel_sparse {

    // A = the sparse matrix of the triplets, duplicates summed in their order: the matrix of setFromTriplets, bit for bit,
    // built by the OpenMP threads (columns counted and filled by chunks of triplets, then sorted and summed per column)
    template<typename T>
    void set_from_triplets(Eigen::SparseMatrix<T>& A, const std::vector<Eigen::Triplet<T>>& trips) {
        typedef typename Eigen::SparseMatrix<T>::StorageIndex index_type;
        const size_t N = trips.size(), nc = size_t(A.cols());
        int nt = 1;
#ifdef _OPENMP
        nt = omp_get_max_threads();
#endif
        if (nt == 1 || N < (size_t(1) << 20)) {
            A.setFromTriplets(trips.begin(), trips.end());
            return;
        }
        // entries of each column per chunk of triplets (chunk t: threads in order)
        std::vector<std::vector<size_t>> cnt(nt, std::vector<size_t>(nc, 0));
        #pragma omp parallel num_threads(nt)
        {
            int t = 0;
#ifdef _OPENMP
            t = omp_get_thread_num();
#endif
            const size_t k0 = N * size_t(t) / size_t(nt), k1 = N * size_t(t + 1) / size_t(nt);
            for (size_t k = k0; k < k1; k++) cnt[t][size_t(trips[k].col())]++;
        }
        // positions: columns in order, the chunks in order inside a column
        std::vector<size_t> start(nc + 1, 0);
        for (size_t c = 0; c < nc; c++) {
            size_t s = 0;
            for (int t = 0; t < nt; t++) s += cnt[t][c];
            start[c + 1] = start[c] + s;
        }
        #pragma omp parallel for schedule(static)
        for (long c = 0; c < (long)nc; c++) {
            size_t p = start[c];
            for (int t = 0; t < nt; t++) { size_t n = cnt[t][c]; cnt[t][c] = p; p += n; }
        }
        std::vector<index_type> rows(N);
        std::vector<T> vals(N);
        #pragma omp parallel num_threads(nt)
        {
            int t = 0;
#ifdef _OPENMP
            t = omp_get_thread_num();
#endif
            const size_t k0 = N * size_t(t) / size_t(nt), k1 = N * size_t(t + 1) / size_t(nt);
            for (size_t k = k0; k < k1; k++) {
                const size_t p = cnt[t][size_t(trips[k].col())]++;
                rows[p] = index_type(trips[k].row());
                vals[p] = trips[k].value();
            }
        }
        std::vector<std::vector<size_t>>().swap(cnt);
        // per column: stable sort by row (order of the triplets kept), duplicates summed in that order
        std::vector<size_t> nnz(nc + 1, 0);
        #pragma omp parallel for schedule(dynamic, 256)
        for (long c = 0; c < (long)nc; c++) {
            const size_t b = start[c], e = start[c + 1];
            for (size_t i = b + 1; i < e; i++) {         // insertion sort: few entries per column, stable
                const index_type r = rows[i];
                const T v = vals[i];
                size_t j = i;
                while (j > b && rows[j - 1] > r) { rows[j] = rows[j - 1]; vals[j] = vals[j - 1]; j--; }
                rows[j] = r;
                vals[j] = v;
            }
            size_t w = b;
            for (size_t i = b; i < e; i++) {
                if (w > b && rows[w - 1] == rows[i]) vals[w - 1] += vals[i];
                else { rows[w] = rows[i]; vals[w] = vals[i]; w++; }
            }
            nnz[c + 1] = w - b;
        }
        for (size_t c = 0; c < nc; c++) nnz[c + 1] += nnz[c];
        A.setZero();
        A.makeCompressed();
        A.resizeNonZeros(Eigen::Index(nnz[nc]));
        index_type* outer = A.outerIndexPtr();
        index_type* inner = A.innerIndexPtr();
        T* values = A.valuePtr();
        #pragma omp parallel for schedule(dynamic, 256)
        for (long c = 0; c < (long)nc; c++) {
            outer[c] = index_type(nnz[c]);
            std::copy(rows.begin() + start[c], rows.begin() + start[c] + (nnz[c + 1] - nnz[c]), inner + nnz[c]);
            std::copy(vals.begin() + start[c], vals.begin() + start[c] + (nnz[c + 1] - nnz[c]), values + nnz[c]);
        }
        outer[nc] = index_type(nnz[nc]);
    }

}

#endif /* parallel_sparse_hpp */
