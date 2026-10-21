//
//  mf_operator.hpp
//
//  Matrix-free products of the explicit second-order operator (solver.matrix_free with formulation = "second_order",
//  no elasto-acoustic interface): with the cell unknowns u and the algebraic face unknowns a of second_order_system,
//      r = F_a - K_au u,    q = K_aa p,    z = D^-1 r,    u'' = M^-1 (F_u - K_uu u - K_ua a),
//  i.e. the products of eval_fast and of the conjugate gradients of solve_faces_fast. Each cell K keeps
//      P_K = K_au (faces of K x K),  L_K = K_aa on the faces of K,  D_K = -M_K^-1 K_uu,  E_K = -M_K^-1 K_ua,  M_K^-1,
//  where L_K holds the entries of K_aa between distinct faces of K (each belongs to the only cell having both faces,
//  checked at the construction) and the diagonal blocks K_aa(f, f) of the faces f whose first cell is K; each face block
//  keeps D^-1(f, f). As in mf_erk_operator: identical cell blocks
//  (structured meshes) are stored once (relative tolerance share_tolerance, 0: equal bit for bit, < 0: never shared),
//  their entries below share_tolerance x the largest one dropped (round-off of the assembly). The cells of a class are
//  processed in batches, their unknowns gathered (precomputed indices) as the columns of a matrix (one row per cell): a
//  small block is applied entry by entry (an AXPY over the cells of the batch, vectorized), a large and dense one by a
//  matrix-matrix product. Their contributions to the faces are written in slots laid out as the batch matrices and
//  summed per face unknown. OpenMP over the batches, then over the face unknowns or blocks: the sums do not depend on
//  the number of threads.
//

#pragma once
#ifndef second_order_mf_operator_hpp
#define second_order_mf_operator_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

template<typename T>
class mf_second_order_operator {
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dense_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;
    typedef Eigen::SparseMatrix<T>                           sparse_type;

    // K: stiffness, cell unknowns [0, nu) then algebraic face unknowns [nu, N); Minv: inverse of the cell mass (nu x nu);
    // Dinv: preconditioner of the faces (n_a x n_a, block diagonal; empty: none); cell_sizes, face_sizes: consecutive
    // diagonal blocks of the cells and of the faces
    mf_second_order_operator(const sparse_type& K, size_t nu, const sparse_type& Minv, const sparse_type& Dinv,
                             const std::vector<size_t>& cell_sizes, const std::vector<size_t>& face_sizes,
                             T share_tolerance = T(1e-10), int max_batch = 64)
        : m_nu(nu), m_na(size_t(K.rows()) - nu), m_tol(share_tolerance), m_max_batch(std::max(1, max_batch)) {
        build_layout(K, cell_sizes, face_sizes);
        build_cells(K, Minv);
        build_faces(Dinv);
        build_batches();
    }

    // r = F_a - K_au u
    void face_rhs(const T* Fa, const T* u, T* r) {
        #pragma omp parallel
        {
            dense_type U;
            #pragma omp for schedule(dynamic, 1)
            for (long b = 0; b < long(m_batches.size()); b++) {
                const batch& bt = m_batches[b];
                const cell_class& C = m_classes[bt.cls];
                gather_cells(bt, C.n, u, U);
                Eigen::Map<dense_type> S(m_slots.data() + bt.f0, bt.count, C.m);
                C.P.apply(U, S, true);
            }
            #pragma omp for schedule(static)
            for (long d = 0; d < long(m_na); d++) r[d] = Fa[d] - slot_sum(d);
        }
    }

    // q = K_aa p; with dot: returns p.q (sums over fixed chunks of faces: independent of the number of threads)
    T face_product(const T* p, T* q, bool dot = false) {
        m_part.assign(m_chunks.size() - 1, T(0));
        #pragma omp parallel
        product_team(p, q, dot ? m_part.data() : nullptr);
        T pq = 0;
        if (dot) for (T v : m_part) pq += v;
        return pq;
    }

    // z = D^-1 r
    void face_precond(const T* r, T* z) const {
        if (m_Di.empty()) throw std::runtime_error("mf_second_order_operator: no preconditioner of the faces");
        #pragma omp parallel for schedule(static)
        for (long f = 0; f < long(m_faces.size()); f++) {
            const face_item& fi = m_faces[f];
            T zz = 0, rz = 0;
            precond_face(fi.s, m_Di.data() + fi.o, r + fi.a0, z + fi.a0, zz, rz);
        }
    }

    // Conjugate gradients of solve_faces_fast on K_aa a = r preconditioned by D^-1 (same iterates up to round-off): a
    // holds the initial guess on entry, stop when |D^-1 res| <= tol |a|. One OpenMP parallel region (work-sharing loops,
    // the scalars by one thread), fused passes over the face unknowns, sums over fixed chunks of faces (independent of
    // the number of threads). Returns the number of iterations (as solve_faces_fast: 1 when the initial guess already
    // satisfies the test), 0 if not converged within max_iter.
    size_t pcg(const T* r, T* a, T tol, size_t max_iter) {
        if (m_Di.empty()) throw std::runtime_error("mf_second_order_operator: no preconditioner of the faces");
        const long nch = long(m_chunks.size()) - 1;
        m_res.resize(m_na);
        m_z.resize(m_na);
        m_p.resize(m_na);
        m_q.resize(m_na);
        T *res = m_res.data(), *z = m_z.data(), *p = m_p.data(), *q = m_q.data();
        std::vector<T> part(3 * static_cast<size_t>(nch)), ppq(static_cast<size_t>(nch));
        T zz = 0, aa = 0, rz = 0, alpha = 0, beta = 0;
        bool done = false;
        size_t result = 0;
        // res = r - q (first) or res -= alpha q, a += alpha p; z = D^-1 res; sums |z|^2, |a|^2, res.z per chunk
        auto update = [&](bool first) {
            #pragma omp for schedule(static)
            for (long c = 0; c < nch; c++) {
                const size_t d0 = m_chunks[c], d1 = m_chunks[c + 1];
                if (first)
                    for (size_t d = d0; d < d1; d++) res[d] = r[d] - q[d];
                else
                    for (size_t d = d0; d < d1; d++) { a[d] += alpha * p[d]; res[d] -= alpha * q[d]; }
                T s0 = 0, s1 = 0, s2 = 0;
                for (size_t d = d0; d < d1; d++) s1 += a[d] * a[d];
                for (int f = m_cf[c]; f < m_cf[c + 1]; f++) {
                    const face_item& fi = m_faces[f];
                    precond_face(fi.s, m_Di.data() + fi.o, res + fi.a0, z + fi.a0, s0, s2);
                }
                part[3 * c] = s0;
                part[3 * c + 1] = s1;
                part[3 * c + 2] = s2;
            }
        };
        #pragma omp parallel
        {
            product_team(a, q, nullptr);
            update(true);
            #pragma omp single
            {
                for (long c = 0; c < nch; c++) { zz += part[3 * c]; aa += part[3 * c + 1]; rz += part[3 * c + 2]; }
                done = std::sqrt(zz) <= tol * std::sqrt(aa);
                if (done) result = 1;
            }
            if (done) {
                #pragma omp for schedule(static)
                for (long d = 0; d < long(m_na); d++) a[d] += z[d];
            }
            else {
                #pragma omp for schedule(static)
                for (long d = 0; d < long(m_na); d++) p[d] = z[d];
                for (size_t m = 1; m <= max_iter; m++) {
                    product_team(p, q, ppq.data());
                    #pragma omp single
                    {
                        T pq = 0;
                        for (T v : ppq) pq += v;
                        alpha = rz / pq;
                    }
                    update(false);
                    #pragma omp single
                    {
                        T s0 = 0, s1 = 0, s2 = 0;
                        for (long c = 0; c < nch; c++) { s0 += part[3 * c]; s1 += part[3 * c + 1]; s2 += part[3 * c + 2]; }
                        zz = s0;
                        aa = s1;
                        done = std::sqrt(zz) <= tol * std::sqrt(aa);
                        if (done) result = m;
                        beta = s2 / rz;
                        rz = s2;
                    }
                    if (done) break;
                    #pragma omp for schedule(static)
                    for (long d = 0; d < long(m_na); d++) p[d] = z[d] + beta * p[d];
                }
            }
        }
        return result;
    }

    // out = M^-1 (F_u - K_uu u - K_ua a) = M^-1 F_u + D u + E a
    void cell_rates(const T* Fu, const T* u, const T* a, T* out) {
        #pragma omp parallel
        {
            dense_type U, A, F, O;
            #pragma omp for schedule(dynamic, 1)
            for (long b = 0; b < long(m_batches.size()); b++) {
                const batch& bt = m_batches[b];
                const cell_class& C = m_classes[bt.cls];
                gather_cells(bt, C.n, u, U);
                gather_cells(bt, C.n, Fu, F);
                O.resize(bt.count, C.n);
                C.D.apply(U, O, true);
                C.Mi.apply(F, O);
                if (C.m > 0) {
                    gather_faces(bt, C.m, a, A);
                    C.E.apply(A, O);
                }
                const int* idx = m_bc.data() + bt.c0;
                const T* o = O.data();
                for (Eigen::Index t = 0; t < O.size(); t++) out[idx[t]] = o[t];
            }
        }
    }

    std::string summary() const {
        size_t bytes = 0;
        size_t entries = 0, kept = 0;
        for (const cell_class& C : m_classes)
            for (const block_op* B : {&C.P, &C.L, &C.D, &C.E, &C.Mi}) {
                bytes += B->bytes();
                entries += size_t(B->rows) * B->cols;
                kept += B->nonzeros();
            }
        bytes += sizeof(T) * m_Di.size();
        return std::to_string(m_cells.size()) + " cells (" + std::to_string(m_classes.size()) + " distinct blocks), "
             + std::to_string(m_faces.size()) + " face blocks, " + std::to_string(entries ? 100 * kept / entries : 0)
             + " % of the entries kept, " + std::to_string(m_batches.size()) + " batches of cells, "
             + std::to_string(bytes >> 20) + " MB";
    }

private:
    // block A (rows x cols) applied to a batch, O (cells x rows) += X (cells x cols) A^T: by a matrix-matrix product when
    // large (at least 16 x 16) and more than half full, entry by entry otherwise (an AXPY over the cells of the batch)
    struct block_op {
        int rows = 0, cols = 0;
        bool dense = true;
        dense_type At;
        std::vector<int> ri, ci;
        std::vector<T> v;
        block_op() = default;
        block_op(const dense_type& A, T drop) : rows(int(A.rows())), cols(int(A.cols())) {
            const T mx = A.size() ? A.cwiseAbs().maxCoeff() : T(0);
            dense_type B = A;
            size_t nnz = 0;
            for (Eigen::Index j = 0; j < B.cols(); j++)
                for (Eigen::Index i = 0; i < B.rows(); i++) {
                    if (std::abs(B(i, j)) <= drop * mx) B(i, j) = T(0);
                    if (B(i, j) != T(0)) nnz++;
                }
            dense = 2 * nnz > size_t(B.size()) && rows >= 16 && cols >= 16;
            if (dense) { At = B.transpose(); return; }
            for (Eigen::Index j = 0; j < B.cols(); j++)
                for (Eigen::Index i = 0; i < B.rows(); i++)
                    if (B(i, j) != T(0)) { ri.push_back(int(i)); ci.push_back(int(j)); v.push_back(B(i, j)); }
        }
        // first: O = X A^T, else O += X A^T (O: X.rows() x rows, contiguous columns)
        template<typename Out>
        void apply(const dense_type& X, Out& O, bool first = false) const {
            if (dense) {
                if (first) O.noalias() = X * At;
                else       O.noalias() += X * At;
                return;
            }
            if (first) O.setZero();
            const Eigen::Index nb = X.rows();
            for (size_t k = 0; k < v.size(); k++) {
                T* __restrict o = O.data() + Eigen::Index(ri[k]) * nb;
                const T* __restrict x = X.data() + Eigen::Index(ci[k]) * nb;
                const T c = v[k];
                for (Eigen::Index i = 0; i < nb; i++) o[i] += c * x[i];
            }
        }
        size_t nonzeros() const { return dense ? size_t(At.size()) : v.size(); }
        size_t bytes() const { return dense ? sizeof(T) * size_t(At.size()) : v.size() * (sizeof(T) + 2 * sizeof(int)); }
    };
    struct cell_blocks { int n, m; dense_type P, L, D, E, Mi; };    // construction (comparisons)
    struct cell_class { int n, m; block_op P, L, D, E, Mi; };
    struct cell_item { int r0, n, m, cls; size_t j0; };
    struct batch { int cls, first, count; size_t f0, c0; };   // f0, c0: offsets in m_bf (and the slots), m_bc
    struct face_item { int a0, s; size_t o; };

    size_t m_nu, m_na;
    T m_tol;
    int m_max_batch;
    std::vector<cell_item> m_cells;
    std::vector<int> m_J;                   // face unknowns of the cells (algebraic numbering), cell c: [j0, j0 + m)
    std::vector<int> m_face_of;             // face block of each algebraic unknown
    std::vector<int> m_owner;               // first cell of each face block
    std::vector<cell_class> m_classes;
    std::vector<int> m_order;               // cells sorted by class
    std::vector<batch> m_batches;
    std::vector<face_item> m_faces;
    std::vector<T> m_Di;                    // diagonal blocks of D^-1 (s x s per face, row major)
    std::vector<int> m_bf, m_bc;            // batch b: face unknowns m_bf[f0 + l count + k] = J_k[l] and cell unknowns
                                            //   m_bc[c0 + j count + k] = r0_k + j of its cells k (batch matrix layout)
    std::vector<int> m_sp, m_sl;            // slots of the algebraic unknowns (CSR): slot t = entry t of m_bf
    std::vector<int> m_s1, m_s2;            // the same when every unknown has at most two slots (m_s2 = -1: one slot)
    std::vector<T> m_slots;

    // sum of the slots of the face unknown d
    T slot_sum(long d) const {
        if (!m_s1.empty()) {
            const int b = m_s2[d];
            return b < 0 ? m_slots[m_s1[d]] : m_slots[m_s1[d]] + m_slots[b];
        }
        T s = 0;
        for (int k = m_sp[d]; k < m_sp[d + 1]; k++) s += m_slots[m_sl[k]];
        return s;
    }
    std::vector<size_t> m_chunks;           // chunks of face unknowns (whole face blocks) of the sums
    std::vector<int> m_cf;                  //   and of face blocks
    std::vector<T> m_part, m_res, m_z, m_p, m_q;

    // q = K_aa p by the threads of the current team (work-sharing loops; called inside a parallel region, or alone),
    // pq[c] = p.q on the chunks c when pq is given
    void product_team(const T* p, T* q, T* pq) {
        const long nch = long(m_chunks.size()) - 1;
        dense_type A;
        #pragma omp for schedule(dynamic, 1)
        for (long b = 0; b < long(m_batches.size()); b++) {
            const batch& bt = m_batches[b];
            const cell_class& C = m_classes[bt.cls];
            gather_faces(bt, C.m, p, A);
            Eigen::Map<dense_type> S(m_slots.data() + bt.f0, bt.count, C.m);
            C.L.apply(A, S, true);
        }
        #pragma omp for schedule(static)
        for (long c = 0; c < nch; c++) {
            T s = 0;
            for (size_t d = m_chunks[c]; d < m_chunks[c + 1]; d++) {
                const T v = slot_sum(long(d));
                q[d] = v;
                s += p[d] * v;
            }
            if (pq) pq[c] = s;
        }
    }

    // z = Di r on a face block (Di row major, S x S), zz += |z|^2, rz += r.z
    template<int S>
    static void precond_face_fixed(const T* Di, const T* r, T* z, T& zz, T& rz) {
        for (int i = 0; i < S; i++) {
            T v = 0;
            for (int j = 0; j < S; j++) v += Di[i * S + j] * r[j];
            z[i] = v;
            zz += v * v;
            rz += r[i] * v;
        }
    }
    static void precond_face(int s, const T* Di, const T* r, T* z, T& zz, T& rz) {
        switch (s) {
            case 1: precond_face_fixed<1>(Di, r, z, zz, rz); return;
            case 2: precond_face_fixed<2>(Di, r, z, zz, rz); return;
            case 3: precond_face_fixed<3>(Di, r, z, zz, rz); return;
            case 4: precond_face_fixed<4>(Di, r, z, zz, rz); return;
            case 5: precond_face_fixed<5>(Di, r, z, zz, rz); return;
            case 6: precond_face_fixed<6>(Di, r, z, zz, rz); return;
            default:
                for (int i = 0; i < s; i++) {
                    T v = 0;
                    for (int j = 0; j < s; j++) v += Di[i * s + j] * r[j];
                    z[i] = v;
                    zz += v * v;
                    rz += r[i] * v;
                }
        }
    }

    // batch matrices: one row per cell of the batch, one column per unknown
    void gather_cells(const batch& bt, int n, const T* x, dense_type& X) const {
        X.resize(bt.count, n);
        const int* idx = m_bc.data() + bt.c0;
        T* o = X.data();
        for (Eigen::Index t = 0; t < X.size(); t++) o[t] = x[idx[t]];
    }
    void gather_faces(const batch& bt, int m, const T* x, dense_type& X) const {
        X.resize(bt.count, m);
        const int* idx = m_bf.data() + bt.f0;
        T* o = X.data();
        for (Eigen::Index t = 0; t < X.size(); t++) o[t] = x[idx[t]];
    }

    // dense block A(r0 .. r0+n-1, c0 .. c0+m-1)
    static dense_type dense_block(const sparse_type& A, Eigen::Index r0, Eigen::Index n, Eigen::Index c0, Eigen::Index m) {
        dense_type B = dense_type::Zero(n, m);
        for (Eigen::Index j = 0; j < m; j++)
            for (typename sparse_type::InnerIterator it(A, c0 + j); it; ++it)
                if (it.row() >= r0 && it.row() < r0 + n) B(it.row() - r0, j) = it.value();
        return B;
    }

    // cells, face blocks, face unknowns of each cell: the unknowns of the face blocks with an entry of K in its columns;
    // K_uu must be block diagonal
    void build_layout(const sparse_type& K, const std::vector<size_t>& cell_sizes, const std::vector<size_t>& face_sizes) {
        size_t r = 0;
        for (size_t s : cell_sizes) { m_cells.push_back(cell_item{int(r), int(s), 0, -1, 0}); r += s; }
        size_t a = 0, o = 0;
        m_face_of.assign(m_na, -1);
        for (size_t s : face_sizes) {
            m_faces.push_back(face_item{int(a), int(s), o});
            for (size_t i = 0; i < s; i++) m_face_of[a + i] = int(m_faces.size()) - 1;
            a += s;
            o += s * s;
        }
        if (r != m_nu || a != m_na)
            throw std::runtime_error("mf_second_order_operator: block sizes do not match the matrix");
        std::vector<std::vector<int>> J(m_cells.size());
        bool block_ok = true;
        #pragma omp parallel for schedule(dynamic, 64) reduction(&& : block_ok)
        for (long c = 0; c < long(m_cells.size()); c++) {
            const cell_item& cb = m_cells[c];
            std::vector<int> fb;
            for (int l = 0; l < cb.n; l++)
                for (typename sparse_type::InnerIterator it(K, cb.r0 + l); it; ++it) {
                    if (size_t(it.row()) >= m_nu) fb.push_back(m_face_of[it.row() - m_nu]);
                    else if (it.value() != T(0) && (it.row() < cb.r0 || it.row() >= cb.r0 + cb.n)) block_ok = false;
                }
            std::sort(fb.begin(), fb.end());
            fb.erase(std::unique(fb.begin(), fb.end()), fb.end());
            for (int f : fb)
                for (int i = 0; i < m_faces[f].s; i++) J[c].push_back(m_faces[f].a0 + i);
        }
        if (!block_ok)
            throw std::runtime_error("mf_second_order_operator: K_uu is not block diagonal");
        for (size_t c = 0; c < m_cells.size(); c++) {
            m_cells[c].j0 = m_J.size();
            m_cells[c].m = int(J[c].size());
            m_J.insert(m_J.end(), J[c].begin(), J[c].end());
        }
        // first cell of each face block (its diagonal block of K_aa goes to that cell)
        m_owner.assign(m_faces.size(), -1);
        for (size_t c = 0; c < m_cells.size(); c++)
            for (int l = 0; l < m_cells[c].m; l++) {
                int& o = m_owner[m_face_of[m_J[m_cells[c].j0 + l]]];
                if (o < 0) o = int(c);
            }
        for (int o : m_owner)
            if (o < 0) throw std::runtime_error("mf_second_order_operator: a face without cell");
    }

    // blocks of cell c; claim: number of cells using each entry of K (entries of K_aa between distinct faces)
    cell_blocks make_blocks(const sparse_type& K, const sparse_type& Minv, size_t c, std::vector<int>& claim) const {
        const cell_item& cb = m_cells[c];
        const int n = cb.n, m = cb.m;
        const int* J = m_J.data() + cb.j0;
        auto local = [&](int d) {
            const int* e = J + m;
            const int* it = std::lower_bound(J, e, d);
            return (it != e && *it == d) ? int(it - J) : -1;
        };
        dense_type Auu = dense_type::Zero(n, n), Aau = dense_type::Zero(m, n), Aua = dense_type::Zero(n, m);
        dense_type L = dense_type::Zero(m, m);
        for (int j = 0; j < n; j++)
            for (typename sparse_type::InnerIterator it(K, cb.r0 + j); it; ++it) {
                if (size_t(it.row()) < m_nu) Auu(it.row() - cb.r0, j) = it.value();
                else Aau(local(int(it.row() - m_nu)), j) = it.value();
            }
        for (int l = 0; l < m; l++)
            for (typename sparse_type::InnerIterator it(K, Eigen::Index(m_nu) + J[l]); it; ++it) {
                if (it.row() >= cb.r0 && it.row() < cb.r0 + n) Aua(it.row() - cb.r0, l) = it.value();
                else if (size_t(it.row()) >= m_nu) {
                    const int d = int(it.row() - m_nu);
                    const int l2 = local(d);
                    if (l2 < 0) continue;
                    if (m_face_of[d] == m_face_of[J[l]]) {      // diagonal block of a face: in its first cell
                        if (m_owner[m_face_of[d]] == int(c)) L(l2, l) = it.value();
                        continue;
                    }
                    L(l2, l) = it.value();
                    #pragma omp atomic
                    claim[&it.value() - K.valuePtr()]++;
                }
            }
        cell_blocks C;
        C.n = n;
        C.m = m;
        C.Mi = dense_block(Minv, cb.r0, n, cb.r0, n);
        C.D = -(C.Mi * Auu);
        C.E = -(C.Mi * Aua);
        C.P = Aau;
        C.L = L;
        return C;
    }

    static bool close_block(const dense_type& A, const dense_type& B, T tol) {
        if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
        if (A.size() == 0) return true;
        if (tol == 0) return A == B;
        const T m = A.cwiseAbs().maxCoeff();
        return (A - B).cwiseAbs().maxCoeff() <= tol * m;
    }
    static bool close(const cell_blocks& A, const cell_blocks& B, T tol) {
        return A.n == B.n && A.m == B.m && close_block(A.P, B.P, tol) && close_block(A.L, B.L, tol)
            && close_block(A.D, B.D, tol) && close_block(A.E, B.E, tol) && close_block(A.Mi, B.Mi, tol);
    }
    // entries rounded to a quantum much coarser than tol (equal blocks almost always hash alike)
    size_t hash(const cell_blocks& C) const {
        size_t h = std::hash<long long>()(C.n * 1000003LL + C.m);
        const T q = std::max(T(1e-4), 100 * m_tol);
        for (const dense_type* B : {&C.P, &C.L, &C.D, &C.E, &C.Mi}) {
            if (B->size() == 0) continue;
            T mx = B->cwiseAbs().maxCoeff();
            if (mx == 0) mx = 1;
            for (Eigen::Index i = 0; i < B->size(); i++) {
                const long long v = std::llround(B->data()[i] / mx / q);
                h ^= std::hash<long long>()(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            }
        }
        return h;
    }

    // cell blocks by chunks (in parallel), classes of identical blocks (serial; at most 16 candidates of the same hash),
    // then their operators (entries below max(share_tolerance, 0) x the largest one dropped); every entry of K_aa between
    // distinct faces must belong to exactly one cell
    void build_cells(const sparse_type& K, const sparse_type& Minv) {
        std::vector<int> claim(size_t(K.nonZeros()), 0);
        std::unordered_multimap<size_t, int> index;
        const size_t chunk = 4096;
        std::vector<cell_blocks> buf, reps;
        for (size_t c0 = 0; c0 < m_cells.size(); c0 += chunk) {
            const size_t nc = std::min(chunk, m_cells.size() - c0);
            buf.assign(nc, cell_blocks());
            std::vector<size_t> h(nc, 0);
            #pragma omp parallel for schedule(dynamic, 16)
            for (long i = 0; i < long(nc); i++) {
                buf[i] = make_blocks(K, Minv, c0 + i, claim);
                if (m_tol >= 0) h[i] = hash(buf[i]);
            }
            for (size_t i = 0; i < nc; i++) {
                int id = -1;
                if (m_tol >= 0) {
                    auto range = index.equal_range(h[i]);
                    int probes = 0;
                    for (auto it = range.first; it != range.second && probes < 16 && id < 0; ++it, ++probes)
                        if (close(reps[it->second], buf[i], m_tol)) id = it->second;
                }
                if (id < 0) {
                    reps.push_back(std::move(buf[i]));
                    id = int(reps.size()) - 1;
                    if (m_tol >= 0) index.emplace(h[i], id);
                }
                m_cells[c0 + i].cls = id;
            }
        }
        for (size_t d = 0; d < m_na; d++)
            for (typename sparse_type::InnerIterator it(K, Eigen::Index(m_nu + d)); it; ++it) {
                if (size_t(it.row()) < m_nu || m_face_of[it.row() - m_nu] == m_face_of[d] || it.value() == T(0)) continue;
                if (claim[&it.value() - K.valuePtr()] != 1)
                    throw std::runtime_error("mf_second_order_operator: an entry of K_aa between two faces does not belong "
                                             "to exactly one cell");
            }
        const T drop = std::max(T(0), m_tol);
        m_classes.resize(reps.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (long k = 0; k < long(reps.size()); k++) {
            const cell_blocks& R = reps[k];
            cell_class& C = m_classes[k];
            C.n = R.n;
            C.m = R.m;
            C.P = block_op(R.P, drop);
            C.L = block_op(R.L, drop);
            C.D = block_op(R.D, drop);
            C.E = block_op(R.E, drop);
            C.Mi = block_op(R.Mi, drop);
        }
    }

    // preconditioner blocks, rows of the face unknowns, chunks of the sums
    void build_faces(const sparse_type& Dinv) {
        size_t total = 0;
        for (const face_item& f : m_faces) total += size_t(f.s) * f.s;
        const bool precond = Dinv.rows() == Eigen::Index(m_na) && m_na > 0;
        if (precond) m_Di.assign(total, T(0));
        const size_t nch = std::max<size_t>(1, std::min<size_t>(256, m_faces.size()));
        m_chunks.assign(1, 0);
        m_cf.assign(1, 0);
        for (size_t k = 1; k <= nch; k++) {
            const size_t f = k * m_faces.size() / nch;
            m_chunks.push_back(f < m_faces.size() ? size_t(m_faces[f].a0) : m_na);
            m_cf.push_back(int(f));
        }
        if (!precond) return;
        #pragma omp parallel for schedule(static)
        for (long f = 0; f < long(m_faces.size()); f++) {
            const face_item& fi = m_faces[f];
            typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major;
            Eigen::Map<row_major>(m_Di.data() + fi.o, fi.s, fi.s) = dense_block(Dinv, fi.a0, fi.s, fi.a0, fi.s);
        }
    }

    // cells sorted by class (increasing order in a class), batches of at most m_max_batch cells of a class
    void build_batches() {
        std::vector<int> count(m_classes.size() + 1, 0);
        for (const cell_item& cb : m_cells) count[cb.cls + 1]++;
        for (size_t k = 0; k < m_classes.size(); k++) count[k + 1] += count[k];
        m_order.assign(m_cells.size(), 0);
        std::vector<int> fill(count.begin(), count.end() - 1);
        for (size_t c = 0; c < m_cells.size(); c++) m_order[fill[m_cells[c].cls]++] = int(c);
        size_t f0 = 0, c0 = 0;
        for (size_t k = 0; k < m_classes.size(); k++)
            for (int i = count[k]; i < count[k + 1]; i += m_max_batch) {
                const int nb = std::min(m_max_batch, count[k + 1] - i);
                m_batches.push_back(batch{int(k), i, nb, f0, c0});
                f0 += size_t(nb) * m_classes[k].m;
                c0 += size_t(nb) * m_classes[k].n;
            }
        // gather tables in the layout of the batch matrices (column l of the batch: unknown l of its cells)
        m_bf.assign(f0, 0);
        m_bc.assign(c0, 0);
        for (const batch& bt : m_batches) {
            const int m = m_classes[bt.cls].m, n = m_classes[bt.cls].n;
            for (int k = 0; k < bt.count; k++) {
                const cell_item& cb = m_cells[m_order[bt.first + k]];
                for (int l = 0; l < m; l++) m_bf[bt.f0 + size_t(l) * bt.count + k] = m_J[cb.j0 + l];
                for (int j = 0; j < n; j++) m_bc[bt.c0 + size_t(j) * bt.count + k] = cb.r0 + j;
            }
        }
        // slots of each algebraic unknown (CSR): the entries of m_bf holding it
        m_sp.assign(m_na + 1, 0);
        for (int d : m_bf) m_sp[d + 1]++;
        for (size_t d = 0; d < m_na; d++) m_sp[d + 1] += m_sp[d];
        m_sl.assign(m_bf.size(), 0);
        std::vector<int> pos(m_sp.begin(), m_sp.end() - 1);
        for (size_t t = 0; t < m_bf.size(); t++) m_sl[pos[m_bf[t]]++] = int(t);
        m_slots.assign(m_bf.size(), T(0));
        bool two = true;
        for (size_t d = 0; d < m_na && two; d++) two = m_sp[d + 1] - m_sp[d] >= 1 && m_sp[d + 1] - m_sp[d] <= 2;
        if (two) {
            m_s1.assign(m_na, 0);
            m_s2.assign(m_na, -1);
            for (size_t d = 0; d < m_na; d++) {
                m_s1[d] = m_sl[m_sp[d]];
                if (m_sp[d + 1] - m_sp[d] == 2) m_s2[d] = m_sl[m_sp[d] + 1];
            }
        }
    }
};

#endif /* second_order_mf_operator_hpp */
