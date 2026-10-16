//
//  mf_erk_operator.hpp
//
//  Matrix-free explicit operator of the first-order HHO elasto-acoustic scheme (same operator as
//  erk_coupling_hho_scheme::erk_weight): with the cell unknowns y_c, the face unknowns y_f and the force F,
//      B_F(y, F) = [k_c ; k_f],   k_c = Mc^-1 (F_c - Kcc y_c - Kcf y_f),   k_f = -Sff^-1 Kfc k_c
//  where Mc, Kcc and the least-squares face block Sff are block diagonal (one block per cell, per face, or per pair of
//  elastic and acoustic faces of an interface). Each cell K keeps the blocks D_K = -Mc^-1 Kcc (n x n), E_K = -Mc^-1 Kcf
//  (n x face unknowns of K), M_K = Mc^-1 and H_K, the columns of K in the blocks G_R = -Sff^-1 Kfc of its faces R. One
//  application: per cell, k_K = D_K y_K + E_K y_f[J_K] (+ M_K F_K) and its contributions H_K k_K to its faces (written
//  in slots of the cell, while k_K is in cache); per face, k_R = sum of the contributions of its cells.
//
//  Identical blocks (structured meshes: same material, cell size and boundary configuration) are stored once: two blocks
//  are compared in the units of the energy norm (cell unknowns scaled by the square root of the diagonal of the mass,
//  face unknowns by the size of their trace map; the stresses and the velocities differ by orders of magnitude) and
//  shared when they agree to share_tolerance relative to their largest entry (0: bitwise equality, < 0: no sharing).
//  The blocks are mostly zeros (orthogonality of the bases) up to the round-off of the assembly: in these units, the
//  entries below share_tolerance x the largest entry are dropped (none when share_tolerance <= 0), and a block at most
//  half full is kept as a sparse (CSR) block, restricted to its nonzero columns. The cells with the same blocks form
//  batches of up to max_batch cells, whose unknowns are gathered as the rows of a matrix (one column per cell): each
//  entry of a sparse block is an AXPY over the cells of the batch (vectorized), a dense block a matrix-matrix product.
//  OpenMP over the batches of cells, then over the faces.
//
//  Integrators (results equal to those of erk_coupling_hho_scheme up to round-off, called by hho_first_order):
//    erk  ERK(s), stage combinations and update fused into the kernel (no intermediate vectors)
//    lts  two-level LTS-RK4 (erk_weight_LTS_coarse + erk_weight_LTS_fine)
//    li   LI-SDIRK(3,4) (erk_weight_LTS_coarse_global + dirk_weight_LTS_fine)
//  For lts and li, with P the projector on the fine unknowns, the coarse Taylor coefficients
//      w_i = B_{F_i}((I-P) a_i, (I-P) F_i),   a_0 = y, a_{i+1} = B_{F_i}(a_i, F_i)   (w_3 without force)
//  are computed as w_i = a_{i+1} - B_{F_i}(P a_i, P F_i): four applications on the whole mesh, the corrections only on
//  the region S where B P does not vanish (fine cells, cells with a fine face, faces of these cells). The fine substeps
//  (lts) and the implicit stages (li, solve restricted to the fine unknowns) only change the unknowns of S; elsewhere
//  the stages reduce to the Taylor polynomial, whose RK4 (Simpson) or SDIRK(3,4) quadrature is applied in closed form
//  (with the last product of the coarse coefficients: one pass over the unknowns per step, as a stage of ERK).
//

#pragma once
#ifndef mf_erk_operator_hpp
#define mf_erk_operator_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "runge_kutta.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

template<typename T>
class mf_erk_operator {
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dense_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;
    typedef Eigen::SparseMatrix<T>                           sparse_type;
    typedef std::function<vector_type(T)>                    force_type;
    // F(t) on a set of cell unknowns: local_force(t, dofs, out) writes F(t)[d] in out[d] for d in dofs (mlts)
    typedef std::function<void(T, const std::vector<int>&, T*)> local_force_type;
    typedef std::function<void(T, const vector_type&)>       observer_type;
    typedef std::function<bool(T, const vector_type&)>       monitor_type;      // false: stop

    // K: stiffness with the coupling terms (cell unknowns, then face unknowns), M: mass (cell block used); cells, faces:
    // layout of the unknowns, {number of cells, unknowns per cell} and {number of faces, unknowns per face} in the order
    // of the unknowns (elastic, acoustic). The inverses of the cell mass blocks and of the face blocks of K (one face, or
    // the elastic and acoustic faces of an interface: the connected blocks) are computed here (dense LU).
    mf_erk_operator(const sparse_type& K, const sparse_type& M, const std::vector<std::pair<size_t, size_t>>& cells,
                    const std::vector<std::pair<size_t, size_t>>& faces, T share_tolerance = T(1e-10), int max_batch = 64)
        : m_tol(share_tolerance), m_drop(std::max(T(0), share_tolerance)), m_max_batch(std::max(1, max_batch)) {
        m_nc = 0; m_nf = 0;
        for (auto [nb, bs] : cells) m_nc += nb * bs;
        for (auto [nb, bs] : faces) m_nf += nb * bs;
        if (m_nc + m_nf != size_t(K.rows()) || K.rows() != K.cols() || M.rows() < Eigen::Index(m_nc))
            throw std::runtime_error("matrix-free operator: layout inconsistent with the matrices");
        build_cell_structure(K, M, cells);
        build_faces(K, faces);
        build_cell_blocks(K, M);
        build_ops();
        m_all = make_set(std::vector<char>(m_cells.size(), 1), std::vector<char>(m_faces.size(), 1));
        m_side.assign(m_n_slots, T(0));
    }

    size_t n_c_dof() const { return m_nc; }
    size_t n_f_dof() const { return m_nf; }
    size_t n_dof() const { return m_nc + m_nf; }

    // Compute blocks of the cells (D, E, M, H) stored in single precision (solver.matrix_free_precision = "single"):
    // half the memory traffic of the products when the blocks are not shared (unstructured meshes), products still
    // accumulated in double; the operator is then perturbed at the 1e-7 relative level. The mass blocks (energy) and
    // the blocks of the fine operator of li stay in double.
    // Fine system of li with the face unknowns inside the fine region (all their cells fine) given by their static
    // constraint u_f = G u_c (G = -Sff^-1 Kfc, the relation the face unknowns keep from consistent initial data) instead
    // of their evolution equation u_f - gamma dt G (D u + E u) = r_f, and the stage slopes on the fine unknowns from the
    // stage solve, k_i = (U_i - R_i) / (gamma dt), R_i = U^n + dt sum_{l<i} a_il k_l (time.li_system = "constraint"):
    // the same stages in exact arithmetic; with a very stiff fine region (1/h stabilization on strongly graded meshes,
    // entries ~ h_min^-2) the evolution rows make the system ill conditioned and B P U_i multiplies the round-off of
    // U_i by its largest entries
    void set_li_constraint(bool constraint) { m_li_constraint = constraint; }
    // Coarse values in the implicit stages of li (time.li_coupling = "stage"): instead of the Taylor polynomial Y of the
    // coarse solution at t_n + c_i dt, its stage values with the SDIRK tableau, Y_i = y_n + dt sum_{l<=i} a_il Y'(c_l dt),
    // i.e. the weights dt sum_l a_il (c_l dt)^(j-1) / (j-1)! of the coefficients w_j: the coarse stage values then have
    // the stage defects of the fine ones, as in the globally implicit SDIRK scheme, and the interface does not act on the
    // stiff fine system as a source term (order reduction of the Prothero-Robinson type). The final update is unchanged:
    // SDIRK(3,4) satisfies the order conditions sum_i b_i a_il c_l = 1/6 and sum_i b_i a_il c_l^2 = 1/12, so that the
    // closed-form integral of the coarse polynomial is the same.
    void set_li_stage_coupling(bool stage) { m_li_stage_coupling = stage; }
    // Static condensation of the fine stage systems of li (time.li_condensation): the fine cell unknowns are eliminated
    // cell by cell (their block I - gamma dt D_K is local) and the Schur complement on the fine faces is factorized
    void set_li_condensation(bool condensation) { m_li_condensation = condensation; }
    // Own coefficients of the MLTS levels reusing the chain (time.mlts_reuse): w_i = B_{Q F_i}(Q a_i) equals a_{i+1} =
    // B_{F_i}(a_i) on the deep unknowns of the band (cells of the band whose faces are all band faces, faces whose cells
    // are all deep), so that the masked products are applied only on the boundary layer of the band, plus one product
    // B a_3 on the deep cells: about 4 instead of 7 products per unknown of a thick band (the coarsest level)
    void set_mlts_reuse(bool reuse) { m_mlts_reuse = reuse; }
    // Products on a set with more unknowns than this run on the threads, the others on the calling thread
    void set_parallel_threshold(size_t n) { m_par_threshold = n; }
    // Work report (output.profile): unknowns of the sets of LTS / MLTS / LI, products per macro step in units of the
    // unknowns, wall time of the parts of the integrators (LI: also the size of the fine system and of its factors)
    void set_profile(bool profile) { m_profile = profile; }

    // Solver of the fine stage systems of li, the fine unknowns in rank order (cell unknowns first): sparse LU of the
    // rows-equilibrated fine operator, or, with condensation, A = [A_CC A_CF; A_FC A_FF] with A_CC block diagonal (one
    // dense block per fine cell, inverted), S = A_FF - A_FC A_CC^-1 A_CF on the fine faces (rows equilibrated, sparse LU)
    struct li_fine_solver {
        bool condensed = false;
        int nc = 0;
        Eigen::SparseLU<sparse_type, Eigen::COLAMDOrdering<int>> lu;
        vector_type row_scale;
        sparse_type Ainv, ACF, AFC;
        size_t n_rows = 0, nnz_matrix = 0;
        static vector_type equilibrate(sparse_type& A) {
            vector_type sc = vector_type::Zero(A.rows());
            for (int k = 0; k < A.outerSize(); ++k)
                for (typename sparse_type::InnerIterator it(A, k); it; ++it)
                    sc(it.row()) = std::max(sc(it.row()), std::abs(it.value()));
            for (Eigen::Index i = 0; i < sc.size(); i++) sc(i) = sc(i) > 0 ? T(1) / sc(i) : T(1);
            A = sc.asDiagonal() * A;
            A.makeCompressed();
            return sc;
        }
        bool factor(sparse_type A, int n_cell, const std::vector<std::pair<int, int>>& cells) {
            if (!condensed) {
                row_scale = equilibrate(A);
                n_rows = size_t(A.rows()); nnz_matrix = size_t(A.nonZeros());
                lu.analyzePattern(A);
                lu.factorize(A);
                return lu.info() == Eigen::Success;
            }
            nc = n_cell;
            const int nf = int(A.rows()) - nc;
            sparse_type ACC = A.topLeftCorner(nc, nc);
            std::vector<Eigen::Triplet<T>> trips;
            for (auto [r0, n] : cells) {
                const dense_type blk = dense_type(ACC.block(r0, r0, n, n));
                const dense_type inv = Eigen::PartialPivLU<dense_type>(blk).inverse();
                for (int j = 0; j < n; j++)
                    for (int i = 0; i < n; i++)
                        if (inv(i, j) != 0) trips.emplace_back(r0 + i, r0 + j, inv(i, j));
            }
            Ainv.resize(nc, nc);
            Ainv.setFromTriplets(trips.begin(), trips.end());
            ACF = A.topRightCorner(nc, nf);
            AFC = A.bottomLeftCorner(nf, nc);
            sparse_type S = sparse_type(A.bottomRightCorner(nf, nf)) - AFC * (Ainv * ACF);
            S.prune(T(0));
            row_scale = equilibrate(S);
            n_rows = size_t(S.rows()); nnz_matrix = size_t(S.nonZeros() + Ainv.nonZeros() + ACF.nonZeros() + AFC.nonZeros());
            lu.analyzePattern(S);
            lu.factorize(S);
            return lu.info() == Eigen::Success;
        }
        vector_type solve(const vector_type& r) const {
            if (!condensed) return lu.solve(vector_type(row_scale.cwiseProduct(r)));
            const int nf = int(r.size()) - nc;
            const vector_type rc = r.head(nc);
            const vector_type g = r.tail(nf) - AFC * (Ainv * rc);
            const vector_type uf = lu.solve(vector_type(row_scale.cwiseProduct(g)));
            vector_type u(r.size());
            u.head(nc) = Ainv * vector_type(rc - ACF * uf);
            u.tail(nf) = uf;
            return u;
        }
        double nnz_lu() const { return double(lu.nnzL() + lu.nnzU()); }
    };

    void set_single_precision() {
        for (auto& o : m_cell_ops)
            for (op_block* b : {&o.A, &o.H}) {
                if (b->single) continue;
                b->valf.assign(b->val.begin(), b->val.end());
                b->densef = b->dense.template cast<float>();
                std::vector<T>().swap(b->val);
                b->dense.resize(0, 0);
                b->single = true;
            }
    }

    // sizes, for the logs
    std::string summary() const {
        size_t bytes = 0, kept = 0, full = 0, sparse = 0;
        for (auto* s : {&m_D, &m_E, &m_M, &m_G, &m_H, &m_mass})
            for (auto& b : *s) bytes += sizeof(T) * b.size();
        for (auto& o : m_cell_ops) {
            for (auto* b : {&o.A, &o.H}) bytes += sizeof(int) * b->idx.size() + b->value_bytes();
            sparse += o.A.sparse;
        }
        for (auto& b : m_mass_ops) bytes += sizeof(int) * b.idx.size() + b.value_bytes();
        bytes += sizeof(int) * (m_J.size() + m_JF.size() + m_face_cells.size()) + sizeof(size_t) * m_cell_slots.size()
               + sizeof(cell_block) * m_cells.size() + sizeof(face_block) * m_faces.size() + sizeof(T) * m_side.size();
        for (auto& c : m_cells) {
            const cell_op& o = m_cell_ops[c.op];
            kept += o.A.stored() + o.H.stored();
            full += size_t(c.n) * (c.n + c.nj) + size_t(o.H.rows) * c.n;
        }
        return std::to_string(m_cells.size()) + " cells (" + std::to_string(m_D.size()) + " D, " + std::to_string(m_E.size()) + " E, "
             + std::to_string(m_n_M) + " M, " + std::to_string(m_n_H) + " H distinct blocks), " + std::to_string(m_faces.size())
             + " face blocks (" + std::to_string(m_G.size()) + " distinct G), " + std::to_string(m_cell_ops.size()) + " classes ("
             + std::to_string(sparse) + " sparse), " + std::to_string(int(std::round(100.0 * kept / std::max<size_t>(1, full))))
             + " % of the entries kept, " + std::to_string(m_all.cells.size()) + " batches of cells, "
             + std::to_string(bytes / (1024 * 1024)) + " MB";
    }

    // out = B_F(y, F) (F: cell part of the force, nullptr for none)
    void apply(const vector_type& y, const T* Fc, vector_type& out) const {
        out.resize(n_dof());
        T* o = out.data();
        auto copy = [o](int r0, int len, const T* v) { std::copy(v, v + len, o + r0); };
        kernel(m_all, y.data(), Fc, false, copy, copy);
    }

    // Faces from the cells: x_f = -Sff^-1 Kfc x_c (erk_coupling_hho_scheme::refresh_faces_unknowns)
    void refresh_faces(vector_type& x) const {
        T* xd = x.data();
        kernel(m_all, xd, nullptr, true, [](int, int, const T*) {}, [xd](int r0, int len, const T* v) { std::copy(v, v + len, xd + r0); });
    }

    // {acoustic, elastic} parts of the discrete energy 1/2 x^T M x (cells whose unknowns precede n_elastic_dof: elastic),
    // summed in the order of the batches (the same result for any number of threads)
    std::array<T, 2> energy_parts(const vector_type& x, size_t n_elastic_dof) const {
        std::vector<T> eb(m_all.cells.size());
        const T* xd = x.data();
        eigen_serial guard;
        #pragma omp parallel
        {
            scratch sb = buffers();
            #pragma omp for schedule(dynamic, 2)
            for (long b = 0; b < (long)m_all.cells.size(); b++) {
                const batch& bt = m_all.cells[b];
                const int* items = m_all.cell_items.data() + bt.first;
                const int n = m_cells[items[0]].n, m = bt.count, MP = m;
                const op_block& A = m_mass_ops[m_mass_id[items[0]]];
                T* X = sb.in;
                for (int q = 0; q < m; q++) {
                    const T* xc = xd + m_cells[items[q]].r0;
                    for (int l = 0; l < n; l++) X[size_t(l) * MP + q] = xc[l];
                }
                op_times(A, false, X, MP, sb.c);
                T e = 0;
                for (int l = 0; l < n; l++)
                    for (int q = 0; q < m; q++) e += X[size_t(l) * MP + q] * sb.c[size_t(l) * MP + q];
                eb[b] = 0.5 * e;
            }
        }
        T ee = 0, ea = 0;
        for (size_t b = 0; b < eb.size(); b++)
            (size_t(m_cells[m_all.cell_items[m_all.cells[b].first]].r0) < n_elastic_dof ? ee : ea) += eb[b];
        return {ea, ee};
    }

    // ERK(s) (a, b, c: Butcher tableau), nt steps up to tf; zero_F: F = 0 (eval_F not called); monitor after each step,
    // stops when it returns false; returns the number of steps done
    size_t erk(vector_type& x, size_t nt, T tf, const dense_type& a, const vector_type& b, const vector_type& c,
               const force_type& eval_F, bool zero_F, const monitor_type& monitor) const {
        const size_t s = b.size(), n = n_dof();
        const T dt = tf / T(nt);
        // stage i reuses k_j (j < i - 1) when a(i, j) != 0: kept only then
        std::vector<char> keep(s, 0);
        for (size_t i = 0; i < s; i++)
            for (size_t j = 0; j + 1 < i; j++)
                if (a(i, j) != 0.0) keep[j] = 1;
        std::vector<vector_type> K(s);
        for (size_t i = 0; i < s; i++)
            if (keep[i]) K[i] = vector_type::Zero(n);
        // y: argument of the stage, written in place for the next stage (each cell reads its unknowns and its faces
        // before its own values are written; the faces are written after all the cells)
        vector_type xn = x, y = x;
        T* xd = x.data(); T* xnd = xn.data(); T* yd = y.data();
        vector_type F;
        for (size_t it = 0; it < nt; it++) {
            const T tn = dt * it;
            for (size_t i = 0; i < s; i++) {
                const T* Fc = nullptr;
                if (!zero_F) { F = eval_F(tn + c(i) * dt); Fc = F.data(); }
                const bool last = i + 1 == s;
                const T cb = dt * b(i);
                // y of the next stage: xn + dt sum_{j <= i} a(i+1, j) k_j (k_i from the kernel, the others kept)
                std::array<const T*, 8> kp{}; std::array<T, 8> kw{}; int nk = 0;
                T ci = 0;
                if (!last) {
                    ci = dt * a(i + 1, i);
                    for (size_t j = 0; j < i; j++)
                        if (a(i + 1, j) != 0.0) { kp[nk] = K[j].data(); kw[nk] = dt * a(i + 1, j); nk++; }
                }
                T* Ki = keep[i] ? K[i].data() : nullptr;
                auto post = [=](int r0, int len, const T* v) {
                    for (int q = 0; q < len; q++) {
                        const int d = r0 + q;
                        const T k = v[q];
                        xd[d] += cb * k;
                        if (Ki) Ki[d] = k;
                        if (last) { yd[d] = xd[d]; xnd[d] = xd[d]; }
                        else {
                            T yy = xnd[d] + ci * k;
                            for (int l = 0; l < nk; l++) yy += kw[l] * kp[l][d];
                            yd[d] = yy;
                        }
                    }
                };
                kernel(m_all, yd, Fc, false, post, post);
            }
            if (monitor && !monitor(dt * (it + 1), x)) return it + 1;
        }
        return nt;
    }

    // Two-level LTS-RK4, p substeps of dt / p in the fine region (fine: 1 on the fine unknowns, diagonal of Pfine)
    void lts(vector_type& x, size_t nt, T tf, size_t p, const std::vector<char>& fine, const force_type& eval_F, bool zero_F,
             const observer_type& observe) const {
        const size_t n = n_dof();
        const T dt = tf / T(nt), dtau = dt / T(p);
        local_data L = make_local(fine);
        std::array<vector_type, 4> W;
        for (auto& w : W) w = vector_type::Zero(n);
        vector_type vP = vector_type::Zero(n), acc = vector_type::Zero(n), xs = vector_type::Zero(n);
        std::array<vector_type, 3> Fco, PF;     // force coefficients F0, F1, F2 (cells), P F at the substep times
        for (auto& f : PF) f = vector_type::Zero(m_nc);
        // closed-form quadrature outside S: x += sum_j alpha_j w_j, alpha_j = sum_m dtau/6 (tm^j + 4 tmh^j + tm1^j) / j!
        std::array<T, 4> alpha{};
        for (size_t m = 0; m < p; m++) {
            T tm = m * dtau, tmh = tm + 0.5 * dtau, tm1 = tm + dtau, fact = 1;
            for (int j = 0; j < 4; j++) {
                if (j > 0) fact *= j;
                alpha[j] += dtau / 6.0 * (std::pow(tm, j) + 4.0 * std::pow(tmh, j) + std::pow(tm1, j)) / fact;
            }
        }
        const T* fmask = L.fine_mask.data();
        T* xsd = xs.data(); T* vPd = vP.data(); T* accd = acc.data();
        std::array<const T*, 4> w{W[0].data(), W[1].data(), W[2].data(), W[3].data()};
        // work report and wall times (set_profile): coarse coefficients, fine substeps
        double t_coarse = 0, t_fine = 0;
        auto now = [] { return std::chrono::steady_clock::now(); };
        auto since = [](std::chrono::steady_clock::time_point c0) {
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
        };
        if (m_profile) {
            const double N = double(n), S = double(L.S.n_dof), Wk = 4 * N + 4 * (double(p) + 1) * S;
            std::cout << "      LTS work per macro step (unknowns x products): 4 x " << N << " + 4 x " << S << " + "
                      << 4 * p << " x " << S << " = " << Wk << " = " << Wk / (4 * N) << " ERK4 steps (fine unknowns "
                      << L.fine_dofs.size() << ", region S " << S << ")" << std::endl;
        }
        for (size_t it = 0; it < nt; it++) {
            const T tn = dt * it;
            const auto c0 = m_profile ? now() : std::chrono::steady_clock::time_point();
            // fine substeps on S (xs, from x^n), the other unknowns in closed form (with the last coarse product)
            L.copy_S(x.data(), xsd);
            if (p != 1) {
                if (!zero_F) force_coefficients(eval_F, tn, dt, Fco);
                coarse_weights(x, xs, zero_F ? nullptr : &Fco, L, W, vP, PF[0], &alpha);
            }
            if (m_profile) t_coarse += since(c0);
            const auto c1 = m_profile ? now() : std::chrono::steady_clock::time_point();
            L.copy_fine(xsd, vPd);
            for (size_t m = 0; m < p; m++) {
                const T tm = m * dtau;
                const std::array<T, 3> times{tm, tm + 0.5 * dtau, tm + dtau};
                if (!zero_F)
                    for (int q = 0; q < 3; q++) L.fine_cells(eval_F(tn + times[q]).data(), PF[q].data());
                const std::array<T, 4> tau{times[0], times[1], times[1], times[2]};
                const std::array<int, 4> fi{0, 1, 1, 2};
                for (int st = 0; st < 4; st++) {
                    const std::array<double, 4> tw = runge_kutta::taylor_weights(tau[st]);
                    const T t1 = tw[1], t2 = tw[2], t3 = tw[3];
                    const T cn = st < 2 ? 0.5 * dtau : dtau;
                    auto post = [=](int r0, int len, const T* v) {
                        for (int q = 0; q < len; q++) {
                            const int d = r0 + q;
                            const T k = v[q] + w[0][d] + t1 * w[1][d] + t2 * w[2][d] + t3 * w[3][d];
                            if (st == 0)      { accd[d] = k;  vPd[d] = fmask[d] * (xsd[d] + cn * k); }
                            else if (st < 3)  { accd[d] += 2.0 * k; vPd[d] = fmask[d] * (xsd[d] + cn * k); }
                            else              { xsd[d] += dtau / 6.0 * (accd[d] + k); vPd[d] = fmask[d] * xsd[d]; }
                        }
                    };
                    kernel(L.S, vPd, zero_F ? nullptr : PF[fi[st]].data(), false, post, post, L.parallel);
                }
            }
            L.copy_S(xsd, x.data());
            if (m_profile) t_fine += since(c1);
            if (observe) observe(dt * (it + 1), x);
        }
        if (m_profile)
            std::cout << "      LTS time (s): coarse coefficients " << t_coarse << ", fine substeps " << t_fine << std::endl;
    }

    // Multilevel LTS-RK4 (local_time_stepping::integrate_mlts, Algorithm 4 of the preprint LTS_HHO) on the levels
    // level[d] = 0 (coarsest) .. L-1 (finest) of the unknowns, ratio p[l] between the steps of levels l and l+1, the own
    // Taylor coefficients of a level computed once per window of the level (refresh_each_substep: at each of its
    // sub-windows). Local: the coefficients of level l are computed on its band and the three rings of cells around it
    // (a_1, a_2, a_3 on shrinking regions, then B_{Q F_i}(Q a_i)), the unknowns of level L-1 and their neighbours (T)
    // take the RK4 steps, and every other unknown is advanced by the exact integral of its polynomial (the RK4 quadrature
    // of a cubic) at the end of each window of the finest level that reads it or changes it. Results equal to those of
    // integrate_mlts up to round-off. local_F (optional): the force on the cells read by a level (its chain) or by the
    // terminal steps only, instead of eval_F on the whole mesh at each window and each fine substep.
    void mlts(vector_type& x, size_t nt, T tf, const std::vector<int>& level, int L, const std::vector<size_t>& p,
              const force_type& eval_F, bool zero_F, const observer_type& observe, bool refresh_each_substep,
              const local_force_type& local_F = {}) const {
        // the chain reuse and the work report / timers are in a separate copy (mlts_ext): this path stays the reference
        // one (with -ffast-math, any change of its code can change the compiled floating-point operations)
        if (m_mlts_reuse || m_profile) { mlts_ext(x, nt, tf, level, L, p, eval_F, zero_F, observe, refresh_each_substep, local_F); return; }
        if (level.size() != n_dof() || L < 2 || p.size() + 1 < size_t(L))
            throw std::runtime_error("matrix-free mlts: levels of the wrong size, a single level or missing ratios");
        mlts_data D = make_mlts(level, L);
        const T dt = tf / T(nt);
        const size_t n = n_dof();
        // scratch: chain a_1..a_3, masked inputs, forces
        std::array<vector_type, 3> A;
        for (auto& a : A) a = vector_type::Zero(n);
        vector_type QA = vector_type::Zero(n), QF = vector_type::Zero(m_nc), PZ = vector_type::Zero(n);
        std::array<vector_type, 3> Fco, PFt, Fs;
        for (auto& f : PFt) f = vector_type::Zero(m_nc);
        if (local_F && !zero_F)
            for (int q = 0; q < 3; q++) { Fco[q] = vector_type::Zero(m_nc); Fs[q] = vector_type::Zero(m_nc); }
        const size_t nT = D.tdofs.size();
        std::array<vector_type, 4> Wt;
        for (auto& w : Wt) w = vector_type::Zero(nT);
        vector_type xs(nT), acc(nT), kv(nT), zs(nT);
        T* xd = x.data();

        // own coefficients of level l over [t0, t0 + h] from the current state
        auto precompute = [&](int l, T t0, T h) {
            mlts_level& V = D.lv[l];
            if (V.empty) return;
            if (!zero_F) {
                if (local_F) {
                    // F at t0, t0 + h/2, t0 + h on the cells of the chain, its Taylor coefficients (quadratic_in_time)
                    const std::array<T, 3> ts{t0, t0 + 0.5 * h, t0 + h};
                    for (int q = 0; q < 3; q++) local_F(ts[q], V.force_cells, Fs[q].data());
                    for (int d : V.force_cells) {
                        const T fn = Fs[0][d], fh = Fs[1][d], fe = Fs[2][d];
                        Fco[0][d] = fn;
                        Fco[1][d] = (-3 * fn + 4 * fh - fe) / h;
                        Fco[2][d] = (4 * fn - 8 * fh + 4 * fe) / (h * h);
                    }
                }
                else force_coefficients(eval_F, t0, h, Fco);
            }
            const T* in = xd;
            for (int i = 0; i < 3; i++) {
                T* o = A[i].data();
                auto copy = [o](int r0, int len, const T* v) { std::copy(v, v + len, o + r0); };
                kernel(V.chain[i], in, zero_F ? nullptr : Fco[i].data(), false, copy, copy, V.chain[i].n_dof > m_par_threshold);
                in = o;
            }
            for (int i = 0; i < 4; i++) {
                const T* ai = i == 0 ? xd : A[i - 1].data();
                for (int d : V.band) QA[d] = ai[d];
                const T* qf = nullptr;
                if (!zero_F && i < 3) {
                    for (int d : V.band_cells) QF[d] = Fco[i][d];
                    qf = QF.data();
                }
                T* o = V.own[i].data();
                const std::vector<int>& wd = V.wdofs;
                auto store = [o, &wd](int r0, int len, const T* v) {
                    const size_t q = std::lower_bound(wd.begin(), wd.end(), r0) - wd.begin();
                    std::copy(v, v + len, o + q);
                };
                kernel(V.wset, QA.data(), qf, false, store, store, V.wset.n_dof > m_par_threshold);
                for (int d : V.band) QA[d] = 0;
                if (qf) for (int d : V.band_cells) QF[d] = 0;
            }
            V.t_own = t0;
        };
        // exact integral over [a, b] of the polynomials of the levels k <= u on the unknowns of group u
        auto close = [&](int u, T a, T b) {
            for (int k = 0; k <= u; k++) {
                const mlts_level& V = D.lv[k];
                if (V.empty) continue;
                const T ea = a - V.t_own, eb = b - V.t_own;
                std::array<T, 4> c;
                T pa = ea, pb = eb, fact = 1;
                for (int j = 0; j < 4; j++) {
                    fact *= T(j + 1);
                    c[j] = (pb - pa) / fact;
                    pa *= ea; pb *= eb;
                }
                const T *w0 = V.own[0].data(), *w1 = V.own[1].data(), *w2 = V.own[2].data(), *w3 = V.own[3].data();
                for (auto [d, q] : D.closed[u][k])
                    xd[d] += c[0] * w0[q] + c[1] * w1[q] + c[2] * w2[q] + c[3] * w3[q];
            }
        };
        // one RK4 step of the unknowns of T over [t0, t0 + h]: k = B_{P F}(P z) + sum of the polynomials
        auto terminal = [&](T t0, T h) {
            for (auto& w : Wt) w.setZero();
            for (int k = 0; k + 1 < L; k++) {
                const mlts_level& V = D.lv[k];
                if (V.empty) continue;
                const T e = t0 - V.t_own, e2 = e * e / 2, e3 = e * e * e / 6;
                const T *w0 = V.own[0].data(), *w1 = V.own[1].data(), *w2 = V.own[2].data(), *w3 = V.own[3].data();
                for (auto [i, q] : D.term_poly[k]) {
                    Wt[0][i] += w0[q] + e * w1[q] + e2 * w2[q] + e3 * w3[q];
                    Wt[1][i] += w1[q] + e * w2[q] + e2 * w3[q];
                    Wt[2][i] += w2[q] + e * w3[q];
                    Wt[3][i] += w3[q];
                }
            }
            if (!zero_F) {
                const std::array<T, 3> ts{t0, t0 + 0.5 * h, t0 + h};
                for (int q = 0; q < 3; q++) {
                    if (local_F) {
                        local_F(ts[q], D.term_band_cells, PFt[q].data());
                        continue;
                    }
                    vector_type F = eval_F(ts[q]);
                    for (int d : D.term_band_cells) PFt[q][d] = F[d];
                }
            }
            for (size_t i = 0; i < nT; i++) xs[i] = xd[D.tdofs[i]];
            const std::array<T, 4> tau{0, 0.5 * h, 0.5 * h, h};
            const std::array<int, 4> fi{0, 1, 1, 2};
            const std::vector<int>& td = D.tdofs;
            T* kvd = kv.data();
            auto store = [kvd, &td](int r0, int len, const T* v) {
                const size_t q = std::lower_bound(td.begin(), td.end(), r0) - td.begin();
                std::copy(v, v + len, kvd + q);
            };
            for (int st = 0; st < 4; st++) {
                // stage state z on the band of level L-1
                if (st == 0) for (size_t q = 0; q < D.term_band_pos.size(); q++) PZ[D.term_band[q]] = xs[D.term_band_pos[q]];
                else {
                    const T cn = st < 3 ? 0.5 * h : h;
                    for (size_t q = 0; q < D.term_band_pos.size(); q++) {
                        const int i = D.term_band_pos[q];
                        PZ[D.term_band[q]] = xs[i] + cn * zs[i];
                    }
                }
                kernel(D.term, PZ.data(), zero_F ? nullptr : PFt[fi[st]].data(), false, store, store, D.term.n_dof > m_par_threshold);
                const std::array<double, 4> tw = runge_kutta::taylor_weights(tau[st]);
                for (size_t i = 0; i < nT; i++) {
                    const T k = kv[i] + Wt[0][i] + tw[1] * Wt[1][i] + tw[2] * Wt[2][i] + tw[3] * Wt[3][i];
                    if (st == 0) acc[i] = k;
                    else if (st < 3) acc[i] += 2 * k;
                    else acc[i] += k;
                    zs[i] = k;
                }
            }
            for (size_t i = 0; i < nT; i++) xd[D.tdofs[i]] = xs[i] + h / 6 * acc[i];
            for (int d : D.term_band) PZ[d] = 0;
        };
        std::function<void(int, T, T)> step = [&](int l, T t0, T h) {
            if (l == L - 1) { terminal(t0, h); return; }
            const T hs = h / T(p[l]);
            if (!refresh_each_substep) precompute(l, t0, h);
            for (size_t m = 0; m < p[l]; m++) {
                const T ts = t0 + T(m) * hs;
                if (refresh_each_substep) precompute(l, ts, hs);
                step(l + 1, ts, hs);
                if (refresh_each_substep) close(l, ts, ts + hs);
            }
            if (!refresh_each_substep) close(l, t0, t0 + h);
        };
        for (size_t it = 0; it < nt; it++) {
            step(0, dt * T(it), dt);
            if (observe) observe(dt * T(it + 1), x);
        }
    }

    // mlts with the options set_mlts_reuse and set_profile (same scheme)
    void mlts_ext(vector_type& x, size_t nt, T tf, const std::vector<int>& level, int L, const std::vector<size_t>& p,
              const force_type& eval_F, bool zero_F, const observer_type& observe, bool refresh_each_substep,
              const local_force_type& local_F = {}) const {
        if (level.size() != n_dof() || L < 2 || p.size() + 1 < size_t(L))
            throw std::runtime_error("matrix-free mlts: levels of the wrong size, a single level or missing ratios");
        mlts_data D = make_mlts(level, L);
        const T dt = tf / T(nt);
        const size_t n = n_dof();
        // scratch: chain a_1..a_3, masked inputs, forces
        std::array<vector_type, 3> A;
        for (auto& a : A) a = vector_type::Zero(n);
        vector_type QA = vector_type::Zero(n), QF = vector_type::Zero(m_nc), PZ = vector_type::Zero(n);
        std::array<vector_type, 3> Fco, PFt, Fs;
        for (auto& f : PFt) f = vector_type::Zero(m_nc);
        if (local_F && !zero_F)
            for (int q = 0; q < 3; q++) { Fco[q] = vector_type::Zero(m_nc); Fs[q] = vector_type::Zero(m_nc); }
        const size_t nT = D.tdofs.size();
        std::array<vector_type, 4> Wt;
        for (auto& w : Wt) w = vector_type::Zero(nT);
        vector_type xs(nT), acc(nT), kv(nT), zs(nT);
        T* xd = x.data();

        // own coefficients of level l over [t0, t0 + h] from the current state
        auto precompute = [&](int l, T t0, T h) {
            mlts_level& V = D.lv[l];
            if (V.empty) return;
            if (!zero_F) {
                if (local_F) {
                    // F at t0, t0 + h/2, t0 + h on the cells of the chain, its Taylor coefficients (quadratic_in_time)
                    const std::array<T, 3> ts{t0, t0 + 0.5 * h, t0 + h};
                    for (int q = 0; q < 3; q++) local_F(ts[q], V.force_cells, Fs[q].data());
                    for (int d : V.force_cells) {
                        const T fn = Fs[0][d], fh = Fs[1][d], fe = Fs[2][d];
                        Fco[0][d] = fn;
                        Fco[1][d] = (-3 * fn + 4 * fh - fe) / h;
                        Fco[2][d] = (4 * fn - 8 * fh + 4 * fe) / (h * h);
                    }
                }
                else force_coefficients(eval_F, t0, h, Fco);
            }
            const T* in = xd;
            for (int i = 0; i < 3; i++) {
                T* o = A[i].data();
                auto copy = [o](int r0, int len, const T* v) { std::copy(v, v + len, o + r0); };
                kernel(V.chain[i], in, zero_F ? nullptr : Fco[i].data(), false, copy, copy, V.chain[i].n_dof > m_par_threshold);
                in = o;
            }
            for (int i = 0; i < 4; i++) {
                const T* ai = i == 0 ? xd : A[i - 1].data();
                for (int d : V.band) QA[d] = ai[d];
                const T* qf = nullptr;
                if (!zero_F && i < 3) {
                    for (int d : V.band_cells) QF[d] = Fco[i][d];
                    qf = QF.data();
                }
                T* o = V.own[i].data();
                const std::vector<int>& wd = V.wdofs;
                auto store = [o, &wd](int r0, int len, const T* v) {
                    const size_t q = std::lower_bound(wd.begin(), wd.end(), r0) - wd.begin();
                    std::copy(v, v + len, o + q);
                };
                if (m_mlts_reuse) {
                    // deep unknowns: a_{i+1} from the chain, B a_3 (without force) for i = 3; then the boundary layer
                    if (i < 3) { const T* an = A[i].data(); for (auto [d, q] : V.deep_pos) o[q] = an[d]; }
                    else kernel(V.deepset, A[2].data(), nullptr, false, store, store, V.deepset.n_dof > m_par_threshold);
                    kernel(V.eset, QA.data(), qf, false, store, store, V.eset.n_dof > m_par_threshold);
                }
                else kernel(V.wset, QA.data(), qf, false, store, store, V.wset.n_dof > m_par_threshold);
                for (int d : V.band) QA[d] = 0;
                if (qf) for (int d : V.band_cells) QF[d] = 0;
            }
            V.t_own = t0;
        };
        // exact integral over [a, b] of the polynomials of the levels k <= u on the unknowns of group u
        auto close = [&](int u, T a, T b) {
            for (int k = 0; k <= u; k++) {
                const mlts_level& V = D.lv[k];
                if (V.empty) continue;
                const T ea = a - V.t_own, eb = b - V.t_own;
                std::array<T, 4> c;
                T pa = ea, pb = eb, fact = 1;
                for (int j = 0; j < 4; j++) {
                    fact *= T(j + 1);
                    c[j] = (pb - pa) / fact;
                    pa *= ea; pb *= eb;
                }
                const T *w0 = V.own[0].data(), *w1 = V.own[1].data(), *w2 = V.own[2].data(), *w3 = V.own[3].data();
                for (auto [d, q] : D.closed[u][k])
                    xd[d] += c[0] * w0[q] + c[1] * w1[q] + c[2] * w2[q] + c[3] * w3[q];
            }
        };
        // one RK4 step of the unknowns of T over [t0, t0 + h]: k = B_{P F}(P z) + sum of the polynomials
        auto terminal = [&](T t0, T h) {
            for (auto& w : Wt) w.setZero();
            for (int k = 0; k + 1 < L; k++) {
                const mlts_level& V = D.lv[k];
                if (V.empty) continue;
                const T e = t0 - V.t_own, e2 = e * e / 2, e3 = e * e * e / 6;
                const T *w0 = V.own[0].data(), *w1 = V.own[1].data(), *w2 = V.own[2].data(), *w3 = V.own[3].data();
                for (auto [i, q] : D.term_poly[k]) {
                    Wt[0][i] += w0[q] + e * w1[q] + e2 * w2[q] + e3 * w3[q];
                    Wt[1][i] += w1[q] + e * w2[q] + e2 * w3[q];
                    Wt[2][i] += w2[q] + e * w3[q];
                    Wt[3][i] += w3[q];
                }
            }
            if (!zero_F) {
                const std::array<T, 3> ts{t0, t0 + 0.5 * h, t0 + h};
                for (int q = 0; q < 3; q++) {
                    if (local_F) {
                        local_F(ts[q], D.term_band_cells, PFt[q].data());
                        continue;
                    }
                    vector_type F = eval_F(ts[q]);
                    for (int d : D.term_band_cells) PFt[q][d] = F[d];
                }
            }
            for (size_t i = 0; i < nT; i++) xs[i] = xd[D.tdofs[i]];
            const std::array<T, 4> tau{0, 0.5 * h, 0.5 * h, h};
            const std::array<int, 4> fi{0, 1, 1, 2};
            const std::vector<int>& td = D.tdofs;
            T* kvd = kv.data();
            auto store = [kvd, &td](int r0, int len, const T* v) {
                const size_t q = std::lower_bound(td.begin(), td.end(), r0) - td.begin();
                std::copy(v, v + len, kvd + q);
            };
            for (int st = 0; st < 4; st++) {
                // stage state z on the band of level L-1
                if (st == 0) for (size_t q = 0; q < D.term_band_pos.size(); q++) PZ[D.term_band[q]] = xs[D.term_band_pos[q]];
                else {
                    const T cn = st < 3 ? 0.5 * h : h;
                    for (size_t q = 0; q < D.term_band_pos.size(); q++) {
                        const int i = D.term_band_pos[q];
                        PZ[D.term_band[q]] = xs[i] + cn * zs[i];
                    }
                }
                kernel(D.term, PZ.data(), zero_F ? nullptr : PFt[fi[st]].data(), false, store, store, D.term.n_dof > m_par_threshold);
                const std::array<double, 4> tw = runge_kutta::taylor_weights(tau[st]);
                for (size_t i = 0; i < nT; i++) {
                    const T k = kv[i] + Wt[0][i] + tw[1] * Wt[1][i] + tw[2] * Wt[2][i] + tw[3] * Wt[3][i];
                    if (st == 0) acc[i] = k;
                    else if (st < 3) acc[i] += 2 * k;
                    else acc[i] += k;
                    zs[i] = k;
                }
            }
            for (size_t i = 0; i < nT; i++) xd[D.tdofs[i]] = xs[i] + h / 6 * acc[i];
            for (int d : D.term_band) PZ[d] = 0;
        };
        // wall times of the parts (set_profile): coefficients of each level, terminal steps, closed-form integrals
        std::vector<double> t_pre(L, 0.0);
        double t_term = 0, t_close = 0;
        auto timed = [this](double& acc, auto&& f) {
            if (!m_profile) { f(); return; }
            const auto c0 = std::chrono::steady_clock::now();
            f();
            acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
        };
        if (m_profile) report_mlts(D, L, p, refresh_each_substep);
        std::function<void(int, T, T)> step = [&](int l, T t0, T h) {
            if (l == L - 1) { timed(t_term, [&] { terminal(t0, h); }); return; }
            const T hs = h / T(p[l]);
            if (!refresh_each_substep) timed(t_pre[l], [&] { precompute(l, t0, h); });
            for (size_t m = 0; m < p[l]; m++) {
                const T ts = t0 + T(m) * hs;
                if (refresh_each_substep) timed(t_pre[l], [&] { precompute(l, ts, hs); });
                step(l + 1, ts, hs);
                if (refresh_each_substep) timed(t_close, [&] { close(l, ts, ts + hs); });
            }
            if (!refresh_each_substep) timed(t_close, [&] { close(l, t0, t0 + h); });
        };
        for (size_t it = 0; it < nt; it++) {
            step(0, dt * T(it), dt);
            if (observe) observe(dt * T(it + 1), x);
        }
        if (m_profile) {
            std::cout << "      MLTS time (s): coefficients";
            for (int l = 0; l + 1 < L; l++) std::cout << " " << t_pre[l];
            std::cout << ", terminal " << t_term << ", closed form " << t_close << std::endl;
        }
    }

    // LI-SDIRK(3,4) (a, b, c: SDIRK(3,4) tableau, gamma = a(0,0)): coarse Taylor polynomial, implicit stages on the fine
    // unknowns (I - gamma dt P B P, sparse LU)
    void li(vector_type& x, size_t nt, T tf, const dense_type& a, const vector_type& b, const vector_type& c,
            const std::vector<char>& fine, const force_type& eval_F, bool zero_F, const observer_type& observe) const {
        const size_t n = n_dof(), s = b.size();
        const T dt = tf / T(nt), gamma = a(0, 0);
        local_data L = make_local(fine);
        li_fine_solver solver;
        solver.condensed = m_li_condensation;
        // rows equilibrated (divided by their largest entry): on strongly graded meshes the rows of the smallest cells
        // carry entries gamma dt / h_T^2 (or 1 / h_T) larger by many orders of magnitude than those of the others, and
        // the factorization of the plain matrix loses its accuracy (the same system, solved more accurately)
        std::vector<char> crow;     // constraint rows (li_system = "constraint")
        // wall times (set_profile): assembly and factorization of the fine system, coarse coefficients, stage solves,
        // stage products and updates
        double t_fact = 0, t_coarse = 0, t_solve = 0, t_stage = 0;
        auto now = [] { return std::chrono::steady_clock::now(); };
        auto since = [](std::chrono::steady_clock::time_point c0) {
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
        };
        {
            const auto cf = now();
            if (m_li_constraint) crow = interior_fine_faces(L);
            // only the fine part is assembled: the rows and columns of the fine unknowns (cells first, then faces)
            sparse_type A = fine_operator(L, gamma * dt, m_li_constraint ? &crow : nullptr);
            std::vector<std::pair<int, int>> cells;
            for (const cell_block& cb : m_cells)
                if (L.rank[cb.r0] >= 0) cells.emplace_back(L.rank[cb.r0], cb.n);
            const int n_cell = int(L.fine_cell_dofs.size());
            if (!solver.factor(A, n_cell, cells))
                throw std::runtime_error("matrix-free li: factorization of the fine operator failed");
            t_fact = since(cf);
            // memory of the implicit part: the factorized matrix (and the condensation blocks) and the LU factors
            const double nnz_lu = solver.nnz_lu();
            std::cout << "      LI fine system: " << A.rows() << " unknowns (" << n_cell << " cell, " << A.rows() - n_cell
                      << " face), nnz " << A.nonZeros();
            if (solver.condensed) std::cout << ", condensed on the faces: " << solver.n_rows << " unknowns";
            std::cout << ", nnz(L + U) " << size_t(nnz_lu) << " (~"
                      << (nnz_lu + double(solver.nnz_matrix)) * double(sizeof(T) + sizeof(int)) / 1048576.0 << " MB)" << std::endl;
            if (m_profile) {
                const double N = double(n), S = double(L.S.n_dof), Wk = 4 * N + 4 * S + double(s) * S;
                std::cout << "      LI work per macro step (unknowns x products): 4 x " << N << " + 4 x " << S << " + "
                          << s << " x " << S << " = " << Wk << " = " << Wk / (4 * N) << " ERK4 steps, plus " << s
                          << " solves (fine unknowns " << L.fine_dofs.size() << ", region S " << S << ", solved rows "
                          << solver.n_rows << ", nnz(L + U) " << size_t(nnz_lu) << ")" << std::endl;
            }
        }
        std::array<vector_type, 4> W;
        for (auto& w : W) w = vector_type::Zero(n);
        vector_type vP = vector_type::Zero(n), zero = vector_type::Zero(n), mpf = vector_type::Zero(n);
        vector_type xs = vector_type::Zero(n), PF = vector_type::Zero(m_nc), rhs(L.fine_dofs.size()), u;
        std::vector<vector_type> K(s, vector_type::Zero(n));
        std::array<vector_type, 3> Fco;
        // closed form outside S: x += sum_j beta_j w_j, beta_j = dt sum_i b_i (c_i dt)^j / j!
        std::array<T, 4> beta{};
        for (size_t i = 0; i < s; i++) {
            T fact = 1;
            for (int j = 0; j < 4; j++) {
                if (j > 0) fact *= j;
                beta[j] += dt * b(i) * std::pow(c(i) * dt, j) / fact;
            }
        }
        std::array<const T*, 4> w{W[0].data(), W[1].data(), W[2].data(), W[3].data()};
        T* mpfd = mpf.data(); T* vPd = vP.data();
        vector_type xo = vector_type::Zero(n);  // U^n on S
        for (size_t it = 0; it < nt; it++) {
            const T tn = dt * it;
            // U^n on S, the other unknowns in closed form (with the last coarse product)
            const auto c0 = m_profile ? now() : std::chrono::steady_clock::time_point();
            L.copy_S(x.data(), xo.data());
            if (!zero_F) force_coefficients(eval_F, tn, dt, Fco);
            coarse_weights(x, xo, zero_F ? nullptr : &Fco, L, W, vP, PF, &beta);
            if (m_profile) t_coarse += since(c0);
            const auto c1 = m_profile ? now() : std::chrono::steady_clock::time_point();
            const T* xd = xo.data();
            for (size_t i = 0; i < s; i++) {
                const std::array<double, 4> tw = runge_kutta::taylor_weights(c(i) * dt);
                T t1 = tw[1], t2 = tw[2], t3 = tw[3];
                if (m_li_stage_coupling) {
                    // coarse stage values y_n + dt sum_{l<=i} a_il Y'(c_l dt) (set_li_stage_coupling)
                    t1 = t2 = t3 = 0;
                    for (size_t l = 0; l <= i; l++) {
                        const T cl = c(l) * dt;
                        t1 += dt * a(i, l);
                        t2 += dt * a(i, l) * cl;
                        t3 += dt * a(i, l) * cl * cl / 2;
                    }
                }
                if (!zero_F) {
                    L.fine_cells(eval_F(tn + c(i) * dt).data(), PF.data());
                    // Mc^-1 P F and its faces (B_F(0, P F))
                    auto post = [=](int r0, int len, const T* v) { std::copy(v, v + len, mpfd + r0); };
                    kernel(L.S, zero.data(), PF.data(), false, post, post, L.parallel);
                }
                // right-hand side on the fine unknowns: U^n + dt sum_{l<i} a_il k_l + gamma dt (Mc^-1 P F + w(tau_i))
                const long nfd = (long)L.fine_dofs.size();
                #pragma omp parallel for schedule(static) if (L.parallel)
                for (long q = 0; q < nfd; q++) {
                    const int d = L.fine_dofs[q];
                    T r = xd[d] + gamma * dt * ((zero_F ? T(0) : mpfd[d]) + w[0][d] + t1 * w[1][d] + t2 * w[2][d] + t3 * w[3][d]);
                    for (size_t l = 0; l < i; l++) r += dt * a(i, l) * K[l][d];
                    rhs[q] = crow.empty() || !crow[q] ? r : T(0);
                }
                const auto c2 = m_profile ? now() : std::chrono::steady_clock::time_point();
                u = solver.solve(rhs);
                if (m_profile) t_solve += since(c2);
                #pragma omp parallel for schedule(static) if (L.parallel)
                for (long q = 0; q < nfd; q++) vPd[L.fine_dofs[q]] = u[q];
                // k_i = B_F(P U, P F) + w(tau_i) on S
                T* Kd = K[i].data();
                auto post = [=](int r0, int len, const T* v) {
                    for (int q = 0; q < len; q++) {
                        const int d = r0 + q;
                        Kd[d] = v[q] + w[0][d] + t1 * w[1][d] + t2 * w[2][d] + t3 * w[3][d];
                    }
                };
                kernel(L.S, vPd, zero_F ? nullptr : PF.data(), false, post, post, L.parallel);
                if (m_li_constraint) {
                    #pragma omp parallel for schedule(static) if (L.parallel)
                    for (long q = 0; q < nfd; q++) {
                        const int d = L.fine_dofs[q];
                        T r = xd[d];
                        for (size_t l = 0; l < i; l++) r += dt * a(i, l) * K[l][d];
                        Kd[d] = (u[q] - r) / (gamma * dt);
                    }
                }
            }
            // U^{n+1} = U^n + dt sum_i b_i k_i on S, closed form elsewhere
            T* xsd = xs.data();
            std::vector<const T*> kp(s);
            for (size_t i = 0; i < s; i++) kp[i] = K[i].data();
            L.for_S([&](int d) {
                T v = xd[d];
                for (size_t i = 0; i < s; i++) v += dt * b(i) * kp[i][d];
                xsd[d] = v;
            });
            L.copy_S(xsd, x.data());
            if (m_profile) t_stage += since(c1);
            if (observe) observe(dt * (it + 1), x);
        }
        if (m_profile)
            std::cout << "      LI time (s): fine system " << t_fact << ", coarse coefficients " << t_coarse
                      << ", stage solves " << t_solve << ", stage products and updates " << t_stage - t_solve << std::endl;
    }

private:
    // cell K: unknowns r0 .. r0+n-1, shared blocks D[d], E[e], M[m], H[h], face unknowns J[j0 .. j0+nj-1] (global
    // numbering), slots of its contributions to its faces m_cell_slots[sl0 .. sl0+ns-1] (rows of H in this order),
    // class op
    struct cell_block { int r0, n, d, e, m, h, nj, ns, op; size_t j0, sl0; };
    // face block R: unknowns r0 .. r0+n-1 (global numbering), shared block G[g] (-Sff^-1 Kfc on the cell unknowns
    // JF[j0 .. j0+nj-1]), contributing cells m_face_cells[c0 .. c0+nc-1], the k-th one at slot s0 + k n of m_side
    struct face_block { int r0, n, g, nj, nc; size_t j0, c0, s0; };
    // Compute block, restricted to its columns with a nonzero entry (used, increasing: cell unknowns up to used_cell,
    // face unknowns up to used_noF, then the force), or to all of them (no compression); CSR on these columns (row i:
    // entries [ptr[i], ptr[i+1]), without the force [ptr[i], mid[i])) or dense (row major)
    struct op_block {
        int rows = 0, cols = 0, cols_noF = 0;
        std::vector<int> used;
        int used_cell = 0, used_noF = 0;
        bool sparse = true;
        size_t kept = 0;                                // entries kept without the force
        std::vector<int> ptr, mid, idx;
        std::vector<T> val;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dense;
        // single precision storage (set_single_precision): valf, densef instead of val, dense
        bool single = false;
        std::vector<float> valf;
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> densef;
        size_t stored() const { return sparse ? kept : size_t(rows) * used_noF; }
        size_t value_bytes() const {
            return single ? sizeof(float) * (valf.size() + densef.size()) : sizeof(T) * (val.size() + dense.size());
        }
    };
    // compute blocks of a class of cells: A = [D | E | M] (columns: cell unknowns, face unknowns, force), H (columns
    // not compressed: the rows of the cell values); slot_rows: rows of H of each face slot
    struct cell_op { op_block A, H; std::vector<int> slot_rows; };
    // batch: items[first .. first+count) of a set; contiguous: their unknowns follow each other
    struct batch { int first, count; bool contiguous; };
    // cells (in batches of the same class) and face blocks (in runs) on which the operator is applied; face_slots: bit k
    // set when the k-th contributing cell of the face block is in the set
    struct block_set {
        std::vector<int> cell_items, face_items;
        std::vector<batch> cells, faces;
        std::vector<uint32_t> face_slots;
        size_t n_dof = 0;           // unknowns of the set (cells and faces)
        int cell_chunk = 2;         // batches per OpenMP chunk (more when the batches are small: no shared blocks)
    };

    size_t m_nc, m_nf;
    T m_tol, m_drop;
    bool m_li_constraint = false;       // set_li_constraint
    bool m_li_stage_coupling = false;   // set_li_stage_coupling
    bool m_li_condensation = false;     // set_li_condensation
    bool m_mlts_reuse = false;          // set_mlts_reuse
    bool m_profile = false;             // set_profile
    size_t m_par_threshold = 20000;     // set_parallel_threshold
    int m_max_batch;
    std::vector<cell_block> m_cells;
    std::vector<face_block> m_faces;
    std::vector<int> m_J, m_JF, m_mass_id, m_cell_of_dof, m_face_cells, m_cell_slot_face;
    std::vector<size_t> m_cell_slots;
    size_t m_n_slots = 0;
    std::vector<dense_type> m_D, m_E, m_M, m_G, m_H, m_mass;
    size_t m_n_M = 0, m_n_H = 0;    // distinct M, H blocks (released after build_ops)
    vector_type m_s, m_t;           // scalings of the cell and face unknowns (energy norm)
    std::vector<cell_op> m_cell_ops;
    std::vector<op_block> m_mass_ops;       // mass blocks (energy), columns not compressed
    block_set m_all;
    mutable std::vector<T> m_side;  // contributions of the cells to their faces (one kernel at a time)
    int m_max_n = 0, m_max_cell_nj = 0, m_max_face_nj = 0, m_max_face_n = 0, m_max_H = 0;
    static constexpr int max_face_run = 256;

    // columns of the gathered matrices processed W at a time (vectorized), the remainder one by one
    static constexpr int W = 16;

    // scratch of the calling thread: gathered inputs (rows x batch), cell values, contributions to the faces (rows x
    // batch), values in the order of the unknowns
    struct scratch { T *in, *c, *h, *o; };
    scratch buffers() const {
        static thread_local std::vector<T> b;
        const size_t mp = size_t(m_max_batch);
        const size_t ni = size_t(std::max(2 * m_max_n + m_max_cell_nj, m_max_n)) * mp;
        const size_t nc = size_t(m_max_n) * mp;
        const size_t nh = size_t(std::max(1, m_max_H)) * mp;
        const size_t no = std::max(size_t(m_max_n) * m_max_batch, size_t(m_max_face_n) * max_face_run);
        if (b.size() < ni + nc + nh + no) b.resize(ni + nc + nh + no);
        return scratch{b.data(), b.data() + ni, b.data() + ni + nc, b.data() + ni + nc + nh};
    }

    // Eigen's own parallelization of the products off inside the kernels (they run in the threads of the batches)
    struct eigen_serial {
        int saved;
        eigen_serial() : saved(Eigen::nbThreads()) { Eigen::setNbThreads(1); }
        ~eigen_serial() { Eigen::setNbThreads(saved); }
    };

    // ---------------------------------------------------------------------------------------------- shared blocks
    // A block B with its scalings r, c to the units of the energy norm, S = diag(r) B diag(c): r, c from the diagonal of
    // the cell mass (cell unknowns) and from the trace map G (face unknowns), see build_cell_structure / build_faces.
    // All the entries of S are comparable (inverse times for D, E, G, H), whatever the physical units of the unknowns.
    struct scaled_block {
        dense_type B, S;
        vector_type r, c;
    };

    class block_store {
    public:
        std::vector<scaled_block> blocks;
        T tol;
        explicit block_store(T t) : tol(t) {}
        // B_A - B_B below tol x the largest entry of S_A in the units of the energy norm of A (the scaled blocks of cells
        // of different sizes may be equal, not the blocks)
        static bool close(const scaled_block& A, const scaled_block& B, T tol) {
            if (tol < 0 || A.B.rows() != B.B.rows() || A.B.cols() != B.B.cols()) return false;
            if (A.B.size() == 0) return true;
            if (tol == 0) return A.B == B.B;
            const T m = A.S.cwiseAbs().maxCoeff();
            return (A.r.asDiagonal() * (A.B - B.B) * A.c.asDiagonal()).cwiseAbs().maxCoeff() <= tol * m;
        }
        // id of a stored block equal to B (within tol), B added if none; at most max_probe blocks of the same hash are
        // compared (nearly equal blocks, e.g. the triangles of a frontal mesh, would make the search quadratic)
        int find_or_add(const scaled_block& B) {
            size_t h = 0;
            if (tol >= 0) {
                h = hash(B.S);
                auto range = m_index.equal_range(h);
                int probes = 0;
                for (auto it = range.first; it != range.second && probes < max_probe; ++it, ++probes)
                    if (close(blocks[it->second], B, tol)) return it->second;
            }
            blocks.push_back(B);
            const int id = int(blocks.size()) - 1;
            if (tol >= 0) m_index.emplace(h, id);
            return id;
        }
    private:
        static constexpr int max_probe = 16;
        std::unordered_multimap<size_t, int> m_index;
        // hash of the entries rounded to a quantum much coarser than tol (equal blocks almost always hash alike)
        size_t hash(const dense_type& B) const {
            size_t h = std::hash<long long>()(B.rows() * 1000003LL + B.cols());
            if (B.size() == 0) return h;
            T m = B.cwiseAbs().maxCoeff();
            if (m == 0) m = 1;
            const T q = std::max(T(1e-4), 100 * tol);
            for (Eigen::Index i = 0; i < B.size(); i++) {
                const long long v = std::llround(B.data()[i] / m / q);
                h ^= std::hash<long long>()(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    // Ids of the blocks of items 0..N-1 in Kinds stores: make(i, out) computes the Kinds blocks of item i (in parallel,
    // by chunks); an item takes the id of the previous one when its block is close to it (checked against the stored
    // block in parallel), the other blocks are looked up in the store (serial, few on structured meshes)
    template<int Kinds, typename Make>
    void share(size_t N, size_t item_bytes, Make&& make, std::array<block_store*, Kinds> stores,
               std::array<std::vector<int>*, Kinds> ids) {
        for (auto* id : ids) id->assign(N, -1);
        const size_t chunk = std::max<size_t>(64, std::min<size_t>(4096, (64u << 20) / std::max<size_t>(1, item_bytes)));
        std::vector<std::array<scaled_block, Kinds>> buf;
        std::array<scaled_block, Kinds> last;
        bool has_last = false;
        for (size_t c0 = 0; c0 < N; c0 += chunk) {
            const size_t nc = std::min(chunk, N - c0);
            buf.assign(nc, std::array<scaled_block, Kinds>());
            std::vector<std::array<char, Kinds>> same(nc);
            #pragma omp parallel for schedule(dynamic, 16)
            for (long i = 0; i < (long)nc; i++) make(c0 + i, buf[i]);
            #pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)nc; i++)
                for (int q = 0; q < Kinds; q++) {
                    const scaled_block* prev = i > 0 ? &buf[i - 1][q] : (has_last ? &last[q] : nullptr);
                    same[i][q] = prev && block_store::close(*prev, buf[i][q], stores[q]->tol);
                }
            for (size_t i = 0; i < nc; i++)
                for (int q = 0; q < Kinds; q++)
                    (*ids[q])[c0 + i] = (same[i][q] && c0 + i > 0) ? (*ids[q])[c0 + i - 1] : stores[q]->find_or_add(buf[i][q]);
            // a run of close blocks may drift: each item is checked against the stored block it uses
            std::vector<std::array<char, Kinds>> bad(nc);
            #pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)nc; i++)
                for (int q = 0; q < Kinds; q++)
                    bad[i][q] = same[i][q] && !block_store::close(stores[q]->blocks[(*ids[q])[c0 + i]], buf[i][q], stores[q]->tol);
            for (size_t i = 0; i < nc; i++)
                for (int q = 0; q < Kinds; q++)
                    if (bad[i][q]) (*ids[q])[c0 + i] = stores[q]->find_or_add(buf[i][q]);
            last = buf[nc - 1];
            has_last = true;
        }
    }

    // Items 0..N-1 whose key (key(i, k): the data read by the computation of their blocks, exactly) equals that of an
    // earlier item: rep[i] = the first such item (their blocks are equal bit for bit), i otherwise
    template<typename Key>
    static std::vector<size_t> exact_duplicates(size_t N, Key&& key) {
        auto hash = [](const std::vector<T>& k) {
            return std::hash<std::string_view>()(std::string_view(reinterpret_cast<const char*>(k.data()), k.size() * sizeof(T)));
        };
        std::vector<size_t> h(N), rep(N);
        #pragma omp parallel for schedule(dynamic, 64)
        for (long i = 0; i < (long)N; i++) {
            std::vector<T> k;
            key(size_t(i), k);
            h[i] = hash(k);
        }
        std::unordered_map<size_t, size_t> first;
        for (size_t i = 0; i < N; i++) rep[i] = first.emplace(h[i], i).first->second;
        #pragma omp parallel for schedule(dynamic, 64)
        for (long i = 0; i < (long)N; i++) {
            if (rep[i] == size_t(i)) continue;
            std::vector<T> a, b;
            key(size_t(i), a);
            key(rep[i], b);
            if (a != b) rep[i] = size_t(i);
        }
        return rep;
    }

    // share() on the items without an earlier exact duplicate (exact_duplicates), the others taking the ids of their
    // representative (not with m_tol < 0: one block per item); returns the representatives
    template<int Kinds, typename Make, typename Key>
    std::vector<size_t> share_unique(size_t N, size_t item_bytes, Make&& make, Key&& key, std::array<block_store*, Kinds> stores,
                                     std::array<std::vector<int>*, Kinds> ids) {
        std::vector<size_t> rep(N);
        if (m_tol < 0) {
            for (size_t i = 0; i < N; i++) rep[i] = i;
            share<Kinds>(N, item_bytes, make, stores, ids);
            return rep;
        }
        rep = exact_duplicates(N, key);
        std::vector<size_t> U, pos(N);
        for (size_t i = 0; i < N; i++)
            if (rep[i] == i) { pos[i] = U.size(); U.push_back(i); }
        std::array<std::vector<int>, Kinds> uid;
        std::array<std::vector<int>*, Kinds> uids;
        for (int q = 0; q < Kinds; q++) uids[q] = &uid[q];
        share<Kinds>(U.size(), item_bytes, [&](size_t k, std::array<scaled_block, Kinds>& out) { make(U[k], out); }, stores, uids);
        for (int q = 0; q < Kinds; q++) {
            ids[q]->assign(N, -1);
            for (size_t i = 0; i < N; i++) (*ids[q])[i] = uid[q][pos[rep[i]]];
        }
        return rep;
    }

    // ---------------------------------------------------------------------------------------------- construction
    // B and S = diag(r) B diag(c), entries of S below m_drop x its largest entry set to zero in both (round-off of the
    // assembly)
    scaled_block scale(const dense_type& B, const vector_type& r, const vector_type& c, bool drop) const {
        scaled_block b{B, r.asDiagonal() * B * c.asDiagonal(), r, c};
        if (drop && m_drop > 0 && b.S.size() > 0) {
            const T m = b.S.cwiseAbs().maxCoeff();
            for (Eigen::Index i = 0; i < b.S.size(); i++)
                if (std::abs(b.S.data()[i]) <= m_drop * m) { b.S.data()[i] = 0; b.B.data()[i] = 0; }
        }
        return b;
    }

    // entries of column j of A with rows in [r0, r1): [first, last) in the arrays of A (rows increasing in a column)
    static std::pair<Eigen::Index, Eigen::Index> col_range(const sparse_type& A, Eigen::Index j, Eigen::Index r0, Eigen::Index r1) {
        typedef typename sparse_type::StorageIndex index_type;
        const index_type* inner = A.innerIndexPtr();
        const Eigen::Index b = A.outerIndexPtr()[j], e = A.isCompressed() ? A.outerIndexPtr()[j + 1] : b + A.innerNonZeroPtr()[j];
        const Eigen::Index lo = std::lower_bound(inner + b, inner + e, index_type(r0)) - inner;
        const Eigen::Index hi = std::lower_bound(inner + lo, inner + e, index_type(r1)) - inner;
        return {lo, hi};
    }

    // dense block A(r0 .. r0+n-1, c0 .. c0+m-1)
    static dense_type dense_block(const sparse_type& A, Eigen::Index r0, Eigen::Index n, Eigen::Index c0, Eigen::Index m) {
        dense_type B = dense_type::Zero(n, m);
        for (Eigen::Index j = 0; j < m; j++) {
            auto [lo, hi] = col_range(A, c0 + j, r0, r0 + n);
            for (Eigen::Index p = lo; p < hi; p++) B(A.innerIndexPtr()[p] - r0, j) = A.valuePtr()[p];
        }
        return B;
    }

    // cells, their face unknowns J (rows of the nonzero entries of Kfc in their columns, global numbering; the assembly
    // stores the zeros of the local matrices); scaling of the cell unknowns: square root of the diagonal of the mass; Kcc
    // must be block diagonal
    void build_cell_structure(const sparse_type& K, const sparse_type& M, const std::vector<std::pair<size_t, size_t>>& layout) {
        size_t r = 0;
        for (auto [nb, cbs] : layout) {
            for (size_t i = 0; i < nb; i++, r += cbs)
                m_cells.push_back(cell_block{int(r), int(cbs), -1, -1, -1, -1, 0, 0, -1, 0, 0});
            m_max_n = std::max(m_max_n, int(cbs));
        }
        m_cell_of_dof.assign(m_nc, -1);
        for (size_t c = 0; c < m_cells.size(); c++)
            for (int l = 0; l < m_cells[c].n; l++) m_cell_of_dof[m_cells[c].r0 + l] = int(c);
        m_s = vector_type::Ones(m_nc);
        std::vector<std::vector<int>> J(m_cells.size());
        bool block_ok = true;
        #pragma omp parallel for schedule(dynamic, 64) reduction(&& : block_ok)
        for (long c = 0; c < (long)m_cells.size(); c++) {
            const cell_block& cb = m_cells[c];
            for (int l = 0; l < cb.n; l++) {
                const Eigen::Index j = cb.r0 + l;
                auto [lo, hi] = col_range(M, j, j, j + 1);
                if (hi > lo && M.valuePtr()[lo] > 0) m_s[j] = std::sqrt(M.valuePtr()[lo]);
                auto [c0, c1] = col_range(K, j, 0, m_nc);
                for (Eigen::Index p = c0; p < c1; p++)
                    if (K.valuePtr()[p] != T(0) && (K.innerIndexPtr()[p] < cb.r0 || K.innerIndexPtr()[p] >= cb.r0 + cb.n)) block_ok = false;
                auto [f0, f1] = col_range(K, j, m_nc, m_nc + m_nf);
                for (Eigen::Index p = f0; p < f1; p++)
                    if (K.valuePtr()[p] != T(0)) J[c].push_back(int(K.innerIndexPtr()[p]));
            }
            std::sort(J[c].begin(), J[c].end());
            J[c].erase(std::unique(J[c].begin(), J[c].end()), J[c].end());
        }
        if (!block_ok)
            throw std::runtime_error("matrix-free operator: Kcc is not block diagonal");
        for (size_t c = 0; c < m_cells.size(); c++) {
            m_cells[c].j0 = m_J.size();
            m_cells[c].nj = int(J[c].size());
            m_J.insert(m_J.end(), J[c].begin(), J[c].end());
            m_max_cell_nj = std::max(m_max_cell_nj, int(J[c].size()));
        }
    }

    // Face blocks (one per face, from the layout), grouped by the nonzero entries of the face block of K (an interface:
    // elastic and acoustic faces); G = -Sff^-1 Kfc on the rows of each face block, Sff^-1 inverse of the face block of K of its
    // group, cell unknowns JF: unknowns of the cells around the group with an entry in its rows of Kfc. Scaling t of the
    // face unknowns: a face value t_f x_f = sum_j (t_f G_fj / s_j) (s_j x_j) is of the size of the scaled cell values
    // (rows of t G diag(1/s) of largest entry 1). Contributing cells of each face block and their slots.
    void build_faces(const sparse_type& K, const std::vector<std::pair<size_t, size_t>>& layout) {
        size_t r = m_nc;
        for (auto [nb, fbs] : layout)
            for (size_t i = 0; i < nb; i++, r += fbs) {
                m_faces.push_back(face_block{int(r), int(fbs), -1, 0, 0, 0, 0, 0});
                m_max_face_n = std::max(m_max_face_n, int(fbs));
            }
        const int nfb = int(m_faces.size());
        std::vector<int> block_of(m_nf);
        for (int b = 0; b < nfb; b++)
            for (int q = 0; q < m_faces[b].n; q++) block_of[m_faces[b].r0 - m_nc + q] = b;
        // groups: connected face blocks (union-find)
        std::vector<int> parent(nfb);
        for (int b = 0; b < nfb; b++) parent[b] = b;
        auto find = [&](int b) { while (parent[b] != b) b = parent[b] = parent[parent[b]]; return b; };
        for (Eigen::Index j = m_nc; j < Eigen::Index(m_nc + m_nf); j++) {
            auto [lo, hi] = col_range(K, j, m_nc, m_nc + m_nf);
            for (Eigen::Index p = lo; p < hi; p++) {
                if (K.valuePtr()[p] == T(0)) continue;
                const int a = find(block_of[j - m_nc]), b = find(block_of[K.innerIndexPtr()[p] - m_nc]);
                if (a != b) parent[std::max(a, b)] = std::min(a, b);
            }
        }
        std::vector<int> group(nfb), first_of_group;
        std::vector<std::vector<int>> members;
        for (int b = 0; b < nfb; b++) {
            const int root = find(b);
            if (root == b) { group[b] = int(members.size()); members.emplace_back(); }
            else group[b] = group[root];
            members[group[b]].push_back(b);
        }
        const int ng = int(members.size());
        // cells around each group (rows of their columns in the group), then the unknowns of these cells with an entry
        std::vector<std::vector<int>> group_cells(ng);
        {
            std::vector<std::vector<int>> cell_groups(m_cells.size());
            #pragma omp parallel for schedule(dynamic, 64)
            for (long c = 0; c < (long)m_cells.size(); c++) {
                const cell_block& cb = m_cells[c];
                for (int j = 0; j < cb.nj; j++) cell_groups[c].push_back(group[block_of[m_J[cb.j0 + j] - m_nc]]);
                std::sort(cell_groups[c].begin(), cell_groups[c].end());
                cell_groups[c].erase(std::unique(cell_groups[c].begin(), cell_groups[c].end()), cell_groups[c].end());
            }
            for (size_t c = 0; c < m_cells.size(); c++)
                for (int g : cell_groups[c]) group_cells[g].push_back(int(c));
        }
        std::vector<std::vector<int>> JF(ng);
        #pragma omp parallel for schedule(dynamic, 64)
        for (long g = 0; g < ng; g++)
            for (int c : group_cells[g])
                for (int l = 0; l < m_cells[c].n; l++) {
                    const int j = m_cells[c].r0 + l;
                    bool any = false;
                    for (int b : members[g]) {
                        auto [lo, hi] = col_range(K, j, m_faces[b].r0, m_faces[b].r0 + m_faces[b].n);
                        for (Eigen::Index p = lo; p < hi && !any; p++) any = K.valuePtr()[p] != T(0);
                    }
                    if (any) JF[g].push_back(j);
                }
        for (int b = 0; b < nfb; b++) {
            face_block& f = m_faces[b];
            const std::vector<int>& jf = JF[group[b]];
            f.j0 = m_JF.size();
            f.nj = int(jf.size());
            m_JF.insert(m_JF.end(), jf.begin(), jf.end());
            m_max_face_nj = std::max(m_max_face_nj, f.nj);
            f.c0 = m_face_cells.size();
            for (int j : jf)
                if (f.nc == 0 || m_face_cells.back() != m_cell_of_dof[j]) { m_face_cells.push_back(m_cell_of_dof[j]); f.nc++; }
            if (f.nc > 32)
                throw std::runtime_error("matrix-free operator: more than 32 cells around a face block");
            f.s0 = m_n_slots;
            m_n_slots += size_t(f.n) * f.nc;
        }
        // slots of each cell, in the order of the face blocks
        std::vector<int> count(m_cells.size(), 0);
        for (auto& f : m_faces)
            for (int k = 0; k < f.nc; k++) count[m_face_cells[f.c0 + k]]++;
        size_t pos = 0;
        for (size_t c = 0; c < m_cells.size(); c++) { m_cells[c].sl0 = pos; m_cells[c].ns = 0; pos += count[c]; }
        m_cell_slots.assign(pos, 0);
        m_cell_slot_face.assign(pos, -1);
        for (int b = 0; b < nfb; b++) {
            const face_block& f = m_faces[b];
            for (int k = 0; k < f.nc; k++) {
                cell_block& cb = m_cells[m_face_cells[f.c0 + k]];
                m_cell_slots[cb.sl0 + cb.ns] = f.s0 + size_t(k) * f.n;
                m_cell_slot_face[cb.sl0 + cb.ns] = b;
                cb.ns++;
            }
        }
        m_t = vector_type::Ones(m_nf);
        // inputs of the block of face b: face block Sff of K of its group (row0: its rows), columns J of the face
        auto face_inputs = [&](size_t b, dense_type& Sff, dense_type& Kb, int& row0) {
            const face_block& f = m_faces[b];
            const std::vector<int>& mem = members[group[b]];
            std::vector<int> S;
            row0 = 0;
            for (int m : mem) {
                if (m == int(b)) row0 = int(S.size());
                for (int q = 0; q < m_faces[m].n; q++) S.push_back(m_faces[m].r0 + q);
            }
            const int ns = int(S.size());
            Sff = dense_type::Zero(ns, ns);
            Kb = dense_type::Zero(ns, f.nj);
            int c0 = 0;
            for (int mc : mem) {
                int r0 = 0;
                for (int mr : mem) {
                    Sff.block(r0, c0, m_faces[mr].n, m_faces[mc].n) = dense_block(K, m_faces[mr].r0, m_faces[mr].n, m_faces[mc].r0, m_faces[mc].n);
                    r0 += m_faces[mr].n;
                }
                c0 += m_faces[mc].n;
            }
            const int* Jf = m_JF.data() + f.j0;
            for (int j = 0; j < f.nj; j++) {
                int r0 = 0;
                for (int mr : mem) {
                    auto [lo, hi] = col_range(K, Jf[j], m_faces[mr].r0, m_faces[mr].r0 + m_faces[mr].n);
                    for (Eigen::Index p = lo; p < hi; p++) Kb(r0 + K.innerIndexPtr()[p] - m_faces[mr].r0, j) = K.valuePtr()[p];
                    r0 += m_faces[mr].n;
                }
            }
        };
        auto make = [&](size_t b, std::array<scaled_block, 1>& out) {
            const face_block& f = m_faces[b];
            dense_type Sff, Kb;
            int row0 = 0;
            face_inputs(b, Sff, Kb, row0);
            const int* Jf = m_JF.data() + f.j0;
            const dense_type Sinv = Sff.partialPivLu().inverse();
            const dense_type G = -(Sinv.middleRows(row0, f.n) * Kb);
            vector_type is(f.nj), t(f.n);
            for (int j = 0; j < f.nj; j++) is[j] = 1 / m_s[Jf[j]];
            for (int q = 0; q < f.n; q++) {
                const T mx = f.nj ? (G.row(q).transpose().cwiseProduct(is)).cwiseAbs().maxCoeff() : T(0);
                t[q] = mx > 0 ? 1 / mx : T(1);
                m_t[f.r0 - m_nc + q] = t[q];
            }
            out[0] = scale(G, t, is, true);
        };
        auto key = [&](size_t b, std::vector<T>& k) {
            const face_block& f = m_faces[b];
            dense_type Sff, Kb;
            int row0 = 0;
            face_inputs(b, Sff, Kb, row0);
            k.assign({T(f.n), T(f.nj), T(row0), T(Sff.rows())});
            k.insert(k.end(), Sff.data(), Sff.data() + Sff.size());
            k.insert(k.end(), Kb.data(), Kb.data() + Kb.size());
            for (int j = 0; j < f.nj; j++) k.push_back(m_s[m_JF[f.j0 + j]]);
        };
        std::vector<int> g;
        block_store sG(m_tol);
        const size_t bytes = sizeof(T) * size_t(2 * m_max_face_n) * size_t(m_max_face_nj + 2 * m_max_face_n);
        std::vector<size_t> rep = share_unique<1>(m_faces.size(), bytes, make, key, {&sG}, {&g});
        for (int b = 0; b < nfb; b++) {
            m_faces[b].g = g[b];
            if (rep[b] != size_t(b))     // scaling of an exact duplicate: that of its representative
                m_t.segment(m_faces[b].r0 - m_nc, m_faces[b].n) = m_t.segment(m_faces[rep[b]].r0 - m_nc, m_faces[rep[b]].n);
        }
        for (auto& sb : sG.blocks) m_G.push_back(sb.B);
    }

    // blocks of the cells, scaled for the comparisons: D = -Mc^-1 Kcc by s, 1/s; E = -Mc^-1 Kcf by s, 1/t; M = Mc^-1 by
    // s, s; the mass by 1/s, 1/s; H (rows of G of its faces on the columns of the cell) by t, 1/s. Kcf must vanish out of
    // the face unknowns J of the cells (checked).
    void build_cell_blocks(const sparse_type& K, const sparse_type& M) {
        for (auto& c : m_cells) {
            int rows = 0;
            for (int k = 0; k < c.ns; k++) rows += m_faces[m_cell_slot_face[c.sl0 + k]].n;
            m_max_H = std::max(m_max_H, rows);
        }
        // nonzero entries of Kcf: all, read through J
        long kcf_all = 0, kcf_read = 0;
        #pragma omp parallel for schedule(static) reduction(+ : kcf_all)
        for (long j = long(m_nc); j < long(m_nc + m_nf); j++) {
            auto [lo, hi] = col_range(K, j, 0, m_nc);
            for (Eigen::Index p = lo; p < hi; p++) kcf_all += K.valuePtr()[p] != T(0);
        }
        // inputs of the blocks of cell ci: mass, cell block and face columns J of K; nonzero entries of these columns read
        auto cell_inputs = [&](size_t ci, dense_type& Mc, dense_type& Kc, dense_type& Kf) {
            const cell_block& c = m_cells[ci];
            Mc = dense_block(M, c.r0, c.n, c.r0, c.n);
            Kc = dense_block(K, c.r0, c.n, c.r0, c.n);
            Kf = dense_type::Zero(c.n, c.nj);
            const int* Jk = m_J.data() + c.j0;
            long read = 0;
            for (int j = 0; j < c.nj; j++) {
                auto [lo, hi] = col_range(K, Jk[j], c.r0, c.r0 + c.n);
                for (Eigen::Index p = lo; p < hi; p++) {
                    Kf(K.innerIndexPtr()[p] - c.r0, j) = K.valuePtr()[p];
                    read += K.valuePtr()[p] != T(0);
                }
            }
            return read;
        };
        #pragma omp parallel for schedule(dynamic, 64) reduction(+ : kcf_read)
        for (long ci = 0; ci < (long)m_cells.size(); ci++) {
            dense_type Mc, Kc, Kf;
            kcf_read += cell_inputs(size_t(ci), Mc, Kc, Kf);
        }
        auto make = [&](size_t ci, std::array<scaled_block, 5>& out) {
            const cell_block& c = m_cells[ci];
            dense_type Mc, Kc, Kf;
            cell_inputs(ci, Mc, Kc, Kf);
            const dense_type Minv = Mc.partialPivLu().inverse();
            const int* Jk = m_J.data() + c.j0;
            const vector_type s = m_s.segment(c.r0, c.n);
            vector_type tf(c.nj);
            for (int j = 0; j < c.nj; j++) tf[j] = m_t[Jk[j] - m_nc];
            out[0] = scale(-(Minv * Kc), s, s.cwiseInverse(), true);
            out[1] = scale(-(Minv * Kf), s, tf.cwiseInverse(), true);
            out[2] = scale(Minv, s, s, true);
            out[3] = scale(Mc, s.cwiseInverse(), s.cwiseInverse(), false);
            // H: rows of G of each face block of the cell (in the order of its slots), columns of the cell
            int rows = 0;
            for (int k = 0; k < c.ns; k++) rows += m_faces[m_cell_slot_face[c.sl0 + k]].n;
            dense_type H = dense_type::Zero(rows, c.n);
            vector_type th(rows);
            int row = 0;
            for (int k = 0; k < c.ns; k++) {
                const face_block& f = m_faces[m_cell_slot_face[c.sl0 + k]];
                const dense_type& G = m_G[f.g];
                for (int j = 0; j < f.nj; j++) {
                    const int d = m_JF[f.j0 + j];
                    if (d >= c.r0 && d < c.r0 + c.n) H.block(row, d - c.r0, f.n, 1) = G.col(j);
                }
                th.segment(row, f.n) = m_t.segment(f.r0 - m_nc, f.n);
                row += f.n;
            }
            out[4] = scale(H, th, s.cwiseInverse(), true);
        };
        // the data read by make: blocks of K and M, scalings, slots (G blocks and columns of the faces of the cell)
        auto key = [&](size_t ci, std::vector<T>& k) {
            const cell_block& c = m_cells[ci];
            dense_type Mc, Kc, Kf;
            cell_inputs(ci, Mc, Kc, Kf);
            k.assign({T(c.n), T(c.nj), T(c.ns)});
            k.insert(k.end(), Mc.data(), Mc.data() + Mc.size());
            k.insert(k.end(), Kc.data(), Kc.data() + Kc.size());
            k.insert(k.end(), Kf.data(), Kf.data() + Kf.size());
            for (int l = 0; l < c.n; l++) k.push_back(m_s[c.r0 + l]);
            for (int j = 0; j < c.nj; j++) k.push_back(m_t[m_J[c.j0 + j] - m_nc]);
            for (int q = 0; q < c.ns; q++) {
                const face_block& f = m_faces[m_cell_slot_face[c.sl0 + q]];
                k.push_back(T(f.g));
                k.push_back(T(f.n));
                k.push_back(T(f.nj));
                for (int j = 0; j < f.nj; j++) {
                    const int d = m_JF[f.j0 + j];
                    k.push_back(d >= c.r0 && d < c.r0 + c.n ? T(d - c.r0) : T(-1));
                }
                for (int l = 0; l < f.n; l++) k.push_back(m_t[f.r0 - m_nc + l]);
            }
        };
        std::vector<int> d, e, m, h;
        block_store sD(m_tol), sE(m_tol), sM(m_tol), sMass(m_tol), sH(m_tol);
        const size_t bytes = sizeof(T) * size_t(m_max_n) * (4 * m_max_n + m_max_cell_nj + m_max_H);
        share_unique<5>(m_cells.size(), bytes, make, key, {&sD, &sE, &sM, &sMass, &sH}, {&d, &e, &m, &m_mass_id, &h});
        if (kcf_read != kcf_all)
            throw std::runtime_error("matrix-free operator: Kcf has entries out of the faces of the cells");
        for (size_t c = 0; c < m_cells.size(); c++) {
            m_cells[c].d = d[c]; m_cells[c].e = e[c]; m_cells[c].m = m[c]; m_cells[c].h = h[c];
        }
        for (auto& b : sD.blocks) m_D.push_back(b.B);
        for (auto& b : sE.blocks) m_E.push_back(b.B);
        for (auto& b : sM.blocks) m_M.push_back(b.B);
        for (auto& b : sMass.blocks) m_mass.push_back(b.B);
        for (auto& b : sH.blocks) m_H.push_back(b.B);
    }

    // Compute block of the given parts (side by side; the last one is the force when force), restricted to its nonzero
    // columns when compress; CSR when at most half of the entries without the force are nonzero
    static op_block make_op(const std::vector<const dense_type*>& parts, bool force, int n_cell_cols, bool compress) {
        op_block A;
        A.rows = int(parts[0]->rows());
        for (auto* p : parts) A.cols += int(p->cols());
        A.cols_noF = force ? A.cols - int(parts.back()->cols()) : A.cols;
        dense_type full = dense_type::Zero(A.rows, A.cols);
        int c0 = 0;
        for (auto* p : parts) {
            full.middleCols(c0, p->cols()) = *p;
            c0 += int(p->cols());
        }
        std::vector<int> pos(A.cols, -1);
        for (int j = 0; j < A.cols; j++)
            if (!compress || (full.col(j).array() != T(0)).any()) {
                pos[j] = int(A.used.size());
                A.used.push_back(j);
                if (j < n_cell_cols) A.used_cell++;
                if (j < A.cols_noF) A.used_noF++;
            }
        A.kept = (full.leftCols(A.cols_noF).array() != T(0)).count();
        A.sparse = 2 * A.kept <= size_t(A.rows) * A.used_noF;
        if (!A.sparse) {
            A.dense.resize(A.rows, A.used.size());
            for (size_t u = 0; u < A.used.size(); u++) A.dense.col(u) = full.col(A.used[u]);
            return A;
        }
        A.ptr.push_back(0);
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                if (j == A.cols_noF) A.mid.push_back(int(A.idx.size()));
                if (full(i, j) != T(0)) { A.idx.push_back(pos[j]); A.val.push_back(full(i, j)); }
            }
            if (A.cols_noF == A.cols) A.mid.push_back(int(A.idx.size()));
            A.ptr.push_back(int(A.idx.size()));
        }
        return A;
    }

    // one class per set of blocks (D, E, M, H) and layout of the slots; the mass blocks. The dense M, H and mass blocks
    // are then released (only in the compute blocks; D, E, G stay for the fine operator of li)
    void build_ops() {
        for (auto& B : m_mass) m_mass_ops.push_back(make_op({&B}, false, int(B.cols()), false));
        std::map<std::vector<int>, int> cls;
        for (auto& c : m_cells) {
            std::vector<int> key{c.d, c.e, c.m, c.h, c.n, c.nj};
            for (int k = 0; k < c.ns; k++) key.push_back(m_faces[m_cell_slot_face[c.sl0 + k]].n);
            auto it = cls.emplace(key, int(m_cell_ops.size())).first;
            if (it->second == int(m_cell_ops.size())) {
                cell_op o;
                o.A = make_op({&m_D[c.d], &m_E[c.e], &m_M[c.m]}, true, c.n, true);
                o.H = make_op({&m_H[c.h]}, false, c.n, false);
                o.slot_rows.assign(key.begin() + 6, key.end());
                m_cell_ops.push_back(std::move(o));
            }
            c.op = it->second;
        }
        m_n_M = m_M.size(); m_n_H = m_H.size();
        std::vector<dense_type>().swap(m_M);
        std::vector<dense_type>().swap(m_H);
        std::vector<dense_type>().swap(m_mass);
    }

    // Cells in batches (same class, increasing order in each class; a batch stays contiguous while its cells follow
    // each other: cut at a gap once a quarter full, gathered otherwise), face blocks in runs of consecutive blocks
    block_set make_set(const std::vector<char>& cell_in, const std::vector<char>& face_in) const {
        block_set set;
        std::vector<std::vector<int>> members(m_cell_ops.size());
        for (size_t K = 0; K < m_cells.size(); K++)
            if (cell_in[K]) { members[m_cells[K].op].push_back(int(K)); set.n_dof += m_cells[K].n; }
        for (auto& cl : members) {
            for (size_t q = 0; q < cl.size(); q++) {
                const int K = cl[q];
                if (q > 0) {
                    batch& b = set.cells.back();
                    const bool next = set.cell_items.back() + 1 == K;
                    if (b.count < m_max_batch && (!b.contiguous || next || b.count < m_max_batch / 4)) {
                        if (b.contiguous && !next) b.contiguous = false;
                        set.cell_items.push_back(K);
                        b.count++;
                        continue;
                    }
                }
                set.cells.push_back(batch{int(set.cell_items.size()), 1, true});
                set.cell_items.push_back(K);
            }
        }
        // about 256 cells per chunk
        set.cell_chunk = std::max(2, int(256 * set.cells.size() / std::max<size_t>(1, set.cell_items.size())));
        for (size_t R = 0; R < m_faces.size(); R++) {
            if (!face_in[R]) continue;
            const face_block& f = m_faces[R];
            set.n_dof += f.n;
            uint32_t mask = 0;
            for (int k = 0; k < f.nc; k++)
                if (cell_in[m_face_cells[f.c0 + k]]) mask |= uint32_t(1) << k;
            if (!set.faces.empty()) {
                batch& b = set.faces.back();
                if (b.contiguous && b.count < max_face_run && set.face_items.back() + 1 == int(R)) {
                    set.face_items.push_back(int(R));
                    set.face_slots.push_back(mask);
                    b.count++;
                    continue;
                }
            }
            set.faces.push_back(batch{int(set.face_items.size()), 1, true});
            set.face_items.push_back(int(R));
            set.face_slots.push_back(mask);
        }
        return set;
    }

    // ---------------------------------------------------------------------------------------------- kernels
    // C (rows x MP) = A X, X (used columns x MP, used_noF without the force), C and X row major
    static void op_times(const op_block& A, bool force, const T* X, int MP, T* C) {
        if (A.single) {
            op_times_single(A, force, X, MP, C);
            return;
        }
        if (!A.sparse) {
            typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rm_type;
            const int nc = force ? int(A.used.size()) : A.used_noF;
            if (MP == 1) {      // one cell (cells without shared blocks): row dot products, no product dispatch
                const T* a = A.dense.data();
                const Eigen::Index ld = A.dense.cols();
                for (int i = 0; i < A.rows; i++) {
                    T acc = 0;
                    const T* ai = a + i * ld;
                    #pragma omp simd reduction(+ : acc)
                    for (int j = 0; j < nc; j++) acc += ai[j] * X[j];
                    C[i] = acc;
                }
                return;
            }
            Eigen::Map<rm_type>(C, A.rows, MP).noalias() = A.dense.leftCols(nc) * Eigen::Map<const rm_type>(X, nc, MP);
            return;
        }
        for (int i = 0; i < A.rows; i++) {
            const int b = A.ptr[i], e = force ? A.ptr[i + 1] : A.mid[i];
            T* ci = C + size_t(i) * MP;
            int q0 = 0;
            for (; q0 + W <= MP; q0 += W) {
                T acc[W] = {};
                for (int p = b; p < e; p++) {
                    const T a = A.val[p];
                    const T* x = X + size_t(A.idx[p]) * MP + q0;
                    #pragma omp simd
                    for (int k = 0; k < W; k++) acc[k] += a * x[k];
                }
                for (int k = 0; k < W; k++) ci[q0 + k] = acc[k];
            }
            if (q0 < MP) {
                const int rem = MP - q0;
                T acc[W] = {};
                for (int p = b; p < e; p++) {
                    const T a = A.val[p];
                    const T* x = X + size_t(A.idx[p]) * MP + q0;
                    for (int k = 0; k < rem; k++) acc[k] += a * x[k];
                }
                for (int k = 0; k < rem; k++) ci[q0 + k] = acc[k];
            }
        }
    }

    // op_times with the single precision blocks (products accumulated in double)
    static void op_times_single(const op_block& A, bool force, const T* X, int MP, T* C) {
        if (!A.sparse) {
            const int nc = force ? int(A.used.size()) : A.used_noF;
            const float* a = A.densef.data();
            const Eigen::Index ld = A.densef.cols();
            if (MP == 1) {      // one cell: row dot products
                for (int i = 0; i < A.rows; i++) {
                    T acc = 0;
                    const float* ai = a + i * ld;
                    #pragma omp simd reduction(+ : acc)
                    for (int j = 0; j < nc; j++) acc += T(ai[j]) * X[j];
                    C[i] = acc;
                }
                return;
            }
            for (int i = 0; i < A.rows; i++) {
                const float* ai = a + i * ld;
                T* ci = C + size_t(i) * MP;
                for (int q = 0; q < MP; q++) ci[q] = 0;
                for (int j = 0; j < nc; j++) {
                    const T aij = T(ai[j]);
                    const T* x = X + size_t(j) * MP;
                    #pragma omp simd
                    for (int q = 0; q < MP; q++) ci[q] += aij * x[q];
                }
            }
            return;
        }
        for (int i = 0; i < A.rows; i++) {
            const int b = A.ptr[i], e = force ? A.ptr[i + 1] : A.mid[i];
            T* ci = C + size_t(i) * MP;
            int q0 = 0;
            for (; q0 + W <= MP; q0 += W) {
                T acc[W] = {};
                for (int p = b; p < e; p++) {
                    const T a = T(A.valf[p]);
                    const T* x = X + size_t(A.idx[p]) * MP + q0;
                    #pragma omp simd
                    for (int k = 0; k < W; k++) acc[k] += a * x[k];
                }
                for (int k = 0; k < W; k++) ci[q0 + k] = acc[k];
            }
            if (q0 < MP) {
                const int rem = MP - q0;
                T acc[W] = {};
                for (int p = b; p < e; p++) {
                    const T a = T(A.valf[p]);
                    const T* x = X + size_t(A.idx[p]) * MP + q0;
                    for (int k = 0; k < rem; k++) acc[k] += a * x[k];
                }
                for (int k = 0; k < rem; k++) ci[q0 + k] = acc[k];
            }
        }
    }

    // Cells of the set: k_K = D_K y_K + E_K y_f[J_K] (+ M_K F_K) (identity: k_K = y_K), post(r0, len, values) on each run
    // of cells (not for identity), contributions H_K k_K to the slots of the faces
    template<typename Post>
    void cell_loop(const block_set& set, const T* y, const T* Fc, bool identity, scratch sb, Post& post) const {
        const bool force = Fc != nullptr && !identity;
        T* side = m_side.data();
        #pragma omp for schedule(dynamic, set.cell_chunk)
        for (long b = 0; b < (long)set.cells.size(); b++) {
            const batch& bt = set.cells[b];
            const int* items = set.cell_items.data() + bt.first;
            const cell_block& c0 = m_cells[items[0]];
            const cell_op& op = m_cell_ops[c0.op];
            const int n = c0.n, m = bt.count, nj = c0.nj, MP = m;
            T* C = sb.c;
            if (identity) {
                for (int q = 0; q < m; q++) {
                    const T* yc = y + m_cells[items[q]].r0;
                    for (int l = 0; l < n; l++) C[size_t(l) * MP + q] = yc[l];
                }
            }
            else {
                // rows of X: the used unknowns of the cells, of their faces, force; one column per cell
                const op_block& A = op.A;
                T* X = sb.in;
                const int* U = A.used.data();
                const int nu = force ? int(A.used.size()) : A.used_noF;
                for (int q = 0; q < m; q++) {
                    const cell_block& c = m_cells[items[q]];
                    const T* yc = y + c.r0;
                    const int* J = m_J.data() + c.j0;
                    int r = 0;
                    for (; r < A.used_cell; r++) X[size_t(r) * MP + q] = yc[U[r]];
                    for (; r < A.used_noF; r++) X[size_t(r) * MP + q] = y[J[U[r] - n]];
                    if (force)
                        for (const T* f = Fc + c.r0 - n - nj; r < nu; r++) X[size_t(r) * MP + q] = f[U[r]];
                }
                op_times(A, force, X, MP, C);
            }
            // contributions to the faces
            if (op.H.rows > 0) {
                op_times(op.H, false, C, MP, sb.h);
                const int ns = int(op.slot_rows.size());
                for (int q = 0; q < m; q++) {
                    const size_t* sl = m_cell_slots.data() + m_cells[items[q]].sl0;
                    for (int k = 0, row = 0; k < ns; row += op.slot_rows[k], k++) {
                        T* dst = side + sl[k];
                        for (int l = 0; l < op.slot_rows[k]; l++) dst[l] = sb.h[size_t(row + l) * MP + q];
                    }
                }
            }
            if (identity) continue;
            // values in the order of the unknowns
            for (int q = 0; q < m; q++)
                for (int l = 0; l < n; l++) sb.o[size_t(q) * n + l] = C[size_t(l) * MP + q];
            if (bt.contiguous)
                post(c0.r0, n * m, sb.o);
            else
                for (int q = 0; q < m; q++) post(m_cells[items[q]].r0, n, sb.o + size_t(q) * n);
        }
    }

    // Face blocks of the set: sum of the contributions of their cells in the set, post(r0, len, values) on each run
    template<typename Post>
    void face_loop(const block_set& set, scratch sb, Post&& post) const {
        const T* side = m_side.data();
        #pragma omp for schedule(dynamic, 4)
        for (long b = 0; b < (long)set.faces.size(); b++) {
            const batch& bt = set.faces[b];
            const int* items = set.face_items.data() + bt.first;
            const uint32_t* masks = set.face_slots.data() + bt.first;
            size_t off = 0;
            for (int q = 0; q < bt.count; q++) {
                const face_block& f = m_faces[items[q]];
                T* o = sb.o + off;
                std::fill(o, o + f.n, T(0));
                for (int k = 0; k < f.nc; k++)
                    if (masks[q] >> k & 1u) {
                        const T* s = side + f.s0 + size_t(k) * f.n;
                        for (int l = 0; l < f.n; l++) o[l] += s[l];
                    }
                off += f.n;
            }
            post(m_faces[items[0]].r0, int(off), sb.o);
        }
    }

    // B_F(y, F) on a set (identity: x_f from x_c), post_c / post_f on each run of cell / face values
    template<typename PostC, typename PostF>
    void kernel(const block_set& set, const T* y, const T* Fc, bool identity, PostC&& post_c, PostF&& post_f,
                bool parallel = true) const {
        eigen_serial guard;
        #pragma omp parallel if (parallel)
        {
            scratch sb = buffers();
            cell_loop(set, y, Fc, identity, sb, post_c);
            face_loop(set, sb, post_f);
        }
    }

    // ---------------------------------------------------------------------------------------------- MLTS
    // Level l < L-1: band (unknowns of the level, its cell unknowns), kernels of a_1, a_2, a_3 (cells X_1 ⊃ X_2 ⊃ X_3, the
    // band and up to three rings of neighbours through the faces) and of w_i = B_{Q F_i}(Q a_i) (cells Z: band cells and
    // cells with a face of the band, faces of Z), support wdofs of w (increasing), own coefficients on wdofs
    // set_mlts_reuse: eset (boundary layer: cells of Z not deep and deep cells with a face not deep, faces of FZ not
    // deep), deepset (deep cells and faces), deep_pos (deep unknown, position in wdofs)
    struct mlts_level {
        bool empty = true;
        std::vector<int> band, band_cells, wdofs, force_cells;     // force_cells: cell unknowns of X_1
        block_set chain[3], wset, eset, deepset;
        std::vector<std::pair<int, int>> deep_pos;
        std::array<vector_type, 4> own;
        T t_own = 0;
    };
    // T: support of B P (level L-1 band, cells with a face of the band, their faces); closed[u][k]: (unknown, position in
    // the support of level k) of the unknowns of group u (finest level that reads or changes them: u) for k <= u;
    // term_poly[k]: (position in tdofs, position in the support of level k)
    struct mlts_data {
        std::vector<mlts_level> lv;
        block_set term;
        std::vector<int> tdofs, term_band, term_band_pos, term_band_cells;
        std::vector<std::vector<std::vector<std::pair<int, int>>>> closed;
        std::vector<std::vector<std::pair<int, int>>> term_poly;
    };

    mlts_data make_mlts(const std::vector<int>& level, int L) const {
        const size_t ncb = m_cells.size(), nfb = m_faces.size();
        std::vector<int> cl(ncb), fl(nfb);
        for (size_t K = 0; K < ncb; K++) cl[K] = level[m_cells[K].r0];
        for (size_t R = 0; R < nfb; R++) fl[R] = level[m_faces[R].r0];
        auto faces_of = [&](size_t K, auto&& f) { const cell_block& c = m_cells[K]; for (int k = 0; k < c.ns; k++) f(size_t(m_cell_slot_face[c.sl0 + k])); };
        auto cells_of = [&](size_t R, auto&& f) { const face_block& fb = m_faces[R]; for (int k = 0; k < fb.nc; k++) f(size_t(m_face_cells[fb.c0 + k])); };
        // faces of a set of cells; cells of a set of faces
        auto faces_of_cells = [&](const std::vector<char>& cin) {
            std::vector<char> fin(nfb, 0);
            for (size_t K = 0; K < ncb; K++) if (cin[K]) faces_of(K, [&](size_t R) { fin[R] = 1; });
            return fin;
        };
        auto cells_of_faces = [&](const std::vector<char>& fin, std::vector<char> cin) {
            for (size_t R = 0; R < nfb; R++) if (fin[R]) cells_of(R, [&](size_t K) { cin[K] = 1; });
            return cin;
        };
        auto dofs_of = [&](const std::vector<char>& cin, const std::vector<char>& fin) {
            std::vector<int> d;
            for (size_t K = 0; K < ncb; K++) if (cin[K]) for (int l = 0; l < m_cells[K].n; l++) d.push_back(m_cells[K].r0 + l);
            for (size_t R = 0; R < nfb; R++) if (fin[R]) for (int l = 0; l < m_faces[R].n; l++) d.push_back(m_faces[R].r0 + l);
            std::sort(d.begin(), d.end());
            return d;
        };
        mlts_data D;
        D.lv.resize(L - 1);
        std::vector<int> u(n_dof(), -1);
        for (int l = 0; l + 1 < L; l++) {
            mlts_level& V = D.lv[l];
            std::vector<char> C(ncb, 0), Fb(nfb, 0);
            for (size_t K = 0; K < ncb; K++) C[K] = cl[K] == l;
            for (size_t R = 0; R < nfb; R++) Fb[R] = fl[R] == l;
            V.band = dofs_of(C, Fb);
            if (V.band.empty()) continue;
            V.empty = false;
            for (int d : V.band) if (size_t(d) < m_nc) V.band_cells.push_back(d);
            // chain regions: a_3 on (X_3, band faces), a_2 on (X_2, faces(X_3)), a_1 on (X_1, faces(X_2))
            std::vector<char> X3 = cells_of_faces(Fb, C);
            std::vector<char> F3 = faces_of_cells(X3), X2 = cells_of_faces(F3, X3);
            std::vector<char> F2 = faces_of_cells(X2), X1 = cells_of_faces(F2, X2);
            std::vector<char> F1 = faces_of_cells(X1);
            V.chain[0] = make_set(X1, F2);
            V.force_cells = dofs_of(X1, std::vector<char>(nfb, 0));
            V.chain[1] = make_set(X2, F3);
            V.chain[2] = make_set(X3, Fb);
            // w: cells Z (band cells and cells with a band face), faces of Z
            std::vector<char> Z = cells_of_faces(Fb, C), FZ = faces_of_cells(Z);
            V.wset = make_set(Z, FZ);
            V.wdofs = dofs_of(Z, FZ);
            if (m_mlts_reuse) {
                std::vector<char> deep(ncb, 0), fdeep(nfb, 0);
                for (size_t K = 0; K < ncb; K++)
                    if (C[K]) { bool all = true; faces_of(K, [&](size_t R) { all = all && Fb[R]; }); deep[K] = all; }
                for (size_t R = 0; R < nfb; R++)
                    if (Fb[R]) { bool all = true; cells_of(R, [&](size_t K) { all = all && deep[K]; }); fdeep[R] = all; }
                std::vector<char> E(ncb, 0), FE(nfb, 0);
                for (size_t R = 0; R < nfb; R++) FE[R] = FZ[R] && !fdeep[R];
                for (size_t K = 0; K < ncb; K++) E[K] = Z[K] && !deep[K];
                for (size_t R = 0; R < nfb; R++) if (FE[R]) cells_of(R, [&](size_t K) { if (deep[K]) E[K] = 1; });
                V.eset = make_set(E, FE);
                V.deepset = make_set(deep, fdeep);
                for (int d : dofs_of(deep, fdeep))
                    V.deep_pos.push_back({d, int(std::lower_bound(V.wdofs.begin(), V.wdofs.end(), d) - V.wdofs.begin())});
            }
            for (auto& w : V.own) w = vector_type::Zero(V.wdofs.size());
            for (int d : V.wdofs) u[d] = std::max(u[d], l);
            for (int d : dofs_of(X1, F1)) u[d] = std::max(u[d], l);     // read by the chain
        }
        // terminal support T
        {
            std::vector<char> C(ncb, 0), Fb(nfb, 0);
            for (size_t K = 0; K < ncb; K++) C[K] = cl[K] == L - 1;
            for (size_t R = 0; R < nfb; R++) Fb[R] = fl[R] == L - 1;
            D.term_band = dofs_of(C, Fb);
            for (int d : D.term_band) if (size_t(d) < m_nc) D.term_band_cells.push_back(d);
            std::vector<char> Z = cells_of_faces(Fb, C), FZ = faces_of_cells(Z);
            D.term = make_set(Z, FZ);
            D.tdofs = dofs_of(Z, FZ);
            for (int d : D.term_band) D.term_band_pos.push_back(int(std::lower_bound(D.tdofs.begin(), D.tdofs.end(), d) - D.tdofs.begin()));
            for (int d : D.tdofs) u[d] = L - 1;
        }
        for (size_t d = 0; d < n_dof(); d++)
            if (u[d] < 0) throw std::runtime_error("matrix-free mlts: an unknown outside every level");
        D.closed.assign(L - 1, std::vector<std::vector<std::pair<int, int>>>(L - 1));
        D.term_poly.assign(L - 1, {});
        for (int k = 0; k + 1 < L; k++) {
            const mlts_level& V = D.lv[k];
            for (size_t q = 0; q < V.wdofs.size(); q++) {
                const int d = V.wdofs[q];
                if (u[d] == L - 1)
                    D.term_poly[k].push_back({int(std::lower_bound(D.tdofs.begin(), D.tdofs.end(), d) - D.tdofs.begin()), int(q)});
                else
                    D.closed[u[d]][k].push_back({d, int(q)});
            }
        }
        return D;
    }

    // Products of one macro step of mlts in units of the unknowns (a product on a set costs its unknowns): per level
    // its chain and own coefficients times its windows per macro step, the terminal steps; and in ERK4 steps on the mesh
    void report_mlts(const mlts_data& D, int L, const std::vector<size_t>& p, bool refresh_each_substep) const {
        double P = 1, W = 0;
        std::cout << "      MLTS work per macro step (unknowns x products x windows):";
        for (int l = 0; l + 1 < L; l++) {
            const mlts_level& V = D.lv[l];
            const double mult = refresh_each_substep ? P * double(p[l]) : P;
            double w = 0;
            if (!V.empty) {
                w = double(V.chain[0].n_dof + V.chain[1].n_dof + V.chain[2].n_dof);
                w += m_mlts_reuse ? 4.0 * double(V.eset.n_dof) + double(V.deepset.n_dof) : 4.0 * double(V.wset.n_dof);
            }
            std::cout << " [" << l << ": band " << V.band.size() << ", " << w << " x " << mult << "]";
            W += w * mult;
            P *= double(p[l]);
        }
        const double wt = 4.0 * double(D.term.n_dof) * P;
        W += wt;
        std::cout << " [terminal: band " << D.term_band.size() << ", " << 4.0 * double(D.term.n_dof) << " x " << P << "]"
                  << std::endl << "      MLTS work: " << W << " per macro step = " << W / (4.0 * double(n_dof()))
                  << " ERK4 steps (" << n_dof() << " unknowns), terminal " << 100.0 * wt / W << " %" << std::endl;
    }

    // ---------------------------------------------------------------------------------------------- LTS / LI
    // Region S where B P does not vanish, fine unknowns
    struct local_data {
        block_set S;
        std::vector<std::pair<int, int>> S_ranges;      // unknowns of S, [begin, end)
        std::vector<int> fine_dofs;                     // fine unknowns, increasing
        std::vector<int> fine_cell_dofs;                // fine cell unknowns
        std::vector<int> rank;                          // rank of an unknown in fine_dofs, -1 if not fine
        vector_type fine_mask;                          // 1 on the fine unknowns, 0 elsewhere
        bool parallel = true;                           // S large enough for the threads

        template<typename F>
        void for_S(F&& f) const {
            #pragma omp parallel for schedule(dynamic, 1) if (parallel)
            for (long r = 0; r < (long)S_ranges.size(); r++)
                for (int d = S_ranges[r].first; d < S_ranges[r].second; d++) f(d);
        }
        void copy_S(const T* from, T* to) const { for_S([=](int d) { to[d] = from[d]; }); }
        void copy_fine(const T* from, T* to) const {
            const long nf = (long)fine_dofs.size();
            #pragma omp parallel for schedule(static) if (parallel)
            for (long q = 0; q < nf; q++) to[fine_dofs[q]] = from[fine_dofs[q]];
        }
        // P F on the cells (zero elsewhere, kept)
        void fine_cells(const T* F, T* PF) const { for (int d : fine_cell_dofs) PF[d] = F[d]; }
    };

    local_data make_local(const std::vector<char>& fine) const {
        if (fine.size() != n_dof())
            throw std::runtime_error("matrix-free operator: fine projector of the wrong size");
        local_data L;
        L.rank.assign(n_dof(), -1);
        L.fine_mask = vector_type::Zero(n_dof());
        for (size_t d = 0; d < n_dof(); d++)
            if (fine[d]) {
                L.rank[d] = int(L.fine_dofs.size());
                L.fine_dofs.push_back(int(d));
                L.fine_mask[d] = 1;
                if (d < m_nc) L.fine_cell_dofs.push_back(int(d));
            }
        // S: fine cells, cells with a fine face unknown; face blocks with a fine unknown or a cell of S
        std::vector<char> cell_in(m_cells.size(), 0), face_in(m_faces.size(), 0);
        for (size_t K = 0; K < m_cells.size(); K++) {
            const cell_block& c = m_cells[K];
            bool in = false;
            for (int l = 0; l < c.n && !in; l++) in = fine[c.r0 + l];
            for (int j = 0; j < c.nj && !in; j++) in = fine[m_J[c.j0 + j]];
            cell_in[K] = in;
        }
        for (size_t R = 0; R < m_faces.size(); R++) {
            const face_block& f = m_faces[R];
            bool in = false;
            for (int q = 0; q < f.n && !in; q++) in = fine[f.r0 + q];
            for (int k = 0; k < f.nc && !in; k++) in = cell_in[m_face_cells[f.c0 + k]];
            face_in[R] = in;
        }
        L.S = make_set(cell_in, face_in);
        for (size_t K = 0; K < m_cells.size(); K++)
            if (cell_in[K]) add_range(L.S_ranges, m_cells[K].r0, m_cells[K].r0 + m_cells[K].n);
        for (size_t R = 0; R < m_faces.size(); R++)
            if (face_in[R]) add_range(L.S_ranges, m_faces[R].r0, m_faces[R].r0 + m_faces[R].n);
        // ranges of moderate length for the threads
        std::vector<std::pair<int, int>> split;
        for (auto [b, e] : L.S_ranges)
            for (int s = b; s < e; s += 4096) split.push_back({s, std::min(e, s + 4096)});
        L.S_ranges.swap(split);
        L.parallel = L.S.n_dof > m_par_threshold;
        return L;
    }

    static void add_range(std::vector<std::pair<int, int>>& r, int b, int e) {
        if (!r.empty() && r.back().second == b) r.back().second = e;
        else r.push_back({b, e});
    }

    // F0 = F(tn), F1 = (-3 F(tn) + 4 F(tn + dt/2) - F(tn + dt)) / dt, F2 = 4 (F(tn) - 2 F(tn + dt/2) + F(tn + dt)) / dt^2 (cells)
    void force_coefficients(const force_type& eval_F, T tn, T dt, std::array<vector_type, 3>& Fco) const {
        vector_type Fn = eval_F(tn).head(m_nc), Fh = eval_F(tn + 0.5 * dt).head(m_nc), F1 = eval_F(tn + dt).head(m_nc);
        Fco[0] = Fn;
        runge_kutta::quadratic_in_time(Fn, Fh, F1, dt, Fco[1], Fco[2]);
    }

    // Coarse Taylor coefficients w_0..w_3 (W): a_1 = B_{F0}(y, F0), a_2 = B_{F1}(a_1), a_3 = B_{F2}(a_2), B a_3 on the whole
    // mesh, then w_i -= B_{F_i}(P a_i, P F_i) on S (descending i: P a_i read before a_i becomes w_{i-1}; P y read in
    // y_S, a copy of y on S). closed: y += sum_j closed_j w_j with the last product (right off S; S is then reset)
    void coarse_weights(vector_type& y, const vector_type& y_S, const std::array<vector_type, 3>* Fco, const local_data& L,
                        std::array<vector_type, 4>& W, vector_type& vP, vector_type& PF, const std::array<T, 4>* closed) const {
        for (int i = 0; i < 4; i++) {
            const T* in = i == 0 ? y.data() : W[i - 1].data();
            T* o = W[i].data();
            const T* Fc = (Fco && i < 3) ? (*Fco)[i].data() : nullptr;
            if (i == 3 && closed) {
                T* yd = y.data();
                const T *w0 = W[0].data(), *w1 = W[1].data(), *w2 = W[2].data();
                const std::array<T, 4> cf = *closed;
                auto post = [=](int r0, int len, const T* v) {
                    for (int q = 0; q < len; q++) {
                        const int d = r0 + q;
                        o[d] = v[q];
                        yd[d] += cf[0] * w0[d] + cf[1] * w1[d] + cf[2] * w2[d] + cf[3] * v[q];
                    }
                };
                kernel(m_all, in, Fc, false, post, post);
            }
            else {
                auto copy = [o](int r0, int len, const T* v) { std::copy(v, v + len, o + r0); };
                kernel(m_all, in, Fc, false, copy, copy);
            }
        }
        T* vPd = vP.data();
        for (int i = 3; i >= 0; i--) {
            L.copy_fine(i == 0 ? y_S.data() : W[i - 1].data(), vPd);
            const T* Fc = nullptr;
            if (Fco && i < 3) { L.fine_cells((*Fco)[i].data(), PF.data()); Fc = PF.data(); }
            T* o = W[i].data();
            auto post = [o](int r0, int len, const T* v) { for (int q = 0; q < len; q++) o[r0 + q] -= v[q]; };
            kernel(L.S, vPd, Fc, false, post, post, L.parallel);
        }
    }

    // I - scale P B P on the fine unknowns (numbered by their rank): rows of a fine cell K: D_K (columns of K), E_K
    // (fine face columns); rows of a fine face unknown: G_R times the rows of D_K', E_K' of its neighbouring cells K'
    // fine face unknowns (rank order) whose cells are all fine: the interior faces of the fine region
    std::vector<char> interior_fine_faces(const local_data& L) const {
        std::vector<char> in(L.fine_dofs.size(), 0);
        for (size_t R = 0; R < m_faces.size(); R++) {
            const face_block& f = m_faces[R];
            if (f.n == 0 || L.rank[f.r0] < 0) continue;
            bool all = true;
            for (int k = 0; k < f.nc && all; k++) all = L.rank[m_cells[m_face_cells[f.c0 + k]].r0] >= 0;
            if (all)
                for (int q = 0; q < f.n; q++) in[L.rank[f.r0 + q]] = 1;
        }
        return in;
    }

    sparse_type fine_operator(const local_data& L, T scale, const std::vector<char>* constraint = nullptr) const {
        std::vector<Eigen::Triplet<T>> trips;
        const std::vector<int>& rank = L.rank;
        for (size_t K = 0; K < m_cells.size(); K++) {
            const cell_block& c = m_cells[K];
            const dense_type& D = m_D[c.d];
            const dense_type& E = m_E[c.e];
            for (int l = 0; l < c.n; l++) {
                const int row = rank[c.r0 + l];
                if (row < 0) continue;
                for (int l2 = 0; l2 < c.n; l2++)
                    if (rank[c.r0 + l2] >= 0 && D(l, l2) != 0) trips.emplace_back(row, rank[c.r0 + l2], -scale * D(l, l2));
                for (int j = 0; j < c.nj; j++) {
                    const int col = rank[m_J[c.j0 + j]];
                    if (col >= 0 && E(l, j) != 0) trips.emplace_back(row, col, -scale * E(l, j));
                }
            }
        }
        for (size_t R = 0; R < m_faces.size(); R++) {
            const face_block& f = m_faces[R];
            bool any = false;
            for (int q = 0; q < f.n; q++) any = any || rank[f.r0 + q] >= 0;
            if (!any) continue;
            const dense_type& G = m_G[f.g];
            const int* JF = m_JF.data() + f.j0;
            if (constraint && (*constraint)[rank[f.r0]]) {
                // static constraint u_f - G u_c[JF] = 0 (identity added below)
                for (int q = 0; q < f.n; q++)
                    for (int j = 0; j < f.nj; j++)
                        if (G(q, j) != 0) trips.emplace_back(rank[f.r0 + q], rank[JF[j]], -G(q, j));
                continue;
            }
            // columns of G grouped by neighbouring cell
            for (int j0 = 0; j0 < f.nj; ) {
                const int K = m_cell_of_dof[JF[j0]];
                int j1 = j0;
                while (j1 < f.nj && m_cell_of_dof[JF[j1]] == K) j1++;
                const cell_block& c = m_cells[K];
                const dense_type& D = m_D[c.d];
                const dense_type& E = m_E[c.e];
                dense_type Dr(j1 - j0, c.n), Er(j1 - j0, c.nj);
                for (int j = j0; j < j1; j++) {
                    Dr.row(j - j0) = D.row(JF[j] - c.r0);
                    Er.row(j - j0) = E.row(JF[j] - c.r0);
                }
                const dense_type GD = G.middleCols(j0, j1 - j0) * Dr, GE = G.middleCols(j0, j1 - j0) * Er;
                for (int q = 0; q < f.n; q++) {
                    const int row = rank[f.r0 + q];
                    if (row < 0) continue;
                    for (int l = 0; l < c.n; l++)
                        if (rank[c.r0 + l] >= 0 && GD(q, l) != 0) trips.emplace_back(row, rank[c.r0 + l], -scale * GD(q, l));
                    for (int j = 0; j < c.nj; j++) {
                        const int col = rank[m_J[c.j0 + j]];
                        if (col >= 0 && GE(q, j) != 0) trips.emplace_back(row, col, -scale * GE(q, j));
                    }
                }
                j0 = j1;
            }
        }
        const int nfd = int(L.fine_dofs.size());
        for (int q = 0; q < nfd; q++) trips.emplace_back(q, q, T(1));
        sparse_type A(nfd, nfd);
        A.setFromTriplets(trips.begin(), trips.end());
        A.makeCompressed();
        return A;
    }
};

#endif /* mf_erk_operator_hpp */
