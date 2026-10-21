//
//  system.hpp
//
//  Second-order formulation (HHOSecondOrder.hpp): the first-order system y = (u_T, v, d) of the second-order
//  equations (class second_order_system: interface splitting of the faces into dynamic and algebraic ones, block
//  inverses, algebraic faces by the splitting of Steins et al. or a factorization, right-hand side f(t, y), stage
//  systems of the implicit integrators, spectral estimates, local evaluations of B for faces.localize).
//

#pragma once
#ifndef second_order_system_hpp
#define second_order_system_hpp

#include "../mesh/domain.hpp"
#include "../physics/manufactured_solutions.hpp"
#include "../studies/convergence.hpp"
#include "../solvers/linear_solver.hpp"
#include "mf_operator.hpp"
#include <algorithm>
#include <memory>

namespace hho_second_order {

    using namespace coupling_common;
    typedef elastoacoustic_two_fields_assembler<mesh_type>  assembler_type;
    typedef SparseMatrix<RealType>                          sparse_type;
    typedef Matrix<RealType, Dynamic, Dynamic>              dense_type;
    typedef SparseMatrix<RealType, RowMajor>                sparse_rm_type;

    // Inverse of a block-diagonal sparse matrix with consecutive diagonal blocks of the given sizes. off_block returns
    // the Frobenius norm of the entries outside the blocks (0 if the matrix is block-diagonal).
    inline sparse_type block_diagonal_inverse(const sparse_type& A, const std::vector<size_t>& sizes, RealType& off_block) {
        std::vector<Triplet<RealType>> triplets;
        std::vector<size_t> block_of(A.rows());
        size_t ofs = 0;
        for (size_t b = 0; b < sizes.size(); ofs += sizes[b], b++)
            for (size_t i = 0; i < sizes[b]; i++)
                block_of[ofs + i] = b;
        RealType off2 = 0.0;
        for (int k = 0; k < A.outerSize(); ++k)
            for (sparse_type::InnerIterator it(A, k); it; ++it)
                if (block_of[it.row()] != block_of[it.col()])
                    off2 += it.value() * it.value();
        ofs = 0;
        for (size_t s : sizes) {
            dense_type blk = A.block(ofs, ofs, s, s).toDense();
            dense_type inv = blk.fullPivLu().inverse();
            for (size_t i = 0; i < s; i++)
                for (size_t j = 0; j < s; j++)
                    triplets.push_back(Triplet<RealType>(ofs + i, ofs + j, inv(i, j)));
            ofs += s;
        }
        assert(ofs == size_t(A.rows()));
        off_block = std::sqrt(off2);
        sparse_type Ainv(A.rows(), A.cols());
        Ainv.setFromTriplets(triplets.begin(), triplets.end());
        return Ainv;
    }

    // Sparse vector on [0, n) for the local evaluations (faces.localize): values, membership stamps, list of its entries
    struct sparse_acc {
        std::vector<RealType> v;
        std::vector<unsigned> mark;
        std::vector<int>      idx;
        unsigned              stamp = 1;
        void resize(size_t n) { v.assign(n, 0.0); mark.assign(n, 0); idx.clear(); stamp = 1; }
        void clear() {
            idx.clear();
            if (++stamp == 0) { std::fill(mark.begin(), mark.end(), 0u); stamp = 1; }
        }
        bool has(int i) const { return mark[i] == stamp; }
        RealType get(int i) const { return has(i) ? v[i] : 0.0; }
        RealType& at(int i) {
            if (mark[i] != stamp) { mark[i] = stamp; v[i] = 0.0; idx.push_back(i); }
            return v[i];
        }
        void sort() { std::sort(idx.begin(), idx.end()); }
    };

    // out += A x on the columns cols of A (column major, increasing), the rows kept by keep: same order of the
    // operations as the product of Eigen on the whole vector
    template<typename X, typename Keep>
    inline void scatter_product(const sparse_type& A, const std::vector<int>& cols, X&& x, sparse_acc& out, Keep&& keep) {
        for (int j : cols) {
            const RealType xj = x(j);
            for (sparse_type::InnerIterator it(A, j); it; ++it)
                if (keep(int(it.row()))) out.at(int(it.row())) += it.value() * xj;
        }
    }

    // First-order form of the semi-discrete problem, unknowns reordered as w = Q x = (u_T, d, a), Q orthogonal
    // (identity except the local QR bases on the elastic interface faces).
    class second_order_system {
    public:
        size_t n_c = 0, n_d = 0, n_a = 0, N = 0;
        sparse_type Q, Kw, Mcc, Cdd, Minv, Cinv;
        bool splitting_Q = true;
        RealType tol = 1.0e-14;
        size_t max_iter = 20000;
        // splitting statistics
        size_t n_solves = 0, n_iterations = 0, n_failures = 0, max_iterations = 0;

        // K: stiffness (stabilization weight gamma), M: mass, C: coupling, D: sigma* with weight gamma, all in the
        // assembler numbering
        second_order_system(mesh_type& msh, assembler_type& as, const sparse_type& K, const sparse_type& D, bool splitting) {

            splitting_Q = splitting;
            auto storage = msh.backend_storage();
            const size_t vfbs = as.get_e_face_basis_data(), sfbs = as.get_a_face_basis_data();
            n_c = as.get_e_n_cells_dof() + as.get_a_n_cells_dof();
            N = size_t(K.rows());
            const size_t e_face_ofs = n_c, a_face_ofs = n_c + as.get_e_face_dof();
            std::vector<size_t> e_compress = as.get_e_compress(), a_compress = as.get_a_compress();
            std::vector<size_t> e_expand = as.get_e_expand(), a_expand = as.get_a_expand();
            auto interfaces = as.get_interfaces();
            n_d = 2 * sfbs * interfaces.size();
            n_a = N - n_c - n_d;

            std::vector<Triplet<RealType>> q;
            for (size_t i = 0; i < n_c; i++)
                q.push_back(Triplet<RealType>(i, i, 1.0));

            // dynamic unknowns: per interface face, [normal part of the elastic unknowns, acoustic unknowns]
            size_t row_d = n_c, row_a = n_c + n_d;
            std::map<size_t, dense_type> Qf;     // local QR basis per interface face
            for (auto& chunk : interfaces) {
                auto& face = storage->edges[chunk.first];
                auto& e_cell = storage->surfaces[chunk.second.first];
                auto& a_cell = storage->surfaces[chunk.second.second];
                dense_type B = as.e_interface_operator(msh, face, e_cell, a_cell);
                HouseholderQR<dense_type> qr(B);
                dense_type Qloc = qr.householderQ() * dense_type::Identity(vfbs, vfbs);
                dense_type R = (Qloc.transpose() * B).topRows(sfbs);
                RealType rmin = R.diagonal().cwiseAbs().minCoeff();
                if (rmin < 1.0e-12 * B.norm())
                    throw std::runtime_error("second_order_system: rank-deficient interface coupling block");
                Qf[chunk.first] = Qloc;
                size_t e_ofs = e_face_ofs + e_compress.at(chunk.first) * vfbs;
                size_t a_ofs = a_face_ofs + a_compress.at(chunk.first) * sfbs;
                for (size_t i = 0; i < sfbs; i++)
                    for (size_t j = 0; j < vfbs; j++)
                        q.push_back(Triplet<RealType>(row_d + i, e_ofs + j, Qloc(j, i)));
                for (size_t i = 0; i < sfbs; i++)
                    q.push_back(Triplet<RealType>(row_d + sfbs + i, a_ofs + i, 1.0));
                row_d += 2 * sfbs;
                m_d_blocks.push_back(2 * sfbs);
            }
            // algebraic unknowns: elastic faces (tangential part on the interface), then acoustic faces off the interface
            for (size_t cf = 0; cf < as.get_elastic_faces(); cf++) {
                size_t fc_id = e_expand[cf];
                size_t e_ofs = e_face_ofs + cf * vfbs;
                if (interfaces.count(fc_id)) {
                    const dense_type& Qloc = Qf.at(fc_id);
                    for (size_t i = sfbs; i < vfbs; i++)
                        for (size_t j = 0; j < vfbs; j++)
                            q.push_back(Triplet<RealType>(row_a + i - sfbs, e_ofs + j, Qloc(j, i)));
                    row_a += vfbs - sfbs;
                    m_a_blocks.push_back(vfbs - sfbs);
                }
                else {
                    for (size_t i = 0; i < vfbs; i++)
                        q.push_back(Triplet<RealType>(row_a + i, e_ofs + i, 1.0));
                    row_a += vfbs;
                    m_a_blocks.push_back(vfbs);
                }
            }
            for (size_t cf = 0; cf < as.get_acoustic_faces(); cf++) {
                if (interfaces.count(a_expand[cf]))
                    continue;
                for (size_t i = 0; i < sfbs; i++)
                    q.push_back(Triplet<RealType>(row_a + i, a_face_ofs + cf * sfbs + i, 1.0));
                row_a += sfbs;
                m_a_blocks.push_back(sfbs);
            }
            assert(row_d == n_c + n_d && row_a == N);
            Q = sparse_type(N, N);
            Q.setFromTriplets(q.begin(), q.end());

            Kw = Q * K * Q.transpose();
            Kw.prune(0.0);
            Mcc = as.MASS.topLeftCorner(n_c, n_c);
            sparse_type Cw = Q * as.COUPLING * Q.transpose();
            Cdd = Cw.block(n_c, n_c, n_d, n_d);
            RealType leak2 = 0.0;
            for (int k = 0; k < Cw.outerSize(); ++k)
                for (sparse_type::InnerIterator it(Cw, k); it; ++it) {
                    bool in_dd = size_t(it.row()) >= n_c && size_t(it.row()) < n_c + n_d && size_t(it.col()) >= n_c && size_t(it.col()) < n_c + n_d;
                    if (!in_dd) leak2 += it.value() * it.value();
                }
            coupling_leak = std::sqrt(leak2) / std::max(Cw.norm(), 1.0e-300);

            std::vector<size_t> c_blocks;
            for (size_t i = 0; i < as.get_elastic_cells(); i++) c_blocks.push_back(as.get_e_cell_basis_data());
            for (size_t i = 0; i < as.get_acoustic_cells(); i++) c_blocks.push_back(as.get_a_cell_basis_data());
            m_cell_blocks = {{as.get_elastic_cells(), as.get_e_cell_basis_data()}, {as.get_acoustic_cells(), as.get_a_cell_basis_data()}};
            m_c_sizes = c_blocks;
            RealType off;
            Minv = block_diagonal_inverse(Mcc, c_blocks, off);
            Cinv = block_diagonal_inverse(Cdd, m_d_blocks, off);

            Kaa = Kw.block(n_c + n_d, n_c + n_d, n_a, n_a);
            if (splitting_Q) {
                sparse_type Dw = Q * D * Q.transpose();
                sparse_type Daa = Dw.block(n_c + n_d, n_c + n_d, n_a, n_a);
                Dinv = block_diagonal_inverse(Daa, m_a_blocks, off);
                sigma_off_block = off / std::max(Daa.norm(), 1.0e-300);
                Raa = Kaa - Daa;
            }
            else {
                m_Kaa_ldlt.compute(Kaa);
                if (m_Kaa_ldlt.info() != Success)
                    throw std::runtime_error("second_order_system: factorization of K_aa failed");
            }
        }
        second_order_system(const second_order_system&) = delete;
        second_order_system& operator=(const second_order_system&) = delete;

        RealType coupling_leak = 0.0;     // |C outside the dynamic block| / |C|, 0 up to round-off
        RealType sigma_off_block = 0.0;   // |sigma*_aa outside its diagonal blocks| / |sigma*_aa|, 0 up to round-off

        size_t state_size() const { return 2 * n_c + n_d; }

        // Fast mode (solver.fast): eval on row-major copies of the row blocks of K_w (OpenMP products of Eigen), the
        // splitting accelerated by preconditioned conjugate gradients (preconditioner D = gamma sigma*_aa: same
        // solution and same stopping quantity |D^{-1} (r - K_aa a)| <= tol |a|, the increment of one splitting step),
        // Pardiso for the direct faces, stage systems condensed on the faces and solved by Pardiso (linear_solver).
        // Results equal up to round-off (algebraic faces: up to faces.tolerance).
        void set_fast(bool fast) {
            m_fast = fast;
            if (!fast) return;
            const size_t nx = n_c + n_d;
            m_Kax = Kw.block(nx, 0, n_a, nx);
            m_Kxx = Kw.topLeftCorner(nx, nx);
            m_Kxa = Kw.block(0, nx, nx, n_a);
            m_Minv_rm = Minv;
            m_Cinv_rm = Cinv;
            if (splitting_Q) {
                m_Kaa_rm = Kaa;
                m_Dinv_rm = Dinv;
                m_Raa_rm = Raa;
            }
#ifdef HAVE_INTEL_MKL
            else {
                m_Kaa_pardiso = std::make_unique<PardisoLDLT<sparse_type>>();
                m_Kaa_pardiso->pardisoParameterArray()[24] = 2;    // parallel forward / backward substitution
                m_Kaa_pardiso->compute(Kaa);
                if (m_Kaa_pardiso->info() != Success)
                    throw std::runtime_error("second_order_system: factorization of K_aa failed");
            }
#endif
        }
        bool fast() const { return m_fast; }
        const std::vector<size_t>& algebraic_blocks() const { return m_a_blocks; }

        // Matrix-free products of eval_fast and of the conjugate gradients of the faces (solver.matrix_free, fast mode,
        // no elasto-acoustic interface; mf_operator.hpp): results equal up to round-off (blocks shared to
        // share_tolerance)
        void set_matrix_free(RealType share_tolerance, int max_batch) {
            if (!m_fast)
                throw std::runtime_error("second_order_system: the matrix-free products require the fast mode");
            if (n_d > 0)
                throw std::runtime_error("second_order_system: matrix-free products without elasto-acoustic interface only");
            m_mf = std::make_unique<mf_second_order_operator<RealType>>(Kw, n_c, Minv, splitting_Q ? Dinv : sparse_type(),
                                                                         m_c_sizes, m_a_blocks, share_tolerance, max_batch);
        }
        bool matrix_free() const { return bool(m_mf); }

        // Matrix-free rates of the velocities (set_matrix_free, n_d = 0): vdot = M^-1 (F_u - K_uu u - K_ua a), the faces a
        // solved from the cell values u (n_c values) with the warm start a_ws
        void mf_rates(const vector_type& Fw, const RealType* u, vector_type& a_ws, RealType* vdot) {
            m_mf_r.resize(Eigen::Index(n_a));
            m_mf->face_rhs(Fw.data() + n_c, u, m_mf_r.data());
            solve_faces(m_mf_r, a_ws);
            m_mf->cell_rates(Fw.data(), u, a_ws.data(), vdot);
        }
        std::string matrix_free_summary() const { return m_mf ? m_mf->summary() : std::string(); }

        // Algebraic faces: K_aa a = r, a holds the initial guess of the splitting on entry
        void solve_faces(const vector_type& r, vector_type& a) {
            if (m_fast) {
                solve_faces_fast(r, a);
                return;
            }
            if (!splitting_Q) {
                a = m_Kaa_ldlt.solve(r);
                return;
            }
            if (size_t(a.size()) != n_a) a = vector_type::Zero(n_a);
            n_solves++;
            if (r.squaredNorm() == 0.0) {
                a.setZero();
                return;
            }
            for (size_t m = 1; m <= max_iter; m++) {
                vector_type a_new = Dinv * (r - Raa * a);
                RealType inc = (a_new - a).norm();
                a = a_new;
                if (inc <= tol * a.norm()) {
                    n_iterations += m;
                    max_iterations = std::max(max_iterations, m);
                    return;
                }
            }
            n_iterations += max_iter;
            max_iterations = max_iter;
            n_failures++;
        }

        // f(t, y) for the forcing Fw = Q F(t); a_ws: warm start of the algebraic faces
        vector_type eval(const vector_type& Fw, const vector_type& y, vector_type& a_ws) {
            if (m_fast) return eval_fast(Fw, y, a_ws);
            vector_type z = vector_type::Zero(N);
            z.head(n_c) = y.head(n_c);
            z.segment(n_c, n_d) = y.tail(n_d);
            vector_type r = Fw.tail(n_a) - (Kw * z).tail(n_a);
            solve_faces(r, a_ws);
            z.tail(n_a) = a_ws;
            vector_type res = Fw - Kw * z;
            vector_type ydot(state_size());
            ydot.head(n_c) = y.segment(n_c, n_c);
            ydot.segment(n_c, n_c) = Minv * res.head(n_c);
            ydot.tail(n_d) = Cinv * res.segment(n_c, n_d);
            return ydot;
        }

        // Unknowns (cells and faces) in the assembler numbering: displacement part of y and faces at forcing Fw
        vector_type unknowns(const vector_type& Fw, const vector_type& y) {
            vector_type z = vector_type::Zero(N);
            z.head(n_c) = y.head(n_c);
            z.segment(n_c, n_d) = y.tail(n_d);
            vector_type r = Fw.tail(n_a) - (Kw * z).tail(n_a);
            vector_type a = vector_type::Zero(n_a);
            solve_faces(r, a);
            z.tail(n_a) = a;
            return Q.transpose() * z;
        }

        // Stage systems Z - h f(F, P Z) = r, P = diag(P_T, P_T, P_d) 0/1 (all ones: implicit RK stage). The unknowns
        // (h v, d, a) solve (K diag(P_T, P_d, 1) + M / h^2 + C_dd / h) X = F + (M r_v / h, C_dd r_d / h, 0) - K (P_T r_u, 0, 0)
        void factorize_stage(RealType h, const vector_type& Py) {
            m_h = h;
            m_Pc = Py.head(n_c);
            vector_type Pw = vector_type::Ones(N);
            Pw.head(n_c) = m_Pc;
            Pw.segment(n_c, n_d) = Py.tail(n_d);
            sparse_type A = Kw * Pw.asDiagonal();
            std::vector<Triplet<RealType>> t;
            for (int k = 0; k < Mcc.outerSize(); ++k)
                for (sparse_type::InnerIterator it(Mcc, k); it; ++it)
                    t.push_back(Triplet<RealType>(it.row(), it.col(), it.value() / (h * h)));
            for (int k = 0; k < Cdd.outerSize(); ++k)
                for (sparse_type::InnerIterator it(Cdd, k); it; ++it)
                    t.push_back(Triplet<RealType>(n_c + it.row(), n_c + it.col(), it.value() / h));
            sparse_type MC(N, N);
            MC.setFromTriplets(t.begin(), t.end());
            A += MC;
            A.makeCompressed();
            if (m_fast) {
                // cells (block-diagonal K_TT P_T + M / h^2) condensed, faces by Pardiso
                m_stage_ls = std::make_unique<linear_solver<RealType>>();
                m_stage_ls->set_Kg(A, N - n_c);
                m_stage_ls->set_fast(true);
                m_stage_ls->condense_equations(m_cell_blocks);
                m_stage_ls->release_condensed_blocks();
                m_stage_ls->factorize();
                return;
            }
            m_stage_lu.analyzePattern(A);
            m_stage_lu.factorize(A);
            if (m_stage_lu.info() != Success)
                throw std::runtime_error("second_order_system: factorization of the stage matrix failed");
        }

        vector_type stage_solve(const vector_type& Fw, const vector_type& r) {
            vector_type rhs = Fw;
            rhs.head(n_c) += Mcc * r.segment(n_c, n_c) / m_h;
            rhs.segment(n_c, n_d) += Cdd * r.tail(n_d) / m_h;
            vector_type z0 = vector_type::Zero(N);
            z0.head(n_c) = m_Pc.cwiseProduct(r.head(n_c));
            rhs -= Kw * z0;
            vector_type X = m_fast ? m_stage_ls->solve(rhs) : vector_type(m_stage_lu.solve(rhs));
            vector_type Z(state_size());
            Z.head(n_c) = r.head(n_c) + m_Pc.cwiseProduct(X.head(n_c));
            Z.segment(n_c, n_c) = X.head(n_c) / m_h;
            Z.tail(n_d) = X.segment(n_c, n_d);
            return Z;
        }

        // Spectral radius of the splitting iteration D^{-1} (K_aa - D) (power iteration; fast mode: stopped when stable),
        // the last value kept in last_splitting_radius
        RealType last_splitting_radius = 0.0;
        RealType splitting_spectral_radius(size_t n_it = 300) {
            last_splitting_radius = splitting_spectral_radius_estimate(n_it);
            return last_splitting_radius;
        }
        RealType splitting_spectral_radius_estimate(size_t n_it) {
            if (!splitting_Q) return 0.0;
            vector_type x = vector_type::Ones(n_a) + 0.5 * vector_type::Random(n_a);
            x /= x.norm();
            RealType rho = 0.0;
            for (size_t it = 0; it < n_it; it++) {
                vector_type y = m_fast ? vector_type(m_Dinv_rm * (m_Raa_rm * x)) : vector_type(Dinv * (Raa * x));
                vector_type z = m_fast ? vector_type(m_Dinv_rm * (m_Raa_rm * y)) : vector_type(Dinv * (Raa * y));
                RealType rho_old = rho;
                rho = std::sqrt(z.norm());          // two steps: robust to a +-rho pair
                x = z / z.norm();
                if (m_fast && it >= 20 && std::abs(rho - rho_old) <= 1.0e-6 * rho) break;   // fast mode: estimate stable
            }
            return rho;
        }

        // Largest frequency omega_max of y' = B y (eigenvalues +-i omega), power iteration on B^2 (fast mode: stopped when
        // stable to 1e-6)
        RealType max_frequency(size_t n_it = 100) {
            vector_type y = vector_type::Ones(state_size()) + 0.5 * vector_type::Random(state_size());
            y /= y.norm();
            vector_type zero = vector_type::Zero(N), a_ws;
            RealType lambda = 0.0;
            for (size_t it = 0; it < n_it; it++) {
                vector_type z = eval(zero, eval(zero, y, a_ws), a_ws);
                RealType lambda_old = lambda;
                lambda = z.norm();
                y = z / lambda;
                if (m_fast && it >= 20 && std::abs(lambda - lambda_old) <= 1.0e-6 * lambda) break;   // fast mode: estimate stable
            }
            return std::sqrt(lambda);
        }

        // ------------------------------------------------------------------------------------ local evaluations
        // faces.localize: B y without force on a state read on a set of cells (u_T, v_T) and of interface blocks (d) only,
        // the algebraic faces from zero (splitting, or the conjugate gradients of the fast mode) on the faces reached,
        // one ring of faces more per iteration, restricted to the faces of cap if given (truncation), the products on
        // the rows reached only. Without cap, equal to eval(0, y masked to the sets, a = 0) up to round-off.
        // Blocks: cells c (unknowns [c_ofs, c_ofs + c_size) of u_T), interface blocks b (of d), algebraic blocks of a.
        struct local_graph {
            std::vector<size_t> c_ofs, c_size, d_ofs, d_size, a_ofs, a_size;
            std::vector<int>    c_of, d_of, a_of;                    // block of an unknown of u_T, d, a
            std::vector<std::vector<int>> c_a, c_d, a_c, d_c;       // algebraic / interface blocks of a cell, cells of a block
        };
        size_t local_evaluations = 0, local_faces = 0, local_rows = 0;   // statistics: evaluations, faces and rows reached

        bool local_ready() const { return m_local_ready; }
        const local_graph& graph() const { return m_g; }

        void build_local() {
            if (m_local_ready) return;
            if (!splitting_Q)
                throw std::runtime_error("faces.localize: requires faces.solver = \"splitting\"");
            const size_t nx = n_c + n_d;
            local_graph& g = m_g;
            for (auto& cb : m_cell_blocks)
                for (size_t i = 0; i < cb.first; i++) { g.c_ofs.push_back(g.c_ofs.empty() ? 0 : g.c_ofs.back() + g.c_size.back()); g.c_size.push_back(cb.second); }
            for (size_t s : m_d_blocks) { g.d_ofs.push_back(g.d_ofs.empty() ? 0 : g.d_ofs.back() + g.d_size.back()); g.d_size.push_back(s); }
            for (size_t s : m_a_blocks) { g.a_ofs.push_back(g.a_ofs.empty() ? 0 : g.a_ofs.back() + g.a_size.back()); g.a_size.push_back(s); }
            g.c_of.assign(n_c, -1); g.d_of.assign(n_d, -1); g.a_of.assign(n_a, -1);
            for (size_t c = 0; c < g.c_ofs.size(); c++) for (size_t k = 0; k < g.c_size[c]; k++) g.c_of[g.c_ofs[c] + k] = int(c);
            for (size_t b = 0; b < g.d_ofs.size(); b++) for (size_t k = 0; k < g.d_size[b]; k++) g.d_of[g.d_ofs[b] + k] = int(b);
            for (size_t b = 0; b < g.a_ofs.size(); b++) for (size_t k = 0; k < g.a_size[b]; k++) g.a_of[g.a_ofs[b] + k] = int(b);
            m_Kax_c = Kw.block(nx, 0, n_a, nx);
            m_Kxx_c = Kw.topLeftCorner(nx, nx);
            m_Kxa_c = Kw.block(0, nx, nx, n_a);
            // cells of the algebraic blocks (rows of K_xa) and of the interface blocks (rows of K_xx in the d columns)
            g.a_c.assign(g.a_ofs.size(), {}); g.d_c.assign(g.d_ofs.size(), {});
            for (int j = 0; j < m_Kxa_c.outerSize(); ++j)
                for (sparse_type::InnerIterator it(m_Kxa_c, j); it; ++it)
                    if (size_t(it.row()) < n_c) g.a_c[g.a_of[j]].push_back(g.c_of[it.row()]);
            for (int j = int(n_c); j < int(nx); ++j)
                for (sparse_type::InnerIterator it(m_Kxx_c, j); it; ++it)
                    if (size_t(it.row()) < n_c) g.d_c[g.d_of[j - n_c]].push_back(g.c_of[it.row()]);
            auto unique = [](std::vector<int>& v) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); };
            g.c_a.assign(g.c_ofs.size(), {}); g.c_d.assign(g.c_ofs.size(), {});
            for (size_t b = 0; b < g.a_c.size(); b++) { unique(g.a_c[b]); for (int c : g.a_c[b]) g.c_a[c].push_back(int(b)); }
            for (size_t b = 0; b < g.d_c.size(); b++) { unique(g.d_c[b]); for (int c : g.d_c[b]) g.c_d[c].push_back(int(b)); }
            for (auto& v : g.c_a) unique(v);
            for (auto& v : g.c_d) unique(v);
            for (sparse_acc* acc : {&m_lr, &m_la, &m_lan, &m_lt, &m_lz, &m_lp, &m_lq}) acc->resize(n_a);
            m_lres.resize(nx);
            m_local_ready = true;
        }

        // Cells within r rings of the cells in (neighbours: cells sharing an algebraic or an interface face)
        std::vector<char> grow_cells(std::vector<char> in, size_t r) const {
            const local_graph& g = m_g;
            std::vector<int> front;
            for (size_t c = 0; c < in.size(); c++) if (in[c]) front.push_back(int(c));
            std::vector<char> fa(g.a_ofs.size(), 0), fd(g.d_ofs.size(), 0);
            for (size_t ring = 0; ring < r && !front.empty(); ring++) {
                std::vector<int> next;
                for (int c : front) {
                    for (int b : g.c_a[c]) if (!fa[b]) { fa[b] = 1; for (int c2 : g.a_c[b]) if (!in[c2]) { in[c2] = 1; next.push_back(c2); } }
                    for (int b : g.c_d[c]) if (!fd[b]) { fd[b] = 1; for (int c2 : g.d_c[b]) if (!in[c2]) { in[c2] = 1; next.push_back(c2); } }
                }
                front.swap(next);
            }
            return in;
        }
        // cells (u_T part of a 0/1 state mask), interface blocks (d part of the mask)
        std::vector<char> mask_cells(const vector_type& mask) const {
            std::vector<char> in(m_g.c_ofs.size(), 0);
            for (size_t c = 0; c < in.size(); c++) in[c] = mask(Eigen::Index(m_g.c_ofs[c])) > 0.5;
            return in;
        }
        std::vector<int> mask_dblocks(const vector_type& mask) const {
            std::vector<int> out;
            for (size_t b = 0; b < m_g.d_ofs.size(); b++) if (mask(Eigen::Index(2 * n_c + m_g.d_ofs[b])) > 0.5) out.push_back(int(b));
            return out;
        }
        // interface blocks with a cell in, algebraic blocks with a cell in (cap of a truncated evaluation)
        std::vector<int> cells_dblocks(const std::vector<char>& in) const {
            std::vector<int> out;
            for (size_t b = 0; b < m_g.d_c.size(); b++)
                for (int c : m_g.d_c[b]) if (in[c]) { out.push_back(int(b)); break; }
            return out;
        }
        std::vector<char> cells_cap(const std::vector<char>& in) const {
            std::vector<char> cap(m_g.a_ofs.size(), 0);
            for (size_t b = 0; b < m_g.a_c.size(); b++)
                for (int c : m_g.a_c[b]) if (in[c]) { cap[b] = 1; break; }
            return cap;
        }
        static std::vector<int> list_of(const std::vector<char>& in) {
            std::vector<int> out;
            for (size_t c = 0; c < in.size(); c++) if (in[c]) out.push_back(int(c));
            return out;
        }

        // B y (no force) reading y on the cells and the interface blocks given (increasing), faces capped to cap if any
        vector_type eval_local(const vector_type& y, const std::vector<int>& cells, const std::vector<int>& dblocks,
                               const std::vector<char>* cap) {
            const local_graph& g = m_g;
            std::vector<int>& xc = m_lxc;
            xc.clear();
            for (int c : cells) for (size_t k = 0; k < g.c_size[c]; k++) xc.push_back(int(g.c_ofs[c] + k));
            for (int b : dblocks) for (size_t k = 0; k < g.d_size[b]; k++) xc.push_back(int(n_c + g.d_ofs[b] + k));
            auto xval = [&](int j) { return j < int(n_c) ? y(j) : y(Eigen::Index(2 * n_c) + (j - Eigen::Index(n_c))); };
            auto keep_a = [&](int i) { return !cap || (*cap)[g.a_of[i]]; };
            // r = 0 - K_ax x
            sparse_acc& r = m_lr;
            r.clear();
            scatter_product(m_Kax_c, xc, xval, r, keep_a);
            for (int i : r.idx) r.v[i] = 0.0 - r.v[i];
            r.sort();
            sparse_acc& a = solve_faces_local(r, keep_a);
            // res = 0 - (K_xx x + K_xa a), rows reached
            sparse_acc& res = m_lres;
            res.clear();
            scatter_product(m_Kxx_c, xc, xval, res, [](int) { return true; });
            a.sort();
            scatter_product(m_Kxa_c, a.idx, [&](int j) { return a.v[j]; }, res, [](int) { return true; });
            res.sort();
            vector_type out = vector_type::Zero(state_size());
            for (int c : cells) for (size_t k = 0; k < g.c_size[c]; k++) { const Eigen::Index i = Eigen::Index(g.c_ofs[c] + k); out(i) = y(Eigen::Index(n_c) + i); }
            for (int j : res.idx) {
                const RealType rj = 0.0 - res.v[j];
                if (j < int(n_c))
                    for (sparse_type::InnerIterator it(Minv, j); it; ++it) out(Eigen::Index(n_c) + it.row()) += it.value() * rj;
                else
                    for (sparse_type::InnerIterator it(Cinv, j - int(n_c)); it; ++it) out(Eigen::Index(2 * n_c) + it.row()) += it.value() * rj;
            }
            local_evaluations++;
            local_faces += a.idx.size();
            local_rows += res.idx.size();
            return out;
        }

        // Local block (faces.projected): B~ y on the cells and interface blocks given, the algebraic faces of cap solved
        // exactly (sparse LDL^T of K_aa restricted to them, factorized once), zero on the other faces; the products on
        // submatrices restricted to the unknowns read, the faces solved and the rows reached (whole cell and interface
        // blocks)
        struct local_block {
            std::vector<int> xi, ci, rc, rd, vc;     // x read, a solved, cell rows (u_T), d rows (d numbering), velocities copied
            sparse_type Kax, Kxx, Kxa, Minv, Cinv;
            SimplicialLDLT<sparse_type> ldlt;
            size_t faces = 0;                         // algebraic blocks solved
        };
        std::unique_ptr<local_block> make_local_block(const std::vector<int>& cells, const std::vector<int>& dblocks,
                                                      const std::vector<char>& cap) {
            const local_graph& g = m_g;
            const size_t nx = n_c + n_d;
            auto lb = std::make_unique<local_block>();
            for (int c : cells) for (size_t k = 0; k < g.c_size[c]; k++) { lb->xi.push_back(int(g.c_ofs[c] + k)); lb->vc.push_back(int(g.c_ofs[c] + k)); }
            for (int b : dblocks) for (size_t k = 0; k < g.d_size[b]; k++) lb->xi.push_back(int(n_c + g.d_ofs[b] + k));
            for (size_t b = 0; b < cap.size(); b++)
                if (cap[b]) { lb->faces++; for (size_t k = 0; k < g.a_size[b]; k++) lb->ci.push_back(int(g.a_ofs[b] + k)); }
            // rows reached by K_xx (columns read) and K_xa (faces solved): whole cell / interface blocks
            std::vector<char> rcb(g.c_ofs.size(), 0), rdb(g.d_ofs.size(), 0);
            auto mark = [&](int i) { if (i < int(n_c)) rcb[g.c_of[i]] = 1; else rdb[g.d_of[i - int(n_c)]] = 1; };
            for (int j : lb->xi) for (sparse_type::InnerIterator it(m_Kxx_c, j); it; ++it) mark(int(it.row()));
            for (int j : lb->ci) for (sparse_type::InnerIterator it(m_Kxa_c, j); it; ++it) mark(int(it.row()));
            for (size_t c = 0; c < rcb.size(); c++) if (rcb[c]) for (size_t k = 0; k < g.c_size[c]; k++) lb->rc.push_back(int(g.c_ofs[c] + k));
            for (size_t b = 0; b < rdb.size(); b++) if (rdb[b]) for (size_t k = 0; k < g.d_size[b]; k++) lb->rd.push_back(int(g.d_ofs[b] + k));
            // local numberings
            std::vector<int> xmap(nx, -1), amap(n_a, -1), rmap(nx, -1);
            for (size_t k = 0; k < lb->xi.size(); k++) xmap[lb->xi[k]] = int(k);
            for (size_t k = 0; k < lb->ci.size(); k++) amap[lb->ci[k]] = int(k);
            for (size_t k = 0; k < lb->rc.size(); k++) rmap[lb->rc[k]] = int(k);
            for (size_t k = 0; k < lb->rd.size(); k++) rmap[n_c + lb->rd[k]] = int(lb->rc.size() + k);
            const int nX = int(lb->xi.size()), nC = int(lb->ci.size()), nR = int(lb->rc.size() + lb->rd.size());
            auto extract = [](const sparse_type& A, const std::vector<int>& cols, const std::vector<int>& rowmap, int nrows, const std::vector<int>* colmap) {
                std::vector<Triplet<RealType>> t;
                for (size_t k = 0; k < cols.size(); k++)
                    for (sparse_type::InnerIterator it(A, cols[k]); it; ++it) {
                        const int r = rowmap[it.row()];
                        if (r >= 0) t.emplace_back(r, colmap ? (*colmap)[cols[k]] : int(k), it.value());
                    }
                sparse_type B(nrows, Eigen::Index(cols.size()));
                B.setFromTriplets(t.begin(), t.end());
                return B;
            };
            lb->Kax = extract(m_Kax_c, lb->xi, amap, nC, nullptr);
            lb->Kxx = extract(m_Kxx_c, lb->xi, rmap, nR, nullptr);
            lb->Kxa = extract(m_Kxa_c, lb->ci, rmap, nR, nullptr);
            sparse_type Kcc = extract(Kaa, lb->ci, amap, nC, nullptr);
            std::vector<int> cmap(n_c, -1), dmap(n_d, -1);
            for (size_t k = 0; k < lb->rc.size(); k++) cmap[lb->rc[k]] = int(k);
            for (size_t k = 0; k < lb->rd.size(); k++) dmap[lb->rd[k]] = int(k);
            lb->Minv = extract(Minv, lb->rc, cmap, int(lb->rc.size()), nullptr);
            lb->Cinv = extract(Cinv, lb->rd, dmap, int(lb->rd.size()), nullptr);
            if (nC > 0) {
                lb->ldlt.compute(Kcc);
                if (lb->ldlt.info() != Success)
                    throw std::runtime_error("faces.projected: factorization of the local face block failed");
            }
            (void)nX;
            return lb;
        }

        // Matrix-free operator of a local block (faces.projected with solver.matrix_free, n_d = 0, splitting): the cells of
        // its rows rc (their cell unknowns outside the block read as zero) and its faces ci, with the sub-blocks of K_w,
        // M^-1 and D^-1; same operator as apply_local_block, faces by conjugate gradients instead of the factorization
        std::unique_ptr<mf_second_order_operator<RealType>> make_local_mf(const local_block& lb, RealType share_tolerance,
                                                                         int max_batch) const {
            if (n_d > 0 || !splitting_Q)
                throw std::runtime_error("faces.projected with solver.matrix_free: no elasto-acoustic interface, splitting only");
            const local_graph& g = m_g;
            const Eigen::Index nr = Eigen::Index(lb.rc.size()), nf = Eigen::Index(lb.ci.size());
            std::vector<Triplet<RealType>> t, tc, ta;
            for (Eigen::Index k = 0; k < nr; k++) {
                t.emplace_back(k, lb.rc[k], 1.0);
                tc.emplace_back(k, lb.rc[k], 1.0);
            }
            for (Eigen::Index k = 0; k < nf; k++) {
                t.emplace_back(nr + k, Eigen::Index(n_c + n_d) + lb.ci[k], 1.0);
                ta.emplace_back(k, lb.ci[k], 1.0);
            }
            sparse_type P(nr + nf, Eigen::Index(N)), Pc(nr, Eigen::Index(n_c)), Pa(nf, Eigen::Index(n_a));
            P.setFromTriplets(t.begin(), t.end());
            Pc.setFromTriplets(tc.begin(), tc.end());
            Pa.setFromTriplets(ta.begin(), ta.end());
            sparse_type Ks = P * Kw * P.transpose(), Ms = Pc * Minv * Pc.transpose(), Ds = Pa * Dinv * Pa.transpose();
            std::vector<size_t> cs, fs;
            for (Eigen::Index k = 0; k < nr; k++)
                if (k == 0 || g.c_of[lb.rc[k]] != g.c_of[lb.rc[k - 1]]) cs.push_back(1); else cs.back()++;
            for (Eigen::Index k = 0; k < nf; k++)
                if (k == 0 || g.a_of[lb.ci[k]] != g.a_of[lb.ci[k - 1]]) fs.push_back(1); else fs.back()++;
            return std::make_unique<mf_second_order_operator<RealType>>(Ks, size_t(nr), Ms, Ds, cs, fs, share_tolerance, max_batch);
        }
        vector_type apply_local_block(local_block& lb, const vector_type& y) const {
            vector_type x(Eigen::Index(lb.xi.size()));
            for (size_t k = 0; k < lb.xi.size(); k++) {
                const int j = lb.xi[k];
                x(Eigen::Index(k)) = j < int(n_c) ? y(j) : y(Eigen::Index(2 * n_c) + (j - Eigen::Index(n_c)));
            }
            vector_type res = -(lb.Kxx * x);
            if (!lb.ci.empty()) {
                const vector_type a = lb.ldlt.solve(-(lb.Kax * x));
                res -= lb.Kxa * a;
            }
            vector_type out = vector_type::Zero(state_size());
            for (int i : lb.vc) out(i) = y(Eigen::Index(n_c) + i);
            const Eigen::Index nrc = Eigen::Index(lb.rc.size());
            const vector_type vc = lb.Minv * res.head(nrc), vd = lb.Cinv * res.tail(res.size() - nrc);
            for (size_t k = 0; k < lb.rc.size(); k++) out(Eigen::Index(n_c) + lb.rc[k]) = vc(Eigen::Index(k));
            for (size_t k = 0; k < lb.rd.size(); k++) out(Eigen::Index(2 * n_c) + lb.rd[k]) = vd(Eigen::Index(k));
            return out;
        }

    private:
        std::vector<size_t> m_d_blocks, m_a_blocks;
        sparse_type Kaa, Raa, Dinv;
        SimplicialLDLT<sparse_type> m_Kaa_ldlt;
        SparseLU<sparse_type> m_stage_lu;
        RealType m_h = 0.0;
        vector_type m_Pc;

        // f(t, y) of the fast mode: x = (u_T, d), r = F_a - K_ax x, faces, residuals of the rows (T, d) only
        vector_type eval_fast(const vector_type& Fw, const vector_type& y, vector_type& a_ws) {
            if (m_mf) {     // n_d = 0: x = u_T
                vector_type ydot(state_size());
                ydot.head(n_c) = y.segment(n_c, n_c);
                mf_rates(Fw, y.data(), a_ws, ydot.data() + n_c);
                return ydot;
            }
            const size_t nx = n_c + n_d;
            vector_type x(nx);
            x.head(n_c) = y.head(n_c);
            x.tail(n_d) = y.tail(n_d);
            vector_type r = Fw.tail(n_a) - m_Kax * x;
            solve_faces(r, a_ws);
            vector_type res = Fw.head(nx) - m_Kxx * x;
            if (n_a > 0) res -= m_Kxa * a_ws;
            vector_type ydot(state_size());
            ydot.head(n_c) = y.segment(n_c, n_c);
            ydot.segment(n_c, n_c) = m_Minv_rm * res.head(n_c);
            ydot.tail(n_d) = m_Cinv_rm * res.tail(n_d);
            return ydot;
        }

        // Algebraic faces of the fast mode: preconditioned conjugate gradients from the initial guess a (splitting), or
        // Pardiso (direct)
        void solve_faces_fast(const vector_type& r, vector_type& a) {
            if (!splitting_Q) {
#ifdef HAVE_INTEL_MKL
                a = m_Kaa_pardiso->solve(r);
#else
                a = m_Kaa_ldlt.solve(r);
#endif
                return;
            }
            if (size_t(a.size()) != n_a) a = vector_type::Zero(n_a);
            n_solves++;
            if (r.squaredNorm() == 0.0) {
                a.setZero();
                return;
            }
            if (m_mf) {     // matrix-free conjugate gradients (fused OpenMP passes, same iterates up to round-off)
                const size_t m = m_mf->pcg(r.data(), a.data(), tol, max_iter);
                n_iterations += m ? m : max_iter;
                max_iterations = m ? std::max(max_iterations, m) : max_iter;
                if (!m) n_failures++;
                return;
            }
            vector_type res = r - m_Kaa_rm * a;
            vector_type z = m_Dinv_rm * res;
            if (z.norm() <= tol * a.norm()) {      // converged: one splitting step, as the splitting
                a += z;
                n_iterations += 1;
                max_iterations = std::max(max_iterations, size_t(1));
                return;
            }
            vector_type p = z, q;
            RealType rz = res.dot(z);
            for (size_t m = 1; m <= max_iter; m++) {
                q = m_Kaa_rm * p;
                RealType alpha = rz / p.dot(q);
                a += alpha * p;
                res -= alpha * q;
                z = m_Dinv_rm * res;
                if (z.norm() <= tol * a.norm()) {
                    n_iterations += m;
                    max_iterations = std::max(max_iterations, m);
                    return;
                }
                RealType rz_new = res.dot(z);
                p = z + (rz_new / rz) * p;
                rz = rz_new;
            }
            n_iterations += max_iter;
            max_iterations = max_iter;
            n_failures++;
        }

        // K_aa a = r from a = 0 on the faces kept by keep (local evaluations): the splitting of solve_faces, or the
        // conjugate gradients of solve_faces_fast in the fast mode, on sparse vectors
        template<typename Keep>
        sparse_acc& solve_faces_local(const sparse_acc& r, Keep&& keep) {
            sparse_acc& a = m_la;
            a.clear();
            n_solves++;
            RealType r2 = 0.0;
            for (int i : r.idx) r2 += r.v[i] * r.v[i];
            if (r2 == 0.0) return a;
            auto norm = [](const sparse_acc& u) { RealType s = 0.0; for (int i : u.idx) s += u.v[i] * u.v[i]; return std::sqrt(s); };
            // z = D^{-1} u on the entries of u (increasing)
            auto apply_dinv = [&](const sparse_acc& u, sparse_acc& z) {
                z.clear();
                scatter_product(Dinv, u.idx, [&](int j) { return u.v[j]; }, z, [](int) { return true; });
                z.sort();
            };
            if (!m_fast) {
                sparse_acc& an = m_lan;
                sparse_acc& t = m_lt;
                for (size_t m = 1; m <= max_iter; m++) {
                    // t = r - R_aa a, a_new = D^{-1} t
                    t.clear();
                    scatter_product(Raa, a.idx, [&](int j) { return a.v[j]; }, t, keep);
                    for (int i : t.idx) t.v[i] = r.get(i) - t.v[i];
                    for (int i : r.idx) if (!t.has(i)) t.at(i) = r.v[i] - 0.0;
                    t.sort();
                    apply_dinv(t, an);
                    RealType inc2 = 0.0;
                    for (int i : an.idx) { const RealType d = an.v[i] - a.get(i); inc2 += d * d; }
                    for (int i : a.idx) if (!an.has(i)) inc2 += a.v[i] * a.v[i];
                    // a = a_new
                    a.clear();
                    for (int i : an.idx) a.at(i) = an.v[i];
                    if (std::sqrt(inc2) <= tol * norm(a)) {
                        n_iterations += m;
                        max_iterations = std::max(max_iterations, m);
                        return a;
                    }
                }
                n_iterations += max_iter;
                max_iterations = max_iter;
                n_failures++;
                return a;
            }
            // conjugate gradients preconditioned by D, from a = 0 (res = r)
            sparse_acc& res = m_lt;
            sparse_acc& z = m_lz;
            sparse_acc& p = m_lp;
            sparse_acc& q = m_lq;
            res.clear();
            for (int i : r.idx) res.at(i) = r.v[i];
            res.sort();
            apply_dinv(res, z);
            if (norm(z) <= tol * norm(a)) {
                for (int i : z.idx) a.at(i) += z.v[i];
                n_iterations += 1;
                max_iterations = std::max(max_iterations, size_t(1));
                return a;
            }
            auto dot = [](const sparse_acc& u, const sparse_acc& w) { RealType s = 0.0; for (int i : u.idx) s += u.v[i] * w.get(i); return s; };
            p.clear();
            for (int i : z.idx) p.at(i) = z.v[i];
            RealType rz = dot(res, z);
            for (size_t m = 1; m <= max_iter; m++) {
                p.sort();
                q.clear();
                scatter_product(Kaa, p.idx, [&](int j) { return p.v[j]; }, q, keep);
                const RealType alpha = rz / dot(p, q);
                for (int i : p.idx) a.at(i) += alpha * p.v[i];
                for (int i : q.idx) res.at(i) -= alpha * q.v[i];
                res.sort();
                apply_dinv(res, z);
                if (norm(z) <= tol * norm(a)) {
                    n_iterations += m;
                    max_iterations = std::max(max_iterations, m);
                    return a;
                }
                const RealType rz_new = dot(res, z), beta = rz_new / rz;
                for (int i : p.idx) p.v[i] *= beta;
                for (int i : z.idx) p.at(i) += z.v[i];
                rz = rz_new;
            }
            n_iterations += max_iter;
            max_iterations = max_iter;
            n_failures++;
            return a;
        }

        // local evaluations (build_local)
        bool m_local_ready = false;
        local_graph m_g;
        sparse_type m_Kax_c, m_Kxx_c, m_Kxa_c;
        sparse_acc m_lr, m_la, m_lan, m_lt, m_lz, m_lp, m_lq, m_lres;
        std::vector<int> m_lxc;

        // fast mode (set_fast); declared last
        bool m_fast = false;
        std::vector<std::pair<size_t, size_t>> m_cell_blocks;     // (cells, cell basis size): elastic, acoustic
        std::vector<size_t> m_c_sizes;                            // cell basis size of each cell
        std::unique_ptr<mf_second_order_operator<RealType>> m_mf; // matrix-free products (set_matrix_free)
        vector_type m_mf_r;                                       // right-hand side of the faces (mf_rates)
        sparse_rm_type m_Kax, m_Kxx, m_Kxa, m_Minv_rm, m_Cinv_rm, m_Kaa_rm, m_Dinv_rm, m_Raa_rm;
#ifdef HAVE_INTEL_MKL
        std::unique_ptr<PardisoLDLT<sparse_type>> m_Kaa_pardiso;
#endif
        std::unique_ptr<linear_solver<RealType>> m_stage_ls;
    };

}

#endif /* second_order_system_hpp */
