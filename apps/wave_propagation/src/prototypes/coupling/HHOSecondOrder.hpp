//
//  HHOSecondOrder.hpp
//
//  Elasto-acoustic coupling, second-order formulation (elastoacoustic_two_fields_assembler: displacement in the
//  solid, velocity potential in the fluid), M u'' + C u' + K u = F(t). Every setting comes from the Lua input file
//  (common/lua_parameters.hpp, template share/second_order.lua). The Runge-Kutta integrators solve a first-order
//  system: the auxiliary variable is the cell velocity v = u_T',
//      u_T' = v,     M v' = F_T - K_TT u_T - K_TF u_F,     C_FF u_F' = F_F - K_FT u_T - K_FF u_F,
//  M being the cell mass. The coupling C = [[0, B], [-B^T, 0]] only acts on the interface faces: a QR factorization of
//  B on each interface face splits the elastic face unknowns into a normal part (range of B) and a tangential part.
//  The normal part and the acoustic interface unknowns, d, are dynamic (C_dd invertible); the other face unknowns, a
//  (all the other faces, tangential part on the interface), are algebraic, K_aa a = F_a - K_aT u_T - K_ad d.
//  State y = (u_T, v, d), y' = f(t, y) = B y + g(t).
//    space : face degree k, equal (k,k) or mixed (k+1,k) order, least-squares (plain least squares / Lehrenfeld-
//            Schoberl) or HHO stabilization, weight in 1/h or O(1), weight gamma (number, or "auto" = weight_factor
//            gamma*, gamma* the critical weight of the splitting computed cell by cell, Steins et al. Lemma 4.3)
//    faces : explicit integrators (erk, lts) and li: algebraic faces from the splitting of Steins et al. (M2AN 57,
//            2023), D a^{m+1} = r - (K_aa - D) a^m with D = gamma sigma*, the block-diagonal face part of the
//            stabilization (fully explicit, contraction gamma*/gamma, relative increment faces.tolerance), or from a
//            sparse factorization of K_aa ("direct"); sdirk solves them within its stage systems
//    time  : ERK(s), SDIRK(s,s+1), two-level LTS-RK4, LI-SDIRK(4,3), final time, number of time steps
//    mesh, problem (manufactured solution), study (one run, space-time convergence, time convergence), output
//  The SDIRK stages and the LI fine solves, (I - a_ii dt B P) Y = R, are solved with a, u_T eliminated, through the
//  sparse system (K + M / (a_ii dt)^2 + C / (a_ii dt)) restricted by the projector P (P = I for SDIRK). The explicit
//  integrators print omega_max (power iteration) and warn when dt is above their stability limit.
//  Steins' optimal explicit setting: mixed order, least squares, 1/h, weight "auto" with weight_factor in [3, 5],
//  splitting; implicit: mixed order, least squares (Lehrenfeld-Schoberl), 1/h, weight 1. With equal order and the
//  HHO stabilization, alpha = lambda_max(S_aa, sigma*_aa) - 1 >= 1 for k = 1, 3 in the solid: the splitting cannot
//  converge (use faces.solver = "direct").
//

#pragma once
#ifndef HHOSecondOrder_hpp
#define HHOSecondOrder_hpp

#include "coupling_common.hpp"

namespace hho_second_order {

    using namespace coupling_common;
    typedef elastoacoustic_two_fields_assembler<mesh_type>  assembler_type;
    typedef SparseMatrix<RealType>                          sparse_type;
    typedef Matrix<RealType, Dynamic, Dynamic>              dense_type;

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
        const std::vector<size_t>& algebraic_blocks() const { return m_a_blocks; }

        // Algebraic faces: K_aa a = r, a holds the initial guess of the splitting on entry
        void solve_faces(const vector_type& r, vector_type& a) {
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
            vector_type X = m_stage_lu.solve(rhs);
            vector_type Z(state_size());
            Z.head(n_c) = r.head(n_c) + m_Pc.cwiseProduct(X.head(n_c));
            Z.segment(n_c, n_c) = X.head(n_c) / m_h;
            Z.tail(n_d) = X.segment(n_c, n_d);
            return Z;
        }

        // Spectral radius of the splitting iteration D^{-1} (K_aa - D) (power iteration)
        RealType splitting_spectral_radius(size_t n_it = 300) {
            if (!splitting_Q) return 0.0;
            vector_type x = vector_type::Ones(n_a) + 0.5 * vector_type::Random(n_a);
            x /= x.norm();
            RealType rho = 0.0;
            for (size_t it = 0; it < n_it; it++) {
                vector_type y = Dinv * (Raa * x);
                vector_type z = Dinv * (Raa * y);
                rho = std::sqrt(z.norm());          // two steps: robust to a +-rho pair
                x = z / z.norm();
            }
            return rho;
        }

        // Largest frequency omega_max of y' = B y (eigenvalues +-i omega), power iteration on B^2
        RealType max_frequency(size_t n_it = 100) {
            vector_type y = vector_type::Ones(state_size()) + 0.5 * vector_type::Random(state_size());
            y /= y.norm();
            vector_type zero = vector_type::Zero(N), a_ws;
            RealType lambda = 0.0;
            for (size_t it = 0; it < n_it; it++) {
                vector_type z = eval(zero, eval(zero, y, a_ws), a_ws);
                lambda = z.norm();
                y = z / lambda;
            }
            return std::sqrt(lambda);
        }

    private:
        std::vector<size_t> m_d_blocks, m_a_blocks;
        sparse_type Kaa, Raa, Dinv;
        SimplicialLDLT<sparse_type> m_Kaa_ldlt;
        SparseLU<sparse_type> m_stage_lu;
        RealType m_h = 0.0;
        vector_type m_Pc;
    };

    // Coupled problem on one mesh: second-order HHO operators with the stabilization weight gamma, first-order system,
    // LTS fine region. The manufactured solution and the Dirichlet data follow m_t, which their functions capture by
    // reference: the object is neither copied nor moved.
    class problem {
    public:
        mesh_type msh;
        disk::hho_degree_info hho_di;
        domain_data dd;
        scal_vec_analytic_functions functions;
        RealType m_t = 0.0;
        std::function<disk::static_vector<RealType, 2>(const point_type&)>    u_fun, v_fun, a_fun, f_fun, s_flux_fun;
        std::function<disk::static_matrix<RealType, 2, 2>(const point_type&)> flux_fun;
        std::function<RealType(const point_type&)>                            s_u_fun, s_v_fun, s_a_fun, s_f_fun;
        sparse_type K, D1;          // stiffness (weight gamma), sigma* (unit weight)
        RealType gamma = 1.0;
        RealType omega_max = 0.0;   // largest frequency of the semi-discrete problem (problem_summary)
        size_t p = 1;               // fine substeps of the LTS integrator
        size_t n_fine_cells = 0;
        vector_type Py;             // fine/coarse projector on the state (0/1)

        // HHO setting of prm.space, stabilization weight gamma_in, algebraic faces from the splitting or direct
        problem(const wave_parameters& prm, size_t level, bool splitting, RealType gamma_in)
            : hho_di(prm.mixed_order() ? prm.space.k + 1 : prm.space.k, prm.space.k), gamma(gamma_in) {

            build_mesh(msh, prm.mesh, level);
            dd = classify_domain(msh, prm.problem.medium);

            functions.set_function_type(solution_type(prm.problem.solution));
            u_fun      = functions.Evaluate_u(m_t);
            v_fun      = functions.Evaluate_v(m_t);
            a_fun      = functions.Evaluate_a(m_t);
            f_fun      = functions.Evaluate_f(m_t);
            flux_fun   = functions.Evaluate_sigma(m_t);
            s_u_fun    = functions.Evaluate_s_u(m_t);
            s_v_fun    = functions.Evaluate_s_v(m_t);
            s_a_fun    = functions.Evaluate_s_a(m_t);
            s_f_fun    = functions.Evaluate_s_f(m_t);
            s_flux_fun = functions.Evaluate_s_q(m_t);

            // second-order formulation: Dirichlet data on the displacement and on the potential, those of the
            // solution or homogeneous (problem.dirichlet = "zero", the initial data remaining those of the solution)
            std::function<disk::static_vector<RealType, 2>(const point_type&)> e_data = u_fun;
            std::function<RealType(const point_type&)> a_data = s_u_fun;
            if (prm.problem.dirichlet == "zero") {
                e_data = [](const point_type&) -> disk::static_vector<RealType, 2> { return disk::static_vector<RealType, 2>::Zero(); };
                a_data = [](const point_type&) -> RealType { return 0.0; };
            }
            e_boundary_type e_bnd(msh);
            a_boundary_type a_bnd(msh);
            e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, e_data);
            a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, a_data);

            m_assembler = std::make_unique<assembler_type>(msh, hho_di, e_bnd, a_bnd, dd.e_material, dd.a_material);
            assembler_type& as = *m_assembler;
            as.set_interface_cell_indexes(dd.interface_cell_pair_indexes);
            if (prm.least_squares())
                as.set_least_squares_stabilization();
            else
                as.set_hho_stabilization();
            as.set_stabilization_scaling(prm.scaled_stabilization());
            as.set_stabilization_size(prm.stabilization_size());
            as.assemble_mass(msh);
            as.set_legacy_interface_quadrature(prm.space.interface_quadrature == "legacy");
            as.assemble_coupling_terms(msh);
            D1 = as.face_stabilization_mass(msh);

            as.set_stabilization_weight(gamma);
            as.assemble(msh, f_fun, s_f_fun);
            K = as.LHS;
            m_sys = std::make_unique<second_order_system>(msh, as, K, gamma * D1, splitting);
            m_sys->tol = prm.faces.tolerance;
            m_sys->max_iter = prm.faces.max_iterations;

            // LTS fine region (as in HHOFirstOrder): cells smaller than h_c and one vertex layer
            p = fine_substeps(prm, dd);
            RealType h_c = fine_size(prm, dd);
            std::vector<bool> e_cell_fine, a_cell_fine, e_face_fine, a_face_fine;
            as.mark_fine_region(msh, h_c, 1, e_cell_fine, a_cell_fine, e_face_fine, a_face_fine);
            second_order_system& s = *m_sys;
            vector_type Pc = vector_type::Zero(s.n_c);
            const size_t n_e_cbs = as.get_e_cell_basis_data(), n_a_cbs = as.get_a_cell_basis_data();
            for (size_t i = 0; i < e_cell_fine.size(); i++)
                if (e_cell_fine[i]) { Pc.segment(i * n_e_cbs, n_e_cbs).setOnes(); n_fine_cells++; }
            for (size_t i = 0; i < a_cell_fine.size(); i++)
                if (a_cell_fine[i]) { Pc.segment(as.get_e_n_cells_dof() + i * n_a_cbs, n_a_cbs).setOnes(); n_fine_cells++; }
            vector_type Pd = vector_type::Zero(s.n_d);
            std::vector<size_t> e_compress = as.get_e_compress(), a_compress = as.get_a_compress();
            size_t ofs = 0, n_d_face = s.n_d / std::max<size_t>(as.get_interfaces().size(), 1);
            for (auto& chunk : as.get_interfaces()) {
                if (e_face_fine[e_compress.at(chunk.first)] || a_face_fine[a_compress.at(chunk.first)])
                    Pd.segment(ofs, n_d_face).setOnes();
                ofs += n_d_face;
            }
            Py = vector_type(s.state_size());
            Py << Pc, Pc, Pd;
        }
        problem(const problem&) = delete;
        problem& operator=(const problem&) = delete;

        assembler_type&      assembler() { return *m_assembler; }
        second_order_system& system()    { return *m_sys; }

        // Critical weight gamma* of the splitting (Steins et al., M2AN 2023, Lemma 4.3 and Table 1), cell by cell without
        // boundary conditions: with the cell unknowns set to zero, the largest eigenvalue of B_T x = mu (2 sigma*_T - S_T) x
        // on the face unknowns of T (B: consistent part of the local stiffness, S: stabilization, sigma*: its face part,
        // unit weights; S_T = sigma*_T on the faces for the Lehrenfeld-Schoberl stabilization). gamma* = max over the
        // cells only depends on the cell shapes, not on h; D^{-1} (K_aa - D) has a spectral radius <= gamma* / gamma.
        RealType critical_weight() {
            assembler_type& as = assembler();
            RealType gamma_star = 0.0;
            auto local = [&](const dense_type& K0, const dense_type& K1, const dense_type& Dff) {
                size_t n_f = Dff.rows();
                dense_type B = K0.bottomRightCorner(n_f, n_f);
                dense_type A = 2.0 * Dff - (K1 - K0).bottomRightCorner(n_f, n_f);
                B = (0.5 * (B + B.transpose())).eval();
                A = (0.5 * (A + A.transpose())).eval();
                if (LLT<dense_type>(A).info() != Success)
                    throw std::runtime_error("critical_weight: 2 sigma*_T - S_T is not positive definite");
                GeneralizedSelfAdjointEigenSolver<dense_type> ges(B, A, EigenvaluesOnly);
                gamma_star = std::max(gamma_star, ges.eigenvalues().maxCoeff());
            };
            auto storage = msh.backend_storage();
            for (auto& e_chunk : as.get_e_material_data()) {
                auto& cell = storage->surfaces[e_chunk.first];
                as.set_stabilization_weight(0.0);
                dense_type K0 = as.e_laplacian_operator(e_chunk.second, msh, cell);
                as.set_stabilization_weight(1.0);
                dense_type K1 = as.e_laplacian_operator(e_chunk.second, msh, cell);
                local(K0, K1, as.e_face_stabilization_mass(e_chunk.second, msh, cell));
            }
            for (auto& a_chunk : as.get_a_material_data()) {
                auto& cell = storage->surfaces[a_chunk.first];
                as.set_stabilization_weight(0.0);
                dense_type K0 = as.a_laplacian_operator(a_chunk.second, msh, cell);
                as.set_stabilization_weight(1.0);
                dense_type K1 = as.a_laplacian_operator(a_chunk.second, msh, cell);
                local(K0, K1, as.a_face_stabilization_mass(a_chunk.second, msh, cell));
            }
            as.set_stabilization_weight(gamma);
            return gamma_star;
        }

        // Q F(t): sources and Dirichlet contributions at time t, in the reordered unknowns
        vector_type eval_F(RealType t) {
            return system().Q * eval_F_unknowns(t);
        }

        // F(t) in the assembler numbering
        vector_type eval_F_unknowns(RealType t) {
            m_t = t;
            assembler().assemble_rhs(msh, f_fun, s_f_fun);
            return assembler().RHS;
        }

        // Newmark initial data in the assembler numbering: projections of the exact displacement, velocity and
        // acceleration at t = 0 on the cells and the faces (as in IHHOSecondOrder)
        std::array<vector_type, 3> newmark_initial_state() {
            m_t = 0.0;
            std::array<vector_type, 3> x;
            assembler().project_over_cells(msh, x[0], u_fun, s_u_fun);
            assembler().project_over_faces(msh, x[0], u_fun, s_u_fun);
            assembler().project_over_cells(msh, x[1], v_fun, s_v_fun);
            assembler().project_over_faces(msh, x[1], v_fun, s_v_fun);
            assembler().project_over_cells(msh, x[2], a_fun, s_a_fun);
            assembler().project_over_faces(msh, x[2], a_fun, s_a_fun);
            return x;
        }

        // State y = (u_T, v, d) of the unknowns u (cells and faces) and of the cell velocities of v
        vector_type state(const vector_type& u, const vector_type& v) {
            const size_t n_c = system().n_c, n_d = system().n_d;
            vector_type uw = system().Q * u;
            vector_type y(system().state_size());
            y << uw.head(n_c), v.head(n_c), uw.segment(n_c, n_d);
            return y;
        }

        // Projection of the exact solution at t = 0: cells (u_T, v), interface faces d
        vector_type initial_state() {
            m_t = 0.0;
            vector_type u, v;
            assembler().project_over_cells(msh, u, u_fun, s_u_fun);
            assembler().project_over_faces(msh, u, u_fun, s_u_fun);
            assembler().project_over_cells(msh, v, v_fun, s_v_fun);
            vector_type uw = system().Q * u;
            const size_t n_c = system().n_c, n_d = system().n_d;
            vector_type y(system().state_size());
            y << uw.head(n_c), v.head(n_c), uw.segment(n_c, n_d);
            return y;
        }

        // {h_min, u, sigma, phi, grad phi, v = u_t, phi_t} errors at time t
        std::array<RealType, 7> errors(const vector_type& y, RealType t) {
            vector_type x = system().unknowns(eval_F(t), y);
            m_t = t;
            auto e = postprocessor<mesh_type>::errors_two_fields_elastoacoustic(msh, hho_di, assembler(), x, u_fun, flux_fun, s_u_fun, s_flux_fun);
            vector_type xv = vector_type::Zero(x.size());
            xv.head(system().n_c) = y.segment(system().n_c, system().n_c);
            auto ev = postprocessor<mesh_type>::errors_two_fields_elastoacoustic(msh, hho_di, assembler(), xv, v_fun, flux_fun, s_v_fun, s_flux_fun);
            return {e[0], e[1], e[2], e[3], e[4], ev[1], ev[3]};
        }

        // Mass norm sqrt(du^T M du + dv^T M dv) (cell unknowns) of the difference of two states
        RealType mass_distance(const vector_type& y1, const vector_type& y2) {
            const size_t n_c = system().n_c;
            vector_type du = y1.head(n_c) - y2.head(n_c), dv = y1.segment(n_c, n_c) - y2.segment(n_c, n_c);
            return std::sqrt(std::abs(du.dot(system().Mcc * du) + dv.dot(system().Mcc * dv)));
        }

        // Energy norm sqrt(dv^T M dv + dx^T K dx) of the difference of two states at time t
        RealType energy_distance(const vector_type& y1, const vector_type& y2, RealType t) {
            vector_type F = eval_F(t);
            vector_type dx = system().unknowns(F, y1) - system().unknowns(F, y2);
            const size_t n_c = system().n_c;
            vector_type dv = y1.segment(n_c, n_c) - y2.segment(n_c, n_c);
            return std::sqrt(std::abs(dv.dot(system().Mcc * dv) + dx.dot(K * dx)));
        }

        // Discrete energy (69) of Burman-Duran-Ern-Steins (JSC 2021) with delta = 0 at time t of the unknowns u, v
        // (assembler numbering): 1/2 v^T M v + 1/2 sum_T u_T^T A_T u_T, the local vectors u_T being completed on the
        // Dirichlet faces by the projection of the data at time t (gather_dof_data, as in the original code of the
        // paper). u^T K u is the same sum with zero Dirichlet values: only the cells with a boundary face are
        // corrected, with their local operators computed once. Equal to 1/2 v^T M v + 1/2 u^T K u for homogeneous data.
        RealType energy(const vector_type& u, const vector_type& v, RealType t) {
            assembler_type& as = assembler();
            auto storage = msh.backend_storage();
            if (!m_boundary_ready) {
                auto on_boundary = [&](const auto& cell) {
                    for (auto& fc : faces(msh, cell))
                        if (msh.is_boundary(fc)) return true;
                    return false;
                };
                for (auto& e_chunk : as.get_e_material_data()) {
                    auto& cell = storage->surfaces[e_chunk.first];
                    if (on_boundary(cell))
                        m_boundary_cells.push_back({e_chunk.first, true, as.e_laplacian_operator(e_chunk.second, msh, cell)});
                }
                for (auto& a_chunk : as.get_a_material_data()) {
                    auto& cell = storage->surfaces[a_chunk.first];
                    if (on_boundary(cell))
                        m_boundary_cells.push_back({a_chunk.first, false, as.a_laplacian_operator(a_chunk.second, msh, cell)});
                }
                m_boundary_ready = true;
            }
            m_t = t;
            const vector_type zero = vector_type::Zero(u.size());
            RealType E = v.dot(as.MASS * v) + u.dot(K * u);
            for (auto& b : m_boundary_cells) {
                auto& cell = storage->surfaces[b.cell];
                vector_type uT = b.elastic ? as.gather_e_dof_data(b.cell, msh, cell, u) : as.gather_a_dof_data(b.cell, msh, cell, u);
                vector_type gT = b.elastic ? as.gather_e_dof_data(b.cell, msh, cell, zero) : as.gather_a_dof_data(b.cell, msh, cell, zero);
                vector_type u0 = uT - gT;
                E += uT.dot(b.A * uT) - u0.dot(b.A * u0);
            }
            return 0.5 * E;
        }

    private:
        struct boundary_cell { size_t cell; bool elastic; dense_type A; };
        std::unique_ptr<assembler_type>      m_assembler;
        std::unique_ptr<second_order_system> m_sys;
        std::vector<boundary_cell>           m_boundary_cells;
        bool                                 m_boundary_ready = false;
    };

    // ERK(s)
    inline vector_type integrate_erk(problem& pb, vector_type y, size_t nt, RealType tf, size_t s) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        second_order_system& sys = pb.system();
        RealType dt = tf / nt;
        std::vector<vector_type> k(s);
        vector_type a_ws;
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            vector_type y_n = y;
            for (size_t i = 0; i < s; i++) {
                vector_type Y = y_n;
                for (size_t j = 0; j < i; j++)
                    Y += a(i,j) * dt * k[j];
                k[i] = sys.eval(pb.eval_F(tn + c(i) * dt), Y, a_ws);
                y += dt * b(i) * k[i];
            }
        }
        return y;
    }

    // SDIRK(s, s+1): stage Y_i = R_i + a_ii dt f(t_i, Y_i), one factorization since a_ii is constant
    inline vector_type integrate_sdirk(problem& pb, vector_type y, size_t nt, RealType tf, size_t s) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(int(s), a, b, c);
        second_order_system& sys = pb.system();
        RealType dt = tf / nt;
        RealType h = a(0,0) * dt;
        sys.factorize_stage(h, vector_type::Ones(sys.state_size()));
        std::vector<vector_type> k(s);
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            vector_type y_n = y;
            for (size_t i = 0; i < s; i++) {
                vector_type R = y_n;
                for (size_t j = 0; j < i; j++)
                    R += a(i,j) * dt * k[j];
                vector_type Y = sys.stage_solve(pb.eval_F(tn + c(i) * dt), R);
                k[i] = (Y - R) / h;
                y += dt * b(i) * k[i];
            }
        }
        return y;
    }

    // Coarse Taylor coefficients of the LTS integrators (two-level Algorithm 3 of the first-order drivers) for
    // y' = B y + g(t), g quadratic on [t_n, t_n + dt], Q = I - P:
    //   w0 = B Q y + Q g0,  w1 = B Q (B y + g0) + Q g1,  w2 = B Q (B^2 y + B g0 + g1) + Q g2,
    //   w3 = B Q (B^3 y + B^2 g0 + B g1 + g2),   g(t_n + tau) = g0 + tau g1 + tau^2/2 g2.
    // ws: warm starts of the algebraic faces of the 13 evaluations.
    struct lts_coarse {
        std::vector<vector_type> w = std::vector<vector_type>(4);
        std::vector<vector_type> ws = std::vector<vector_type>(13);

        void compute(problem& pb, const vector_type& y, RealType tn, RealType dt) {
            second_order_system& sys = pb.system();
            vector_type zero = vector_type::Zero(sys.N);
            vector_type ys = vector_type::Zero(sys.state_size());
            vector_type Qp = vector_type::Ones(sys.state_size()) - pb.Py;
            auto B = [&](const vector_type& z, size_t i) { return sys.eval(zero, z, ws[i]); };
            vector_type g0  = sys.eval(pb.eval_F(tn), ys, ws[0]);
            vector_type g12 = sys.eval(pb.eval_F(tn + 0.5 * dt), ys, ws[1]);
            vector_type g1e = sys.eval(pb.eval_F(tn + dt), ys, ws[2]);
            vector_type g1 = (-3.0 * g0 + 4.0 * g12 - g1e) / dt;
            vector_type g2 = (4.0 * g0 - 8.0 * g12 + 4.0 * g1e) / (dt * dt);
            vector_type By = B(y, 3), B2y = B(By, 4), B3y = B(B2y, 5);
            vector_type Bg0 = B(g0, 6), B2g0 = B(Bg0, 7), Bg1 = B(g1, 8);
            w[0] = B(Qp.cwiseProduct(y), 9) + Qp.cwiseProduct(g0);
            w[1] = B(Qp.cwiseProduct(By + g0), 10) + Qp.cwiseProduct(g1);
            w[2] = B(Qp.cwiseProduct(B2y + Bg0 + g1), 11) + Qp.cwiseProduct(g2);
            w[3] = B(Qp.cwiseProduct(B3y + B2g0 + Bg1 + g2), 12);
        }

        vector_type taylor(RealType tau) const {
            return w[0] + tau * w[1] + (tau * tau / 2.0) * w[2] + (tau * tau * tau / 6.0) * w[3];
        }
    };

    // Two-level LTS-RK4: p fine RK4 substeps of z' = B P z + P g(t) + w(tau)
    inline vector_type integrate_lts(problem& pb, vector_type y, size_t nt, RealType tf) {
        second_order_system& sys = pb.system();
        RealType dt = tf / nt;
        size_t p = pb.p;
        RealType dtau = dt / p;
        bool all_fine = pb.Py.minCoeff() > 0.5;
        lts_coarse coarse;
        vector_type ys = vector_type::Zero(sys.state_size());
        std::vector<vector_type> ws(2);
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            if (all_fine)
                for (auto& wi : coarse.w) wi = vector_type::Zero(sys.state_size());
            else
                coarse.compute(pb, y, tn, dt);
            // k = B P z + P g(t) + w(tau) = f(F(t), P z) - (I - P) g(t) + w(tau), g(t) = f(F(t), 0)
            auto fine_stage = [&](const vector_type& z, RealType tau) {
                vector_type F = pb.eval_F(tn + tau);
                vector_type k = sys.eval(F, pb.Py.cwiseProduct(z), ws[0]);
                if (!all_fine)
                    k -= (vector_type::Ones(sys.state_size()) - pb.Py).cwiseProduct(sys.eval(F, ys, ws[1]));
                return vector_type(k + coarse.taylor(tau));
            };
            for (size_t m = 0; m < p; m++) {
                RealType tm = m * dtau;
                vector_type k1 = fine_stage(y, tm);
                vector_type k2 = fine_stage(y + 0.5 * dtau * k1, tm + 0.5 * dtau);
                vector_type k3 = fine_stage(y + 0.5 * dtau * k2, tm + 0.5 * dtau);
                vector_type k4 = fine_stage(y + dtau * k3, tm + dtau);
                y += dtau * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            }
        }
        return y;
    }

    // LI-SDIRK(4,3): coarse ERK4 Taylor coefficients, one SDIRK(3,4) step of z' = B P z + P g(t) + w(tau) whose
    // stages (I - a_ii dt B P) Z = R + a_ii dt (P g + w) are solved with the projected stage matrix
    inline vector_type integrate_li(problem& pb, vector_type y, size_t nt, RealType tf) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(3, a, b, c);
        second_order_system& sys = pb.system();
        RealType dt = tf / nt;
        RealType h = a(0,0) * dt;
        sys.factorize_stage(h, pb.Py);
        lts_coarse coarse;
        vector_type ys = vector_type::Zero(sys.state_size()), zero = vector_type::Zero(sys.N), ws;
        std::vector<vector_type> k(3);
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            coarse.compute(pb, y, tn, dt);
            vector_type y_n = y;
            for (size_t i = 0; i < 3; i++) {
                vector_type R = y_n;
                for (size_t j = 0; j < i; j++)
                    R += a(i,j) * dt * k[j];
                vector_type g = sys.eval(pb.eval_F(tn + c(i) * dt), ys, ws);
                vector_type rhs = R + h * (pb.Py.cwiseProduct(g) + coarse.taylor(c(i) * dt));
                vector_type Z = sys.stage_solve(zero, rhs);
                k[i] = (Z - R) / h;
                y += dt * b(i) * k[i];
            }
        }
        return y;
    }

    // Newmark(beta, gamma) on M u'' + C u' + K u = F (cells and faces, M zero on the faces): predictor u*, v*, solve
    // (M + gamma dt C + beta dt^2 K) a = F(t_{n+1}) - K u* - C v*, corrector. Discrete energy (69) of [JSC21]
    // (problem::energy, Dirichlet faces = data) written to energy_os at each step if given (exactly conserved for
    // beta = 1/4, gamma = 1/2 without sources). Returns the state y = (u_T, v_T, d) at the final time.
    inline vector_type integrate_newmark(problem& pb, const wave_parameters& prm, size_t nt, std::ostream* energy_os = nullptr) {
        assembler_type& as = pb.assembler();
        const RealType tf = prm.time.final_time, dt = tf / nt;
        const RealType beta = prm.time.newmark_beta, gam = prm.time.newmark_gamma;
        const sparse_type& K = pb.K;
        const sparse_type& M = as.MASS;
        const sparse_type& C = as.COUPLING;
        sparse_type A = M + gam * dt * C + beta * dt * dt * K;
        linear_solver<RealType> analysis;
        if (prm.solver.static_condensation) {
            analysis.set_Kg(A, as.get_n_face_dof());
            std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(2);
            vec_cell_basis_data[0] = std::make_pair(as.get_e_material_data().size(), as.get_e_cell_basis_data());
            vec_cell_basis_data[1] = std::make_pair(as.get_a_material_data().size(), as.get_a_cell_basis_data());
            analysis.condense_equations(vec_cell_basis_data);
        }
        else
            analysis.set_Kg(A);
        analysis.factorize();

        auto x0 = pb.newmark_initial_state();
        vector_type u = x0[0], v = x0[1], a = x0[2];
        if (energy_os)
            *energy_os << std::setprecision(16) << std::scientific << 0.0 << " " << pb.energy(u, v, 0.0) << std::endl;
        for (size_t it = 0; it < nt; it++) {
            vector_type F = pb.eval_F_unknowns(dt * (it + 1));
            u += dt * v + 0.5 * dt * dt * (1.0 - 2.0 * beta) * a;
            v += dt * (1.0 - gam) * a;
            vector_type rhs = F - K * u - C * v;
            a = analysis.solve(rhs);
            u += beta * dt * dt * a;
            v += gam * dt * a;
            if (energy_os)
                *energy_os << dt * (it + 1) << " " << pb.energy(u, v, dt * (it + 1)) << std::endl;
        }
        return pb.state(u, v);
    }

    inline vector_type integrate(problem& pb, const wave_parameters& prm, vector_type y, size_t nt) {
        RealType tf = prm.time.final_time;
        if (prm.time.scheme == "newmark") return integrate_newmark(pb, prm, nt);
        if (prm.time.scheme == "erk")   return integrate_erk(pb, y, nt, tf, prm.stages());
        if (prm.time.scheme == "sdirk") return integrate_sdirk(pb, y, nt, tf, prm.stages());
        if (prm.time.scheme == "lts")   return integrate_lts(pb, y, nt, tf);
        return integrate_li(pb, y, nt, tf);
    }

    // Mesh, discretization and splitting data; omega_max (explicit integrators) and the spectral radius of the splitting
    // come from power iterations
    inline std::string problem_summary(problem& pb, bool explicit_scheme) {
        second_order_system& sys = pb.system();
        std::ostringstream os;
        os << "      h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min << "  gamma = " << pb.gamma;
        os << "  unknowns: cells " << sys.n_c << ", dynamic faces " << sys.n_d << ", algebraic faces " << sys.n_a
           << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << std::endl;
        if (explicit_scheme) {
            RealType omega = sys.max_frequency();
            pb.omega_max = omega;
            os << "      omega_max ~ " << omega << " (ERK4 stable for dt < ~" << 2.8 / omega << ")  ";
        }
        else
            os << "      ";
        if (sys.splitting_Q)
            os << "splitting spectral radius = " << sys.splitting_spectral_radius() << "  ";
        os << "[checks: |C off d|/|C| = " << sys.coupling_leak << ", |sigma* off blocks|/|sigma*| = " << sys.sigma_off_block << "]" << std::endl;
        sys.n_solves = sys.n_iterations = sys.n_failures = sys.max_iterations = 0;
        return os.str();
    }

    // Explicit integrators: warning when dt omega_max exceeds the imaginary-axis stability bound (ERK(3): sqrt(3),
    // ERK(4) and the fine substeps of LTS-RK4: 2 sqrt(2); ERK(1), ERK(2) are unstable for any dt on this axis and only
    // converge for small dt omega_max)
    inline void check_time_step(const problem& pb, const wave_parameters& prm, RealType dt) {
        RealType bound = 0.0;
        const std::string& scheme = prm.time.scheme;
        if (scheme == "erk" && prm.stages() == 3) bound = std::sqrt(3.0);
        if ((scheme == "erk" && prm.stages() == 4) || scheme == "lts") bound = 2.0 * std::sqrt(2.0);
        RealType z = dt * pb.omega_max / (scheme == "lts" ? RealType(pb.p) : 1.0);
        if (bound > 0.0 && z > bound)
            std::cout << bold << red << "      WARNING: dt omega_max = " << z << " > " << bound
                      << " (stability limit of the explicit integrator): increase time.n_steps" << reset << std::endl;
    }

    inline void print_splitting_stats(second_order_system& sys, std::ostream& os) {
        if (!sys.splitting_Q || sys.n_solves == 0) return;
        os << "      splitting: " << sys.n_solves << " face solves, " << RealType(sys.n_iterations) / sys.n_solves
           << " iterations on average, max " << sys.max_iterations << ", not converged " << sys.n_failures << std::endl;
        sys.n_solves = sys.n_iterations = sys.n_failures = sys.max_iterations = 0;
    }

}

void HHOSecondOrder(const wave_parameters& prm);

void HHOSecondOrder(const wave_parameters& prm) {

    using namespace hho_second_order;

    // algebraic faces of the explicit (erk, lts) and locally implicit (li) integrators: splitting or direct; sdirk
    // solves them within its stage systems
    bool splitting = prm.explicit_scheme() && prm.faces.solver == "splitting";
    std::cout << std::endl << bold << red << "   SECOND-ORDER ELASTO-ACOUSTIC COUPLING, " << prm.scheme_label()
              << ", " << prm.space_label() << (splitting ? ", faces: Steins splitting" : ", faces: direct") << reset << std::endl;
    std::string fname = prm.output_file();
    std::ofstream log(fname);
    log << "second-order formulation" << std::endl;
    prm.print(log);

    // stabilization weight: space.weight, or weight_factor gamma* (gamma* computed cell by cell on the first mesh of
    // the study, it only depends on the cell shapes)
    RealType gamma = prm.space.weight;
    if (prm.space.weight_auto) {
        problem probe(prm, prm.mesh.level, false, 1.0);
        RealType gamma_star = probe.critical_weight();
        gamma = prm.space.weight_factor * gamma_star;
        std::ostringstream os;
        os << "      gamma* = " << gamma_star << " (cells of the mesh l = " << prm.mesh.level << "), gamma = "
           << prm.space.weight_factor << " gamma* = " << gamma << std::endl;
        std::cout << os.str();
        log << os.str();
    }

    // errors at the final time: elastic displacement, stress (reconstruction) and velocity, acoustic potential, its
    // reconstructed gradient and its time derivative, L2 norms
    error_columns all = columns(prm, {"elastic_u", "elastic_sigma", "elastic_v"}, {"acoustic_phi", "acoustic_grad_phi", "acoustic_phi_t"});
    auto all_errors = [](const std::array<RealType, 7>& e) { return std::vector<RealType>{e[1], e[2], e[5], e[3], e[4], e[6]}; };
    const bool explicit_scheme = prm.explicit_scheme();

    // Newmark energy history of the last mesh of a space study / of the finest run of a time study
    std::ofstream energy_os;
    if (prm.output.energy && prm.time.scheme == "newmark")
        energy_os.open(prm.energy_file());
    auto run = [&](problem& pb, size_t nt, bool record) {
        if (prm.time.scheme == "newmark" && record && energy_os.is_open()) {
            energy_os.seekp(0);
            std::ostringstream hist;
            vector_type y = integrate_newmark(pb, prm, nt, &hist);
            // relative energy variation max_n |E^n - E^0| / E^0
            std::istringstream is(hist.str());
            RealType t, E, E0 = -1.0, var = 0.0;
            while (is >> t >> E) {
                if (E0 < 0.0) E0 = E;
                var = std::max(var, std::abs(E - E0) / E0);
            }
            energy_os << "# t E (discrete energy (69) of JSC21, Dirichlet faces = data), max |E - E0| / E0 = " << var << std::endl << hist.str();
            std::cout << bold << cyan << "      Newmark energy: max |E^n - E^0| / E^0 = " << var << reset << std::endl;
            return y;
        }
        return integrate(pb, prm, pb.initial_state(), nt);
    };

    const RealType tf = prm.time.final_time;
    std::vector<conv_row> rows;
    std::vector<std::string> names;
    if (prm.study.type != "time") {
        names = all.names;
        std::vector<size_t> levels = study_levels(prm);
        for (size_t l : levels) {
            problem pb(prm, l, splitting, gamma);
            size_t nt = study_steps(prm, l, pb.dd.h_min);
            std::cout << bold << cyan << "      l = " << l << "  nt = " << nt << reset << std::endl;
            std::string summary = problem_summary(pb, explicit_scheme);
            std::cout << summary;
            log << summary;
            if (explicit_scheme) check_time_step(pb, prm, tf / nt);
            vector_type y = run(pb, nt, l == levels.back());
            print_splitting_stats(pb.system(), std::cout);
            rows.push_back({pb.dd.h_max, tf / nt, all.select(all_errors(pb.errors(y, tf))), pb.dd.h_min});
            print_table(names, rows, false, std::cout);
        }
        print_table(names, rows, false, log);
    }
    else {
        // against the exact solution only, or against a reference solution of the same integrator (mass and energy
        // norms of the difference) with the exact errors of the displacement / potential
        bool reference = prm.study.reference == "numerical";
        error_columns cols = reference ? columns(prm, {"elastic_u"}, {"acoustic_phi"}) : all;
        names = cols.names;
        if (reference) names.insert(names.begin(), {"mass_norm", "energy_norm"});
        problem pb(prm, prm.mesh.level, splitting, gamma);
        std::string summary = problem_summary(pb, explicit_scheme);
        std::cout << summary;
        log << summary;
        size_t n0 = prm.time.n_steps;
        size_t n_levels = prm.study.time_levels;
        if (explicit_scheme) check_time_step(pb, prm, tf / n0);
        vector_type y_ref;
        if (reference) {
            size_t nt_ref = prm.study.reference_factor * n0 * (size_t(1) << (n_levels - 1));
            y_ref = run(pb, nt_ref, false);
            std::cout << bold << cyan << "      reference nt = " << nt_ref << reset << std::endl;
            print_splitting_stats(pb.system(), std::cout);
        }
        for (size_t i = 0; i < n_levels; i++) {
            size_t nt = n0 * (size_t(1) << i);
            vector_type y = run(pb, nt, i + 1 == n_levels);
            auto e = pb.errors(y, tf);
            std::vector<RealType> errors = reference ? cols.select({e[1], e[3]}) : cols.select(all_errors(e));
            if (reference) errors.insert(errors.begin(), {pb.mass_distance(y, y_ref), pb.energy_distance(y, y_ref, tf)});
            rows.push_back({pb.dd.h_max, tf / nt, errors, pb.dd.h_min});
            print_table(names, rows, true, std::cout);
        }
        print_splitting_stats(pb.system(), std::cout);
        print_table(names, rows, true, log);
    }
    write_table(prm.table_file(), names, rows);
    std::cout << bold << red << "   --> " << fname << ", " << prm.table_file() << reset << std::endl << std::endl;
}

#endif /* HHOSecondOrder_hpp */
