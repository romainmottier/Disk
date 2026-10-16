//
//  problem.hpp
//
//  Second-order formulation (HHOSecondOrder.hpp): the coupled problem on one mesh (class problem: assembly,
//  stabilization weight, initial state, F(t), energy, errors).
//

#pragma once
#ifndef second_order_problem_hpp
#define second_order_problem_hpp

#include "system.hpp"

namespace hho_second_order {

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
            as.set_parallel_assembly(prm.solver.fast);     // fast mode: local operators once per geometric class
            as.set_class_reuse(prm.solver.fast);           //   of cells (same matrices up to round-off)
            as.assemble_mass(msh);
            as.set_legacy_interface_quadrature(prm.space.interface_quadrature == "legacy");
            as.assemble_coupling_terms(msh);
            D1 = as.face_stabilization_mass(msh);

            as.set_stabilization_weight(gamma);
            as.set_cache_rhs_data(prm.solver.fast);        // F(t) at each stage: time-independent data kept
            as.assemble(msh, f_fun, s_f_fun);
            K = as.LHS;
            m_sys = std::make_unique<second_order_system>(msh, as, K, gamma * D1, splitting);
            m_sys->tol = prm.faces.tolerance;
            m_sys->max_iter = prm.faces.max_iterations;
            m_sys->set_fast(prm.solver.fast);

            // LTS fine region (as in HHOFirstOrder): cells smaller than h_c and one vertex layer
            p = fine_substeps(prm, dd);
            RealType h_c = fine_size(prm, dd);
            std::vector<bool> e_cell_fine, a_cell_fine, e_face_fine, a_face_fine;
            as.mark_fine_region(msh, h_c, prm.time.fine_layers, e_cell_fine, a_cell_fine, e_face_fine, a_face_fine);
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

        // Multilevel LTS (time.scheme = "mlts"): level of each unknown of the state y = (u_T, v_T, d) from the levels of
        // the cells (coupling_common::mlts_cell_levels): the cell unknowns at the level of their cell, the interface
        // unknowns at the level of the finer of their two cells. Returns the number of levels.
        size_t mlts_state_levels(const wave_parameters& prm, std::vector<int>& level) {
            std::vector<int> cl;
            const size_t L = coupling_common::mlts_cell_levels(msh, dd, prm, cl);
            assembler_type& as = assembler();
            second_order_system& s = system();
            const size_t n_e_cbs = as.get_e_cell_basis_data(), n_a_cbs = as.get_a_cell_basis_data();
            std::vector<int> lc(s.n_c, 0), ld(s.n_d, 0);
            for (auto& ci : as.get_e_cell_index())
                for (size_t k = 0; k < n_e_cbs; k++) lc[ci.second * n_e_cbs + k] = cl[ci.first];
            for (auto& ci : as.get_a_cell_index())
                for (size_t k = 0; k < n_a_cbs; k++) lc[as.get_e_n_cells_dof() + ci.second * n_a_cbs + k] = cl[ci.first];
            size_t ofs = 0, n_d_face = s.n_d / std::max<size_t>(as.get_interfaces().size(), 1);
            for (auto& chunk : as.get_interfaces()) {
                const int l = std::max(cl[chunk.second.first], cl[chunk.second.second]);
                for (size_t k = 0; k < n_d_face; k++) ld[ofs + k] = l;
                ofs += n_d_face;
            }
            level.clear();
            level.insert(level.end(), lc.begin(), lc.end());
            level.insert(level.end(), lc.begin(), lc.end());
            level.insert(level.end(), ld.begin(), ld.end());
            return L;
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

        // Random state (uniform in [-1, 1], fixed seed): initial data of the CFL study exciting all the modes
        // (study.cfl_random)
        vector_type random_state() {
            vector_type y = initial_state();
            std::mt19937 gen(12345);
            std::uniform_real_distribution<RealType> u(-1.0, 1.0);
            for (Eigen::Index i = 0; i < y.size(); i++) y(i) = u(gen);
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

}

#endif /* second_order_problem_hpp */
