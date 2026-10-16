//
//  problem.hpp
//
//  First-order formulation (HHOFirstOrder.hpp): the coupled problem on one mesh (class problem: assembly of M, K and
//  F(t), initial state, energy, errors, ERK and matrix-free operators, LTS projectors) and the elastic problem of the
//  three-field formulation (class problem_three_fields).
//

#pragma once
#ifndef first_order_problem_hpp
#define first_order_problem_hpp

#include "../mesh/domain.hpp"
#include "../physics/manufactured_solutions.hpp"
#include "../studies/convergence.hpp"
#include "../assembly/elastodynamic_three_fields_assembler.hpp"
#include <chrono>
#include "../time_schemes/dirk_coupling_hho_scheme.hpp"
#include "../time_schemes/mf_erk_operator.hpp"

namespace hho_first_order {

    using namespace coupling_common;
    typedef elastoacoustic_four_fields_assembler<mesh_type>  assembler_type;
    typedef erk_coupling_hho_scheme<RealType>               erk_type;
    typedef mf_erk_operator<RealType>                       mf_type;

    // Coupled problem on one mesh: assembled M, K (stiffness + coupling), F(t), the ERK operator with its
    // block-diagonal inverses and the LTS fine/coarse projectors. The manufactured solution and the Dirichlet data
    // follow m_t, which their functions capture by reference: the object is neither copied nor moved.
    // zero_data (cfl, spectrum studies), problem.solution = "ricker", problem.source / dirichlet = "zero": no source
    // and / or homogeneous Dirichlet data; without both, F does not depend on time and is assembled once.
    class problem {
    public:
        mesh_type msh;
        disk::hho_degree_info hho_di;
        domain_data dd;
        std::vector<int> tags;      // cell tags (mesh.type = "polygons")
        scal_vec_analytic_functions functions;
        RealType m_t = 0.0;
        std::function<disk::static_vector<RealType, 2>(const point_type&)>    v_fun, f_fun, s_flux_fun, e_src;
        std::function<disk::static_matrix<RealType, 2, 2>(const point_type&)> flux_fun;
        std::function<RealType(const point_type&)>                            s_v_fun, s_f_fun, a_src;
        size_t p = 1;               // fine substeps of the LTS integrator
        RealType h_c = 0.0;         // cells smaller than h_c are fine
        size_t n_fine_cells = 0;
        bool ls_Q = true;           // least-squares stabilization: block-diagonal face block
        bool ricker_Q = false;      // initial Ricker wavelet (problem.ricker), no source, homogeneous data
        bool constant_F = false;    // no source and homogeneous Dirichlet data: F assembled once
        // F(t) = sum_j F_coeffs(t)[j] F_basis[j] (solver.matrix_free and problem.solution = "lshape_corner": the
        // Dirichlet part and the source part, separable in time), set by the first eval_F; eval_F_local evaluates F
        // on a few unknowns only (the matrix-free MLTS)
        std::vector<vector_type> F_basis;
        std::function<std::vector<RealType>(RealType)> F_coeffs;
        bool separable_F = false;
        bool lshape_sin = false;    // lshape_corner_sin: F(t) = w cos(w t) F_D - w^2 sin(w t) F_S
        bool consistent_faces = true;   // initial faces from the cells (K_ff x_F = -K_FT x_T), else projected
        bool m_matrix_free = false;     // solver.matrix_free: matrix-free explicit operator (mf()), never the ERK operator
        RealType m_mf_tolerance = 1e-10;
        int m_mf_batch = 64;
        bool m_mf_single = false;       // solver.matrix_free_precision = "single"
        bool m_mf_li_constraint = false;    // time.li_system = "constraint" (li)
        bool m_mf_li_stage = false;         // time.li_coupling = "stage" (li)
        bool m_mf_li_condensation = false;  // time.li_condensation (li)
        bool m_mf_mlts_reuse = false;       // time.mlts_reuse (mlts)
        bool m_mf_profile = false;          // output.profile (lts, mlts)
        size_t m_mf_par_threshold = 20000;  // solver.parallel_threshold
        wave_parameters::problem_data::ricker_data ricker;

        problem(const wave_parameters& prm, size_t level, bool zero_data = false)
            : hho_di(prm.mixed_order() ? prm.space.k + 1 : prm.space.k, prm.space.k), ls_Q(prm.least_squares()),
              m_matrix_free(prm.solver.matrix_free && prm.study.type != "spectrum"), m_mf_tolerance(prm.solver.matrix_free_tolerance),
              m_mf_batch(int(prm.solver.matrix_free_batch)), m_mf_single(prm.solver.matrix_free_precision == "single"),
              m_mf_li_constraint(prm.time.li_system == "constraint"),
              m_mf_li_stage(prm.time.li_coupling == "stage"), m_mf_li_condensation(prm.time.li_condensation),
              m_mf_mlts_reuse(prm.time.mlts_reuse), m_mf_profile(prm.output.profile),
              m_mf_par_threshold(prm.solver.parallel_threshold) {

            build_mesh(msh, prm.mesh, level, &tags);
            dd = classify_domain(msh, prm, tags);
            ricker_Q = prm.problem.solution == "ricker";
            ricker = prm.problem.ricker;
            const bool zero_bc  = zero_data || ricker_Q || prm.problem.dirichlet == "zero";
            const bool zero_src = zero_data || ricker_Q || prm.problem.source == "zero";
            constant_F = zero_bc && zero_src;
            // as the original Ricker drivers (SDIRK): initial faces = projections (zero)
            consistent_faces = !(ricker_Q && prm.time.scheme == "sdirk");

            functions.set_function_type(solution_type(prm.problem.solution));
            functions.set_paper_frequencies(prm.problem.omega, prm.problem.theta);
            scal_vec_analytic_functions::lshape_cx = prm.problem.lshape_corner[0];
            scal_vec_analytic_functions::lshape_cy = prm.problem.lshape_corner[1];
            scal_vec_analytic_functions::lshape_nu = prm.problem.lshape_exponent;
            v_fun      = functions.Evaluate_v(m_t);
            f_fun      = functions.Evaluate_f(m_t);
            flux_fun   = functions.Evaluate_sigma(m_t);
            s_v_fun    = functions.Evaluate_s_v(m_t);
            s_f_fun    = functions.Evaluate_s_f(m_t);
            s_flux_fun = functions.Evaluate_s_q(m_t);
            auto zero_v = [](const point_type&) -> disk::static_vector<RealType, 2> { return disk::static_vector<RealType, 2>::Zero(); };
            auto zero_s = [](const point_type&) -> RealType { return 0.0; };
            e_src = zero_src ? std::function<disk::static_vector<RealType, 2>(const point_type&)>(zero_v) : f_fun;
            a_src = zero_src ? std::function<RealType(const point_type&)>(zero_s) : s_f_fun;

            // first-order formulation: Dirichlet data on the elastic velocity and on the acoustic scalar unknown
            e_boundary_type e_bnd(msh);
            a_boundary_type a_bnd(msh);
            if (zero_bc) {
                e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, zero_v);
                a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, zero_s);
            }
            else {
                e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, v_fun);
                a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_v_fun);
            }

            m_assembler = std::make_unique<assembler_type>(msh, hho_di, e_bnd, a_bnd, dd.e_material, dd.a_material);
            assembler_type& as = *m_assembler;
            as.set_interface_cell_indexes(dd.interface_cell_pair_indexes);
            if (ls_Q)
                as.set_least_squares_stabilization();
            else
                as.set_hho_stabilization();
            if (prm.scaled_stabilization())
                as.set_scaled_stabilization();
            if (prm.space.scaling == "graded")
                as.set_stabilization_graded(dd.h_max);
            if (prm.space.scaling == "hybrid")
                as.set_stabilization_hybrid(prm.time.fine_size > 0.0 ? RealType(prm.time.fine_size) : RealType(0.75) * dd.h_max);   // the implicit region of li
            as.set_stabilization_size(prm.stabilization_size());
            as.set_stabilization_weights(prm.elastic_weight(), prm.acoustic_weight());
            as.set_parallel_assembly(m_matrix_free);    // (the parallel copy of the local computations may round
                                                        //   differently under -ffast-math: not by default)
            as.set_skip_zero_entries(m_matrix_free);
            as.set_class_reuse(m_matrix_free);          // local operators once per geometric class (round-off)
            as.set_cache_rhs_data(m_matrix_free && !constant_F);   // F(t) at each stage: time-independent data kept
            as.assemble_mass(msh);
            as.set_legacy_interface_quadrature(prm.space.interface_quadrature == "legacy");
            as.assemble_coupling_terms(msh);
            as.assemble(msh, e_src, a_src, ls_Q);       // the "explicit" branch of the assembler is least squares
            as.LHS += as.COUPLING;
            as.LHS_STAB = SparseMatrix<RealType>();     // copy of the stabilization blocks, not used here (memory)

            p = fine_substeps(prm, dd);
            h_c = fine_size(prm, dd);
            if (prm.time.fine_measure == "width") as.set_fine_measure_width(true);
            as.assemble_P(msh, h_c, prm.time.fine_layers);
            for (bool fine : as.cell_is_fine_silo)
                n_fine_cells += fine;
            separable_F = m_matrix_free && !constant_F && (prm.problem.solution == "lshape_corner" || prm.problem.solution == "lshape_corner_sin");
            lshape_sin = prm.problem.solution == "lshape_corner_sin";
        }
        problem(const problem&) = delete;
        problem& operator=(const problem&) = delete;

        assembler_type& assembler() { return *m_assembler; }

        // ERK operator (blocks of K and M, block-diagonal inverses of the cell mass and of the least-squares face block),
        // built at its first use: not for sdirk without consistent initial faces (memory of the large meshes)
        erk_type& erk() {
            if (!m_erk) {
                assembler_type& as = *m_assembler;
                m_erk = std::make_unique<erk_type>(as.LHS, as.RHS, as.MASS, as.COUPLING,
                    as.get_e_n_cells_dof(), as.get_a_n_cells_dof(), as.get_e_face_dof(), as.get_a_face_dof());
                m_erk->Mcc_inverse(as.get_elastic_cells(), as.get_acoustic_cells(), as.get_e_cell_basis_data(), as.get_a_cell_basis_data());
                if (ls_Q) {
                    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
                    m_erk->Sff_inverse(as.get_elastic_faces(), as.get_acoustic_faces(), as.get_e_face_basis_data(), as.get_a_face_basis_data(),
                        as.get_e_compress(), as.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, dd.interface_face_indexes);
                }
            }
            return *m_erk;
        }

        // Matrix-free explicit operator (solver.matrix_free), built at its first use from the assembled matrices (no ERK
        // operator): the initial faces and the energy use it
        mf_type& mf() {
            if (!m_mf) {
                auto t0 = std::chrono::steady_clock::now();
                assembler_type& as = assembler();
                m_mf = std::make_unique<mf_type>(as.LHS, as.MASS, cell_blocks(),
                    std::vector<std::pair<size_t, size_t>>{{as.get_elastic_faces(), as.get_e_face_basis_data()},
                                                           {as.get_acoustic_faces(), as.get_a_face_basis_data()}},
                    m_mf_tolerance, m_mf_batch);
                if (m_mf_single) m_mf->set_single_precision();
                if (m_mf_li_constraint) m_mf->set_li_constraint(true);
                if (m_mf_li_stage) m_mf->set_li_stage_coupling(true);
                if (m_mf_li_condensation) m_mf->set_li_condensation(true);
                if (m_mf_mlts_reuse) m_mf->set_mlts_reuse(true);
                if (m_mf_profile) m_mf->set_profile(true);
                if (m_mf_par_threshold != 20000) m_mf->set_parallel_threshold(m_mf_par_threshold);
                std::chrono::duration<double> el = std::chrono::steady_clock::now() - t0;
                std::cout << bold << cyan << "      matrix-free operator: " << m_mf->summary() << ", built in " << el.count()
                          << " s" << reset << std::endl;
            }
            return *m_mf;
        }

        // Multilevel LTS (time.scheme = "mlts"): level of each cell (indexed by mesh cell id), coupling_common::
        // mlts_cell_levels. Returns the number of levels.
        size_t mlts_cell_levels(const wave_parameters& prm, std::vector<int>& level) {
            return coupling_common::mlts_cell_levels(msh, dd, prm, level);
        }

        // 1 on the fine unknowns of the LTS / LI integrators (diagonal of Pfine)
        std::vector<char> fine_unknowns() {
            const SparseMatrix<RealType>& P = assembler().Pfine;
            std::vector<char> fine(P.rows(), 0);
            for (int j = 0; j < P.outerSize(); j++)
                for (SparseMatrix<RealType>::InnerIterator it(P, j); it; ++it)
                    if (it.row() == it.col() && it.value() > 0.5) fine[it.row()] = 1;
            return fine;
        }

        // Cell blocks (number of cells, unknowns per cell) and face unknowns, for the static condensation
        std::vector<std::pair<size_t, size_t>> cell_blocks() {
            return {std::make_pair(assembler().get_e_material_data().size(), assembler().get_e_cell_basis_data()),
                    std::make_pair(assembler().get_a_material_data().size(), assembler().get_a_cell_basis_data())};
        }
        size_t n_face_dof() { return assembler().get_n_face_dof(); }

        // Ricker wavelet amplitude theta exp(-pi^2 r^2 / lambda^2), lambda = velocity / frequency
        RealType ricker_amplitude(const point_type& pt) const {
            RealType dx = pt.x() - ricker.center[0], dy = pt.y() - ricker.center[1];
            RealType lambda = ricker.velocity / ricker.frequency;
            return ricker.amplitude * std::exp(-M_PI * M_PI * (dx * dx + dy * dy) / (lambda * lambda));
        }

        // Converging tone burst of a focused transducer (problem.ricker, shape = "focused"): around the focus F = center,
        // r = |x - F|, pressure A W(theta) exp(-((r - r_c) / s)^2) cos(2 pi (r - r_c) / lambda), lambda = velocity /
        // frequency, r_c = radius - cycles lambda / 2, s = cycles lambda / 4, W(theta) = exp(-(theta / aperture)^8) the
        // window in the angle between x - F and the axis; acoustic flux of the one-way wave travelling towards F,
        // m = -p / (density velocity) (F - x) / r (convention of the flux of this code, checked on a converging wave)
        static constexpr double focused_flux_sign = -1.0;  // flux of a wave travelling along d: -p d / (rho c)
        RealType focused_pressure(const point_type& pt) const {
            const RealType dx = pt.x() - ricker.center[0], dy = pt.y() - ricker.center[1];
            const RealType r = std::sqrt(dx * dx + dy * dy);
            if (r == 0) return 0;
            const RealType a = ricker.axis * M_PI / 180.0, al = ricker.aperture * M_PI / 180.0;
            const RealType c = std::max(RealType(-1), std::min(RealType(1), (dx * std::cos(a) + dy * std::sin(a)) / r));
            const RealType th = std::acos(c) / al;
            const RealType lambda = ricker.velocity / ricker.frequency;
            const RealType rc = ricker.radius - 0.5 * ricker.cycles * lambda, sg = 0.25 * ricker.cycles * lambda;
            const RealType u = (r - rc) / sg;
            return ricker.amplitude * std::exp(-std::pow(th, 8)) * std::exp(-u * u) * std::cos(2 * M_PI * (r - rc) / lambda);
        }

        // Projection of the exact solution (or of the Ricker wavelet) at t = 0, faces consistent with the cells
        // (K_ff x_F = -K_FT x_T) unless consistent_faces is false
        vector_type initial_state() {
            m_t = 0.0;
            vector_type x;
            if (ricker_Q) {
                auto zero_v = [](const point_type&) -> disk::static_vector<RealType, 2> { return disk::static_vector<RealType, 2>::Zero(); };
                auto zero_m = [](const point_type&) -> disk::static_matrix<RealType, 2, 2> { return disk::static_matrix<RealType, 2, 2>::Zero(); };
                auto zero_s = [](const point_type&) -> RealType { return 0.0; };
                auto wave_v = [this](const point_type& pt) -> disk::static_vector<RealType, 2> {
                    RealType w = ricker_amplitude(pt);
                    return disk::static_vector<RealType, 2>{w * (pt.x() - ricker.center[0]), w * (pt.y() - ricker.center[1])};
                };
                auto wave_s = [this](const point_type& pt) -> RealType { return ricker_amplitude(pt); };
                std::function<disk::static_vector<RealType, 2>(const point_type&)> e_v = zero_v, a_q = zero_v;
                std::function<RealType(const point_type&)> a_p = zero_s;
                if (ricker.shape == "focused") {
                    a_p = [this](const point_type& pt) -> RealType { return focused_pressure(pt); };
                    a_q = [this](const point_type& pt) -> disk::static_vector<RealType, 2> {
                        const RealType dx = ricker.center[0] - pt.x(), dy = ricker.center[1] - pt.y();
                        const RealType r = std::sqrt(dx * dx + dy * dy);
                        if (r == 0) return disk::static_vector<RealType, 2>::Zero();
                        const RealType w = focused_pressure(pt) / (ricker.density * ricker.velocity * r);
                        return disk::static_vector<RealType, 2>{focused_flux_sign * w * dx, focused_flux_sign * w * dy};
                    };
                }
                else if (ricker.field == "acoustic_flux")    a_q = wave_v;
                else if (ricker.field == "acoustic_pressure") a_p = wave_s;
                else                                          e_v = wave_v;
                assembler().project_over_cells(msh, x, e_v, zero_m, a_p, a_q);
                assembler().project_over_faces(msh, x, e_v, a_p);
            }
            else {
                assembler().project_over_cells(msh, x, v_fun, flux_fun, s_v_fun, s_flux_fun);
                assembler().project_over_faces(msh, x, v_fun, s_v_fun);
            }
            if (!consistent_faces)
                return x;
            if (ls_Q && m_matrix_free)
                mf().refresh_faces(x);
            else if (ls_Q)
                erk().refresh_faces_unknowns(x);
            else {
                const size_t n_c = assembler().get_e_n_cells_dof() + assembler().get_a_n_cells_dof();
                const size_t n_f = assembler().get_n_face_dof();
                SparseMatrix<RealType> Kff = assembler().LHS.block(n_c, n_c, n_f, n_f);
                SparseMatrix<RealType> Kfc = assembler().LHS.block(n_c, 0, n_f, n_c);
                SparseLU<SparseMatrix<RealType>> lu(Kff);
                x.tail(n_f) = lu.solve(vector_type(-(Kfc * x.head(n_c))));
            }
            return x;
        }

        // Random cell unknowns (uniform in [-1, 1], fixed seed) with consistent faces: initial data of the CFL study
        // exciting all the modes (study.cfl_random)
        vector_type random_state() {
            vector_type x = initial_state();
            const size_t n_c = assembler().get_e_n_cells_dof() + assembler().get_a_n_cells_dof();
            std::mt19937 gen(12345);
            std::uniform_real_distribution<RealType> u(-1.0, 1.0);
            for (size_t i = 0; i < n_c; i++) x(i) = u(gen);
            consistent_face_unknowns(x);
            return x;
        }

        // Face unknowns of x from its cell unknowns (static condensation)
        void consistent_face_unknowns(vector_type& x) {
            if (ls_Q && m_matrix_free)
                mf().refresh_faces(x);
            else if (ls_Q)
                erk().refresh_faces_unknowns(x);
            else {
                const size_t n_c = assembler().get_e_n_cells_dof() + assembler().get_a_n_cells_dof();
                const size_t n_f = assembler().get_n_face_dof();
                SparseMatrix<RealType> Kff = assembler().LHS.block(n_c, n_c, n_f, n_f);
                SparseMatrix<RealType> Kfc = assembler().LHS.block(n_c, 0, n_f, n_c);
                SparseLU<SparseMatrix<RealType>> lu(Kff);
                x.tail(n_f) = lu.solve(vector_type(-(Kfc * x.head(n_c))));
            }
        }

        // F(t): sources and Dirichlet contributions at time t (assembled once when it does not depend on time)
        vector_type eval_F(RealType t) {
            if (separable_F) {
                prepare_separable_F();
                std::vector<RealType> a = F_coeffs(t);
                vector_type F = a[0] * F_basis[0];
                for (size_t j = 1; j < F_basis.size(); j++) F += a[j] * F_basis[j];
                return F;
            }
            m_t = t;
            if (constant_F) {
                if (!m_F_ready) {
                    assembler().assemble_rhs(msh, e_src, a_src, ls_Q);
                    m_F = assembler().RHS;
                    m_F_ready = true;
                }
                return m_F;
            }
            assembler().assemble_rhs(msh, e_src, a_src, ls_Q);
            return assembler().RHS;
        }

        // F_basis, F_coeffs of the separable F: F(t) = a_D(t) F_D + a_S(t) F_S (Dirichlet part, source part), with
        // lshape_corner a = (-w sin(w t), -w^2 cos(w t)), lshape_corner_sin a = (w cos(w t), -w^2 sin(w t)); F_D, F_S
        // from F(0) and F(pi / (2 w)) assembled once (2 x 2 system)
        void prepare_separable_F() {
            if (!F_basis.empty()) return;
            const RealType w = scal_vec_analytic_functions::lshape_w, t1 = 0.0, t2 = 0.5 * M_PI / w;
            if (lshape_sin)
                F_coeffs = [w](RealType t) { return std::vector<RealType>{w * std::cos(w * t), -w * w * std::sin(w * t)}; };
            else
                F_coeffs = [w](RealType t) { return std::vector<RealType>{-w * std::sin(w * t), -w * w * std::cos(w * t)}; };
            m_t = t1;
            assembler().assemble_rhs(msh, e_src, a_src, ls_Q);
            vector_type F1 = assembler().RHS;
            m_t = t2;
            assembler().assemble_rhs(msh, e_src, a_src, ls_Q);
            vector_type F2 = assembler().RHS;
            const std::vector<RealType> a1 = F_coeffs(t1), a2 = F_coeffs(t2);
            const RealType det = a1[0] * a2[1] - a1[1] * a2[0];
            vector_type FD = (a2[1] * F1 - a1[1] * F2) / det;
            vector_type FS = (a1[0] * F2 - a2[0] * F1) / det;
            F_basis = {FD, FS};
        }

        // F(t) on the unknowns dofs (written in out[d]): separable F only
        void eval_F_local(RealType t, const std::vector<int>& dofs, RealType* out) {
            prepare_separable_F();
            std::vector<RealType> a = F_coeffs(t);
            for (int d : dofs) {
                RealType v = a[0] * F_basis[0][d];
                for (size_t j = 1; j < F_basis.size(); j++) v += a[j] * F_basis[j][d];
                out[d] = v;
            }
        }

        // Discrete energy 1/2 x^T M x (cell unknowns: velocity and stress, scalar and flux)
        RealType energy(const vector_type& x) {
            if (m_mf) {
                auto e = mf_energy_parts(x);
                return e[0] + e[1];
            }
            return 0.5 * x.dot(assembler().MASS * x);
        }

        // energy() is the sum of energy_parts() (one pass for both: matrix-free operator built)
        bool energy_from_parts() const { return bool(m_mf); }

        // {acoustic, elastic} parts of the discrete energy (M is block diagonal, elastic cells first)
        std::array<RealType, 2> energy_parts(const vector_type& x) {
            if (m_mf) return mf_energy_parts(x);
            vector_type Mx = assembler().MASS * x;
            const size_t ne = assembler().get_e_n_cells_dof(), na = assembler().get_a_n_cells_dof();
            return {0.5 * x.segment(ne, na).dot(Mx.segment(ne, na)), 0.5 * x.head(ne).dot(Mx.head(ne))};
        }

        // {h_min, elastic velocity, elastic stress, acoustic scalar, acoustic flux} L2 errors at time t
        std::array<RealType, 5> errors(vector_type& x, RealType t) {
            m_t = t;
            return postprocessor<mesh_type>::errors_four_fields_elastoacoustic(msh, hho_di, assembler(), x, v_fun, flux_fun, s_v_fun, s_flux_fun);
        }

        // Components of the errors in the norms of Mottier-Ern-Khot-Guillot (6.3) at time t (output.errors = "energy"),
        // names in energy_error_names():
        //   hho_p, hho_v      ||Pi p - p_T||_{L2(1/kappa)}, ||Pi v - v_T||_{L2(rho)}: cell (HHO) unknowns against the L2
        //                     projections of the exact solution (as errors_four_fields_elastoacoustic)
        //   dg_m, dg_s        ||m - m_T||_{L2(rho)}, ||s - s_T||_{L2(C^-1)}: dual (dG) unknowns, by quadrature
        //   dg_m_projection, dg_s_projection    the same against the L2 projections of m (gradients of P^{k'+1}) and
        //                     of s (P^k_sym), as compute_errors_four_fields_elastoacoustic_energy_norm
        //   dg_s_plain, dg_s_projection_plain   ||s - s_T||, ||Pi s - s_T|| in the plain L2 norm
        // with C^-1 s = (s - lambda / (2 (lambda + mu)) tr(s) I) / (2 mu) (plane strain)
        static const std::vector<std::string>& energy_error_names() {
            static const std::vector<std::string> n = {"hho_p", "hho_v", "dg_m", "dg_s", "dg_m_projection", "dg_s_projection",
                                                       "dg_s_plain", "dg_s_projection_plain"};
            return n;
        }
        std::vector<RealType> energy_errors(vector_type& x, RealType t) {
            m_t = t;
            typedef Matrix<RealType, Dynamic, Dynamic> dense;
            auto& as = assembler();
            const auto& di = hho_di;
            const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
            const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
            const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
            const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
            const size_t e_cdof = n_ten + n_vec, a_cdof = n_flux + n_scal;
            const size_t n_e = as.get_e_material_data().size();
            auto storage = msh.backend_storage();
            std::vector<RealType> e(8, 0.0);
            auto cinv_norm2 = [](const disk::static_matrix<RealType, 2, 2>& d, RealType mu, RealType lambda) {
                RealType tr = d(0, 0) + d(1, 1);
                return (d.squaredNorm() - lambda / (2.0 * (lambda + mu)) * tr * tr) / (2.0 * mu);
            };
            size_t ei = 0;
            for (auto& ec : as.get_e_material_data()) {
                auto& cell = storage->surfaces[ec.first];
                const RealType rho = ec.second.rho(), mu = rho * ec.second.vs() * ec.second.vs();
                const RealType lambda = rho * ec.second.vp() * ec.second.vp() - 2.0 * mu;
                vector_type xs = x.segment(ei * e_cdof, n_ten), xv = x.segment(ei * e_cdof + n_ten, n_vec);
                auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                dense mv = disk::make_mass_matrix(msh, cell, vb);
                vector_type dv = vector_type(mv.llt().solve(disk::make_rhs(msh, cell, vb, v_fun))) - xv;
                e[1] += rho * dv.dot(mv * dv);
                auto tb = disk::make_sym_matrix_monomial_basis(msh, cell, di.grad_degree());
                dense mt = disk::make_mass_matrix(msh, cell, tb);
                vector_type ps = mt.llt().solve(disk::make_rhs(msh, cell, tb, flux_fun));
                for (auto& qp : disk::integrate(msh, cell, 2 * (di.cell_degree() + 1))) {
                    auto phi = tb.eval_functions(qp.point());
                    disk::static_matrix<RealType, 2, 2> sh = disk::eval(xs, phi), sp = disk::eval(ps, phi);
                    disk::static_matrix<RealType, 2, 2> d = flux_fun(qp.point()) - sh, dp = sp - sh;
                    e[3] += qp.weight() * cinv_norm2(d, mu, lambda);
                    e[5] += qp.weight() * cinv_norm2(dp, mu, lambda);
                    e[6] += qp.weight() * d.squaredNorm();
                    e[7] += qp.weight() * dp.squaredNorm();
                }
                ei++;
            }
            size_t ai = 0;
            for (auto& ac : as.get_a_material_data()) {
                auto& cell = storage->surfaces[ac.first];
                const RealType rho = ac.second.rho(), kappa = rho * ac.second.vp() * ac.second.vp();
                const size_t off = n_e * e_cdof + ai * a_cdof;
                vector_type xq = x.segment(off, n_flux), xp = x.segment(off + n_flux, n_scal);
                auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                dense mp = disk::make_mass_matrix(msh, cell, cb);
                vector_type dp = vector_type(mp.llt().solve(disk::make_rhs(msh, cell, cb, s_v_fun))) - xp;
                e[0] += dp.dot(mp * dp) / kappa;
                // flux space: gradients of the reconstruction basis without the constant
                auto rb = disk::make_scalar_monomial_basis(msh, cell, di.reconstruction_degree());
                dense mq = dense(disk::make_stiffness_matrix(msh, cell, rb)).block(1, 1, n_flux, n_flux);
                vector_type rq = vector_type::Zero(n_flux);
                auto qps = disk::integrate(msh, cell, 2 * (di.cell_degree() + 1));
                for (auto& qp : qps) {
                    auto dphi = rb.eval_gradients(qp.point());
                    auto m = s_flux_fun(qp.point());
                    for (size_t i = 0; i < n_flux; i++)
                        rq(i) += qp.weight() * (m(0) * dphi(i + 1, 0) + m(1) * dphi(i + 1, 1));
                }
                vector_type pq = mq.llt().solve(rq);
                for (auto& qp : qps) {
                    auto dphi = rb.eval_gradients(qp.point());
                    RealType hx = 0.0, hy = 0.0, px = 0.0, py = 0.0;
                    for (size_t i = 0; i < n_flux; i++) {
                        hx += xq(i) * dphi(i + 1, 0); hy += xq(i) * dphi(i + 1, 1);
                        px += pq(i) * dphi(i + 1, 0); py += pq(i) * dphi(i + 1, 1);
                    }
                    auto m = s_flux_fun(qp.point());
                    e[2] += qp.weight() * rho * ((m(0) - hx) * (m(0) - hx) + (m(1) - hy) * (m(1) - hy));
                    e[4] += qp.weight() * rho * ((px - hx) * (px - hx) + (py - hy) * (py - hy));
                }
                ai++;
            }
            for (auto& v : e) v = std::sqrt(v);
            return e;
        }

    private:
        // energy parts by the matrix-free operator, kept out of line: the code of the other branch of energy and
        // energy_parts (and so their last bits) stays as without the matrix-free operator
        __attribute__((noinline)) std::array<RealType, 2> mf_energy_parts(const vector_type& x) {
            return m_mf->energy_parts(x, assembler().get_e_n_cells_dof());
        }

        std::unique_ptr<assembler_type> m_assembler;
        std::unique_ptr<erk_type>       m_erk;
        std::unique_ptr<mf_type>        m_mf;
        vector_type                     m_F;
        bool                            m_F_ready = false;
    };

    // Elastic medium with the formulation of Burman-Duran-Ern-Steins 2021: elastodynamic_three_fields_assembler, the
    // assembler of the paper (unknowns per cell: stress 2 mu eps, hydrostatic part lambda div u, velocity; Dirichlet data
    // on the velocity; rho = 1, vp = sqrt(3), vs = 1), integrated with sdirk only.
    class problem_three_fields {
    public:
        typedef elastodynamic_three_fields_assembler<mesh_type> assembler3_type;
        mesh_type msh;
        disk::hho_degree_info hho_di;
        domain_data dd;
        scal_vec_analytic_functions functions;
        RealType m_t = 0.0;
        std::function<disk::static_vector<RealType, 2>(const point_type&)>    v_fun, f_fun;
        std::function<disk::static_matrix<RealType, 2, 2>(const point_type&)> flux_fun;
        size_t p = 1;
        size_t n_fine_cells = 0;

        problem_three_fields(const wave_parameters& prm, size_t level)
            : hho_di(prm.mixed_order() ? prm.space.k + 1 : prm.space.k, prm.space.k) {
            build_mesh(msh, prm.mesh, level);
            dd = classify_domain(msh, "elastic");
            functions.set_function_type(solution_type(prm.problem.solution));
            v_fun    = functions.Evaluate_v(m_t);
            f_fun    = functions.Evaluate_f(m_t);
            flux_fun = functions.Evaluate_sigma(m_t);
            e_boundary_type bnd(msh);
            bnd.addDirichletEverywhere(v_fun);
            m_assembler = std::make_unique<assembler3_type>(msh, hho_di, bnd);
            assembler3_type& as = *m_assembler;
            as.load_material_data(msh);
            if (prm.mixed_order())
                as.set_hdg_stabilization();
            else
                as.set_hho_stabilization();
            // weight in 1/h: face diameter h_F if set_scaled_stabilization, cell diameter h_T otherwise
            if (prm.stabilization_size() == 2)
                as.set_scaled_stabilization();
            as.assemble_mass(msh);
            as.assemble(msh, f_fun);
            n_fine_cells = msh.cells_size();
        }
        problem_three_fields(const problem_three_fields&) = delete;
        problem_three_fields& operator=(const problem_three_fields&) = delete;

        assembler3_type& assembler() { return *m_assembler; }

        std::vector<std::pair<size_t, size_t>> cell_blocks() {
            return {std::make_pair(msh.cells_size(), assembler().get_cell_basis_data())};
        }
        size_t n_face_dof() { return assembler().get_n_face_dof(); }

        // Projection of the exact velocity and stress at t = 0 (cells and faces), as in the original code of the paper
        vector_type initial_state() {
            m_t = 0.0;
            vector_type x;
            assembler().project_over_cells(msh, x, v_fun, flux_fun);
            assembler().project_over_faces(msh, x, v_fun);
            return x;
        }

        vector_type eval_F(RealType t) {
            m_t = t;
            assembler().assemble_rhs(msh, f_fun);
            return assembler().RHS;
        }

        RealType energy(const vector_type& x) {
            return 0.5 * x.dot(assembler().MASS * x);
        }

        // {h_min, elastic velocity, elastic stress, -, -} L2 errors at time t (layout of problem::errors)
        std::array<RealType, 5> errors(vector_type& x, RealType t) {
            m_t = t;
            auto e = postprocessor<mesh_type>::errors_three_fields_vectorial(msh, hho_di, x, v_fun, flux_fun);
            return {dd.h_min, e[0], e[1], 0.0, 0.0};
        }

    private:
        std::unique_ptr<assembler3_type> m_assembler;
    };

}

#endif /* first_order_problem_hpp */
