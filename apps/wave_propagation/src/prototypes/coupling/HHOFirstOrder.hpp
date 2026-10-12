//
//  HHOFirstOrder.hpp
//
//  Elasto-acoustic coupling, first-order formulation (velocity-stress in the solid, velocity-pressure in the fluid,
//  elastoacoustic_four_fields_assembler): M y' + K y = F(t). Every setting comes from the Lua input file
//  (common/lua_parameters.hpp, template share/first_order.lua):
//    space : face degree k, equal (k,k) or mixed (k+1,k) order, least-squares stabilization (plain least squares in
//            equal order, Lehrenfeld-Schoberl in mixed order) or HHO stabilization, weight O(1) or in 1/h
//    time  : ERK(s), SDIRK(s,s+1), two-level LTS-RK4 (substeps in the locally refined region), LI-SDIRK(4,3) (coarse
//            ERK4 Taylor polynomial, SDIRK(3,4) in the fine region), final time, number of time steps
//    mesh, problem (manufactured solution), study (one run, space-time convergence, time convergence), output
//  The explicit integrators (erk, lts) and li eliminate the faces with the block-diagonal face block of the
//  least-squares stabilization (erk_coupling_hho_scheme::Sff_inverse); sdirk accepts both stabilizations.
//  Papers' settings (MEG, MEKG): explicit = equal order, plain least squares, O(1); implicit = mixed order,
//  Lehrenfeld-Schoberl, 1/h.
//  space.elastic_formulation = "three_fields" (elastic medium, sdirk): the formulation and the assembler of
//  Burman-Duran-Ern-Steins 2021 (elastodynamic_three_fields_assembler: stress 2 mu eps, hydrostatic part lambda div u
//  and velocity), instead of the velocity-stress block of the coupled assembler.
//

#pragma once
#ifndef HHOFirstOrder_hpp
#define HHOFirstOrder_hpp

#include "coupling_common.hpp"
#include "../../common/elastodynamic_three_fields_assembler.hpp"
#include <chrono>
#include "../LTS_ERK_SDIRK/dirk_coupling_hho_scheme.hpp"

namespace hho_first_order {

    using namespace coupling_common;
    typedef elastoacoustic_four_fields_assembler<mesh_type>  assembler_type;
    typedef erk_coupling_hho_scheme<RealType>               erk_type;

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
        bool consistent_faces = true;   // initial faces from the cells (K_ff x_F = -K_FT x_T), else projected
        wave_parameters::problem_data::ricker_data ricker;

        problem(const wave_parameters& prm, size_t level, bool zero_data = false)
            : hho_di(prm.mixed_order() ? prm.space.k + 1 : prm.space.k, prm.space.k), ls_Q(prm.least_squares()) {

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
            as.set_stabilization_size(prm.stabilization_size());
            as.set_stabilization_weights(prm.elastic_weight(), prm.acoustic_weight());
            as.assemble_mass(msh);
            as.set_legacy_interface_quadrature(prm.space.interface_quadrature == "legacy");
            as.assemble_coupling_terms(msh);
            as.assemble(msh, e_src, a_src, ls_Q);       // the "explicit" branch of the assembler is least squares
            as.LHS += as.COUPLING;
            as.LHS_STAB = SparseMatrix<RealType>();     // copy of the stabilization blocks, not used here (memory)

            p = fine_substeps(prm, dd);
            h_c = fine_size(prm, dd);
            as.assemble_P(msh, h_c);
            for (bool fine : as.cell_is_fine_silo)
                n_fine_cells += fine;
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
                if (ricker.field == "acoustic_flux")         a_q = wave_v;
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
            if (ls_Q)
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

        // F(t): sources and Dirichlet contributions at time t (assembled once when it does not depend on time)
        vector_type eval_F(RealType t) {
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

        // Discrete energy 1/2 x^T M x (cell unknowns: velocity and stress, scalar and flux)
        RealType energy(const vector_type& x) {
            return 0.5 * x.dot(assembler().MASS * x);
        }

        // {acoustic, elastic} parts of the discrete energy (M is block diagonal, elastic cells first)
        std::array<RealType, 2> energy_parts(const vector_type& x) {
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
        std::unique_ptr<assembler_type> m_assembler;
        std::unique_ptr<erk_type>       m_erk;
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

    // Point and segment geometry of the sensors: is pt in the (closed) polygon of a cell, within a tolerance eps?
    inline bool cell_contains(const mesh_type& msh, const typename mesh_type::cell_type& cell, const point_type& pt, RealType eps) {
        auto pts = points(msh, cell);
        const size_t n = pts.size();
        bool inside = false;
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            const auto& a = pts[j];
            const auto& b = pts[i];
            // on the edge [a, b]
            RealType ex = b.x() - a.x(), ey = b.y() - a.y();
            RealType len2 = ex * ex + ey * ey;
            RealType s = len2 > 0.0 ? ((pt.x() - a.x()) * ex + (pt.y() - a.y()) * ey) / len2 : 0.0;
            s = std::min(1.0, std::max(0.0, s));
            RealType dx = a.x() + s * ex - pt.x(), dy = a.y() + s * ey - pt.y();
            if (dx * dx + dy * dy <= eps * eps)
                return true;
            if (((a.y() > pt.y()) != (b.y() > pt.y())) &&
                (pt.x() < (b.x() - a.x()) * (pt.y() - a.y()) / (b.y() - a.y()) + a.x()))
                inside = !inside;
        }
        return inside;
    }

    inline bool face_contains(const mesh_type& msh, const typename mesh_type::face_type& face, const point_type& pt, RealType eps) {
        auto pts = points(msh, face);
        const auto& a = pts[0];
        const auto& b = pts[1];
        RealType ex = b.x() - a.x(), ey = b.y() - a.y();
        RealType len2 = ex * ex + ey * ey;
        RealType s = len2 > 0.0 ? ((pt.x() - a.x()) * ex + (pt.y() - a.y()) * ey) / len2 : 0.0;
        if (s < -1.0e-12 || s > 1.0 + 1.0e-12) return false;
        RealType dx = a.x() + s * ex - pt.x(), dy = a.y() + s * ey - pt.y();
        return dx * dx + dy * dy <= eps * eps;
    }

    // Point sensors (output.sensors): one CSV file <table>_<name>.csv per sensor, a line "t,values" every
    // output.sensor_every steps. Fields: cell values at the point (acoustic_pressure p_T, acoustic_flux m_T,
    // elastic_velocity v_T, elastic_stress s_T: xx, xy, yy) in a cell of the corresponding medium containing the point,
    // and at a point of the interface the face unknowns interface_pressure p_F (acoustic face) and interface_velocity
    // v_F (elastic face) of the interface face containing it.
    class sensors {
        struct item {
            std::string name;
            point_type pt;
            std::vector<std::string> fields;
            long e_cell = -1, a_cell = -1;      // cells of the media containing the point
            long e_face = -1, a_face = -1;      // position of the interface face among the faces of the cell
            std::ofstream os;
        };
        std::vector<std::unique_ptr<item>> m_items;
        size_t m_every = 1;

    public:
        sensors(problem& pb, const wave_parameters& prm, const std::string& base) : m_every(prm.output.sensor_every) {
            auto& msh = pb.msh;
            const RealType eps = 1.0e-8 * std::max(RealType(1.0), pb.dd.h_max);
            for (auto& sn : prm.output.sensors) {
                auto it = std::make_unique<item>();
                it->name = sn.name;
                it->pt = point_type(sn.point[0], sn.point[1]);
                it->fields = sn.fields;
                for (auto& cell : msh) {
                    auto ci = msh.lookup(cell);
                    if (!cell_contains(msh, cell, it->pt, eps)) continue;
                    if (pb.dd.e_material.count(ci) && it->e_cell < 0) it->e_cell = long(ci);
                    if (pb.dd.a_material.count(ci) && it->a_cell < 0) it->a_cell = long(ci);
                }
                auto interface_face = [&](long ci) -> long {
                    auto& cell = msh.backend_storage()->surfaces[size_t(ci)];
                    auto fcs = faces(msh, cell);
                    for (size_t j = 0; j < fcs.size(); j++)
                        if (pb.dd.interface_face_indexes.count(msh.lookup(fcs[j])) && face_contains(msh, fcs[j], it->pt, eps))
                            return long(j);
                    return -1;
                };
                std::string header = "t";
                for (auto& f : it->fields) {
                    bool acoustic = f == "acoustic_pressure" || f == "acoustic_flux" || f == "interface_pressure";
                    long ci = acoustic ? it->a_cell : it->e_cell;
                    if (ci < 0)
                        throw std::runtime_error("sensor " + it->name + ": no " + std::string(acoustic ? "acoustic" : "elastic")
                                                 + " cell contains the point for the field " + f);
                    if (f == "interface_pressure") {
                        it->a_face = interface_face(it->a_cell);
                        if (it->a_face < 0) throw std::runtime_error("sensor " + it->name + ": the point is not on the interface");
                    }
                    if (f == "interface_velocity") {
                        it->e_face = interface_face(it->e_cell);
                        if (it->e_face < 0) throw std::runtime_error("sensor " + it->name + ": the point is not on the interface");
                    }
                    if (f == "acoustic_pressure")  header += ",p";
                    if (f == "acoustic_flux")      header += ",mx,my";
                    if (f == "elastic_velocity")   header += ",vx,vy";
                    if (f == "elastic_stress")     header += ",sxx,sxy,syy";
                    if (f == "interface_pressure") header += ",pF";
                    if (f == "interface_velocity") header += ",vFx,vFy";
                }
                it->os.open(base + "_" + it->name + ".csv");
                it->os << "# sensor " << it->name << " (" << sn.point[0] << ", " << sn.point[1] << ")" << std::endl
                       << header << std::endl << std::setprecision(12) << std::scientific;
                m_items.push_back(std::move(it));
            }
        }

        bool empty() const { return m_items.empty(); }

        void record(problem& pb, size_t step, RealType t, const vector_type& x) {
            if (step % m_every != 0) return;
            auto& msh = pb.msh;
            auto& as = pb.assembler();
            const auto& di = pb.hho_di;
            const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
            const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
            const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
            const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
            for (auto& it : m_items) {
                it->os << t;
                vector_type xe, xa;
                if (it->e_cell >= 0) xe = as.gather_e_dof_data(size_t(it->e_cell), msh, msh.backend_storage()->surfaces[size_t(it->e_cell)], x);
                if (it->a_cell >= 0) xa = as.gather_a_dof_data(size_t(it->a_cell), msh, msh.backend_storage()->surfaces[size_t(it->a_cell)], x);
                for (auto& f : it->fields) {
                    if (f == "acoustic_pressure") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                        it->os << "," << xa.segment(n_flux, n_scal).dot(cb.eval_functions(it->pt));
                    }
                    else if (f == "acoustic_flux") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto rb = disk::make_scalar_monomial_basis(msh, cell, di.reconstruction_degree());
                        auto dphi = rb.eval_gradients(it->pt);
                        RealType mx = 0.0, my = 0.0;
                        for (long i = 1; i < dphi.rows(); i++) { mx += xa(i - 1) * dphi(i, 0); my += xa(i - 1) * dphi(i, 1); }
                        it->os << "," << mx << "," << my;
                    }
                    else if (f == "elastic_velocity") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                        auto v = disk::eval(vector_type(xe.segment(n_ten, n_vec)), vb.eval_functions(it->pt));
                        it->os << "," << v(0) << "," << v(1);
                    }
                    else if (f == "elastic_stress") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto tb = disk::make_sym_matrix_monomial_basis(msh, cell, di.grad_degree());
                        auto st = disk::eval(vector_type(xe.head(n_ten)), tb.eval_functions(it->pt));
                        it->os << "," << st(0, 0) << "," << st(0, 1) << "," << st(1, 1);
                    }
                    else if (f == "interface_pressure") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto fc = faces(msh, cell)[size_t(it->a_face)];
                        auto fb = disk::make_scalar_monomial_basis(msh, fc, di.face_degree());
                        const size_t n_fbs = disk::scalar_basis_size(di.face_degree(), 1);
                        const size_t a_cbs = as.get_a_cell_basis_data();
                        it->os << "," << xa.segment(a_cbs + size_t(it->a_face) * n_fbs, n_fbs).dot(fb.eval_functions(it->pt));
                    }
                    else if (f == "interface_velocity") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto fc = faces(msh, cell)[size_t(it->e_face)];
                        auto fb = disk::make_vector_monomial_basis(msh, fc, di.face_degree());
                        const size_t n_fbs = disk::vector_basis_size(di.face_degree(), 1, 2);
                        const size_t e_cbs = as.get_e_cell_basis_data();
                        auto v = disk::eval(vector_type(xe.segment(e_cbs + size_t(it->e_face) * n_fbs, n_fbs)), fb.eval_functions(it->pt));
                        it->os << "," << v(0) << "," << v(1);
                    }
                }
                it->os << "\n";
            }
        }
    };

    // Snapshots (output.snapshots): at the time step nearest to each requested time, the file
    // <table>_snapshot_<i>.dat with one line per cell, the acoustic pressure (acoustic cells) or the norm of the elastic
    // velocity (elastic cells) of the cell polynomial at the barycenter, then at the vertices of the cell; the cells
    // (medium, number of vertices, coordinates) once in <table>_cells.dat.
    class snapshots {
        std::vector<size_t> m_steps;
        std::vector<RealType> m_times;
        std::string m_base;
    public:
        snapshots(problem& pb, const wave_parameters& prm, const std::string& base, RealType dt) : m_base(base) {
            if (prm.output.snapshots.empty()) return;
            for (auto t : prm.output.snapshots) {
                m_steps.push_back(size_t(std::llround(t / dt)));
                m_times.push_back(t);
            }
            std::ofstream os(base + "_cells.dat");
            os << "# cells: medium (0 elastic, 1 acoustic), number of vertices, x_1 y_1 ... x_n y_n" << std::endl;
            os << std::setprecision(12);
            for (auto& cell : pb.msh) {
                auto ci = pb.msh.lookup(cell);
                auto pts = points(pb.msh, cell);
                os << (pb.dd.a_material.count(ci) ? 1 : 0) << " " << pts.size();
                for (auto& q : pts) os << " " << q.x() << " " << q.y();
                os << "\n";
            }
        }

        bool empty() const { return m_steps.empty(); }

        void maybe_write(problem& pb, size_t step, RealType t, const vector_type& x) {
            for (size_t i = 0; i < m_steps.size(); i++) {
                if (m_steps[i] != step) continue;
                auto& msh = pb.msh;
                auto& as = pb.assembler();
                const auto& di = pb.hho_di;
                const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
                const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
                const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
                const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
                std::ofstream os(m_base + "_snapshot_" + std::to_string(i) + ".dat");
                os << "# t = " << std::setprecision(12) << t << " (requested " << m_times[i] << "): acoustic pressure "
                   << "(acoustic cells), |elastic velocity| (elastic cells) at the barycenter, then at the vertices" << std::endl;
                os << std::setprecision(8) << std::scientific;
                for (auto& cell : msh) {
                    auto ci = msh.lookup(cell);
                    auto pts = points(msh, cell);
                    std::vector<point_type> at(1, barycenter(msh, cell));
                    at.insert(at.end(), pts.begin(), pts.end());
                    if (pb.dd.a_material.count(ci)) {
                        vector_type xa = as.gather_a_dof_data(ci, msh, cell, x);
                        vector_type xp = xa.segment(n_flux, n_scal);
                        auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                        for (size_t j = 0; j < at.size(); j++)
                            os << (j ? " " : "") << xp.dot(cb.eval_functions(at[j]));
                    }
                    else {
                        vector_type xe = as.gather_e_dof_data(ci, msh, cell, x);
                        vector_type xv = xe.segment(n_ten, n_vec);
                        auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                        for (size_t j = 0; j < at.size(); j++) {
                            auto v = disk::eval(xv, vb.eval_functions(at[j]));
                            os << (j ? " " : "") << std::sqrt(v(0) * v(0) + v(1) * v(1));
                        }
                    }
                    os << "\n";
                }
            }
        }
    };

    // Called after each time step with the time and the state (energy history)
    typedef std::function<void(RealType, const vector_type&)> observer_type;

    // ERK(s)
    inline vector_type integrate_erk(problem& pb, vector_type x, size_t nt, RealType tf, size_t s, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        RealType dt = tf / nt;
        std::vector<vector_type> k(s);
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            vector_type x_n = x;
            for (size_t i = 0; i < s; i++) {
                vector_type y = x_n;
                for (size_t j = 0; j < i; j++)
                    y += a(i,j) * dt * k[j];
                vector_type F = pb.eval_F(tn + c(i) * dt);
                pb.erk().SetFg(F);
                pb.erk().erk_weight(y, k[i]);
                x += dt * b(i) * k[i];
            }
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

    // SDIRK(s, s+1): (M + a_ii dt K) k_i = F(t_i) - K y_i, one factorization since a_ii is constant
    // iterative_tol > 0: BiCGSTAB with incomplete LU (relative residual iterative_tol) instead of the direct solver
    template<typename Problem>
    inline vector_type integrate_sdirk(Problem& pb, vector_type x, size_t nt, RealType tf, size_t s, bool sc_Q,
                                       const observer_type& observe = {}, RealType iterative_tol = -1.0) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(int(s), a, b, c);
        RealType dt = tf / nt;
        auto& as = pb.assembler();
        // reference mode: no copies of K and M (memory of the large meshes), same arithmetic
        dirk_hho_scheme<RealType> dirk_an(as.LHS, as.RHS, as.MASS, typename dirk_hho_scheme<RealType>::by_reference_t{});
        if (sc_Q)
            dirk_an.set_static_condensation_data(pb.cell_blocks(), pb.n_face_dof());
        RealType scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        dirk_an.ComposeMatrix();
        if (iterative_tol > 0.0)
            dirk_an.setIterativeSolver(iterative_tol);
        dirk_an.DecomposeMatrix();
        std::vector<vector_type> k(s);
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            vector_type x_n = x;
            for (size_t i = 0; i < s; i++) {
                vector_type y = x_n;
                for (size_t j = 0; j < i; j++)
                    y += a(i,j) * dt * k[j];
                vector_type F = pb.eval_F(tn + c(i) * dt);
                dirk_an.SetFg(F);
                dirk_an.irk_weight(y, k[i], dt, a(i,i), true);
                x += dt * b(i) * k[i];
            }
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

    // Two-level LTS-RK4: coarse Taylor coefficients w once per step, p fine RK4 substeps
    inline vector_type integrate_lts(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        RealType dt = tf / nt;
        size_t p = pb.p;
        RealType dtau = dt / p;
        assembler_type& as = pb.assembler();
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            std::vector<vector_type> w(4, vector_type::Zero(x.rows()));
            pb.erk().ZeroFc();
            if (p != 1) {
                vector_type Fn = pb.eval_F(tn), Fn12 = pb.eval_F(tn + 0.5*dt), Fn1 = pb.eval_F(tn + dt);
                pb.erk().ZeroFc();
                pb.erk().erk_weight_LTS_coarse(x, as.Pcoarse, w, Fn, Fn12, Fn1, dt);
            }
            for (size_t m = 0; m < p; m++) {
                RealType tm = m * dtau;
                vector_type Fm = pb.eval_F(tn + tm), Fmh = pb.eval_F(tn + tm + 0.5*dtau), Fm1 = pb.eval_F(tn + tm + dtau);
                pb.erk().erk_weight_LTS_fine(x, as.Pfine, w, Fm, Fmh, Fm1, tm, dtau);
            }
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

    // LI-SDIRK(3,4): coarse ERK4 Taylor coefficients, one implicit SDIRK(3,4) step in the fine region
    inline vector_type integrate_li(problem& pb, vector_type x, size_t nt, RealType tf, const observer_type& observe = {}) {
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        dirk_butcher_tableau::sdirk_tables(3, a, b, c);
        RealType gamma = a(0,0);
        RealType dt = tf / nt;
        assembler_type& as = pb.assembler();
        dirk_coupling_hho_scheme<RealType> dirk_op;
        dirk_op.compose(pb.erk(), as.Pfine, gamma, dt);
        std::function<vector_type(RealType)> eval_F = [&](RealType t) { return pb.eval_F(t); };
        for (size_t it = 0; it < nt; it++) {
            RealType tn = dt * it;
            std::vector<vector_type> w(4, vector_type::Zero(x.rows()));
            pb.erk().ZeroFc();
            vector_type Fn = pb.eval_F(tn), Fn12 = pb.eval_F(tn + 0.5*dt), Fn1 = pb.eval_F(tn + dt);
            pb.erk().ZeroFc();
            pb.erk().erk_weight_LTS_coarse_global(x, as.Pcoarse, w, Fn, Fn12, Fn1, dt);
            dirk_weight_LTS_fine<RealType>(pb.erk(), dirk_op, x, as.Pfine, w, tn, dt, eval_F, a, b, c, gamma);
            if (observe) observe(dt * (it + 1), x);
        }
        return x;
    }

    // iterative solver tolerance of the mesh level (space study: multiplied by solver.tolerance_growth per level)
    inline RealType level_tolerance(const wave_parameters& prm, size_t level) {
        if (prm.study.type != "space" || prm.solver.tolerance_growth == 1.0) return prm.solver.tolerance;
        return prm.solver.tolerance * std::pow(prm.solver.tolerance_growth, RealType(level - prm.mesh.level));
    }

    template<typename Problem>
    inline vector_type integrate(Problem& pb, const wave_parameters& prm, vector_type x, size_t nt, const observer_type& observe = {},
                                 size_t level = 0) {
        RealType tf = prm.time.final_time;
        if (prm.time.scheme == "sdirk")
            return integrate_sdirk(pb, x, nt, tf, prm.stages(), prm.solver.static_condensation, observe,
                                   prm.solver.iterative ? level_tolerance(prm, std::max(level, prm.mesh.level)) : -1.0);
        if constexpr (std::is_same_v<Problem, problem>) {
            if (prm.time.scheme == "erk")   return integrate_erk(pb, x, nt, tf, prm.stages(), observe);
            if (prm.time.scheme == "lts")   return integrate_lts(pb, x, nt, tf, observe);
            return integrate_li(pb, x, nt, tf, observe);
        }
        throw std::runtime_error("space.elastic_formulation = \"three_fields\": time.scheme = \"sdirk\" only");
    }

    // Base name of the outputs: the data table without its extension
    inline std::string output_base(const wave_parameters& prm) {
        std::string f = prm.table_file();
        return f.substr(0, f.size() - 4);
    }

    // Integration with the outputs of the run: history of the discrete energy 1/2 x^T M x (file prm.energy_file():
    // "# t E" header with max |E^n - E^0| / E^0, one line per step; acoustic and elastic parts E_acoustic E_elastic in
    // two more columns for a coupled domain), sensors and snapshots
    template<typename Problem>
    inline vector_type integrate_observed(Problem& pb, const wave_parameters& prm, vector_type x, size_t nt, size_t level = 0) {
        const RealType dt = prm.time.final_time / RealType(nt);
        const std::string base = output_base(prm);
        const bool energy_Q = prm.output.energy;
        std::ostringstream hist;
        hist << std::setprecision(16) << std::scientific;
        const RealType E0 = pb.energy(x);
        RealType var = 0.0;
        // acoustic and elastic parts when both media are present
        bool parts_Q = false;
        if constexpr (std::is_same_v<Problem, problem>)
            parts_Q = !pb.dd.e_material.empty() && !pb.dd.a_material.empty();
        auto energy_line = [&](RealType t, const vector_type& y) {
            RealType E = pb.energy(y);
            var = std::max(var, std::abs(E - E0) / E0);
            hist << t << " " << E;
            if constexpr (std::is_same_v<Problem, problem>) {
                if (parts_Q) {
                    auto parts = pb.energy_parts(y);
                    hist << " " << parts[0] << " " << parts[1];
                }
            }
            hist << "\n";
        };
        std::unique_ptr<sensors> sens;
        std::unique_ptr<snapshots> snaps;
        if constexpr (std::is_same_v<Problem, problem>) {
            if (!prm.output.sensors.empty())   sens  = std::make_unique<sensors>(pb, prm, base);
            if (!prm.output.snapshots.empty()) snaps = std::make_unique<snapshots>(pb, prm, base, dt);
            if (sens)  sens->record(pb, 0, 0.0, x);
            if (snaps) snaps->maybe_write(pb, 0, 0.0, x);
        }
        if (energy_Q) energy_line(0.0, x);
        size_t step = 0;
        x = integrate(pb, prm, x, nt, [&](RealType t, const vector_type& y) {
            step++;
            if (energy_Q) energy_line(t, y);
            if constexpr (std::is_same_v<Problem, problem>) {
                if (sens)  sens->record(pb, step, t, y);
                if (snaps) snaps->maybe_write(pb, step, t, y);
            }
        }, level);
        if (energy_Q) {
            std::ofstream os(prm.energy_file());
            os << "# t E" << (parts_Q ? " E_acoustic E_elastic" : "")
               << " (discrete energy 1/2 x^T M x), max |E - E0| / E0 = " << var << std::endl << hist.str();
            std::cout << bold << cyan << "      energy: max |E^n - E^0| / E^0 = " << var << reset << std::endl;
        }
        return x;
    }

    // CFL study (Mottier-Ern-Guillot, Algorithm 1): ERK(s) with the initial data of the solution, homogeneous Dirichlet
    // data and no source; starting from time.n_steps (doubled until stable, then decreased by 10 % while stable), the
    // number of steps N is decreased by cfl_decrement N + 1 after each stable run, until a run where the discrete
    // energy grows by more than
    // cfl_threshold of the initial or of the previous energy (or is not finite); cfl_absolute: variation in absolute
    // value instead of growth. Table: last stable and first unstable N, CFL* = c_max dt / h with dt = T / (last stable
    // N), h = h_min and h_max (cell diameters).
    inline void cfl_study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        problem pb(prm, prm.mesh.level, true);
        const vector_type x0 = pb.initial_state();
        const RealType E0 = pb.energy(x0), tf = prm.time.final_time, thr = prm.study.cfl_threshold;
        const size_t s = prm.stages();
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        const vector_type F = pb.eval_F(0.0);
        auto stable_run = [&](size_t N) {
            const RealType dt = tf / RealType(N);
            vector_type x = x0;
            RealType E_prev = E0;
            std::vector<vector_type> k(s);
            for (size_t it = 0; it < N; it++) {
                vector_type x_n = x;
                for (size_t i = 0; i < s; i++) {
                    vector_type y = x_n;
                    for (size_t j = 0; j < i; j++)
                        y += a(i, j) * dt * k[j];
                    vector_type Fi = F;
                    pb.erk().SetFg(Fi);
                    pb.erk().erk_weight(y, k[i]);
                    x += dt * b(i) * k[i];
                }
                RealType E = pb.energy(x);
                RealType d_prev = (E - E_prev) / E_prev, d_0 = (E - E0) / E0;
                if (prm.study.cfl_absolute) { d_prev = std::abs(d_prev); d_0 = std::abs(d_0); }
                if (!std::isfinite(E) || d_prev > thr || d_0 > thr)
                    return false;
                E_prev = E;
            }
            return true;
        };
        size_t N = prm.time.n_steps, stable = 0, unstable = 0;
        for (size_t tries = 0; tries < 12 && !stable_run(N); tries++)
            N *= 2;
        // coarse approach of the limit (steps of 10 %), then the decrements of Algorithm 1 from the last stable N
        for (size_t M = N - (N / 10 + 1); M > 0 && stable_run(M); M -= M / 10 + 1)
            N = M;
        for (size_t run = 0; run < prm.study.cfl_runs && N > 0; run++) {
            bool ok = stable_run(N);
            std::cout << bold << cyan << "      ERK(" << s << ") N = " << N << ": " << (ok ? "stable" : "unstable") << reset << std::endl;
            if (!ok) { unstable = N; break; }
            stable = N;
            N -= size_t(RealType(N) * prm.study.cfl_decrement) + 1;
        }
        if (stable == 0)
            throw std::runtime_error("cfl study: no stable number of steps found from time.n_steps");
        const RealType dt = tf / RealType(stable);
        names = {"stable_steps", "unstable_steps", "cfl_h_min", "cfl_h_max"};
        rows.push_back({pb.dd.h_max, dt, {RealType(stable), RealType(unstable), pb.dd.c_max * dt / pb.dd.h_min,
                                          pb.dd.c_max * dt / pb.dd.h_max}, pb.dd.h_min});
        std::cout << bold << cyan << "      c_max = " << pb.dd.c_max << "  h_min = " << pb.dd.h_min << "  h_max = " << pb.dd.h_max << reset << std::endl;
        print_table(names, rows, false, std::cout);
        print_table(names, rows, false, log);
    }

    // Spectral study (Mottier-Ern-Khot-Guillot, Section 6.1): largest eigenvalue gamma of
    // K_schur^T M_T^-1 K_schur X = gamma M_T X, K_schur = K_TT - K_TF K_FF^-1 K_FT the Schur complement of the face block
    // of the stiffness (coupling included), M_T the cell mass matrix; spectral radius rho = sqrt(gamma), normalized by
    // sqrt(number of cells). Dense computation, on each mesh of the study.
    inline void spectrum_study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        typedef Matrix<RealType, Dynamic, Dynamic> dense_type;
        names = {"rho", "rho_normalized", "cells"};
        for (size_t l : study_levels(prm)) {
            problem pb(prm, l, true);
            auto& as = pb.assembler();
            const long n_c = long(as.get_e_n_cells_dof() + as.get_a_n_cells_dof()), n_f = long(as.get_n_face_dof());
            dense_type K = dense_type(as.LHS), M = dense_type(as.MASS);
            dense_type S = K.topLeftCorner(n_c, n_c);
            if (n_f > 0)
                S -= K.topRightCorner(n_c, n_f) * K.bottomRightCorner(n_f, n_f).partialPivLu().solve(dense_type(K.bottomLeftCorner(n_f, n_c)));
            dense_type Mc = M.topLeftCorner(n_c, n_c);
            dense_type A = S.transpose() * Mc.llt().solve(S);
            A = 0.5 * (A + A.transpose()).eval();
            Eigen::GeneralizedSelfAdjointEigenSolver<dense_type> es(A, Mc, Eigen::EigenvaluesOnly);
            const RealType gamma = es.eigenvalues().maxCoeff();
            const RealType cells = RealType(pb.msh.cells_size());
            rows.push_back({pb.dd.h_max, 0.0, {std::sqrt(gamma), std::sqrt(gamma) / std::sqrt(cells), cells}, pb.dd.h_min});
            std::cout << bold << cyan << "      l = " << l << "  cells = " << cells << "  cell dofs = " << n_c << "  face dofs = " << n_f << reset << std::endl;
            print_table(names, rows, false, std::cout);
        }
        print_table(names, rows, false, log);
    }

}

namespace hho_first_order {

    // Space (or single run) and time studies of the first-order formulation for a problem type (problem: coupled
    // assembler; problem_three_fields: three-field formulation of the elastic medium): convergence rows and names
    template<typename Problem>
    inline void study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        error_columns all = columns(prm, {"elastic_v", "elastic_sigma"}, {"acoustic_p", "acoustic_q"});
        auto all_errors = [](const std::array<RealType, 5>& e) { return std::vector<RealType>{e[1], e[2], e[3], e[4]}; };
        const RealType tf = prm.time.final_time;
        // no error for the Ricker wavelet; energy history, sensors and snapshots of the last (or only) mesh
        const bool ricker = prm.problem.solution == "ricker";
        const bool observed = prm.output.energy || !prm.output.sensors.empty() || !prm.output.snapshots.empty();
        // components of the errors in the energy norms (output.errors = "energy", problem::energy_error_names)
        const bool energy_norms = prm.output.errors == "energy";
        auto final_errors = [&](Problem& pb, vector_type& x) -> std::vector<RealType> {
            if constexpr (std::is_same_v<Problem, problem>) {
                if (energy_norms)
                    return pb.energy_errors(x, tf);
            }
            return all.select(all_errors(pb.errors(x, tf)));
        };
        if (prm.study.type != "time") {
            names = ricker ? std::vector<std::string>{} : energy_norms ? problem::energy_error_names() : all.names;
            if (prm.output.cpu) { names.push_back("cpu"); names.push_back("dofs"); names.push_back("face_dofs"); }
            std::vector<size_t> levels = study_levels(prm);
            for (size_t l : levels) {
                auto t0 = std::chrono::steady_clock::now();
                Problem pb(prm, l);
                size_t nt = study_steps(prm, l, pb.dd.h_min);
                vector_type x = observed && l == levels.back() ? integrate_observed(pb, prm, pb.initial_state(), nt, l)
                                                               : integrate(pb, prm, pb.initial_state(), nt, {}, l);
                RealType cpu = std::chrono::duration<RealType>(std::chrono::steady_clock::now() - t0).count();
                std::vector<RealType> e = ricker ? std::vector<RealType>{} : final_errors(pb, x);
                if (prm.output.cpu) {
                    e.push_back(cpu);
                    e.push_back(RealType(pb.assembler().LHS.rows()));     // unknowns before / after the static condensation
                    e.push_back(RealType(pb.n_face_dof()));               //   of the cell unknowns
                }
                rows.push_back({pb.dd.h_max, tf / nt, e, pb.dd.h_min});
                std::cout << bold << cyan << "      l = " << l << "  h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min
                          << "  nt = " << nt << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << reset << std::endl;
                print_table(names, rows, false, std::cout);
            }
            print_table(names, rows, false, log);
        }
        else {
            // against the exact solution only, or against a reference solution of the same integrator (mass norm of the
            // difference) with the exact errors of the velocity / scalar unknowns
            bool reference = prm.study.reference == "numerical";
            error_columns cols = reference ? columns(prm, {"elastic_v"}, {"acoustic_p"}) : all;
            names = energy_norms ? problem::energy_error_names() : cols.names;
            if (reference) names.insert(names.begin(), "mass_norm");
            Problem pb(prm, prm.mesh.level);
            size_t n0 = prm.time.n_steps;
            size_t n_levels = prm.study.time_levels;
            vector_type x_ref;
            if (reference) {
                size_t nt_ref = prm.study.reference_factor * n0 * (size_t(1) << (n_levels - 1));
                x_ref = integrate(pb, prm, pb.initial_state(), nt_ref);
                std::cout << bold << cyan << "      reference nt = " << nt_ref << reset << std::endl;
            }
            std::cout << bold << cyan << "      h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min
                      << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << reset << std::endl;
            for (size_t i = 0; i < n_levels; i++) {
                size_t nt = n0 * (size_t(1) << i);
                vector_type x = integrate(pb, prm, pb.initial_state(), nt);
                auto e = pb.errors(x, tf);
                std::vector<RealType> errors = energy_norms ? final_errors(pb, x) : reference ? cols.select({e[1], e[3]}) : cols.select(all_errors(e));
                if (reference) errors.insert(errors.begin(), mass_norm(pb.assembler().MASS, x - x_ref));
                rows.push_back({pb.dd.h_max, tf / nt, errors, pb.dd.h_min});
                print_table(names, rows, true, std::cout);
            }
            print_table(names, rows, true, log);
        }
    }

}

void HHOFirstOrder(const wave_parameters& prm);

void HHOFirstOrder(const wave_parameters& prm) {

    using namespace hho_first_order;

    std::cout << std::endl << bold << red << "   FIRST-ORDER ELASTO-ACOUSTIC COUPLING, " << prm.scheme_label()
              << ", " << prm.space_label() << reset << std::endl;
    std::string fname = prm.output_file();
    std::ofstream log(fname);
    log << "first-order formulation" << std::endl;
    prm.print(log);

    // errors at the final time: elastic velocity and stress, acoustic scalar (time derivative of the potential) and
    // flux (gradient of the potential), L2 norms on the cell unknowns
    std::vector<conv_row> rows;
    std::vector<std::string> names;
    if (prm.study.type == "cfl")
        cfl_study(prm, log, rows, names);
    else if (prm.study.type == "spectrum")
        spectrum_study(prm, log, rows, names);
    else if (prm.space.elastic_formulation == "three_fields")
        study<problem_three_fields>(prm, log, rows, names);
    else
        study<problem>(prm, log, rows, names);
    write_table(prm.table_file(), names, rows);
    std::cout << bold << red << "   --> " << fname << ", " << prm.table_file() << reset << std::endl << std::endl;
}

#endif /* HHOFirstOrder_hpp */
