//
//  LI_SDIRK_LTS_conv_test.hpp
//
//  Empirical TEMPORAL convergence test for LI-SDIRK(4,3): fixed mesh
//  (space discretization held constant, k=3 by default), sweeping ONLY
//  the macro time step dt (halved 4 times), measuring the L2 error
//  against the manufactured solution at final time T=1, and reporting
//  the observed convergence rate (log2 of successive error ratios).
//
//  Runs three cases with the SAME mesh/matrices, differing only in the
//  fine-region projector P:
//    "MIXED" : the real interface-based split from assembler.assemble_P
//              (the actual LI-SDIRK(4,3) scheme).
//    "P=0"   : Pfine == 0 everywhere -- should degenerate to the coarse
//              ERK4 Taylor precompute alone driving the whole domain
//              (dirk_weight_LTS_fine's implicit solve becomes trivial,
//              A = I, k_i = w(tau_i) for all stages).
//    "P=I"   : Pfine == Identity everywhere (and correspondingly
//              Pcoarse == 0, so erk_weight_LTS_coarse_global returns
//              w == 0) -- should degenerate to a standalone whole-domain
//              SDIRK(3,4) scheme, order exactly 4.
//
//  Mirrors prototypes/LTS/ERK4_LTS_conv_test.hpp's mesh/assembly/BC/
//  manufactured-solution setup (copied essentially verbatim -- this
//  prototype is about the TIME INTEGRATOR, not the spatial problem), but
//  builds the mesh and the erk_coupling_hho_scheme ONCE and sweeps dt
//  (not the mesh refinement level l).
//
//  Uses erk_weight_LTS_coarse_global's Pcoarse-taking overload (an
//  EXISTING, unmodified "stateless twin" of erk_weight_LTS_coarse already
//  present in erk_coupling_hho_scheme.hpp -- see its own header comment
//  there) instead of erk_weight_LTS_coarse itself, because this driver
//  deliberately swaps Pcoarse across the three P-cases on the SAME
//  erk_coupling_hho_scheme object; erk_weight_LTS_coarse caches its
//  active-dof list on first use and is only safe for a single fixed
//  Pcoarse for the object's lifetime (documented in its own file, see
//  the comment block directly above erk_weight_LTS_coarse_global). The
//  primary single-run driver (LI_SDIRK_LTS.hpp) uses erk_weight_LTS_coarse
//  itself, unchanged, exactly as ERK4_LTS.hpp does, since it only ever
//  uses one fixed Pcoarse for its whole run.
//
//  Created by Claude (Sonnet 5), for Romain Mottier -- 2026-09-14.
//

#include "dirk_coupling_hho_scheme.hpp"
#include <regex>

void LI_SDIRK_LTS_conv_test(int argc, char **argv);

namespace li_sdirk_conv_detail {

    // Pulls the two "L2-norm error = <value>" occurrences (elastic, then
    // acoustic) out of the text the postprocessor writes, and returns the
    // combined L2 error sqrt(e_L2^2 + a_L2^2).
    inline double combined_L2_error(const std::string & log_text) {
        std::regex re("L2-norm error = ([0-9eE+\\-\\.naninfNANIF]+)");
        auto begin = std::sregex_iterator(log_text.begin(), log_text.end(), re);
        auto end   = std::sregex_iterator();
        std::vector<double> vals;
        for (auto it = begin; it != end; ++it) {
            vals.push_back(std::stod((*it)[1].str()));
        }
        if (vals.size() < 2) {
            std::cout << bold << red << "   [LI-SDIRK conv test] WARNING: could not parse two L2 errors from log." << reset << std::endl;
            return -1.0;
        }
        return std::sqrt(vals[0]*vals[0] + vals[1]*vals[1]);
    }

}

void LI_SDIRK_LTS_conv_test(int argc, char **argv){

    std::cout << std::endl << bold << red << "   LI-SDIRK(4,3) LTS -- TEMPORAL CONVERGENCE TEST" << reset << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;

    // #############################################################################################
    // ############################## Mesh generation (ONCE, fixed) ####################################
    // #############################################################################################

    mesh_type msh;
    bool local_refinement = true;
    {
        RealType lx = 2.0, ly = 1.0;
        size_t nx = 4, ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-1.0, 0.0);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        if (local_refinement) {
            cartesian_2d_mesh_builder<RealType> mesh_builder2(lx,ly,nx,ny);
            mesh_builder2.refine_mesh(sim_data.m_n_divs);
            mesh_builder2.set_translation_data(-1.0, 0.0);
            mesh_builder2.build_mesh();
            typename mesh_type::point_type pt1(-0.5,0.5);
            typename mesh_type::point_type pt2(0.5,0.5);
            std::set<size_t> cell_indexes1 = postprocessor<mesh_type>::find_cells(pt1, msh, true);
            std::set<size_t> cell_indexes2 = postprocessor<mesh_type>::find_cells(pt2, msh, true);
            std::vector<size_t> vec;
            vec.insert(vec.end(), cell_indexes1.begin(), cell_indexes1.end());
            vec.insert(vec.end(), cell_indexes2.begin(), cell_indexes2.end());
            mesh_builder2.refine_cells(vec, sim_data.m_substeps_Q);
            mesh_builder2.move_to_mesh_storage(msh);
        }
    }

    RealType h_max = 1e-5, h_min = 10;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min) h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    int p = static_cast<int>(std::round(h_max / h_min));
    auto h_c = (p == 1) ? 1.25*h_max : 0.75*h_max;
    std::cout << bold << cyan << "      h_max = " << h_max << ", h_min = " << h_min << ", h_max/h_min = " << p << reset << std::endl;

    // #############################################################################################
    // ############################## Manufactured solution / HHO / materials / BCs ####################
    // #############################################################################################

    RealType ti = 0.0, tf = 1.0, t = ti;
    scal_vec_analytic_functions functions;
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
    auto u_fun    = functions.Evaluate_u(t);
    auto v_fun    = functions.Evaluate_v(t);
    auto f_fun    = functions.Evaluate_f(t);
    auto flux_fun = functions.Evaluate_sigma(t);
    auto s_u_fun    = functions.Evaluate_s_u(t);
    auto s_v_fun    = functions.Evaluate_s_v(t);
    auto s_f_fun    = functions.Evaluate_s_f(t);
    auto s_flux_fun = functions.Evaluate_s_q(t);

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
        return elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0);
    };
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };

    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    RealType eps = 1.0e-10;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
        mesh_type::point_type bar = barycenter(msh, *face_it);
        auto fc_id = msh.lookup(*face_it);
        if (std::fabs(bar.x()) < eps) interface_face_indexes.insert(fc_id);
    }
    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        if (bar.x() > 0) a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
        else e_material.insert(std::make_pair(cell_ind, elastic_mat_fun(bar)));
        for (auto face : faces(msh,cell)) {
            auto fc_id = msh.lookup(face);
            if (interface_face_indexes.count(fc_id)) {
                if (bar.x() > 0) interface_cell_pair_indexes[fc_id].second = cell_ind;
                else interface_cell_pair_indexes[fc_id].first = cell_ind;
            }
        }
    }
    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        mesh_type::point_type bar = barycenter(msh, *face_it);
        auto fc_id = msh.lookup(*face_it);
        if (interface_face_indexes.count(fc_id)) {
            if (bar.y() > 0) acoustic_internal_faces.insert(fc_id);
            else elastic_internal_faces.insert(fc_id);
        }
    }
    size_t bc_elastic_id = 0, bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
        mesh_type::point_type bar = barycenter(msh, *face_it);
        auto fc_id = msh.lookup(*face_it);
        if (bar.x() > 0) {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_acoustic_id, true};
            acoustic_bc_face_indexes.insert(fc_id);
        } else {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_elastic_id, true};
            elastic_bc_face_indexes.insert(fc_id);
        }
    }
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);

    // #############################################################################################
    // ###################################### Assembly (ONCE) ##########################################
    // #############################################################################################

    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_coupling_stabilization();
    if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);

    Matrix<RealType, Dynamic, 1> x_dof0;
    assembler.project_over_cells(msh, x_dof0, v_fun, flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof0, v_fun, s_v_fun);

    assembler.assemble(msh, f_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING;

    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();

    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    erk_an.refresh_faces_unknowns(x_dof0);

    size_t n_dof = assembler.RHS.rows();

    Matrix<double, Dynamic, Dynamic> a_tab;
    Matrix<double, Dynamic, 1> b_tab;
    Matrix<double, Dynamic, 1> c_tab;
    dirk_butcher_tableau::sdirk_tables(3, a_tab, b_tab, c_tab);
    RealType gamma = a_tab(0,0);

    assembler.assemble_P(msh, h_c, 1);
    SparseMatrix<RealType> Pfine_mixed   = assembler.Pfine;
    SparseMatrix<RealType> Pcoarse_mixed = assembler.Pcoarse;

    SparseMatrix<RealType> Pfine_zero(n_dof, n_dof);   // all-zero
    SparseMatrix<RealType> Pcoarse_full(n_dof, n_dof);
    Pcoarse_full.setIdentity();

    SparseMatrix<RealType> Pfine_full(n_dof, n_dof);
    Pfine_full.setIdentity();
    SparseMatrix<RealType> Pcoarse_zero(n_dof, n_dof); // all-zero

    // NOTE on per-case nt0: the coarse-region "B" operator (erk_weight /
    // the Taylor-freeze recursion inside erk_weight_LTS_coarse_global) is
    // an EXPLICIT stencil -- same as in the production ERK4_LTS.hpp driver,
    // where the macro dt is already tuned (there via -n, typically several
    // hundred steps) to respect the coarse region's own CFL limit
    // (dt < ~C*h_coarse/v_max). Whenever P selects real coarse cells
    // (MIXED) or the WHOLE domain is treated as coarse (P=0), dt must stay
    // inside that explicit stability limit or the scheme diverges outright
    // -- this is a property of the (already-validated) explicit coarse
    // precompute itself, not something introduced by this prototype. The
    // implicit SDIRK-only case (P=I) has no such restriction (A-stable),
    // so it can use a much coarser starting nt0 and still show the full
    // asymptotic order-4 range before hitting the spatial-error floor.
    // Explicit per-case nt lists, calibrated by direct probing of THIS
    // mesh/problem (see report): the coarse-region "B" operator is an
    // EXPLICIT stencil (same as the already-validated production
    // ERK4_LTS.hpp driver, where dt is likewise tuned via -n to respect
    // the coarse region's own CFL limit dt < ~C*h_coarse/v_max). Whenever
    // P selects real coarse cells (MIXED) or the WHOLE domain is treated
    // as coarse (P=0), dt must stay inside that explicit stability limit
    // or the scheme diverges outright -- a property of the explicit
    // coarse precompute itself, not introduced by this prototype. P=I
    // (pure implicit SDIRK, A-stable) has no such restriction and can
    // use much larger dt.
    //   MIXED : stable from nt>=~230 on this mesh; window before the
    //           k=3 spatial-error floor (~5.7e-6) is nt in [250,400].
    //   P=0   : whole-domain coarse CFL only stabilizes around nt>=~400,
    //           by which point the temporal error is ALREADY at/below
    //           the spatial floor -- no clean pre-floor window exists at
    //           this mesh resolution (reported as-is, see write-up).
    //   P=I   : A-stable, clean order-4 window at nt in [10,80].
    struct PCase { std::string label; SparseMatrix<RealType> Pfine; SparseMatrix<RealType> Pcoarse; std::vector<size_t> nt_list; };
    std::vector<PCase> cases = {
        {"MIXED (real interface split)", Pfine_mixed, Pcoarse_mixed, {250, 275, 300, 325, 350, 375, 400}},
        {"P=0  (pure coarse precompute)", Pfine_zero,  Pcoarse_full, {420, 460, 500, 700, 900, 1200}},
        {"P=I  (pure implicit SDIRK)   ", Pfine_full,  Pcoarse_zero, {10, 20, 40, 80, 160}},
    };

    for (auto & pc : cases) {

        std::cout << std::endl << bold << red << "   ================================================================" << reset << std::endl;
        std::cout << bold << red << "   LI-SDIRK(4,3) TEMPORAL CONVERGENCE -- case: " << pc.label << reset << std::endl;
        std::cout << bold << red << "   ================================================================" << reset << std::endl;

        std::vector<RealType> dts, errs;

        for (size_t nt : pc.nt_list) {

            RealType dt = (tf - ti) / RealType(nt);

            dirk_coupling_hho_scheme<RealType> dirk_op;
            dirk_op.compose(erk_an, pc.Pfine, gamma, dt);

            Matrix<RealType, Dynamic, 1> x_dof = x_dof0;

            auto eval_F = [&](RealType t_abs) -> Matrix<RealType, Dynamic, 1> {
                t = t_abs;
                auto v_fun_t   = functions.Evaluate_v(t);
                auto f_fun_t   = functions.Evaluate_f(t);
                auto s_v_fun_t = functions.Evaluate_s_v(t);
                auto s_f_fun_t = functions.Evaluate_s_f(t);
                assembler.get_e_bc_conditions().updateDirichletFunction(v_fun_t, 0);
                assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
                assembler.assemble_rhs(msh, f_fun_t, s_f_fun_t, true);
                return assembler.RHS;
            };
            std::function<Matrix<RealType, Dynamic, 1>(RealType)> eval_F_std = eval_F;

            for (size_t it = 1; it <= nt; it++) {

                RealType tn = dt*(it-1) + ti;
                auto x_dof_n = x_dof;

                std::vector<Matrix<RealType, Dynamic, 1>> w(4);
                for (int i = 0; i < 4; ++i) { w[i].resize(x_dof.rows()); w[i].setZero(); }

                erk_an.ZeroFc();
                RealType tn12 = tn + 0.5*dt, tn1 = tn + dt;
                Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
                Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn12);
                Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn1);
                erk_an.ZeroFc();
                erk_an.erk_weight_LTS_coarse_global(x_dof_n, pc.Pcoarse, w, Fn, Fn12, Fn1, dt);

                dirk_weight_LTS_fine<RealType>(erk_an, dirk_op, x_dof_n, pc.Pfine, w, tn, dt, eval_F_std, a_tab, b_tab, c_tab, gamma);

                x_dof = x_dof_n;
                t = tn + dt;
            }

            t = tf;
            auto v_fun_T      = functions.Evaluate_v(t);
            auto flux_fun_T   = functions.Evaluate_sigma(t);
            auto s_v_fun_T    = functions.Evaluate_s_v(t);
            auto s_flux_fun_T = functions.Evaluate_s_q(t);

            std::ostringstream err_log;
            postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler, x_dof, v_fun_T, flux_fun_T, s_v_fun_T, s_flux_fun_T, err_log);

            double err = li_sdirk_conv_detail::combined_L2_error(err_log.str());
            dts.push_back(dt);
            errs.push_back(err);

            std::cout << bold << cyan << "      nt = " << std::setw(6) << nt
                      << "   dt = " << std::setprecision(8) << dt
                      << "   L2 error = " << std::setprecision(10) << err << reset << std::endl;
        }

        std::cout << std::endl << bold << yellow << "   ---- Convergence table (" << pc.label << ") ----" << reset << std::endl;
        std::cout << std::setw(14) << "dt" << std::setw(20) << "L2 error" << std::setw(14) << "order" << std::endl;
        for (size_t i = 0; i < dts.size(); i++) {
            std::cout << std::setw(14) << std::setprecision(6) << dts[i]
                      << std::setw(20) << std::setprecision(10) << errs[i];
            if (i > 0 && errs[i] > 0 && errs[i-1] > 0) {
                // General log-ratio order estimate, robust to non-power-of-2
                // nt spacing: order = ln(e_prev/e_cur) / ln(dt_prev/dt_cur)
                double order = std::log(errs[i-1]/errs[i]) / std::log(dts[i-1]/dts[i]);
                std::cout << std::setw(14) << std::setprecision(4) << order;
            } else {
                std::cout << std::setw(14) << "--";
            }
            std::cout << std::endl;
        }

        // Least-squares log-log slope over ALL reported points: a single,
        // more robust "observed order" summary complementing the
        // pointwise pairwise orders above.
        {
            size_t N = dts.size();
            double mx = 0, my = 0;
            for (size_t i = 0; i < N; i++) { mx += std::log(dts[i]); my += std::log(errs[i]); }
            mx /= N; my /= N;
            double sxy = 0, sxx = 0;
            for (size_t i = 0; i < N; i++) {
                double dx = std::log(dts[i]) - mx, dy = std::log(errs[i]) - my;
                sxy += dx*dy; sxx += dx*dx;
            }
            double slope = sxy / sxx;
            std::cout << bold << green << "   Least-squares log-log slope (observed order, all points) = " << slope << reset << std::endl;
        }
    }

    std::cout << std::endl << bold << red << "   LI-SDIRK(4,3) TEMPORAL CONVERGENCE TEST -- DONE" << reset << std::endl << std::endl;
}
