//
//  LI_SDIRK_LTS.hpp
//
//  LI-SDIRK(4,3): coarse region (I-P)U advanced by the EXISTING explicit
//  LTS-ERK "coarse precompute" (erk_weight_LTS_coarse, frozen degree-3
//  Taylor polynomial over the macro step dt, UNCHANGED); fine region P*U
//  advanced by a SINGLE implicit SDIRK(3,4) stage recursion of size dt
//  (dirk_weight_LTS_fine, new -- see dirk_coupling_hho_scheme.hpp), instead
//  of the p explicit RK4 sub-steps of size dt/p used by erk_weight_LTS_fine.
//
//  Mirrors prototypes/LTS/ERK4_LTS.hpp's overall structure (mesh, elasto-
//  acoustic assembler, BCs, manufactured solution -- copied essentially
//  verbatim, this prototype is about the TIME INTEGRATOR only), with the
//  p/dtau fine sub-stepping loop replaced by a 3-stage SDIRK call and NO
//  sub-stepping at all (single global dt everywhere).
//
//  Created by Claude (Sonnet 5), for Romain Mottier -- 2026-09-14.
//

#include "../../time_schemes/dirk_coupling_hho_scheme.hpp"

void LI_SDIRK_LTS(int argc, char **argv);

void LI_SDIRK_LTS(int argc, char **argv){

    // #############################################################################################
    // ############################## Simulation paramaters ########################################
    // #############################################################################################

    std::cout << std::endl << bold << red << "   LI-SDIRK(4,3) LTS TEST" << reset << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu, tcit;

    // #############################################################################################
    // ############################## Mesh generation ##############################################
    // #############################################################################################

    tc.tic();
    cpu.tic();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;
    bool local_refinement = true;

    {
        RealType lx = 2.0;
        RealType ly = 1.0;
        size_t nx = 4;
        size_t ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-1.0, 0.0);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        if (local_refinement) {
            RealType lx = 2.0;
            RealType ly = 1.0;
            size_t nx = 4;
            size_t ny = 2;
            cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
            mesh_builder.refine_mesh(sim_data.m_n_divs);
            mesh_builder.set_translation_data(-1.0, 0.0);
            mesh_builder.build_mesh();
            typename mesh_type::point_type pt1(-0.5,0.5);
            typename mesh_type::point_type pt2(0.5,0.5);
            std::set<size_t> cell_indexes1 = postprocessor<mesh_type>::find_cells(pt1, msh, true);
            std::set<size_t> cell_indexes2 = postprocessor<mesh_type>::find_cells(pt2, msh, true);
            std::vector<size_t> vec;
            vec.insert(vec.end(), cell_indexes1.begin(), cell_indexes1.end());
            vec.insert(vec.end(), cell_indexes2.begin(), cell_indexes2.end());
            auto n_loc_ref_lvl = sim_data.m_substeps_Q;
            mesh_builder.refine_cells(vec, n_loc_ref_lvl);
            mesh_builder.move_to_mesh_storage(msh);
        }
    }
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;
    RealType h_max = 1e-5;
    RealType h_min = 10;
    for (auto & cell : msh ) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min) h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    int p = static_cast<int>(std::round(h_max / h_min));
    auto h_c = 0.75*h_max;
    if (p == 1) h_c = 1.25*h_max;
    std::cout << bold << cyan << "      h_max = " << h_max << reset << std::endl;
    std::cout << bold << cyan << "      h_min = " << h_min << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << p << reset << std::endl << std::endl;

    // #############################################################################################
    // ################################ Time controls ##############################################
    // #############################################################################################

    size_t nt = std::max<size_t>(1, sim_data.m_nt_divs);
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = (tf-ti)/nt;
    RealType t = ti;

    // #############################################################################################
    // ############################## Manufactured solution #########################################
    // #############################################################################################

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

    // #############################################################################################
    // ################################## HHO setting ###############################################
    // #############################################################################################

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    // #############################################################################################
    // ################################ Material data ################################################
    // #############################################################################################

    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
        RealType rho = 1.0, vp = std::sqrt(3.0), vs = 1.0;
        return elastic_material_data<RealType>(rho,vp,vs);
    };
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        RealType rho = 1.0, vp = 1.0;
        return acoustic_material_data<RealType>(rho,vp);
    };

    // #############################################################################################
    // ############################## Boundary conditions ############################################
    // #############################################################################################

    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    RealType eps = 1.0e-10;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.x()) < eps) { interface_face_indexes.insert(fc_id); continue; }
    }
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        if (bar.x() > 0) a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
        else e_material.insert(std::make_pair(cell_ind, elastic_mat_fun(bar)));
        auto cell_faces = faces(msh,cell);
        for (auto face :cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) {
                if (bar.x() > 0) interface_cell_pair_indexes[fc_id].second = cell_ind;
                else interface_cell_pair_indexes[fc_id].first = cell_ind;
            }
        }
    }
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
        if (is_member_Q) {
            if (bar.y() > 0) acoustic_internal_faces.insert(fc_id);
            else elastic_internal_faces.insert(fc_id);
        }
    }

    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (bar.x() > 0) {
            disk::boundary_descriptor bi{bc_acoustic_id, true};
            msh.backend_storage()->boundary_info.at(fc_id) = bi;
            acoustic_bc_face_indexes.insert(fc_id);
        } else {
            disk::boundary_descriptor bi{bc_elastic_id, true};
            msh.backend_storage()->boundary_info.at(fc_id) = bi;
            elastic_bc_face_indexes.insert(fc_id);
        }
    }
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);

    // #############################################################################################
    // ###################################### Assembly ################################################
    // #############################################################################################

    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_coupling_stabilization();
    if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);

    // #############################################################################################
    // ###################### Projecting initial data #################################################
    // #############################################################################################

    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);

    // #############################################################################################
    // ###################################### Solving ##################################################
    // #############################################################################################

    assembler.assemble(msh, f_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING;

    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();

    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    erk_an.refresh_faces_unknowns(x_dof);

    // #############################################################################################
    // ############################### SDIRK(3,4) Butcher tableau ######################################
    // #############################################################################################

    Matrix<double, Dynamic, Dynamic> a_tab;
    Matrix<double, Dynamic, 1> b_tab;
    Matrix<double, Dynamic, 1> c_tab;
    dirk_butcher_tableau::sdirk_tables(3, a_tab, b_tab, c_tab);
    RealType gamma = a_tab(0,0);

    // #############################################################################################
    // ##################################### Fine-region projector #####################################
    // #############################################################################################

    assembler.assemble_P(msh, h_c, 1);

    dirk_coupling_hho_scheme<RealType> dirk_op;
    dirk_op.compose(erk_an, assembler.Pfine, gamma, dt);

    // #############################################################################################
    // ##################################### Preprocessor ###############################################
    // #############################################################################################

    std::ostringstream filename;
    filename << "li_sdirk_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << ".txt";
    std::ofstream simulation_log(filename.str());
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Space step = " << h_max << std::endl;
    simulation_log.flush();

    // #############################################################################################
    // ##################################### Time marching ###############################################
    // #############################################################################################

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

    std::cout << bold << red << "   TIME MARCHING SCHEME (LI-SDIRK(4,3)): " << reset << std::endl;
    for (size_t it = 1; it <= nt; it++) {

        tcit.tic();
        RealType tn = dt*(it-1) + ti;
        auto x_dof_n = x_dof;

        std::vector<Matrix<RealType, Dynamic, 1>> w(4);
        for (int i = 0; i < 4; ++i) { w[i].resize(x_dof.rows()); w[i].setZero(); }

        erk_an.ZeroFc();
        RealType tn12 = tn + 0.5*dt;
        RealType tn1  = tn + dt;
        Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
        Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn12);
        Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn1);
        erk_an.ZeroFc();
        erk_an.erk_weight_LTS_coarse(x_dof_n, assembler.Pcoarse, w, Fn, Fn12, Fn1, dt);

        std::function<Matrix<RealType, Dynamic, 1>(RealType)> eval_F_std = eval_F;
        dirk_weight_LTS_fine<RealType>(erk_an, dirk_op, x_dof_n, assembler.Pfine, w, tn, dt, eval_F_std, a_tab, b_tab, c_tab, gamma);

        x_dof = x_dof_n;
        t = tn + dt;
        tcit.toc();

        if (it == nt) {
            t = tn + dt;
            auto v_fun_T      = functions.Evaluate_v(t);
            auto flux_fun_T   = functions.Evaluate_sigma(t);
            auto s_v_fun_T    = functions.Evaluate_s_v(t);
            auto s_flux_fun_T = functions.Evaluate_s_q(t);
            std::cout << std::endl;
            postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler, x_dof, v_fun_T, flux_fun_T, s_v_fun_T, s_flux_fun_T, simulation_log);
            postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic_energy_norm(msh, hho_di, assembler, x_dof, v_fun_T, flux_fun_T, s_v_fun_T, s_flux_fun_T, simulation_log);
        }
    }

    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << std::endl << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;

}
