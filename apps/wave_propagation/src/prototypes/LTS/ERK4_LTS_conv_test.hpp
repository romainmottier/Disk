

//  Created by Romain Mottier

// ../wave_propagation -s0 -k3 -r0 -c0 -p0 -l4 -n300 -i0 -f0 -e0

void ERK4_LTS_conv_test(int argc, char **argv);

void ERK4_LTS_conv_test(int argc, char **argv){
  
    // #############################################################################################
    // ############################## Simulation paramaters ######################################## 
    // #############################################################################################
    
    std::cout << std::endl << bold << red << "   EXPLICIT ACOUSTIC CONV TEST" << reset << std::endl;
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu, tcit;
    
    // #############################################################################################
    // ############################## Mesh generation ##############################################
    // #############################################################################################
    
    cpu.tic();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;
    
    if (sim_data.m_polygonal_mesh_Q) {
        size_t l = sim_data.m_n_divs;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        bool use_poly_mesh = false; 
        bool use_simp_mesh = false; 
        if (use_poly_mesh) {
            for (int i = 0; i <= 9; ++i) {
                mesh_files.push_back("../../meshes/conv_test/poly/poly_" + std::to_string(32 * (1 << i)) + ".txt");
            }
        } 
        else if (use_simp_mesh) {
            std::vector<double> conv_vals = {1.0, 0.35, 0.15, 0.07, 0.035, 0.026, 0.017, 0.0125, 0.0085, 0.005};
            for (int i = 0; i < conv_vals.size(); ++i) {
                mesh_files.push_back(
                    "../../meshes/conv_test/simplices/unstructured/l" + std::to_string(i) + "_conv_test_" + std::to_string(conv_vals[i]) + ".txt");
            }
        }
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    else {
        RealType lx = 2.0;  
        RealType ly = 1.0;          
        size_t nx = 4;
        size_t ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-1.0, 0.0);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
    }
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;
    RealType h_max = 1e-5;
    RealType h_min = 10;
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min) {
            h_min = h_l;
        }
        else if (h_l > h_max) {
            h_max = h_l;
        }
    }
    auto h_c = 0.75*h_max;
    std::cout << bold << cyan << "      h_max = " << h_max << reset << std::endl;
    std::cout << bold << cyan << "      h_min = " << h_min << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << h_max/h_min << reset << std::endl << std::endl;
    
    // #############################################################################################
    // ################################ Time controls ##############################################
    // #############################################################################################
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) 
        nt  = sim_data.m_nt_divs;
    
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = (tf-ti)/nt;
    RealType t = ti;
    
    // #############################################################################################
    // ############################## Manufactured solution ########################################
    // #############################################################################################

    scal_vec_analytic_functions functions;
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionCubicInTimeAcoustic);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuarticInTimeAcoustic);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpaceAcoustic);
    
    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_matrix<double,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        disk::static_matrix<double, 2, 2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
        return sigma;
    };

    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
        disk::static_vector<double, 2> f{0,0};
        return f;
    };

    // Acoustic analytical functions
    auto s_u_fun    = functions.Evaluate_s_u(t);
    auto s_v_fun    = functions.Evaluate_s_v(t);
    auto s_a_fun    = functions.Evaluate_s_a(t);
    auto s_f_fun    = functions.Evaluate_s_f(t);
    auto s_flux_fun = functions.Evaluate_s_q(t);
    
    // #############################################################################################
    // ################################## HHO setting ##############################################
    // #############################################################################################

    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);
    
    // #############################################################################################
    // ################################ Material data ##############################################
    // #############################################################################################
    
    // Classify cells per material data and bc faces
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        RealType rho, vp, vs;
        rho = 1.0;            // Solid mass density
        vp  = std::sqrt(3.0); // Seismic compressional velocity vp
        vs  = 1.0;            // Seismic shear velocity vs
        elastic_material_data<RealType> material(rho,vp,vs);
        return material;
    };
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        RealType rho, vp;
        rho = 1.0; // Fluid mass density
        vp  = 1.0; // Seismic compressional velocity vp
        acoustic_material_data<RealType> material(rho,vp);
        return material;
    };
    
    // #############################################################################################
    // ############################## Boundary conditions ##########################################
    // #############################################################################################

    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;

    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        // Assigning the material properties
        acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
        a_material.insert(std::make_pair(cell_ind,material));
    }
    
    // Internal faces structure 
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        disk::boundary_descriptor bi{bc_acoustic_id, true};
        msh.backend_storage()->boundary_info.at(fc_id) = bi;
        acoustic_bc_face_indexes.insert(fc_id);
    }

    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);
    
    // #############################################################################################
    // ###################################### Assembly #############################################
    // #############################################################################################

    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_coupling_stabilization();
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }    
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);
    
    // #############################################################################################
    // ###################### Projecting initial data ##############################################
    // #############################################################################################

    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun);
    
    // #############################################################################################
    // ###################################### Solving ##############################################
    // #############################################################################################

    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    // ERK schemes
    int s = 4;
    erk_butcher_tableau::erk_tables(s, a, b, c);
    assembler.assemble(msh, null_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING; 
    assembler.assemble_P(msh, 1);

    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();
    
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);//assembler.get_interfaces());
    erk_an.refresh_faces_unknowns(x_dof);
    
    // ##################################################
    // ################################################## Preprocessor
    // ##################################################  
    
    std::ostringstream filename;
    filename << "explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
    std::string filename_str = filename.str();
    std::ofstream simulation_log(filename_str);
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Space step = " << h_max << std::endl;
    simulation_log.flush();
    std::cout << std::endl << std::endl;

    size_t it = 0;
    std::ostringstream filename_silo;
    filename_silo << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
    std::string silo_file_name = filename_silo.str();
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);

    // ##################################################
    // ################################################## Time marching
    // ##################################################
    
    assembler.assemble_P(msh, h_c);
    assembler.assemble_P_bis(msh, h_c);
    size_t nb_silo_files = 25;
    size_t step_interval = std::max(size_t(1), nt / nb_silo_files);
    std::cout << bold << red << "   TIME MARCHING SCHEME: " << reset << std::endl;
    auto p = std::pow(2, sim_data.m_substeps_Q);
    auto dtau = dt / p;
    for(size_t it = 1; it <= nt; it++) {

        //////////////////////////////////////////////////////////////////////////
        tcit.tic();
        RealType tn = dt*(it-1)+ti;
        if (it % step_interval == 0 || it == nt) {
            std::cout << bold << cyan << "      Time step number " << it << ": t = " << t << reset << std::endl;
        }

        ////////////////////////////////////////////////////////////////////////// PRECOMPUTATIONS: ERK ON THE GLOBAL DOFS 
        size_t n_dof = x_dof.rows();
        std::vector<Matrix<RealType, Dynamic, 1>> w(4), yn(4), k(4);
        for (int i = 0; i < 4; ++i) {
            w[i].resize(n_dof);  w[i].setZero();
            yn[i].resize(n_dof); yn[i].setZero();
            k[i].resize(n_dof);  k[i].setZero();
        }
        auto x_dof_n = x_dof;
        
        // Manufactured solution + BC at tn, tn+1/2, tn+1
        auto tn12 = tn + 0.5*dt;
        auto tn1 = tn + dt;
        auto v_fun_n = functions.Evaluate_v(tn);
        auto f_fun_n = functions.Evaluate_f(tn);
        auto s_v_fun_n = functions.Evaluate_s_v(tn);
        auto s_f_fun_n = functions.Evaluate_s_f(tn);
        assembler.get_e_bc_conditions().updateDirichletFunction(v_fun_n, 0);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_n, 0);
        assembler.assemble_rhs(msh, f_fun_n, s_f_fun_n, false);
        Matrix<RealType, Dynamic, 1> Fn = assembler.RHS;
        // Manufactured solution + BC at tn+1/2
        auto v_fun_n12 = functions.Evaluate_v(tn12);
        auto f_fun_n12 = functions.Evaluate_f(tn12);
        auto s_v_fun12 = functions.Evaluate_s_v(tn12);
        auto s_f_fun_n12 = functions.Evaluate_s_f(tn12);
        assembler.get_e_bc_conditions().updateDirichletFunction(v_fun_n12, 0);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun12, 0);
        assembler.assemble_rhs(msh, f_fun_n12, s_f_fun_n12, false);        
        Matrix<RealType, Dynamic, 1> Fn12 = assembler.RHS;
        // Manufactured solution + BC at tn+1
        auto v_fun_n1 = functions.Evaluate_v(tn1);
        auto f_fun_n1   = functions.Evaluate_f(tn1);
        auto s_v_fun_n1 = functions.Evaluate_s_v(tn1);
        auto s_f_fun_n1 = functions.Evaluate_s_f(tn1);
        assembler.get_e_bc_conditions().updateDirichletFunction(v_fun_n1, 0);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_n1, 0);
        assembler.assemble_rhs(msh, f_fun_n1, s_f_fun_n1, false);
        Matrix<RealType, Dynamic, 1> Fn1 = assembler.RHS;
        
        erk_an.compute_wi_with_F(x_dof_n, assembler.Pcoarse, w, Fn, Fn12, Fn1, dt);

        ////////////////////////////////////////////////////////////////////////// LOOP OVER THE SUBSTEPS: ERK4 ON THE LOCAL DOFS WITH INJECTION OF THE GLOBAL DOFS
        // Butcher tableau offsets for RK4: 0, 1/2, 1/2, 1
        const std::array<double, 4> c = {0.0, 0.5, 0.5, 1.0};
        // RK4 stage increments
        const std::array<double, 4> a = {0.0, 0.5, 0.5, 1.0};
        for (int m = 0; m < p; m++) {
            for (int s = 0; s < 4; ++s) {
                if (s == 0) {
                    yn[s] = assembler.Pfine * x_dof_n;
                }
                else {
                    yn[s] = assembler.Pfine * (x_dof_n + dtau * a[s] * k[s-1]);
                }            
                
                erk_an.erk_weight(yn[s], k[s]);   // erk weight
                double t  = (m + c[s]) * dtau;
                double t2 = t * t;
                double t3 = t * t2;
                k[s] += w[0] + t*w[1] + t2*w[2]/2.0 + t3*w[3]/6.0;
            }

            // FINAL UPDATE
            x_dof_n += dtau * (k[0] + 2.0*k[1] + 2.0*k[2] + k[3]) / 6.0;

        }
        x_dof = x_dof_n;
        t += dt;
        if (sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::ostringstream filename;
            filename << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << "_";
            std::string silo_file_name = filename.str();
            postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false, h_c);
        }
        tcit.toc();
        if (sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::cout << bold << yellow << "         Iteration completed in " << tcit << " seconds" << reset << std::endl;
        }
    }
    
    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << std::endl << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;

}