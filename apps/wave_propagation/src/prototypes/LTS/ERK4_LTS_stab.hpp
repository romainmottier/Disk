
//  Created by Romain Mottier
// WITHOUT LOCAL REFINEMENT: ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l4 -n220 -p0 -f1 -e0
// WITH LOCAL REFINEMENT:    ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l4 -n220 -p3 -f1 -e0

void ERK4_LTS_stab(int argc, char **argv);

void ERK4_LTS_stab(int argc, char **argv){

    std::cout << std::endl << bold << red << "   ERK4-LTS STABILITY ANALYSIS" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    cpu.tic();

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
    
    if (sim_data.m_polygonal_mesh_Q) {       
        size_t l = sim_data.m_n_divs;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        {   // Simplicial meshes
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l0_conv_test_1.0.txt");    // l = 0
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l1_conv_test_0.35.txt");   // l = 1 
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l2_conv_test_0.15.txt");   // l = 2
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l3_conv_test_0.07.txt");   // l = 3 
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l4_conv_test_0.035.txt");  // l = 4
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l5_conv_test_0.026.txt");  // l = 5 
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l6_conv_test_0.017.txt");  // l = 6
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l7_conv_test_0.0125.txt"); // l = 7 
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l8_conv_test_0.0085.txt"); // l = 8
            // mesh_files.push_back("../../../../../meshes/conv_test/simplices/unstructured/l9_conv_test_0.005.txt");  // l = 9 
        }  
        {   // Polyhedral meshes
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_32.txt");     // -l 0
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_64.txt");     // -l 1
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_128.txt");    // -l 2
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_256.txt");    // -l 3
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_512.txt");    // -l 4
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_1024.txt");   // -l 5 
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_2048.txt");   // -l 6
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_4096.txt");   // -l 7 
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_8192.txt");   // -l 8
            mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_16384.txt");  // -l 9
        }      
        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    else {
        RealType lx = 2.0;  
        RealType ly = 1.0;          
        size_t nx = 6;
        size_t ny = 3;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-0.5, -0.5);
        mesh_builder.build_mesh();
        std::vector<size_t> cells_to_refine = {7,10};
        mesh_builder.refine_cells(cells_to_refine, sim_data.m_substeps_Q);
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
    auto p = h_max/h_min;
    auto h_c = 0.75*h_max;
    if (p == 1) {
        h_c = 1.25*h_max;
    }
    std::cout << bold << cyan << "      h_max = " << h_max << reset << std::endl;
    std::cout << bold << cyan << "      h_min = " << h_min << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << p << reset << std::endl << std::endl;

    // ##################################################
    // ################################################## Time
    // ##################################################

    size_t   nt = sim_data.m_nt_divs;
    RealType ti = 0.0, tf = 0.25;
    RealType dt = (tf - ti) / nt;

    // ##################################################
    // ################################################## HHO
    // ##################################################

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    // ##################################################
    // ################################################## Materials
    // ##################################################

    auto elastic_mat_fun = [](const typename mesh_type::point_type&) -> elastic_material_data<RealType> {
        return elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0);
    };
    auto acoustic_mat_fun = [](const typename mesh_type::point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };

    // #############################################################################################
    // ############################## Manufactured solution ########################################
    // #############################################################################################
    
    scal_vec_analytic_functions functions;
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
    
    // Elastic analytical functions
    auto u_fun    = functions.Evaluate_u(ti);
    auto v_fun    = functions.Evaluate_v(ti);
    auto a_fun    = functions.Evaluate_a(ti);
    auto f_fun    = functions.Evaluate_f(ti);
    auto flux_fun = functions.Evaluate_sigma(ti);
    
    // Acoustic analytical functions
    auto s_u_fun    = functions.Evaluate_s_u(ti);
    auto s_v_fun    = functions.Evaluate_s_v(ti);
    auto s_a_fun    = functions.Evaluate_s_a(ti);
    auto s_f_fun    = functions.Evaluate_s_f(ti);
    auto s_flux_fun = functions.Evaluate_s_q(ti);

    // ##################################################
    // ################################################## Structure
    // ##################################################
    
    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    RealType eps = 1.0e-10;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.x()) < eps) {
            interface_face_indexes.insert(fc_id);
            continue;
        } 
    }
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        // Assigning the material properties
        if (bar.x() > 0) {
            acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
            a_material.insert(std::make_pair(cell_ind,material));
        }
        else {
            elastic_material_data<RealType> material = elastic_mat_fun(bar);
            e_material.insert(std::make_pair(cell_ind,material));
        }
        // Detection of faces on the interfaces
        auto cell_faces = faces(msh,cell);
        for (auto face :cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) {
                if (bar.x() > 0) 
                interface_cell_pair_indexes[fc_id].second = cell_ind;
                else 
                interface_cell_pair_indexes[fc_id].first = cell_ind;
            }
        }
    }
    // Internal faces structure 
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);      
        bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
        if (is_member_Q) {
            if (bar.y() > 0) {
                acoustic_internal_faces.insert(fc_id);
            }
            else {
                elastic_internal_faces.insert(fc_id);
            }
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
        }
        else {
            disk::boundary_descriptor bi{bc_elastic_id, true};
            msh.backend_storage()->boundary_info.at(fc_id) = bi;
            elastic_bc_face_indexes.insert(fc_id);
        }   
    }
    // Detect interface elastic - acoustic
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);
    

    // ##################################################
    // ################################################## Assembly
    // ##################################################

    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);

    // Pulse acoustique comme condition initiale (meme que simulation de base)
    auto v_fun_adi_acoustic = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type& pt) -> disk::static_vector<double,2> {
        double x = pt.x(), y = pt.y(), xc = 0.0, yc = 0.1, fc = 10.0;
        double lp = std::sqrt(1.0)/fc;
        double r  = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        double wave = 10.0 / std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI);
        return disk::static_vector<double,2>{wave*(x-xc), wave*(y-yc)};
    };

    // #############################################################################################
    // ###################### Projecting initial data ##############################################
    // #############################################################################################
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);
    

    assembler.assemble(msh, f_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING;

    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
        assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(),
        assembler.get_e_face_dof(),    assembler.get_a_face_dof());
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
                       assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
                       assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
                       assembler.get_e_compress(), assembler.get_a_compress(),
                       elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    erk_an.refresh_faces_unknowns(x_dof);
    assembler.assemble_P(msh, h_c);

    if (sim_data.m_render_silo_files_Q) {
        std::ostringstream sn;
        sn << "silo_stab_l_" << sim_data.m_n_divs << "_k_" << sim_data.m_k_degree << "_p_" << p << "_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(sn.str(), 0, msh, hho_di, x_dof, e_material, a_material, false, h_c);
    }

    // ##################################################
    // ################################################## Log
    // ##################################################

    std::ostringstream fname;
    fname << "stab_l_" << sim_data.m_n_divs << "_n_" << nt
          << "_k_" << sim_data.m_k_degree << "_p_" << p << ".txt";
    std::ofstream log(fname.str());
    sim_data.write_simulation_data(log);
    log << "nt=" << nt << " dt=" << dt << " p=" << p
        << " n_dof=" << x_dof.rows() << "\n";
    log.flush();

    // ##################################################
    // ################################################## Stability sweep
    // ##################################################

    // C est de taille n_c x n_c.
    // Les faces ne sont pas des inconnues independantes — elles sont
    // determinees par les cellules via condensation statique dans erk_weight.
    //
    // Pour construire la colonne i de C :
    //   - e_i = vecteur de taille n_dof, nul partout sauf position i (cellule)
    //   - les faces de e_i restent a zero : erk_weight les recalcule en interne
    //   - on applique un pas LTS-RK4 exactement comme dans la simulation
    //   - on extrait les n_c premieres composantes (cellules) du resultat

    const int n_dof = static_cast<int>(x_dof.rows());
    const int n_c   = static_cast<int>(erk_an.n_c_dof());
    const int n_f   = n_dof - n_c;

    std::cout << bold << cyan << "      n_dof=" << n_dof << "  n_c=" << n_c << "  n_f=" << n_f
              << "  dt=" << dt << reset << std::endl;

    const double dt_min = dt;
    const double dt_max = 10.0 * dt;
    const int    n_pts  = 100;
    const double ddt    = (dt_max - dt_min) / static_cast<double>(n_pts - 1);

    std::cout << bold << red << "\n   STABILITY SWEEP (p=" << p << ")" << reset << std::endl;
    log << std::setw(20) << "dt" << std::setw(22) << "rho" << std::setw(10) << "stable\n";

    erk_an.ZeroFc();
    double dt_max_stable = -1.0;

    for (int s = 0; s < n_pts; ++s) {

        const double dt_s   = dt_min + s * ddt;
        const double dtau_s = dt_s / static_cast<double>(p);

        // Construction de C (n_c x n_c)
        tc.tic();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_c, n_c);

        for (int i = 0; i < n_c; ++i) {

            // Vecteur de base cellule i — faces restent a zero
            // erk_weight recalcule les faces via condensation statique
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            e_i(i) = 1.0;

            // Phase 1 : predicateur coarse — identique a la simulation
            std::vector<Matrix<RealType, Dynamic, 1>> w(4);
            for (int j = 0; j < 4; ++j) { w[j].resize(n_dof); w[j].setZero(); }
            erk_an.erk_weight_LTS_coarse_old(e_i, assembler.Pcoarse, w);

            // Phase 2 : p substeps RK4 fins — copie exacte de la boucle simulation
            Matrix<RealType, Dynamic, 1> x = e_i;
            for (int m = 0; m < p; ++m) {
                Matrix<RealType, Dynamic, 1> yn1, yn2, yn3, yn4, k1, k2, k3, k4;
                const double tm   =  m        * dtau_s;
                const double tm12 = (m + 0.5) * dtau_s;
                const double tm1  = (m + 1.0) * dtau_s;

                auto Taylor = [&](double tau) {
                    return w[0] + tau*w[1] + (tau*tau/2.0)*w[2] + (tau*tau*tau/6.0)*w[3];
                };

                yn1 = assembler.Pfine * x;
                erk_an.erk_weight(yn1, k1);
                k1 += Taylor(tm);

                yn2 = assembler.Pfine * (x + dtau_s/2.0*k1);
                erk_an.erk_weight(yn2, k2);
                k2 += Taylor(tm12);

                yn3 = assembler.Pfine * (x + dtau_s/2.0*k2);
                erk_an.erk_weight(yn3, k3);
                k3 += Taylor(tm12);

                yn4 = assembler.Pfine * (x + dtau_s*k3);
                erk_an.erk_weight(yn4, k4);
                k4 += Taylor(tm1);

                x += dtau_s/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
            }

            // Colonne i = partie cellules du resultat
            C.col(i) = x.head(n_c);
        }
        tc.toc();

        std::cout << bold << cyan << "      C(" << n_c << "x" << n_c << ") built in " << tc
                  << "   norm=" << std::setprecision(6) << C.norm() << reset << std::endl;

        // rho(C) exact via EigenSolver dense
        Eigen::EigenSolver<Eigen::MatrixXd> es(C);
        double rho = -1.0;
        if (es.info() == Eigen::Success) {
            rho = es.eigenvalues().cwiseAbs().maxCoeff();
            if (rho <= 1.0) dt_max_stable = dt_s;
        }
        if (rho <= 1.0) {
            dt_max_stable = dt_s;
            std::cout << "   --> dt_max_stable updated = " << dt_max_stable << std::endl;
        }
        std::cout << bold << cyan
                  << "      dt=" << std::setw(14) << std::setprecision(8)  << dt_s
                  << "   rho=" << std::setw(22) << std::setprecision(15) << rho;
        std::cout << (rho >= 0 && rho <= 1.0 ? "  [stable]" : "  [UNSTABLE]");
        std::cout << reset << std::endl;

        log << std::setw(20) << std::setprecision(10) << dt_s
            << std::setw(22) << std::setprecision(15) << rho
            << std::setw(10) << (rho >= 0 && rho <= 1.0 ? "yes" : "no") << "\n";

            

    }

    std::cout << bold << red << "\n   dt_max_stable (p=" << p << ") = "
              << std::setprecision(10) << dt_max_stable << reset << std::endl;
    log << "dt_max_stable=" << dt_max_stable << "\n";
    log.flush();

    cpu.toc();
    log << "CPU=" << cpu << "\n";
    std::cout << bold << red << "   CPU: " << cpu << reset << std::endl << std::endl;
}
