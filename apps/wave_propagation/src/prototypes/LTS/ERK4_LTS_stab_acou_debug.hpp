
//  Created by Romain Mottier
// Debug version of ERK4_LTS_stab_acou :
// Compares C_lts (erk_weight_LTS_fine) vs C_std (Butcher + erk_weight)
// to identify why rho(C_lts) != rho(C_std) for p=1.

void ERK4_LTS_stab_acou_debug(int argc, char **argv);

void ERK4_LTS_stab_acou_debug(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS STABILITY ANALYSIS ACOUSTIC [DEBUG]" << std::endl << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    cpu.tic();

    // =========================================================================
    // Mesh generation
    // =========================================================================
    tc.tic();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;

    mesh_type msh;

    if (sim_data.m_polygonal_mesh_Q)
    {
        size_t l = sim_data.m_n_divs;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_32.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_64.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_128.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_256.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_512.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_1024.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_2048.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_4096.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_8192.txt");
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_16384.txt");
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    else {
        RealType lx = 1.0, ly = 1.0;
        int    L  = sim_data.m_substeps_Q;
        size_t nx = 18, ny = 18;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.build_mesh();
        auto is_fine_zone = [](const typename mesh_type::point_type& p) {
            return p.x() > 1.0/3.0 && p.x() < 2.0/3.0
                && p.y() > 1.0/3.0 && p.y() < 2.0/3.0;
        };
        if (L > 0) mesh_builder.refine_with_protection(is_fine_zone, L);
        mesh_builder.move_to_mesh_storage(msh);
    }
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;

    RealType h_max = 1.0e-5, h_min = 10.0;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min)      h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    int p = static_cast<int>(std::round(h_max / h_min));
    auto h_c = (p == 1) ? 1.25 * h_max : 0.75 * h_max;

    std::cout << bold << cyan << "      h_max       = " << h_max << reset << std::endl;
    std::cout << bold << cyan << "      h_min       = " << h_min           << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << p     << reset << std::endl << std::endl;

    // =========================================================================
    // Time parameters
    // =========================================================================
    size_t   nt = sim_data.m_nt_divs;
    RealType ti = 0.0, tf = 0.25;
    RealType dt = (tf - ti) / nt;

    // =========================================================================
    // HHO discretisation degree
    // =========================================================================
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    // =========================================================================
    // Material / boundary
    // =========================================================================
    auto acoustic_mat_fun = [](const typename mesh_type::point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };
    auto null_s_fun = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> double {
        return 0.0;
    };
    auto null_fun = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> disk::static_vector<double,2> {
        return disk::static_vector<double,2>{0,0};
    };
    auto null_flux_fun = [](const typename disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> disk::static_matrix<double,2,2> {
        return disk::static_matrix<double,2,2>::Zero(2,2);
    };

    std::map<size_t, elastic_material_data<RealType>>  e_material;
    std::map<size_t, acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t, std::pair<size_t,size_t>> interface_cell_pair_indexes;

    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        a_material.insert({cell_ind, acoustic_mat_fun(barycenter(msh, cell))});
    }
    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        auto fc_id = msh.lookup(*face_it);
        if (interface_face_indexes.count(fc_id)) acoustic_internal_faces.insert(fc_id);
    }
    size_t bc_elastic_id = 0, bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto fc_id = msh.lookup(*face_it);
        msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_acoustic_id, true};
        acoustic_bc_face_indexes.insert(fc_id);
    }
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);

    // =========================================================================
    // HHO assembly
    // =========================================================================
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_coupling_stabilization();
    if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, null_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
    assembler.assemble(msh, null_fun, null_s_fun, true);
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
    assembler.assemble_P(msh, h_c, 1);

    // =========================================================================
    // Log file
    // =========================================================================
    std::ostringstream fname;
    fname << "debug_l_" << sim_data.m_n_divs << "_n_" << nt
          << "_k_" << sim_data.m_k_degree << "_p_" << p << ".txt";
    std::ofstream log(fname.str());
    sim_data.write_simulation_data(log);
    log << "nt=" << nt << " dt=" << dt << " p=" << p << " n_dof=" << x_dof.rows() << "\n";
    log.flush();

    const int n_dof = static_cast<int>(x_dof.rows());
    const int n_c   = static_cast<int>(erk_an.n_c_dof());
    const int n_f   = n_dof - n_c;

    std::cout << bold << red << "   DISCRETIZATION" << reset << std::endl;
    std::cout << bold << cyan
              << "      n_dof=" << n_dof << "  n_c=" << n_c << "  n_f=" << n_f
              << reset << std::endl;

    Matrix<RealType, Dynamic, 1> F_zero = Matrix<RealType, Dynamic, 1>::Zero(n_dof);

    // Butcher tableau RK4 classique
    Matrix<RealType, Dynamic, Dynamic> a_bt;
    Matrix<RealType, Dynamic, 1>       b_bt, c_bt;
    const int ss = 4;
    erk_butcher_tableau::erk_tables(ss, a_bt, b_bt, c_bt);

    // =========================================================================
    // DEBUG : compare C_lts vs C_std pour n_diag valeurs de dt
    // =========================================================================
    // C_lts : erk_weight_LTS_coarse (si p!=1) + p * erk_weight_LTS_fine
    // C_std : boucle Butcher RK4 + erk_weight (sans Pfine)
    //
    // Pour p=1, Pfine=I, w=0 : les deux doivent donner ||C_lts - C_std|| ~ 0
    // Si ce n'est pas le cas, le bug est dans erk_weight_LTS_fine pour p=1.
    // =========================================================================

    const int    n_diag      = 5;
    const double dt_min_diag = dt;
    const double dt_max_diag = 2.0 * dt;
    const double ddt_diag    = (dt_max_diag - dt_min_diag) / static_cast<double>(n_diag - 1);

    std::cout << bold << red << "\n   DEBUG : C_lts vs C_std" << reset << std::endl;
    log << "\n--- DEBUG : C_lts vs C_std ---\n";
    log << std::setw(12) << "dt"
        << std::setw(26) << "||C_lts-C_std||/||C_lts||"
        << std::setw(16) << "rho(C_lts)"
        << std::setw(16) << "rho(C_std)" << "\n";

    for (int id = 0; id < n_diag; ++id) {

        const double dt_d   = dt_min_diag + id * ddt_diag;
        const double dtau_d = dt_d / p;

        // ---- C_lts : chemin LTS ----
        Eigen::MatrixXd C_lts = Eigen::MatrixXd::Zero(n_c, n_c);
        for (int i = 0; i < n_c; ++i) {
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            e_i(i) = 1.0;
            erk_an.refresh_faces_unknowns(e_i);
            std::vector<Matrix<RealType, Dynamic, 1>> w(4);
            for (int j = 0; j < 4; ++j) { w[j].resize(n_dof); w[j].setZero(); }
            erk_an.ZeroFc();
            if (p != 1)
                erk_an.erk_weight_LTS_coarse(e_i, assembler.Pcoarse, w, F_zero, F_zero, F_zero, dt_d);
            for (int m = 0; m < p; ++m) {
                RealType tm = m * dtau_d;
                erk_an.erk_weight_LTS_fine(e_i, assembler.Pfine, w, F_zero, F_zero, F_zero, tm, dtau_d);
            }
            C_lts.col(i) = e_i.head(n_c);
        }

        // ---- C_std : Butcher + erk_weight ----
        Eigen::MatrixXd C_std = Eigen::MatrixXd::Zero(n_c, n_c);
        for (int i = 0; i < n_c; ++i) {
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            e_i(i) = 1.0;
            erk_an.refresh_faces_unknowns(e_i);
            Matrix<RealType, Dynamic, 1> y0 = e_i;
            Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, ss);
            Matrix<RealType, Dynamic, 1> yn, ki;
            erk_an.ZeroFc();
            for (int ii = 0; ii < ss; ++ii) {
                yn = y0;
                for (int j = 0; j < ss - 1; ++j)
                    yn += a_bt(ii,j) * dt_d * k.block(0, j, n_dof, 1);
                erk_an.erk_weight(yn, ki);
                e_i += dt_d * b_bt(ii,0) * ki;
                k.block(0, ii, n_dof, 1) = ki;
            }
            C_std.col(i) = e_i.head(n_c);
        }

        // ---- Comparaison ----
        double diff_norm  = (C_lts - C_std).norm();
        double C_lts_norm = C_lts.norm();
        double rel_diff   = diff_norm / (C_lts_norm + 1e-300);

        Eigen::EigenSolver<Eigen::MatrixXd> es_lts(C_lts), es_std(C_std);
        double rho_lts = (es_lts.info() == Eigen::Success)
                         ? es_lts.eigenvalues().cwiseAbs().maxCoeff() : -1.0;
        double rho_std = (es_std.info() == Eigen::Success)
                         ? es_std.eigenvalues().cwiseAbs().maxCoeff() : -1.0;

        std::cout << bold << cyan
                  << "      dt=" << std::setprecision(5) << dt_d
                  << "  ||C_lts-C_std||/||C_lts||=" << std::setprecision(4) << rel_diff
                  << "  rho_lts=" << std::setprecision(8) << rho_lts
                  << "  rho_std=" << std::setprecision(8) << rho_std
                  << reset << std::endl;

        log << std::setw(12) << std::setprecision(6)  << dt_d
            << std::setw(26) << std::setprecision(6)  << rel_diff
            << std::setw(16) << std::setprecision(10) << rho_lts
            << std::setw(16) << std::setprecision(10) << rho_std << "\n";
        log.flush();
    }

    cpu.toc();
    log << "CPU=" << cpu << "\n";
    std::cout << bold << red << "   CPU: " << cpu << reset << std::endl << std::endl;
}
