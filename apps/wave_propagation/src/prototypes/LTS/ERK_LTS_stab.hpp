
//  Created by Romain Mottier
// WITH LOCAL REFINEMENT:    ../../../wave_propagation -k0 -s0 -r0 -c0 -m0 -l0 -n11 -p1 -f1 -e0

void ERK_LTS_stab(int argc, char **argv);

void ERK_LTS_stab(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS STABILITY ANALYSIS ERK(s)" << std::endl << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    cpu.tic();
    const bool save_eigenvalues_Q = false;

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

        // Polyhedral meshes
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_32.txt");    // -l 0
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_64.txt");    // -l 1
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_128.txt");   // -l 2
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_256.txt");   // -l 3
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_512.txt");   // -l 4
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_1024.txt");  // -l 5
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_2048.txt");  // -l 6
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_4096.txt");  // -l 7
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_8192.txt");  // -l 8
        mesh_files.push_back("../../../../../meshes/conv_test/poly/poly_16384.txt"); // -l 9

        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    // else {
    //     RealType lx = 1.0;
    //     RealType ly = 1.0;
    //     int    L    = sim_data.m_substeps_Q;
    //     size_t nx   = 18;
    //     size_t ny   = 18;
    //     cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
    //     mesh_builder.build_mesh();
    //     auto is_fine_zone = [](const typename mesh_type::point_type& p) {
    //         return p.x() > 1.0/3.0 && p.x() < 2.0/3.0
    //         && p.y() > 1.0/3.0 && p.y() < 2.0/3.0;
    //     };
    //     if (L > 0) {
    //         mesh_builder.refine_with_protection(is_fine_zone, L);
    //     }
    //     mesh_builder.move_to_mesh_storage(msh);
    // }
    else {
        RealType lx = 1.0;
        RealType ly = 1.0;
        int    L    = sim_data.m_substeps_Q;
        size_t nx   = 18;
        size_t ny   = 18;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.build_mesh();
        auto is_fine_zone = [](const typename mesh_type::point_type& p) {
            return p.x() > 1.0/3.0 && p.x() < 2.0/3.0
            && p.y() > 1.0/3.0 && p.y() < 2.0/3.0;
        };
        if (L > 0) {
            mesh_builder.refine_with_protection(is_fine_zone, L);
        }
        mesh_builder.move_to_mesh_storage(msh);
    }
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;
    RealType h_max = 1.0e-5;
    RealType h_min = 10.0;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min) {
            h_min = h_l;
        }
        else if (h_l > h_max) {
            h_max = h_l;
        }
    }
    // auto p   = h_max / h_min;
    auto p   = 1;
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
    // Material parameters
    // =========================================================================
    auto elastic_mat_fun = [](const typename mesh_type::point_type&) -> elastic_material_data<RealType> {
        return elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0);
    };
    auto acoustic_mat_fun = [](const typename mesh_type::point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };

    auto null_s_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> double {
      return 0.0;
    }; 

    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
      disk::static_vector<double, 2> f{0,0};
      return f;
    };
    
    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_matrix<double,2,2> {
      double x,y;
      x = pt.x();
      y = pt.y();
      disk::static_matrix<double, 2, 2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
      return sigma;
    };

    // =========================================================================
    // Domain decomposition: elastic (x<0) / acoustic (x>0) / interface (x=0)
    // =========================================================================
    std::map<size_t, elastic_material_data<RealType>>  e_material;
    std::map<size_t, acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t, std::pair<size_t, size_t>> interface_cell_pair_indexes;

    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
    }
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id    = msh.lookup(face);
        bool is_iface = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
        if (is_iface) {
            acoustic_internal_faces.insert(fc_id);
        }
    }
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        disk::boundary_descriptor bi{bc_acoustic_id, true};
        msh.backend_storage()->boundary_info.at(fc_id) = bi;
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
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, null_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);

    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    int s = 4;
    erk_butcher_tableau::erk_tables(s, a, b, c);

    assembler.assemble(msh, null_fun, null_s_fun, true);
    assembler.LHS += assembler.COUPLING;
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(), assembler.get_e_face_dof(), assembler.get_a_face_dof());
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    erk_an.refresh_faces_unknowns(x_dof);
    assembler.assemble_P(msh, h_c, 1);

    if (sim_data.m_render_silo_files_Q) {
        std::ostringstream sn;
        sn << "silo_stab_l_" << sim_data.m_n_divs << "_k_" << sim_data.m_k_degree << "_p_" << p << "_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(sn.str(), 0, msh, hho_di, x_dof, e_material, a_material, false, h_c, assembler.cell_is_fine_silo);
    }
 
    // =========================================================================
    // Log file
    // =========================================================================
    std::ostringstream fname;
    fname << "stab_l_" << sim_data.m_n_divs << "_n_" << nt << "_k_" << sim_data.m_k_degree << "_p_" << p << ".txt";
    std::ofstream log(fname.str());
    sim_data.write_simulation_data(log);
    log << "nt=" << nt << " dt=" << dt << " p=" << p << " n_dof=" << x_dof.rows() << "\n";
    log.flush();

    // =========================================================================
    // Stability sweep
    // =========================================================================
    //
    // We build the amplification matrix C of size (n_c x n_c), where n_c is
    // the number of cell DOFs. Face unknowns are NOT independent — they are
    // determined by static condensation (18b): K_F U_F = -K_FT U_T.
    // Each column i of C is obtained by:
    //   1. Setting e_i = canonical basis vector on cell dofs; faces are set
    //      to their condensed values via refresh_faces_unknowns (admissible state).
    //   2. Applying one full LTS-RK4 step (coarse predictor + p fine substeps)
    //      with source term f = 0 (valid for stability analysis).
    //   3. Extracting the first n_c components (cell part) of the result.
    //
    // The scheme is stable for a given dt iff rho(C) <= 1.
    // The sweep scans dt in [dt_min, dt_max] with n_pts points and stops
    // as soon as rho > rho_stop.
    //
    // TEST 1 is run once at the first dt point to validate that C*x matches
    // a direct LTS step on a random admissible state (should be ~machine epsilon).
    //
    // If save_eigenvalues_Q is true, eigenvalue files are written for each dt
    // near the stability boundary for complex-plane post-processing.

    const int n_dof = static_cast<int>(x_dof.rows());
    const int n_c   = static_cast<int>(erk_an.n_c_dof());
    const int n_f   = n_dof - n_c;

    // Display tolerance: rho in (1, 1+rho_eps] is shown as boundary noise
    const double rho_eps  = 1.0e-5;

    // Sweep hard stop threshold
    const double rho_stop = 1.01;

    std::cout << bold << red << "   DISCRETIZATION" << reset << std::endl;
    std::cout << bold << cyan << "      n_dof=" << n_dof << "  n_c=" << n_c << "  n_f=" << n_f << reset << std::endl;
    std::cout << bold << cyan << "      rho_eps=" << rho_eps << "  rho_stop=" << rho_stop << reset << std::endl;
    std::cout << bold << cyan << "      save_eigenvalues_Q=" << (save_eigenvalues_Q ? "true" : "false") << reset << std::endl;

    const double dt_min = dt;
    const double dt_max = 2.0 * dt;
    const int    n_pts  = 50;
    const double ddt    = (dt_max - dt_min) / static_cast<double>(n_pts - 1);

    erk_an.ZeroFc();

    std::cout << bold << red << "\n   STABILITY SWEEP (p=" << p << ")" << reset << std::endl;
    log << std::setw(20) << "dt" << std::setw(22) << "rho" << std::setw(10) << "stable\n";

    double dt_max_stable  = -1.0;
    bool   test1_done_Q   = false;   // TEST 1 is run only once (first point)

    Matrix<RealType, Dynamic, 1> F_zero = Matrix<RealType, Dynamic, 1>::Zero(n_dof);

    for (int s = 0; s < n_pts; ++s) {

        const double dt_s   = dt_min + s * ddt;
        const double dtau_s = dt_s / p;

        std::cout << bold << cyan << "      dt = " << std::setprecision(5) << dt_s << ":" << reset;

        // -----------------------------------------------------------------
        // Build the amplification matrix C column by column
        // -----------------------------------------------------------------
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_c, n_c);

        for (int i = 0; i < n_c; ++i) {

            // Canonical basis vector: cell i = 1; faces set by static condensation
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            e_i(i) = 1.0;
            erk_an.refresh_faces_unknowns(e_i);

            Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<RealType, Dynamic, 1> Fg, Fg_c;            
            Matrix<RealType, Dynamic, 1> yn, ki;
            for (int ii = 0; ii < s; ii++) {
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(ii,j) * dt * k.block(0, j, n_dof, 1);
                }
                erk_an.erk_weight(yn, ki);
                // Accumulated solution
                e_i += dt*b(ii,0)*ki;
                k.block(0, ii, n_dof, 1) = ki;      
            }
            
            C.col(i) = e_i.head(n_c);
            // std::cout << C.col(i) << std::endl;
        }

        // -----------------------------------------------------------------
        // TEST 1 : run once at the first sweep point to validate C
        // C * x_rand must equal a direct LTS step on the same admissible x_rand
        // -----------------------------------------------------------------
        if (!test1_done_Q) {
            Matrix<RealType, Dynamic, 1> x_rand_c = Matrix<RealType, Dynamic, 1>::Random(n_c);
            Matrix<RealType, Dynamic, 1> Cx = C * x_rand_c;

            Matrix<RealType, Dynamic, 1> x_rand = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            x_rand.head(n_c) = x_rand_c;
            erk_an.refresh_faces_unknowns(x_rand);  

            Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<RealType, Dynamic, 1> Fg, Fg_c;            
            Matrix<RealType, Dynamic, 1> yn, ki;
            for (int ii = 0; ii < s; ii++) {
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(ii,j) * dt * k.block(0, j, n_dof, 1);
                }
                erk_an.erk_weight(yn, ki);
                x_rand += dt*b(ii,0)*ki;
                k.block(0, ii, n_dof, 1) = ki;      
            }

            double err1 = (Cx - x_rand.head(n_c)).norm();
            std::cout << bold << yellow << "   ||C*x - step(x)|| = " << std::setprecision(6) << err1 << reset;
            test1_done_Q = true;
        }

        // -----------------------------------------------------------------
        // Compute spectral radius rho(C) = max |lambda_i(C)|
        // -----------------------------------------------------------------
        Eigen::EigenSolver<Eigen::MatrixXd> es(C);
        double rho = -1.0;

        if (es.info() == Eigen::Success) {
            rho = es.eigenvalues().cwiseAbs().maxCoeff();

            // Save eigenvalue files only if requested and near/above stability boundary
            if (save_eigenvalues_Q && rho > 1.0 - 0.05) {
                std::ostringstream ev_fname;
                ev_fname << "eigenvalues_dt_" << std::setprecision(6) << dt_s << ".txt";
                std::ofstream ev_file(ev_fname.str());
                ev_file << "# dt=" << dt_s << "  rho=" << rho << "\n";
                ev_file << "# real  imag\n";
                for (int j = 0; j < es.eigenvalues().size(); ++j) {
                    ev_file << std::setprecision(15)
                            << es.eigenvalues()(j).real() << "  "
                            << es.eigenvalues()(j).imag() << "\n";
                }
                ev_file.close();
            }
        }
        else {
            std::cout << bold << red << "   --> EigenSolver FAILED at dt=" << dt_s << reset << std::endl;
        }

        // Display: cyan = stable (rho <= 1 + rho_eps), red = unstable
        bool is_stable = (rho >= 0.0 && rho <= 1.0 + rho_eps);
        if (is_stable) {
            std::cout << bold << cyan;
            dt_max_stable = dt_s;
        }
        else {
            std::cout << bold << red;
        }

        std::cout << "   rho = " << std::setprecision(10) << rho << reset << std::endl;

        log << std::setw(20) << std::setprecision(10) << dt_s
            << std::setw(22) << std::setprecision(15) << rho
            << std::setw(10) << (is_stable ? "yes" : "no") << "\n";
        log.flush();

        // -----------------------------------------------------------------
        // If unstable beyond threshold: TEST 2 — propagate the dominant
        // eigenvector and confirm divergence, then stop the sweep
        // -----------------------------------------------------------------
        if (rho > rho_stop) {

            int idx_max = 0;
            es.eigenvalues().cwiseAbs().maxCoeff(&idx_max);

            // Use real part of dominant eigenvector; fall back to imag if trivial
            Matrix<RealType, Dynamic, 1> v_unstable = es.eigenvectors().col(idx_max).real();
            if (v_unstable.norm() < 1.0e-10)
                v_unstable = es.eigenvectors().col(idx_max).imag();
            v_unstable.normalize();

            // Build admissible initial state from eigenvector (cell part only)
            Matrix<RealType, Dynamic, 1> x_eig = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            x_eig.head(n_c) = v_unstable;
            erk_an.refresh_faces_unknowns(x_eig);   // ← admissible state

            const int N_check = 200;
            std::cout << bold << red << "\n   TEST 2 : propagation du vecteur propre instable\n" << reset;

            for (int step = 0; step < N_check; ++step) {
                std::vector<Matrix<RealType, Dynamic, 1>> w_eig(4);
                for (int j = 0; j < 4; ++j) { w_eig[j].resize(n_dof); w_eig[j].setZero(); }

                erk_an.ZeroFc();
                if (p != 1) {
                    erk_an.erk_weight_LTS_coarse(x_eig, assembler.Pcoarse, w_eig, F_zero, F_zero, F_zero, dt_s);
                }
                for (int m = 0; m < p; ++m) {
                    RealType tm = m * dtau_s;
                    erk_an.erk_weight_LTS_fine(x_eig, assembler.Pfine, w_eig, F_zero, F_zero, F_zero, tm, dtau_s);
                }

                double norm_eig = x_eig.head(n_c).norm();
                if (step % 20 == 0) {
                    std::cout << bold << cyan
                              << "      step " << std::setw(3) << step
                              << "   ||x_c|| = " << std::setprecision(8) << norm_eig
                              << reset << std::endl;
                }
                if (norm_eig > 1.0e10 || std::isnan(norm_eig)) {
                    std::cout << bold << red << "      --> DIVERGENCE confirmée" << reset << std::endl;
                    break;
                }
            }

            break;   // stop sweep
        }

    }

    // =========================================================================
    // Summary
    // =========================================================================

    std::cout << bold << red
              << "\n   dt_max_stable (p=" << p << ") = "
              << std::setprecision(12) << dt_max_stable
              << reset << std::endl;

    log << "dt_max_stable=" << dt_max_stable << "\n";
    log.flush();

    cpu.toc();
    log << "CPU=" << cpu << "\n";
    std::cout << bold << red << "   CPU: " << cpu << reset << std::endl << std::endl;
}
