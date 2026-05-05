
//  Created by Romain Mottier
// WITHOUT LOCAL REFINEMENT: ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l4 -n220 -p0 -f1 -e0
// WITH LOCAL REFINEMENT:    ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l4 -n220 -p3 -f1 -e0

void ERK4_LTS_stab(int argc, char **argv);

void ERK4_LTS_stab(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS STABILITY ANALYSIS" << std::endl << std::endl;

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

    if (sim_data.m_polygonal_mesh_Q) {
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
    else {
        RealType lx = 2.0;
        RealType ly = 1.0;
        size_t   nx = 6;
        size_t   ny = 3;

        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-0.5, -0.5);
        mesh_builder.build_mesh();

        std::vector<size_t> cells_to_refine = {7, 10};
        mesh_builder.refine_cells(cells_to_refine, sim_data.m_substeps_Q);
        mesh_builder.move_to_mesh_storage(msh);
    }

    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;

    // Compute mesh size extrema and refinement ratio p = h_max / h_min
    RealType h_max = 1.0e-5;
    RealType h_min = 10.0;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min) h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }

    // p is the local refinement ratio; used to set h_c (the coarse threshold)
    auto p   = h_max / h_min;
    auto h_c = (p == 1) ? 1.25 * h_max : 0.75 * h_max;

    std::cout << bold << cyan << "      h_max       = " << h_max         << reset << std::endl;
    std::cout << bold << cyan << "      h_min       = " << h_min                  << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << p              << reset << std::endl << std::endl;

    // =========================================================================
    // Time parameters
    // =========================================================================

    // nt is taken directly from the command-line argument -n;
    // dt = (tf - ti) / nt is the coarse time step whose stability we study
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

    // Elastic domain (x < 0): rho=1, lambda=sqrt(3), mu=1
    auto elastic_mat_fun = [](const typename mesh_type::point_type&) -> elastic_material_data<RealType> {
        return elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0);
    };
    // Acoustic domain (x > 0): rho=1, c=1
    auto acoustic_mat_fun = [](const typename mesh_type::point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };

    // =========================================================================
    // Analytical functions (used only for projection of initial data)
    // =========================================================================

    scal_vec_analytic_functions functions;
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);

    auto u_fun      = functions.Evaluate_u(ti);
    auto v_fun      = functions.Evaluate_v(ti);
    auto a_fun      = functions.Evaluate_a(ti);
    auto f_fun      = functions.Evaluate_f(ti);
    auto flux_fun   = functions.Evaluate_sigma(ti);

    auto s_u_fun    = functions.Evaluate_s_u(ti);
    auto s_v_fun    = functions.Evaluate_s_v(ti);
    auto s_a_fun    = functions.Evaluate_s_a(ti);
    auto s_f_fun    = functions.Evaluate_s_f(ti);
    auto s_flux_fun = functions.Evaluate_s_q(ti);

    // =========================================================================
    // Domain decomposition: elastic (x<0) / acoustic (x>0) / interface (x=0)
    // =========================================================================

    std::map<size_t, elastic_material_data<RealType>>  e_material;
    std::map<size_t, acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t, std::pair<size_t, size_t>> interface_cell_pair_indexes;

    const RealType eps = 1.0e-10;

    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.x()) < eps) {
            interface_face_indexes.insert(fc_id);
            continue;
        }
    }

    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        if (bar.x() > 0) {
            a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
        }
        else {
            e_material.insert(std::make_pair(cell_ind, elastic_mat_fun(bar)));
        }
        // For each face of this cell, check if it lies on the interface
        auto cell_faces = faces(msh, cell);
        for (auto face : cell_faces) {
            auto fc_id     = msh.lookup(face);
            bool is_iface  = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_iface) {
                if (bar.x() > 0)
                    interface_cell_pair_indexes[fc_id].second = cell_ind;
                else
                    interface_cell_pair_indexes[fc_id].first  = cell_ind;
            }
        }
    }

    // Classify internal interface faces as elastic-side or acoustic-side
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id   = msh.lookup(face);
        bool is_iface = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
        if (is_iface) {
            if (bar.y() > 0)
                acoustic_internal_faces.insert(fc_id);
            else
                elastic_internal_faces.insert(fc_id);
        }
    }

    // Assign Dirichlet boundary conditions
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
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

    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  u_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);

    // =========================================================================
    // HHO assembly
    // =========================================================================

    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);

    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);

    // Project initial data onto cell and face unknowns
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);

    // Build the global stiffness operator (LHS = stiffness + coupling)
    assembler.assemble(msh, f_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING;

    // Build the ERK scheme object and invert mass / Schur complement on cell/face blocks
    erk_coupling_hho_scheme<RealType> erk_an(
        assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
        assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(),
        assembler.get_e_face_dof(),    assembler.get_a_face_dof());

    erk_an.Mcc_inverse(
        assembler.get_elastic_cells(),  assembler.get_acoustic_cells(),
        assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());

    erk_an.Sff_inverse(
        assembler.get_elastic_faces(),  assembler.get_acoustic_faces(),
        assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
        assembler.get_e_compress(),     assembler.get_a_compress(),
        elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);

    erk_an.refresh_faces_unknowns(x_dof);

    // Build the coarse/fine projection matrices P (threshold h_c)
    assembler.assemble_P(msh, h_c);

    // Optional: dump initial state to Silo for visualisation
    if (sim_data.m_render_silo_files_Q) {
        std::ostringstream sn;
        sn << "silo_stab_l_" << sim_data.m_n_divs
           << "_k_" << sim_data.m_k_degree
           << "_p_" << p << "_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(
            sn.str(), 0, msh, hho_di, x_dof, e_material, a_material, false, h_c);
    }

    // =========================================================================
    // Log file
    // =========================================================================

    std::ostringstream fname;
    fname << "stab_l_" << sim_data.m_n_divs
          << "_n_"     << nt
          << "_k_"     << sim_data.m_k_degree
          << "_p_"     << p << ".txt";
    std::ofstream log(fname.str());
    sim_data.write_simulation_data(log);
    log << "nt=" << nt << " dt=" << dt << " p=" << p
        << " n_dof=" << x_dof.rows() << "\n";
    log.flush();

    // =========================================================================
    // Stability sweep
    // =========================================================================
    //
    // We build the amplification matrix C of size (n_c x n_c), where n_c is
    // the number of cell DOFs. Face unknowns are NOT independent — they are
    // eliminated via static condensation inside erk_weight, so they do not
    // appear as columns of C.
    //
    // Column i of C is obtained by:
    //   1. Setting e_i = canonical basis vector of size n_dof (zero everywhere
    //      except at cell position i); face entries are left at zero.
    //   2. Applying one full LTS-RK4 step (coarse predictor + p fine substeps)
    //      with source term f = 0 (valid for stability analysis).
    //   3. Extracting the first n_c components (cell part) of the result.
    //
    // The scheme is stable for a given dt iff the spectral radius rho(C) <= 1.
    // We sweep dt from dt_min to dt_max and stop as soon as rho > 1.25.

    const int n_dof = static_cast<int>(x_dof.rows());
    const int n_c   = static_cast<int>(erk_an.n_c_dof());
    const int n_f   = n_dof - n_c;

    std::cout << bold << cyan
              << "      n_dof=" << n_dof << "  n_c=" << n_c << "  n_f=" << n_f
              << "  dt=" << dt << reset << std::endl;

    const double dt_min = dt;
    const double dt_max = 10.0 * dt;
    const int    n_pts  = 100;
    const double ddt    = (dt_max - dt_min) / static_cast<double>(n_pts - 1);

    std::cout << bold << red << "\n   STABILITY SWEEP (p=" << p << ")" << reset << std::endl;
    log << std::setw(20) << "dt" << std::setw(22) << "rho" << std::setw(10) << "stable\n";

    // Zero the face correction accumulator (no source term for stability analysis)
    erk_an.ZeroFc();
    double dt_max_stable = -1.0;

    for (int s = 0; s < n_pts; ++s) {
        const double dt_s   = dt_min + s * ddt;
        const double dtau_s = dt_s / static_cast<double>(p);

        // -----------------------------------------------------------------
        // Build the amplification matrix C column by column
        // -----------------------------------------------------------------
        tc.tic();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_c, n_c);

        for (int i = 0; i < n_c; ++i) {
            // Canonical basis vector: cell i = 1, all faces = 0
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            e_i(i) = 1.0;

            // --- Coarse predictor (no source term) ---
            // Computes the Taylor coefficients w[0..3] used to interpolate
            // the coarse-grid right-hand side across the fine substeps.
            // The "_old" variant is used because there is no source term (f=0),
            // which makes it faster than the general version.
            std::vector<Matrix<RealType, Dynamic, 1>> w(4);
            for (int j = 0; j < 4; ++j) { w[j].resize(n_dof); w[j].setZero(); }
            erk_an.erk_weight_LTS_coarse_old(e_i, assembler.Pcoarse, w);

            // --- Fine substeps: p RK4 steps of size dtau_s ---
            // At each substep m, the coarse contribution is reconstructed
            // from the Taylor polynomial w(tau) = w0 + tau*w1 + tau^2/2*w2 + tau^3/6*w3
            Matrix<RealType, Dynamic, 1> x = e_i;
            for (int m = 0; m < p; ++m) {
                const double tm   =  m        * dtau_s;
                const double tm12 = (m + 0.5) * dtau_s;
                const double tm1  = (m + 1.0) * dtau_s;

                // Taylor interpolation of the coarse predictor at time tau
                auto Taylor = [&](double tau) {
                    return w[0] + tau * w[1]
                                + (tau * tau / 2.0)       * w[2]
                                + (tau * tau * tau / 6.0) * w[3];
                };

                Matrix<RealType, Dynamic, 1> k1, k2, k3, k4;
                Matrix<RealType, Dynamic, 1> yn;

                // RK4 stage 1
                yn = assembler.Pfine * x;
                erk_an.erk_weight(yn, k1);
                k1 += Taylor(tm);

                // RK4 stage 2
                yn = assembler.Pfine * (x + dtau_s / 2.0 * k1);
                erk_an.erk_weight(yn, k2);
                k2 += Taylor(tm12);

                // RK4 stage 3
                yn = assembler.Pfine * (x + dtau_s / 2.0 * k2);
                erk_an.erk_weight(yn, k3);
                k3 += Taylor(tm12);

                // RK4 stage 4
                yn = assembler.Pfine * (x + dtau_s * k3);
                erk_an.erk_weight(yn, k4);
                k4 += Taylor(tm1);

                // Update solution
                x += dtau_s / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            }

            // Column i = cell part of the output vector
            C.col(i) = x.head(n_c);
        }
        tc.toc();

        std::cout << bold << cyan
                  << "      C(" << n_c << "x" << n_c << ") built in " << tc
                  << "   norm=" << std::setprecision(6) << C.norm()
                  << reset << std::endl;

        // -----------------------------------------------------------------
        // Compute spectral radius rho(C) = max |lambda_i(C)|
        // The scheme is stable iff rho(C) <= 1.
        // -----------------------------------------------------------------
        Eigen::EigenSolver<Eigen::MatrixXd> es(C);
        double rho = -1.0;

        if (es.info() == Eigen::Success) {
            rho = es.eigenvalues().cwiseAbs().maxCoeff();
            if (rho <= 1.0) {
                dt_max_stable = dt_s;
                std::cout << "   --> dt_max_stable updated = " << dt_max_stable << std::endl;
            }
        }
        else {
            std::cout << bold << red
                      << "   --> EigenSolver FAILED at dt=" << dt_s
                      << reset << std::endl;
        }

        std::cout << bold << cyan
                  << "      dt=" << std::setw(14) << std::setprecision(8)  << dt_s
                  << "   rho=" << std::setw(22) << std::setprecision(15) << rho
                  << (rho >= 0 && rho <= 1.0 ? "  [stable]" : "  [UNSTABLE]")
                  << reset << std::endl;

        log << std::setw(20) << std::setprecision(10) << dt_s
            << std::setw(22) << std::setprecision(15) << rho
            << std::setw(10) << (rho >= 0 && rho <= 1.0 ? "yes" : "no") << "\n";
        log.flush();

        // Save all eigenvalues (real + imag parts) to a separate file for
        // post-processing and visualisation in the complex plane
        if (es.info() == Eigen::Success) {
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

        // Stop the sweep as soon as the scheme becomes significantly unstable
        if (rho > 1.25) {
            std::cout << bold << red
                      << "   --> rho > 1.25, stopping sweep."
                      << reset << std::endl;
            break;
        }

    } // end stability sweep

    // =========================================================================
    // Summary
    // =========================================================================

    std::cout << bold << red
              << "\n   dt_max_stable (p=" << p << ") = "
              << std::setprecision(10) << dt_max_stable
              << reset << std::endl;

    log << "dt_max_stable=" << dt_max_stable << "\n";
    log.flush();

    cpu.toc();
    log << "CPU=" << cpu << "\n";
    std::cout << bold << red << "   CPU: " << cpu << reset << std::endl << std::endl;
}
