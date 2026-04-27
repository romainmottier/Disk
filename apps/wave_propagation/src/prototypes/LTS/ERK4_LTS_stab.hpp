
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

    // ##################################################
    // ################################################## Mesh
    // ##################################################

    tc.tic();
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;
    mesh_type msh;

    {
        RealType lx = 1, ly = 1;
        size_t nx = 2, ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-0.5, -0.5);
        mesh_builder.build_mesh();
        std::vector<size_t> cells_to_refine = {2589, 2590, 2591, 2592, 2593, 2594,
                                               2525, 2526, 2527, 2528, 2529, 2530,
                                               2461, 2462, 2463, 2464, 2465, 2466,
                                               2397, 2398, 2399, 2400, 2401, 2402,
                                               2333, 2334, 2335, 2336, 2337, 2338};
        mesh_builder.refine_cells(cells_to_refine, sim_data.m_substeps_Q);
        mesh_builder.move_to_mesh_storage(msh);
    }
    tc.toc();
    std::cout << bold << red << "   MESH GENERATION : " << tc << " seconds" << reset << std::endl;

    RealType h_max = 1e-5, h_min = 10;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min)      h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    const int p = static_cast<int>(std::round(h_max / h_min));
    RealType h_c = (p == 1) ? 1.1 * h_max : 0.75 * h_max;

    std::cout << bold << cyan << "      h_max = " << h_max << "   h_min = " << h_min
              << "   p = " << p << reset << std::endl << std::endl;

    // ##################################################
    // ################################################## Time
    // ##################################################

    const size_t   nt = sim_data.m_nt_divs;
    const RealType ti = 0.0, tf = 0.25;
    const RealType dt = (tf - ti) / nt;

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

    // ##################################################
    // ################################################## Structure
    // ##################################################

    std::map<size_t, elastic_material_data<RealType>>  e_material;
    std::map<size_t, acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t, std::pair<size_t,size_t>> interface_cell_pair_indexes;
    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;

    const RealType eps = 1e-10, y_interface = 0.0;

    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        auto fc_id = msh.lookup(*face_it);
        if (std::fabs(barycenter(msh, *face_it).y() - y_interface) < eps)
            interface_face_indexes.insert(fc_id);
    }
    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        auto bar = barycenter(msh, cell);
        if (bar.y() > y_interface) a_material.insert({cell_ind, acoustic_mat_fun(bar)});
        else                       e_material.insert({cell_ind, elastic_mat_fun(bar)});
        for (auto face : faces(msh, cell)) {
            auto fc_id = msh.lookup(face);
            if (interface_face_indexes.count(fc_id)) {
                if (bar.y() > y_interface) interface_cell_pair_indexes[fc_id].second = cell_ind;
                else                       interface_cell_pair_indexes[fc_id].first  = cell_ind;
            }
        }
    }
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        auto fc_id = msh.lookup(*face_it);
        if (!interface_face_indexes.count(fc_id)) {
            if (barycenter(msh, *face_it).y() > y_interface) acoustic_internal_faces.insert(fc_id);
            else                                              elastic_internal_faces.insert(fc_id);
        }
    }

    size_t bc_elastic_id = 0, bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto fc_id = msh.lookup(*face_it);
        auto bar   = barycenter(msh, *face_it);
        if (bar.y() > y_interface) {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_acoustic_id, true};
            acoustic_bc_face_indexes.insert(fc_id);
        } else {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_elastic_id, true};
            elastic_bc_face_indexes.insert(fc_id);
        }
    }

    auto null_s_fun = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> double {
        return 0.0;
    };
    auto null_fun = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> disk::static_vector<double,2> {
        return disk::static_vector<double,2>{0,0};
    };
    auto null_flux_fun = [](const disk::mesh<double,2,disk::generic_mesh_storage<double,2>>::point_type&) -> disk::static_matrix<double,2,2> {
        return disk::static_matrix<double,2,2>::Zero();
    };

    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);

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

    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_fun_adi_acoustic);
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
    assembler.assemble_P(msh, h_c);

    if (sim_data.m_render_silo_files_Q) {
        std::ostringstream sn;
        sn << "silo_stab_l_" << sim_data.m_n_divs << "_k_" << sim_data.m_k_degree << "_p_" << p << "_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(
            sn.str(), 0, msh, hho_di, x_dof, e_material, a_material, false, h_c);
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
