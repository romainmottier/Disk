
//  Contributions by Omar Durán and Romain Mottier

// Heterogeneous acoustic pulse (vp = 5 below y = 0.5, vp = 1 above): first-order HHO discretization of
// elastoacoustic_four_fields_assembler restricted to acoustic cells, with three explicit time integrators:
//   AcousticHeterogeneousPulse         : ERK(4)
//   AcousticLTSEulerHeterogeneousPulse : two-level LTS with explicit Euler substeps
//   AcousticHeterogeneousPulse_LTS_RK4 : two-level LTS-RK4
//
// ../wave_propagation -k3 -s0 -r0 -c0 -m1 -l0 -n750 -p1 -f1 -e0
// ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n650 -p1 -f1 -e0
// ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l6 -n1350 -p1 -f1 -e0
// ../wave_propagation -k 3 -s 0 -r 0 -c 0 -p 0 -l 6 -n 9500 -f 1 -e 0    (LTS-Euler)

enum class acoustic_pulse_scheme { erk, lts_euler, lts_rk4 };

void acoustic_heterogeneous_pulse(int argc, char **argv, acoustic_pulse_scheme scheme) {

    // ######################################################################
    // ###################################################################### Simulation paramaters
    // ######################################################################

    if (scheme == acoustic_pulse_scheme::erk)
        std::cout << std::endl << bold << red << "   LTS EXPLICIT PULSE - Acoustic" << std::endl << std::endl;
    else if (scheme == acoustic_pulse_scheme::lts_euler)
        std::cout << std::endl << bold << red << "   LTS EULER PULSE - Acoustic" << std::endl << std::endl;
    else
        std::cout << std::endl << bold << red << "   LTS RK4 PULSE - Acoustic" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, simulation_tc;
    simulation_tc.tic();

    // ##################################################
    // ################################################## Mesh generation
    // ##################################################

    tc.tic();
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;

    if (sim_data.m_polygonal_mesh_Q) {
        auto validate_l = [](size_t l) -> size_t {
            if ((0 <= l) && (l < 15) ) {
                return l;
            }
            else {
                std::cout << std::endl << std::endl;
                std::cout << "Warning:: Only few polygonal meshes available.";
                std::cout << std::endl << std::endl;
                return 4;
            }
        };

        size_t l = validate_l(sim_data.m_n_divs);
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;

        if (scheme == acoustic_pulse_scheme::lts_rk4) {
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p1.txt");    // l = 0
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p2.txt");    // l = 1
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p3.txt");    // l = 2
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p4.txt");    // l = 3
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p5.txt");    // l = 4
        }
        else {
            mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square.txt");    // l = 0
            // mesh_files.push_back("../../meshes/pulse/simplices/simplex_l3_0.21.txt");   // l = 1
            // mesh_files.push_back("../../meshes/pulse/simplices/simplex_l4_0.096.txt");  // l = 2
            // mesh_files.push_back("../../meshes/pulse/simplices/simplex_l5_0.0485.txt"); // l = 3
            // mesh_files.push_back("../../meshes/pulse/simplices/simplex_l6_0.024.txt");  // l = 4

            // mesh_files.push_back("../../meshes/pulse/poly/poly_l2.txt");   // -l 0
            // mesh_files.push_back("../../meshes/pulse/poly/poly_l3.txt");   // -l 1
            // mesh_files.push_back("../../meshes/pulse/poly/poly_l4.txt");   // -l 2
            // mesh_files.push_back("../../meshes/pulse/poly/poly_l5.txt");   // -l 3
            // mesh_files.push_back("../../meshes/pulse/poly/poly_l6.txt");   // -l 4
            // mesh_files.push_back("../../meshes/pulse/poly/poly_l7.txt");   // -l 5
        }

        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        if (scheme == acoustic_pulse_scheme::erk) {
            mesh_builder.build_bassin();
        }
        else {
            mesh_builder.build_mesh();
        }
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    else {
        RealType lx = 1;
        RealType ly = 1;
        size_t nx = 2;
        size_t ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-0.0, -0.0);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
    }

    RealType h = 10;
    for (auto & cell : msh ) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h) {
            h = h_l;
        }
        if (scheme == acoustic_pulse_scheme::erk) {
            auto fcs = faces(msh, cell);
            for (auto& fc : fcs) {
                RealType h_f = measure(msh, fc);  // Longueur de la face
                std::cout << "Face id: " << msh.lookup(fc) << " - length: " << h_f << std::endl;
            }
        }
    }

    tc.toc();
    std::cout << std::endl << std::endl;
    std::cout << bold << red << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl << std::endl;

    // ######################################################################
    // ###################################################################### Time controls
    // ######################################################################

    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt = sim_data.m_nt_divs;
    }

    RealType ti = 0.0;
    RealType tf = 0.25;
    RealType dt = (tf-ti)/nt;
    RealType t  = ti;

    // ######################################################################
    // ###################################################################### HHO setting
    // ######################################################################

    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    // ##################################################
    // ################################################## Material data
    // ##################################################

    auto null_fun = [](const mesh_type::point_type& pt) -> RealType {
        return 0.0;
    };

    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        return {0,0};
    };

    // elastic data, unused: the whole domain is acoustic
    auto e_null_fun = [](const typename mesh_type::point_type& pt) -> disk::static_vector<RealType, 2> {
        return disk::static_vector<RealType, 2>::Zero();
    };

    auto e_null_flux_fun = [](const typename mesh_type::point_type& pt) -> disk::static_matrix<RealType, 2, 2> {
        return disk::static_matrix<RealType, 2, 2>::Zero();
    };

    auto a_null_flux_fun = [](const typename mesh_type::point_type& pt) -> disk::static_vector<RealType, 2> {
        return disk::static_vector<RealType, 2>::Zero();
    };

    auto vel_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y,xc,yc,r,wave;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 0.65;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = 0.1*(-4*std::sqrt(10.0/3.0)*(-1 + 1600.0*r*r))/(std::exp(800*r*r)*std::pow(M_PI,0.25));
            return wave;
    };

    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    a_bnd.addDirichletEverywhere(null_fun);
    tc.tic();

    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        RealType rho, vp;
        rho = 1.0;
        if (y < 0.5) {
            vp = 5.0;
        }
        else{
            vp = 1.0;
        }
        acoustic_material_data<RealType> material(rho, vp); // rho, seismic compressional velocity vp
        return material;
    };

    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(barycenter(msh, cell))));
    }

    // ##################################################
    // ################################################## Solving a primal HHO mixed problem
    // ##################################################

    std::cout << bold << red << "   ASSEMBLY: " << reset << std::endl;
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);

    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "      Assembler generation: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "      Mass Assembly completed: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble(msh, e_null_fun, null_fun, true);
    tc.toc();
    std::cout << bold << cyan << "      Stiffness and rhs assembly completed: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    std::set<size_t> elastic_internal_faces, acoustic_internal_faces, interface_face_indexes;
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
        assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(), assembler.get_e_face_dof(), assembler.get_a_face_dof());
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
        assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    // Explicit schemes use the plain least-squares stabilization: the face-face block is block-diagonal
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
        assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
        assembler.get_e_compress(), assembler.get_a_compress(),
        elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    tc.toc();
    std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;

    // ######################################################################
    // ###################################################################### Projecting initial data
    // ######################################################################

    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, e_null_fun, e_null_flux_fun, vel_fun, a_null_flux_fun);
    assembler.project_over_faces(msh, x_dof, e_null_fun, vel_fun);
    erk_an.refresh_faces_unknowns(x_dof);
    if (sim_data.m_render_silo_files_Q) {
        std::string silo_file_name = "e_inhomogeneous_scalar_mixed_";
        postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, 0, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
    }

    if (scheme == acoustic_pulse_scheme::erk) {

        // ##################################################
        // ################################################## Time marching: ERK - HHO
        // ##################################################

        std::ofstream simulation_log("inhomogeneous_acoustic_two_fields_explicit.txt");

        std::ofstream sensor_top_log("top_sensor_e_acoustic_two_fields.csv");
        std::ofstream sensor_bot_log("bot_sensor_e_acoustic_two_fields.csv");
        typename mesh_type::point_type top_pt(0.5, 2.0/3.0);
        typename mesh_type::point_type bot_pt(0.5, 1.0/3.0);
        std::pair<typename mesh_type::point_type,size_t> top_pt_cell = std::make_pair(top_pt, -1);
        std::pair<typename mesh_type::point_type,size_t> bot_pt_cell = std::make_pair(bot_pt, -1);

        postprocessor<mesh_type>::record_data_acoustic_two_fields(0, top_pt_cell, msh, hho_di, x_dof, sensor_top_log);
        postprocessor<mesh_type>::record_data_acoustic_two_fields(0, bot_pt_cell, msh, hho_di, x_dof, sensor_bot_log);

        // Solving a first order equation HDG/HHO propagation problem
        int s = 4;
        Matrix<RealType, Dynamic, Dynamic> a;
        Matrix<RealType, Dynamic, 1> b;
        Matrix<RealType, Dynamic, 1> c;
        erk_butcher_tableau::erk_tables(s, a, b, c);

        Matrix<RealType, Dynamic, 1> x_dof_n;
        timecounter erk_tc;
        erk_tc.tic();

        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            RealType tn = dt*(it-1)+ti;
            // ERK step
            tc.tic();
            {
                size_t n_dof = x_dof.rows();
                Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
                Matrix<RealType, Dynamic, 1> yn, ki;

                x_dof_n = x_dof;
                for (int i = 0; i < s; i++) {
                    yn = x_dof;
                    for (int j = 0; j < s - 1; j++) {
                        yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                    }
                    erk_an.erk_weight(yn, ki);
                    x_dof_n += dt*b(i,0)*ki;
                    k.block(0, i, n_dof, 1) = ki;
                }
            }
            tc.toc();
            std::cout << bold << cyan << "ERK step completed: " << tc << " seconds" << reset << std::endl;
            x_dof = x_dof_n;

            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "e_inhomogeneous_scalar_mixed_";
                postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
            }

            postprocessor<mesh_type>::record_data_acoustic_two_fields(it, top_pt_cell, msh, hho_di, x_dof, sensor_top_log);
            postprocessor<mesh_type>::record_data_acoustic_two_fields(it, bot_pt_cell, msh, hho_di, x_dof, sensor_bot_log);

        }
        erk_tc.toc();
        simulation_log << "Simulation time : " << erk_tc << " seconds" << std::endl;
        simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
        simulation_log << "Number of ERK steps =  " << s << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
        return;
    }

    // ##################################################
    // ################################################## Time marching: LTS - HHO
    // ##################################################

    // LTS-Euler:
    // ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n4500 -p1 -f1 -e0 CPU TIME: 19.4503
    // ../wave_propagation -k3 -s0 -r0 -c0 -m1 -l0 -n4500 -p2 -f1 -e0 CPU TIME: 2303.61 ~ 38.39min
    // Fine cells: diameter below 0.015625, without protection layer
    assembler.assemble_P(msh, 0.015625, 0);
    size_t nb_silo_files = 25;
    size_t step_interval = std::max(size_t(1), nt / nb_silo_files);
    std::cout << std::endl;
    std::cout << bold << red << "   TIME MARCHING SCHEME: " << reset << std::endl;
    auto p = sim_data.m_substeps_Q;
    auto dtau = dt / p;
    for(size_t it = 1; it <= nt; it++) {
        //////////////////////////////////////////////////////////////////////////
        tcit.tic();
        if (it % step_interval == 0 || it == nt) {
            std::cout << bold << cyan << "      Time step number " << it << ": t = " << t << reset << std::endl;
        }
        //////////////////////////////////////////////////////////////////////////
        if (scheme == acoustic_pulse_scheme::lts_euler) {
            Matrix<RealType, Dynamic, 1> k;
            auto yn = x_dof;
            erk_an.erk_euler_LTS(yn, k, dtau, p, assembler.IminusP_cell, assembler.Pfacecoarse, assembler.P_cell, assembler.Pfacefine);
            x_dof = yn;
        }
        else {
            std::vector<Matrix<RealType, Dynamic, 1>> w(4);
            size_t n_dof = x_dof.rows();
            for (int i = 0; i < 4; ++i) {
                w[i].resize(n_dof);
                w[i].setZero();
            }
            Matrix<RealType, Dynamic, 1> yn1(n_dof), yn2(n_dof), yn3(n_dof), yn4(n_dof);
            Matrix<RealType, Dynamic, 1> k1(n_dof),  k2(n_dof),  k3(n_dof),  k4(n_dof);
            auto x_dof_n = x_dof;
            erk_an.compute_wi(x_dof_n, assembler.IminusP_cell, assembler.Pfacecoarse, w);
            //////////////////////////////////////////////////////////////////////////
            for (int m = 0; m < p; m++) {
                // k1
                yn1 = assembler.Pfine * x_dof_n;
                erk_an.erk_weight(yn1, k1);
                k1 += w[0] + m*dtau*w[1] + m*m*dtau*dtau*w[2]/2.0 + m*m*m*dtau*dtau*dtau*w[3]/6.0;
                // k2
                yn2 = assembler.Pfine * (x_dof_n+dtau*k1/2.0);
                erk_an.erk_weight(yn2, k2);
                k2 += w[0] + (m+0.5)*dtau*w[1] + (m+0.5)*(m+0.5)*dtau*dtau*w[2]/2.0 + (m+0.5)*(m+0.5)*(m+0.5)*dtau*dtau*dtau*w[3]/6.0;
                // k3
                yn3 = assembler.Pfine * (x_dof_n+dtau*k2/2.0);
                erk_an.erk_weight(yn3, k3);
                k3 += w[0] + (m+0.5)*dtau*w[1] + (m+0.5)*(m+0.5)*dtau*dtau*w[2]/2.0 + (m+0.5)*(m+0.5)*(m+0.5)*dtau*dtau*dtau*w[3]/6.0;
                // k4
                yn4 = assembler.Pfine * (x_dof_n+dtau*k3);
                erk_an.erk_weight(yn4, k4);
                k4 += w[0] + (m+1.0)*dtau*w[1] + (m+1.0)*(m+1.0)*dtau*dtau*w[2]/2.0 + (m+1.0)*(m+1.0)*(m+1.0)*dtau*dtau*dtau*w[3]/6.0;
                // FINAL UPDATE
                x_dof_n += dtau*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
            }
            x_dof = x_dof_n;
        }
        t += dt;
        if (sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::string silo_file_name = "ricker_LTS_euler_";
            postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
        }
        tcit.toc();
        if (scheme == acoustic_pulse_scheme::lts_rk4 && sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::cout << bold << yellow << "         Iteration completed in " << tcit << " seconds" << reset << std::endl;
        }
    }
    simulation_tc.toc();
    std::cout << std::endl << bold << red << "   CPU TIME: " << simulation_tc << std::endl << std::endl;

}

void AcousticHeterogeneousPulse(int argc, char **argv);
void AcousticHeterogeneousPulse(int argc, char **argv) {
    acoustic_heterogeneous_pulse(argc, argv, acoustic_pulse_scheme::erk);
}

void AcousticLTSEulerHeterogeneousPulse(int argc, char **argv);
void AcousticLTSEulerHeterogeneousPulse(int argc, char **argv) {
    acoustic_heterogeneous_pulse(argc, argv, acoustic_pulse_scheme::lts_euler);
}

void AcousticHeterogeneousPulse_LTS_RK4(int argc, char **argv);
void AcousticHeterogeneousPulse_LTS_RK4(int argc, char **argv) {
    acoustic_heterogeneous_pulse(argc, argv, acoustic_pulse_scheme::lts_rk4);
}
