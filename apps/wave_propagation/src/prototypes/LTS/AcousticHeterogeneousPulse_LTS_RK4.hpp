
//  Contributions by Omar Durán and Romain Mottier

// ../wave_propagation -k3 -s0 -r0 -c0 -m1 -l0 -n750 -p1 -f1 -e0
// ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n650 -p1 -f1 -e0
// ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l6 -n1350 -p1 -f1 -e0

void AcousticHeterogeneousPulse_LTS_RK4(int argc, char **argv);
void AcousticHeterogeneousPulse_LTS_RK4(int argc, char **argv) {
    
    // ######################################################################
    // ###################################################################### Simulation paramaters 
    // ######################################################################

    std::cout << std::endl << bold << red << "   LTS RK4 PULSE - Acoustic" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, simulation_tc;
    simulation_tc.tic();;
//     
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
        
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p1.txt");    // l = 0
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p2.txt");    // l = 1
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p3.txt");    // l = 2
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p4.txt");    // l = 3
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_p5.txt");    // l = 4

        // mesh_files.push_back("/home/romain/GitHub/MESHES_DISK/nonconform_square.txt");       // l = 0

        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
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
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        RealType h_l = diameter(msh, cell);
        if (h_l < h) {
            h = h_l;
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
            RealType x,y;
            x = pt.x();
            y = pt.y();
            return 0.0;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        return {0,0};
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
    
    a_boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(2);
        RealType rho, vp;
        rho = 1.0;
        if (y < 0.5) {
            vp = 5.0;
        }
        else{
            vp = 1.0;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
        return mat_data;
    };

    // ##################################################
    // ################################################## Solving a primal HHO mixed problem 
    // ##################################################

    std::cout << bold << red << "   ASSEMBLY: " << reset << std::endl;
    auto assembler = acoustic_two_fields_assembler_LTS<mesh_type>(msh, hho_di, bnd);
    assembler.load_material_data(msh, acoustic_mat_fun);
    
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
    assembler.assemble(msh, null_fun, true);
    tc.toc();
    std::cout << bold << cyan << "      Stiffness and rhs assembly completed: " << tc << " seconds" << reset << std::endl;
    size_t n_face_dof = assembler.get_n_face_dof();
    tc.tic();
    erk_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS,n_face_dof);
    erk_an.Mcc_inverse(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
    erk_an.Sff_inverse(std::make_pair(assembler.get_n_faces(), assembler.get_face_basis_data()));
    tc.toc();
    std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;
    
    // ######################################################################
    // ###################################################################### Projecting initial data 
    // ######################################################################
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, vel_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vel_fun);
    erk_an.refresh_faces_unknowns(x_dof);
    if (sim_data.m_render_silo_files_Q) {
        std::string silo_file_name = "e_inhomogeneous_scalar_mixed_";
        postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, 0, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
    }

    // ##################################################
    // ################################################## Time marching: LTS - EULER - HHO
    // ##################################################
    
    // ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n4500 -p1 -f1 -e0 CPU TIME: 19.4503
    // ../wave_propagation -k3 -s0 -r0 -c0 -m1 -l0 -n4500 -p2 -f1 -e0 CPU TIME: 2303.61 ~ 38.39min
    assembler.assemble_P(msh, 0.015625);
    assembler.assemble_P_bis(msh, 0.015625);
    size_t nb_silo_files = 25;
    size_t step_interval = std::max(size_t(1), nt / nb_silo_files);
    std::cout << std::endl;
    std::cout << bold << red << "   TIME MARCHING SCHEME: " << reset << std::endl;
    auto p = sim_data.m_substeps_Q;
    auto dtau = dt / p;
    for(size_t it = 1; it <= nt; it++) {
        //////////////////////////////////////////////////////////////////////////
        tcit.tic();
        RealType tn = dt*(it-1)+ti;
        if (it % step_interval == 0 || it == nt) {
            std::cout << bold << cyan << "      Time step number " << it << ": t = " << t << reset << std::endl;
        }
        //////////////////////////////////////////////////////////////////////////
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
        t += dt;
        if (sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::string silo_file_name = "ricker_LTS_euler_";
            postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
        }
        tcit.toc();
        if (sim_data.m_render_silo_files_Q && (it % step_interval == 0 || it == nt)) {
            std::cout << bold << yellow << "         Iteration completed in " << tcit << " seconds" << reset << std::endl;
        }
    }
    simulation_tc.toc();
    std::cout << std::endl << bold << red << "   CPU TIME: " << simulation_tc << std::endl << std::endl;

}

