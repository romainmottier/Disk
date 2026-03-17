

//  Created by Romain Mottier
// ../wave_propagation -k 2 -s 0 -r 0 -c 0 -p 0 -l 6 -n 2500 -f 1 -e 0
// WITHOUT LOCAL REFINEMENT: ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n220 -p1 -f1 -e0
// WITH LOCAL REFINEMENT LVL 3:../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n220 -p3 -f1 -e0
void HeterogeneousERK4_LTS_HHO_FirstOrder(int argc, char **argv);

void HeterogeneousERK4_LTS_HHO_FirstOrder(int argc, char **argv){
    
    // ######################################################################
    // ###################################################################### Simulation paramaters 
    // ######################################################################
    
    std::cout << std::endl << bold << red << "   RK4 - LTS - PULSE - COUPLING" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, cpu;
    cpu.tic();
    
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
        
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p1.txt");    // l = 0
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p2.txt");    // l = 1
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p3.txt");    // l = 2
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p4.txt");    // l = 3
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p5.txt");    // l = 4
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
        mesh_builder.set_translation_data(-0.5, -0.5);
        mesh_builder.build_mesh();
        std::vector<size_t> cells_to_refine = {2589, 2590, 2591, 2592, 2593, 2594,
                                               2525, 2526, 2527, 2528, 2529, 2530,
                                               2461, 2462, 2463, 2464, 2465, 2466,
                                               2397, 2398, 2399, 2400, 2401, 2402,
                                               2333, 2334, 2335, 2336, 2337, 2338};
        mesh_builder.refine_cells(cells_to_refine, 7);
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
    
    auto water_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> acoustic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp;
        rho = 1.0;            
        vp  = 1.0;     
        acoustic_material_data<double> material(rho,vp);
        return material;
    };
    
    auto granit_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> elastic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp, vs;
        rho = 2.624390244; 
        vp  = 4.0;    
        vs  = 2.0; 
        elastic_material_data<double> material(rho,vp,vs);
        return material;
    };

    // ###################################################################### 
    // ###################################################################### Structure setting 
    // ###################################################################### 

    std::map<size_t,elastic_material_data<RealType>>  e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    
    RealType eps = 1.0e-10;
    RealType y_interface = 0.0;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.y() - y_interface) < eps) {
            interface_face_indexes.insert(fc_id);
            continue;
        }
    }

    
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        
        // Assigning the material properties
        if (bar.y() > y_interface) {
            // acoustic_material_data<RealType> material = water_mat_fun_adi(bar); 
            acoustic_material_data<RealType> material = acoustic_mat_fun(bar); 
            a_material.insert(std::make_pair(cell_ind,material));
        }
        else {
            // elastic_material_data<RealType> material = granit_mat_fun_adi(bar); 
            elastic_material_data<RealType> material = elastic_mat_fun(bar); 
            e_material.insert(std::make_pair(cell_ind,material));
        }
        
        // Detection of faces on the interfaces
        auto cell_faces = faces(msh,cell);
        for (auto face :cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) {
                if (bar.y() > y_interface) {
                    interface_cell_pair_indexes[fc_id].second = cell_ind;
                }
                else {
                    interface_cell_pair_indexes[fc_id].first = cell_ind;
                }
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
        }
        else {
            if (bar.y() > y_interface) {
                acoustic_internal_faces.insert(fc_id);
            }
            else {
                elastic_internal_faces.insert(fc_id);
            }
        }
    }
    
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (bar.y() > y_interface) {
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

    // Boundary condition
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);

    // ##################################################
    // ################################################## Solving a primal HHO mixed problem 
    // ##################################################

    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }

    tc.toc();
    std::cout << bold << red << "   ASSEMBLY 1 : " << std::endl;
    std::cout << bold << cyan << "      Assembler generation : ";
    std::cout << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "      Mass Assembly : ";
    std::cout << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_coupling_terms(msh);
    tc.toc();
    std::cout << bold << cyan << "      Coupling assembly : ";
    std::cout << tc << " seconds" << reset << std::endl << std::endl;    
  
    // ######################################################################
    // ###################################################################### Projecting initial data 
    // ######################################################################
    
    auto v_fun_adi_acoustic = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
        double x,y,xc,yc,r,wave,vx,vy,c,lp, fc, vp;
        x    = pt.x();
        y    = pt.y();
        xc   = 0.0;
        yc   = 0.1; // 0.1;
        fc   = 10.0;
        c    = 10;
        vp   = std::sqrt(1.0);
        lp   = vp/fc;
        r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
        vx   = wave*(x-xc);
        vy   = wave*(y-yc);
        disk::static_vector<double, 2> v{vx,vy};
        return v;
    };
    
    Matrix<RealType, Dynamic, 1> x_dof;
    // Acoustic pulse intialized in pressure 
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_fun_adi_acoustic);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
    // Elastic pulse intialized in pressure 
    // assembler.project_over_cells(msh, x_dof, v_fun, null_flux_fun, null_s_fun, null_fun);
    // assembler.project_over_faces(msh, x_dof, v_fun, null_s_fun);
  
    // ##################################################
    // ################################################## Solving a first order equation HDG/HHO propagation problem
    // ##################################################

    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    // ERK(s) schemes
    std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
    std::cout << bold << cyan << "      First stiffness assembly completed: ";
    tc.tic();
    assembler.assemble(msh, null_fun, null_s_fun, true);
    tc.toc();
    std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
    assembler.LHS += assembler.COUPLING; 
    
    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();
    
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);//assembler.get_interfaces());
    
    tc.toc();
    std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;
    tc.tic();
    erk_an.refresh_faces_unknowns(x_dof);
    tc.toc();
    std::cout << bold << cyan << "      Inverse of Sff + Coupling in: " << tc << " seconds" << reset << std::endl;
    
    // ##################################################
    // ################################################## Preprocessor
    // ##################################################
    
    std::ostringstream filename;
    filename << "Explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << ".txt";
    std::string filename_str = filename.str();
    std::ofstream simulation_log(filename_str);
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << 4 << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log.flush();

    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::ostringstream filename;
        filename << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << "_";
        std::string silo_file_name = filename.str();
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false, h_c);
    }

    // ##################################################
    // ################################################## Sensors
    // ##################################################

    bool e_side_Q = true;
    bool a_side_Q = false;

    std::ostringstream filename_acou;
    filename_acou << "A_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << ".csv";
    std::string filename_acou_str = filename_acou.str();
    std::ofstream Acoustic_sensor_1_log(filename_acou_str);
    typename mesh_type::point_type Acoustic_s1_pt(-0.15,  0.1);
    std::pair<typename mesh_type::point_type,size_t> Acoustic_s1_pt_cell  = std::make_pair(Acoustic_s1_pt, -1);

    std::ostringstream filename_int;
    filename_int <<  "I_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << ".csv";
    std::string filename_int_str = filename_int.str();
    std::ofstream Interface_sensor_1_log(filename_int_str);    
    typename mesh_type::point_type Interface_s1_pt(-0.15, 0.0);
    std::pair<typename mesh_type::point_type,size_t> Interface_s1_pt_cell = std::make_pair(Interface_s1_pt, -1);

    std::ostringstream filename_ela;
    filename_ela <<  "E_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << 4 << ".csv";
    std::string filename_ela_str = filename_ela.str();
    std::ofstream Elastic_sensor_1_log(filename_ela_str);
    typename mesh_type::point_type Elastic_s1_pt(-0.15,  -0.1);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s1_pt_cell = std::make_pair(Elastic_s1_pt, -1);

    bool sensors = false;
    if (sensors) {
        postprocessor<mesh_type>::record_acoustic_data_elasto_acoustic_four_fields(0, Acoustic_s1_pt_cell, msh, hho_di, assembler, x_dof, a_side_Q, Acoustic_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Interface_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Interface_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
    }

    std::cout << std::endl;

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
        ////////////////////////////////////////////////////////////////////////// LOOP OVER THE SUBSTEPS: ERK4 ON THE LOCAL DOFS WITH INJECTION OF THE GLOBAL DOFS
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





























