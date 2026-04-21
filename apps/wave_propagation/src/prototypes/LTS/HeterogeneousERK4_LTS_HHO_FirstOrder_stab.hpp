

//  Created by Romain Mottier
// ../wave_propagation -k3 -s0 -r0 -c0 -m1 -l0 -n4500 -p1 -f1 -e0
// WITHOUT LOCAL REFINEMENT: ../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n220 -p1 -f1 -e0
// WITH LOCAL REFINEMENT LVL 3:../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n220 -p3 -f1 -e0
// ../../../wave_propagation -k3 -s0 -r0 -c0 -m0 -l5 -n225 -p5 -f1 -e0

void HeterogeneousERK4_LTS_HHO_FirstOrder_stab(int argc, char **argv);

void HeterogeneousERK4_LTS_HHO_FirstOrder_stab(int argc, char **argv){
    
    // ######################################################################
    // ###################################################################### Simulation parameters
    // ######################################################################
    
    std::cout << std::endl << bold << red << "   RK4 - LTS - PULSE - COUPLING" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    cpu.tic();
    
    // ######################################################################
    // ###################################################################### Mesh generation
    // ######################################################################
    
    tc.tic();
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;
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
        
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p5.txt");
        mesh_files.push_back("/home/romain/GitHub/Disk/meshes/nonconform_square_coupling_p5.txt");
        
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
        // std::vector<size_t> cells_to_refine = {2589, 2590, 2591, 2592, 2593, 2594,
        //                                        2525, 2526, 2527, 2528, 2529, 2530,
        //                                        2461, 2462, 2463, 2464, 2465, 2466,
        //                                        2397, 2398, 2399, 2400, 2401, 2402,
        //                                        2333, 2334, 2335, 2336, 2337, 2338};
        // mesh_builder.refine_cells(cells_to_refine, sim_data.m_substeps_Q);
        mesh_builder.move_to_mesh_storage(msh);
    }
    
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;

    RealType h_max = 1e-5;
    RealType h_min = 10;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min)      h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    const int p_stab = static_cast<int>(std::round(h_max / h_min));
    auto h_c = 0.75 * h_max;
    if (p_stab == 1) {
        h_c = 1.1 * h_max;
    }

    std::cout << bold << cyan << "      h_max = "            << h_max       << reset << std::endl;
    std::cout << bold << cyan << "      h_min = "            << h_min       << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = "      << h_max/h_min << std::endl;
    std::cout << bold << cyan << "      p_stab (rounded) = " << p_stab      << reset << std::endl << std::endl;

    // ######################################################################
    // ###################################################################### Time controls
    // ######################################################################

    const size_t   nt = sim_data.m_nt_divs;
    const RealType ti = 0.0;
    const RealType tf = 0.25;
    const RealType dt = (tf - ti) / static_cast<RealType>(nt);
    
    // ######################################################################
    // ###################################################################### HHO setting
    // ######################################################################
    
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);
    
    // ##################################################
    // ################################################## Material data
    // ##################################################
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
        RealType rho = 1.0;
        RealType vp  = std::sqrt(3.0);
        RealType vs  = 1.0;
        return elastic_material_data<RealType>(rho,vp,vs);
    };
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        RealType rho = 1.0;
        RealType vp  = 1.0;
        return acoustic_material_data<RealType>(rho,vp);
    };

    auto water_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> acoustic_material_data<double> {
        return acoustic_material_data<double>(1.0, 1.0);
    };
    
    auto granit_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> elastic_material_data<double> {
        return elastic_material_data<double>(2.624390244, 4.0, 2.0);
    };
    
    // ##################################################
    // ################################################## Structure setting
    // ##################################################
    
    std::map<size_t,elastic_material_data<RealType>>  e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    
    RealType eps         = 1.0e-10;
    RealType y_interface = 0.0;
    
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.y() - y_interface) < eps) {
            interface_face_indexes.insert(fc_id);
        }
    }
    
    for (auto & cell : msh) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        
        if (bar.y() > y_interface) {
            a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
        } else {
            e_material.insert(std::make_pair(cell_ind, elastic_mat_fun(bar)));
        }
        
        auto cell_faces = faces(msh, cell);
        for (auto face : cell_faces) {
            auto fc_id = msh.lookup(face);
            if (interface_face_indexes.find(fc_id) != interface_face_indexes.end()) {
                if (bar.y() > y_interface)
                    interface_cell_pair_indexes[fc_id].second = cell_ind;
                else
                    interface_cell_pair_indexes[fc_id].first  = cell_ind;
            }
        }
    }
    
    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (interface_face_indexes.find(fc_id) == interface_face_indexes.end()) {
            if (bar.y() > y_interface) acoustic_internal_faces.insert(fc_id);
            else                       elastic_internal_faces.insert(fc_id);
        }
    }
    
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (bar.y() > y_interface) {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_acoustic_id, true};
            acoustic_bc_face_indexes.insert(fc_id);
        } else {
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_elastic_id, true};
            elastic_bc_face_indexes.insert(fc_id);
        }
    }
    
    auto null_s_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type&) -> double {
        return 0.0;
    };
    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type&) -> disk::static_vector<double, 2> {
        return disk::static_vector<double, 2>{0, 0};
    };
    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type&) -> disk::static_matrix<double,2,2> {
        return disk::static_matrix<double,2,2>::Zero(2,2);
    };
    
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
    
    // ######################################################################
    // ###################################################################### Assembly
    // ######################################################################
    
    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << red  << "   ASSEMBLY 1 : " << std::endl;
    std::cout << bold << cyan << "      Assembler generation : " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "      Mass Assembly : " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_coupling_terms(msh);
    tc.toc();
    std::cout << bold << cyan << "      Coupling assembly : " << tc << " seconds" << reset << std::endl << std::endl;
    
    // ######################################################################
    // ###################################################################### Initial condition
    // ######################################################################
    
    auto v_fun_adi_acoustic = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
        double x    = pt.x(),  y  = pt.y();
        double xc   = 0.0,     yc = 0.1;
        double fc   = 10.0,    vp = std::sqrt(1.0);
        double lp   = vp/fc;
        double r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        double wave = 10.0 / std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI);
        return disk::static_vector<double, 2>{wave*(x-xc), wave*(y-yc)};
    };
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_fun_adi_acoustic);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
    
    // ######################################################################
    // ###################################################################### ERK scheme setup
    // ######################################################################
    
    Matrix<RealType, Dynamic, Dynamic> a_mat;
    Matrix<RealType, Dynamic, 1>       b_vec, c_vec;
    
    std::cout << bold << red  << "   ASSEMBLY 2 : " << std::endl;
    std::cout << bold << cyan << "      First stiffness assembly completed: ";
    tc.tic();
    assembler.assemble(msh, null_fun, null_s_fun, true);
    tc.toc();
    std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
    assembler.LHS += assembler.COUPLING;
    
    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs        = assembler.get_e_face_dof();
    size_t a_face_dofs        = assembler.get_a_face_dof();
    
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
                                              elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
                       assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
                       assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
                       assembler.get_e_compress(), assembler.get_a_compress(),
                       elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    
    tc.toc();
    std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;
    tc.tic();
    erk_an.refresh_faces_unknowns(x_dof);
    tc.toc();
    std::cout << bold << cyan << "      Inverse of Sff + Coupling in: " << tc << " seconds" << reset << std::endl;
    
    // ######################################################################
    // ###################################################################### Log file
    // ######################################################################
    
    std::ostringstream filename;
    filename << "Explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs
             << "_k_" << sim_data.m_k_degree << "_s_" << 4 << ".txt";
    std::ofstream simulation_log(filename.str());
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << 4      << std::endl;
    simulation_log << "Number of time steps =  " << nt     << std::endl;
    simulation_log << "Step size =  "             << dt     << std::endl;
    simulation_log << "p_stab =  "                << p_stab << std::endl;
    simulation_log << "Number of equations : "    << assembler.RHS.rows() << std::endl;
    simulation_log.flush();

    // ######################################################################
    // ###################################################################### Projection matrix P + silo (mesh visualisation)
    // ######################################################################

    assembler.assemble_P(msh, h_c);

    // Export silo de la condition initiale pour visualiser le maillage
    // et le partitionnement coarse/fine avant le sweep de stabilité.
    if (sim_data.m_render_silo_files_Q) {
        std::ostringstream silo_name;
        silo_name << "silo_stab_l_" << sim_data.m_n_divs
                  << "_k_" << sim_data.m_k_degree
                  << "_p_" << p_stab << "_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic_LTS(
                silo_name.str(), 0, msh, hho_di, x_dof, e_material, a_material, false, h_c);
        std::cout << bold << cyan << "      Silo file written: "
                  << silo_name.str() << "0" << reset << std::endl;
    }
    
    // ######################################################################
    // ###################################################################### Stability analysis
    // ######################################################################

    std::cout << bold << red << "\n   STABILITY ANALYSIS (LTS-RK4, p=" << p_stab << ")" << reset << std::endl;
    simulation_log << "\n============================================\n";
    simulation_log << "STABILITY ANALYSIS (LTS-RK4, p=" << p_stab << ")\n";
    simulation_log << "============================================\n";
    simulation_log << std::setw(20) << "dt"
                   << std::setw(22) << "rho(C_LTS)"
                   << std::setw(12) << "stable" << "\n";

    const double dt_stab_min  = 1e-3;
    const double dt_stab_max  = 2.0 * dt;
    const int    nb_dt_points = 60;

    const int    n_dof  = static_cast<int>(x_dof.rows());
    const int    n_eigs = 10;
    const int    n_cv   = std::min(n_dof, std::max(30, 20 * n_eigs));

    const double ddt = (dt_stab_max - dt_stab_min) / static_cast<double>(nb_dt_points - 1);

    auto apply_C_LTS = [&](const Matrix<RealType, Dynamic, 1> & x_in, double dtau)
        -> Matrix<RealType, Dynamic, 1>
    {
        std::vector<Matrix<RealType, Dynamic, 1>> w(4);
        for (int i = 0; i < 4; ++i) { w[i].resize(x_dof.rows()); w[i].setZero(); }
        erk_an.ZeroFc();
        erk_an.erk_weight_LTS_coarse_old(x_in, assembler.Pcoarse, w);

        Matrix<RealType, Dynamic, 1> x = x_in;
        for (int m = 0; m < p_stab; ++m) {
            const double tm   =  m        * dtau;
            const double tm12 = (m + 0.5) * dtau;
            const double tm1  = (m + 1.0) * dtau;

            auto coarse_at = [&](double t) -> Matrix<RealType, Dynamic, 1> {
                return w[0] + t*w[1] + (t*t/2.0)*w[2] + (t*t*t/6.0)*w[3];
            };

            Matrix<RealType, Dynamic, 1> yn1, yn2, yn3, yn4;
            Matrix<RealType, Dynamic, 1> k1, k2, k3, k4;

            yn1 = assembler.Pfine * x;
            erk_an.erk_weight(yn1, k1);
            k1 += coarse_at(tm);

            yn2 = assembler.Pfine * (x + dtau/2.0*k1);
            erk_an.erk_weight(yn2, k2);
            k2 += coarse_at(tm12);

            yn3 = assembler.Pfine * (x + dtau/2.0*k2);
            erk_an.erk_weight(yn3, k3);
            k3 += coarse_at(tm12);

            yn4 = assembler.Pfine * (x + dtau*k3);
            erk_an.erk_weight(yn4, k4);
            k4 += coarse_at(tm1);

            x += dtau/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
        }
        return x;
    };

    double dt_max_stable = -1.0;

    for (int s = 0; s < nb_dt_points; ++s) {

        const double dt_s = dt_stab_min + s * ddt;
        const double dtau = dt_s / static_cast<double>(p_stab);

        tc.tic();
        Eigen::SparseMatrix<double> C_LTS(n_dof, n_dof);
        {
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(n_dof);
            Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);
            for (int i = 0; i < n_dof; ++i) {
                e_i(i) = 1.0;
                Matrix<RealType, Dynamic, 1> col = apply_C_LTS(e_i, dtau);
                for (int j = 0; j < n_dof; ++j) {
                    if (std::abs(col(j)) > 1e-15)
                        triplets.emplace_back(j, i, col(j));
                }
                e_i(i) = 0.0;
            }
            C_LTS.setFromTriplets(triplets.begin(), triplets.end());
        }
        tc.toc();
        std::cout << bold << cyan << "      C_LTS built in " << tc << " s" << reset << std::endl;

        Spectra::SparseGenMatProd<double> op(C_LTS);
        Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>> eigs(op, n_eigs, n_cv);
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-12);

        const bool ok  = (eigs.info() == Spectra::CompInfo::Successful);
        double     rho = -1.0;
        if (ok) {
            rho = eigs.eigenvalues().cwiseAbs().maxCoeff();
            if (rho <= 1.0) dt_max_stable = dt_s;
        }

        // Affichage haute précision (15 chiffres) pour distinguer
        // rho = 1.0000001 de rho = 0.9999999, invisibles avec setprecision(6).
        std::cout << bold << cyan
                  << "      dt = " << std::setw(14) << std::setprecision(8)  << dt_s
                  << "   rho = "   << std::setw(22) << std::setprecision(15) << rho;
        if (!ok)            std::cout << "  [Spectra failed]";
        if (ok && rho<=1.0) std::cout << "  [stable]";
        else                std::cout << "  [UNSTABLE]";
        std::cout << reset << std::endl;

        simulation_log << std::setw(20) << std::setprecision(10) << dt_s
                       << std::setw(22) << std::setprecision(15) << rho
                       << std::setw(12) << (ok && rho<=1.0 ? "yes" : "no") << "\n";
    }

    std::cout << bold << red
              << "   --> dt_max_stable (p=" << p_stab << ") = "
              << std::setprecision(10) << dt_max_stable
              << reset << std::endl;
    simulation_log << "dt_max_stable (p=" << p_stab << ") = "
                   << std::setprecision(10) << dt_max_stable << "\n";

    simulation_log.flush();

    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << std::endl << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
}
