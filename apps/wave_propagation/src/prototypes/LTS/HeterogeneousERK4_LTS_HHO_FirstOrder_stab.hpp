

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
        
        mesh_files.push_back("/home/mottie0000/Github/Diskpp/meshes/nonconform_square_coupling_p5.txt");    // l = 4
        mesh_files.push_back("/home/romain/GitHub/Disk/meshes/nonconform_square_coupling_p5.txt");          // l = 4

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
        mesh_builder.refine_cells(cells_to_refine, sim_data.m_substeps_Q);
        mesh_builder.move_to_mesh_storage(msh);
    }
    
    tc.toc();
    std::cout << bold << red << std::endl << std::endl << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl;
    RealType h_max = 1e-5;
    RealType h_min = 10;
    for (auto & cell : msh ) {
        RealType h_l = diameter(msh, cell);
        if (h_l < h_min)      h_min = h_l;
        else if (h_l > h_max) h_max = h_l;
    }
    auto h_c = 0.75*h_max;
    auto p = h_max/h_min;
    std::cout << bold << cyan << "      h_max = "       << h_max << reset << std::endl;
    std::cout << bold << cyan << "      h_min = "       << h_min << std::endl;
    std::cout << bold << cyan << "      h_max/h_min = " << p     << reset << std::endl << std::endl;

    // ######################################################################
    // ###################################################################### Time controls
    // ######################################################################
    // dt is used as reference upper bound for the stability sweep (dt_stab_max = 2*dt).

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
        
        auto cell_faces = faces(msh,cell);
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
    std::cout << bold << red   << "   ASSEMBLY 1 : " << std::endl;
    std::cout << bold << cyan  << "      Assembler generation : " << tc << " seconds" << reset << std::endl;
    
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
    simulation_log << "Number of ERK steps =  " << 4       << std::endl;
    simulation_log << "Number of time steps =  " << nt      << std::endl;
    simulation_log << "Step size =  "             << dt      << std::endl;
    simulation_log << "Number of equations : "    << assembler.RHS.rows() << std::endl;
    simulation_log.flush();

    // ######################################################################
    // ###################################################################### Projection matrix P
    // ######################################################################

    assembler.assemble_P(msh, h_c);

    // ######################################################################
    // ###################################################################### Stability analysis
    // ######################################################################
    //
    // OBJECTIF : déterminer le rayon spectral rho(C_LTS) de la matrice
    // d'amplification du schéma LTS-RK4, en fonction de dt et de p.
    //
    // CRITERE DE STABILITE : le schéma est stable si et seulement si
    //   rho(C_LTS) <= 1
    // car y_{n+1} = C_LTS * y_n => ||y_n|| ~ rho^n ||y_0||.
    // Si rho > 1, la solution croît exponentiellement : instabilité.
    //
    // PRINCIPE GENERAL :
    // Le schéma LTS-RK4 avec F=0 est linéaire en y_n, donc il existe
    // une matrice C_LTS telle que y_{n+1} = C_LTS * y_n pour tout y_n.
    // On ne dispose pas de C_LTS explicitement, mais on sait calculer
    // le produit C_LTS * x pour n'importe quel vecteur x : c'est
    // exactement ce que fait un pas de temps du schéma LTS-RK4 avec F=0.
    // On exploite cela pour construire C_LTS colonne par colonne.

    std::cout << bold << red << "\n   STABILITY ANALYSIS (LTS-RK4)" << reset << std::endl;
    simulation_log << "\n============================================\n";
    simulation_log << "STABILITY ANALYSIS (LTS-RK4)\n";
    simulation_log << "============================================\n";

    // Valeurs de p à tester :
    //   p=1  => RK4 standard sans LTS (référence)
    //   p>1  => LTS-RK4(p) avec p sous-pas locaux
    const std::vector<int> p_values = {1, 2, 3, 5};

    // Plage de dt pour le sweep de stabilité.
    // On part de dt_stab_min (très stable) et on monte jusqu'à 2*dt_sim
    // pour détecter le seuil d'instabilité.
    const double dt_stab_min  = 1e-4;
    const double dt_stab_max  = 2.0 * dt;   // 2x le dt de simulation
    const int    nb_dt_points = 60;          // résolution du sweep

    const int    n_dof  = static_cast<int>(x_dof.rows()); // taille du système
    const int    n_eigs = 1;   // on cherche uniquement la valeur propre dominante
    // n_cv : taille de l'espace de Krylov pour Spectra.
    // Doit vérifier n_eigs < n_cv <= n_dof.
    // Plus n_cv est grand, plus la convergence est robuste mais coûteuse.
    const int    n_cv   = std::min(n_dof, std::max(30, 6 * n_eigs));

    const double ddt = (dt_stab_max - dt_stab_min) / static_cast<double>(nb_dt_points - 1);

    // ----------------------------------------------------------------
    // LAMBDA apply_C_LTS
    // ----------------------------------------------------------------
    // Applique un pas LTS-RK4 complet (F=0) au vecteur x_in et retourne
    // le résultat x_out = C_LTS * x_in.
    //
    // C'est le même algorithme que la boucle temporelle, mais :
    //   - F = 0 partout (on étudie la stabilité, pas la précision)
    //   - on travaille sur un vecteur quelconque x_in (pas forcément
    //     une condition initiale physique)
    //
    // Paramètres :
    //   x_in   : vecteur d'entrée (un vecteur de base e_i en pratique)
    //   dtau   : pas de temps local = dt / p
    //   p_stab : nombre de sous-pas fins
    // ----------------------------------------------------------------
    auto apply_C_LTS = [&](const Matrix<RealType, Dynamic, 1> & x_in,
                            double dtau, int p_stab)
        -> Matrix<RealType, Dynamic, 1>
    {
        // ------ Phase 1 : prédicateur coarse ------
        // Calcule les vecteurs w[0]..w[3] = contributions coarse gelées.
        // w[i] = B(I-P) B^i x_in  (avec F=0)
        // Ces vecteurs encodent les dérivées successives de la solution
        // sur la région coarse, évaluées une seule fois au début du pas.
        std::vector<Matrix<RealType, Dynamic, 1>> w(4);
        for (int i = 0; i < 4; ++i) { w[i].resize(n_dof); w[i].setZero(); }
        erk_an.erk_weight_LTS_coarse_old(x_in, assembler.Pcoarse, w);

        // ------ Phase 2 : sous-pas fins ------
        // On avance la solution sur p_stab sous-pas de taille dtau.
        // A chaque sous-pas m, on effectue un RK4 classique sur la
        // partie fine (projection P), synchronisé avec la contribution
        // coarse via le polynôme de Taylor en t.
        Matrix<RealType, Dynamic, 1> x = x_in;

        for (int m = 0; m < p_stab; ++m) {

            // Temps locaux pour les 4 stades RK4 :
            //   c1 = 0   => t = m*dtau
            //   c2 = 1/2 => t = (m+1/2)*dtau
            //   c3 = 1/2 => t = (m+1/2)*dtau
            //   c4 = 1   => t = (m+1)*dtau
            const double tm   =  m        * dtau;
            const double tm12 = (m + 0.5) * dtau;
            const double tm1  = (m + 1.0) * dtau;

            // Polynôme de Taylor coarse évalué au temps t :
            //   sum_{j=0}^{3} t^j/j! * w[j]
            // C'est la contribution coarse interpolée au bon instant
            // pour chaque stade RK4 (synchronisation coarse/fine).
            auto coarse_at = [&](double t) -> Matrix<RealType, Dynamic, 1> {
                return w[0] + t*w[1] + (t*t/2.0)*w[2] + (t*t*t/6.0)*w[3];
            };

            // Les yn sont matérialisés explicitement (ne pas passer
            // une expression temporaire à erk_weight qui prend un
            // argument par référence non-const).
            Matrix<RealType, Dynamic, 1> yn1, yn2, yn3, yn4;
            Matrix<RealType, Dynamic, 1> k1, k2, k3, k4;

            // -- Stade 1 : évaluation en t = m*dtau --
            // yn1 = P * x  (partie fine de x)
            // k1  = B(yn1) + coarse_at(tm)  (fine + coarse synchronisés)
            yn1 = assembler.Pfine * x;
            erk_an.erk_weight(yn1, k1);   // k1_fine = B * yn1
            k1 += coarse_at(tm);           // synchronisation coarse/fine

            // -- Stade 2 : évaluation en t = (m+1/2)*dtau --
            // On avance x d'un demi-pas avec k1, puis on prend la partie fine
            yn2 = assembler.Pfine * (x + dtau/2.0*k1);
            erk_an.erk_weight(yn2, k2);
            k2 += coarse_at(tm12);

            // -- Stade 3 : évaluation en t = (m+1/2)*dtau --
            // Même instant que le stade 2, mais avancé avec k2
            yn3 = assembler.Pfine * (x + dtau/2.0*k2);
            erk_an.erk_weight(yn3, k3);
            k3 += coarse_at(tm12);

            // -- Stade 4 : évaluation en t = (m+1)*dtau --
            // On avance x d'un pas entier avec k3
            yn4 = assembler.Pfine * (x + dtau*k3);
            erk_an.erk_weight(yn4, k4);
            k4 += coarse_at(tm1);

            // -- Mise à jour RK4 classique --
            // x_{m+1} = x_m + dtau/6 * (k1 + 2*k2 + 2*k3 + k4)
            x += dtau/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
        }

        // x contient maintenant C_LTS * x_in
        return x;
    };

    // ----------------------------------------------------------------
    // BOUCLE EXTERNE : sur les valeurs de p
    // Pour chaque p, on trace la courbe rho(C_LTS) en fonction de dt
    // et on identifie le dt_max stable.
    // ----------------------------------------------------------------
    for (int p_stab : p_values) {

        std::cout << bold << red
                  << "\n   === p = " << p_stab << " ===" << reset << std::endl;
        simulation_log << "\n--- p = " << p_stab << " ---\n";
        simulation_log << std::setw(18) << "dt"
                       << std::setw(18) << "rho(C_LTS)"
                       << std::setw(12) << "stable" << "\n";

        // dt_max_stable : le plus grand dt testé tel que rho <= 1.
        // Initialisé à -1 pour détecter le cas où aucun dt n'est stable.
        double dt_max_stable = -1.0;

        // ----------------------------------------------------------------
        // BOUCLE INTERNE : sweep sur dt
        // ----------------------------------------------------------------
        for (int s = 0; s < nb_dt_points; ++s) {

            const double dt_s = dt_stab_min + s * ddt;  // dt courant
            const double dtau = dt_s / p_stab;           // pas de temps local

            // ------------------------------------------------------------
            // CONSTRUCTION EXPLICITE DE C_LTS
            //
            // Idée : le schéma LTS-RK4 est linéaire en y_n (F=0), donc
            // il existe une matrice C_LTS telle que y_{n+1} = C_LTS * y_n.
            //
            // Pour construire C_LTS, on utilise l'identité :
            //   C_LTS * e_i = i-ème colonne de C_LTS
            // où e_i = (0,...,0,1,0,...,0)^T est le i-ème vecteur de base.
            //
            // Preuve : par définition du produit matrice-vecteur,
            //   (C_LTS * e_i)_j = sum_k [C_LTS]_{jk} * (e_i)_k
            //                    = sum_k [C_LTS]_{jk} * delta_{ki}
            //                    = [C_LTS]_{ji}
            // donc C_LTS * e_i donne exactement la i-ème colonne.
            //
            // On applique donc apply_C_LTS à chaque e_i pour obtenir
            // la i-ème colonne, et on assemble C_LTS colonne par colonne.
            //
            // Coût : n_dof applications de apply_C_LTS
            //        = n_dof pas de temps LTS-RK4
            // ------------------------------------------------------------
            tc.tic();
            Eigen::SparseMatrix<double> C_LTS(n_dof, n_dof);
            {
                // On stocke C_LTS en format triplet (ligne, colonne, valeur)
                // pour construire la SparseMatrix efficacement à la fin.
                // On ignore les entrées < 1e-15 pour conserver la creusité :
                // C_LTS est creuse car le schéma ne couple que les DDL voisins.
                std::vector<Eigen::Triplet<double>> triplets;
                triplets.reserve(n_dof); // estimation basse, sera étendu si besoin

                // e_i est réutilisé à chaque itération pour éviter des
                // réallocations inutiles : on met 1 à la position i,
                // on appelle apply_C_LTS, puis on remet 0.
                Matrix<RealType, Dynamic, 1> e_i = Matrix<RealType, Dynamic, 1>::Zero(n_dof);

                for (int i = 0; i < n_dof; ++i) {

                    // Construire e_i : vecteur nul sauf à la position i
                    e_i(i) = 1.0;

                    // Appliquer un pas LTS-RK4 à e_i :
                    // col = C_LTS * e_i = i-ème colonne de C_LTS
                    Matrix<RealType, Dynamic, 1> col = apply_C_LTS(e_i, dtau, p_stab);

                    // Stocker les entrées non nulles de la colonne i
                    // sous forme de triplets (j, i, col(j))
                    for (int j = 0; j < n_dof; ++j) {
                        if (std::abs(col(j)) > 1e-15)
                            triplets.emplace_back(j, i, col(j));
                    }

                    // Remettre e_i à zéro pour la prochaine itération
                    e_i(i) = 0.0;
                }

                // Assembler la SparseMatrix à partir des triplets
                C_LTS.setFromTriplets(triplets.begin(), triplets.end());
            }
            tc.toc();
            std::cout << bold << cyan << "      C_LTS built in " << tc << " s" << reset << std::endl;

            // ------------------------------------------------------------
            // CALCUL DU RAYON SPECTRAL rho(C_LTS)
            //
            // On utilise Spectra::GenEigsSolver (valeurs propres générales,
            // C_LTS est non-symétrique) avec SparseGenMatProd comme wrapper,
            // exactement comme dans compute_eigenvalues existant.
            //
            // On ne cherche que la valeur propre de plus grand module
            // (SortRule::LargestMagn), ce qui suffit pour rho.
            //
            // rho = max_i |lambda_i|  (les lambda_i sont complexes en général)
            //
            // Critère de stabilité :
            //   rho <= 1  =>  stable
            //   rho >  1  =>  instable (la solution diverge)
            // ------------------------------------------------------------
            Spectra::SparseGenMatProd<double> op(C_LTS);
            Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>>
                    eigs(op, n_eigs, n_cv);
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestMagn);

            const bool ok  = (eigs.info() == Spectra::CompInfo::Successful);
            double     rho = -1.0;
            if (ok) {
                // eigenvalues() retourne un vecteur complexe
                // cwiseAbs() prend le module de chaque entrée
                // maxCoeff() retourne le maximum => rho(C_LTS)
                rho = eigs.eigenvalues().cwiseAbs().maxCoeff();
                // Mettre à jour dt_max_stable tant que le schéma est stable
                if (rho <= 1.0) dt_max_stable = dt_s;
            }

            std::cout << bold << cyan
                      << "      dt = " << std::setw(12) << std::setprecision(6) << dt_s
                      << "   rho = "   << std::setw(12) << std::setprecision(6) << rho;
            if (!ok)            std::cout << "  [Spectra failed]";
            if (ok && rho<=1.0) std::cout << "  [stable]";
            else                std::cout << "  [UNSTABLE]";
            std::cout << reset << std::endl;

            simulation_log << std::setw(18) << dt_s
                           << std::setw(18) << rho
                           << std::setw(12) << (ok && rho<=1.0 ? "yes" : "no") << "\n";
        }

        // Le plus grand dt testé tel que rho(C_LTS) <= 1
        std::cout << bold << red
                  << "   --> dt_max_stable (p=" << p_stab << ") = " << dt_max_stable
                  << reset << std::endl;
        simulation_log << "dt_max_stable (p=" << p_stab << ") = " << dt_max_stable << "\n";
    }
    simulation_log.flush();

    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << std::endl << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
}
