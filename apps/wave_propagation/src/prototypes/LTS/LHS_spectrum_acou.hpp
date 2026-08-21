
#include <filesystem>

//  Created by Romain Mottier
//
// Spectral analysis of the acoustic stiffness matrix (assembler.LHS) on a
// series of coarse meshes. LHS is square, size n_dof = n_c + n_f, with
// ordering [cell dofs | face dofs] (n_c = get_e_n_cells_dof() +
// get_a_n_cells_dof(), i.e. the *combined* elastic+acoustic cell block —
// here purely acoustic (V_T flux, P_T pressure) since only a_material is
// populated). The top-left n_c x n_c block LHS.topLeftCorner(n_c,n_c) is
// exactly the "cell-cell" block of the picture (zero corner + G/-G^T +
// stabilization restricted to P_T-P_T).
//
// LHS is NOT symmetric: reconstruction part is skew-symmetric (zero on the
// V_T-V_T corner, +/-G off-diagonal), stabilization part is symmetric but
// confined to the (P_T,P_F) sub-block. So we use the general EigenSolver
// and get complex eigenvalues/eigenvectors in general.
//
// For each eigenvalue we tag how much of its eigenvector's energy sits in
// the cell block vs the face block:
//     frac_cell = ||v.head(n_c)||^2 / (||v.head(n_c)||^2 + ||v.tail(n_f)||^2)
// frac_cell = 1  -> eigenvector purely supported on cell dofs (cell-cell block)
// frac_cell = 0  -> eigenvector purely supported on face dofs
// This lets the companion python plot color each point of the full LHS
// spectrum by its cell/face localization, highlighting the cell-cell block.
//
// In addition to the per-level spectrum files, a summary file
// "lhs_scaling_summary.txt" collects, for each mesh level: h_max, rho(LHS)
// (full spectrum), rho restricted to cell-dominated eigenvectors
// (frac_cell > 0.9), rho restricted to face-dominated eigenvectors
// (frac_cell < 0.1), and rho(S_SCHUR) (see below). This is meant to check,
// via a log-log fit (done in the companion python script), whether these
// spectral radii scale like h, 1/h or 1/h^2 as the mesh is refined.
//
// We also compute the Schur complement of LHS with respect to the face
// dofs, i.e. the n_c x n_c matrix obtained by eliminating faces via static
// condensation:
//     S_SCHUR = Kcc - Kcf * Sff_inv * Kfc
// (same formula/blocks as erk_coupling_hho_scheme::compute_eigenvalues_bis,
// reusing Kcc()/Kcf()/Kfc()/SffInv() populated by Sff_inverse()). Unlike the
// raw LHS.topLeftCorner(n_c,n_c) block, S_SCHUR is the operator that
// actually governs the cell unknowns once faces are condensed out — the
// physically relevant "cell-cell" operator used by the real explicit
// schemes in this codebase.
//
// Output layout, all under a "spectrum/" subdirectory of the current
// working directory:
//     spectrum/lhs_scaling_summary.txt   <- one row per level, both matrices
//     spectrum/dense/level_<lvl>.txt     <- full LHS spectrum (+ frac_cell)
//     spectrum/condensed/level_<lvl>.txt <- S_SCHUR spectrum (n_c x n_c)
//
// Usage example: ../../../wave_propagation -k0 -s0 -r0 -c0 -m0 -l0 -n11 -p1 -f1 -e0

void LHS_spectrum_acou(int argc, char **argv);

void LHS_spectrum_acou(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   LHS SPECTRUM ANALYSIS (ACOUSTIC)" << std::endl << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    const std::string out_dir       = "spectrum";
    const std::string dense_dir     = out_dir + "/dense";
    const std::string condensed_dir = out_dir + "/condensed";
    std::filesystem::create_directories(dense_dir);
    std::filesystem::create_directories(condensed_dir);

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    auto acoustic_mat_fun = [](const typename mesh_type::point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };
    auto null_s_fun = [](const mesh_type::point_type&) -> double {
        return 0.0;
    };
    auto null_fun = [](const mesh_type::point_type&) -> disk::static_vector<double, 2> {
        disk::static_vector<double, 2> f{0, 0};
        return f;
    };
    auto null_flux_fun = [](const mesh_type::point_type&) -> disk::static_matrix<double, 2, 2> {
        return disk::static_matrix<double, 2, 2>::Zero(2, 2);
    };

    // Series of coarse cartesian meshes (unit square), small enough for a
    // full dense eigendecomposition (Eigen::EigenSolver, O(n_dof^3)), spread
    // over close to a decade in h for a reliable log-log scaling fit.
    std::vector<std::pair<size_t, size_t>> mesh_levels = { {4, 4}, {6, 6}, {8, 8}, {12, 12}, {16, 16}, {24, 24} };

    std::ofstream summary(out_dir + "/lhs_scaling_summary.txt");
    summary << "# level  nx  ny  h_max  n_dof  n_c  n_f  rho_full  rho_cell  rho_face  rho_schur\n";

    for (size_t lvl = 0; lvl < mesh_levels.size(); ++lvl) {

        size_t nx = mesh_levels[lvl].first;
        size_t ny = mesh_levels[lvl].second;

        std::cout << bold << red << "\n   MESH LEVEL " << lvl << "  (nx=" << nx << ", ny=" << ny << ")" << reset << std::endl;

        RealType lx = 1.0, ly = 1.0;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.build_mesh();
        mesh_type msh;
        mesh_builder.move_to_mesh_storage(msh);

        RealType h_max = 0.0;
        for (auto & cell : msh) {
            RealType h_l = diameter(msh, cell);
            if (h_l > h_max) h_max = h_l;
        }

        // Pure acoustic domain: no elastic cells, no interface.
        std::map<size_t, acoustic_material_data<RealType>> a_material;
        std::map<size_t, elastic_material_data<RealType>>  e_material;
        std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
        std::map<size_t, std::pair<size_t, size_t>> interface_cell_pair_indexes;

        for (auto & cell : msh) {
            auto cell_ind = msh.lookup(cell);
            a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(barycenter(msh, cell))));
        }

        size_t bc_elastic_id  = 0;
        size_t bc_acoustic_id = 1;
        for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
            auto face  = *face_it;
            auto fc_id = msh.lookup(face);
            disk::boundary_descriptor bi{bc_acoustic_id, true};
            msh.backend_storage()->boundary_info.at(fc_id) = bi;
            acoustic_bc_face_indexes.insert(fc_id);
        }

        e_boundary_type e_bnd(msh);
        a_boundary_type a_bnd(msh);
        e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
        a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);

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
        assembler.assemble(msh, null_fun, null_s_fun, true);
        assembler.LHS += assembler.COUPLING;

        const int n_c   = static_cast<int>(assembler.get_e_n_cells_dof() + assembler.get_a_n_cells_dof());
        const int n_dof = static_cast<int>(assembler.LHS.rows());
        const int n_f   = n_dof - n_c;

        std::cout << bold << cyan << "      n_dof=" << n_dof << "  n_c=" << n_c << "  n_f=" << n_f << reset << std::endl;

        // -----------------------------------------------------------------
        // Schur complement w.r.t. face dofs: S_SCHUR = Kcc - Kcf*Sff_inv*Kfc
        // -----------------------------------------------------------------
        std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
        erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
            assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(), assembler.get_e_face_dof(), assembler.get_a_face_dof());
        erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
            assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
            assembler.get_e_compress(), assembler.get_a_compress(),
            elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);

        SparseMatrix<RealType> S_SCHUR = erk_an.Kcc() - erk_an.Kcf() * erk_an.SffInv() * erk_an.Kfc();
        Eigen::MatrixXd S_dense = S_SCHUR.toDense();
        Eigen::EigenSolver<Eigen::MatrixXd> es_schur(S_dense);

        double rho_schur = -1.0;
        if (es_schur.info() == Eigen::Success) {
            rho_schur = es_schur.eigenvalues().cwiseAbs().maxCoeff();
        } else {
            std::cout << bold << red << "      --> Schur EigenSolver FAILED at level " << lvl << reset << std::endl;
        }
        std::cout << bold << cyan << "      rho(S_SCHUR) = " << std::setprecision(10) << rho_schur << reset << std::endl;

        std::ostringstream schur_fname;
        schur_fname << condensed_dir << "/level_" << lvl << ".txt";
        std::ofstream schur_out(schur_fname.str());
        schur_out << "# level=" << lvl << " nx=" << nx << " ny=" << ny << " k=" << sim_data.m_k_degree
                  << " h_max=" << std::setprecision(15) << h_max << " n_c=" << n_c << "\n";
        schur_out << "# real  imag\n";
        if (es_schur.info() == Eigen::Success) {
            for (int j = 0; j < n_c; ++j) {
                schur_out << std::setprecision(15)
                          << es_schur.eigenvalues()(j).real() << "  "
                          << es_schur.eigenvalues()(j).imag() << "\n";
            }
        }
        schur_out.close();
        std::cout << bold << cyan << "      --> " << schur_fname.str() << reset << std::endl;

        // -----------------------------------------------------------------
        // Full dense eigendecomposition of LHS
        // -----------------------------------------------------------------
        Eigen::MatrixXd LHS_dense = assembler.LHS.toDense();
        Eigen::EigenSolver<Eigen::MatrixXd> es(LHS_dense);

        if (es.info() != Eigen::Success) {
            std::cout << bold << red << "      --> EigenSolver FAILED at level " << lvl << reset << std::endl;
            continue;
        }

        double rho = es.eigenvalues().cwiseAbs().maxCoeff();
        std::cout << bold << cyan << "      h_max=" << h_max << "  rho(LHS) = " << std::setprecision(10) << rho << reset << std::endl;

        // -----------------------------------------------------------------
        // Per-eigenvalue cell/face localization of the eigenvector, and
        // block-restricted spectral radii (cell-dominated / face-dominated)
        // -----------------------------------------------------------------
        std::ostringstream fname;
        fname << dense_dir << "/level_" << lvl << ".txt";
        std::ofstream out(fname.str());
        out << "# level=" << lvl << " nx=" << nx << " ny=" << ny << " k=" << sim_data.m_k_degree
            << " h_max=" << std::setprecision(15) << h_max
            << " n_dof=" << n_dof << " n_c=" << n_c << " n_f=" << n_f << "\n";
        out << "# real  imag  frac_cell\n";

        double rho_cell = 0.0, rho_face = 0.0;
        for (int j = 0; j < n_dof; ++j) {
            Eigen::VectorXcd v = es.eigenvectors().col(j);
            double norm_c2    = v.head(n_c).squaredNorm();
            double norm_f2    = v.tail(n_f).squaredNorm();
            double frac_cell  = norm_c2 / (norm_c2 + norm_f2 + 1.0e-300);
            double lambda_abs = std::abs(es.eigenvalues()(j));

            if (frac_cell > 0.9) rho_cell = std::max(rho_cell, lambda_abs);
            if (frac_cell < 0.1) rho_face = std::max(rho_face, lambda_abs);

            out << std::setprecision(15)
                << es.eigenvalues()(j).real() << "  "
                << es.eigenvalues()(j).imag() << "  "
                << frac_cell << "\n";
        }
        out.close();

        summary << std::setw(6)  << lvl
                << std::setw(6)  << nx
                << std::setw(6)  << ny
                << std::setw(20) << std::setprecision(15) << h_max
                << std::setw(10) << n_dof
                << std::setw(10) << n_c
                << std::setw(10) << n_f
                << std::setw(22) << std::setprecision(15) << rho
                << std::setw(22) << std::setprecision(15) << rho_cell
                << std::setw(22) << std::setprecision(15) << rho_face
                << std::setw(22) << std::setprecision(15) << rho_schur << "\n";
        summary.flush();

        std::cout << bold << cyan << "      rho_cell(frac>0.9)=" << rho_cell
                  << "  rho_face(frac<0.1)=" << rho_face << reset << std::endl;
        std::cout << bold << cyan << "      --> " << fname.str() << reset << std::endl;
    }

    summary.close();
    std::cout << bold << red << "\n   --> " << out_dir << "/lhs_scaling_summary.txt" << reset << std::endl;
    std::cout << std::endl;
}
