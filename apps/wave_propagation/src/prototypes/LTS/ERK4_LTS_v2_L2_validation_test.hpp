
//  Created by Romain Mottier
//
#include <filesystem>
//
// Validation gate for the new multi-level LTS building blocks
// (LTS_subblock_set, build_LTS_subblock, erk_weight_LTS_coarse_v2,
// erk_weight_LTS_fine_v2 in erk_coupling_hho_scheme.hpp): at L=2 (one
// "coarse" band = assembler.Pcoarse, one "fine" band = assembler.Pfine,
// exactly the existing partition), the new v2 functions must reproduce
// the trusted, already-validated ERK4_LTS_Lshape_MMS_conv_test.hpp result
// bit-for-bit (up to roundoff) -- N=1, k=3: L2 error = 2.76864549890122e-4.
// A mismatch here means a real bug in the v2 restriction/masking, to be
// fixed BEFORE attempting genuine multi-level (L>2).
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_LTS_v2_L2_validation_test(int argc, char **argv);

void ERK4_LTS_v2_L2_validation_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS v2 BUILDING-BLOCK VALIDATION (L=2, N=1)" << std::endl << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;
    typedef typename mesh_type::point_type point_type;

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    auto acoustic_mat_fun = [](const point_type&) -> acoustic_material_data<RealType> {
        return acoustic_material_data<RealType>(1.0, 1.0);
    };
    auto null_s_fun = [](const point_type&) -> RealType { return 0.0; };
    auto null_fun = [](const point_type&) -> disk::static_vector<double, 2> {
        return disk::static_vector<double, 2>{0, 0};
    };
    auto null_flux_fun = [](const point_type&) -> disk::static_matrix<double, 2, 2> {
        return disk::static_matrix<double, 2, 2>::Zero(2, 2);
    };

    const RealType corner_omega = 1.5 * M_PI;
    const RealType nu           = M_PI / corner_omega;
    const RealType w            = 2.0;

    auto polar = [](const point_type& pt) -> std::pair<RealType, RealType> {
        RealType dx = pt.x() - 0.5, dy = pt.y() - 0.5;
        RealType r = std::sqrt(dx * dx + dy * dy);
        RealType theta = std::atan2(dy, dx) - 0.5 * M_PI;
        if (theta < 0.0) theta += 2.0 * M_PI;
        return {r, theta};
    };
    auto phi = [&](RealType r, RealType theta) -> RealType {
        if (r < 1.0e-12) return 0.0;
        return std::pow(r, nu) * std::sin(nu * theta);
    };
    auto grad_phi = [&](RealType r, RealType theta) -> disk::static_vector<double, 2> {
        disk::static_vector<double, 2> g{0.0, 0.0};
        if (r < 1.0e-10) return g;
        RealType A = nu * std::pow(r, nu - 1.0) * std::sin(nu * theta);
        RealType B = nu * std::pow(r, nu - 1.0) * std::cos(nu * theta);
        RealType c = std::cos(theta), s = std::sin(theta);
        g[0] = -A * s - B * c;
        g[1] =  A * c - B * s;
        return g;
    };
    auto p_exact_at = [&](RealType t) {
        return [&, t](const point_type& pt) -> RealType {
            auto [r, theta] = polar(pt);
            return -w * std::sin(w * t) * phi(r, theta);
        };
    };
    auto v_exact_at = [&](RealType t) {
        return [&, t](const point_type& pt) -> disk::static_vector<double, 2> {
            auto [r, theta] = polar(pt);
            return std::cos(w * t) * grad_phi(r, theta);
        };
    };

    const std::string mesh_dir = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/";
    const int N = 1;
    const size_t mesh_k = sim_data.m_k_degree;

    std::ostringstream mesh_fname;
    mesh_fname << mesh_dir << "lshape_graded_k" << mesh_k << "_N" << N << ".txt";

    mesh_type msh;
    {
        polygon_2d_mesh_reader<RealType> mesh_builder;
        mesh_builder.set_poly_mesh_file(mesh_fname.str());
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
    }

    RealType h_max = 0.0, h_min = 1.0e10;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l > h_max) h_max = h_l;
        if (h_l < h_min) h_min = h_l;
    }
    int p = static_cast<int>(std::round(h_max / h_min));
    if (p < 1) p = 1;
    RealType h_c = (p == 1) ? 1.25 * h_max : 0.75 * h_max;

    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p=" << p << reset << std::endl;

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
    if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();

    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);

    Matrix<RealType, Dynamic, 1> x_dof;
    auto v0_fun = v_exact_at(0.0);
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v0_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);

    assembler.assemble(msh, null_fun, null_s_fun, true);
    assembler.LHS += assembler.COUPLING;

    std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
        assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(), assembler.get_e_face_dof(), assembler.get_a_face_dof());
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
        assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
        assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
        assembler.get_e_compress(), assembler.get_a_compress(),
        elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
    erk_an.refresh_faces_unknowns(x_dof);

    // nb_layer=0 (no buffer dilation), unlike the production nb_layer=3 used
    // by ERK4_LTS_Lshape_MMS_conv_test.hpp -- with nb_layer=3 on this small
    // N=1 mesh the dilation eats the ENTIRE coarse region (verified: 36/228
    // cells are geometrically >= h_c, but 3 rings of vertex-adjacency
    // dilation from the other 192 fine cells consumes all of them), so the
    // "trusted" 2.76864549890122e-4 reference never actually exercised a
    // nonempty coarse band -- worthless as a test of the coarse role /
    // halo / masking logic specifically. nb_layer=0 keeps a genuine
    // ~36-cell coarse region so this validation actually stresses the new
    // code's hardest part.
    assembler.assemble_P(msh, h_c, 0);

    // ---------------- Build the L=2 LTS_subblock_set pair from Pcoarse/Pfine ----------------
    auto extract_active = [&](const SparseMatrix<RealType>& P, size_t n_c_dof) {
        std::vector<int> c_idx, f_idx;
        for (int k = 0; k < P.outerSize(); ++k)
            for (typename SparseMatrix<RealType>::InnerIterator it(P, k); it; ++it) {
                int i = (int)it.row();
                if (i < (int)n_c_dof) c_idx.push_back(i);
                else f_idx.push_back(i - (int)n_c_dof);
            }
        return std::make_pair(c_idx, f_idx);
    };
    auto [coarse_c, coarse_f] = extract_active(assembler.Pcoarse, erk_an.n_c_dof());
    auto [fine_c,   fine_f]   = extract_active(assembler.Pfine,   erk_an.n_c_dof());
    std::cout << bold << cyan << "      band0 (coarse): " << coarse_c.size() << " cells, " << coarse_f.size() << " faces" << reset << std::endl;
    std::cout << bold << cyan << "      band1 (fine):   " << fine_c.size()   << " cells, " << fine_f.size()   << " faces" << reset << std::endl;

    auto blocks_coarse = erk_an.build_LTS_subblock(coarse_c, coarse_f, 6);
    auto blocks_fine   = erk_an.build_LTS_subblock(fine_c,   fine_f,   6);

    // ---------------- Time discretization (identical to the trusted MMS test) ----------------
    RealType ti = 0.0, tf = 0.05;
    const RealType cfl_factor = 0.002;
    RealType dt_target = cfl_factor * h_max;
    size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_target)));
    RealType dt = tf / nt;
    RealType dtau = dt / p;

    auto phi_field_fun = [&](const point_type& pt) -> RealType {
        auto [r, theta] = polar(pt);
        return phi(r, theta);
    };
    Matrix<RealType, Dynamic, 1> RHS_bc_unit, RHS_src_unit;
    {
        assembler.get_a_bc_conditions().updateDirichletFunction(phi_field_fun, 0);
        assembler.assemble_rhs(msh, null_fun, null_s_fun, true);
        RHS_bc_unit = assembler.RHS;

        assembler.get_a_bc_conditions().updateDirichletFunction(null_s_fun, 0);
        assembler.assemble_rhs(msh, null_fun, phi_field_fun, true);
        RHS_src_unit = assembler.RHS;
    }
    auto eval_F = [&](RealType t_abs) -> Matrix<RealType, Dynamic, 1> {
        RealType g = -w * std::sin(w * t_abs);
        RealType f = -w * w * std::cos(w * t_abs);
        return g * RHS_bc_unit + f * RHS_src_unit;
    };

    // ---------------- DIAGNOSTIC: direct wk[] comparison, coarse role only ----------------
    {
        Matrix<RealType, Dynamic, 1> Fn   = eval_F(ti);
        Matrix<RealType, Dynamic, 1> Fn12 = eval_F(ti + 0.5 * dt);
        Matrix<RealType, Dynamic, 1> Fn1  = eval_F(ti + dt);

        std::vector<Matrix<RealType, Dynamic, 1>> wk_orig(4), wk_v2(4);
        for (int j = 0; j < 4; ++j) { wk_orig[j].resize(x_dof.rows()); wk_orig[j].setZero(); }
        for (int j = 0; j < 4; ++j) { wk_v2[j].resize(x_dof.rows()); wk_v2[j].setZero(); }

        erk_an.ZeroFc();
        erk_an.erk_weight_LTS_coarse(x_dof, assembler.Pcoarse, wk_orig, Fn, Fn12, Fn1, dt);
        erk_an.erk_weight_LTS_coarse_v2(x_dof, blocks_coarse, wk_v2, Fn, Fn12, Fn1, dt);

        std::cout << bold << red << "\n   [DIAGNOSTIC] wk[] comparison (coarse role only, first macro-step):" << reset << std::endl;
        for (int j = 0; j < 4; ++j) {
            Matrix<RealType, Dynamic, 1> diff = wk_orig[j] - wk_v2[j];
            RealType max_diff = diff.cwiseAbs().maxCoeff();
            RealType max_orig = wk_orig[j].cwiseAbs().maxCoeff();

            // Where is the worst mismatch, cell or face?
            int worst_idx = 0; RealType worst_val = 0;
            for (int i = 0; i < diff.rows(); ++i) if (std::abs(diff(i)) > worst_val) { worst_val = std::abs(diff(i)); worst_idx = i; }
            bool worst_is_face = (size_t)worst_idx >= erk_an.n_c_dof();
            bool worst_in_coarse_own = false, worst_in_coarse_halo = false;
            if (!worst_is_face) {
                for (size_t i = 0; i < blocks_coarse.n_own_c; ++i) if (blocks_coarse.active_c[i] == worst_idx) worst_in_coarse_own = true;
                for (size_t i = blocks_coarse.n_own_c; i < blocks_coarse.active_c.size(); ++i) if (blocks_coarse.active_c[i] == worst_idx) worst_in_coarse_halo = true;
            } else {
                int local_f = worst_idx - (int)erk_an.n_c_dof();
                for (size_t i = 0; i < blocks_coarse.n_own_f; ++i) if (blocks_coarse.active_f[i] == local_f) worst_in_coarse_own = true;
                for (size_t i = blocks_coarse.n_own_f; i < blocks_coarse.active_f.size(); ++i) if (blocks_coarse.active_f[i] == local_f) worst_in_coarse_halo = true;
            }
            std::cout << bold << cyan << "      w[" << j << "]: max|diff|=" << std::setprecision(6) << max_diff
                      << "  max|orig|=" << max_orig
                      << "  worst_idx=" << worst_idx << (worst_is_face ? " (face)" : " (cell)")
                      << "  in_coarse_own=" << worst_in_coarse_own << "  in_coarse_halo=" << worst_in_coarse_halo
                      << reset << std::endl;
        }
    }

    // ---------------- Time marching: ORIGINAL (trusted) vs NEW v2, side by side ----------------
    Matrix<RealType, Dynamic, 1> x_dof_orig = x_dof;
    Matrix<RealType, Dynamic, 1> x_dof_v2   = x_dof;

    for (size_t it = 1; it <= nt; it++) {

        RealType tn = ti + dt * (it - 1);

        // ---- original (unrestricted, trusted) algorithm ----
        {
            auto x_dof_n = x_dof_orig;
            std::vector<Matrix<RealType, Dynamic, 1>> wk(4);
            for (int j = 0; j < 4; ++j) { wk[j].resize(x_dof_orig.rows()); wk[j].setZero(); }
            erk_an.ZeroFc();
            if (p != 1) {
                Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
                Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn + 0.5 * dt);
                Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn + dt);
                erk_an.ZeroFc();
                erk_an.erk_weight_LTS_coarse(x_dof_n, assembler.Pcoarse, wk, Fn, Fn12, Fn1, dt);
            }
            for (int m = 0; m < p; m++) {
                RealType tm  =  m        * dtau;
                RealType tmh = (m + 0.5) * dtau;
                RealType tm1 = (m + 1.0) * dtau;
                Matrix<RealType, Dynamic, 1> Fm  = eval_F(tn + tm);
                Matrix<RealType, Dynamic, 1> Fmh = eval_F(tn + tmh);
                Matrix<RealType, Dynamic, 1> Fm1 = eval_F(tn + tm1);
                erk_an.erk_weight_LTS_fine(x_dof_n, assembler.Pfine, wk, Fm, Fmh, Fm1, tm, dtau);
            }
            x_dof_orig = x_dof_n;
        }

        // ---- new v2 (restricted, multi-level-ready) algorithm ----
        {
            auto x_dof_n = x_dof_v2;
            std::vector<Matrix<RealType, Dynamic, 1>> wk(4);
            for (int j = 0; j < 4; ++j) { wk[j].resize(x_dof_v2.rows()); wk[j].setZero(); }
            if (p != 1) {
                Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
                Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn + 0.5 * dt);
                Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn + dt);
                erk_an.erk_weight_LTS_coarse_v2(x_dof_n, blocks_coarse, wk, Fn, Fn12, Fn1, dt);
                erk_an.erk_weight_LTS_coarse_advance_uncovered(x_dof_n, blocks_coarse, wk,
                    blocks_fine.active_c, blocks_fine.active_f, dt);
            }
            for (int m = 0; m < p; m++) {
                RealType tm  =  m        * dtau;
                RealType tmh = (m + 0.5) * dtau;
                RealType tm1 = (m + 1.0) * dtau;
                Matrix<RealType, Dynamic, 1> Fm  = eval_F(tn + tm);
                Matrix<RealType, Dynamic, 1> Fmh = eval_F(tn + tmh);
                Matrix<RealType, Dynamic, 1> Fm1 = eval_F(tn + tm1);
                erk_an.erk_weight_LTS_fine_v2(x_dof_n, blocks_fine, wk, Fm, Fmh, Fm1, tm, dtau);
            }
            x_dof_v2 = x_dof_n;
        }

        if (it == 1) {
            Matrix<RealType, Dynamic, 1> diff = (x_dof_orig - x_dof_v2).cwiseAbs();
            RealType d = diff.maxCoeff();
            std::cout << bold << cyan << "      [step " << it << "] max|orig-v2| = " << std::setprecision(6) << d << reset << std::endl;

            std::set<int> band0_own(coarse_c.begin(), coarse_c.end());
            std::set<int> band0_own_f(coarse_f.begin(), coarse_f.end());
            std::set<int> band1_own(fine_c.begin(), fine_c.end());
            std::set<int> band1_own_f(fine_f.begin(), fine_f.end());
            std::set<int> band1_halo_c(blocks_fine.active_c.begin() + blocks_fine.n_own_c, blocks_fine.active_c.end());
            std::set<int> band1_halo_f(blocks_fine.active_f.begin() + blocks_fine.n_own_f, blocks_fine.active_f.end());
            std::cout << bold << yellow << "      band1 halo size: " << band1_halo_c.size() << " cells, " << band1_halo_f.size() << " faces" << reset << std::endl;

            // Categorize every dof by its |diff| magnitude and band membership
            std::map<std::string, RealType> max_by_cat;
            std::map<std::string, int> count_by_cat;
            for (int i = 0; i < diff.rows(); ++i) {
                if (diff(i) < 1e-12) continue;
                bool is_face = (size_t)i >= erk_an.n_c_dof();
                int local = is_face ? (i - (int)erk_an.n_c_dof()) : i;
                std::string cat;
                if (!is_face) {
                    if (band0_own.count(local)) cat = "band0_own_cell";
                    else if (band1_own.count(local)) cat = "band1_own_cell";
                    else cat = "UNCLASSIFIED_cell";
                } else {
                    if (band0_own_f.count(local)) cat = "band0_own_face";
                    else if (band1_own_f.count(local)) cat = "band1_own_face";
                    else cat = "UNCLASSIFIED_face";
                }
                if (diff(i) > max_by_cat[cat]) max_by_cat[cat] = diff(i);
                count_by_cat[cat]++;
            }
            for (auto& [cat, mx] : max_by_cat) {
                std::cout << bold << cyan << "         " << cat << ": count=" << count_by_cat[cat] << "  max|diff|=" << std::setprecision(6) << mx << reset << std::endl;
            }
        }
    }

    // ---------------- Direct DOF-vector comparison (strongest possible check) ----------------
    RealType max_abs_diff = (x_dof_orig - x_dof_v2).cwiseAbs().maxCoeff();
    RealType max_abs_val  = x_dof_orig.cwiseAbs().maxCoeff();
    std::cout << bold << red << "\n   max|x_dof_orig - x_dof_v2|                 = "
              << std::setprecision(6) << max_abs_diff << reset << std::endl;
    std::cout << bold << red << "   max|x_dof_orig| (scale reference)          = "
              << std::setprecision(6) << max_abs_val << reset << std::endl;

    // ---------------- Error at t=tf vs closed-form solution, both runs ----------------
    size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
    size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
    size_t cell_dof   = n_scal_dof + n_vec_dof;

    auto compute_l2 = [&](const Matrix<RealType,Dynamic,1>& x_dof_final) -> RealType {
        auto p_exact_tf = p_exact_at(tf);
        RealType l2_error_sq = 0.0;
        size_t cell_i = 0;
        for (auto & cell : msh) {
            Matrix<RealType, Dynamic, 1> p_dof = x_dof_final.block(cell_i * cell_dof + n_vec_dof, 0, n_scal_dof, 1);
            auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
            Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
            Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, p_exact_tf);
            Matrix<RealType, Dynamic, 1> proj = mass.llt().solve(rhs);
            Matrix<RealType, Dynamic, 1> diff = proj - p_dof;
            l2_error_sq += diff.dot(mass * diff);
            cell_i++;
        }
        return std::sqrt(l2_error_sq);
    };

    RealType l2_orig = compute_l2(x_dof_orig);
    RealType l2_v2   = compute_l2(x_dof_v2);

    std::cout << bold << red   << "\n   L2_error(pressure, t=" << tf << ") [original, unrestricted] = " << std::setprecision(15) << l2_orig << reset << std::endl;
    std::cout << bold << yellow << "   L2_error(pressure, t=" << tf << ") [new v2, restricted]     = " << std::setprecision(15) << l2_v2   << reset << std::endl;

    bool gate_pass = (max_abs_diff < 1e-9 * std::max<RealType>(1.0, max_abs_val));
    std::cout << bold << (gate_pass ? green : red)
              << "   " << (gate_pass ? "MATCH (gate passes -- v2 building blocks are exact)" : "MISMATCH (gate fails -- do not proceed to L>2)")
              << reset << std::endl << std::endl;
}
