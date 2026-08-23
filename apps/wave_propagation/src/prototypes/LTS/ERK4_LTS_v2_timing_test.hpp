
//  Created by Romain Mottier
//
#include <filesystem>
//
// Timing experiment for the validated, correctly cost-restricted 2-level
// LTS building blocks (erk_weight_LTS_coarse_v2 / erk_weight_LTS_fine_v2 /
// erk_weight_LTS_coarse_advance_uncovered in erk_coupling_hho_scheme.hpp,
// validated bit-exact against the original at N=1 in
// ERK4_LTS_v2_L2_validation_test.hpp), run on N=2 (p=1024) with a
// DELIBERATELY TIGHT h_c so the "fine" band covers only the last couple of
// octaves (near h_min), not almost the whole mesh -- unlike the production
// h_c=0.75*h_max used elsewhere. The number of fine RK4 sub-steps (p) is a
// property of h_min/h_max and does NOT change with this choice, but the
// PER-CALL cost (now genuinely restricted to the small fine band + halo,
// not the whole mesh) should drop sharply, giving a real wall-clock
// speedup even without genuine multi-level (L>2) recursion.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_LTS_v2_timing_test(int argc, char **argv);

void ERK4_LTS_v2_timing_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS v2 TIMING EXPERIMENT (N=2, tight h_c)" << std::endl << std::endl;

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
    const int N = 2;
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

    // Deliberately tight h_c -- isolate only the last couple of octaves
    // (near h_min) as "fine", instead of the production 0.75*h_max.
    RealType h_c = 16.0 * h_min;

    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p=" << p << "  h_c=" << h_c << " (4*h_min)" << reset << std::endl;

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

    // NOTE: deliberately NOT using assembler.assemble_P() here -- it rounds
    // both the cell diameter and hc to 3 decimal places
    // (elastoacoustic_four_fields_assembler.hpp:1658, "round(diameter*1000)/1000"),
    // which silently zeroes out any h_c below ~5e-4. Every graded mesh in
    // this session has h_min well under that (8.6e-5 at N=2, 1.3e-6 at
    // N=3), so a genuinely tight h_c (needed for real cost restriction)
    // always rounds to 0 there, classifying every cell "coarse" -- this
    // was confirmed directly (band1 came back with 0 dofs). Classifying
    // cells here instead, with a direct (unrounded) comparison, since
    // build_LTS_subblock only needs plain index lists, not assemble_P's
    // sparse Pcoarse/Pfine matrices.
    size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
    size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
    size_t cell_dof   = n_scal_dof + n_vec_dof;
    size_t a_fbs      = disk::scalar_basis_size(hho_di.face_degree(), mesh_type::dimension - 1);

    std::vector<bool> cell_is_fine(msh.cells_size(), false);
    {
        size_t ci = 0;
        for (auto & cell : msh) { if (diameter(msh, cell) < h_c) cell_is_fine[ci] = true; ci++; }
    }
    std::vector<bool> face_is_fine(msh.faces_size(), false);
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            if (cell_is_fine[ci]) for (auto & face : faces(msh, cell)) face_is_fine[msh.lookup(face)] = true;
            ci++;
        }
    }
    // NOTE: face dof offsets are NOT simply mesh_face_index*a_fbs -- the
    // assembler builds a COMPRESSED face index that EXCLUDES Dirichlet
    // faces entirely (elastoacoustic_four_fields_assembler.hpp:120-154,
    // "a_egdes" is built by skipping is_a_dirichlet(face) faces, then
    // m_a_compress_indexes maps mesh face id -> compressed offset). Using
    // the raw mesh face index instead produced out-of-range accesses
    // (confirmed via an Eigen assertion in a -UNDEBUG debug rebuild) --
    // must go through get_a_compress() and skip Dirichlet faces exactly
    // like the assembler does.
    std::vector<size_t> a_compress = assembler.get_a_compress();
    std::vector<int> coarse_c, fine_c, coarse_f, fine_f;
    for (size_t i = 0; i < msh.cells_size(); ++i)
        for (size_t d = 0; d < cell_dof; ++d)
            (cell_is_fine[i] ? fine_c : coarse_c).push_back((int)(i*cell_dof + d));
    for (size_t j = 0; j < msh.faces_size(); ++j) {
        if (acoustic_bc_face_indexes.count(j)) continue;  // Dirichlet face: no dof slot at all
        size_t compressed = a_compress.at(j);
        for (size_t d = 0; d < a_fbs; ++d)
            (face_is_fine[j] ? fine_f : coarse_f).push_back((int)(compressed*a_fbs + d));
    }
    std::cout << bold << cyan << "      band0 (coarse): " << coarse_c.size() << " cell-dofs" << reset << std::endl;
    std::cout << bold << cyan << "      band1 (fine):   " << fine_c.size()   << " cell-dofs (vs "
              << msh.cells_size() << " mesh cells total)" << reset << std::endl;

    timecounter tc;
    tc.tic();
    auto blocks_coarse = erk_an.build_LTS_subblock(coarse_c, coarse_f, 6);
    auto blocks_fine   = erk_an.build_LTS_subblock(fine_c,   fine_f,   6);
    tc.toc();
    std::cout << bold << cyan << "      band1 halo reach: " << (blocks_fine.active_c.size() - blocks_fine.n_own_c)
              << " cell-dofs beyond own  (block-build time: " << tc << " s)" << reset << std::endl;

    RealType ti = 0.0, tf = 0.05;
    const RealType cfl_factor = 0.002;
    RealType dt_target = cfl_factor * h_max;
    size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_target)));
    RealType dt = tf / nt;
    RealType dtau = dt / p;
    std::cout << bold << cyan << "      nt=" << nt << "  p=" << p << "  (total fine sub-steps = " << nt*p << ")" << reset << std::endl;

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

    timecounter tc_run;
    tc_run.tic();

    for (size_t it = 1; it <= nt; it++) {

        RealType tn = ti + dt * (it - 1);
        auto x_dof_n = x_dof;

        std::vector<Matrix<RealType, Dynamic, 1>> wk(4);
        for (int j = 0; j < 4; ++j) { wk[j].resize(x_dof.rows()); wk[j].setZero(); }

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

        x_dof = x_dof_n;

        if (it % std::max<size_t>(1, nt/10) == 0)
            std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
    }
    tc_run.toc();

    auto p_exact_tf = p_exact_at(tf);
    RealType l2_error_sq = 0.0;
    size_t cell_i = 0;
    for (auto & cell : msh) {
        Matrix<RealType, Dynamic, 1> p_dof = x_dof.block(cell_i * cell_dof + n_vec_dof, 0, n_scal_dof, 1);
        auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
        Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
        Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, p_exact_tf);
        Matrix<RealType, Dynamic, 1> proj = mass.llt().solve(rhs);
        Matrix<RealType, Dynamic, 1> diff = proj - p_dof;
        l2_error_sq += diff.dot(mass * diff);
        cell_i++;
    }
    RealType l2_error = std::sqrt(l2_error_sq);

    std::cout << bold << red << "\n   TOTAL TIME-MARCHING WALL TIME: " << tc_run << " s" << reset << std::endl;
    std::cout << bold << red << "   L2_error(pressure, t=" << tf << ") = " << std::setprecision(10) << l2_error << reset << std::endl << std::endl;
}
