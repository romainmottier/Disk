
//  Created by Romain Mottier
//
#include <filesystem>
#include <functional>
//
// Genuine multi-level (L>2) LTS-RK4 driver, built on the already-validated
// restricted building blocks (erk_weight_LTS_coarse_v2 / erk_weight_LTS_fine_v2
// / erk_weight_LTS_coarse_advance_uncovered in erk_coupling_hho_scheme.hpp,
// bit-exact validated at L=2 in ERK4_LTS_v2_L2_validation_test.hpp).
//
// KEY FIX vs the earlier (failed) "just shrink h_c" 2-level experiment
// (ERK4_LTS_v2_timing_test.hpp): a single coarse/fine split can only
// legitimately span ~1-2 octaves before the cubic-in-time Taylor polynomial
// wk[] stops being a valid approximation (confirmed empirically: NaN even
// at h_c=16*h_min while dt stayed tied to h_max). Here, EVERY level i gets
// its OWN dt_i, scaled to ITS OWN characteristic cell size (not the global
// h_max), so no level's cubic-Taylor fit is ever asked to represent more
// than a couple of octaves' worth of dynamics.
//
// Recursive structure per macro-step call recurse(0, tn, dt, x_dof, zero):
//   at level i (i<L-1): compute this level's OWN "coarse role" wk_i
//   (erk_weight_LTS_coarse_v2, using level i's own dt_i, NOT the global
//   macro dt), combine with the ancestor Taylor contribution passed in
//   (coefficient-wise sum -- both live in the SAME local-time coordinate
//   system, tau=0 <-> t_start), advance band i's own dofs NOT covered by
//   band i+1's halo via the closed-form integral of the COMBINED
//   polynomial (erk_weight_LTS_coarse_advance_uncovered), then loop p_i
//   times, each time Taylor-shifting the combined polynomial to that
//   sub-interval's own local origin and either calling erk_weight_LTS_fine_v2
//   directly (base case, i+1==L-1) or recursing (i+1<L-1) -- band i+1's own
//   advancement (both its "deep" cells via ITS OWN closed-form step, and its
//   "boundary" cells via whatever deeper level's fine_v2 eventually scatters
//   into them) is entirely handled INSIDE that recursive call, so there is
//   no separate/redundant fine_v2 call for band i+1 here.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_MLTS_timing_test(int argc, char **argv);

void ERK4_MLTS_timing_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   MULTI-LEVEL LTS TIMING/STABILITY TEST (N=2)" << std::endl << std::endl;

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
    int p_global = static_cast<int>(std::round(h_max / h_min));
    if (p_global < 1) p_global = 1;

    // ---------------- Multi-level bands: thresholds spanning ~2 octaves each ----------------
    std::vector<RealType> thresholds; // h_c_1 > h_c_2 > ... (L-1 thresholds)
    {
        RealType t = h_max / 4.0; // band0 = top 2 octaves
        while (t > 2.0 * h_min) {
            thresholds.push_back(t);
            t /= 4.0; // each subsequent band also ~2 octaves
        }
    }
    const int L = (int)thresholds.size() + 1;
    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << "  L=" << L << " levels" << reset << std::endl;
    for (size_t i = 0; i < thresholds.size(); ++i)
        std::cout << bold << cyan << "        threshold[" << i << "] = " << thresholds[i] << reset << std::endl;

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

    // ---------------- Build L bands (own_c/own_f) via successive threshold differencing ----------------
    size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
    size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
    size_t cell_dof   = n_scal_dof + n_vec_dof;
    size_t a_fbs      = disk::scalar_basis_size(hho_di.face_degree(), mesh_type::dimension - 1);
    std::vector<size_t> a_compress = assembler.get_a_compress();

    std::vector<int> cell_band(msh.cells_size());
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            RealType h = diameter(msh, cell);
            int band = L - 1;
            for (int i = 0; i < (int)thresholds.size(); ++i) if (h >= thresholds[i]) { band = i; break; }
            cell_band[ci] = band;
            ci++;
        }
    }
    // A face's band = the FINEST (highest band index) among its adjacent
    // cells -- matches assemble_P's union rule ("a face is fine if ANY
    // adjacent cell is fine"), generalized to L bands.
    std::vector<int> face_band(msh.faces_size(), 0);
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            for (auto & face : faces(msh, cell)) {
                auto fid = msh.lookup(face);
                if (cell_band[ci] > face_band[fid]) face_band[fid] = cell_band[ci];
            }
            ci++;
        }
    }

    std::vector<std::vector<int>> own_c(L), own_f(L);
    for (size_t i = 0; i < msh.cells_size(); ++i) {
        int b = cell_band[i];
        for (size_t d = 0; d < cell_dof; ++d) own_c[b].push_back((int)(i*cell_dof + d));
    }
    for (size_t j = 0; j < msh.faces_size(); ++j) {
        if (acoustic_bc_face_indexes.count(j)) continue;
        int b = face_band[j];
        size_t compressed = a_compress.at(j);
        for (size_t d = 0; d < a_fbs; ++d) own_f[b].push_back((int)(compressed*a_fbs + d));
    }

    std::vector<erk_coupling_hho_scheme<RealType>::LTS_subblock_set> blocks(L);
    for (int i = 0; i < L; ++i) {
        blocks[i] = erk_an.build_LTS_subblock(own_c[i], own_f[i], 6);
        std::cout << bold << cyan << "      band" << i << ": own_c=" << own_c[i].size()
                  << " own_f=" << own_f[i].size()
                  << " halo_c=" << (blocks[i].active_c.size() - blocks[i].n_own_c) << reset << std::endl;
    }

    // ---------------- Per-level scales/dt/ratios ----------------
    const RealType cfl_factor = 0.002;
    std::vector<RealType> scale(L), dt_level(L);
    scale[0] = h_max;
    for (int i = 1; i < L - 1; ++i) scale[i] = thresholds[i-1];
    scale[L-1] = h_min;
    for (int i = 0; i < L; ++i) dt_level[i] = cfl_factor * scale[i];
    std::vector<int> p_level(L, 1);
    for (int i = 0; i < L - 1; ++i) {
        p_level[i] = std::max(1, (int)std::round(dt_level[i] / dt_level[i+1]));
        std::cout << bold << cyan << "      level " << i << ": scale=" << scale[i] << " dt=" << dt_level[i] << " p=" << p_level[i] << reset << std::endl;
    }

    RealType ti = 0.0, tf = 0.05;
    RealType dt_macro = dt_level[0];
    size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_macro)));
    dt_macro = tf / nt;
    dt_level[0] = dt_macro; // re-snap level 0 to exactly divide [ti,tf]
    std::cout << bold << cyan << "      nt=" << nt << reset << std::endl;

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

    auto taylor_shift = [&](const std::vector<Matrix<RealType,Dynamic,1>>& wsrc, RealType a) {
        std::vector<Matrix<RealType,Dynamic,1>> out(4);
        RealType a2 = a*a, a3 = a2*a;
        out[0] = wsrc[0] + a*wsrc[1] + (a2/2.0)*wsrc[2] + (a3/6.0)*wsrc[3];
        out[1] = wsrc[1] + a*wsrc[2] + (a2/2.0)*wsrc[3];
        out[2] = wsrc[2] + a*wsrc[3];
        out[3] = wsrc[3];
        return out;
    };

    size_t n_dof = x_dof.rows();

    std::function<void(int, RealType, RealType, Matrix<RealType,Dynamic,1>&, const std::vector<Matrix<RealType,Dynamic,1>>&)> recurse;
    recurse = [&](int i, RealType t_start, RealType len, Matrix<RealType,Dynamic,1>& xn,
                  const std::vector<Matrix<RealType,Dynamic,1>>& ancestor_w) {

        std::vector<Matrix<RealType,Dynamic,1>> combined = ancestor_w;

        if (i < L - 1) {
            std::vector<Matrix<RealType,Dynamic,1>> own_w(4);
            for (int j = 0; j < 4; ++j) { own_w[j] = Matrix<RealType,Dynamic,1>::Zero(n_dof); }
            Matrix<RealType, Dynamic, 1> Fn   = eval_F(t_start);
            Matrix<RealType, Dynamic, 1> Fn12 = eval_F(t_start + 0.5*len);
            Matrix<RealType, Dynamic, 1> Fn1  = eval_F(t_start + len);
            erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w, Fn, Fn12, Fn1, len);
            for (int j = 0; j < 4; ++j) combined[j] = combined[j] + own_w[j];

            // Coverage must be checked against the TERMINAL (base-case)
            // band's halo, not the immediate next band's: only the
            // terminal band actually scatters via erk_weight_LTS_fine_v2
            // (a real write); every intermediate band only recurses, so
            // its own halo is read-only B-chain context, never a scatter
            // target. Using blocks[i+1] here for i<L-2 was a real bug
            // (found and fixed via ERK4_MLTS_L3_validation_test.hpp: it
            // dropped the L2-vs-reference error ratio from ~109x to ~9x
            // at a small controlled p_global=32 case).
            erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, blocks[i], combined,
                blocks[L-1].active_c, blocks[L-1].active_f, len);
        }

        if (i == L - 1) return;

        int p_i = p_level[i];
        RealType dtau_i = len / p_i;
        for (int m = 0; m < p_i; ++m) {
            RealType tm = m * dtau_i;
            auto shifted = taylor_shift(combined, tm);
            if (i + 1 == L - 1) {
                Matrix<RealType, Dynamic, 1> Fm  = eval_F(t_start + tm);
                Matrix<RealType, Dynamic, 1> Fmh = eval_F(t_start + tm + 0.5*dtau_i);
                Matrix<RealType, Dynamic, 1> Fm1 = eval_F(t_start + tm + dtau_i);
                erk_an.erk_weight_LTS_fine_v2(xn, blocks[i+1], shifted, Fm, Fmh, Fm1, 0.0, dtau_i);
            } else {
                recurse(i + 1, t_start + tm, dtau_i, xn, shifted);
            }
        }
    };

    timecounter tc_run;
    tc_run.tic();
    std::vector<Matrix<RealType,Dynamic,1>> zero_w(4);
    for (int j = 0; j < 4; ++j) zero_w[j] = Matrix<RealType,Dynamic,1>::Zero(n_dof);

    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt_macro * (it - 1);
        recurse(0, tn, dt_macro, x_dof, zero_w);
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
    std::cout << bold << red << "   L2_error(pressure, t=" << tf << ") = " << std::setprecision(10) << l2_error << reset << std::endl;
    std::cout << bold << yellow << "   (for reference, the original unrestricted 2-level scheme gave 1.14498093168055e-05 at N=2)" << reset << std::endl << std::endl;
}
