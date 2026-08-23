
//  Created by Romain Mottier
//
#include <filesystem>
#include <functional>
//
// Rigorous validation of the GENERIC recursive driver at L=3 (the first
// case that actually exercises multi-level ancestor-accumulation: with
// L=2 there is only ever ONE "coarse role" contributor, so the
// Taylor-shift-and-combine logic is never exercised beyond a single
// term -- see ERK4_MLTS_L2_validation_test.hpp, which is bit-exact but
// structurally cannot catch an L>=3 accumulation bug).
//
// Reference: the ORIGINAL, unrestricted 2-level scheme with Pcoarse=band0,
// Pfine={band1+band2 combined}. This is NOT the same algorithm as genuine
// L=3 (the 2-level reference gives band1 real RK4 dynamics at the tiny
// global dtau, while L=3 additionally approximates band1's OWN dynamics
// via a second cubic-Taylor fit) -- so an exact match is not expected.
// What IS expected, if the L=3 recursion has no bug: a SMALL, plausible
// discrepancy (comparable in magnitude to the accuracy loss already
// budgeted for by design), not a catastrophic (order-of-magnitude or
// blow-up) one. This test uses a deliberately CHEAP, SMALL p_global
// (N=1, p_global=32, split into 3 small bands) specifically so that
// "expected approximation accumulation" would predict only a mild
// discrepancy -- making a real implementation bug easy to distinguish
// from expected behaviour.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_MLTS_L3_validation_test(int argc, char **argv);

void ERK4_MLTS_L3_validation_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   MULTI-LEVEL DRIVER VALIDATION (L=3, N=1)" << std::endl << std::endl;

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
    int p_global = static_cast<int>(std::round(h_max / h_min));
    if (p_global < 1) p_global = 1;

    // L=3: two thresholds splitting the p_global=32 (5 octave) range into
    // 3 roughly-equal parts (~1.7 octaves each -- deliberately narrow
    // enough to exercise real multi-level accumulation, but the ABSOLUTE
    // p_global stays tiny so any bug shows up as a gross, unmissable
    // discrepancy rather than being masked by legitimate small
    // approximation error).
    // p_global=32=2^5; pick threshold[0] giving EXACT integer p0*p1=p_global
    // (p0=4, p1=8, since scale[1]=threshold[0]=h_max/4=8*h_min and
    // scale[2]=h_min by construction) so the comparison isn't confounded
    // by a quantization mismatch between the recursive scheme's total
    // fine-scale resolution and the reference's. threshold[1] (the actual
    // geometric band1/band2 cell-classification boundary) doesn't feed
    // into the p0/p1 calculation (that uses scale[2]=h_min directly), so
    // it can be picked independently -- just above h_min.
    std::vector<RealType> thresholds = { h_max / 4.0, 2.0 * h_min };
    const int L = 3;

    std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
              << "  h_max=" << h_max << "  h_min=" << h_min
              << "  p_global=" << p_global << reset << std::endl;

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
    std::vector<int> face_band(msh.faces_size(), 0);
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            for (auto & face : faces(msh, cell)) { auto fid = msh.lookup(face); if (cell_band[ci] > face_band[fid]) face_band[fid] = cell_band[ci]; }
            ci++;
        }
    }
    std::vector<std::vector<int>> own_c(L), own_f(L);
    for (size_t i = 0; i < msh.cells_size(); ++i) { int b = cell_band[i]; for (size_t d = 0; d < cell_dof; ++d) own_c[b].push_back((int)(i*cell_dof + d)); }
    for (size_t j = 0; j < msh.faces_size(); ++j) {
        if (acoustic_bc_face_indexes.count(j)) continue;
        int b = face_band[j]; size_t compressed = a_compress.at(j);
        for (size_t d = 0; d < a_fbs; ++d) own_f[b].push_back((int)(compressed*a_fbs + d));
    }
    std::vector<erk_coupling_hho_scheme<RealType>::LTS_subblock_set> blocks(L);
    for (int i = 0; i < L; ++i) {
        blocks[i] = erk_an.build_LTS_subblock(own_c[i], own_f[i], 6);
        std::cout << bold << cyan << "      band" << i << ": own_c=" << own_c[i].size() << " own_f=" << own_f[i].size()
                  << " halo_c=" << (blocks[i].active_c.size() - blocks[i].n_own_c) << reset << std::endl;
    }

    const RealType cfl_factor = 0.002;
    std::vector<RealType> scale(L), dt_level(L);
    scale[0] = h_max; scale[1] = thresholds[0]; scale[2] = h_min;
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
    dt_level[0] = dt_macro;
    std::cout << bold << cyan << "      nt=" << nt << "  (p0*p1=" << p_level[0]*p_level[1] << " vs p_global=" << p_global << ")" << reset << std::endl;

    auto phi_field_fun = [&](const point_type& pt) -> RealType { auto [r, theta] = polar(pt); return phi(r, theta); };
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
        RealType g = -w * std::sin(w * t_abs), f = -w * w * std::cos(w * t_abs);
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
    recurse = [&](int i, RealType t_start, RealType len, Matrix<RealType,Dynamic,1>& xn, const std::vector<Matrix<RealType,Dynamic,1>>& ancestor_w) {
        std::vector<Matrix<RealType,Dynamic,1>> combined = ancestor_w;
        if (i < L - 1) {
            std::vector<Matrix<RealType,Dynamic,1>> own_w(4);
            for (int j = 0; j < 4; ++j) own_w[j] = Matrix<RealType,Dynamic,1>::Zero(n_dof);
            Matrix<RealType, Dynamic, 1> Fn = eval_F(t_start), Fn12 = eval_F(t_start + 0.5*len), Fn1 = eval_F(t_start + len);
            erk_an.erk_weight_LTS_coarse_v2(xn, blocks[i], own_w, Fn, Fn12, Fn1, len);
            for (int j = 0; j < 4; ++j) combined[j] = combined[j] + own_w[j];
            // Coverage must be checked against the TERMINAL (base-case)
            // band's halo, not the immediate next band's: only the
            // terminal band actually scatters via erk_weight_LTS_fine_v2
            // (a real write); every intermediate band i+1<L-1 only
            // recurses -- its own halo is read-only context for ITS OWN
            // B-chain, never a scatter target. Using blocks[i+1] here for
            // i<L-2 was the bug: it treated "reachable by the next band's
            // halo" as "will be written by someone", which is only true
            // when that next band IS the terminal one.
            erk_an.erk_weight_LTS_coarse_advance_uncovered(xn, blocks[i], combined, blocks[L-1].active_c, blocks[L-1].active_f, len);
        }
        if (i == L - 1) return;
        int p_i = p_level[i];
        RealType dtau_i = len / p_i;
        for (int m = 0; m < p_i; ++m) {
            RealType tm = m * dtau_i;
            auto shifted = taylor_shift(combined, tm);
            if (i + 1 == L - 1) {
                Matrix<RealType, Dynamic, 1> Fm = eval_F(t_start + tm), Fmh = eval_F(t_start + tm + 0.5*dtau_i), Fm1 = eval_F(t_start + tm + dtau_i);
                erk_an.erk_weight_LTS_fine_v2(xn, blocks[i+1], shifted, Fm, Fmh, Fm1, 0.0, dtau_i);
            } else {
                recurse(i + 1, t_start + tm, dtau_i, xn, shifted);
            }
        }
    };

    // ---------------- reference: ORIGINAL 2-level scheme, Pcoarse=band0, Pfine={band1+band2} ----------------
    Matrix<RealType, Dynamic, 1> x_dof_ref  = x_dof;
    Matrix<RealType, Dynamic, 1> x_dof_mlts = x_dof;
    std::vector<Matrix<RealType,Dynamic,1>> zero_w(4);
    for (int j = 0; j < 4; ++j) zero_w[j] = Matrix<RealType,Dynamic,1>::Zero(n_dof);

    std::vector<Eigen::Triplet<RealType>> tc, tf_trip;
    for (int g : own_c[0]) tc.emplace_back(g, g, 1.0);
    for (int g : own_f[0]) tc.emplace_back((int)erk_an.n_c_dof()+g, (int)erk_an.n_c_dof()+g, 1.0);
    for (int g : own_c[1]) tf_trip.emplace_back(g, g, 1.0);
    for (int g : own_f[1]) tf_trip.emplace_back((int)erk_an.n_c_dof()+g, (int)erk_an.n_c_dof()+g, 1.0);
    for (int g : own_c[2]) tf_trip.emplace_back(g, g, 1.0);
    for (int g : own_f[2]) tf_trip.emplace_back((int)erk_an.n_c_dof()+g, (int)erk_an.n_c_dof()+g, 1.0);
    SparseMatrix<RealType> Pcoarse_ref(n_dof, n_dof), Pfine_ref(n_dof, n_dof);
    Pcoarse_ref.setFromTriplets(tc.begin(), tc.end());
    Pfine_ref.setFromTriplets(tf_trip.begin(), tf_trip.end());
    RealType dtau_ref = dt_macro / p_global;

    for (size_t it = 1; it <= nt; it++) {
        RealType tn = ti + dt_macro * (it - 1);

        {
            auto x_dof_n = x_dof_ref;
            std::vector<Matrix<RealType, Dynamic, 1>> wk(4);
            for (int j = 0; j < 4; ++j) { wk[j].resize(n_dof); wk[j].setZero(); }
            erk_an.ZeroFc();
            if (p_global != 1) {
                Matrix<RealType, Dynamic, 1> Fn = eval_F(tn), Fn12 = eval_F(tn + 0.5*dt_macro), Fn1 = eval_F(tn + dt_macro);
                erk_an.ZeroFc();
                erk_an.erk_weight_LTS_coarse(x_dof_n, Pcoarse_ref, wk, Fn, Fn12, Fn1, dt_macro);
                for (int m = 0; m < p_global; m++) {
                    RealType tm = m*dtau_ref, tmh = (m+0.5)*dtau_ref, tm1 = (m+1.0)*dtau_ref;
                    Matrix<RealType, Dynamic, 1> Fm = eval_F(tn+tm), Fmh = eval_F(tn+tmh), Fm1 = eval_F(tn+tm1);
                    erk_an.erk_weight_LTS_fine(x_dof_n, Pfine_ref, wk, Fm, Fmh, Fm1, tm, dtau_ref);
                }
            }
            x_dof_ref = x_dof_n;
        }

        recurse(0, tn, dt_macro, x_dof_mlts, zero_w);
    }

    RealType max_abs_diff = (x_dof_ref - x_dof_mlts).cwiseAbs().maxCoeff();
    RealType max_abs_val  = x_dof_ref.cwiseAbs().maxCoeff();

    {
        Matrix<RealType,Dynamic,1> diff = (x_dof_ref - x_dof_mlts).cwiseAbs();
        std::vector<RealType> max_by_band(L, 0.0);
        std::vector<std::set<int>> band_own_c(L), band_own_f(L);
        for (int b = 0; b < L; ++b) { band_own_c[b] = std::set<int>(own_c[b].begin(), own_c[b].end()); band_own_f[b] = std::set<int>(own_f[b].begin(), own_f[b].end()); }
        for (int i = 0; i < diff.rows(); ++i) {
            if (diff(i) < 1e-12) continue;
            bool is_face = (size_t)i >= erk_an.n_c_dof();
            int local = is_face ? (i - (int)erk_an.n_c_dof()) : i;
            for (int b = 0; b < L; ++b) {
                bool in = is_face ? band_own_f[b].count(local) : band_own_c[b].count(local);
                if (in) { if (diff(i) > max_by_band[b]) max_by_band[b] = diff(i); break; }
            }
        }
        for (int b = 0; b < L; ++b)
            std::cout << bold << cyan << "      max|diff| in band" << b << " own dofs: " << std::setprecision(6) << max_by_band[b] << reset << std::endl;
    }

    auto compute_l2 = [&](const Matrix<RealType,Dynamic,1>& xf) -> RealType {
        auto p_exact_tf = p_exact_at(tf);
        RealType l2_error_sq = 0.0; size_t cell_i = 0;
        for (auto & cell : msh) {
            Matrix<RealType, Dynamic, 1> p_dof = xf.block(cell_i * cell_dof + n_vec_dof, 0, n_scal_dof, 1);
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
    RealType l2_ref  = compute_l2(x_dof_ref);
    RealType l2_mlts = compute_l2(x_dof_mlts);

    std::cout << bold << red << "\n   max|x_dof_ref - x_dof_mlts_L3| = " << std::setprecision(6) << max_abs_diff << reset << std::endl;
    std::cout << bold << red << "   max|x_dof_ref| (scale)         = " << std::setprecision(6) << max_abs_val << reset << std::endl;
    std::cout << bold << red   << "\n   L2_error [2-level reference, Pfine=band1+band2] = " << std::setprecision(10) << l2_ref  << reset << std::endl;
    std::cout << bold << yellow << "   L2_error [generic recursive, L=3]                = " << std::setprecision(10) << l2_mlts << reset << std::endl;
    RealType ratio = l2_mlts / std::max(l2_ref, (RealType)1e-300);
    std::cout << bold << (ratio < 10.0 ? green : red)
              << "   ratio = " << ratio << "  (" << (ratio < 10.0 ? "plausible -- no gross bug" : "SUSPICIOUS -- likely a real bug in L>=3 accumulation")
              << reset << std::endl << std::endl;
}
