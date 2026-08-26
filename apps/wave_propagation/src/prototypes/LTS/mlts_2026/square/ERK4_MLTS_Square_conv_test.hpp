
//  Created by Romain Mottier
//
#include <filesystem>
#include <functional>
//
// Convergence-order check for the multi-level restricted driver
// (ERK4_MLTS_Lshape_MehlinRestricted_conv_test.hpp's own machinery,
// verbatim: build_LTS_subblock_Pell + erk_weight_LTS_coarse_v2/_fine_v2
// + erk_weight_LTS_coarse_advance_uncovered_full), applied to the SAME
// proven graded-quadtree mesh generator as the L-shape drivers (see
// mesh_generation/square_graded/generate_square_graded_mesh.py, a
// direct adaptation of generate_lshape_graded_mesh.py with the notch
// removed -- same quadtree/hanging-node machinery, same poly_2d file
// format, same polygon_2d_mesh_reader path), grading toward the domain
// CENTER (an ordinary point, no singularity) instead of the L-shape's
// re-entrant corner, with a SMOOTH sinusoidal manufactured solution
// (scal_vec_analytic_functions::reproduction_acoustic: u = t^2 sin(pi x)
// sin(pi y), C-infinity everywhere, homogeneous Dirichlet on [0,1]^2 by
// construction) instead of the singular corner solution.
//
// PURPOSE: isolate whether the restricted drivers' residual (~5x on the
// L-shape corner problem, even with the P_ell-exact submatrix + overlap
// fix) comes from the algorithm itself or from the corner singularity's
// steep gradients interacting badly with a TRUNCATED local reach. Same
// band/threshold/CFL/recursion machinery, same p_global sequence
// (1, 32, 1024, 32768 for N=0..3, by construction -- the mesh generator
// is literally the same script) as the L-shape drivers; the ONLY things
// that change are the mesh's geometry (graded but non-singular) and the
// solution (smooth). If order ~k+1=4 is recovered here, the L-shape
// residual is specific to the singularity; if it is NOT recovered here
// either, the restricted algorithm itself still has a genuine bug
// unrelated to the singularity.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_MLTS_Square_conv_test(int argc, char **argv);

void ERK4_MLTS_Square_conv_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   MULTI-LEVEL LTS SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP" << std::endl << std::endl;

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

    const std::string mesh_dir = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/square_graded/meshes/";
    const size_t mesh_k = sim_data.m_k_degree;
    std::vector<int> level_N = {0, 1, 2, 3};
    if (const char* env_n = std::getenv("MLTS_N_MAX")) {
        int nmax = std::atoi(env_n);
        level_N.clear();
        for (int n = 0; n <= nmax; ++n) level_N.push_back(n);
    }

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);
    std::ostringstream conv_fname;
    conv_fname << out_dir << "/square_mlts_convergence_k_" << mesh_k << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# N  h_max  h_min  p_global  L_levels  L2_error_pressure  wall_time_s\n";

    for (int N : level_N) {

        std::cout << bold << red << "\n   ================ MESH LEVEL N=" << N << " ================" << reset << std::endl;

        std::ostringstream mesh_fname;
        mesh_fname << mesh_dir << "square_graded_k" << mesh_k << "_N" << N << ".txt";

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

        std::vector<RealType> thresholds;
        if (const char* env_l = std::getenv("MLTS_FORCE_L")) {
            int Lforce = std::atoi(env_l);
            for (int k = 1; k < Lforce; ++k)
                thresholds.push_back(h_max / std::pow(2.0, k));
        } else {
            RealType t_th = h_max / 4.0;
            while (t_th > 2.0 * h_min) { thresholds.push_back(t_th); t_th /= 4.0; }
        }
        const int L = (int)thresholds.size() + 1;
        std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p_global=" << p_global << "  L=" << L << " levels" << reset << std::endl;

        // ##################################################
        // Single-material (acoustic-only) setup, identical to the
        // L-shape drivers.
        // ##################################################

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

        // ##################################################
        // Smooth manufactured solution: u = t^2 sin(pi x) sin(pi y),
        // homogeneous Dirichlet on [0,1]^2 by construction (already
        // validated infrastructure, used e.g. in IAcoustic_conv_test.hpp).
        // `t` MUST be a single, shared, mutable outer variable:
        // Evaluate_s_*(double& t) captures t BY REFERENCE in the
        // returned closure, so every closure must be (re)created right
        // after `t` is set to the value it should evaluate at, and used
        // immediately -- the same pattern ERK4_LTS_conv_test.hpp relies on.
        // ##################################################

        RealType t = 0.0;
        scal_vec_analytic_functions functions;
        functions.set_function_type(scal_vec_analytic_functions::EFunctionType::reproduction_acoustic);

        auto s_u_fun0    = functions.Evaluate_s_u(t);
        auto s_v_fun0    = functions.Evaluate_s_v(t);
        auto s_flux_fun0 = functions.Evaluate_s_q(t);

        e_boundary_type e_bnd(msh);
        a_boundary_type a_bnd(msh);
        e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
        a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun0);

        auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
        assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
        assembler.set_coupling_stabilization();
        if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();

        assembler.assemble_mass(msh);
        assembler.assemble_coupling_terms(msh);

        Matrix<RealType, Dynamic, 1> x_dof;
        assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun0, s_flux_fun0);
        assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun0);

        assembler.assemble(msh, null_fun, s_v_fun0, true);
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

        // HYBRID architecture (see ERK4_MLTS_Lshape_HybridTerminal_conv_test.hpp
        // for the original derivation on the L-shape problem, re-validated
        // here on the smooth square problem): ancestor levels (i < L-1)
        // keep the GLOBAL exact computation (erk_weight_LTS_coarse_global,
        // full unrestricted Kcc/Kcf/Kfc) -- these are visited relatively
        // few times (product of ancestor p_i's), so their O(n_dof) cost
        // per call is affordable, and being GLOBAL they carry no
        // truncation error at all, regardless of L. Only the TERMINAL
        // level (visited nt*p_global times -- the dominant cost driver)
        // is restricted to a small own+halo sub-block via
        // erk_weight_LTS_fine_v2. This sidesteps the P_ell/overlap
        // truncation issue entirely for every non-terminal band: the
        // multi-level RESTRICTED approach's residual (confirmed via a
        // forced-L sweep on this exact problem: L=2 gives ~5e-8, L>=3
        // jumps by 3-4 orders of magnitude even with tuned overlap) was
        // traced to intermediate ("middle") bands needing a coarse-role
        // reach beyond what any FIXED-width halo/overlap can give --
        // going GLOBAL for those bands removes the truncation outright.
        std::vector<std::vector<int>> active_dofs_band(L);
        for (int b = 0; b < L; ++b) {
            active_dofs_band[b] = own_c[b];
            for (int g : own_f[b]) active_dofs_band[b].push_back((int)erk_an.n_c_dof() + g);
            std::sort(active_dofs_band[b].begin(), active_dofs_band[b].end());
        }
        const int terminal_halo = std::getenv("MLTS_TERMINAL_HALO") ? std::atoi(std::getenv("MLTS_TERMINAL_HALO")) : 8;
        auto blocks_L1 = erk_an.build_LTS_subblock(own_c[L-1], own_f[L-1], terminal_halo);
        std::set<int> term_cov_c(blocks_L1.active_c.begin(), blocks_L1.active_c.end());
        std::set<int> term_cov_f(blocks_L1.active_f.begin(), blocks_L1.active_f.end());
        for (int b = 0; b < L; ++b) {
            std::cout << bold << cyan << "      band" << b << ": own_c=" << own_c[b].size()
                      << " own_f=" << own_f[b].size() << reset << std::endl;
        }

        RealType cfl_factor = 0.05;
        if (const char* env_cfl = std::getenv("MLTS_CFL_FACTOR")) cfl_factor = std::atof(env_cfl);
        std::vector<RealType> scale(L), dt_level(L);
        scale[0] = h_max;
        for (int i = 1; i < L - 1; ++i) scale[i] = thresholds[i-1];
        scale[L-1] = h_min;
        for (int i = 0; i < L; ++i) dt_level[i] = cfl_factor * scale[i];
        std::vector<int> p_level(L, 1);
        for (int i = 0; i < L - 1; ++i)
            p_level[i] = std::max(1, (int)std::round(dt_level[i] / dt_level[i+1]));

        RealType ti = 0.0, tf = 0.1;
        RealType dt_macro = dt_level[0];
        size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_macro)));
        dt_macro = tf / nt;
        dt_level[0] = dt_macro;
        std::cout << bold << cyan << "      nt=" << nt << reset << std::endl;

        auto eval_F = [&](RealType t_abs) -> Matrix<RealType, Dynamic, 1> {
            t = t_abs;
            auto s_v_fun_t = functions.Evaluate_s_v(t);
            auto s_f_fun_t = functions.Evaluate_s_f(t);
            assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
            assembler.assemble_rhs(msh, null_fun, s_f_fun_t, true);
            return assembler.RHS;
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
                      const std::vector<Matrix<RealType,Dynamic,1>>& w_accum) {
            if (i == L - 1) {
                Matrix<RealType, Dynamic, 1> Fm = eval_F(t_start), Fmh = eval_F(t_start + 0.5*len), Fm1 = eval_F(t_start + len);
                erk_an.erk_weight_LTS_fine_v2(xn, blocks_L1, w_accum, Fm, Fmh, Fm1, 0.0, len);
                return;
            }
            int p_i = p_level[i];
            RealType dtau_i = len / p_i;
            RealType dt2 = dtau_i*dtau_i, dt3 = dt2*dtau_i, dt4 = dt3*dtau_i;
            for (int m = 0; m < p_i; ++m) {
                RealType tm = m * dtau_i;
                auto shifted = taylor_shift(w_accum, tm);
                std::vector<Matrix<RealType,Dynamic,1>> own_w_i(4);
                for (int j = 0; j < 4; ++j) own_w_i[j] = Matrix<RealType,Dynamic,1>::Zero(n_dof);
                Matrix<RealType, Dynamic, 1> Fn = eval_F(t_start + tm), Fn12 = eval_F(t_start + tm + 0.5*dtau_i), Fn1 = eval_F(t_start + tm + dtau_i);
                erk_an.erk_weight_LTS_coarse_global(xn, active_dofs_band[i], own_w_i, Fn, Fn12, Fn1, dtau_i);
                std::vector<Matrix<RealType,Dynamic,1>> combined(4);
                for (int j = 0; j < 4; ++j) combined[j] = shifted[j] + own_w_i[j];
                // Read-before-write ordering: recurse into finer levels
                // BEFORE writing band i's own closed-form advance into xn.
                // In the exact global algorithm xn is never written by any
                // ancestor level before its descendants have read it --
                // writing here first would let band (i+1)'s own_w_{i+1}
                // computation (erk_weight_LTS_coarse_global reads the FULL
                // current xn) see a xn already contaminated by band i's
                // direct write, a genuine ordering bug (see the analogous
                // fix in ERK4_MLTS_Lshape_RestrictedExact_conv_test.hpp).
                recurse(i + 1, t_start + tm, dtau_i, xn, combined);
                for (int g : own_c[i]) {
                    if (term_cov_c.count(g)) continue;
                    xn(g) += own_w_i[0](g)*dtau_i + own_w_i[1](g)*dt2/2.0 + own_w_i[2](g)*dt3/6.0 + own_w_i[3](g)*dt4/24.0;
                }
                for (int g : own_f[i]) {
                    if (term_cov_f.count(g)) continue;
                    int gg = (int)erk_an.n_c_dof() + g;
                    xn(gg) += own_w_i[0](gg)*dtau_i + own_w_i[1](gg)*dt2/2.0 + own_w_i[2](gg)*dt3/6.0 + own_w_i[3](gg)*dt4/24.0;
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

        t = tf;
        auto s_v_fun_T = functions.Evaluate_s_v(t);
        RealType l2_error_sq = 0.0;
        size_t cell_i = 0;
        for (auto & cell : msh) {
            Matrix<RealType, Dynamic, 1> p_dof = x_dof.block(cell_i * cell_dof + n_vec_dof, 0, n_scal_dof, 1);
            auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
            Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
            Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, s_v_fun_T);
            Matrix<RealType, Dynamic, 1> proj = mass.llt().solve(rhs);
            Matrix<RealType, Dynamic, 1> diff = proj - p_dof;
            l2_error_sq += diff.dot(mass * diff);
            cell_i++;
        }
        RealType l2_error = std::sqrt(l2_error_sq);

        std::cout << bold << red << "      wall time: " << tc_run << " s   L2_error(pressure) = "
                  << std::setprecision(10) << l2_error << reset << std::endl;

        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << p_global << " " << L
                  << " " << l2_error << " " << tc_run.elapsed() << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}
