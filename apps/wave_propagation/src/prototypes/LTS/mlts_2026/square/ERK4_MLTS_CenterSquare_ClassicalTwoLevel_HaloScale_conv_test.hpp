
//  Created by Romain Mottier
//
#include <filesystem>
#include <functional>
#include <sys/resource.h>
//
// Classical-vs-multilevel EXACT comparison: same GENUINE classical
// two-level LTS-RK4 machinery as ERK4_MLTS_CenterSquare_ClassicalTwoLevel_
// conv_test.hpp (assemble_P threshold split + plain erk_weight_LTS_coarse/
// fine, no band recursion), but reading the SAME mesh family
// (centersquarehaloscale_graded, generate_haloscale_family.py) used by the
// multilevel ERK4_MLTS_CenterSquare_Adaptive_conv_test.hpp -- so both
// drivers run on IDENTICAL meshes (same h_max, same p_global per N),
// making wall/CPU time and accuracy directly comparable, not just
// order-of-magnitude.
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

static inline double cpu_seconds_now_classical() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec
         + (double)ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
}

void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_conv_test(int argc, char **argv);

void ERK4_MLTS_CenterSquare_ClassicalTwoLevel_HaloScale_conv_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   CLASSICAL TWO-LEVEL LTS-RK4 SQUARE SMOOTH-SOLUTION CONVERGENCE SWEEP" << std::endl << std::endl;

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
    auto null_fun = [](const point_type&) -> disk::static_vector<double, 2> {
        return disk::static_vector<double, 2>{0, 0};
    };
    auto null_flux_fun = [](const point_type&) -> disk::static_matrix<double, 2, 2> {
        return disk::static_matrix<double, 2, 2>::Zero(2, 2);
    };

    const std::string mesh_dir = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/center_square_fixed/meshes/";
    const size_t mesh_k = sim_data.m_k_degree;
    std::vector<int> level_N = {0, 1, 2, 3, 4};
    if (const char* env_n = std::getenv("MLTS_N_MAX")) {
        int nmax = std::atoi(env_n);
        level_N.clear();
        for (int n = 0; n <= nmax; ++n) level_N.push_back(n);
    }

    // nb_layer: vertex-sharing dilation passes around the fine region in
    // assemble_P -- ERK4_LTS_Lshape_conv_test.hpp's comment explains a
    // first-order mixed acoustic system needs more than the default 1 as p
    // grows; the haloscale family reaches p_global=1024, so default to 3
    // for safety (overridable).
    const int nb_layer = std::getenv("MLTS_NB_LAYER") ? std::atoi(std::getenv("MLTS_NB_LAYER")) : 3;

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);
    std::ostringstream conv_fname;
    conv_fname << out_dir << "/centersquare_classicaltwolevel_haloscale_convergence_k_" << mesh_k << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# N  h_max  h_min  p_global  L_levels  L2_error_pressure  wall_time_s  L2_error_velocity  cpu_time_s\n";

    for (int N : level_N) {

        std::cout << bold << red << "\n   ================ MESH LEVEL N=" << N << " (classical two-level, haloscale mesh) ================" << reset << std::endl;

        // ---------------- Mesh ----------------
        std::ostringstream mesh_fname;
        mesh_fname << mesh_dir << "centersquarehaloscale_graded_k" << mesh_k << "_N" << N << ".txt";

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
        RealType h_c = (p_global == 1) ? 1.25 * h_max : 0.75 * h_max;

        std::cout << bold << cyan << "      n_cells=" << msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p_global=" << p_global << reset << std::endl;

        // ---------------- Materials & boundary conditions ----------------
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

        // ---------------- Assembly ----------------
        auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
        assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
        assembler.set_coupling_stabilization();
        if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();

        assembler.assemble_mass(msh);
        assembler.assemble_coupling_terms(msh);

        // ---------------- Initial condition ----------------
        Matrix<RealType, Dynamic, 1> x_dof;
        assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun0, s_flux_fun0);
        assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun0);

        assembler.assemble(msh, null_fun, s_v_fun0, true);
        assembler.LHS += assembler.COUPLING;

        // ---------------- ERK operator setup ----------------
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

        // Classical two-level P projectors: threshold cell diameter vs h_c,
        // NOT the multilevel Pell/band index bookkeeping used elsewhere in
        // this investigation.
        assembler.assemble_P(msh, h_c, nb_layer);

        // ---------------- Time discretization ----------------
        RealType cfl_factor = 0.05;
        if (const char* env_cfl = std::getenv("MLTS_CFL_FACTOR")) cfl_factor = std::atof(env_cfl);
        RealType ti = 0.0, tf = 0.1;
        RealType dt_macro = cfl_factor * h_max;
        size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_macro)));
        dt_macro = tf / nt;
        RealType dtau = dt_macro / p_global;
        std::cout << bold << cyan << "      nt=" << nt << reset << std::endl;

        // ---------------- Source term ----------------
        auto eval_F = [&](RealType t_abs) -> Matrix<RealType, Dynamic, 1> {
            t = t_abs;
            auto s_v_fun_t = functions.Evaluate_s_v(t);
            auto s_f_fun_t = functions.Evaluate_s_f(t);
            assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
            assembler.assemble_rhs(msh, null_fun, s_f_fun_t, true);
            return assembler.RHS;
        };

        size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
        size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
        size_t cell_dof   = n_scal_dof + n_vec_dof;
        size_t n_dof      = x_dof.rows();

        // ---------------- Time marching ----------------
        timecounter tc_run;
        tc_run.tic();
        RealType cpu_start = cpu_seconds_now_classical();
        for (size_t it = 1; it <= nt; it++) {
            RealType tn = ti + dt_macro * (it - 1);
            auto x_dof_n = x_dof;

            std::vector<Matrix<RealType,Dynamic,1>> w(4);
            for (int j = 0; j < 4; ++j) { w[j].resize(n_dof); w[j].setZero(); }

            erk_an.ZeroFc();
            if (p_global != 1) {
                RealType tn12 = tn + 0.5*dt_macro;
                RealType tn1  = tn + dt_macro;
                Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
                Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn12);
                Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn1);
                erk_an.ZeroFc();
                erk_an.erk_weight_LTS_coarse(x_dof_n, assembler.Pcoarse, w, Fn, Fn12, Fn1, dt_macro);
            }

            for (int m = 0; m < p_global; m++) {
                RealType tm  =  m      * dtau;
                RealType tmh = (m+0.5) * dtau;
                RealType tm1 = (m+1.0) * dtau;
                Matrix<RealType, Dynamic, 1> Fm  = eval_F(tn + tm);
                Matrix<RealType, Dynamic, 1> Fmh = eval_F(tn + tmh);
                Matrix<RealType, Dynamic, 1> Fm1 = eval_F(tn + tm1);
                erk_an.erk_weight_LTS_fine(x_dof_n, assembler.Pfine, w, Fm, Fmh, Fm1, tm, dtau);
            }

            x_dof = x_dof_n;
            if (it % std::max<size_t>(1, nt/10) == 0)
                std::cout << bold << yellow << "         step " << it << "/" << nt << reset << std::endl;
        }
        tc_run.toc();
        RealType cpu_time = cpu_seconds_now_classical() - cpu_start;

        // ---------------- Error computation (pressure + velocity) ----------------
        t = tf;
        auto s_v_fun_T    = functions.Evaluate_s_v(t);
        auto s_flux_fun_T = functions.Evaluate_s_q(t);

        // Pressure (scalar cell unknown, degree = cell_degree): mass-matrix
        // L2 projection of the exact field vs. the HHO cell dofs.
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

        // Velocity (the DUAL variable, at the START of each cell's dof
        // block, size n_vec_dof): NOT a plain vector-monomial field -- its
        // n_vec_dof coefficients are the gradients of a scalar basis at
        // hho_di.reconstruction_degree() (row 0, the constant, is skipped
        // since its gradient is zero). Computed by direct quadrature
        // (same convention as postprocessor.hpp's four-fields error
        // routines), not a mass-matrix projection.
        RealType v_l2_error_sq = 0.0;
        cell_i = 0;
        for (auto & cell : msh) {
            Matrix<RealType, Dynamic, 1> v_dof = x_dof.block(cell_i * cell_dof, 0, n_vec_dof, 1);
            auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
            auto int_rule = integrate(msh, cell, 2 * (hho_di.cell_degree() + 1));
            for (auto & qp : int_rule) {
                auto t_dphi = rec_basis.eval_gradients(qp.point());
                Matrix<RealType, 1, 2> vh = Matrix<RealType, 1, 2>::Zero();
                for (size_t i = 1; i < t_dphi.rows(); i++)
                    vh = vh + v_dof(i - 1) * t_dphi.block(i, 0, 1, 2);
                disk::static_vector<RealType, 2> v_exact = s_flux_fun_T(qp.point());
                Matrix<RealType, 1, 2> v_exact_row; v_exact_row << v_exact(0), v_exact(1);
                v_l2_error_sq += qp.weight() * (v_exact_row - vh).squaredNorm();
            }
            cell_i++;
        }
        RealType v_l2_error = std::sqrt(v_l2_error_sq);

        std::cout << bold << red << "      wall time: " << tc_run << " s   cpu time: " << cpu_time
                  << " s   L2_error(pressure) = " << std::setprecision(10) << l2_error
                  << "   L2_error(velocity) = " << v_l2_error << reset << std::endl;

        // Columns: N h_max h_min p_global L_levels L2_error_pressure
        // wall_time_s L2_error_velocity cpu_time_s -- cpu_time_s is the
        // exact process CPU time (user+system, via getrusage), immune to
        // wall-clock contention from unrelated processes.
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << p_global << " " << 2
                  << " " << l2_error << " " << tc_run.elapsed() << " " << v_l2_error << " " << cpu_time << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}
