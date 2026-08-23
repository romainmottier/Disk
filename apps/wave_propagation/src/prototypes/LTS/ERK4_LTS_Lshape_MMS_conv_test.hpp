
//  Created by Romain Mottier
//
#include <filesystem>
//
// Method-of-manufactured-solutions convergence test for the actual
// time-dependent LTS-RK4 + HHO scheme (erk_coupling_hho_scheme, "Algorithm
// 3") on the graded L-shaped mesh with reentrant corner -- unlike
// EllipticLshape_conv_test.hpp (a static Poisson solve, no time-stepping at
// all) this exercises the real explicit LTS time integrator, and unlike
// ERK4_LTS_Lshape_conv_test.hpp (compares against a numerically-computed
// reference mesh solution) the error is measured against a CLOSED-FORM
// exact solution, so no reference mesh is needed and every level's error is
// exact rather than "vs finest available mesh".
//
// The pure acoustic wave equation u_tt - Delta u = f has no known closed
// form on the L-shaped domain in general (the domain's eigenmodes are not
// analytic), but a manufactured solution is easy to build by tacking a
// smooth time profile onto the static corner-singular harmonic function
// phi(r,theta) = r^{2/3} sin(2*theta/3) (same one used in
// EllipticLshape_conv_test.hpp, Delta phi = 0, vanishes on the two notch
// edges, r,theta polar coordinates centered at the reentrant corner
// (0.5,0.5), theta measured from the edge toward (0.5,1), interior angle
// omega_corner = 3*pi/2):
//
//   U(x,y,t) = cos(w t) * phi(r,theta)                (not itself a state
//                                                       variable, shown only
//                                                       for reference)
//   P(x,y,t) = dU/dt   = -w*sin(w t) * phi(r,theta)    (our acoustic
//                                                       "pressure" state, = udot)
//   V(x,y,t) = grad U  =  cos(w t) * grad(phi)         (our acoustic
//                                                       "flux/velocity" state)
//   S(x,y,t) = d2U/dt2 - Delta U = -w^2*cos(w t)*phi(r,theta) - 0
//                                                       (source injected into
//                                                       the pressure-equation
//                                                       RHS at every RK
//                                                       substage, following
//                                                       exactly the
//                                                       eval_F/assemble_rhs
//                                                       pattern already used
//                                                       for manufactured
//                                                       solutions in
//                                                       ERK4_LTS_conv_test.hpp)
//
// Since P(x,y,0) = 0 identically (sin(0)=0), the initial condition is simply
// p(0)=0, v(0)=grad(phi) -- no auxiliary Poisson solve needed for the IC,
// unlike ERK4_LTS_Lshape_conv_test.hpp.
//
// Convergence metric: L2 error of the pressure field p at t=tf against the
// closed-form P(x,y,tf), using the same mass-matrix-weighted L2-projection
// technique as postprocessor::compute_errors_one_field /
// EllipticLshape_conv_test.hpp (project the exact function onto each cell's
// own scalar basis, compare to that cell's actual pressure DOFs).
//
// Usage example: ../../../wave_propagation -k3 -s0 -r0 -c0 -f0 -e0

void ERK4_LTS_Lshape_MMS_conv_test(int argc, char **argv);

void ERK4_LTS_Lshape_MMS_conv_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ERK4-LTS L-SHAPE MMS CONV TEST (reentrant corner)" << std::endl << std::endl;

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

    // ---------------- corner-singular manufactured solution ----------------
    const RealType corner_omega = 1.5 * M_PI;      // 270 deg reentrant angle
    const RealType nu           = M_PI / corner_omega; // = 2/3
    const RealType w            = 2.0;             // time-oscillation frequency

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
        RealType A = nu * std::pow(r, nu - 1.0) * std::sin(nu * theta); // du/dr
        RealType B = nu * std::pow(r, nu - 1.0) * std::cos(nu * theta); // (1/r) du/dtheta
        RealType c = std::cos(theta), s = std::sin(theta);
        g[0] = -A * s - B * c;
        g[1] =  A * c - B * s;
        return g;
    };

    // p_exact(t), s_exact(t): closures over the running time, updated by
    // eval_F at every RK substage (same idiom as ERK4_LTS_conv_test.hpp).
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

    // -----------------------------------------------------------------
    // Mesh sequence (same coupled-refinement quadtree family as
    // EllipticLshape_conv_test.hpp / ERK4_LTS_Lshape_conv_test.hpp).
    // -----------------------------------------------------------------
    const std::string mesh_dir = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/";
    std::vector<int> level_N = {0, 1, 2};
    const size_t mesh_k = sim_data.m_k_degree;

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);

    std::ostringstream conv_fname;
    conv_fname << out_dir << "/lshape_mms_convergence_k_" << mesh_k << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# N  h_max  h_min  L2_error_pressure\n";

    for (int N : level_N) {

        std::cout << bold << red << "\n   MESH LEVEL N=" << N << reset << std::endl;

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

        // ---------------- Pure acoustic domain ----------------
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
        a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun); // p_exact(0) = 0

        auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
        assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
        assembler.set_coupling_stabilization();
        if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();

        assembler.assemble_mass(msh);
        assembler.assemble_coupling_terms(msh);

        // ---------------- Initial condition: p(0)=0, v(0)=grad(phi) ----------------
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

        assembler.assemble_P(msh, h_c, 3);

        // ---------------- Time discretization ----------------
        RealType ti = 0.0, tf = 0.05;
        const RealType cfl_factor = 0.002;
        RealType dt_target = cfl_factor * h_max;
        size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt_target)));
        RealType dt = tf / nt;
        RealType dtau = dt / p;

        // ---------------- Source-term lambda (eval_F) ----------------
        // S(x,y,t) = -w^2*cos(w t)*phi(x,y) and the Dirichlet trace
        // p_exact(x,y,t) = -w*sin(w t)*phi(x,y) are both separable in time
        // (a fixed spatial profile phi(x,y) times a scalar time factor), and
        // assemble_rhs is linear in both its source argument and the
        // currently-set Dirichlet function. So instead of calling the O(cells)
        // assemble_rhs at every one of the ~nt*(3+3*p) RK substages (which
        // times out already at N=2, p=1024), precompute the two unit RHS
        // vectors ONCE per level and recombine them with scalar time factors
        // -- exact, and O(dof) per substage instead of O(cells).
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
            RealType g = -w * std::sin(w * t_abs);       // Dirichlet-trace time factor
            RealType f = -w * w * std::cos(w * t_abs);   // source time factor
            return g * RHS_bc_unit + f * RHS_src_unit;
        };

        // ---------------- Time marching ----------------
        for (size_t it = 1; it <= nt; it++) {

            RealType tn = ti + dt * (it - 1);
            auto x_dof_n = x_dof;

            std::vector<Matrix<RealType, Dynamic, 1>> wk(4);
            for (int j = 0; j < 4; ++j) { wk[j].resize(x_dof.rows()); wk[j].setZero(); }

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

            x_dof = x_dof_n;
        }

        // ---------------- Error at t=tf vs closed-form solution ----------------
        size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
        size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
        size_t cell_dof   = n_scal_dof + n_vec_dof;

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

        std::cout << bold << cyan << "      L2_error(pressure, t=" << tf << ") = "
                  << std::setprecision(10) << l2_error << reset << std::endl;
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << l2_error << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}
