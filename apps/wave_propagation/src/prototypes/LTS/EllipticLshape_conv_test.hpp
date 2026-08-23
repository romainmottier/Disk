
//  Created by Romain Mottier
//
// Analytic-solution convergence test on the graded L-shaped mesh, decoupled
// from the wave-equation/LTS-RK4 machinery in ERK4_LTS_Lshape_conv_test.hpp.
//
// Domain Omega = (0,1)^2 \ [0.5,1)x(0.5,1], reentrant corner at H=(0.5,0.5),
// interior angle omega = 3*pi/2 (270 deg). The classical corner-singular
// harmonic function
//     u(r,theta) = r^{pi/omega} sin(pi*theta/omega) = r^{2/3} sin(2*theta/3)
// (r, theta: polar coordinates centered at H, theta measured from the edge
// H-(0.5,1) so that theta=0 and theta=omega are exactly the two edges
// meeting at H) satisfies Delta u = 0 everywhere in Omega except at H
// itself, and vanishes identically on those two edges. Prescribing u's own
// trace as Dirichlet data on the REST of the boundary makes u the *exact*
// solution of -Delta u = 0 with u = u_exact on the whole of dOmega -- this
// is the standard L-shaped-domain corner-singularity benchmark used
// throughout the graded-mesh FEM literature (u in H^{1+2/3-eps}, so a
// quasi-uniform mesh caps the L2 rate at h^{4/3} regardless of polynomial
// degree -- exactly what was measured with the background-fixed corner
// refinement earlier this session; a properly graded mesh should recover
// the full h^{k+1} rate).
//
// Reuses the acoustic_one_field_assembler + linear_solver + a
// compute_errors_one_field-style loop, following the pattern of
// EllipticOneFieldConvergenceTest.hpp -- a single sparse solve per level,
// no time-stepping, no reference mesh needed since the exact solution is
// known in closed form.
//
// Usage example: ../../../wave_propagation -k3

void EllipticLshape_conv_test(int argc, char **argv);

void EllipticLshape_conv_test(int argc, char **argv)
{
    std::cout << std::endl << bold << red << "   ELLIPTIC L-SHAPE CORNER-SINGULARITY CONV TEST" << std::endl << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, true> boundary_type;

    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);

    // ---------------- corner-singular exact solution ----------------
    const RealType omega = 1.5 * M_PI;   // 270 deg reentrant angle
    const RealType nu    = M_PI / omega; // = 2/3

    auto polar = [](const typename mesh_type::point_type& pt) -> std::pair<RealType, RealType> {
        RealType dx = pt.x() - 0.5, dy = pt.y() - 0.5;
        RealType r = std::sqrt(dx * dx + dy * dy);
        RealType theta = std::atan2(dy, dx) - 0.5 * M_PI;
        if (theta < 0.0) theta += 2.0 * M_PI;
        return {r, theta};
    };

    auto exact_scal_fun = [&](const typename mesh_type::point_type& pt) -> RealType {
        auto [r, theta] = polar(pt);
        if (r < 1.0e-12) return 0.0;
        return std::pow(r, nu) * std::sin(nu * theta);
    };

    auto exact_flux_fun = [&](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        auto [r, theta] = polar(pt);
        if (r < 1.0e-10) return {0.0, 0.0};
        RealType A = nu * std::pow(r, nu - 1.0) * std::sin(nu * theta); // du/dr
        RealType B = nu * std::pow(r, nu - 1.0) * std::cos(nu * theta); // (1/r) du/dtheta
        RealType c = std::cos(theta), s = std::sin(theta);
        // grad u = A*(-sin(theta), cos(theta)) + B*(-cos(theta), -sin(theta))
        std::vector<RealType> g(2);
        g[0] = -A * s - B * c;
        g[1] =  A * c - B * s;
        return g;
    };

    auto null_rhs_fun = [](const typename mesh_type::point_type&) -> RealType { return 0.0; };
    acoustic_material_data<RealType> unit_material(1.0, 1.0);

    // ---------------- mesh sequence: same generator as the LTS test ----------------
    const std::string mesh_dir = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/lshape_graded/meshes/";
    std::vector<int> level_N = {0, 1, 2, 3, 4};
    const size_t mesh_k = sim_data.m_k_degree;

    std::ostringstream fname;
    fname << "elliptic_lshape_convergence_k_" << mesh_k << ".txt";
    std::ofstream error_file(fname.str());

    for (int N : level_N) {

        std::ostringstream mesh_fname;
        mesh_fname << mesh_dir << "lshape_graded_k" << mesh_k << "_N" << N << ".txt";

        mesh_type msh;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        mesh_builder.set_poly_mesh_file(mesh_fname.str());
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);

        RealType h_max = 0.0;
        for (auto & cell : msh) {
            RealType h_l = diameter(msh, cell);
            if (h_l > h_max) h_max = h_l;
        }

        std::cout << bold << red << "\n   MESH LEVEL N=" << N << reset << std::endl;
        std::cout << bold << cyan << "      n_cells=" << msh.cells_size() << "  h_max=" << h_max << reset << std::endl;

        boundary_type bnd(msh);
        bnd.addDirichletEverywhere(exact_scal_fun);

        auto assembler = acoustic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
        assembler.load_material_data(msh, unit_material);
        assembler.assemble(msh, null_rhs_fun);
        assembler.apply_bc(msh);

        linear_solver<RealType> analysis(assembler.LHS);
        analysis.factorize();
        Matrix<RealType, Dynamic, 1> x_dof = analysis.solve(assembler.RHS);

        error_file << "N=" << N << " h_max=" << std::setprecision(15) << h_max << "\n";
        postprocessor<mesh_type>::compute_errors_one_field(msh, hho_di, assembler, x_dof, exact_scal_fun, exact_flux_fun, error_file);
    }

    error_file.close();
    std::cout << bold << red << "\n   --> " << fname.str() << reset << std::endl << std::endl;
}
