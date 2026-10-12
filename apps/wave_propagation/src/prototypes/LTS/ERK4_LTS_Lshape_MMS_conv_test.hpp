//  Created by Romain Mottier
//
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
#pragma once

#include "lts_common.hpp"

// Exact solution: lts::lshape_corner_solution (its right-hand side is the
// precomputed pair of unit vectors -- the Dirichlet trace and the source are
// both a fixed spatial profile phi times a scalar time factor and
// assemble_rhs is linear in both, so the O(cells) assemble_rhs is not
// called at every one of the ~nt*(3+3*p) RK substages, which timed out
// already at N=2, p=1024). Classical two-level LTS-RK4 split by
// assemble_P(h_c = 0.75 h_max, nb_layer = 3).
void ERK4_LTS_Lshape_MMS_conv_test(int argc, char **argv)
{
    using namespace lts;
    std::cout << std::endl << bold << red << "   ERK4-LTS L-SHAPE MMS CONV TEST (reentrant corner)" << std::endl << std::endl;

    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);

    // Mesh sequence (same coupled-refinement quadtree family as
    // EllipticLshape_conv_test.hpp / ERK4_LTS_Lshape_conv_test.hpp).
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
        mesh_fname << mesh_root << "lshape_graded/meshes/lshape_graded_k" << mesh_k << "_N" << N << ".txt";
        acoustic_problem<lshape_corner_solution> pb(mesh_fname.str());
        const RealType h_max = pb.h_max, h_min = pb.h_min;
        const int p = pb.p_global;
        RealType h_c = (p == 1) ? 1.25 * h_max : 0.75 * h_max;
        std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p=" << p << reset << std::endl;

        // Pure acoustic domain, p(0)=0, v(0)=grad(phi).
        pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
        assembler_type& assembler = pb.assembler();
        assembler.assemble_P(pb.msh, h_c, 3);

        RealType ti = 0.0, tf = 0.05;
        RealType dt = 0.002 * h_max;
        size_t nt = macro_steps(tf, dt);
        for (size_t it = 1; it <= nt; it++)
            two_level_reference_step(pb, assembler.Pcoarse, assembler.Pfine, p, pb.x_dof, ti + dt * (it - 1), dt, true);

        // Error at t=tf vs the closed-form solution.
        RealType l2_error = pb.pressure_error(pb.x_dof, tf);
        std::cout << bold << cyan << "      L2_error(pressure, t=" << tf << ") = "
                  << std::setprecision(10) << l2_error << reset << std::endl;
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << l2_error << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}
