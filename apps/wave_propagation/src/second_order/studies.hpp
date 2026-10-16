//
//  studies.hpp
//
//  Second-order formulation (HHOSecondOrder.hpp): summaries and checks of a run (stabilization weight, frequencies,
//  stability limit of the explicit integrators, statistics of the splitting); the studies themselves are in
//  HHOSecondOrder().
//

#pragma once
#ifndef second_order_studies_hpp
#define second_order_studies_hpp

#include "time_schemes.hpp"

namespace hho_second_order {

    // Mesh, discretization and splitting data; omega_max (explicit integrators) and the spectral radius of the splitting
    // come from power iterations
    inline std::string problem_summary(problem& pb, bool explicit_scheme) {
        second_order_system& sys = pb.system();
        std::ostringstream os;
        os << "      h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min << "  gamma = " << pb.gamma;
        os << "  unknowns: cells " << sys.n_c << ", dynamic faces " << sys.n_d << ", algebraic faces " << sys.n_a
           << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << std::endl;
        if (explicit_scheme) {
            RealType omega = sys.max_frequency();
            pb.omega_max = omega;
            os << "      omega_max ~ " << omega << " (ERK4 stable for dt < ~" << 2.8 / omega << ")  ";
        }
        else
            os << "      ";
        if (sys.splitting_Q)
            os << "splitting spectral radius = " << sys.splitting_spectral_radius() << "  ";
        os << "[checks: |C off d|/|C| = " << sys.coupling_leak << ", |sigma* off blocks|/|sigma*| = " << sys.sigma_off_block << "]" << std::endl;
        sys.n_solves = sys.n_iterations = sys.n_failures = sys.max_iterations = 0;
        return os.str();
    }

    // Explicit integrators: warning when dt omega_max exceeds the imaginary-axis stability bound (ERK(3): sqrt(3),
    // ERK(4) and the fine substeps of LTS-RK4: 2 sqrt(2); ERK(1), ERK(2) are unstable for any dt on this axis and only
    // converge for small dt omega_max)
    inline void check_time_step(const problem& pb, const wave_parameters& prm, RealType dt) {
        RealType bound = 0.0;
        const std::string& scheme = prm.time.scheme;
        if (scheme == "erk" && prm.stages() == 3) bound = std::sqrt(3.0);
        if ((scheme == "erk" && prm.stages() == 4) || scheme == "lts") bound = 2.0 * std::sqrt(2.0);
        RealType z = dt * pb.omega_max / (scheme == "lts" ? RealType(pb.p) : 1.0);
        if (bound > 0.0 && z > bound)
            std::cout << bold << red << "      WARNING: dt omega_max = " << z << " > " << bound
                      << " (stability limit of the explicit integrator): increase time.n_steps" << reset << std::endl;
    }

    inline void print_splitting_stats(second_order_system& sys, std::ostream& os) {
        if (!sys.splitting_Q || sys.n_solves == 0) return;
        os << "      splitting: " << sys.n_solves << " face solves, " << RealType(sys.n_iterations) / sys.n_solves
           << " iterations on average, max " << sys.max_iterations << ", not converged " << sys.n_failures << std::endl;
        sys.n_solves = sys.n_iterations = sys.n_failures = sys.max_iterations = 0;
    }

}

#endif /* second_order_studies_hpp */
