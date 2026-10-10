/* Convergence test of the unfitted HHO method on hexahedra */
#include "unfitted_HHO_run.hpp"

void run_hexahedra(const simulation_options& opts, bool agglomeration)
{
    run_mesh<disk::cartesian_mesh<double, 3>>(opts, agglomeration);
}
