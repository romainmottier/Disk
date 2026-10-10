/* Convergence test of the unfitted HHO method on tetrahedra */
#include "unfitted_HHO_run.hpp"

void run_tetrahedra(const simulation_options& opts, bool agglomeration)
{
    run_mesh<disk::simplicial_mesh<double, 3>>(opts, agglomeration);
}
