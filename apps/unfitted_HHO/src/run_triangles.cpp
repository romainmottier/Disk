/* Convergence test of the unfitted HHO method on triangles */
#include "unfitted_HHO_run.hpp"

void run_triangles(const simulation_options& opts, bool agglomeration)
{
    run_mesh<disk::simplicial_mesh<double, 2>>(opts, agglomeration);
}
