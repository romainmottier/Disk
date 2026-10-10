/* Convergence test of the unfitted HHO method on polygons */
#include "unfitted_HHO_run.hpp"

void run_polygons(const simulation_options& opts, bool agglomeration)
{
    Elliptic_Interface<disk::generic_mesh<double, 2>>(opts, agglomeration);
}
