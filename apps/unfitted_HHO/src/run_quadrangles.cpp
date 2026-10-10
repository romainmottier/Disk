/* Convergence test of the unfitted HHO method on quadrangles */
#include "unfitted_HHO_run.hpp"

void run_quadrangles(const simulation_options& opts, bool agglomeration)
{
    Elliptic_Interface<disk::cartesian_mesh<double, 2>>(opts, agglomeration);
}
