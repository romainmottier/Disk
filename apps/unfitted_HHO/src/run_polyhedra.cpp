/* Convergence test of the unfitted HHO method on polyhedra */
#include "unfitted_HHO_run.hpp"

void run_polyhedra(const simulation_options& opts, bool agglomeration)
{
    Elliptic_Interface<disk::generic_mesh<double, 3>>(opts, agglomeration);
}
