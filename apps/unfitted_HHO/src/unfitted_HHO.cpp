/*
 * Unfitted HHO methods stabilized by polynomial extension for elliptic
 * interface problems (E. Burman, A. Ern, R. Mottier).
 *
 * Prototypes: -p0 stabilization by polynomial extension
 *             -p1 stabilization by cell agglomeration [12]
 * Meshes:     -M0 quadrangles (paper), 1 triangles, 2 polygons (hexagons),
 *             -M3 hexahedra, 4 tetrahedra, 5 polyhedra (hexagonal prisms)
 * Local refinement: -x n (default 3) refines at most n times the cells that
 *             the interface does not cut once, if any (polytopal meshes); -x 0
 *             keeps the meshes (paper)
 * Example: convergence test of Fig. 6 (circle, k = 0..3, l = 0..4)
 *   ./unfitted_HHO -p0 -i0 -t0 -k3 -l4 -r8
 */

#include <iostream>

#include "unfitted_HHO_run.hpp"

int main(int argc, char **argv)
{
    DBSetDeprecateWarnings(0);

    /* all the options are parsed at once: a first pass looking only for -p
     * would make getopt permute the values given as separate arguments */
    simulation_options opts;
    simulation_data(argc, argv, opts);

    if (opts.prototype > 1) {
        std::cout << "Unknown prototype " << opts.prototype << std::endl;
        return 1;
    }
    bool agglomeration = (opts.prototype == 1);

    std::cout << std::endl << "            " << bold << red;
    if (agglomeration)
        std::cout << "ELLIPTIC INTERFACE PROBLEM STABILIZED WITH CELL AGGLOMERATION" << reset << std::endl << std::endl;
    else
        std::cout << "ELLIPTIC INTERFACE PROBLEM STABILIZED WITH POLYNOMIAL EXTENSION" << reset << std::endl << std::endl;

    switch (mesh_kind(opts.mesh)) {
        case mesh_kind::quadrangles: run_quadrangles(opts, agglomeration); break;
        case mesh_kind::triangles:   run_triangles(opts,   agglomeration); break;
        case mesh_kind::polygons:    run_polygons(opts,    agglomeration); break;
        case mesh_kind::hexahedra:   run_hexahedra(opts,   agglomeration); break;
        case mesh_kind::tetrahedra:  run_tetrahedra(opts,  agglomeration); break;
        case mesh_kind::polyhedra:   run_polyhedra(opts,   agglomeration); break;
        default:
            std::cout << "Unknown mesh " << opts.mesh << std::endl;
            return 1;
    }

    return 0;
}
