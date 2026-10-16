/*
 * Unfitted HHO methods stabilized by polynomial extension for elliptic
 * interface problems (E. Burman, A. Ern, R. Mottier).
 *
 * Prototypes: -p0 stabilization by polynomial extension
 *             -p1 stabilization by cell agglomeration [12]
 *             -p2 boundary value correction on a polygonal interface, cell
 *                 agglomeration (E. Burman, common/bvc.hpp)
 *             -p3 boundary value correction, stabilization by polynomial
 *                 extension (common/bvc_extension.hpp)
 * Meshes:     -M0 quadrangles (paper), 1 triangles, 2 polygons (hexagons),
 *             -M3 hexahedra, 4 tetrahedra, 5 polyhedra (hexagonal prisms)
 * Local refinement: -x n (default 3) refines at most n times the cells that
 *             the interface does not cut once, if any (polytopal meshes); -x 0
 *             keeps the meshes (paper)
 * Example: convergence test of Fig. 6 (circle, k = 0..3, l = 0..4)
 *   ./unfitted_HHO -p0 -i0 -t0 -k3 -l4 -r8
 * Boundary value correction, Table 2 of the paper (N = 8..128):
 *   ./unfitted_HHO -p2 -i0 -t0 -m4 -N8 -K0 -k3 -l4 -x0
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

    if (opts.prototype > 3) {
        std::cout << "Unknown prototype " << opts.prototype << std::endl;
        return 1;
    }
    bool agglomeration = (opts.prototype == 1);

    std::cout << std::endl << "            " << bold << red;
    if (opts.prototype == 2)
        std::cout << "ELLIPTIC INTERFACE PROBLEM, BOUNDARY VALUE CORRECTION ON A POLYGONAL INTERFACE" << reset
                  << std::endl << std::endl;
    else if (opts.prototype == 3)
        std::cout << "ELLIPTIC INTERFACE PROBLEM, BOUNDARY VALUE CORRECTION STABILIZED WITH POLYNOMIAL EXTENSION"
                  << reset << std::endl << std::endl;
    else if (agglomeration)
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
