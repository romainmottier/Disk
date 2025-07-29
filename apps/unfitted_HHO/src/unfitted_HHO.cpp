
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <list>
#include <regex>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/output/cuthho_output.hpp"
#include "diskpp/methods/hho"

using RealType = double;
typedef disk::cartesian_mesh<RealType, 2> mesh_type;
#include "diskpp/mesh/cut_mesh.hpp"
#include "common/preprocessor.hpp"
#include "common/test_cases.hpp"
#include "common/cuthho_assemblers.hpp"
#include "common/operators.hpp"
#include "common/postprocessor.hpp"

#define reset   "\033[0m"
#define red     "\033[31m"
#define bold    "\033[1m"
#define yellow  "\033[33m"
#define blue    "\033[34m"

typedef disk::BoundaryConditions<mesh_type, true> boundary_type;

// ../unfitted_HHO -p0 -l4 -i0 -k3 -n0 -r4 -c0 -s0 -t0 -f1 
void Elliptic_Poly_Ext(int argc, char **argv);
void Elliptic_Poly_Ext(int argc, char **argv) {

    timecounter tc, tck, tcl;
    tc.tic();

    using T = double;    
    
    // SIMULATION DATA & MESH PREPROCESSING
    preprocessor<mesh_type> preprocess; 
    preprocessor<mesh_type>::SimulationOptions opts;
    preprocess.simulation_data(argc, argv, opts); 
    auto level_set_function = preprocess.make_level_set_function<double>(opts.lv_set);
    mesh_type msh = preprocess.msh_generation(*level_set_function, opts.l_divs, opts.int_refsteps);
    make_polynomial_extension(msh, level_set_function);
    preprocess.msh_debug(msh, *level_set_function, opts.l_divs, opts.debug);

    // TEST CASE 
    preprocessor<mesh_type>::params<double> mat_params;
    preprocessor<mesh_type>::update_mat_prop(mat_params, 1.0, 1.0);
    auto test_case = select_test_case<mesh_type>(opts.test_case, level_set_function.get());

    // HHO DISCRTEIZATION
    disk::hho_degree_info hdi(opts.degree+1, opts.degree); // Mixed order discretization
    auto method = make_gradrec_interface_method<double, mesh_type, decltype(test_case)>(test_case);
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(test_case.bcs_fun);
    auto assembler = unfitted_elliptic_interface_assembler<mesh_type>(msh, hdi, bnd);


}


int main(int argc, char **argv) {

    DBSetDeprecateWarnings(0);
    using RealType = double;
    
    int ch, prototype = 0;
    while ((ch = getopt(argc, argv, "p:")) != -1) {
        if (ch == 'p') {
            prototype = atoi(optarg);
            break;  // stop après avoir trouvé -p
        }
    }

    // POSTPRO HHO SOLUTION
    std::string error_file_txt = "solution_error_file.txt";
    std::ofstream error_file(error_file_txt);
    postprocessor<mesh_type>::write_conv_sol(error_file_txt);

    if (prototype == 0) {
        std::cout << std::endl << "            " << bold << red;
        std::cout << "ELLIPTIC INTERFACE PROBLEM STABILIZED WITH POLYNOMIAL EXTENSION" << reset << std::endl << std::endl;
        Elliptic_Poly_Ext(argc, argv);
        return 0;
    }

}
