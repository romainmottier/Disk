//
//  wave_propagation.cpp
//
//  Elasto-acoustic wave propagation with HHO methods. Two prototypes, entirely configured by a Lua input file
//  (parameters/lua_parameters.hpp; annotated templates share/first_order.lua and share/second_order.lua):
//    formulation = "first_order"  : HHOFirstOrder  (velocity-stress / velocity-pressure, four-field assembler)
//    formulation = "second_order" : HHOSecondOrder (displacement / velocity potential, two-field assembler)
//  The input file sets the whole discretization: polynomial degree, equal or mixed order, stabilization and its
//  scaling, time integrator (ERK, SDIRK, LTS, LI), time step, mesh and local refinement, manufactured solution, study.
//      ./wave_propagation input.lua ["lua statements" ...]
//  The statements after the file are executed after it, e.g. "time.scheme = 'sdirk'; time.n_steps = 20".
//  The former research prototypes are kept in wave_propagation_legacy.cpp (CMake option WAVE_PROPAGATION_LEGACY).
//

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

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
using namespace Eigen;

#include "diskpp/common/timecounter.hpp"
#include "diskpp/methods/hho"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/output/silo.hpp"

// application common sources
#include "output/display_settings.hpp"
#include "mesh/fitted_geometry_builders.hpp"
#include "solvers/linear_solver.hpp"
#include "physics/acoustic_material_data.hpp"
#include "physics/elastic_material_data.hpp"
#include "physics/scal_vec_analytic_functions.hpp"
#include "output/postprocessor.hpp"
#include "parameters/lua_parameters.hpp"

// RK schemes
#include "time_schemes/dirk_hho_scheme.hpp"
#include "time_schemes/dirk_butcher_tableau.hpp"
#include "time_schemes/erk_butcher_tableau.hpp"
#include "time_schemes/erk_coupling_hho_scheme.hpp"

// the two prototypes
#include "first_order/HHOFirstOrder.hpp"
#include "second_order/HHOSecondOrder.hpp"

int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " input.lua [\"lua statements\" ...]" << std::endl
                  << "       templates: share/first_order.lua, share/second_order.lua" << std::endl;
        return 1;
    }

    try {
        std::vector<std::string> statements(argv + 2, argv + argc);
        wave_parameters prm = load_wave_parameters(argv[1], statements);
        std::cout << bold << red << "   WAVE PROPAGATION, input " << argv[1] << reset << std::endl;
        prm.print(std::cout);
        if (prm.formulation == "first_order")
            HHOFirstOrder(prm);
        else
            HHOSecondOrder(prm);
    }
    catch (const std::exception& e) {
        std::cerr << "wave_propagation: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
