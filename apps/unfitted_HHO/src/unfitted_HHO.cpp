
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

#include "diskpp/mesh/cut_mesh.hpp"
#include "common/cutmesh_generation.hpp"

using RealType = double;

void Elliptic_Poly_Ext(int argc, char **argv);
void Elliptic_Poly_Ext(int argc, char **argv) {
    
    timecounter tc, tck, tcl;
    tc.tic();

    using T = double;    

    // ################################################## Mesh generation 
    // ########## Level set function
    RealType radius = 1.0/3.0;  
    auto level_set_function = circle_level_set<RealType>(radius, 0.5, 0.5);          
    // auto level_set_function = flower_level_set<RealType>(radius, 0.5, 0.5, 8, 0.03);  
    // auto level_set_function = flower_level_set<RealType>(radius, 0.5, 0.5, 6, 0.045);  
    
    mesh_type msh = MeshGeneration(level_set_function);
    
}


int main(int argc, char **argv) {
    Elliptic_Poly_Ext(argc, argv);
    return 0;
}
