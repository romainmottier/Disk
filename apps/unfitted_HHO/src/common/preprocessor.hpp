#pragma once
#ifndef preprocessor_hpp
#define preprocessor_hpp

#include <iostream>
#include <memory>
#include <unistd.h>  // for getopt
#include <stdexcept>
#include <vector>
#include <array>
#include "diskpp/mesh/mesh.hpp"

#define reset   "\033[0m"
#define red     "\033[31m"
#define bold    "\033[1m"
#define yellow  "\033[33m"
#define blue    "\033[34m"

template<typename Mesh>
class preprocessor {

public:

    struct SimulationOptions {
        size_t prototype = 0;
        size_t l_divs = 0;
        size_t lv_set = 0;
        size_t degree = 0;
        size_t nt_divs = 0;
        size_t int_refsteps = 0;
        size_t sc = 0;
        size_t iterative_solver = 0;
        size_t test_case = 0;
        size_t debug = 0;
    };

    void simulation_data(int argc, char** argv, SimulationOptions& opts) const {
        int ch;
        while ((ch = getopt(argc, argv, "p:l:i:k:n:r:c:s:t:f:")) != -1) {
            switch (ch) {
                case 'p': opts.prototype = std::atoi(optarg); break;
                case 'l': opts.l_divs = std::atoi(optarg); break;
                case 'i': opts.lv_set = std::atoi(optarg); break;
                case 'k': opts.degree = std::atoi(optarg); break;
                case 'n': opts.nt_divs = std::atoi(optarg); break;
                case 'r': opts.int_refsteps = std::atoi(optarg); break;
                case 'c': opts.sc = std::atoi(optarg); break;
                case 's': opts.iterative_solver = std::atoi(optarg); break;
                case 't': opts.test_case = std::atoi(optarg); break;
                case 'f': opts.debug = std::atoi(optarg); break;
            }
        }

        std::cout << "    " << bold << red << "SIMULATION DATA:" << reset << blue << bold << std::endl;
        std::cout << "        Prototype                  -p : " << opts.prototype << std::endl;
        std::cout << "        Space refinement level     -l : " << opts.l_divs << std::endl;
        std::cout << "        Polynomial degree          -k : " << opts.degree << "     (Face unknowns)" << std::endl;
        std::cout << "        Time refinement level      -n : " << opts.nt_divs << std::endl;
        std::cout << "        Interface refinement level -r : " << opts.int_refsteps << std::endl;
        std::cout << "        Static condensation        -c : " << opts.sc << std::endl;
        std::cout << "        Iterative solver           -s : " << opts.iterative_solver << std::endl;
        std::cout << "        Debug & Silo files         -f : " << opts.debug << std::endl;

    }

    template<typename RealType>
    std::unique_ptr<level_set<RealType>> make_level_set_function(size_t lv_set) const {

        RealType radius = 1.0/3.0;

        switch (lv_set) {
            case 0:
                std::cout << "        Level set function         -i : Circle with R = " << radius << "\n";
                return std::make_unique<circle_level_set<RealType>>(radius, 0.5, 0.5);
            case 1:
                std::cout << "        Level set function         -i : Flower\n";
                return std::make_unique<flower_level_set<RealType>>(radius, 0.5, 0.5, 8, 0.03);
            default:
                throw std::runtime_error("Unsupported level set type: " + std::to_string(lv_set));
        }
    }

    inline void run_cmd(const std::string& cmd) {
        int ret = std::system(cmd.c_str());
    }

    mesh_type 
    msh_generation(level_set<RealType> & level_set_function, size_t l_divs, size_t int_refsteps){
        
        mesh_type msh;
        auto mesher = disk::make_simple_mesher(msh);
        for (auto nr = 0; nr < l_divs; nr++)
            mesher.refine();
        
        detect_node_position(msh, level_set_function); 
        detect_cut_faces(msh, level_set_function); 
        detect_cut_cells(msh, level_set_function);
        detect_cut_type(msh, level_set_function);
        make_neighbors_info_cartesian(msh);
        refine_interface(msh, level_set_function, int_refsteps);
        
        return msh;
        
    }
    
    void 
    msh_debug(mesh_type msh, level_set<RealType> & level_set_function, size_t l_divs, size_t debug){
        
        if (debug) {
            // print_polynomial_extension(msh);
            output_mesh_info(msh, level_set_function);
        }

        std::string mesh_info = "cuthho_meshinfo_l" + std::to_string(l_divs) + ".silo";
        std::string command = "mv cuthho_meshinfo.silo " + mesh_info;
        std::string interface = "interface_l" + std::to_string(l_divs) + ".3D";
        std::string command2 = "mv interface.3D " + interface;
        std::string pairing1 = "agglo_five_l" + std::to_string(l_divs) + ".okc";
        std::string command3 = "mv agglo_five.okc " + pairing1;
        std::string pairing2 = "agglo_four_l" + std::to_string(l_divs) + ".okc";
        std::string command4 = "mv agglo_four.okc " + pairing2;
        run_cmd(command);
        run_cmd(command2);
        run_cmd(command3);
        run_cmd(command4);
        
        return;
        
    }
    
    template<typename T>
    struct params {
        T kappa_1, kappa_2;
        T c_1, c_2;
        params() : kappa_1(1.0), kappa_2(1.0), c_1(1.0), c_2(1.0) {}
        params(T kap1, T kap2) : kappa_1(kap1), kappa_2(kap2), c_1(1.0), c_2(1.0) {}
    };
    
    template<typename T>
    static void update_mat_prop(params<T>& parms, T new_kappa1, T new_kappa2) {
        parms.kappa_1 = new_kappa1;
        parms.kappa_2 = new_kappa2;
        std::cout << "        Contrast                      : " << parms.kappa_2/parms.kappa_1 << std::endl;
    }

};

#endif // preprocessor_hpp
