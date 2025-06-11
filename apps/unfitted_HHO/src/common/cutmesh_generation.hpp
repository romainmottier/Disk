
#ifndef CutMesh_hpp
#define CutMesh_hpp

using RealType = double;
typedef disk::cartesian_mesh<RealType, 2> mesh_type;

mesh_type 
MeshGeneration(level_set<RealType> & level_set_function){
    
    mesh_type msh;
    auto mesher = disk::make_simple_mesher(msh);
    for (auto nr = 0; nr < 2; nr++)
        mesher.refine();
    
    detect_node_position(msh, level_set_function); 
    detect_cut_faces(msh, level_set_function); 
    detect_cut_cells(msh, level_set_function);
    detect_cut_type(msh, level_set_function);
    make_neighbors_info_cartesian(msh);
    refine_interface(msh, level_set_function, 4);
    // make_polynomial_extension(msh, level_set_function);
    
    bool dump_debug = true;
    if (dump_debug) {
        output_mesh_info(msh, level_set_function);
        // std::string mesh_info = "cuthho_meshinfo_l" + std::to_string(l) + ".silo";
        // std::string command = "mv cuthho_meshinfo.silo " + mesh_info;
        // std::string interface = "interface_l" + std::to_string(l) + ".3D";
        // std::string command2 = "mv interface.3D " + interface;
        // std::system(command.c_str());
        // std::system(command2.c_str());
    }
    
    return msh;

}


#endif
