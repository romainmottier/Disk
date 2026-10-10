/*
 *       /\        Matteo Cicuttin (C) 2017,2018; Guillaume Delay 2018,2019
 *      /__\       matteo.cicuttin@enpc.fr        guillaume.delay@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is ProtoN, a library for fast Prototyping of
 *  /_\/_\/_\/_\   Numerical methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <fstream>
#include <vector>

#include "silo.hpp"
#include "diskpp/mesh/cut_mesh.hpp"

/* Silo file (cut cells, cut types, level set, nodes) and text files
 * (intersection points and interface: "X Y Z val") describing a cut mesh */
template<disk::CutMesh Mesh, typename Function>
void
output_mesh_info(Mesh& msh, const Function& level_set_function) {

    using RealType = typename Mesh::coordinate_type;

    /************** OPEN SILO DATABASE **************/
    disk::silo_database silo;
    silo.create("cuthho_meshinfo.silo");
    auto background = msh.background();
    silo.add_mesh(background, "mesh");

    /************** MAKE A SILO VARIABLE FOR CELL POSITIONING **************/
    std::vector<RealType> cut_cell_markers;
    std::vector<RealType> cell_indexes;
    size_t cell_ind = 0;
    for (auto& cl : msh) {
        if (locate(msh, cl) == disk::location::IN_POSITIVE_SIDE)
            cut_cell_markers.push_back(1.0);
        else if (locate(msh, cl) == disk::location::IN_NEGATIVE_SIDE)
            cut_cell_markers.push_back(-1.0);
        else if (locate(msh, cl) == disk::location::ON_INTERFACE)
            cut_cell_markers.push_back(0.0);
        else
            throw std::logic_error("shouldn't have arrived here...");
        cell_indexes.push_back(cell_ind);
        cell_ind++;
    }
    silo.add_variable("mesh", "cut_cells", cut_cell_markers.data(), cut_cell_markers.size(), disk::zonal_variable_t);
    silo.add_variable("mesh", "cell_index", cell_indexes.data(), cell_indexes.size(), disk::zonal_variable_t);

    /************** MAKE A SILO VARIABLE FOR LEVEL SET FUNCTION **************/
    std::vector<RealType> level_set_vals;
    for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        level_set_vals.push_back( level_set_function(*itor) );
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), disk::nodal_variable_t);

    /************** MAKE A SILO VARIABLE FOR NODE POSITIONING **************/
    std::vector<RealType> node_pos;
    for (auto& nd_data : msh.cut_node_data)
        node_pos.push_back(nd_data.loc == disk::location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), disk::nodal_variable_t);

    std::vector<RealType> cell_set;
    for (auto& cl : msh) {

        RealType r = 0.0;

        switch (cut_info(msh, cl).cut) {

            case disk::cut_type::UNDEF:
                r = 0.0;
                break;

            case disk::cut_type::T_OK:
                r = 1.0;
                break;

            case disk::cut_type::T_KO_NEG:
                r = 2.0;
                break;

            case disk::cut_type::T_KO_POS:
                r = 3.0;
                break;

        }

        cell_set.push_back( r );
    }
    silo.add_variable("mesh", "cut_type", cell_set.data(), cell_set.size(), disk::zonal_variable_t);

    silo.close();

    auto write_point = [](std::ofstream& ofs, const typename Mesh::point_type& p) {
        ofs << p[0] << "   " << p[1] << "   " << ((Mesh::dimension == 3) ? p[2] : 0.0) << "     0.0" << std::endl;
    };

    /*************  MAKE AN OUTPUT FOR THE INTERSECTION POINTS *************/
    std::ofstream points_file("int_points.3D", std::ios::out | std::ios::trunc);
    if (points_file) {
        points_file << "X   Y   Z   val" << std::endl;
        if constexpr (Mesh::dimension == 2) {
            for (auto& fc : faces(msh)) {
                const auto& fd = cut_info(msh, fc);
                if (fd.loc == disk::location::ON_INTERFACE)
                    write_point(points_file, fd.intersection_point);
            }
        }
        else {
            for (auto& ei : msh.edge_intersections)
                write_point(points_file, ei.second);
        }
        points_file.close();
    }
    else
        std::cerr << "Points_file has not been opened" << std::endl;

    /*************  MAKE AN OUTPUT FOR THE INTERFACE *************/
    std::ofstream interface_file("interface.3D", std::ios::out | std::ios::trunc);
    if (interface_file) {
        interface_file << "X   Y   Z   val" << std::endl;
        for (auto& cl : msh) {
            if (!is_cut(msh, cl))
                continue;
            if constexpr (Mesh::dimension == 2) {
                for (auto& p : cut_info(msh, cl).poly.polyline)
                    write_point(interface_file, p);
            }
            else {
                std::vector<typename Mesh::point_type> itf;
                disk::cutgeo::polyhedron_interface_samples(cut_info(msh, cl).poly, *msh.level_set_function, itf);
                for (auto& p : itf)
                    write_point(interface_file, p);
            }
        }
        interface_file.close();
    }
    else
        std::cerr << "Interface_file has not been opened" << std::endl;
}
