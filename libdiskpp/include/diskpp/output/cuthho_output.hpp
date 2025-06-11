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

#include "silo.hpp"

// #include "level_set.hpp"
// #include "mesh_storage.hpp"
// #include "diskpp/quadratures/quadratures.hpp"


template<typename T>
class postprocess_output_object {

public:
    postprocess_output_object()
    {}

    virtual bool write() = 0;
};

template<typename T>
class silo_output_object : public postprocess_output_object<T>
{

};

template<typename T>
class gnuplot_output_object : public postprocess_output_object<T>
{
    std::string                                 output_filename;
    std::vector< std::pair<disk::point<T,2>, T > >   data;

public:
    gnuplot_output_object(const std::string& filename)
        : output_filename(filename)
    {}

    void add_data(const disk::point<T,2>& pt, const T& val)
    {
        data.push_back( std::make_pair(pt, val) );
    }

    bool write()
    {
        std::ofstream ofs(output_filename);

        for (auto& d : data)
            ofs << d.first.x() << " " << d.first.y() << " " << d.second << std::endl;

        ofs.close();

        return true;
    }
};


template<typename T>
class postprocess_output
{
    std::list< std::shared_ptr< postprocess_output_object<T>> >     postprocess_objects;

public:
    postprocess_output()
    {}

    void add_object( std::shared_ptr<postprocess_output_object<T>> obj )
    {
        postprocess_objects.push_back( obj );
    }

    bool write(void) const
    {
        for (auto& obj : postprocess_objects)
            obj->write();

        return true;
    }
};



template<disk::mesh_2D Mesh, typename Function>
void
output_mesh_info(Mesh& msh, const Function& level_set_function) {

    using RealType = typename Mesh::coordinate_type;

    /************** OPEN SILO DATABASE **************/
    disk::silo_database silo;
    silo.create("cuthho_meshinfo.silo");
    silo.add_mesh(msh, "mesh");

    /************** MAKE A SILO VARIABLE FOR CELL POSITIONING **************/
    std::vector<RealType> cut_cell_markers;
    std::vector<RealType> cell_indexes;
    size_t cell_ind = 0;
    for (auto& cl : cells(msh)) {
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
    for (auto& pt : points(msh))
        level_set_vals.push_back( level_set_function(pt) );
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), disk::nodal_variable_t);

    /************** MAKE A SILO VARIABLE FOR NODE POSITIONING **************/
    std::vector<RealType> node_pos;
    for (auto& n : nodes(msh))
        node_pos.push_back(locate(msh, n) == disk::location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), disk::nodal_variable_t);

    auto storage = msh.backend_storage();
    std::vector<RealType> cell_set;
    for (auto& cl : cells(msh)) {

        RealType r;

        auto cl_id = offset(msh, cl);
        auto& cl_cut = storage -> cut_cell_data[cl_id].cut;

        switch (cl_cut) {

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

    /*************  MAKE AN OUTPUT FOR THE INTERSECTION POINTS *************/
    std::vector<RealType> int_pts_x;
    std::vector<RealType> int_pts_y;

    for (auto& fc : faces(msh)) {

        auto fc_id = offset(msh, fc);
        auto& fc_loc = storage -> cut_face_data[fc_id].loc;
        auto& fc_intersection_point = storage -> cut_face_data[fc_id].intersection_point;

        if (fc_loc != disk::location::ON_INTERFACE) continue;

        RealType x = fc_intersection_point.x();
        RealType y = fc_intersection_point.y();

        int_pts_x.push_back(x);
        int_pts_y.push_back(y);

    }

    std::ofstream points_file("int_points.3D", std::ios::out | std::ios::trunc);

    if(points_file)
    {
        // instructions
        points_file << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_pts_x.size(); i++)
        {
            points_file << int_pts_x[i] << "   " <<  int_pts_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        points_file.close();
    }

    else
        std::cerr << "Points_file has not been opened" << std::endl;


    /*************  MAKE AN OUTPUT FOR THE INTERFACE *************/
    std::vector<RealType> int_x;
    std::vector<RealType> int_y;

    for (auto& cl : cells(msh)) {

        auto cl_id = offset(msh, cl);
        auto& cl_loc = storage -> cut_cell_data[cl_id].loc;
        auto& cl_interface = storage -> cut_cell_data[cl_id].interface;

        if (cl_loc != disk::location::ON_INTERFACE) continue;

        for(size_t i = 0; i < cl_interface.size(); i++) {

            RealType x = cl_interface.at(i).x();
            RealType y = cl_interface.at(i).y();

            int_x.push_back(x);
            int_y.push_back(y);
        }
    }
    std::ofstream interface_file("interface.3D", std::ios::out | std::ios::trunc);

    if(interface_file) {

        // instructions
        interface_file << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_x.size(); i++) {
            interface_file << int_x[i] << "   " <<  int_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        interface_file.close();
    }

    else
        std::cerr << "Interface_file has not been opened" << std::endl;

}




/////// test_info -> for error output
template<typename T>
class test_info {
public:
    test_info()
        {
            H1 = 0.0;
            L2 = 0.0;
            cond = 0.0;
        }
    T H1; // H1-error
    T L2; // L2-error
    T cond; // condition number
};


/////// stokes_test_info -> for stokes problem
template<typename T>
class stokes_test_info {
public:
    stokes_test_info()
        {
            H1_vel = 0.0;
            L2_vel = 0.0;
            L2_p = 0.0;
            cond = 0.0;
        }
    T H1_vel; // H1-error for velocity
    T L2_vel; // L2-error for velocity
    T L2_p;   // L2-error for pressure
    T cond;   // condition number
};

