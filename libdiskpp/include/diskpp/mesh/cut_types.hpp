/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include <vector>
#include <optional>

#include "point.hpp"
#include "diskpp/quadratures/quadrature_point.hpp"

namespace disk {

enum class location {
    IN_NEGATIVE_SIDE,
    IN_POSITIVE_SIDE,
    ON_INTERFACE,
    UNDEF
};

enum class cut_type {
    T_OK,
    T_KO_NEG,
    T_KO_POS,
    UNDEF
};


// CUT CELLS INFOS 
template<typename T, size_t DIM>
struct subcell;

template<typename T>
struct subcell<T, 2>
{};

template<typename T, size_t DIM>
struct cut_cell_info;

template<typename T>
struct cut_cell_info<T, 2> {
    subcell<T, 2>           neg;
    subcell<T, 2>           pos;
    std::optional<size_t>   stabilizing_cell;
    std::vector<size_t>     dependent_ill_cut_neg_cells;
    std::vector<size_t>     dependent_ill_cut_pos_cells;
    std::vector<point<T,2>> interface;
    location                loc;
    cut_type                cut;
    point<T,2>              p0, p1;
    std::vector<disk::quadrature_point<T, 2>> integration_n; // composite integration rules
    std::vector<disk::quadrature_point<T, 2>> integration_p;
    bool                        distorted;
    std::set<size_t>            f_neighbors; // face neighbors
    std::set<size_t>            d_neighbors; // diagonal neighbors
    size_t local_dofs;

    cut_cell_info() :
        loc(location::UNDEF),
        cut(cut_type::UNDEF),
        distorted(false)
    {}
};

// CUT FACES INFOS 
template<typename T, size_t DIM>
struct subface;

template<typename T>
struct subface<T, 2>
{};

template<typename T, size_t DIM>
struct cut_face_info;

template<typename T>
struct cut_face_info<T, 2> {
    subface<T, 2>           neg;
    subface<T, 2>           pos;
    point<T,2>              intersection_point;
    std::vector<point<T,2>> interface_nodes;
    location                loc;
    size_t                  node_inside; 
};

// CUT NODES INFOS 
template<typename T, size_t DIM>
struct cut_node_info;

template<typename T>
struct cut_node_info<T, 2> {
    location loc;
};

} // namespace disk
