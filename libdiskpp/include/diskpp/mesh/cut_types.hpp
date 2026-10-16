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

#include <array>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "point.hpp"

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

/* Representation of the interface inside a cut polygon:
 *  - polyline:   the points of cut_polygon::polyline (2^r segments, paper);
 *  - polynomial: x(t) = p0 + t (p1 - p0) + h(t) nu, h interpolating the heights
 *                of the interface above its chord at the nodes;
 *  - implicit:   the same graph above the chord, the height being the root of
 *                the level set along nu (exact geometry). */
enum class interface_kind {
    polyline,
    polynomial,
    implicit
};

/* A polygon cut once by the interface: a cell of a 2D mesh, or a face of a 3D
 * mesh in the coordinates of its plane. */
template<typename T>
struct cut_polygon {
    std::vector<point<T,2>> vertices;    // counterclockwise
    std::vector<location>   vertex_loc;  // IN_NEGATIVE_SIDE or IN_POSITIVE_SIDE
    std::vector<size_t>     vertex_ids;  // point ids of the vertices in the mesh
    point<T,2>              p0, p1;      // ends of the interface, the negative side on the left of p0 -> p1
    size_t                  edge0 = 0, edge1 = 0; // edges (i, i+1) containing p0 and p1
    bool                    convex = true;
    T                       h = 0;       // diameter

    interface_kind          kind = interface_kind::polyline;
    std::vector<point<T,2>> polyline;    // points of the interface, from p0 to p1
    std::vector<T>          nodes, heights; // polynomial interface
    // points with respect to which the negative and positive parts are
    // star-shaped (centers of the sectors of their quadratures)
    std::array<point<T,2>, 2> centers;

    // implicit interface of a polygon whose parts are not star-shaped: the
    // polygon is split into cut polygons (ordered along the interface, from
    // p0 to p1) and polygons lying on one side
    std::vector<cut_polygon<T>>                                 children;
    std::vector<std::pair<std::vector<point<T,2>>, location>>   uncut_children;
};

/* Planar face of a polyhedron cut by the interface */
template<typename T>
struct cut_polyhedron_face {
    std::vector<point<T,3>> vertices;    // counterclockwise seen from outside
    Eigen::Matrix<T,3,1>    normal;      // outward unit normal
    location                loc = location::UNDEF; // side, or ON_INTERFACE if cut
    point<T,3>              origin;      // frame x = origin + U u of its plane (cut faces)
    Eigen::Matrix<T,3,2>    U;
    cut_polygon<T>          poly;        // cut in the coordinates of the plane (cut faces)
};

/* A convex polyhedron (3D cell or part of a cell) cut once by the interface.
 * The interface inside it is the graph, in the direction nu, of the level set
 * above the region D of the plane (plane_origin, a1, a2) enclosed by the
 * projection of the loop formed by the interface curves of its cut faces.
 * A polyhedron for which D or its parts are not star-shaped is split by a
 * plane into the polyhedra "children"; "linear" polyhedra are integrated by
 * subdivision (fallback). */
template<typename T>
struct cut_polyhedron {
    std::vector<cut_polyhedron_face<T>>  faces;
    location                loc = location::UNDEF; // side, or ON_INTERFACE if cut
    T                       h = 0;       // diameter
    bool                    linear = false;
    bool                    cut_more_than_once = false; // several loops: the mesh does not resolve the interface

    Eigen::Matrix<T,3,1>    nu, a1, a2;  // a1 x a2 = nu, nu points towards the positive side
    point<T,3>              plane_origin;
    // cut faces (index in faces, and whether their interface curve is
    // traversed from p1 to p0) in the counterclockwise order around nu
    std::vector<std::pair<size_t, bool>> loop;
    point<T,2>              D_center;    // center of the sectors of D
    T                       center_height = 0; // height of the interface above D_center
    std::array<point<T,3>, 2> centers;   // centers of the cones of the two parts
    std::vector<point<T,3>> interface;   // sample points of the interface

    // flat interface (cut_polyhedron_options::flat, for the boundary value
    // correction): its triangles, vertices on the level set, oriented from the
    // negative to the positive side (see cutgeo::build_flat_cut_polyhedron)
    bool                    flat = false;
    std::vector<std::array<point<T,3>, 3>> flat_interface;

    std::vector<cut_polyhedron<T>> children;
};

/* Neighbors and pairing of a cell (polynomial extension and agglomeration) */
struct cut_cell_pairing {
    std::set<size_t>        f_neighbors; // cells sharing a face
    std::set<size_t>        d_neighbors; // other cells sharing a vertex

    // polynomial extension: cell paired with an ill-cut cell, and ill-cut
    // cells paired with a cell
    std::optional<size_t>   stabilizing_cell;
    std::vector<size_t>     dependent_ill_cut_neg_cells;
    std::vector<size_t>     dependent_ill_cut_pos_cells;
};

// CUT CELLS INFOS
template<typename T, size_t DIM>
struct cut_cell_info;

template<typename T>
struct cut_cell_info<T, 2> : public cut_cell_pairing {
    location                loc = location::UNDEF;
    cut_type                cut = cut_type::UNDEF;
    cut_polygon<T>          poly;        // geometry of the cut (cut cells)
};

template<typename T>
struct cut_cell_info<T, 3> : public cut_cell_pairing {
    location                loc = location::UNDEF;
    cut_type                cut = cut_type::UNDEF;
    cut_polyhedron<T>       poly;        // geometry of the cut (cut cells)
};

// CUT FACES INFOS
template<typename T, size_t DIM>
struct cut_face_info;

template<typename T>
struct cut_face_info<T, 2> {
    location                loc = location::UNDEF;
    point<T,2>              intersection_point;
};

/* A face of a 3D mesh: the cut is described in the orthonormal frame
 * (origin, U) of its plane */
template<typename T>
struct cut_face_info<T, 3> {
    location                loc = location::UNDEF;
    point<T,3>              origin;
    Eigen::Matrix<T,3,2>    U;
    cut_polygon<T>          poly;
};

// CUT NODES INFOS
template<typename T, size_t DIM>
struct cut_node_info {
    location loc = location::UNDEF;
};

} // namespace disk
