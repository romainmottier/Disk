//
//  fitted_geometry_builder.hpp
//
//  Created by Omar Durán
//  Contributor: Romain Mottier


#pragma once
#ifndef fitted_geometry_builder_hpp
#define fitted_geometry_builder_hpp

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cassert>
#include <thread>
#include <set>
#include "diskpp/geometry/geometry.hpp"
#include <unordered_map>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

template<typename MESH>
class fitted_geometry_builder {
    
protected:
    
    std::string m_log_file = "mesh_log.txt";
    
    size_t m_dimension = 0;
    
    size_t m_n_elements = 0;
    
public:
    
    fitted_geometry_builder() {
        
    }
    
    // set the log file name
    void set_log_file(std::string log_file) {
        m_log_file = log_file;
    }
    
    // get the geometry dimension
    size_t get_dimension() {
        return m_dimension;
    }
    
    // get the number of elements
    size_t get_n_elements() {
        return m_n_elements;
    }
    
    // build the mesh
    virtual bool build_mesh()  = 0;
    
    // move generated mesh data to an external mesh storage
    virtual void move_to_mesh_storage(MESH& msh) = 0;
    
    // Print in log file relevant mesh information
    virtual void print_log_file()    = 0;
    
    virtual ~ fitted_geometry_builder() {
        
    }
    
};

template<typename T>
class cartesian_2d_mesh_builder : public fitted_geometry_builder<disk::generic_mesh<T, 2>>
{
    typedef disk::generic_mesh<T,2>                 mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    
    struct polygon_2d
    {
        std::vector<size_t>                 m_member_nodes;
        std::set<std::array<size_t, 2>>     m_member_edges;
        int                                 m_material          = 0;
        bool                                m_elastic_material  = false;
        int                                 m_refinement_level  = 0;
        bool                                m_to_refine        = false;

        bool operator<(const polygon_2d & other) {
            return m_member_nodes < other.m_member_nodes;
        }
    };
    
    std::vector<point_type>                         points;
    std::vector<node_type>                          vertices;
    std::vector<std::array<size_t, 2>>              facets;
    std::vector<std::array<size_t, 2>>              skeleton_edges;
    std::vector<std::array<size_t, 2>>              boundary_edges;
    std::vector<polygon_2d>                         polygons;
    

    T m_lx = 0.0;
    T m_ly = 0.0;
            
    T m_x_t = 0.0;
    T m_y_t = 0.0;
    
    size_t m_nx = 0;
    size_t m_ny = 0;
    
    void reserve_storage(){
        size_t n_points = (m_nx + 1) * (m_ny + 1);
        points.reserve(n_points);
        vertices.reserve(n_points);
        
        size_t n_edges = 2*m_nx*m_ny + m_nx + m_ny;
        size_t n_skel_edges = 2*m_nx*m_ny + m_nx + m_ny;
        size_t n_bc_edges = n_edges - n_skel_edges;
        skeleton_edges.reserve(n_skel_edges);
        boundary_edges.reserve(n_bc_edges);
        
        size_t n_polygons = m_nx * m_ny;
        polygons.reserve(n_polygons);
    }
            
    void validate_edge(std::array<size_t, 2> & edge){
        assert(edge[0] != edge[1]);
        if (edge[0] > edge[1]){
            std::swap(edge[0], edge[1]);
        }
    }
    
public:

    cartesian_2d_mesh_builder(T lx, T ly, size_t nx, size_t ny) : fitted_geometry_builder<mesh_type>()
    {
        m_lx = lx;
        m_ly = ly;
        m_nx = nx;
        m_ny = ny;
        fitted_geometry_builder<mesh_type>::m_dimension = 2;
    }
    
    // uniform refinement x-direction
    void refine_mesh_x_direction(size_t n_refinements){
        for (unsigned int i = 0; i < n_refinements; i++) {
            m_nx *= 2;
        }
    }
    
    // uniform refinement y-direction
    void refine_mesh_y_direction(size_t n_refinements){
        for (unsigned int i = 0; i < n_refinements; i++) {
            m_ny *= 2;
        }
    }
    
    // uniform refinement
    void refine_mesh(size_t n_refinements){
        refine_mesh_x_direction(n_refinements);
        refine_mesh_y_direction(n_refinements);
    }
    
    // build the mesh
    bool build_mesh(){
        
        reserve_storage();
        
        std::vector<T> range_x(m_nx+1,0.0);
        T dx = m_lx/T(m_nx);
        for (unsigned int i = 0; i <= m_nx; i++) {
            range_x[i] = i*dx;
        }
        
        std::vector<T> range_y(m_ny+1,0.0);
        T dy = m_ly/T(m_ny);
        for (unsigned int i = 0; i <= m_ny; i++) {
            range_y[i] = i*dy;
        }
        
        size_t node_id = 0;
        for (unsigned int j = 0; j <= m_ny; j++) {
            T yv = range_y[j] + m_y_t;
            for (unsigned int i = 0; i <= m_nx; i++) {
                T xv = range_x[i] + m_x_t;
                point_type point(xv, yv);
                points.push_back(point);
                vertices.push_back(node_type(disk::point_identifier<2>(node_id)));
                node_id++;
            }
        }
        
        size_t edge_id = 0;
        for (size_t j = 0; j < m_ny; j++) {
            for (size_t i = 0; i < m_nx; i++) {
                
                size_t id_0 = i + j * (m_nx + 1);
                size_t id_1 = id_0 + 1;
                size_t id_2 = i + (m_nx + 1) + 1 + j* (m_nx + 1);
                size_t id_3 = id_2 - 1;
                
                // Adding edges: Cases to avoid edge duplicity
                if (i==0 && j==0) {
                    
                    std::array<size_t, 2> e0 = {id_0,id_1};
                    validate_edge(e0);
                    facets.push_back( e0 );
                    boundary_edges.push_back( e0 );
                    edge_id++;
                    
                    std::array<size_t, 2> e1 = {id_1,id_2};
                    validate_edge(e1);
                    facets.push_back( e1 );
                    if(j == m_ny - 1) boundary_edges.push_back( e1 );
                    edge_id++;
                    
                    std::array<size_t, 2> e2 = {id_2,id_3};
                    validate_edge(e2);
                    facets.push_back( e2 );
                    edge_id++;
                    
                    std::array<size_t, 2> e3 = {id_3,id_0};
                    validate_edge(e3);
                    facets.push_back( e3 );
                    boundary_edges.push_back( e3 );
                    edge_id++;
                }
                
                if ((i>0 && i < m_nx) && j==0) {
                    std::array<size_t, 2> e0 = {id_0,id_1};
                    validate_edge(e0);
                    facets.push_back( e0 );
                    boundary_edges.push_back( e0 );
                    edge_id++;
                    
                    std::array<size_t, 2> e1 = {id_1,id_2};
                    validate_edge(e1);
                    facets.push_back( e1 );
                    if(i == m_nx - 1) boundary_edges.push_back( e1 );
                    edge_id++;
                    
                    std::array<size_t, 2> e2 = {id_2,id_3};
                    validate_edge(e2);
                    facets.push_back( e2 );
                    if(j == m_ny - 1) boundary_edges.push_back( e2 );
                    edge_id++;
                }
                if (i==0 && j>0) {
                    
                    std::array<size_t, 2> e1 = {id_1,id_2};
                    validate_edge(e1);
                    facets.push_back( e1 );
                    edge_id++;
                    
                    std::array<size_t, 2> e2 = {id_2,id_3};
                    validate_edge(e2);
                    facets.push_back( e2 );
                    if(j == m_ny - 1) boundary_edges.push_back( e2 );
                    edge_id++;
                    
                    std::array<size_t, 2> e3 = {id_3,id_0};
                    validate_edge(e3);
                    boundary_edges.push_back( e3 );
                    facets.push_back( e3 );
                    edge_id++;
                }
                
                if (i>0 && j>0) {
                    
                    std::array<size_t, 2> e1 = {id_1,id_2};
                    validate_edge(e1);
                    facets.push_back( e1 );
                    if(i == m_nx - 1) boundary_edges.push_back( e1 );
                    edge_id++;
                    
                    std::array<size_t, 2> e2 = {id_2,id_3};
                    validate_edge(e2);
                    facets.push_back( e2 );
                    if(j == m_ny - 1) boundary_edges.push_back( e2 );
                    edge_id++;
                    
                }
            }
        }
        
        size_t surface_id = 0;
        for (size_t j = 0; j < m_ny; j++) {
            for (size_t i = 0; i < m_nx; i++) {
                
                size_t id_0 = i + j * (m_nx + 1);
                size_t id_1 = id_0 + 1;
                size_t id_2 = i + (m_nx + 1) + 1 + j* (m_nx + 1);
                size_t id_3 = id_2 - 1;
                
                polygon_2d polygon;
                polygon.m_member_nodes = {id_0,id_1,id_2,id_3};
                std::array<size_t, 2> e0 = {id_0,id_1};
                validate_edge(e0);
                std::array<size_t, 2> e1 = {id_1,id_2};
                validate_edge(e1);
                std::array<size_t, 2> e2 = {id_2,id_3};
                validate_edge(e2);
                std::array<size_t, 2> e3 = {id_3,id_0};
                validate_edge(e3);
                
                polygon.m_member_edges = {e0,e1,e2,e3};
                polygons.push_back( polygon );
                surface_id++;
                
            }
        }
        
        // std::cout << bold << red << std::endl << std::endl;
        // std::cout << "Debug mesh" << std::endl;
        // std::cout << reset << "Points: " << vertices.size() << std::endl;
        // std::cout << reset << "Polygons: " << polygons.size() << std::endl;
        // std::cout << reset << "Boundaries: " << boundary_edges.size() << std::endl;
        // std::cout << reset << "Facets: " << facets.size() << std::endl;
        return true;
    }
    
    void remove_duplicate_points() {
        std::unordered_map<size_t, size_t> point_mapping;  // Map old point indices to new point indices
        std::vector<point_type> unique_points;  // Unique points
        
        for (size_t i = 0; i < points.size(); ++i) {
            if (point_mapping.find(i) == point_mapping.end()) {
                // This is a new unique point
                point_mapping[i] = unique_points.size();
                unique_points.push_back(points[i]);
            }
        }
        
        // Update cell and face data with unique point indices
        for (auto& polygon : polygons) {
            for (size_t i = 0; i < polygon.m_member_nodes.size(); ++i) {
                size_t old_point_index = polygon.m_member_nodes[i];
                polygon.m_member_nodes[i] = point_mapping[old_point_index];
            }
            
            std::set<std::array<size_t, 2>> updated_edges;
            for (const auto& edge : polygon.m_member_edges) {
                size_t old_point1 = edge[0];
                size_t old_point2 = edge[1];
                size_t new_point1 = point_mapping[old_point1];
                size_t new_point2 = point_mapping[old_point2];
                updated_edges.insert({new_point1, new_point2});
            }
            polygon.m_member_edges = updated_edges;
        }
        
        // Update the internal mesh_reader data with the unique points
        points = unique_points;
        
        // Update the vertices
        vertices.clear();
        for (size_t i = 0; i < unique_points.size(); ++i) {
            vertices.push_back(node_type());
        }
    }
    
void rebuild_member_nodes_from_edges() {
    for (auto& poly : polygons) {
        const auto& edges = poly.m_member_edges;
        if (edges.empty()) continue;

        // Build adjacency: node -> list of connected nodes via edges
        std::unordered_map<size_t, std::vector<size_t>> adj;
        for (const auto& e : edges) {
            adj[e[0]].push_back(e[1]);
            adj[e[1]].push_back(e[0]);
        }

        // Walk the chain starting from the first node of the first edge
        size_t start = edges.begin()->operator[](0);
        std::vector<size_t> ordered;
        ordered.reserve(edges.size());

        size_t prev = std::numeric_limits<size_t>::max();
        size_t cur  = start;

        for (size_t step = 0; step < edges.size(); ++step) {
            ordered.push_back(cur);
            const auto& neighbors = adj[cur];
            size_t next = std::numeric_limits<size_t>::max();
            for (size_t nb : neighbors) {
                if (nb != prev) { next = nb; break; }
            }
            if (next == std::numeric_limits<size_t>::max()) break;
            prev = cur;
            cur  = next;
        }

        poly.m_member_nodes = ordered;
    }
}

void rebuild_all_from_nodes_and_edges() {

    // Step 1: rebuild m_member_edges from consecutive node pairs
    for (auto& poly : polygons) {
        poly.m_member_edges.clear();
        const auto& nodes = poly.m_member_nodes;
        size_t nn = nodes.size();
        for (size_t k = 0; k < nn; ++k) {
            std::array<size_t,2> edge = {nodes[k], nodes[(k+1) % nn]};
            validate_edge(edge);
            poly.m_member_edges.insert(edge);
        }
    }

    // Step 2: rebuild facets as the union of all polygon edges
    std::set<std::array<size_t,2>> facet_set;
    for (auto& poly : polygons)
        for (auto& e : poly.m_member_edges)
            facet_set.insert(e);
    facets.assign(facet_set.begin(), facet_set.end());
}

std::array<size_t, 4> extract_quad_corners(size_t cell_index) const {

    const auto& nodes = polygons[cell_index].m_member_nodes;

    // Compute centroid of all nodes (including hanging nodes)
    T cx = T(0), cy = T(0);
    for (size_t nid : nodes) {
        cx += points[nid].x();
        cy += points[nid].y();
    }
    cx /= T(nodes.size());
    cy /= T(nodes.size());

    // The 4 corners are the farthest node in each of the 4 quadrants.
    // Hanging nodes are always closer to the centroid than true corners,
    // so the farthest node per quadrant is always a true corner.
    //
    // Quadrant mapping:
    //   q=0 : dx<0, dy<0 -> bottom-left
    //   q=1 : dx>0, dy<0 -> bottom-right
    //   q=2 : dx<0, dy>0 -> top-left
    //   q=3 : dx>0, dy>0 -> top-right

    std::array<size_t, 4> corners = {0, 0, 0, 0};
    std::array<T, 4>      best    = {T(-1), T(-1), T(-1), T(-1)};

    for (size_t nid : nodes) {
        T dx   = points[nid].x() - cx;
        T dy   = points[nid].y() - cy;
        T dist = dx*dx + dy*dy;

        int q = (dx >= T(0) ? 1 : 0) + (dy >= T(0) ? 2 : 0);

        if (dist > best[q]) {
            best[q]    = dist;
            corners[q] = nid;
        }
    }

    // Return in CCW order: n0=bot-left, n1=bot-right, n2=top-right, n3=top-left
    return {corners[0],   // bot-left  (dx<0, dy<0) q=0
            corners[1],   // bot-right (dx>0, dy<0) q=1
            corners[3],   // top-right (dx>0, dy>0) q=3
            corners[2]};  // top-left  (dx<0, dy>0) q=2
}

void refine_quad_cell(size_t cell_index,
                      const std::array<size_t, 4>& corners) {

    assert(cell_index < polygons.size());

    const int  parent_level = polygons[cell_index].m_refinement_level;
    const int  parent_mat   = polygons[cell_index].m_material;
    const bool parent_elas  = polygons[cell_index].m_elastic_material;

    // ----------------------------------------------------------------
    // 1. The 4 corners in CCW order (precomputed before any refinement):
    //
    //   n3 (top-left) ---- n2 (top-right)
    //        |                    |
    //   n0 (bot-left) ---- n1 (bot-right)
    // ----------------------------------------------------------------

    const size_t n0 = corners[0]; // bottom-left
    const size_t n1 = corners[1]; // bottom-right
    const size_t n2 = corners[2]; // top-right
    const size_t n3 = corners[3]; // top-left

    // ----------------------------------------------------------------
    // 2. Create the 4 edge midpoints and the barycenter
    //
    //   n3 ---m23--- n2
    //   |      |      |
    //  m03----c-----m12
    //   |      |      |
    //   n0 ---m01--- n1
    // ----------------------------------------------------------------

    auto get_or_create_midpoint = [&](size_t na, size_t nb) -> size_t {
        point_type mid(
            (points[na].x() + points[nb].x()) * T(0.5),
            (points[na].y() + points[nb].y()) * T(0.5)
        );
        constexpr T tol = T(1e-14);
        for (size_t k = 0; k < points.size(); ++k) {
            T dx = points[k].x() - mid.x();
            T dy = points[k].y() - mid.y();
            if (dx*dx + dy*dy < tol*tol)
                return k;
        }
        size_t new_id = points.size();
        points.push_back(mid);
        vertices.push_back(node_type(disk::point_identifier<2>(new_id)));
        return new_id;
    };

    const size_t m01 = get_or_create_midpoint(n0, n1);
    const size_t m12 = get_or_create_midpoint(n1, n2);
    const size_t m23 = get_or_create_midpoint(n2, n3);
    const size_t m03 = get_or_create_midpoint(n0, n3);

    const size_t c = [&]() -> size_t {
        point_type center(
            (points[n0].x() + points[n1].x() + points[n2].x() + points[n3].x()) * T(0.25),
            (points[n0].y() + points[n1].y() + points[n2].y() + points[n3].y()) * T(0.25)
        );
        size_t new_id = points.size();
        points.push_back(center);
        vertices.push_back(node_type(disk::point_identifier<2>(new_id)));
        return new_id;
    }();

    // ----------------------------------------------------------------
    // 3. Update boundary edges
    // ----------------------------------------------------------------

    auto split_boundary_edge = [&](size_t a, size_t mid, size_t b) {
        std::array<size_t,2> full_edge = {a, b};
        validate_edge(full_edge);
        auto it = std::find(boundary_edges.begin(), boundary_edges.end(), full_edge);
        if (it != boundary_edges.end()) {
            boundary_edges.erase(it);
            std::array<size_t,2> e1 = {a, mid}; validate_edge(e1);
            std::array<size_t,2> e2 = {mid, b}; validate_edge(e2);
            boundary_edges.push_back(e1);
            boundary_edges.push_back(e2);
        }
    };

    split_boundary_edge(n0, m01, n1);
    split_boundary_edge(n1, m12, n2);
    split_boundary_edge(n2, m23, n3);
    split_boundary_edge(n3, m03, n0);

    // ----------------------------------------------------------------
    // 4. Build the 4 child cells
    //
    //   Q0 : n0,  m01, c,   m03
    //   Q1 : m01, n1,  m12, c
    //   Q2 : c,   m12, n2,  m23
    //   Q3 : m03, c,   m23, n3
    // ----------------------------------------------------------------

    auto make_quad = [&](size_t a, size_t b, size_t d_, size_t e_) -> polygon_2d {
        polygon_2d q;
        q.m_member_nodes     = {a, b, d_, e_};
        q.m_material         = parent_mat;
        q.m_elastic_material = parent_elas;
        q.m_refinement_level = parent_level + 1;

        auto add_edge = [&](size_t x, size_t y) {
            std::array<size_t,2> edge = {x, y};
            validate_edge(edge);
            q.m_member_edges.insert(edge);
            facets.push_back(edge);
        };
        add_edge(a, b);
        add_edge(b, d_);
        add_edge(d_, e_);
        add_edge(e_, a);
        return q;
    };

    polygon_2d q0 = make_quad(n0,  m01, c,   m03);
    polygon_2d q1 = make_quad(m01, n1,  m12, c  );
    polygon_2d q2 = make_quad(c,   m12, n2,  m23);
    polygon_2d q3 = make_quad(m03, c,   m23, n3 );

    // ----------------------------------------------------------------
    // 5. Erase parent cell and append the 4 children
    // ----------------------------------------------------------------

    polygons.erase(polygons.begin() + cell_index);
    polygons.push_back(q0);
    polygons.push_back(q1);
    polygons.push_back(q2);
    polygons.push_back(q3);
}

// In struct polygon_2d, add:
//   bool m_to_refine = false;

void refine_cells(const std::vector<size_t>& cell_indices,
                  int refinement_level) {

    if (refinement_level <= 0) return;

    // Mark the initial targets
    for (size_t idx : cell_indices) {
        assert(idx < polygons.size());
        polygons[idx].m_to_refine = true;
    }

    for (int pass = 0; pass < refinement_level; ++pass) {

        // Collect current indices of all marked cells
        std::vector<size_t> targets;
        for (size_t i = 0; i < polygons.size(); ++i)
            if (polygons[i].m_to_refine)
                targets.push_back(i);

        // std::cout << "[refine_cells] pass " << pass
        //           << " | polygons=" << polygons.size()
        //           << " | targets=" << targets.size() << std::endl;

        // Sort descending so that erasing cell at index i does not
        // shift the indices of cells not yet processed in this pass.
        std::sort(targets.begin(), targets.end(), std::greater<size_t>());

        // Extract all corners BEFORE any refinement
        std::vector<std::array<size_t,4>> all_corners(targets.size());
        for (size_t i = 0; i < targets.size(); ++i)
            all_corners[i] = extract_quad_corners(targets[i]);

        // Clear flags — children will be re-marked below
        for (size_t idx : targets)
            polygons[idx].m_to_refine = false;

        for (size_t i = 0; i < targets.size(); ++i) {

            // Since targets is sorted descending and we erase from high
            // to low, erasing targets[i] does not shift any targets[j]
            // with j > i (they all have lower indices). No adjustment needed.
            size_t idx = targets[i];

            refine_quad_cell(idx, all_corners[i]);

            // Mark the 4 children for the next pass
            size_t new_end = polygons.size();
            for (size_t ci = new_end - 4; ci < new_end; ++ci)
                polygons[ci].m_to_refine = true;
        }

        // Step 1: initial facet set from all polygon edges
        {
            std::set<std::array<size_t,2>> facet_set;
            for (auto& poly : polygons)
                for (auto& e : poly.m_member_edges)
                    facet_set.insert(e);
            facets.assign(facet_set.begin(), facet_set.end());
        }

        // Step 2: update stale edges in neighbor polygons
        {
            std::set<std::array<size_t,2>> facet_set(facets.begin(),
                                                      facets.end());
            for (auto& poly : polygons) {
                std::set<std::array<size_t,2>> updated_edges;
                for (const auto& e : poly.m_member_edges) {
                    if (facet_set.count(e)) {
                        updated_edges.insert(e);
                    } else {
                        bool found = false;
                        for (const auto& f : facet_set) {
                            size_t shared = std::numeric_limits<size_t>::max();
                            if      (f[0] == e[0]) shared = f[1];
                            else if (f[1] == e[0]) shared = f[0];
                            else continue;
                            if (shared == e[1]) continue;

                            const auto& pa = points[e[0]];
                            const auto& pb = points[e[1]];
                            const auto& pm = points[shared];
                            T ex = pb.x()-pa.x(), ey = pb.y()-pa.y();
                            T fx = pm.x()-pa.x(), fy = pm.y()-pa.y();
                            T cross = ex*fy - ey*fx;
                            T len2  = ex*ex + ey*ey;
                            constexpr T tol = T(1e-10);
                            if (std::abs(cross) > tol*std::sqrt(len2)) continue;
                            T t = (fx*ex + fy*ey) / len2;
                            if (t <= T(0) || t >= T(1)) continue;

                            std::array<size_t,2> h1 = {e[0], shared};
                            std::array<size_t,2> h2 = {shared, e[1]};
                            validate_edge(h1); validate_edge(h2);
                            if (facet_set.count(h1) && facet_set.count(h2)) {
                                updated_edges.insert(h1);
                                updated_edges.insert(h2);
                                found = true;
                                break;
                            }
                        }
                        if (!found) updated_edges.insert(e);
                    }
                }
                poly.m_member_edges = updated_edges;
            }
        }

        // Step 3: rebuild m_member_nodes by chaining m_member_edges
        for (auto& poly : polygons) {
            const auto& edges = poly.m_member_edges;
            if (edges.empty()) continue;

            std::unordered_map<size_t, std::vector<size_t>> adj;
            for (const auto& e : edges) {
                adj[e[0]].push_back(e[1]);
                adj[e[1]].push_back(e[0]);
            }

            size_t start = edges.begin()->operator[](0);
            std::vector<size_t> ordered;
            ordered.reserve(edges.size());
            size_t prev = std::numeric_limits<size_t>::max();
            size_t cur  = start;

            for (size_t step = 0; step < edges.size(); ++step) {
                ordered.push_back(cur);
                const auto& nbrs = adj[cur];
                size_t next = std::numeric_limits<size_t>::max();
                for (size_t nb : nbrs)
                    if (nb != prev) { next = nb; break; }
                if (next == std::numeric_limits<size_t>::max()) break;
                prev = cur;
                cur  = next;
            }
            poly.m_member_nodes = ordered;
        }

        // Step 4: rebuild facets from consistent m_member_edges
        {
            std::set<std::array<size_t,2>> facet_set;
            for (auto& poly : polygons)
                for (auto& e : poly.m_member_edges)
                    facet_set.insert(e);
            facets.assign(facet_set.begin(), facet_set.end());
        }
    }

    // Clear all flags
    for (auto& poly : polygons)
        poly.m_to_refine = false;
}

void move_to_mesh_storage(mesh_type& msh){
                
                auto storage = msh.backend_storage();
                storage->points = std::move(points);
                storage->nodes = std::move(vertices);
                
                std::vector<edge_type> edges;
                edges.reserve(facets.size());
                for (size_t i = 0; i < facets.size(); i++)
                {
                    assert(facets[i][0] < facets[i][1]);
                    auto node1 = typename node_type::id_type(facets[i][0]);
                    auto node2 = typename node_type::id_type(facets[i][1]);
                    auto e = edge_type(node1, node2);
                    edges.push_back(e);
                }
                std::sort(edges.begin(), edges.end());
                
                storage->boundary_info.resize(edges.size());
                for (size_t i = 0; i < boundary_edges.size(); i++)
                {
                    assert(boundary_edges[i][0] < boundary_edges[i][1]);
                    auto node1 = typename node_type::id_type(boundary_edges[i][0]);
                    auto node2 = typename node_type::id_type(boundary_edges[i][1]);
                    auto e = edge_type(node1, node2);
                    auto position = find_element_id(edges.begin(), edges.end(), e);
                    if (position.first == false)
                    {
                        std::cout << "Bad bug at " << __FILE__ << "("
                        << __LINE__ << ")" << std::endl;
                        // -- DIAGNOSTIC --
                        std::cout << "  [diag] missing boundary edge ["
                        << boundary_edges[i][0] << ", "
                        << boundary_edges[i][1] << "]" << std::endl;
                        // -- END DIAGNOSTIC --
                        return;
                    }
                    disk::boundary_descriptor bi{0, true};
                    storage->boundary_info.at(position.second) = bi;
                }
                
                storage->edges = std::move(edges);
                
                std::vector<surface_type> surfaces;
                surfaces.reserve( polygons.size() );
                
                size_t pi = 0; // -- DIAGNOSTIC: polygon index --
                for (auto& p : polygons)
                {
                    std::vector<typename edge_type::id_type> surface_edges;
                    for (auto& e : p.m_member_edges)
                    {
                        assert(e[0] < e[1]);
                        auto n1 = typename node_type::id_type(e[0]);
                        auto n2 = typename node_type::id_type(e[1]);
                        edge_type edge(n1, n2);
                        auto edge_id = find_element_id(storage->edges.begin(),
                        storage->edges.end(), edge);
                        if (!edge_id.first)
                        {
                            std::cout << "Bad bug at " << __FILE__ << "("
                            << __LINE__ << ")" << std::endl;
                            // -- DIAGNOSTIC --
                            std::cout << "  [diag] polygon " << pi
                            << " | missing edge ["
                            << e[0] << ", " << e[1] << "]" << std::endl;
                            std::cout << "  [diag] polygon nodes: ";
                            for (auto n : p.m_member_nodes)
                            std::cout << n << " ";
                            std::cout << std::endl;
                            std::cout << "  [diag] polygon edges: ";
                            for (auto& pe : p.m_member_edges)
                            std::cout << "[" << pe[0] << "," << pe[1] << "] ";
                            std::cout << std::endl;
                            // -- END DIAGNOSTIC --
                            return;
                        }
                        surface_edges.push_back(edge_id.second);
                    }
                    auto surface = surface_type(surface_edges);
                    surface.set_point_ids(p.m_member_nodes.begin(), p.m_member_nodes.end());
                    surfaces.push_back( surface );
                    pi++; // -- DIAGNOSTIC --
                }
                
                std::sort(surfaces.begin(), surfaces.end());
                storage->surfaces = std::move(surfaces);
            }
            
void set_translation_data(T x_t, T y_t){
                m_x_t = x_t;
                m_y_t = y_t;
            }
            
            size_t get_nx(){
                return m_nx;
            }
            
            size_t get_ny(){
                return m_ny;
            }
            
            // Print in log file relevant mesh information
            void print_log_file(){
                fitted_geometry_builder<mesh_type>::m_n_elements = polygons.size();
                std::ofstream file;
        file.open (fitted_geometry_builder<mesh_type>::m_log_file.c_str());
        file << "Number of surfaces : " << polygons.size() << std::endl;
        file << "Number of skeleton edges : " << skeleton_edges.size() << std::endl;
        file << "Number of boundary edges : " << boundary_edges.size() << std::endl;
        file << "Number of vertices : " << vertices.size() << std::endl;
        file.close();
    }
    
};


template<typename T>
class polygon_2d_mesh_reader : 

public fitted_geometry_builder<disk::generic_mesh<T, 2>> {

    typedef disk::generic_mesh<T,2>                 mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;   

public:

    struct polygon_2d {
        std::vector<size_t>                 m_member_nodes;
        std::set<std::array<size_t, 2>>     m_member_edges;
        int                                 m_material;
        bool                                m_elastic_material;

        bool operator<(const polygon_2d & other) {
            return m_member_nodes < other.m_member_nodes;
        }
    };

    // struct FaceInfo {
    //     size_t face_id;
    //     size_t cell1_id = -1;
    //     size_t cell2_id = -1;
    // };

    std::vector<point_type>                         points;
    std::vector<node_type>                          vertices;
    std::vector<std::array<size_t, 2>>              facets;
    std::vector<std::array<size_t, 2>>              interfaces;
    std::vector<std::array<size_t, 2>>              skeleton_edges;
    std::vector<std::array<size_t, 2>>              boundary_edges;
    std::vector<polygon_2d>                         polygons;
    std::string poly_mesh_file;
    std::set<size_t> bc_points;
    
    void clear_storage() {
        points.clear();
        vertices.clear();
        skeleton_edges.clear();
        boundary_edges.clear();
        polygons.clear();
        bc_points.clear();
    }
    
    void reserve_storage(){
        
        std::ifstream input;
        input.open(poly_mesh_file.c_str());
        
        size_t n_points, n_polygons, n_bc_curves, n_bc_edges, n_edges;
        if (input.is_open()) {
            std::string line;
            std::getline(input, line);
            std::stringstream(line) >> n_points >> n_polygons >> n_bc_curves;
            for(size_t id = 0; id < n_points; id++){
                if(std::getline(input, line)){

                }
                else{
                    break;
                }
            }
            
            n_edges = 0;
            n_bc_edges = 0;
            size_t n_polygon_vertices;
            for(size_t surface_id=0; surface_id < n_polygons; surface_id++)
            {
                if(std::getline(input, line)){
                    std::stringstream(line) >> n_polygon_vertices;
                    n_edges += n_polygon_vertices;
                  }
                  else{
                      break;
                  }
            }
            
            bc_points.clear();
            size_t bc_point_id;
            for(size_t bc_id=0; bc_id < n_bc_curves; bc_id++)
            {
                if(std::getline(input, line)){
                    std::stringstream input_line(line);
                    while(!input_line.eof()){
                        input_line >> bc_point_id;
                        bc_point_id--;
                        bc_points.insert(bc_point_id);
                    }
                }
                else{
                  break;
                }
            }
            n_bc_edges = bc_points.size();
            
            points.reserve(n_points);
            vertices.reserve(n_points);
            
            size_t n_skel_edges = n_edges - n_bc_edges;
            skeleton_edges.reserve(n_skel_edges);
            boundary_edges.reserve(n_bc_edges);
            polygons.reserve(n_polygons);
            
        }
    }
            
    void validate_edge(std::array<size_t, 2> & edge){
        assert(edge[0] != edge[1]);
        if (edge[0] > edge[1]){
            std::swap(edge[0], edge[1]);
        }
    }

    polygon_2d_mesh_reader() : fitted_geometry_builder<mesh_type>()
    {
        fitted_geometry_builder<mesh_type>::m_dimension = 2;
    }
    
    void set_poly_mesh_file(std::string mesh_file){
        poly_mesh_file = mesh_file;
    }
    
    // build the mesh
    bool build_mesh(){
        
        clear_storage();
        reserve_storage();
        
        std::ifstream input;
        input.open(poly_mesh_file.c_str());
        
        size_t n_points, n_polygons, n_bc_curves;
        if (input.is_open()) {
            std::string line;
            std::getline(input, line);
            std::stringstream(line) >> n_points >> n_polygons >> n_bc_curves;

            T xv, yv;
            for(size_t id = 0; id < n_points; id++){
                if(std::getline(input, line)){
                    std::stringstream(line) >> xv >> yv;
                    point_type point(xv, yv);
                    points.push_back(point);
                    vertices.push_back(node_type(disk::point_identifier<2>(id)));
                }
                else{
                    break;
                }
            }
            
            size_t n_polygon_vertices, id;
            for(size_t surface_id=0; surface_id < n_polygons; surface_id++) {

                if(std::getline(input, line)){
                    std::stringstream input_line(line);
                    input_line >> n_polygon_vertices;
                    
                    polygon_2d polygon;
                    std::vector<size_t> member_nodes;
                    for (size_t i = 0; i < n_polygon_vertices; i++) {
                        input_line >> id;
                        id--;
                        member_nodes.push_back(id);
                    }
                    polygon.m_member_nodes = member_nodes;
                    assert(member_nodes.size() == n_polygon_vertices);

                    std::set< std::array<size_t, 2> > member_edges;
                    std::array<size_t, 2> edge;
                    for (size_t i = 0; i < member_nodes.size(); i++) {
                        
                        if (i == n_polygon_vertices - 1) {
                            edge = {member_nodes[i],member_nodes[0]};
                        }

                        else {
                            edge = {member_nodes[i],member_nodes[i+1]};
                        }
                        
                        validate_edge(edge);
                        facets.push_back( edge );
                        member_edges.insert(edge);
                        
                        bool is_bc_point_l_Q = bc_points.find(edge.at(0)) != bc_points.end();
                        bool is_bc_point_r_Q = bc_points.find(edge.at(1)) != bc_points.end();
                        if (is_bc_point_l_Q && is_bc_point_r_Q) {
                            boundary_edges.push_back( edge );
                        }
                    }
                    
                    polygon.m_member_edges = member_edges;

                    polygons.push_back( polygon );
                    
                  }

                  else {
                      break;
                  }
            }
        }
          
        // Duplicated facets are eliminated
        std::sort( facets.begin(), facets.end() );
        facets.erase( std::unique( facets.begin(), facets.end() ), facets.end() );

        // std::cout << bold << red << std::endl << std::endl;
        // std::cout << "Debug mesh" << std::endl;
        // std::cout << reset << "Points: " << vertices.size() << std::endl;
        // std::cout << reset << "Polygons: " << polygons.size() << std::endl;
        // std::cout << reset << "Boundaries: " << boundary_edges.size() << std::endl;
        // std::cout << reset << "Facets: " << facets.size() << std::endl;

        return true;
    }
    
       // build the mesh
    bool build_bassin(){
        
        clear_storage();
        reserve_storage();
        
        std::ifstream input;
        input.open(poly_mesh_file.c_str());
        
        size_t n_points, n_polygons, n_bc_curves;
        if (input.is_open()) {
            std::string line;
            std::getline(input, line);
            std::stringstream(line) >> n_points >> n_polygons >> n_bc_curves;

            T xv, yv;
            for(size_t id = 0; id < n_points; id++){
                if(std::getline(input, line)){
                    std::stringstream(line) >> xv >> yv;
                    point_type point(xv, yv);
                    points.push_back(point);
                    vertices.push_back(node_type(disk::point_identifier<2>(id)));
                }
                else{
                    break;
                }
            }
            
            size_t n_polygon_vertices, id;
            for(size_t surface_id=0; surface_id < n_polygons; surface_id++) {

                if(std::getline(input, line)){
                    std::stringstream input_line(line);
                    input_line >> n_polygon_vertices;
                    
                    polygon_2d polygon;
                    std::vector<size_t> member_nodes;
                    int material;
                    for (size_t i = 0; i < n_polygon_vertices; i++) {
                        input_line >> id;
                        id--;
                        member_nodes.push_back(id);
                    }
                    polygon.m_member_nodes = member_nodes;
                    assert(member_nodes.size() == n_polygon_vertices);

                    // FaceInfo face_info;
                    std::set< std::array<size_t, 2> > member_edges;
                    std::array<size_t, 2> edge;
                    for (size_t i = 0; i < member_nodes.size(); i++) {
                        
                        if (i == n_polygon_vertices - 1) {
                            edge = {member_nodes[i],member_nodes[0]};
                        }

                        else {
                            edge = {member_nodes[i],member_nodes[i+1]};
                        }
                        
                        validate_edge(edge);
                        facets.push_back( edge );
                        member_edges.insert(edge);


                        bool is_bc_point_l_Q = bc_points.find(edge.at(0)) != bc_points.end();
                        bool is_bc_point_r_Q = bc_points.find(edge.at(1)) != bc_points.end();
                        if (is_bc_point_l_Q && is_bc_point_r_Q) {
                            boundary_edges.push_back( edge );
                        }
                    }
                    
                    polygon.m_member_edges = member_edges;
                    polygons.push_back( polygon );
                    
                  }
                  else {
                      break;
                  }
            }


            for(size_t bord=0; bord < 4; bord++) {
                if(std::getline(input, line)){
                  std::stringstream input_line(line);
                }
            }


            for(size_t polygon=0; polygon < n_polygons; polygon++) {
                size_t material;
                if(std::getline(input, line)){
                  std::stringstream input_line(line);
                  input_line >> material;
                }
                polygons[polygon].m_material = material;
            }
            
            
        }
        
        // // Duplicated facets are eliminated
        // std::sort( facets.begin(), facets.end() );
        // facets.erase( std::unique( facets.begin(), facets.end() ), facets.end() );
        
        // // ---------------------- DEBUG ----------------------
        // std::cout << "Debug Facets & Cells Attachments\n";
        // std::cout << "Total facets: " << facets.size() << "\n";
        
        // // Construire un map des arêtes → cellules qui les contiennent
        // std::map<std::array<size_t,2>, std::vector<size_t>> edge_to_cells;
        // for(size_t c=0; c<polygons.size(); ++c){
        //     for(const auto& e : polygons[c].m_member_edges){
        //         edge_to_cells[e].push_back(c);
        //     }
        // }
        
        // for(const auto& e : facets){
        //     std::cout << "Facet [" << e[0] << ", " << e[1] << "]\n";
        //     // Coordonnées des points de la face
        //     std::cout << "  Nodes:\n";
        //     std::cout << "    " << e[0] << " : (" << points[e[0]].x() << ", " << points[e[0]].y() << ")\n";
        //     std::cout << "    " << e[1] << " : (" << points[e[1]].x() << ", " << points[e[1]].y() << ")\n";
        
        //     // Cellules qui contiennent cette face
        //     auto it = edge_to_cells.find(e);
        //     if(it != edge_to_cells.end()){
        //         std::cout << "  Attached cells: ";
        //         for(auto cidx : it->second){
        //             std::cout << cidx << " ";
        //         }
        //         std::cout << "\n";
        //     }
        // }
        
        // std::cout << "Debug terminé.\n";
        // // --------------------------------------------------
        
        // std::cout << bold << red << std::endl << std::endl;
        // std::cout << "Debug mesh" << std::endl;
        // std::cout << reset << "Points: " << vertices.size() << std::endl;
        // std::cout << reset << "Polygons: " << polygons.size() << std::endl;
        // std::cout << reset << "Boundaries: " << boundary_edges.size() << std::endl;
        // std::cout << reset << "Facets: " << facets.size() << std::endl;
        
        return true;
    }
    
    void move_to_mesh_storage(mesh_type& msh){
        
        auto storage = msh.backend_storage();
        storage->points = std::move(points);
        storage->nodes = std::move(vertices);
        
        std::vector<edge_type> edges;
        edges.reserve(facets.size());
        for (size_t i = 0; i < facets.size(); i++) {
            assert(facets[i][0] < facets[i][1]);
            auto node1 = typename node_type::id_type(facets[i][0]);
            auto node2 = typename node_type::id_type(facets[i][1]);
            
            auto e = edge_type(node1, node2);
            
            // e.set_point_ids(facets[i].begin(), facets[i].end());
            edges.push_back(e);
        }
        std::sort(edges.begin(), edges.end());
        
        storage->boundary_info.resize(edges.size());
        for (size_t i = 0; i < boundary_edges.size(); i++) {
            assert(boundary_edges[i][0] < boundary_edges[i][1]);
            auto node1 = typename node_type::id_type(boundary_edges[i][0]);
            auto node2 = typename node_type::id_type(boundary_edges[i][1]);
            auto e = edge_type(node1, node2);
            auto position = find_element_id(edges.begin(), edges.end(), e);
            if (position.first == false) {
                std::cout << "Bad bug at " << __FILE__ << "("
                << __LINE__ << ")" << std::endl;
                return;
            }
            disk::boundary_descriptor bi{0, true};
            storage->boundary_info.at(position.second) = bi;
        }
        
        storage->edges = std::move(edges);
        
        std::vector<surface_type> surfaces;
        surfaces.reserve( polygons.size() );
        
        std::vector<int> materials;
        materials.reserve( polygons.size() );
        
        int compteur = 0;
        for (auto& p : polygons) {
            compteur = compteur + 1;
            std::vector<typename edge_type::id_type> surface_edges;
            for (auto& e : p.m_member_edges) {
                assert(e[0] < e[1]);
                auto n1 = typename node_type::id_type(e[0]);
                auto n2 = typename node_type::id_type(e[1]);
                
                edge_type edge(n1, n2);
                auto edge_id = find_element_id(storage->edges.begin(),
                storage->edges.end(), edge);
                if (!edge_id.first)
                {
                    std::cout << "Bad bug at " << __FILE__ << "("
                    << __LINE__ << ")" << std::endl;
                    return;
                }
                
                surface_edges.push_back(edge_id.second);
                
                
            }
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(p.m_member_nodes.begin(), p.m_member_nodes.end());
            surfaces.push_back( surface );
            materials.push_back(p.m_material);
        }
        
        ////////////////////////////// Trie pour les matériaux
        std::vector<std::pair<surface_type, int>> pairedVec;
        for (size_t i = 0; i < surfaces.size(); ++i) {
            pairedVec.push_back(std::make_pair(surfaces[i], materials[i]));
        }
        // Trier les paires en fonction des valeurs (les surfaces dans ce cas)
        std::sort(pairedVec.begin(), pairedVec.end(), [](const auto& left, const auto& right) {
            return left.first < right.first;
        });
        // Réorganiser le vecteur material en fonction des indices triés
        std::vector<int> sortedMaterial;
        for (const auto& pair : pairedVec) {
            sortedMaterial.push_back(pair.second);
        }
        // Mettre à jour le vecteur material
        for (size_t i = 0; i < surfaces.size(); ++i) {
            polygons[i].m_material = sortedMaterial[i];
            if (polygons[i].m_material == 1 || polygons[i].m_material == 2) {
                polygons[i].m_elastic_material = true;
            }
            else {
                polygons[i].m_elastic_material = false;
            }
        }
        //////////////////////////////////////////////////
        
        std::sort(surfaces.begin(), surfaces.end());
        storage->surfaces = std::move(surfaces);
        
    }
    
    // Print in log file relevant mesh information
    void print_log_file(){
        fitted_geometry_builder<mesh_type>::m_n_elements = polygons.size();
        std::ofstream file;
        file.open (fitted_geometry_builder<mesh_type>::m_log_file.c_str());
        file << "Number of polygons : " << polygons.size() << std::endl;
        file << "Number of skeleton edges : " << skeleton_edges.size() << std::endl;
        file << "Number of boundary edges : " << boundary_edges.size() << std::endl;
        file << "Number of vertices : " << vertices.size() << std::endl;
        file.close();
    }
    
    const std::vector<point_type>& getPoints() const {
        return points;
    }
    
    const std::vector<node_type>& getVertices() const {
        return vertices;
    }
    
    const std::vector<std::array<size_t, 2>>& getFacets() const {
        return facets;
    }
    
    const std::vector<std::array<size_t, 2>>& getSkeletonEdges() const {
        return skeleton_edges;
    }
    
    const std::vector<std::array<size_t, 2>>& getBoundaryEdges() const {
        return boundary_edges;
    }
    
    const std::vector<polygon_2d>& getPolygons() const {
        return polygons;
    }
    
    void remove_duplicate_points() {
        std::unordered_map<size_t, size_t> point_mapping;  // Map old point indices to new point indices
        std::vector<point_type> unique_points;  // Unique points
        
        for (size_t i = 0; i < points.size(); ++i) {
            if (point_mapping.find(i) == point_mapping.end()) {
                // This is a new unique point
                point_mapping[i] = unique_points.size();
                unique_points.push_back(points[i]);
            }
        }
        
        // Update cell and face data with unique point indices
        for (auto& polygon : polygons) {
            for (size_t i = 0; i < polygon.m_member_nodes.size(); ++i) {
                size_t old_point_index = polygon.m_member_nodes[i];
                polygon.m_member_nodes[i] = point_mapping[old_point_index];
            }
            
            std::set<std::array<size_t, 2>> updated_edges;
            for (const auto& edge : polygon.m_member_edges) {
                size_t old_point1 = edge[0];
                size_t old_point2 = edge[1];
                size_t new_point1 = point_mapping[old_point1];
                size_t new_point2 = point_mapping[old_point2];
                updated_edges.insert({new_point1, new_point2});
            }
            polygon.m_member_edges = updated_edges;
        }
        
        // Update the internal mesh_reader data with the unique points
        points = unique_points;
        
        // Update the vertices
        vertices.clear();
        for (size_t i = 0; i < unique_points.size(); ++i) {
            vertices.push_back(node_type());
        }
    }
    
    void refine_cells(const std::vector<size_t>& cell_indices, int refinement_level) {
        
        if (refinement_level <= 0) return;
        
        // Work on a sorted copy (descending) so that erasing cell at
        // index i does not shift the indices of cells not yet processed
        // in the same pass.
        std::vector<size_t> targets(cell_indices.begin(), cell_indices.end());
        std::sort(targets.begin(), targets.end(), std::greater<size_t>());
        
        for (int pass = 0; pass < refinement_level; ++pass) {
            
            if (pass == 0) {
                // First pass: refine the originally requested cells
                for (size_t idx : targets)
                refine_quad_cell(idx);
            } else {
                // Subsequent passes: the children of the previous pass
                // are always the last 4 * targets.size() cells appended.
                size_t n        = polygons.size();
                size_t n_new    = targets.size() * 4;
                size_t start    = n - n_new;
                
                // Collect and sort descending so erasure stays safe
                std::vector<size_t> next_targets(n_new);
                std::iota(next_targets.begin(), next_targets.end(), start);
                std::sort(next_targets.begin(), next_targets.end(), std::greater<size_t>());
                
                targets = next_targets;
                for (size_t idx : targets)
                refine_quad_cell(idx);
            }
        }
    }
    
    void refine_quad_cell(size_t cell_index) {
        
        assert(cell_index < polygons.size());
        
        // Copy all needed data upfront before any modification of
        // polygons[] or facets[] that could invalidate references.
        const std::vector<size_t>              parent_nodes = polygons[cell_index].m_member_nodes;
        const std::set<std::array<size_t, 2>>  parent_edges = polygons[cell_index].m_member_edges;
        const int                              parent_level = polygons[cell_index].m_refinement_level;
        const int                              parent_mat   = polygons[cell_index].m_material;
        const bool                             parent_elas  = polygons[cell_index].m_elastic_material;
        
        assert(parent_nodes.size() == 4);
        
        // ----------------------------------------------------------------
        // 1. Retrieve the 4 vertices in order
        // ----------------------------------------------------------------
        const size_t n0 = parent_nodes[0];
        const size_t n1 = parent_nodes[1];
        const size_t n2 = parent_nodes[2];
        const size_t n3 = parent_nodes[3];
        
        // ----------------------------------------------------------------
        // 2. Create the 4 edge midpoints and the barycenter
        // ----------------------------------------------------------------
        
        auto get_or_create_midpoint = [&](size_t na, size_t nb) -> size_t {
            point_type mid(
                (points[na].x() + points[nb].x()) * T(0.5),
                (points[na].y() + points[nb].y()) * T(0.5)
            );
            constexpr T tol = T(1e-14);
            for (size_t k = 0; k < points.size(); ++k) {
                T dx = points[k].x() - mid.x();
                T dy = points[k].y() - mid.y();
                if (dx*dx + dy*dy < tol*tol)
                return k;
            }
            size_t new_id = points.size();
            points.push_back(mid);
            vertices.push_back(node_type(disk::point_identifier<2>(new_id)));
            return new_id;
        };
        
        const size_t m01 = get_or_create_midpoint(n0, n1);
        const size_t m12 = get_or_create_midpoint(n1, n2);
        const size_t m23 = get_or_create_midpoint(n2, n3);
        const size_t m03 = get_or_create_midpoint(n0, n3);
        
        const size_t c = [&]() -> size_t {
            point_type center(
                (points[n0].x() + points[n1].x() + points[n2].x() + points[n3].x()) * T(0.25),
                (points[n0].y() + points[n1].y() + points[n2].y() + points[n3].y()) * T(0.25)
            );
            size_t new_id = points.size();
            points.push_back(center);
            vertices.push_back(node_type(disk::point_identifier<2>(new_id)));
            return new_id;
        }();
        
        // ----------------------------------------------------------------
        // 3. Update boundary edges
        // ----------------------------------------------------------------
        
        auto split_boundary_edge = [&](size_t a, size_t mid, size_t b) {
            std::array<size_t,2> full_edge = {a, b};
            validate_edge(full_edge);
            auto it = std::find(boundary_edges.begin(), boundary_edges.end(), full_edge);
            if (it != boundary_edges.end()) {
                boundary_edges.erase(it);
                std::array<size_t,2> e1 = {a, mid}; validate_edge(e1);
                std::array<size_t,2> e2 = {mid, b}; validate_edge(e2);
                boundary_edges.push_back(e1);
                boundary_edges.push_back(e2);
            }
        };
        
        split_boundary_edge(n0, m01, n1);
        split_boundary_edge(n1, m12, n2);
        split_boundary_edge(n2, m23, n3);
        split_boundary_edge(n3, m03, n0);
        
        // ----------------------------------------------------------------
        // 4. Update neighboring cells
        // ----------------------------------------------------------------
        
        struct EdgeSplit { size_t a, b, mid; };
        std::array<EdgeSplit, 4> splits = {{
            {n0, n1, m01},
            {n1, n2, m12},
            {n2, n3, m23},
            {n0, n3, m03}
        }};
        
        for (auto& [ea, eb, emid] : splits) {
            
            std::array<size_t,2> parent_edge = {ea, eb};
            validate_edge(parent_edge);
            
            for (size_t ci = 0; ci < polygons.size(); ++ci) {
                
                if (ci == cell_index) continue;
                
                polygon_2d& neighbor = polygons[ci];
                
                if (neighbor.m_member_edges.find(parent_edge) == neighbor.m_member_edges.end())
                continue;
                
                neighbor.m_member_edges.erase(parent_edge);
                std::array<size_t,2> e1 = {ea,   emid}; validate_edge(e1);
                std::array<size_t,2> e2 = {emid, eb  }; validate_edge(e2);
                neighbor.m_member_edges.insert(e1);
                neighbor.m_member_edges.insert(e2);
                
                auto& nodes = neighbor.m_member_nodes;
                size_t nn   = nodes.size();
                for (size_t k = 0; k < nn; ++k) {
                    size_t cur  = nodes[k];
                    size_t next = nodes[(k + 1) % nn];
                    if ((cur == ea && next == eb) || (cur == eb && next == ea)) {
                        nodes.insert(nodes.begin() + k + 1, emid);
                        break;
                    }
                }
                break;
            }
        }
        
        // ----------------------------------------------------------------
        // 5. Build the 4 child cells
        // ----------------------------------------------------------------
        
        auto make_quad = [&](size_t a, size_t b, size_t d_, size_t e_) -> polygon_2d {
            polygon_2d q;
            q.m_member_nodes      = {a, b, d_, e_};
            q.m_material          = parent_mat;
            q.m_elastic_material  = parent_elas;
            q.m_refinement_level  = parent_level + 1;
            
            auto add_edge = [&](size_t x, size_t y) {
                std::array<size_t,2> edge = {x, y};
                validate_edge(edge);
                q.m_member_edges.insert(edge);
                facets.push_back(edge);
            };
            add_edge(a, b);
            add_edge(b, d_);
            add_edge(d_, e_);
            add_edge(e_, a);
            return q;
        };
        
        polygon_2d q0 = make_quad(n0,  m01, c,   m03);
        polygon_2d q1 = make_quad(m01, n1,  m12, c  );
        polygon_2d q2 = make_quad(c,   m12, n2,  m23);
        polygon_2d q3 = make_quad(m03, c,   m23, n3 );
        
        // ----------------------------------------------------------------
        // 6. Remove parent facets, erase parent cell, append children
        //    Use the copied parent_edges — poly reference is now unsafe
        // ----------------------------------------------------------------
        
        for (const auto& e : parent_edges) {
            auto it = std::find(facets.begin(), facets.end(), e);
            if (it != facets.end())
            facets.erase(it);
        }
        
        polygons.erase(polygons.begin() + cell_index);
        
        polygons.push_back(q0);
        polygons.push_back(q1);
        polygons.push_back(q2);
        polygons.push_back(q3);
        
        // Remove duplicate facets
        std::sort(facets.begin(), facets.end());
        facets.erase(std::unique(facets.begin(), facets.end()), facets.end());
    }
    
    void refine_cell_poly(size_t cell_index, int nelem = 50, int maxiter = 100) {
        
        if (cell_index >= polygons.size()) {
            std::cerr << "Invalid cell index" << std::endl;
            return;
        }
        
        py::scoped_interpreter guard{};   
        
        py::module bridge = py::module::import("polymesher_bridge");
        
        const auto& poly = polygons[cell_index];
        
        std::vector<std::vector<double>> polygon_points;
        for (auto node_id : poly.m_member_nodes) {
            const auto& pt = points[node_id];
            polygon_points.push_back({pt.x(), pt.y()});
        }
        
        if (polygon_points.front() != polygon_points.back()) {
            polygon_points.push_back(polygon_points.front());
        }
        
        auto result = bridge.attr("generate_mesh")(polygon_points, nelem, maxiter);
        auto result_tuple = result.cast<py::tuple>();  // <-- conversion en tuple
        
        py::array_t<double> Node = result_tuple[0].cast<py::array_t<double>>();
        py::array_t<int> Element = result_tuple[1].cast<py::array_t<int>>();
        
        auto nodes = Node.unchecked<2>();
        auto elems = Element.unchecked<2>();
        
        std::cout << "Refined cell " << cell_index << ":\n";
        std::cout << "Nb nodes = " << nodes.shape(0) << "\n";
        std::cout << "Nb elems = " << elems.shape(0) << "\n";
    }

    
    
    // void refine_cells_with_pypolymesher(const std::vector<size_t>& cell_indices) {
    
    //     // Export des polygones à raffiner
    //     std::ofstream out("cells_to_refine.csv");
    //     for (auto idx : cell_indices) {
    //         const auto& poly = polygons[idx];
    //         for (auto node_id : poly.m_member_nodes) {
    //             const auto& p = points[node_id];
    //             out << p.x() << " " << p.y() << " ";
    //         }
    //         out << std::endl;
    //     }
    //     out.close();
        
    //     // Appel Python
    //     system("python3 /home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/common/pyPolyMesher/refine_cells.py");
        
    //     // Import des cellules raffinées
    //     std::ifstream in("refined_cells.csv");
    //     size_t new_node_id = points.size();
    //     std::string line;
    //     while (std::getline(in, line)) {
    //         std::stringstream ss(line);
    //         std::vector<T> coords;
    //         T x, y;
    //         while (ss >> x >> y) coords.push_back(x), coords.push_back(y);
            
    //         polygon_2d new_poly;
    //         std::set<std::array<size_t,2>> new_edges;
    //         for (size_t i = 0; i < coords.size(); i += 2) {
    //             point_type p(coords[i], coords[i + 1]);
    //             points.push_back(p);
    //             new_poly.m_member_nodes.push_back(new_node_id);
    //             if (i > 0) {
    //                 new_edges.insert({new_node_id - 1, new_node_id});
    //             }
    //             new_node_id++;
    //         }
    //         // fermer le polygone
    //         if(new_poly.m_member_nodes.size() > 2){
    //             new_edges.insert({new_poly.m_member_nodes.back(), new_poly.m_member_nodes[0]});
    //         }
    //         new_poly.m_member_edges = new_edges;
    //         polygons.push_back(new_poly);
    //     }
        
    //     // Recréer les facets
    //     facets.clear();
    //     for (auto& poly : polygons) {
    //         for (auto& e : poly.m_member_edges) {
    //             facets.push_back(e);
    //         }
    //     }
    // }
    
    
};

#endif /* fitted_geometry_builder_hpp */



