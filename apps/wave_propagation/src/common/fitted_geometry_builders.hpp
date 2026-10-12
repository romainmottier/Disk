//
//  fitted_geometry_builder.hpp
//
//  Created by Omar Durán
//  Contributor: Romain Mottier


#pragma once
#ifndef fitted_geometry_builder_hpp
#define fitted_geometry_builder_hpp

#include <queue>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cassert>
#include <thread>
#include <set>
#include "diskpp/geometry/geometry.hpp"
#include <unordered_map>

template<typename MESH>
class fitted_geometry_builder {
    
protected:
    
    size_t m_dimension = 0;
    
    // orient an edge from its smallest to its largest node index
    static void validate_edge(std::array<size_t, 2> & edge){
        assert(edge[0] != edge[1]);
        if (edge[0] > edge[1]){
            std::swap(edge[0], edge[1]);
        }
    }
    
    // Moves points and nodes to the mesh storage, builds the sorted edges from the facets, flags the
    // boundary edges and returns the surfaces in polygon order (sorted and stored by the caller).
    // Returns false if an edge cannot be found.
    template<typename Polygon>
    static bool fill_mesh_storage(MESH& msh, std::vector<typename MESH::point_type>& points, std::vector<typename MESH::node_type>& vertices,
                                  const std::vector<std::array<size_t, 2>>& facets, const std::vector<std::array<size_t, 2>>& boundary_edges,
                                  const std::vector<Polygon>& polygons, std::vector<typename MESH::surface_type>& surfaces){
        
        typedef typename MESH::node_type        node_type;
        typedef typename MESH::edge_type        edge_type;
        typedef typename MESH::surface_type     surface_type;
        
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
        
        // reset: a mesh object can receive a second builder (e.g. local refinement), stale flags would remain
        storage->boundary_info.assign(edges.size(), disk::boundary_descriptor{});
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
                return false;
            }
            disk::boundary_descriptor bi{0, true};
            storage->boundary_info.at(position.second) = bi;
        }
        
        storage->edges = std::move(edges);
        
        surfaces.reserve(polygons.size());
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
                    return false;
                }
                surface_edges.push_back(edge_id.second);
            }
            auto surface = surface_type(surface_edges);
            surface.set_point_ids(p.m_member_nodes.begin(), p.m_member_nodes.end());
            surfaces.push_back(surface);
        }
        return true;
    }
    
public:
    
    fitted_geometry_builder() {
        
    }
    
    // get the geometry dimension
    size_t get_dimension() {
        return m_dimension;
    }
    
    // build the mesh
    virtual bool build_mesh()  = 0;
    
    // move generated mesh data to an external mesh storage
    virtual void move_to_mesh_storage(MESH& msh) = 0;
    
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
    using fitted_geometry_builder<mesh_type>::validate_edge;
    
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
    
void refine_cells(const std::vector<size_t>& cell_indices,
                  int refinement_level) {

    if (refinement_level <= 0) return;

    for (size_t idx : cell_indices) {
        assert(idx < polygons.size());
        polygons[idx].m_to_refine = true;
    }

    for (int pass = 0; pass < refinement_level; ++pass) {

        std::vector<size_t> targets;
        for (size_t i = 0; i < polygons.size(); ++i)
            if (polygons[i].m_to_refine)
                targets.push_back(i);

        // std::cout << "[refine_cells] pass " << pass
        //           << " | polygons=" << polygons.size()
        //           << " | targets=" << targets.size() << std::endl;

        // Extract corners BEFORE any modification
        std::vector<std::array<size_t,4>> all_corners(targets.size());
        for (size_t i = 0; i < targets.size(); ++i)
            all_corners[i] = extract_quad_corners(targets[i]);

        // ----------------------------------------------------------------
        // PHASE 1: compute all midpoints and update neighbor edges FIRST,
        // before any polygon is erased. This ensures cell indices are
        // still valid when we search for neighbors.
        // ----------------------------------------------------------------
        // Store the midpoints for each target so phase 2 can use them
        struct CellSplit {
            size_t n0, n1, n2, n3;
            size_t m01, m12, m23, m03, c;
        };
        std::vector<CellSplit> splits(targets.size());

        for (size_t i = 0; i < targets.size(); ++i) {
            size_t cell_index = targets[i];
            const auto& cors = all_corners[i];

            const size_t n0 = cors[0];
            const size_t n1 = cors[1];
            const size_t n2 = cors[2];
            const size_t n3 = cors[3];

            auto get_or_create_midpoint = [&](size_t na, size_t nb) -> size_t {
                point_type mid(
                    (points[na].x() + points[nb].x()) * T(0.5),
                    (points[na].y() + points[nb].y()) * T(0.5)
                );
                constexpr T tol = T(1e-14);
                for (size_t k = 0; k < points.size(); ++k) {
                    T dx = points[k].x() - mid.x();
                    T dy = points[k].y() - mid.y();
                    if (dx*dx + dy*dy < tol*tol) return k;
                }
                size_t new_id = points.size();
                points.push_back(mid);
                vertices.push_back(node_type(disk::point_identifier<2>(new_id)));
                return new_id;
            };

            size_t m01 = get_or_create_midpoint(n0, n1);
            size_t m12 = get_or_create_midpoint(n1, n2);
            size_t m23 = get_or_create_midpoint(n2, n3);
            size_t m03 = get_or_create_midpoint(n0, n3);
            size_t c;
            {
                point_type center(
                    (points[n0].x()+points[n1].x()+points[n2].x()+points[n3].x())*T(0.25),
                    (points[n0].y()+points[n1].y()+points[n2].y()+points[n3].y())*T(0.25)
                );
                constexpr T tol = T(1e-14);
                c = points.size();
                for (size_t k = 0; k < points.size(); ++k) {
                    T dx = points[k].x() - center.x();
                    T dy = points[k].y() - center.y();
                    if (dx*dx + dy*dy < tol*tol) { c = k; break; }
                }
                if (c == points.size()) {
                    points.push_back(center);
                    vertices.push_back(node_type(disk::point_identifier<2>(c)));
                }
            }

            splits[i] = {n0, n1, n2, n3, m01, m12, m23, m03, c};

            // Update boundary edges
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

            // Update neighbor m_member_edges — cell indices still valid
            // because no polygon has been erased yet in this pass
            struct EdgeSplit { size_t a, b, mid; };
            std::array<EdgeSplit,4> esplits = {{
                {n0,n1,m01},{n1,n2,m12},{n2,n3,m23},{n0,n3,m03}
            }};
            for (auto& [ea, eb, emid] : esplits) {
                std::array<size_t,2> parent_edge = {ea, eb};
                validate_edge(parent_edge);
                std::array<size_t,2> h1 = {ea, emid}; validate_edge(h1);
                std::array<size_t,2> h2 = {emid, eb}; validate_edge(h2);
                for (size_t ci = 0; ci < polygons.size(); ++ci) {
                    if (ci == cell_index) continue;
                    polygon_2d& nb = polygons[ci];
                    if (!nb.m_member_edges.count(parent_edge)) continue;
                    nb.m_member_edges.erase(parent_edge);
                    nb.m_member_edges.insert(h1);
                    nb.m_member_edges.insert(h2);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------
        // PHASE 2: build children and erase parents, descending order
        // ----------------------------------------------------------------
        // Clear flags before erasing (indices will shift)
        for (size_t idx : targets)
            polygons[idx].m_to_refine = false;

        // Process descending to keep indices stable
        // Re-sort targets descending
        std::vector<size_t> sorted_targets = targets;
        std::sort(sorted_targets.begin(), sorted_targets.end(), std::greater<size_t>());

        // Map original index -> split data
        std::map<size_t, size_t> idx_to_split;
        for (size_t i = 0; i < targets.size(); ++i)
            idx_to_split[targets[i]] = i;

        for (size_t orig_idx : sorted_targets) {
            size_t i = idx_to_split[orig_idx];
            const auto& sp = splits[i];
            const int parent_level = polygons[orig_idx].m_refinement_level;
            const int parent_mat   = polygons[orig_idx].m_material;
            const bool parent_elas = polygons[orig_idx].m_elastic_material;

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
                add_edge(a, b); add_edge(b, d_);
                add_edge(d_, e_); add_edge(e_, a);
                return q;
            };

            polygon_2d q0 = make_quad(sp.n0, sp.m01, sp.c,   sp.m03);
            polygon_2d q1 = make_quad(sp.m01, sp.n1, sp.m12, sp.c  );
            polygon_2d q2 = make_quad(sp.c,  sp.m12, sp.n2,  sp.m23);
            polygon_2d q3 = make_quad(sp.m03, sp.c,  sp.m23, sp.n3 );

            polygons.erase(polygons.begin() + orig_idx);
            polygons.push_back(q0);
            polygons.push_back(q1);
            polygons.push_back(q2);
            polygons.push_back(q3);

            // Mark children for next pass
            size_t new_end = polygons.size();
            for (size_t ci = new_end - 4; ci < new_end; ++ci)
                polygons[ci].m_to_refine = true;
        }

        // ----------------------------------------------------------------
        // Rebuild m_member_nodes and facets after all refinements
        // ----------------------------------------------------------------
        for (auto& poly : polygons) {
            const auto& medges = poly.m_member_edges;
            if (medges.empty()) continue;
            std::unordered_map<size_t, std::vector<size_t>> adj;
            for (const auto& e : medges) {
                adj[e[0]].push_back(e[1]);
                adj[e[1]].push_back(e[0]);
            }
            size_t start = std::numeric_limits<size_t>::max();
            for (const auto& e : medges) {
                if (e[0] < start) start = e[0];
                if (e[1] < start) start = e[1];
            }
            std::vector<size_t> ordered;
            ordered.reserve(medges.size());
            size_t prev = std::numeric_limits<size_t>::max();
            size_t cur  = start;
            for (size_t step = 0; step < medges.size(); ++step) {
                ordered.push_back(cur);
                const auto& nbrs = adj[cur];
                size_t next = std::numeric_limits<size_t>::max();
                for (size_t nb : nbrs)
                    if (nb != prev) { next = nb; break; }
                if (next == std::numeric_limits<size_t>::max()) break;
                prev = cur; cur = next;
            }
            poly.m_member_nodes = ordered;
        }

        {
            std::set<std::array<size_t,2>> facet_set;
            for (auto& poly : polygons)
                for (auto& e : poly.m_member_edges)
                    facet_set.insert(e);
            facets.assign(facet_set.begin(), facet_set.end());
        }
    }

    for (auto& poly : polygons)
        poly.m_to_refine = false;
}


std::array<size_t, 4> extract_quad_corners(size_t cell_index) const {

    const auto& nodes = polygons[cell_index].m_member_nodes;

    T cx = T(0), cy = T(0);
    for (size_t nid : nodes) {
        cx += points[nid].x();
        cy += points[nid].y();
    }
    cx /= T(nodes.size());
    cy /= T(nodes.size());

    // The 4 corners are the farthest node in each quadrant.
    // Hanging nodes are always closer to the centroid than true corners.
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
        int q  = (dx >= T(0) ? 1 : 0) + (dy >= T(0) ? 2 : 0);
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


// ------------------------------------------------------------
// Barycentre d'une cellule du builder
// (opère sur polygons/points, avant move_to_mesh_storage)
// ------------------------------------------------------------
point_type cell_barycenter(size_t cell_id) const
{
    const auto& nodes = polygons[cell_id].m_member_nodes;
    T cx = T(0), cy = T(0);
    for (size_t nid : nodes) {
        cx += points[nid].x();
        cy += points[nid].y();
    }
    cx /= T(nodes.size());
    cy /= T(nodes.size());
    return point_type(cx, cy);
}
 
 
// ------------------------------------------------------------
// Voisins d'une cellule du builder
//
// Deux cellules sont voisines si elles partagent au moins
// un noeud ET que ce noeud est strictement sur une arête
// de l'autre cellule (gère le cas grosse/petite cellule).
//
// Inspiré de find_cells() : on parcourt les point_ids de
// chaque cellule, comme dans le style diskpp existant.
// ------------------------------------------------------------
std::vector<size_t> get_neighbors(size_t cell_id) const
{
    // Collecter les noeuds de la cellule courante
    const std::set<size_t> cell_nodes(
        polygons[cell_id].m_member_nodes.begin(),
        polygons[cell_id].m_member_nodes.end()
    );
 
    std::vector<size_t> neighbors;
 
    for (size_t j = 0; j < polygons.size(); ++j) {
        if (j == cell_id) continue;
 
        // Cas 1 : partage une arête complète via m_member_edges
        bool shared_edge = false;
        for (const auto& e : polygons[cell_id].m_member_edges) {
            if (polygons[j].m_member_edges.count(e)) {
                shared_edge = true;
                break;
            }
        }
        if (shared_edge) {
            neighbors.push_back(j);
            continue;
        }
 
        // Cas 2 : grosse cellule voisine d'une petite
        // Un noeud de j est strictement sur une arête de cell_id
        // (hanging node) — même logique que find_cells pour les pt_ids
        bool hanging = false;
        for (const auto& e_big : polygons[cell_id].m_member_edges) {
            const point_type& pa = points[e_big[0]];
            const point_type& pb = points[e_big[1]];
            T dx = pb.x() - pa.x();
            T dy = pb.y() - pa.y();
            T len2 = dx*dx + dy*dy;
            if (len2 < T(1e-28)) continue;
 
            for (size_t nid : polygons[j].m_member_nodes) {
                if (cell_nodes.count(nid)) continue; // déjà coin commun
                const point_type& p = points[nid];
                T t = ((p.x()-pa.x())*dx + (p.y()-pa.y())*dy) / len2;
                if (t < T(1e-12) || t > T(1)-T(1e-12)) continue;
                T px = pa.x() + t*dx - p.x();
                T py = pa.y() + t*dy - p.y();
                if (px*px + py*py < T(1e-24)) {
                    hanging = true;
                    break;
                }
            }
            if (hanging) break;
        }
        if (hanging) neighbors.push_back(j);
    }
 
    return neighbors;
}
 
 
// ------------------------------------------------------------
// Raffinement avec couches de protection (BFS)
//
// Garantit : max 1 hanging node par transition de niveau.
// Les cellules hors de portée du BFS restent au niveau 0.
//
// Étape 1 — BFS depuis la zone fine :
//   required[i] = L pour la zone fine
//   required[voisin] = max(required[voisin], required[cur]-1)
//   propagation s'arrête quand required[cur] <= 1
//   → les cellules non atteintes restent à 0 (coarse)
//
// Étape 2 — L passes de raffinement, r = 1..L :
//   raffiner d'un niveau les cellules avec required >= r
//   resynchroniser required après chaque appel à refine_cells
// ------------------------------------------------------------
void refine_with_protection(
    std::function<bool(const point_type&)> is_fine_zone,
    int L)
{
    if (L <= 0) return;
 
    // ----------------------------------------------------------
    // Étape 1 — BFS
    // ----------------------------------------------------------
    std::vector<int> required(polygons.size(), 0);
 
    std::queue<size_t> q;
    for (size_t i = 0; i < polygons.size(); ++i) {
        if (is_fine_zone(cell_barycenter(i))) {
            required[i] = L;
            q.push(i);
        }
    }
 
    while (!q.empty()) {
        size_t cur = q.front(); q.pop();
        int lvl = required[cur];
        if (lvl <= 1) continue;   // voisins peuvent rester à 0
 
        for (size_t nb : get_neighbors(cur)) {
            if (required[nb] < lvl - 1) {
                required[nb] = lvl - 1;
                q.push(nb);
            }
        }
    }
 
    // ----------------------------------------------------------
    // Étape 2 — L passes de raffinement
    // ----------------------------------------------------------
    for (int r = 1; r <= L; ++r) {
 
        // Collecter les cibles : required >= r et pas encore au niveau r
        std::vector<size_t> targets;
        std::vector<int>    target_req;
        for (size_t i = 0; i < polygons.size(); ++i) {
            if (required[i] >= r && polygons[i].m_refinement_level < r) {
                targets.push_back(i);
                target_req.push_back(required[i]);
            }
        }
        if (targets.empty()) continue;
 
        // Snapshot du required pour les survivants avant raffinement
        // refine_cells va supprimer les cibles (ordre décroissant)
        // et pousser 4 enfants en fin de vecteur pour chacune
        size_t n_before = polygons.size();
        std::set<size_t> target_set(targets.begin(), targets.end());
 
        // Snapshot complet : required de TOUTES les cellules actuelles
        std::vector<int> req_snapshot = required;
 
        refine_cells(targets, 1);
 
        // ----------------------------------------------------------
        // Resynchronisation de required
        //
        // Après refine_cells(targets, 1) :
        //   - les n_before cellules originales moins les |targets|
        //     cibles supprimées donnent (n_before - |targets|)
        //     survivants, conservés en tête dans le même ordre relatif
        //   - 4*|targets| nouveaux enfants sont ajoutés en queue,
        //     dans l'ordre de targets (croissant)
        // ----------------------------------------------------------
        std::vector<int> new_required;
        new_required.reserve(polygons.size());
 
        // Survivants : cellules originales non supprimées
        for (size_t i = 0; i < n_before; ++i) {
            if (!target_set.count(i))
                new_required.push_back(req_snapshot[i]);
        }
 
        // Enfants : 4 enfants par cible, héritent du required parent
        for (size_t k = 0; k < targets.size(); ++k) {
            int parent_req = target_req[k];
            for (int c = 0; c < 4; ++c)
                new_required.push_back(parent_req);
        }
 
        required = std::move(new_required);
    }
}


void move_to_mesh_storage(mesh_type& msh){
    
    std::vector<surface_type> surfaces;
    if (!fitted_geometry_builder<mesh_type>::fill_mesh_storage(msh, points, vertices, facets, boundary_edges, polygons, surfaces))
        return;
    
    std::sort(surfaces.begin(), surfaces.end());
    msh.backend_storage()->surfaces = std::move(surfaces);
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
};


template<typename T>
class polygon_2d_mesh_reader : 

public fitted_geometry_builder<disk::generic_mesh<T, 2>> {

    typedef disk::generic_mesh<T,2>                 mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;   
    using fitted_geometry_builder<mesh_type>::validate_edge;

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

    polygon_2d_mesh_reader() : fitted_geometry_builder<mesh_type>()
    {
        fitted_geometry_builder<mesh_type>::m_dimension = 2;
    }
    
    void set_poly_mesh_file(std::string mesh_file){
        poly_mesh_file = mesh_file;
    }
    
    // build the mesh
    bool build_mesh(){
        return read_poly_mesh(false);
    }
    
    // build the mesh of the basin/brain files, which also carry one material per polygon
    bool build_bassin(){
        return read_poly_mesh(true);
    }
    
    // Reads poly_mesh_file. With materials, the polygons are followed by 4 border lines and one
    // material per polygon, and the duplicated facets are kept (they always were for these files).
    bool read_poly_mesh(bool with_materials){
        
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

            if (with_materials) {
                // skip the 4 border lines
                for(size_t bord=0; bord < 4; bord++) {
                    std::getline(input, line);
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
        }
          
        // Duplicated facets are eliminated
        if (!with_materials) {
            std::sort( facets.begin(), facets.end() );
            facets.erase( std::unique( facets.begin(), facets.end() ), facets.end() );
        }

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
        
        std::vector<surface_type> surfaces;
        if (!fitted_geometry_builder<mesh_type>::fill_mesh_storage(msh, points, vertices, facets, boundary_edges, polygons, surfaces))
            return;
        
        std::vector<int> materials;
        materials.reserve( polygons.size() );
        for (auto& p : polygons) {
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
        msh.backend_storage()->surfaces = std::move(surfaces);
        
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
};

#endif /* fitted_geometry_builder_hpp */



