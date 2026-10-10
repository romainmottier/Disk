/*
 * Structured meshes of the unit square and of the unit cube:
 *  - 2D: Cartesian grid (Section 5), grid of right triangles, hexagonal
 *    meshes of the FVCA5 benchmark (polygons);
 *  - 3D: Cartesian grid of hexahedra, Kuhn subdivision of the cubes into six
 *    tetrahedra, prisms extruded from the hexagonal meshes (polyhedra).
 * All the boundary faces get the boundary id 0, except on the Cartesian
 * grids of the square (0 bottom, 1 right, 2 top, 3 left).
 */

#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <vector>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/geometry/geometry.hpp"

namespace unfitted_meshes {

template<typename T>
void sort_unique(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

/* Structured nx x ny grid of [x0,x1]x[y0,y1] (the meshes of Section 5). Cells
 * are numbered row by row as in ProtoN (cell (i,j) has offset j*nx + i) and
 * their vertices are stored in tensor order {BL, BR, TL, TR}, like the other
 * cartesian meshes of diskpp. Boundary ids: 0 bottom, 1 right, 2 top, 3 left. */
template<typename T>
void make_cartesian_grid(disk::cartesian_mesh<T,2>& msh, size_t nx, size_t ny,
                         T x0 = 0.0, T x1 = 1.0, T y0 = 0.0, T y1 = 1.0)
{
    using mesh_type = disk::cartesian_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;

    if (nx == 0 or ny == 0)
        throw std::invalid_argument("make_cartesian_grid: empty grid");

    auto storage = msh.backend_storage();

    auto hx = (x1 - x0)/nx;
    auto hy = (y1 - y0)/ny;
    auto ptid = [&](size_t i, size_t j) { return disk::point_identifier<2>(j*(nx+1) + i); };

    for (size_t j = 0; j <= ny; j++) {
        for (size_t i = 0; i <= nx; i++) {
            storage->points.push_back( point_type(x0 + i*hx, y0 + j*hy) );
            storage->nodes.push_back( node_type({ptid(i,j)}) );
        }
    }

    using ne_pair = std::pair<edge_type, disk::boundary_descriptor>;
    std::vector<ne_pair> edges;
    for (size_t j = 0; j <= ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            disk::boundary_descriptor bd;
            if (j == 0)  bd = disk::boundary_descriptor(0, true);
            if (j == ny) bd = disk::boundary_descriptor(2, true);
            edges.push_back( std::make_pair(edge_type({ptid(i,j), ptid(i+1,j)}), bd) );
        }
    }
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i <= nx; i++) {
            disk::boundary_descriptor bd;
            if (i == 0)  bd = disk::boundary_descriptor(3, true);
            if (i == nx) bd = disk::boundary_descriptor(1, true);
            edges.push_back( std::make_pair(edge_type({ptid(i,j), ptid(i,j+1)}), bd) );
        }
    }
    std::sort(edges.begin(), edges.end(),
        [](const ne_pair& a, const ne_pair& b) { return a.first < b.first; });
    for (auto& e : edges) {
        storage->edges.push_back(e.first);
        storage->boundary_info.push_back(e.second);
    }

    for (size_t j = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++)
            storage->surfaces.push_back( surface_type({ptid(i,j), ptid(i+1,j), ptid(i,j+1), ptid(i+1,j+1)}) );
    std::sort(storage->surfaces.begin(), storage->surfaces.end());

    storage->subdomain_info.resize( storage->surfaces.size() );
}

/* n x n squares of the unit square, each split into two triangles along the
 * same diagonal */
template<typename T>
void make_triangle_grid(disk::simplicial_mesh<T,2>& msh, size_t n)
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;

    auto storage = msh.backend_storage();
    T h = 1.0/n;
    auto id = [&](size_t i, size_t j) { return disk::point_identifier<2>(j*(n+1) + i); };
    for (size_t j = 0; j <= n; j++)
        for (size_t i = 0; i <= n; i++) {
            storage->points.push_back(point_type(i*h, j*h));
            storage->nodes.push_back(node_type({id(i,j)}));
        }

    std::vector<edge_type> edges;
    std::vector<surface_type> surfaces;
    for (size_t j = 0; j < n; j++)
        for (size_t i = 0; i < n; i++) {
            auto a = id(i,j), b = id(i+1,j), c = id(i,j+1), d = id(i+1,j+1);
            surfaces.push_back(surface_type({a, b, d}));
            surfaces.push_back(surface_type({a, d, c}));
            for (auto e : {edge_type({a,b}), edge_type({b,d}), edge_type({a,d}), edge_type({d,c}), edge_type({a,c})})
                edges.push_back(e);
        }
    sort_unique(edges);
    sort_unique(surfaces);

    auto on_boundary = [&](size_t p, size_t q) {
        size_t ip = p % (n+1), jp = p / (n+1), iq = q % (n+1), jq = q / (n+1);
        return (ip == iq && (ip == 0 || ip == n)) || (jp == jq && (jp == 0 || jp == n));
    };
    storage->edges = edges;
    storage->boundary_info.resize(edges.size());
    for (size_t e = 0; e < edges.size(); e++) {
        auto pids = edges[e].point_ids();
        if (on_boundary(pids[0], pids[1]))
            storage->boundary_info[e] = disk::boundary_descriptor(0, true);
    }
    storage->surfaces = surfaces;
    storage->subdomain_info.resize(surfaces.size());
}

/* Hexagonal meshes of the FVCA5 benchmark (level >= 1), see hexagon_mesh.cpp
 * (the mesher of diskpp can only be included in one translation unit) */
void make_hexagon_mesh(disk::generic_mesh<double,2>& msh, size_t level);

/* n x n x n cubes of the unit cube; the vertices of the cells are stored in
 * tensor order, as by the cartesian loader of diskpp */
template<typename T>
void make_hexahedral_grid(disk::cartesian_mesh<T,3>& msh, size_t n)
{
    using mesh_type = disk::cartesian_mesh<T,3>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;
    using volume_type = typename mesh_type::volume_type;

    auto storage = msh.backend_storage();
    T h = 1.0/n;
    auto id = [&](size_t i, size_t j, size_t k) { return size_t((k*(n+1) + j)*(n+1) + i); };
    for (size_t k = 0; k <= n; k++)
        for (size_t j = 0; j <= n; j++)
            for (size_t i = 0; i <= n; i++) {
                storage->points.push_back(point_type(i*h, j*h, k*h));
                storage->nodes.push_back(node_type({disk::point_identifier<3>(id(i,j,k))}));
            }

    std::vector<edge_type> edges;
    std::vector<surface_type> surfaces;
    std::vector<volume_type> volumes;
    for (size_t k = 0; k < n; k++)
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < n; i++) {
                size_t p0 = id(i,j,k),   p1 = id(i+1,j,k),   p2 = id(i,j+1,k),   p3 = id(i+1,j+1,k);
                size_t p4 = id(i,j,k+1), p5 = id(i+1,j,k+1), p6 = id(i,j+1,k+1), p7 = id(i+1,j+1,k+1);
                for (auto e : {edge_type({p0, p1}), edge_type({p0, p2}), edge_type({p0, p4}), edge_type({p1, p3}),
                               edge_type({p1, p5}), edge_type({p2, p3}), edge_type({p2, p6}), edge_type({p3, p7}),
                               edge_type({p4, p5}), edge_type({p4, p6}), edge_type({p5, p7}), edge_type({p6, p7})})
                    edges.push_back(e);
                surfaces.push_back(surface_type({p0, p2, p6, p4}));
                surfaces.push_back(surface_type({p1, p3, p7, p5}));
                surfaces.push_back(surface_type({p0, p1, p3, p2}));
                surfaces.push_back(surface_type({p4, p5, p7, p6}));
                surfaces.push_back(surface_type({p0, p4, p5, p1}));
                surfaces.push_back(surface_type({p2, p6, p7, p3}));
                volumes.push_back(volume_type({p0, p1, p2, p3, p4, p5, p6, p7}));
            }
    sort_unique(edges);
    sort_unique(surfaces);
    std::sort(volumes.begin(), volumes.end());

    auto ijk = [&](size_t p) { return std::array<size_t,3>{p % (n+1), (p / (n+1)) % (n+1), p / ((n+1)*(n+1))}; };
    storage->edges = edges;
    storage->surfaces = surfaces;
    storage->volumes = volumes;
    storage->boundary_info.resize(surfaces.size());
    for (size_t s = 0; s < surfaces.size(); s++) {
        auto pids = surfaces[s].point_ids();
        for (size_t d = 0; d < 3; d++) {
            bool lo = true, hi = true;
            for (auto p : pids) {
                lo = lo && (ijk(p)[d] == 0);
                hi = hi && (ijk(p)[d] == n);
            }
            if (lo || hi)
                storage->boundary_info[s] = disk::boundary_descriptor(0, true);
        }
    }
    storage->subdomain_info.resize(volumes.size());
}

/* n x n x n cubes of the unit cube, each split into the six tetrahedra of
 * the Kuhn subdivision (along the diagonal (0,0,0)-(1,1,1)) */
template<typename T>
void make_tetrahedral_grid(disk::simplicial_mesh<T,3>& msh, size_t n)
{
    using mesh_type = disk::simplicial_mesh<T,3>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;
    using volume_type = typename mesh_type::volume_type;
    using pid = disk::point_identifier<3>;

    auto storage = msh.backend_storage();
    T h = 1.0/n;
    auto id = [&](size_t i, size_t j, size_t k) { return size_t((k*(n+1) + j)*(n+1) + i); };
    for (size_t k = 0; k <= n; k++)
        for (size_t j = 0; j <= n; j++)
            for (size_t i = 0; i <= n; i++) {
                storage->points.push_back(point_type(i*h, j*h, k*h));
                storage->nodes.push_back(node_type({pid(id(i,j,k))}));
            }

    const std::array<std::array<size_t,3>, 6> perms = {{ {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} }};
    std::vector<edge_type> edges;
    std::vector<surface_type> surfaces;
    std::vector<volume_type> volumes;
    for (size_t k = 0; k < n; k++)
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < n; i++)
                for (auto& pm : perms) {
                    std::array<size_t,3> c = {i, j, k};
                    std::array<pid,4> v;
                    v[0] = pid(id(c[0], c[1], c[2]));
                    for (size_t s = 0; s < 3; s++) {
                        c[pm[s]]++;
                        v[s+1] = pid(id(c[0], c[1], c[2]));
                    }
                    volumes.push_back(volume_type({v[0], v[1], v[2], v[3]}));
                    for (size_t a = 0; a < 4; a++)
                        for (size_t b = a+1; b < 4; b++)
                            edges.push_back(edge_type({v[a], v[b]}));
                    surfaces.push_back(surface_type({v[1], v[2], v[3]}));
                    surfaces.push_back(surface_type({v[0], v[2], v[3]}));
                    surfaces.push_back(surface_type({v[0], v[1], v[3]}));
                    surfaces.push_back(surface_type({v[0], v[1], v[2]}));
                }
    sort_unique(edges);
    sort_unique(surfaces);
    std::sort(volumes.begin(), volumes.end());

    auto ijk = [&](size_t p) { return std::array<size_t,3>{p % (n+1), (p / (n+1)) % (n+1), p / ((n+1)*(n+1))}; };
    storage->edges = edges;
    storage->surfaces = surfaces;
    storage->volumes = volumes;
    storage->boundary_info.resize(surfaces.size());
    for (size_t s = 0; s < surfaces.size(); s++) {
        auto pids = surfaces[s].point_ids();
        for (size_t d = 0; d < 3; d++) {
            bool lo = true, hi = true;
            for (auto p : pids) {
                lo = lo && (ijk(p)[d] == 0);
                hi = hi && (ijk(p)[d] == n);
            }
            if (lo || hi)
                storage->boundary_info[s] = disk::boundary_descriptor(0, true);
        }
    }
    storage->subdomain_info.resize(volumes.size());
}

/* Prisms of the unit cube: the polygons of a 2D mesh of the unit square,
 * extruded in nz layers */
template<typename T>
void make_prism_mesh(disk::generic_mesh<T,3>& msh, const disk::generic_mesh<T,2>& m2, size_t nz)
{
    using mesh_type = disk::generic_mesh<T,3>;
    using point_type = typename mesh_type::point_type;
    using node_type = typename mesh_type::node_type;
    using edge_type = typename mesh_type::edge_type;
    using surface_type = typename mesh_type::surface_type;
    using volume_type = typename mesh_type::volume_type;
    using edge_id = typename edge_type::id_type;
    using surface_id = typename surface_type::id_type;

    auto storage = msh.backend_storage();
    size_t n2 = m2.points_size();
    T hz = 1.0/nz;
    std::vector<typename disk::generic_mesh<T,2>::point_type> pts2(m2.points_begin(), m2.points_end());
    for (size_t k = 0; k <= nz; k++)
        for (size_t p = 0; p < n2; p++) {
            storage->points.push_back(point_type(pts2[p].x(), pts2[p].y(), k*hz));
            storage->nodes.push_back(node_type(disk::point_identifier<3>(k*n2 + p)));
        }
    auto pid = [&](size_t k, size_t p) { return k*n2 + p; };

    /* the polygons of the 2D mesh (vertices counterclockwise) */
    std::vector<std::vector<size_t>> polys;
    for (auto& cl : m2) {
        auto ids = cl.point_ids();
        std::vector<size_t> pv(ids.begin(), ids.end());
        T area = 0;
        for (size_t i = 0; i < pv.size(); i++) {
            const auto& a = pts2[pv[i]];
            const auto& b = pts2[pv[(i+1) % pv.size()]];
            area += a.x()*b.y() - a.y()*b.x();
        }
        if (area < 0)
            std::reverse(pv.begin(), pv.end());
        polys.push_back(pv);
    }

    /* edges */
    std::vector<edge_type> edges;
    for (size_t k = 0; k <= nz; k++)
        for (auto& pv : polys)
            for (size_t i = 0; i < pv.size(); i++)
                edges.push_back(edge_type(typename edge_type::sub_id_type(pid(k, pv[i])),
                                          typename edge_type::sub_id_type(pid(k, pv[(i+1) % pv.size()]))));
    for (size_t k = 0; k < nz; k++)
        for (size_t p = 0; p < n2; p++)
            edges.push_back(edge_type(typename edge_type::sub_id_type(pid(k, p)),
                                      typename edge_type::sub_id_type(pid(k+1, p))));
    sort_unique(edges);
    auto edge_of = [&](size_t a, size_t b) {
        auto e = edge_type(typename edge_type::sub_id_type(a), typename edge_type::sub_id_type(b));
        auto it = std::lower_bound(edges.begin(), edges.end(), e);
        if (it == edges.end() || !(*it == e))
            throw std::logic_error("make_prism_mesh: edge not found");
        return edge_id(std::distance(edges.begin(), it));
    };

    /* faces: (edges, points) */
    struct face_def { std::vector<edge_id> edges; std::vector<size_t> points; };
    std::vector<face_def> fdefs;
    std::map<std::vector<size_t>, size_t> face_index; // sorted points -> face
    auto add_face = [&](const std::vector<size_t>& fp) {
        auto key = fp;
        std::sort(key.begin(), key.end());
        auto it = face_index.find(key);
        if (it != face_index.end())
            return it->second;
        face_def fd;
        fd.points = fp;
        for (size_t i = 0; i < fp.size(); i++)
            fd.edges.push_back(edge_of(fp[i], fp[(i+1) % fp.size()]));
        fdefs.push_back(fd);
        face_index[key] = fdefs.size() - 1;
        return fdefs.size() - 1;
    };
    std::vector<std::vector<size_t>> vol_faces, vol_points;
    for (size_t k = 0; k < nz; k++)
        for (auto& pv : polys) {
            std::vector<size_t> vf, vp;
            std::vector<size_t> bottom, top;
            for (auto p : pv) {
                bottom.push_back(pid(k, p));
                top.push_back(pid(k+1, p));
            }
            std::reverse(bottom.begin(), bottom.end()); /* outward orientation */
            vf.push_back(add_face(bottom));
            vf.push_back(add_face(top));
            for (size_t i = 0; i < pv.size(); i++) {
                size_t a = pv[i], b = pv[(i+1) % pv.size()];
                vf.push_back(add_face({pid(k, a), pid(k, b), pid(k+1, b), pid(k+1, a)}));
            }
            vp.insert(vp.end(), bottom.begin(), bottom.end());
            vp.insert(vp.end(), top.begin(), top.end());
            vol_faces.push_back(vf);
            vol_points.push_back(vp);
        }

    /* sort the faces by their edges (lookup), renumber the faces of the volumes */
    std::vector<size_t> order(fdefs.size());
    for (size_t i = 0; i < order.size(); i++)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return fdefs[a].edges < fdefs[b].edges; });
    std::vector<size_t> new_index(fdefs.size());
    std::vector<surface_type> surfaces;
    for (size_t i = 0; i < order.size(); i++) {
        new_index[order[i]] = i;
        const auto& fd = fdefs[order[i]];
        surface_type s(fd.edges);
        std::vector<disk::point_identifier<3>> sp;
        for (auto p : fd.points)
            sp.push_back(disk::point_identifier<3>(p));
        s.set_point_ids(sp.begin(), sp.end());
        surfaces.push_back(s);
    }

    std::vector<std::pair<volume_type, size_t>> volumes;
    std::vector<size_t> owners(fdefs.size(), 0);
    for (size_t v = 0; v < vol_faces.size(); v++) {
        std::vector<surface_id> vf;
        for (auto f : vol_faces[v]) {
            vf.push_back(surface_id(new_index[f]));
            owners[new_index[f]]++;
        }
        volume_type vol(vf);
        std::vector<disk::point_identifier<3>> vp;
        for (auto p : vol_points[v])
            vp.push_back(disk::point_identifier<3>(p));
        vol.set_point_ids(vp.begin(), vp.end());
        volumes.push_back(std::make_pair(vol, v));
    }
    std::sort(volumes.begin(), volumes.end(),
              [](const std::pair<volume_type, size_t>& a, const std::pair<volume_type, size_t>& b) { return a.first < b.first; });

    storage->edges = edges;
    storage->surfaces = surfaces;
    for (auto& v : volumes)
        storage->volumes.push_back(v.first);
    storage->boundary_info.resize(surfaces.size());
    for (size_t s = 0; s < surfaces.size(); s++)
        if (owners[s] == 1)
            storage->boundary_info[s] = disk::boundary_descriptor(0, true);
    storage->subdomain_info.resize(storage->volumes.size());
}

} // namespace unfitted_meshes
