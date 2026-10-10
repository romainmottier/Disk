/*
 * Local structured refinement of the background meshes where the interface
 * is not resolved, i.e. where a cell is not cut once (the assumption of the
 * method, see disk::cut_once_violations). The marked cells are split as on
 * the next level of the mesh, until all the cells are cut once or have been
 * refined max_levels times:
 *  - polygons: triangles into 4 triangles, the other polygons into
 *    quadrangles joining the middles of their edges to their barycenter
 *    (quadrangles into 4 quadrangles, hexagons into 6 quadrangles);
 *  - tetrahedra into 8 tetrahedra (4 at the vertices, the inner octahedron
 *    split along its shortest diagonal);
 *  - prisms (hexahedra are prisms on quadrangles): base polygon split as
 *    above, height split in 2.
 * The refinement is balanced (two cells sharing a point differ by at most one
 * level); the coarser neighbors of a refined cell keep their faces, split at
 * the new points (hanging nodes), so that the refined mesh is a polytopal
 * mesh (disk::generic_mesh), handled by the HHO methods without any change.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/cut_mesh.hpp"
#include "meshes.hpp"

namespace unfitted_meshes {

enum class refinable_shape { polygon, tetrahedron, prism };

template<typename T, size_t DIM>
class refinable_mesh {

public:
    using point_type = disk::point<T,DIM>;

    struct cell {
        refinable_shape     shape;
        std::vector<size_t> v;      // polygon: vertices counterclockwise; tetrahedron: 4 vertices;
                                    // prism: base polygon, then top polygon (v[n+i] above v[i])
        size_t              level;  // number of refinements
    };

    std::vector<point_type> points;
    std::vector<cell>       cells;

private:
    std::map<std::pair<size_t,size_t>, size_t>                      m_middles;
    std::map<std::vector<size_t>, size_t>                           m_centers;
    std::map<std::vector<size_t>, std::vector<std::vector<size_t>>> m_face_children;  // 3D

    static std::vector<size_t>
    sorted(std::vector<size_t> v)
    {
        std::sort(v.begin(), v.end());
        return v;
    }

    /* middle of the segment [a, b], barycenter of a polygon (shared by all
     * the cells refined) */
    size_t
    middle(size_t a, size_t b)
    {
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = m_middles.find(key);
        if (it != m_middles.end())
            return it->second;
        points.push_back((points[a] + points[b]) / T(2));
        m_middles[key] = points.size() - 1;
        return points.size() - 1;
    }

    size_t
    center(const std::vector<size_t>& p)
    {
        auto key = sorted(p);
        auto it = m_centers.find(key);
        if (it != m_centers.end())
            return it->second;
        point_type c;
        for (auto i : p)
            c = c + points[i];
        points.push_back(c / T(p.size()));
        m_centers[key] = points.size() - 1;
        return points.size() - 1;
    }

    /* polygons of the refinement of a polygon */
    std::vector<std::vector<size_t>>
    split_polygon(const std::vector<size_t>& p)
    {
        size_t n = p.size();
        std::vector<size_t> m(n);
        for (size_t i = 0; i < n; i++)
            m[i] = middle(p[i], p[(i+1) % n]);
        if (n == 3)
            return {{p[0], m[0], m[2]}, {m[0], p[1], m[1]}, {m[2], m[1], p[2]}, {m[0], m[1], m[2]}};
        size_t c = center(p);
        std::vector<std::vector<size_t>> ret;
        for (size_t i = 0; i < n; i++)
            ret.push_back({p[i], m[i], c, m[(i + n - 1) % n]});
        return ret;
    }

    /* faces of a cell of a 3D mesh */
    static std::vector<std::vector<size_t>>
    cell_faces(const cell& c)
    {
        const auto& v = c.v;
        if (c.shape == refinable_shape::tetrahedron)
            return {{v[1], v[2], v[3]}, {v[0], v[2], v[3]}, {v[0], v[1], v[3]}, {v[0], v[1], v[2]}};
        size_t n = v.size() / 2;
        std::vector<std::vector<size_t>> ret;
        ret.push_back(std::vector<size_t>(v.begin(), v.begin() + n));
        ret.push_back(std::vector<size_t>(v.begin() + n, v.end()));
        for (size_t i = 0; i < n; i++)
            ret.push_back({v[i], v[(i+1) % n], v[n + (i+1) % n], v[n + i]});
        return ret;
    }

    /* the neighbors of a refined cell see its faces split */
    void
    split_face(const std::vector<size_t>& f)
    {
        auto key = sorted(f);
        if (m_face_children.find(key) == m_face_children.end())
            m_face_children[key] = split_polygon(f);
    }

    std::vector<cell>
    children(const cell& c)
    {
        std::vector<cell> ret;
        auto make = [&](refinable_shape s, std::vector<size_t> v) { ret.push_back(cell{s, std::move(v), c.level + 1}); };

        if (c.shape == refinable_shape::polygon) {
            for (auto& q : split_polygon(c.v))
                make(refinable_shape::polygon, q);
            return ret;
        }

        if (c.shape == refinable_shape::tetrahedron) {
            const auto& a = c.v;
            size_t m01 = middle(a[0], a[1]), m02 = middle(a[0], a[2]), m03 = middle(a[0], a[3]);
            size_t m12 = middle(a[1], a[2]), m13 = middle(a[1], a[3]), m23 = middle(a[2], a[3]);
            make(refinable_shape::tetrahedron, {a[0], m01, m02, m03});
            make(refinable_shape::tetrahedron, {m01, a[1], m12, m13});
            make(refinable_shape::tetrahedron, {m02, m12, a[2], m23});
            make(refinable_shape::tetrahedron, {m03, m13, m23, a[3]});
            /* inner octahedron: 4 tetrahedra around its shortest diagonal */
            std::array<std::array<size_t,2>, 3> diag = {{ {m01, m23}, {m02, m13}, {m03, m12} }};
            auto length = [&](const std::array<size_t,2>& d) { return (points[d[0]] - points[d[1]]).to_vector().norm(); };
            size_t best = 0;
            for (size_t d = 1; d < 3; d++)
                if (length(diag[d]) < length(diag[best]))
                    best = d;
            /* the other two diagonals: equator e0, e1, e2, e3 (e0, e2 and e1, e3 opposite) */
            const auto& d1 = diag[(best + 1) % 3];
            const auto& d2 = diag[(best + 2) % 3];
            std::array<size_t,4> e = {d1[0], d2[0], d1[1], d2[1]};
            for (size_t i = 0; i < 4; i++)
                make(refinable_shape::tetrahedron, {diag[best][0], diag[best][1], e[i], e[(i+1) % 4]});
            for (auto& f : cell_faces(c))
                split_face(f);
            return ret;
        }

        /* prism: base and top polygons split, height split in 2 */
        size_t n = c.v.size() / 2;
        std::vector<size_t> b(c.v.begin(), c.v.begin() + n), t(c.v.begin() + n, c.v.end());
        std::map<size_t, std::pair<size_t,size_t>> up;  // point of the base -> (middle layer, top)
        for (size_t i = 0; i < n; i++) {
            size_t j = (i+1) % n;
            up[b[i]] = std::make_pair(middle(b[i], t[i]), t[i]);
            up[middle(b[i], b[j])] = std::make_pair(center({b[i], b[j], t[j], t[i]}), middle(t[i], t[j]));
        }
        if (n > 3) {
            size_t cb = center(b), ct = center(t);
            up[cb] = std::make_pair(middle(cb, ct), ct);
        }
        for (auto& q : split_polygon(b)) {
            std::vector<size_t> lower = q, upper;
            for (auto p : q)
                lower.push_back(up.at(p).first);
            for (auto p : q)
                upper.push_back(up.at(p).first);
            for (auto p : q)
                upper.push_back(up.at(p).second);
            make(refinable_shape::prism, lower);
            make(refinable_shape::prism, upper);
        }
        for (auto& f : cell_faces(c))
            split_face(f);
        return ret;
    }

    /* the points of the segment [a, b] but b: a, then the middles inserted by
     * the refinements of the cells containing it */
    void
    edge_points(size_t a, size_t b, std::vector<size_t>& out) const
    {
        auto it = m_middles.find(std::make_pair(std::min(a, b), std::max(a, b)));
        if (it == m_middles.end()) {
            out.push_back(a);
            return;
        }
        edge_points(a, it->second, out);
        edge_points(it->second, b, out);
    }

    /* polygon with the hanging nodes of its edges */
    std::vector<size_t>
    with_hanging_nodes(const std::vector<size_t>& p) const
    {
        std::vector<size_t> ret;
        for (size_t i = 0; i < p.size(); i++)
            edge_points(p[i], p[(i+1) % p.size()], ret);
        return ret;
    }

    /* faces of the mesh covering the face f of a cell */
    void
    leaf_faces(const std::vector<size_t>& f, std::vector<std::vector<size_t>>& out) const
    {
        auto it = m_face_children.find(sorted(f));
        if (it == m_face_children.end()) {
            out.push_back(with_hanging_nodes(f));
            return;
        }
        for (auto& c : it->second)
            leaf_faces(c, out);
    }

public:
    /* refine the cells "marked" (indices in cells) */
    void
    refine(const std::vector<size_t>& marked)
    {
        std::vector<bool> is_marked(cells.size(), false);
        for (auto i : marked)
            is_marked[i] = true;
        std::vector<cell> new_cells;
        for (size_t i = 0; i < cells.size(); i++) {
            if (!is_marked[i]) {
                new_cells.push_back(cells[i]);
                continue;
            }
            auto ch = children(cells[i]);
            new_cells.insert(new_cells.end(), ch.begin(), ch.end());
        }
        cells = std::move(new_cells);
    }

    size_t
    max_level() const
    {
        size_t ret = 0;
        for (auto& c : cells)
            ret = std::max(ret, c.level);
        return ret;
    }

    /* The diskpp mesh of the cells (all the boundary faces get the id 0).
     * Returns, for each cell of msh, its index in cells. */
    std::vector<size_t>
    build(disk::generic_mesh<T,DIM>& msh) const
    {
        using mesh_type = disk::generic_mesh<T,DIM>;
        using node_type = typename mesh_type::node_type;
        using edge_type = typename mesh_type::edge_type;
        using surface_type = typename mesh_type::surface_type;
        using edge_id = typename edge_type::id_type;

        auto storage = msh.backend_storage();
        storage->points = points;
        for (size_t i = 0; i < points.size(); i++)
            storage->nodes.push_back(node_type(disk::point_identifier<DIM>(i)));

        /* boundary of the cells: polygons (2D), faces of the polyhedra (3D) */
        std::vector<std::vector<std::vector<size_t>>> cell_polys(cells.size());
        for (size_t c = 0; c < cells.size(); c++) {
            if constexpr (DIM == 2)
                cell_polys[c].push_back(with_hanging_nodes(cells[c].v));
            else
                for (auto& f : cell_faces(cells[c]))
                    leaf_faces(f, cell_polys[c]);
        }

        std::vector<edge_type> edges;
        for (auto& polys : cell_polys)
            for (auto& p : polys)
                for (size_t i = 0; i < p.size(); i++)
                    edges.push_back(edge_type(typename edge_type::sub_id_type(p[i]),
                                              typename edge_type::sub_id_type(p[(i+1) % p.size()])));
        sort_unique(edges);
        auto edge_of = [&](size_t a, size_t b) {
            auto e = edge_type(typename edge_type::sub_id_type(a), typename edge_type::sub_id_type(b));
            auto it = std::lower_bound(edges.begin(), edges.end(), e);
            if (it == edges.end() || !(*it == e))
                throw std::logic_error("refinable_mesh: edge not found");
            return edge_id(std::distance(edges.begin(), it));
        };

        std::vector<size_t> order(cells.size());

        if constexpr (DIM == 2) {
            std::vector<size_t> owners(edges.size(), 0);
            std::vector<std::pair<surface_type, size_t>> surfaces;
            for (size_t c = 0; c < cells.size(); c++) {
                const auto& p = cell_polys[c][0];
                std::vector<edge_id> se;
                for (size_t i = 0; i < p.size(); i++) {
                    auto e = edge_of(p[i], p[(i+1) % p.size()]);
                    se.push_back(e);
                    owners[e]++;
                }
                surface_type s(se);
                std::vector<disk::point_identifier<2>> sp;
                for (auto i : p)
                    sp.push_back(disk::point_identifier<2>(i));
                s.set_point_ids(sp.begin(), sp.end());
                surfaces.push_back(std::make_pair(s, c));
            }
            std::sort(surfaces.begin(), surfaces.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            storage->edges = edges;
            storage->boundary_info.resize(edges.size());
            for (size_t e = 0; e < edges.size(); e++)
                if (owners[e] == 1)
                    storage->boundary_info[e] = disk::boundary_descriptor(0, true);
            for (size_t k = 0; k < surfaces.size(); k++) {
                storage->surfaces.push_back(surfaces[k].first);
                order[k] = surfaces[k].second;
            }
            storage->subdomain_info.resize(storage->surfaces.size());
        }
        else {
            using volume_type = typename mesh_type::volume_type;
            using surface_id = typename surface_type::id_type;

            /* faces, identified by their points */
            struct face_def { std::vector<edge_id> edges; std::vector<size_t> points; };
            std::vector<face_def> fdefs;
            std::map<std::vector<size_t>, size_t> face_index;
            std::vector<std::vector<size_t>> vol_faces(cells.size());
            for (size_t c = 0; c < cells.size(); c++)
                for (auto& fp : cell_polys[c]) {
                    auto key = sorted(fp);
                    auto it = face_index.find(key);
                    if (it == face_index.end()) {
                        face_def fd;
                        fd.points = fp;
                        for (size_t i = 0; i < fp.size(); i++)
                            fd.edges.push_back(edge_of(fp[i], fp[(i+1) % fp.size()]));
                        fdefs.push_back(fd);
                        it = face_index.insert(std::make_pair(key, fdefs.size() - 1)).first;
                    }
                    vol_faces[c].push_back(it->second);
                }

            std::vector<size_t> forder(fdefs.size());
            for (size_t i = 0; i < forder.size(); i++)
                forder[i] = i;
            std::sort(forder.begin(), forder.end(), [&](size_t a, size_t b) { return fdefs[a].edges < fdefs[b].edges; });
            std::vector<size_t> new_index(fdefs.size());
            std::vector<surface_type> surfaces;
            for (size_t i = 0; i < forder.size(); i++) {
                new_index[forder[i]] = i;
                const auto& fd = fdefs[forder[i]];
                surface_type s(fd.edges);
                std::vector<disk::point_identifier<3>> sp;
                for (auto p : fd.points)
                    sp.push_back(disk::point_identifier<3>(p));
                s.set_point_ids(sp.begin(), sp.end());
                surfaces.push_back(s);
            }

            std::vector<size_t> owners(fdefs.size(), 0);
            std::vector<std::pair<volume_type, size_t>> volumes;
            for (size_t c = 0; c < cells.size(); c++) {
                std::vector<surface_id> vf;
                std::vector<size_t> vp;
                for (auto f : vol_faces[c]) {
                    vf.push_back(surface_id(new_index[f]));
                    owners[new_index[f]]++;
                    vp.insert(vp.end(), fdefs[f].points.begin(), fdefs[f].points.end());
                }
                sort_unique(vp);
                volume_type vol(vf);
                std::vector<disk::point_identifier<3>> vpi;
                for (auto p : vp)
                    vpi.push_back(disk::point_identifier<3>(p));
                vol.set_point_ids(vpi.begin(), vpi.end());
                volumes.push_back(std::make_pair(vol, c));
            }
            std::sort(volumes.begin(), volumes.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });

            storage->edges = edges;
            storage->surfaces = surfaces;
            storage->boundary_info.resize(surfaces.size());
            for (size_t s = 0; s < surfaces.size(); s++)
                if (owners[s] == 1)
                    storage->boundary_info[s] = disk::boundary_descriptor(0, true);
            for (size_t k = 0; k < volumes.size(); k++) {
                storage->volumes.push_back(volumes[k].first);
                order[k] = volumes[k].second;
            }
            storage->subdomain_info.resize(storage->volumes.size());
        }
        return order;
    }
};

/* The background meshes of meshes.hpp as refinable meshes */

/* n x n squares of the unit square, split into two triangles if "triangles" */
template<typename T>
void
refinable_square_grid(refinable_mesh<T,2>& rm, size_t n, bool triangles)
{
    T h = 1.0/n;
    auto id = [&](size_t i, size_t j) { return j*(n+1) + i; };
    for (size_t j = 0; j <= n; j++)
        for (size_t i = 0; i <= n; i++)
            rm.points.push_back(disk::point<T,2>(i*h, j*h));
    for (size_t j = 0; j < n; j++)
        for (size_t i = 0; i < n; i++) {
            auto a = id(i,j), b = id(i+1,j), c = id(i,j+1), d = id(i+1,j+1);
            if (triangles) {
                rm.cells.push_back({refinable_shape::polygon, {a, b, d}, 0});
                rm.cells.push_back({refinable_shape::polygon, {a, d, c}, 0});
            }
            else
                rm.cells.push_back({refinable_shape::polygon, {a, b, d, c}, 0});
        }
}

/* polygons of a 2D mesh, vertices counterclockwise */
template<typename T>
std::vector<std::vector<size_t>>
counterclockwise_polygons(const disk::generic_mesh<T,2>& m2)
{
    std::vector<std::vector<size_t>> ret;
    std::vector<typename disk::generic_mesh<T,2>::point_type> pts(m2.points_begin(), m2.points_end());
    for (auto& cl : m2) {
        auto ids = cl.point_ids();
        std::vector<size_t> pv(ids.begin(), ids.end());
        disk::point<T,2> c;
        for (auto p : pv)
            c = c + pts[p];
        c = c / T(pv.size());
        auto angle = [&](size_t p) { return std::atan2(pts[p].y() - c.y(), pts[p].x() - c.x()); };
        std::sort(pv.begin(), pv.end(), [&](size_t a, size_t b) { return angle(a) < angle(b); });
        ret.push_back(pv);
    }
    return ret;
}

template<typename T>
void
refinable_polygonal_mesh(refinable_mesh<T,2>& rm, const disk::generic_mesh<T,2>& m2)
{
    rm.points.assign(m2.points_begin(), m2.points_end());
    for (auto& pv : counterclockwise_polygons(m2))
        rm.cells.push_back({refinable_shape::polygon, pv, 0});
}

/* n x n x n cubes of the unit cube, or their Kuhn tetrahedra */
template<typename T>
void
refinable_cube_grid(refinable_mesh<T,3>& rm, size_t n, bool tetrahedra)
{
    T h = 1.0/n;
    auto id = [&](size_t i, size_t j, size_t k) { return (k*(n+1) + j)*(n+1) + i; };
    for (size_t k = 0; k <= n; k++)
        for (size_t j = 0; j <= n; j++)
            for (size_t i = 0; i <= n; i++)
                rm.points.push_back(disk::point<T,3>(i*h, j*h, k*h));
    const std::array<std::array<size_t,3>, 6> perms = {{ {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} }};
    for (size_t k = 0; k < n; k++)
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < n; i++) {
                if (!tetrahedra) {
                    rm.cells.push_back({refinable_shape::prism,
                                        {id(i,j,k), id(i+1,j,k), id(i+1,j+1,k), id(i,j+1,k),
                                         id(i,j,k+1), id(i+1,j,k+1), id(i+1,j+1,k+1), id(i,j+1,k+1)}, 0});
                    continue;
                }
                for (auto& pm : perms) {
                    std::array<size_t,3> c = {i, j, k};
                    std::vector<size_t> v = {id(c[0], c[1], c[2])};
                    for (size_t s = 0; s < 3; s++) {
                        c[pm[s]]++;
                        v.push_back(id(c[0], c[1], c[2]));
                    }
                    rm.cells.push_back({refinable_shape::tetrahedron, v, 0});
                }
            }
}

/* prisms on the polygons of a 2D mesh, in nz layers */
template<typename T>
void
refinable_prism_mesh(refinable_mesh<T,3>& rm, const disk::generic_mesh<T,2>& m2, size_t nz)
{
    size_t n2 = m2.points_size();
    T hz = 1.0/nz;
    std::vector<typename disk::generic_mesh<T,2>::point_type> pts2(m2.points_begin(), m2.points_end());
    for (size_t k = 0; k <= nz; k++)
        for (size_t p = 0; p < n2; p++)
            rm.points.push_back(disk::point<T,3>(pts2[p].x(), pts2[p].y(), k*hz));
    for (size_t k = 0; k < nz; k++)
        for (auto& pv : counterclockwise_polygons(m2)) {
            std::vector<size_t> v;
            for (auto p : pv)
                v.push_back(k*n2 + p);
            for (auto p : pv)
                v.push_back((k+1)*n2 + p);
            rm.cells.push_back({refinable_shape::prism, v, 0});
        }
}

/* Statistics of the refinement */
struct refinement_statistics {
    size_t refined = 0;      // cells refined (all the steps)
    size_t levels = 0;       // largest number of refinements of a cell
    size_t unresolved = 0;   // cells still not cut once (max_levels reached)
};

/* Refine the cells of rm that are not cut once by the interface, at most
 * max_levels times each, and build the refined mesh in msh. The refinement
 * is balanced (2:1): two cells sharing a point differ by at most one level,
 * so that an ill-cut cell and the neighbor that stabilizes it have
 * comparable sizes (the polynomial of the neighbor is extended to the
 * ill-cut cell). */
template<typename T, size_t DIM>
refinement_statistics
refine_until_cut_once(refinable_mesh<T,DIM>& rm, disk::generic_mesh<T,DIM>& msh,
                      const level_set<T,DIM>& ls, const disk::cut_geometry_options& geometry,
                      size_t max_levels)
{
    refinement_statistics stats;
    while (true) {
        disk::cut_mesh<disk::generic_mesh<T,DIM>> tmp;
        auto order = rm.build(static_cast<disk::generic_mesh<T,DIM>&>(tmp));
        auto bad = disk::cut_once_violations(tmp, ls, geometry);

        size_t n = rm.cells.size();
        std::vector<bool> marked(n, false);
        std::vector<size_t> queue;
        stats.unresolved = 0;
        for (size_t k = 0; k < bad.size(); k++) {
            if (!bad[k])
                continue;
            if (rm.cells[order[k]].level < max_levels) {
                marked[order[k]] = true;
                queue.push_back(order[k]);
            }
            else
                stats.unresolved++;
        }
        if (queue.empty())
            break;

        /* balance: the cells sharing a point with a refined cell are refined
         * if they would be two levels coarser */
        std::vector<size_t> where(n);
        for (size_t k = 0; k < n; k++)
            where[order[k]] = k;
        std::vector<std::vector<size_t>> point_cells(rm.points.size());
        for (size_t k = 0; k < n; k++)
            for (auto p : tmp[k].point_ids())
                point_cells[p].push_back(order[k]);
        auto new_level = [&](size_t c) { return rm.cells[c].level + (marked[c] ? 1 : 0); };
        while (!queue.empty()) {
            size_t c = queue.back();
            queue.pop_back();
            for (auto p : tmp[where[c]].point_ids())
                for (auto nb : point_cells[p])
                    if (!marked[nb] && new_level(c) > rm.cells[nb].level + 1) {
                        marked[nb] = true;
                        queue.push_back(nb);
                    }
        }

        std::vector<size_t> to_refine;
        for (size_t c = 0; c < n; c++)
            if (marked[c])
                to_refine.push_back(c);
        rm.refine(to_refine);
        stats.refined += to_refine.size();
    }
    stats.levels = rm.max_level();
    rm.build(msh);
    return stats;
}

} // namespace unfitted_meshes
