/*
 * Meshes cut by an interface (unfitted methods).
 *
 * cut_mesh<Mesh> is a mesh of diskpp of any type (cartesian, simplicial or
 * generic, 2D or 3D), together with the description of its cut by the zero
 * level set of a function: location of the nodes, faces and cells, geometry
 * of the cut faces and cells, quadratures of their parts, neighbors and
 * pairing of the ill-cut cells. Each cell is assumed to be cut once.
 *
 * Geometry of the cut (see cut_geometry.hpp):
 *  - 2D: a cut cell is a polygon cut once; the interface is a polyline (2^r
 *    segments, ProtoN), the polynomial interpolation of the heights of the
 *    interface above its chord, or the exact graph of the level set above it
 *    (implicit, default). The parts are integrated by sectors from a point
 *    with respect to which they are star-shaped.
 *  - 3D: the cut faces are polygons cut once, treated as above in their
 *    plane. In a cut cell, the curves of its cut faces form a loop; the
 *    interface is the graph, in the direction nu, of the level set above the
 *    region D of a plane enclosed by the projection of the loop. D is
 *    integrated by sectors from its center, and the interface by lifting
 *    these points onto the zero level set. Each part of the cell is
 *    integrated by cones from a point with respect to which it is
 *    star-shaped, the bases of the cones being its flat faces and the
 *    interface. Cells where these assumptions fail are integrated by
 *    subdivision into tetrahedra.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "mesh.hpp"
#include "level_set.hpp"
#include "cut_types.hpp"
#include "cut_geometry.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/quadratures/quadratures.hpp"

namespace disk {

enum class connectivity_via {
    node,
    edge,
    face,
    undef
};

inline std::ostream&
operator<<(std::ostream& os, connectivity_via cv) {
    switch (cv) {
        case connectivity_via::node:
            os << "node";
            break;
        case connectivity_via::edge:
            os << "edge";
            break;
        case connectivity_via::face:
            os << "face";
            break;
        case connectivity_via::undef:
            os << "undef";
            break;
    }
    return os;
}

template<typename Mesh>
class mesh_connectivity
{
    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using conn_info_t = std::pair<size_t, connectivity_via>;

    std::vector<std::vector<size_t>>        n2cs_;
    std::vector<std::vector<size_t>>        neighbours_;
    std::vector<std::vector<conn_info_t>>   neighbours_via_;

public:
    mesh_connectivity();
    mesh_connectivity(const mesh_type& msh) {
        n2cs_.resize( msh.points_size() );
        for (auto& cl : msh) {
            auto myofs = offset(msh, cl);
            auto ptids = cl.point_ids();
            for (size_t i = 0; i < ptids.size(); i++) {
                n2cs_[ ptids[i] ].push_back(myofs);
            }
        }

        neighbours_.resize( msh.cells_size() );
        for (auto& cl : msh) {
            auto myofs = offset(msh, cl);
            auto ptids = cl.point_ids();
            for (size_t i = 0; i < ptids.size(); i++) {
                const auto& neighs = n2cs_[ptids[i]];
                neighbours_[myofs].insert(neighbours_[myofs].end(),
                    neighs.begin(), neighs.end());
            }

            std::sort(neighbours_[myofs].begin(),
                neighbours_[myofs].end());
        }

        neighbours_via_.resize(msh.cells_size());
        for (size_t ofs = 0; ofs < neighbours_.size(); ofs++) {
            const auto& neighs = neighbours_[ofs];
            size_t i = 0;
            while ( i < neighs.size() ) {
                /* neighs is sorted: count how many nodes are shared with neighs[i] */
                size_t curr_num = neighs[i];
                size_t curr_count = 0;
                while (i < neighs.size() && neighs[i] == curr_num) {
                    curr_count++;
                    i++;
                }
                if (curr_num == ofs)
                    continue;
                conn_info_t ci;
                ci.first = curr_num;
                if (curr_count == 1) {
                    ci.second = connectivity_via::node;
                } else if (curr_count == 2) {
                    ci.second = (Mesh::dimension == 2) ? connectivity_via::face : connectivity_via::edge;
                } else {
                    ci.second = connectivity_via::face;
                }

                neighbours_via_[ofs].push_back(ci);
            }
        }
    }

    std::vector<conn_info_t>
    neighbours_offsets(const mesh_type& msh, const cell_type& cl)
    {
        auto ofs = offset(msh, cl);
        return {neighbours_via_[ofs].begin(), neighbours_via_[ofs].end()};
    }

    std::vector<std::pair<cell_type, connectivity_via>>
    neighbours(const mesh_type& msh, const cell_type& cl)
    {
        auto ofs = offset(msh, cl);
        std::vector<std::pair<cell_type, connectivity_via>> ret;
        for (auto& c : neighbours_via_[ofs])
            ret.push_back( std::make_pair(msh[c.first], c.second) );
        return ret;
    }
};

/* Options of the geometry of the cut */
struct cut_geometry_options {
    interface_kind  kind = interface_kind::implicit; // interface in the cut cells (2D) and faces (3D)
    size_t          degree = 0;             // degree of the polynomial interface
    size_t          refine_levels = 8;      // 2D: polyline of 2^levels segments (and fallback)
    size_t          subdivision_levels = 3; // 3D: subdivision of the cells where the cones fail
};

/* Bounded cache of the quadratures of the cut cells and faces (the most
 * recently used are kept) */
template<typename T, size_t DIM>
class cut_quadrature_cache {
public:
    struct entry {
        std::vector<quadrature_point<T, DIM>>   qps;
        std::vector<Eigen::Matrix<T, DIM, 1>>   normals;
        bool                                    subdivided = false; // computed by subdivision
    };
    using key_type = std::array<size_t, 4>;

private:
    using list_type = std::list<key_type>;
    std::map<key_type, std::pair<entry, typename list_type::iterator>> m_map;
    list_type   m_order;
    size_t      m_capacity;

public:
    explicit cut_quadrature_cache(size_t capacity = 1024) : m_capacity(capacity) {}

    const entry*
    find(const key_type& k) {
        auto it = m_map.find(k);
        if (it == m_map.end())
            return nullptr;
        m_order.splice(m_order.begin(), m_order, it->second.second);
        return &it->second.first;
    }

    const entry&
    insert(const key_type& k, entry&& e) {
        while (m_map.size() >= m_capacity && !m_order.empty()) {
            m_map.erase(m_order.back());
            m_order.pop_back();
        }
        m_order.push_front(k);
        auto& slot = m_map[k];
        slot = std::make_pair(std::move(e), m_order.begin());
        return slot.first;
    }

    void clear() { m_map.clear(); m_order.clear(); }
};

/* A mesh of diskpp cut by the zero level set of a function */
template<typename Mesh>
class cut_mesh : public Mesh {
public:
    using background_mesh_type = Mesh;
    using coordinate_type = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using level_set_type = level_set<coordinate_type, Mesh::dimension>;

    std::vector<cut_cell_info<coordinate_type, Mesh::dimension>> cut_cell_data;
    std::vector<cut_face_info<coordinate_type, Mesh::dimension>> cut_face_data;
    std::vector<cut_node_info<coordinate_type, Mesh::dimension>> cut_node_data; // per point
    std::map<std::pair<size_t, size_t>, point_type>  edge_intersections;       // 3D, edges by point ids
    std::shared_ptr<const level_set_type>           level_set_function;
    cut_geometry_options                            geometry;

    mutable cut_quadrature_cache<coordinate_type, Mesh::dimension> quadrature_cache;

    cut_mesh() = default;
    explicit cut_mesh(const Mesh& msh) : Mesh(msh) {}

    const Mesh& background() const { return *this; }
};

template<typename M>
concept CutMesh = requires(const M& m) {
    typename M::background_mesh_type;
    m.cut_cell_data;
};

/* The geometric queries of diskpp are forwarded to the background mesh: they
 * are overloaded for the mesh types (or generic templates which would be
 * preferred for a derived class) */
template<typename Mesh>
auto points(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl) { return points(msh.background(), cl); }
template<typename Mesh>
auto points(const cut_mesh<Mesh>& msh, const typename Mesh::face_type& fc) { return points(msh.background(), fc); }
template<typename Mesh>
auto faces(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl) { return faces(msh.background(), cl); }
template<typename Mesh>
auto barycenter(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl) { return barycenter(msh.background(), cl); }
template<typename Mesh>
auto barycenter(const cut_mesh<Mesh>& msh, const typename Mesh::face_type& fc) { return barycenter(msh.background(), fc); }
template<typename Mesh>
auto diameter(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl) { return diameter(msh.background(), cl); }
template<typename Mesh>
auto diameter(const cut_mesh<Mesh>& msh, const typename Mesh::face_type& fc) { return diameter(msh.background(), fc); }
template<typename Mesh>
auto measure(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl) { return measure(msh.background(), cl); }
template<typename Mesh>
auto measure(const cut_mesh<Mesh>& msh, const typename Mesh::face_type& fc) { return measure(msh.background(), fc); }
template<typename Mesh>
auto integrate(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl, size_t degree) {
    return integrate(msh.background(), cl, degree);
}
template<typename Mesh>
auto integrate(const cut_mesh<Mesh>& msh, const typename Mesh::face_type& fc, size_t degree) {
    return integrate(msh.background(), fc, degree);
}
template<typename Mesh>
auto normal(const cut_mesh<Mesh>& msh, const typename Mesh::cell_type& cl, const typename Mesh::face_type& fc) {
    return normal(msh.background(), cl, fc);
}

template<typename Mesh, typename Element>
std::vector<typename Mesh::node_type>
nodes(const Mesh& msh, const Element& elem)
{
    auto ptids = elem.point_ids();

    auto storage = msh.backend_storage();

    auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto {
        return *std::next(storage->nodes.begin(), pi);
    };

    std::vector<typename Mesh::node_type> nds(ptids.size());
    std::transform(ptids.begin(), ptids.end(), nds.begin(), ptid_to_point);

    return nds;
}

/* The cut-cell geometry walks along the boundary of the cells, whereas the
 * storage order of the vertices depends on the mesh type (cartesian cells are
 * stored in tensor order {BL, BR, TL, TR}). This returns the permutation of the
 * stored vertices that lists them counter-clockwise, starting from the first
 * stored vertex. */
template<mesh_2D Mesh>
std::vector<size_t>
ccw_permutation(const Mesh& msh, const typename Mesh::cell_type& cl) {

    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    std::vector<size_t> perm(pts.size());
    std::iota(perm.begin(), perm.end(), 0);
    auto angle = [&](size_t i) {
        auto d = pts[i] - bar;
        return std::atan2(d.y(), d.x());
    };
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b) { return angle(a) < angle(b); });
    std::rotate(perm.begin(), std::find(perm.begin(), perm.end(), 0), perm.end());

    return perm;
}

template<mesh_2D Mesh>
std::vector<typename Mesh::node_type>
ccw_nodes(const Mesh& msh, const typename Mesh::cell_type& cl) {

    auto nds  = nodes(msh, cl);
    auto perm = ccw_permutation(msh, cl);
    std::vector<typename Mesh::node_type> ret;
    for (auto p : perm)
        ret.push_back(nds[p]);

    return ret;
}

template<mesh_2D Mesh>
std::vector<typename Mesh::point_type>
ccw_points(const Mesh& msh, const typename Mesh::cell_type& cl) {

    auto pts  = points(msh, cl);
    auto perm = ccw_permutation(msh, cl);
    std::vector<typename Mesh::point_type> ret;
    for (auto p : perm)
        ret.push_back(pts[p]);

    return ret;
}

/* Faces in counter-clockwise order: face i joins the ccw vertices i and i+1. */
template<mesh_2D Mesh>
std::vector<typename Mesh::face_type>
ccw_faces(const Mesh& msh, const typename Mesh::cell_type& cl) {

    using face_type = typename Mesh::face_type;

    auto fcs   = faces(msh, cl);
    auto ptids = cl.point_ids();
    auto perm  = ccw_permutation(msh, cl);

    std::vector<face_type> ret;
    for (size_t i = 0; i < perm.size(); i++) {
        size_t a = ptids[perm[i]];
        size_t b = ptids[perm[(i+1) % perm.size()]];
        auto itor = std::find_if(fcs.begin(), fcs.end(), [&](const face_type& fc) {
            auto fids = fc.point_ids();
            return (size_t(fids[0]) == a && size_t(fids[1]) == b) || (size_t(fids[0]) == b && size_t(fids[1]) == a);
        });
        if (itor == fcs.end())
            throw std::logic_error("ccw_faces: face not found");
        ret.push_back(*itor);
    }

    return ret;
}

/***************************************************************************/
/*                        LOCATION OF THE ELEMENTS                          */
/***************************************************************************/

template<CutMesh Mesh>
location
locate(const Mesh& msh, const typename Mesh::cell_type& cl) {
    auto& cl_loc = msh.cut_cell_data[offset(msh, cl)].loc;
    assert(cl_loc != location::UNDEF);
    return cl_loc;
}

template<CutMesh Mesh>
location
locate(const Mesh& msh, const typename Mesh::face_type& fc) {
    auto& fc_loc = msh.cut_face_data[offset(msh, fc)].loc;
    assert(fc_loc != location::UNDEF);
    return fc_loc;
}

template<CutMesh Mesh>
location
locate(const Mesh& msh, const typename Mesh::node_type& nd) {
    auto& nd_loc = msh.cut_node_data[nd.point_ids()[0]].loc;
    assert(nd_loc != location::UNDEF);
    return nd_loc;
}

template<CutMesh Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::cell_type& cl) {
    return msh.cut_cell_data[offset(msh, cl)].loc == location::ON_INTERFACE;
}

template<CutMesh Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::face_type& fc) {
    return msh.cut_face_data[offset(msh, fc)].loc == location::ON_INTERFACE;
}

template<CutMesh Mesh>
const auto&
cut_info(const Mesh& msh, const typename Mesh::cell_type& cl) {
    return msh.cut_cell_data[offset(msh, cl)];
}

template<CutMesh Mesh>
const auto&
cut_info(const Mesh& msh, const typename Mesh::face_type& fc) {
    return msh.cut_face_data[offset(msh, fc)];
}

/* Level set of a 2D mesh, or restricted to the plane of a cut face of a 3D mesh */
template<CutMesh Mesh>
cutgeo::planar_level_set<typename Mesh::coordinate_type, Mesh::dimension>
cell_level_set(const Mesh& msh) requires (Mesh::dimension == 2) {
    return cutgeo::planar_level_set<typename Mesh::coordinate_type, 2>(*msh.level_set_function);
}

template<CutMesh Mesh>
cutgeo::planar_level_set<typename Mesh::coordinate_type, 3>
face_level_set(const Mesh& msh, const cut_face_info<typename Mesh::coordinate_type, 3>& fd) {
    return cutgeo::planar_level_set<typename Mesh::coordinate_type, 3>(*msh.level_set_function, fd.origin, fd.U);
}

/***************************************************************************/
/*                    GEOMETRY OF THE CUT: 2D (POLYGONS)                    */
/***************************************************************************/

/* Polygon of the subcell "where", the interface being replaced by its chord if
 * chord_only (a convex polygon for a convex cell cut by a straight line) */
template<CutMesh MeshType>
std::vector<typename MeshType::point_type>
collect_triangulation_points(const MeshType& msh, const typename MeshType::cell_type& cl, location where,
                             bool chord_only = false) requires (MeshType::dimension == 2) {
    assert(is_cut(msh, cl));
    return cutgeo::side_polygon(cut_info(msh, cl).poly, where, chord_only);
}

template<CutMesh MeshType>
bool
has_curved_interface(const MeshType& msh, const typename MeshType::cell_type& cl) requires (MeshType::dimension == 2) {
    return is_cut(msh, cl) && cut_info(msh, cl).poly.kind != interface_kind::polyline;
}

/* Center of the fan triangulating the subcell of a quadrangle with a polyline
 * interface (ProtoN) */
template<CutMesh MeshType>
typename MeshType::point_type
tesselation_center(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    auto fcs = ccw_faces(msh, cl);
    auto pts = ccw_points(msh, cl);
    auto nds = ccw_nodes(msh, cl);

    if (fcs.size() != 4)
        throw std::invalid_argument("This works only on quads for now");

    if( !is_cut(msh, cl) )
        throw std::invalid_argument("No tesselation centers for uncut cells");

    // if two consecutive faces are cut
    // return either the common node or the opposite node
    for (size_t i = 0; i < fcs.size(); i++) {
        auto f1 = i;
        auto f2 = (i+1) % fcs.size();
        auto n = (i+1) % fcs.size();

        if (is_cut(msh, fcs[f1]) && is_cut(msh, fcs[f2])) {
            if (locate(msh, nds[n]) == where)
                return pts[n];
            else
                return pts[(n+2)%4];
        }
    }

    // if two opposite faces are cut
    // return the center of one of the other faces
    for (size_t i = 0; i < 2; i++) {
        auto f1 = i;
        auto f2 = i+2;
        auto n = i+1;
        if (is_cut(msh, fcs[f1]) && is_cut(msh, fcs[f2])) {
            if(locate(msh, nds[n]) == where )
                return 0.5*(pts[n] + pts[n+1]);
            else
                return 0.5*(pts[(n+2)%4] + pts[(n+3)%4]);
        }
    }

    // normally the tesselation center is already found
    throw std::logic_error("we shouldn't arrive here !!");
}

/* End points of the part of the face lying in the subdomain "where" */
template<CutMesh MeshType>
std::array<typename MeshType::point_type, 2>
points(const MeshType& msh, const typename MeshType::face_type& fc, location where) requires (MeshType::dimension == 2) {

    if (locate(msh, fc) != where && locate(msh, fc) != location::ON_INTERFACE)
        throw std::invalid_argument("This face has no points where requested");

    auto pts = points(msh, fc);
    std::array<typename MeshType::point_type, 2> ret = {pts[0], pts[1]};
    if (!is_cut(msh, fc))
        return ret;

    auto ip  = cut_info(msh, fc).intersection_point;
    auto nds = nodes(msh, fc);
    if (locate(msh, nds[0]) == where && locate(msh, nds[1]) != where)
        ret[1] = ip;
    else if (locate(msh, nds[0]) != where && locate(msh, nds[1]) == where)
        ret[0] = ip;
    else
        throw std::logic_error("Invalid point configuration");

    return ret;
}

/***************************************************************************/
/*                     GEOMETRY OF THE CUT: 3D (POLYHEDRA)                  */
/***************************************************************************/

/* The cell as a polyhedron: faces with their vertices counterclockwise seen
 * from outside and their outward normals */
template<CutMesh Mesh>
cut_polyhedron<typename Mesh::coordinate_type>
cell_polyhedron(const Mesh& msh, const typename Mesh::cell_type& cl) requires (Mesh::dimension == 3)
{
    using T = typename Mesh::coordinate_type;
    cut_polyhedron<T> P;
    for (auto& fc : faces(msh, cl)) {
        cut_polyhedron_face<T> f;
        auto pts = points(msh, fc);
        f.vertices.assign(pts.begin(), pts.end());
        f.normal = normal(msh, cl, fc);
        if (cutgeo::newell_normal(f.vertices).dot(f.normal) < 0)
            std::reverse(f.vertices.begin(), f.vertices.end());
        P.faces.push_back(f);
    }
    return P;
}

template<CutMesh Mesh>
cutgeo::cut_polyhedron_options
polyhedron_options(const Mesh& msh)
{
    cutgeo::cut_polyhedron_options opt;
    opt.refine_levels = msh.geometry.refine_levels;
    opt.subdivision_levels = msh.geometry.subdivision_levels;
    return opt;
}

/* Quadrature of the part "where" of a cut face of a 3D mesh */
template<CutMesh Mesh>
std::vector<quadrature_point<typename Mesh::coordinate_type, 3>>
face_side_quadrature_3d(const Mesh& msh, const cut_face_info<typename Mesh::coordinate_type, 3>& fd,
                        size_t degree, location where)
{
    using T = typename Mesh::coordinate_type;
    auto phi = face_level_set(msh, fd);
    std::vector<quadrature_point<T,3>> ret;
    for (auto& qp : cutgeo::side_quadrature(fd.poly, phi, degree, where))
        ret.push_back(make_qp(phi.to_global(qp.first), qp.second));
    return ret;
}

/* Quadrature points and normals of the interface of a 3D cut cell (cached) */
template<CutMesh Mesh>
const typename cut_quadrature_cache<typename Mesh::coordinate_type, 3>::entry&
cached_interface_3d(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    using T = typename Mesh::coordinate_type;
    typename cut_quadrature_cache<T,3>::key_type key = {1, offset(msh, cl), 0, degree};
    if (auto e = msh.quadrature_cache.find(key))
        return *e;
    cutgeo::cut_quadratures<T,3> q;
    cutgeo::polyhedron_interface_quadrature(cut_info(msh, cl).poly, *msh.level_set_function, degree,
                                            polyhedron_options(msh), q);
    typename cut_quadrature_cache<T,3>::entry e;
    e.qps = q.itf;
    e.normals = q.normals;
    return msh.quadrature_cache.insert(key, std::move(e));
}

/***************************************************************************/
/*                    QUADRATURES AND MEASURES OF THE PARTS                 */
/***************************************************************************/

/* Quadrature on the interface of a cut cell. The weights are signed so that
 * the interface is oriented outward from the subcell "where": since the
 * interface is oriented with the negative side on its left (2D), or its
 * normal points towards the positive side (3D), the sign only depends on
 * "where". (ProtoN deduced it from the position of the barycenter of the
 * subcell, which gives a wrong sign for non-convex subcells, e.g. the flower
 * on the coarsest meshes.) */
template<CutMesh MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, MeshType::dimension>>
integrate_interface(const MeshType& msh, const typename MeshType::cell_type& cl, size_t degree, location where) {

    using T = typename MeshType::coordinate_type;
    static const size_t DIM = MeshType::dimension;

    assert(is_cut(msh, cl));
    T int_sign = (where == location::IN_POSITIVE_SIDE) ? -1.0 : +1.0;
    std::vector<disk::quadrature_point<T, DIM>> ret;

    if constexpr (DIM == 3) {
        const auto& e = cached_interface_3d(msh, cl, degree);
        for (auto& qp : e.qps)
            ret.push_back( disk::make_qp(qp.point(), int_sign * qp.weight()) );
        return ret;
    }
    else {
        const auto& poly = cut_info(msh, cl).poly;
        if (poly.kind == interface_kind::implicit) {
            auto cq = cutgeo::curve_quadrature(poly, cell_level_set(msh), degree);
            for (auto& qp : cq.qps)
                ret.push_back( disk::make_qp(qp.first, int_sign * qp.second) );
            return ret;
        }
        if (poly.kind == interface_kind::polynomial) {
            /* Gauss rule in the parameter of the curve, ds = |x'(t)| dt */
            auto phi = cell_level_set(msh);
            cutgeo::polygon_curve<T, decltype(phi)> crv(poly, phi);
            for (auto& tw : cutgeo::gauss_legendre_01<T>(crv.nodes_count(degree, false))) {
                point<T,2> x;
                cutgeo::vec2<T> dx;
                crv.eval(tw.first, x, dx);
                ret.push_back( disk::make_qp(x, int_sign * tw.second * dx.norm()) );
            }
            return ret;
        }
        const auto& interface = poly.polyline;
        for (size_t i = 1; i < interface.size(); i++) {
            auto qps = disk::quadrature::gauss_legendre(degree, interface.at(i-1), interface.at(i));
            for (auto& qp : qps)
                ret.push_back( disk::make_qp(qp.point(), int_sign * qp.weight()) );
        }
        return ret;
    }
}

/* Unit normal of the interface, one entry per quadrature point of
 * integrate_interface(). It points from the negative to the positive side. */
template<CutMesh MeshType>
std::vector<Eigen::Matrix<typename MeshType::coordinate_type, MeshType::dimension, 1>>
interface_normals(const MeshType& msh, const typename MeshType::cell_type& cl, size_t degree) {

    using T = typename MeshType::coordinate_type;
    static const size_t DIM = MeshType::dimension;

    assert(is_cut(msh, cl));

    if constexpr (DIM == 3) {
        return cached_interface_3d(msh, cl, degree).normals;
    }
    else {
        const auto& poly = cut_info(msh, cl).poly;
        std::vector<Eigen::Matrix<T, 2, 1>> ret;
        if (poly.kind == interface_kind::implicit)
            return cutgeo::curve_quadrature(poly, cell_level_set(msh), degree).normals;
        if (poly.kind == interface_kind::polynomial) {
            auto phi = cell_level_set(msh);
            cutgeo::polygon_curve<T, decltype(phi)> crv(poly, phi);
            for (auto& tw : cutgeo::gauss_legendre_01<T>(crv.nodes_count(degree, false))) {
                point<T,2> x;
                cutgeo::vec2<T> d;
                crv.eval(tw.first, x, d);
                ret.push_back( Eigen::Matrix<T, 2, 1>(d(1), -d(0)) / d.norm() );
            }
            return ret;
        }
        const auto& interface = poly.polyline;
        for (size_t i = 1; i < interface.size(); i++) {
            auto p0 = interface.at(i-1);
            auto p1 = interface.at(i);
            auto scale = p1 - p0;
            auto meas = scale.to_vector().norm();
            Eigen::Matrix<T, 2, 1> ni;
            ni(0) =  scale.y() / meas;
            ni(1) = -scale.x() / meas;
            auto nqps = disk::quadrature::gauss_legendre(degree, p0, p1).size();
            for (size_t j = 0; j < nqps; j++)
                ret.push_back(ni);
        }
        return ret;
    }
}

template<CutMesh MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, MeshType::dimension>>
integrate(const MeshType& msh, const typename MeshType::face_type& fc, size_t degree, const location& where) {

    using T = typename MeshType::coordinate_type;
    static const size_t DIM = MeshType::dimension;

    std::vector<disk::quadrature_point<T, DIM>> ret;
    if (locate(msh, fc) != where && locate(msh, fc) != location::ON_INTERFACE)
        return ret;

    if (!is_cut(msh, fc)) /* Element is not cut, use std. integration */
        return integrate(msh, fc, degree);

    if constexpr (DIM == 2) {
        auto pts = points(msh, fc, where);
        return disk::quadrature::gauss_legendre(degree, pts[0], pts[1]);
    }
    else {
        typename cut_quadrature_cache<T,3>::key_type key = {2, offset(msh, fc), size_t(where), degree};
        if (auto e = msh.quadrature_cache.find(key))
            return e->qps;
        typename cut_quadrature_cache<T,3>::entry e;
        e.qps = face_side_quadrature_3d(msh, cut_info(msh, fc), degree, where);
        return msh.quadrature_cache.insert(key, std::move(e)).qps;
    }
}

/* Fan triangulation of a subcell of a quadrangle with a polyline interface
 * (ProtoN) */
template<typename T>
struct temp_tri {

    std::array<point<T,2>,3> pts;

    T area() const {
        auto v1 = pts[1] - pts[0];
        auto v2 = pts[2] - pts[0];

        return ( v1.x()*v2.y() - v2.x()*v1.y() ) / 2.0;
        // can be negative
    }
};

template<CutMesh MeshType>
std::vector<temp_tri<typename MeshType::coordinate_type >>
triangulate(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    assert(is_cut(msh, cl));

    auto tp = collect_triangulation_points(msh, cl, where);
    auto bar = tesselation_center(msh, cl, where);

    using T = typename MeshType::coordinate_type;
    std::vector<temp_tri<T>> tris;

    for (size_t i = 0; i < tp.size(); i++) {
        temp_tri<T> t;
        t.pts[0] = bar;
        t.pts[1] = tp[i];
        t.pts[2] = tp[(i+1)%tp.size()];
        tris.push_back(t);
    }

    return tris;
}

template<CutMesh MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, MeshType::dimension>>
integrate(const MeshType& msh, const typename MeshType::cell_type& cl, size_t degree, const location& where) {

    using T = typename MeshType::coordinate_type;
    static const size_t DIM = MeshType::dimension;

    std::vector<disk::quadrature_point<T, DIM>> ret;

    if (locate(msh, cl) != where && locate(msh, cl) != location::ON_INTERFACE )
        return ret;

    if (!is_cut(msh, cl)) // Element is not cut, use std. integration
        return integrate(msh, cl, degree);

    if constexpr (DIM == 3) {
        typename cut_quadrature_cache<T,3>::key_type key = {0, offset(msh, cl), size_t(where), degree};
        if (auto e = msh.quadrature_cache.find(key))
            return e->qps;
        typename cut_quadrature_cache<T,3>::entry e;
        cutgeo::polyhedron_volume_quadrature(cut_info(msh, cl).poly, *msh.level_set_function, degree, where,
                                             polyhedron_options(msh), e.qps);
        return msh.quadrature_cache.insert(key, std::move(e)).qps;
    }
    else {
        const auto& poly = cut_info(msh, cl).poly;
        if (poly.kind != interface_kind::polyline) {
            /* polynomial and implicit interfaces: sectors from the center of
             * the subcell (positive weights, points in the subcell) */
            auto build = [&]() {
                std::vector<disk::quadrature_point<T, 2>> qps;
                for (auto& qp : cutgeo::side_quadrature(poly, cell_level_set(msh), degree, where))
                    qps.push_back( disk::make_qp(qp.first, qp.second) );
                return qps;
            };
            if (poly.kind == interface_kind::polynomial)
                return build();
            typename cut_quadrature_cache<T,2>::key_type key = {0, offset(msh, cl), size_t(where), degree};
            if (auto e = msh.quadrature_cache.find(key))
                return e->qps;
            typename cut_quadrature_cache<T,2>::entry e;
            e.qps = build();
            return msh.quadrature_cache.insert(key, std::move(e)).qps;
        }

        if (poly.vertices.size() != 4) {
            /* polyline in a polygon: ear clipping */
            for (auto& qp : cutgeo::side_quadrature(poly, cell_level_set(msh), degree, where))
                ret.push_back( disk::make_qp(qp.first, qp.second) );
            return ret;
        }

        /* The rule of the reference triangle is mapped onto each triangle of the
         * fan, the weights being scaled by the signed area computed with the
         * cross product, as in ProtoN:
         *  - the subcell need not be star-shaped with respect to the center of
         *    the fan (strongly curved interfaces on coarse meshes): signed
         *    weights keep the decomposition of the polygon exact;
         *  - the subcells can be extremely thin (square interface of Fig. 12),
         *    a case where the Kahan formula used by triangle_gauss() returns NaN. */
        const auto& ref_qps = cutgeo::reference_triangle_rule<T>(degree);
        auto tris = triangulate(msh, cl, where);
        for (auto& tri : tris) {
            auto v1 = tri.pts[1] - tri.pts[0];
            auto v2 = tri.pts[2] - tri.pts[0];
            T det = v1.x()*v2.y() - v1.y()*v2.x(); // twice the signed area
            for (auto& qp : ref_qps) {
                auto xi  = qp.first.x();
                auto eta = qp.first.y();
                ret.push_back( disk::make_qp(tri.pts[0] + v1*xi + v2*eta, qp.second*det) );
            }
        }
        return ret;
    }
}

template<CutMesh MeshType>
typename MeshType::coordinate_type
measure(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    using T = typename MeshType::coordinate_type;

    if (!is_cut(msh, cl)) /* Element is not cut, use std. integration */
        return measure(msh, cl);

    T totmeas = 0.0;
    auto qpsi = integrate(msh, cl, 0, where);
    for (auto& qp : qpsi) {
        totmeas += qp.weight();
    }

    return totmeas;
}

/* Points describing the part "where" of a cell (its vertices if it is not
 * cut): polygon of the subcell with the interface samples (2D); vertices,
 * interface curves of the faces and samples of the interface (3D) */
template<CutMesh MeshType>
std::vector<typename MeshType::point_type>
region_points(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    using T = typename MeshType::coordinate_type;
    using point_type = typename MeshType::point_type;

    if (!is_cut(msh, cl)) {
        auto pts = points(msh, cl);
        return std::vector<point_type>(pts.begin(), pts.end());
    }

    if constexpr (MeshType::dimension == 2) {
        return collect_triangulation_points(msh, cl, where);
    }
    else {
        std::vector<point_type> ret;
        auto pts = points(msh, cl);
        auto ids = cl.point_ids();
        for (size_t i = 0; i < pts.size(); i++)
            if (msh.cut_node_data[ids[i]].loc == where)
                ret.push_back(pts[i]);
        for (auto& fc : faces(msh, cl)) {
            const auto& fd = cut_info(msh, fc);
            if (fd.loc != location::ON_INTERFACE)
                continue;
            cutgeo::planar_level_set<T,3> phi(*msh.level_set_function, fd.origin, fd.U);
            for (auto& p : fd.poly.polyline)
                ret.push_back(phi.to_global(p));
        }
        cutgeo::polyhedron_interface_samples(cut_info(msh, cl).poly, *msh.level_set_function, ret);
        return ret;
    }
}

/* Centroid of the subcell "where" */
template<CutMesh MeshType>
typename MeshType::point_type
barycenter(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    using T = typename MeshType::coordinate_type;
    using point_type = typename MeshType::point_type;

    if (!is_cut(msh, cl))
        return barycenter(msh, cl);

    if constexpr (MeshType::dimension == 2) {
        if (!has_curved_interface(msh, cl)) {
            /* centroid of the polygon approximating the subcell (as in ProtoN) */
            auto tp = collect_triangulation_points(msh, cl, where);
            auto p0 = tp.front();
            point_type ret;
            T den = 0.0;
            for (size_t i = 2; i < tp.size(); i++) {
                auto pprev = tp[i-1] - p0;
                auto pcur  = tp[i] - p0;
                T d = (pprev.x()*pcur.y() - pprev.y()*pcur.x()) / 2.0;
                ret = ret + (pprev + pcur) * d;
                den += d;
            }
            return p0 + ret/(den*3);
        }
    }

    point_type c;
    T meas = 0.0;
    for (auto& qp : integrate(msh, cl, 1, where)) {
        c = c + qp.point() * qp.weight();
        meas += qp.weight();
    }
    return c / meas;
}

template<CutMesh MeshType>
typename MeshType::coordinate_type
diameter(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    if (!is_cut(msh, cl))
        return diameter(msh, cl);

    auto tp = region_points(msh, cl, where);
    typename MeshType::coordinate_type diam = 0.0;
    for (size_t i = 0; i < tp.size(); i++)
        for (size_t j = i+1; j < tp.size(); j++)
            diam = std::max(diam, (tp[i]-tp[j]).to_vector().norm());

    return diam;
}

/* Points describing the part "where" of a face (end points of the
 * subsegment in 2D, polygon of the subface with the samples of the
 * interface in 3D) */
template<CutMesh MeshType>
std::vector<typename MeshType::point_type>
face_region_points(const MeshType& msh, const typename MeshType::face_type& fc, location where) {

    if (!is_cut(msh, fc)) {
        auto pts = points(msh, fc);
        return std::vector<typename MeshType::point_type>(pts.begin(), pts.end());
    }

    if constexpr (MeshType::dimension == 2) {
        auto tp = points(msh, fc, where);
        return {tp[0], tp[1]};
    }
    else {
        using T = typename MeshType::coordinate_type;
        const auto& fd = cut_info(msh, fc);
        cutgeo::planar_level_set<T,3> phi(*msh.level_set_function, fd.origin, fd.U);
        std::vector<typename MeshType::point_type> ret;
        for (auto& p : cutgeo::side_polygon(fd.poly, where))
            ret.push_back(phi.to_global(p));
        return ret;
    }
}

template<CutMesh MeshType>
typename MeshType::point_type
barycenter(const MeshType& msh, const typename MeshType::face_type& fc, location where) {

    using T = typename MeshType::coordinate_type;

    if (!is_cut(msh, fc))
        return barycenter(msh, fc);

    if constexpr (MeshType::dimension == 2) {
        auto tp = points(msh, fc, where);
        return 0.5 * (tp[0] + tp[1]);
    }
    else {
        typename MeshType::point_type c;
        T meas = 0;
        for (auto& qp : integrate(msh, fc, 1, where)) {
            c = c + qp.point() * qp.weight();
            meas += qp.weight();
        }
        return c / meas;
    }
}

template<CutMesh MeshType>
typename MeshType::coordinate_type
diameter(const MeshType& msh, const typename MeshType::face_type& fc, location where) {

    if (!is_cut(msh, fc))
        return diameter(msh, fc);

    auto tp = face_region_points(msh, fc, where);
    typename MeshType::coordinate_type diam = 0.0;
    for (size_t i = 0; i < tp.size(); i++)
        for (size_t j = i+1; j < tp.size(); j++)
            diam = std::max(diam, (tp[i]-tp[j]).to_vector().norm());
    return diam;
}

/***************************************************************************/
/*                         DETECTION OF THE CUT                             */
/***************************************************************************/

/* Statistics of the construction of a cut mesh */
struct cut_mesh_statistics {
    size_t cut_cells = 0;
    size_t clamped = 0;     // cells whose interface left the cell (a face or an edge crossed twice), projected onto its faces
    size_t fallbacks = 0;   // cut cells (2D) or faces (3D) keeping the polyline
    size_t split = 0;       // 3D: cut cells split into polyhedra
    size_t subdivided = 0;  // 3D: cut cells with parts integrated by subdivision
};

template<CutMesh MeshType>
void
detect_node_position(MeshType& msh) {

    const auto& ls = *msh.level_set_function;
    size_t nb_points = msh.points_size();
    msh.cut_node_data.assign(nb_points, cut_node_info<typename MeshType::coordinate_type, MeshType::dimension>());
    size_t i = 0;
    for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++, i++)
        msh.cut_node_data[i].loc = (ls(*itor) < 0) ? location::IN_NEGATIVE_SIDE : location::IN_POSITIVE_SIDE;
}

/* Location of the faces, intersection points (2D) and cut polygons (3D) */
template<CutMesh MeshType>
void
detect_cut_faces(MeshType& msh, cut_mesh_statistics& stats) {

    using T = typename MeshType::coordinate_type;
    using point_type = typename MeshType::point_type;
    static const size_t DIM = MeshType::dimension;

    const auto& ls = *msh.level_set_function;
    msh.cut_face_data.assign(msh.faces_size(), cut_face_info<T, DIM>());
    msh.edge_intersections.clear();

    for (auto& fc : faces(msh)) {
        auto fc_id = offset(msh, fc);
        auto& fd = msh.cut_face_data[fc_id];
        auto pts = points(msh, fc);
        auto ids = fc.point_ids();
        bool has_neg = false, has_pos = false;
        for (auto id : ids) {
            if (msh.cut_node_data[id].loc == location::IN_NEGATIVE_SIDE)
                has_neg = true;
            else
                has_pos = true;
        }
        if (!has_neg) {
            fd.loc = location::IN_POSITIVE_SIDE;
            continue;
        }
        if (!has_pos) {
            fd.loc = location::IN_NEGATIVE_SIDE;
            continue;
        }
        fd.loc = location::ON_INTERFACE;

        if constexpr (DIM == 2) {
            auto threshold = diameter(msh, fc) / 1e20;
            auto phi = cutgeo::planar_level_set<T,2>(ls);
            fd.intersection_point = cutgeo::find_zero_crossing(pts[0], pts[1], phi, threshold);
        }
        else {
            /* the face in the frame of its plane, its vertices counterclockwise */
            cutgeo::polygon_frame(pts, fd.origin, fd.U);
            cutgeo::planar_level_set<T,3> phi(ls, fd.origin, fd.U);
            std::vector<point<T,2>> lv;
            std::vector<location> lloc;
            std::vector<size_t> lid;
            for (size_t i = 0; i < pts.size(); i++) {
                lv.push_back(phi.to_local(pts[i]));
                lloc.push_back(msh.cut_node_data[ids[i]].loc);
                lid.push_back(ids[i]);
            }
            auto edge_point = [&](size_t i) {
                size_t a = lid[i], b = lid[(i+1) % lid.size()];
                auto key = std::make_pair(std::min(a, b), std::max(a, b));
                auto it = msh.edge_intersections.find(key);
                if (it == msh.edge_intersections.end()) {
                    auto x = cutgeo::edge_zero(pts[i], pts[(i+1) % pts.size()], ls);
                    it = msh.edge_intersections.insert(std::make_pair(key, x)).first;
                }
                return phi.to_local(it->second);
            };
            if (!cutgeo::init_cut_polygon(fd.poly, lv, lloc, lid, edge_point))
                throw std::logic_error("A face is cut more than once: the mesh does not resolve the interface");
            if (msh.geometry.kind == interface_kind::polyline || !cutgeo::make_implicit_interface(fd.poly, phi)) {
                if (msh.geometry.kind != interface_kind::polyline)
                    stats.fallbacks++;
                cutgeo::refine_polyline(fd.poly, phi, std::max<size_t>(msh.geometry.refine_levels, 1));
            }
        }
    }
}


/* Location of the cells, geometry of the cut cells (polygons in 2D, loops
 * and cones in 3D) */
template<CutMesh MeshType>
void
detect_cut_cells(MeshType& msh, cut_mesh_statistics& stats) {

    using T = typename MeshType::coordinate_type;
    static const size_t DIM = MeshType::dimension;

    msh.cut_cell_data.assign(msh.cells_size(), cut_cell_info<T, DIM>());
    const auto& ls = *msh.level_set_function;

    for (auto& cl : msh) {
        auto cl_id = offset(msh, cl);
        auto& ccd = msh.cut_cell_data[cl_id];

        bool cut = false;
        for (auto& fc : faces(msh, cl))
            cut = cut || is_cut(msh, fc);
        if (!cut) {
            auto ids = cl.point_ids();
            ccd.loc = msh.cut_node_data[ids[0]].loc;
            continue;
        }
        ccd.loc = location::ON_INTERFACE;
        stats.cut_cells++;

        if constexpr (DIM == 2) {
            auto pts = ccw_points(msh, cl);
            auto nds = ccw_nodes(msh, cl);
            auto fcs = ccw_faces(msh, cl);
            std::vector<location> vloc;
            std::vector<size_t> vid;
            for (auto& nd : nds) {
                vloc.push_back(locate(msh, nd));
                vid.push_back(nd.point_ids()[0]);
            }
            auto edge_point = [&](size_t i) { return cut_info(msh, fcs[i]).intersection_point; };
            if (!cutgeo::init_cut_polygon(ccd.poly, pts, vloc, vid, edge_point))
                throw std::logic_error("invalid number of cuts in cell");

            auto phi = cell_level_set(msh);
            if (cutgeo::refine_polyline(ccd.poly, phi, msh.geometry.refine_levels))
                stats.clamped++;
            bool ok = true;
            if (msh.geometry.kind == interface_kind::polynomial && msh.geometry.degree >= 2)
                ok = cutgeo::make_polynomial_interface(ccd.poly, phi, msh.geometry.degree);
            else if (msh.geometry.kind == interface_kind::implicit)
                ok = cutgeo::make_implicit_interface(ccd.poly, phi);
            if (!ok)
                stats.fallbacks++;
        }
        else {
            ccd.poly = cell_polyhedron(msh, cl);
            cutgeo::build_cut_polyhedron(ccd.poly, ls, polyhedron_options(msh));
            if (ccd.poly.cut_more_than_once)
                throw std::logic_error("A cell is cut more than once: the mesh does not resolve the interface");
            size_t split = 0, linear = 0, projected = 0;
            cutgeo::polyhedron_statistics(ccd.poly, split, linear, projected);
            if (split > 0)
                stats.split++;
            if (linear > 0)
                stats.subdivided++;
            if (projected > 0)
                stats.clamped++;
        }
    }
}

/* Flag the cut cells as well-cut (T_OK) or ill-cut on one side (T_KO_NEG,
 * T_KO_POS). A side is ill-cut if its measure is below threshold*|T|, or, for
 * the quadrangles (ProtoN), if the interface cuts the faces at a relative
 * distance below threshold from the vertices of that side. threshold = 0
 * flags no cell at all. */
template<CutMesh MeshType>
void
detect_cut_type(MeshType& msh, typename MeshType::coordinate_type threshold = 0.3) {

    const auto threshold_cells = threshold;

    for (auto& cl : cells(msh)) {

        auto cl_id = offset(msh, cl);
        auto& cl_cut_type = msh.cut_cell_data[cl_id].cut;
        if (!is_cut(msh, cl)) {
            cl_cut_type = cut_type::T_OK;
            continue;
        }

        // criterion on the measure of the parts. ProtoN flags the cells before
        // refining the interface: with the polyline and polynomial interfaces
        // of 2D meshes, the parts are those of the cell cut by the chord.
        auto part_measure = [&](location where) {
            if constexpr (MeshType::dimension == 2) {
                if (msh.geometry.kind != interface_kind::implicit) {
                    typename MeshType::coordinate_type area = 0;
                    cutgeo::polygon_centroid(cutgeo::side_polygon(cut_info(msh, cl).poly, where, true), &area);
                    return area;
                }
            }
            return measure(msh, cl, where);
        };
        if (part_measure(location::IN_NEGATIVE_SIDE) < threshold_cells*measure(msh, cl)) {
            cl_cut_type = cut_type::T_KO_NEG;
            continue;
        }
        else if (part_measure(location::IN_POSITIVE_SIDE) < threshold_cells*measure(msh, cl)) {
            cl_cut_type = cut_type::T_KO_POS;
            continue;
        }
        cl_cut_type = cut_type::T_OK;

        if constexpr (MeshType::dimension == 2) {
            auto fcs = ccw_faces(msh, cl);
            auto pts = ccw_points(msh, cl);
            auto nds = ccw_nodes(msh, cl);

            if (fcs.size() != 4)
                continue;

            auto intersection = [&](size_t f) {
                return cut_info(msh, fcs[f]).intersection_point;
            };

            /* If it is a quadrilateral we have 6 possible configurations of the
             * element-cut intersection. */

            // the two cut faces share the vertex n
            auto agglo_set_single_node = [&](size_t n) -> void {
                auto f1 = (n == 0) ? fcs.size()-1 : n-1;
                auto f2 = n;

                auto da = (pts[n] - intersection(f1)).to_vector().norm() / measure(msh, fcs[f1]);
                auto db = (pts[n] - intersection(f2)).to_vector().norm() / measure(msh, fcs[f2]);

                assert(da >= 0 && da <= 1);
                assert(db >= 0 && db <= 1);

                if (std::min(da, db) > threshold) {
                    cl_cut_type = cut_type::T_OK;
                    return;
                }

                if (locate(msh, nds[n]) == location::IN_NEGATIVE_SIDE)
                    cl_cut_type = cut_type::T_KO_NEG;
                else
                    cl_cut_type = cut_type::T_KO_POS;
            };

            // the two cut faces f1 and f2 are opposite
            auto agglo_set_double_node = [&](size_t f1, size_t f2) -> void {
                assert( (f1 == 0 && f2 == 2) || (f1 == 1 && f2 == 3) );

                auto n1 = f1;
                auto n2 = (f2+1) % fcs.size();

                auto da = (pts[n1] - intersection(f1)).to_vector().norm() / measure(msh, fcs[f1]);
                auto db = (pts[n2] - intersection(f2)).to_vector().norm() / measure(msh, fcs[f2]);

                auto m1 = std::max(da, db);
                auto m2 = std::max(1-da, 1-db);

                if (std::min(m1, m2) > threshold) {
                    cl_cut_type = cut_type::T_OK;
                    return;
                }

                if (locate(msh, nds[n1]) == location::IN_NEGATIVE_SIDE)
                    cl_cut_type = (m1 <= threshold) ? cut_type::T_KO_NEG : cut_type::T_KO_POS;
                else
                    cl_cut_type = (m2 <= threshold) ? cut_type::T_KO_NEG : cut_type::T_KO_POS;
            };

            for (size_t i = 0; i < fcs.size(); i++) {
                auto f1 = i;
                auto f2 = (i+1) % fcs.size();
                auto n  = (i+1) % fcs.size();
                if (is_cut(msh, fcs[f1]) && is_cut(msh, fcs[f2]))
                    agglo_set_single_node(n);
            }

            if (is_cut(msh, fcs[0]) && is_cut(msh, fcs[2]))
                agglo_set_double_node(0, 2);

            if (is_cut(msh, fcs[1]) && is_cut(msh, fcs[3]))
                agglo_set_double_node(1, 3);
        }
    }
}

/* Delta(T): the cells sharing a face with T (f_neighbors) and the other cells
 * sharing a vertex with T (d_neighbors) */
template<CutMesh Mesh>
void
make_neighbors_info(Mesh& msh) {

    std::vector<std::vector<size_t>> face_cells(msh.faces_size());
    for (auto& cl : msh) {
        auto cl_id = offset(msh, cl);
        for (auto& fc : faces(msh, cl))
            face_cells[offset(msh, fc)].push_back(cl_id);
    }

    mesh_connectivity conn(msh.background());
    for (auto& cl : msh) {
        auto cl_id = offset(msh, cl);
        auto& ccd = msh.cut_cell_data[cl_id];
        ccd.f_neighbors.clear();
        ccd.d_neighbors.clear();
        for (auto& fc : faces(msh, cl))
            for (auto c : face_cells[offset(msh, fc)])
                if (c != cl_id)
                    ccd.f_neighbors.insert(c);
        for (auto& n : conn.neighbours_offsets(msh.background(), cl))
            if (ccd.f_neighbors.find(n.first) == ccd.f_neighbors.end())
                ccd.d_neighbors.insert(n.first);
    }
}

/* Build the description of the cut of the mesh by the level set: location of
 * the nodes, faces and cells, geometry of the cut elements, flagging of the
 * ill-cut cells (parameter threshold), neighbors */
template<CutMesh Mesh>
cut_mesh_statistics
make_cut_mesh(Mesh& msh, std::shared_ptr<const typename Mesh::level_set_type> ls,
              const cut_geometry_options& geometry, typename Mesh::coordinate_type threshold)
{
    cut_mesh_statistics stats;
    msh.level_set_function = ls;
    msh.geometry = geometry;
    if (Mesh::dimension == 3 && geometry.kind == interface_kind::polynomial)
        msh.geometry.kind = interface_kind::implicit;
    msh.quadrature_cache.clear();

    detect_node_position(msh);
    detect_cut_faces(msh, stats);
    detect_cut_cells(msh, stats);
    detect_cut_type(msh, threshold);
    make_neighbors_info(msh);
    return stats;
}

/***************************************************************************/
/*                    PAIRING OF THE ILL-CUT CELLS                          */
/***************************************************************************/

////// find_good_neighbor
template<CutMesh MeshType>
typename MeshType::cell_type
find_good_neighbor(const MeshType& msh, const typename MeshType::cell_type cl, const location where,
                   bool largest_subcell = false) {

    typename MeshType::cell_type best_neigh = cl;

    location other_where;
    if (where == location::IN_NEGATIVE_SIDE)
        other_where = location::IN_POSITIVE_SIDE;
    else
        other_where = location::IN_NEGATIVE_SIDE;

    // look for a neighboring cell to merge with
    typename MeshType::coordinate_type area = 1000000;

    // exact ties (symmetric configurations) are broken in favor of the first
    // candidate, independently of the round-off errors
    auto tol = 1e-12 * measure(msh, cl);

    // measure of the part of a neighbor in other_where. As in ProtoN, the
    // measure of an uncut neighbor is its whole measure, whatever its side,
    // so that the cut neighbors are preferred; with largest_subcell, an uncut
    // neighbor lying in "where" has no part in other_where and is preferred.
    auto other_measure = [&](const typename MeshType::cell_type& cl_n) {
        if (largest_subcell && !is_cut(msh, cl_n))
            return (locate(msh, cl_n) == where) ? typename MeshType::coordinate_type(0) : measure(msh, cl_n);
        return measure(msh, cl_n, other_where);
    };

    const auto& f_neigh = cut_info(msh, cl).f_neighbors;

    for (auto it = f_neigh.begin(); it != f_neigh.end(); ++it) {
        auto cl_n = msh[*it];

        // if cl_n is on the wrong size -> do not consider it
        if (locate(msh, cl_n) != where && locate(msh, cl_n) != location::ON_INTERFACE)
            continue;

        // if cl_n is a small cut of the same size -> do not consider it
        auto cl_n_cut_type = cut_info(msh, cl_n).cut;
        if (where == location::IN_NEGATIVE_SIDE && cl_n_cut_type == cut_type::T_KO_NEG)
            continue;
        if (where == location::IN_POSITIVE_SIDE && cl_n_cut_type == cut_type::T_KO_POS)
            continue;

        // search for the "best" neighbor -> the one with the smallest volume in other_where
        if( other_measure(cl_n) < area - tol )
        {
            area = other_measure(cl_n);
            best_neigh = cl_n;
        }
    }

    if(best_neigh == cl) // no possible face agglomerations
    {
        // look for a diagonal agglomeration
        const auto& d_neigh = cut_info(msh, cl).d_neighbors;

        for (auto it = d_neigh.begin(); it != d_neigh.end(); ++it) {
            auto cl_n = msh[*it];

            // if cl_n is on the wrong size or cut -> do not consider it
            if (locate(msh, cl_n) != where)
                continue;

            // search for the "best" neighbor -> the one with the smallest volume
            if( other_measure(cl_n) < area - tol ) {
                area = other_measure(cl_n);
                best_neigh = cl_n;
            }
        }
    }

    if(best_neigh == cl) {
        // coarse meshes (all the diagonal neighbors are cut): a diagonal
        // neighbor which is not ill-cut on the side "where"
        const auto& d_neigh = cut_info(msh, cl).d_neighbors;
        for (auto it = d_neigh.begin(); it != d_neigh.end(); ++it) {
            auto cl_n = msh[*it];
            if (!is_cut(msh, cl_n))
                continue;
            auto cl_n_cut_type = cut_info(msh, cl_n).cut;
            if (where == location::IN_NEGATIVE_SIDE && cl_n_cut_type == cut_type::T_KO_NEG)
                continue;
            if (where == location::IN_POSITIVE_SIDE && cl_n_cut_type == cut_type::T_KO_POS)
                continue;
            if( other_measure(cl_n) < area - tol ) {
                area = other_measure(cl_n);
                best_neigh = cl_n;
            }
        }
    }

    if(best_neigh == cl) {
        throw std::logic_error("No possible agglomerations !!!");
    }

    return best_neigh;
}

template<typename Mesh>
void
output_agglo_lists(Mesh& msh, std::vector<int> table_neg, std::vector<int> table_pos,
                   std::string file)
{
    // number of arrows
    size_t nb_arrows = 0;
    for(size_t i=0; i<table_neg.size(); i++)
    {
        if(table_neg.at(i) != -1)
            nb_arrows++;
        if(table_pos.at(i) != -1)
            nb_arrows++;
    }

    // initiate the output file
    std::ofstream output(file, std::ios::out | std::ios::trunc);
    if( !output )
        std::cerr << "agglo output file has not been opened" << std::endl;

    output << "5 " << nb_arrows << " 12" << std::endl;
    output << "x" << std::endl;
    output << "y" << std::endl;
    output << "z" << std::endl;
    output << "u" << std::endl;
    output << "v" << std::endl;

    output << "0.  1.  10" << std::endl;
    output << "0.  1.  10" << std::endl;
    output << "0.  0.  10" << std::endl;
    output << "-1  1.  10" << std::endl;
    output << "-1  1.  10" << std::endl;


    // loop on the cells
    size_t cp = 0;
    for (auto& cl : cells(msh)) {
        int TN = table_neg.at(cp);
        if (TN != -1) {
            auto bar_cl = barycenter(msh,cl);
            auto bar_neigh = barycenter(msh,msh[TN]);

            auto vect = bar_neigh - bar_cl;

            output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                   << vect[0] << "   " << vect[1] << std::endl;
        }


        int TP = table_pos.at(cp);
        if( TP != -1)
        {
            auto bar_cl = barycenter(msh,cl);
            auto bar_neigh = barycenter(msh,msh[TP]);

            auto vect = bar_neigh - bar_cl;

            output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                   << vect[0] << "   " << vect[1] << std::endl;
        }
        cp++;
    }

    output.close();
}

template<typename Mesh>
void
output_agglo_lists_step4(Mesh& msh, std::vector<int> table_neg, std::vector<int> table_pos,
                   std::string file) {
    // number of arrows
    size_t nb_arrows = 0;
    for (size_t i=0; i<table_neg.size(); i++) {
      if(table_neg.at(i) != -1)
        nb_arrows++;
    }
    for (size_t i=0; i<table_pos.size(); i++) {
      if(table_pos.at(i) != -1)
        nb_arrows++;
    }

    // initiate the output file
    std::ofstream output(file, std::ios::out | std::ios::trunc);
    if( !output )
        std::cerr << "agglo output file has not been opened" << std::endl;

    output << "5 " << nb_arrows << " 12" << std::endl;
    output << "x" << std::endl;
    output << "y" << std::endl;
    output << "z" << std::endl;
    output << "u" << std::endl;
    output << "v" << std::endl;

    output << "0.  1.  10" << std::endl;
    output << "0.  1.  10" << std::endl;
    output << "0.  0.  10" << std::endl;
    output << "-1  1.  10" << std::endl;
    output << "-1  1.  10" << std::endl;


    // loop on the cells
    size_t cp = 0;
    for (auto& cl : cells(msh)) {
        int TN = table_neg.at(cp);
        if( TN != -1) {
          auto h = diameter(msh, cl);
          auto bar_cl = barycenter(msh, cl);
          auto bar_neigh = barycenter(msh, msh[TN]);
          auto vect = bar_neigh - bar_cl;
          // VERS LE HAUT
          if (vect[0] <= 1e-5 && vect[1] > 0) {
            bar_cl[0] = bar_cl[0] - h/6.0;
          }
          // VERS LE BAS
          if (vect[0] <= 1e-5 && vect[1] < 0) {
            // EN HAUT A GAUCHE
            if (bar_cl[0] <= 0.5 && bar_cl[1] >= 0.5)
                bar_cl[0] = bar_cl[0] + h/6.0;
            else
                bar_cl[0] = bar_cl[0] - h/6.0;
          }
          // VERS LA DROITE
          else if (vect[1] <= 1e-5 && vect[0] > 0){
            // EN BAS A GAUCHE
            if (bar_cl[0] <= 0.5 && bar_cl[1] <= 0.5)
                bar_cl[1] = bar_cl[1] + h/6.0;
            // ELSE
            else
                bar_cl[1] = bar_cl[1] - h/6.0;
          }
          // VERS LA GAUCHE
          else if (vect[1] <= 1e-5 && vect[0] < 0){
            // EN HAUT A GAUCHE
            if (bar_cl[0] >= 0.5 && bar_cl[1] <= 0.5)
                bar_cl[1] = bar_cl[1] + h/6.0;
            else
                bar_cl[1] = bar_cl[1] + h/6.0;
          }
          output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                 << vect[0] << "   " << vect[1] << std::endl;
        }
        int TP = table_pos.at(cp);
        if( TP != -1) {
          auto bar_cl = barycenter(msh,cl);
          auto bar_neigh = barycenter(msh, msh[TP]);
          auto vect = bar_neigh - bar_cl;
          output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                 << vect[0] << "   " << vect[1] << std::endl;
        }
        cp++;
    }
    output.close();
}

/* Steps 1 to 3 of the pairing/agglomeration procedure ([12, Algorithm 1]):
 * each ill-cut side points towards a neighbor (table_neg/table_pos give the
 * target of each cell, -1 if none), chains being avoided for the cells of
 * domain 1. */
template<CutMesh Mesh>
void
make_pairing_tables(Mesh& msh, std::vector<int>& table_neg, std::vector<int>& table_pos, bool output_pairing = false,
                    bool largest_subcell = false) {

    // Initiate lists to store the pairing infos
    size_t nb_cells = msh.cells_size();
    table_neg.assign(nb_cells, -1);
    table_pos.assign(nb_cells, -1);

    // the targets (for the test "already agglomerated" of step 1)
    std::vector<bool> targeted(nb_cells, false);

    // start the process for domain 1, and then domain 2
    for(size_t domain=1; domain < 3; domain++) {

        // loop on the cells
        for (auto& cl : cells(msh)) {
            auto cl_id = offset(msh, cl);
            auto cl_cut_type = msh.cut_cell_data[cl_id].cut;
            location where;
            if(domain == 1)
                where = location::IN_NEGATIVE_SIDE;
            else
                where = location::IN_POSITIVE_SIDE;

            if(cl_cut_type == cut_type::T_OK)
                continue;
            else if(cl_cut_type == cut_type::UNDEF)
                throw std::logic_error("UNDEF agglo_set");
            else if(cl_cut_type == cut_type::T_KO_NEG
                    && where != location::IN_NEGATIVE_SIDE)
                continue;
            else if(cl_cut_type == cut_type::T_KO_POS
                    && where != location::IN_POSITIVE_SIDE)
                continue;

            // if cl is already agglomerated : no need for further agglomerations
            if (targeted[cl_id])
                continue;

            typename Mesh::cell_type best_neigh = find_good_neighbor(msh, cl, where, largest_subcell);

            // prepare agglomeration of cells cl and best_neigh
            size_t offset1 = cl_id;
            size_t offset2 = offset(msh,best_neigh);

            if(where == location::IN_NEGATIVE_SIDE)
                table_neg.at(offset1) = offset2;
            else
                table_pos.at(offset1) = offset2;
            targeted[offset2] = true;
        }

        if(output_pairing && domain == 1)
            output_agglo_lists(msh, table_neg, table_pos, "agglo_one.okc");
        if(output_pairing && domain == 2)
            output_agglo_lists(msh, table_neg, table_pos, "agglo_two.okc");
    }

    //////////////   CHANGE THE AGGLO FOR THE CELLS OF DOMAIN 1 THAT ARE TARGETTED ///////
    std::vector<int> pos_source(nb_cells, -1);  // first cell targeting each cell in table_pos
    for (size_t i = 0; i < nb_cells; i++)
        if (table_pos[i] != -1 && pos_source[table_pos[i]] == -1)
            pos_source[table_pos[i]] = i;
    std::vector<size_t> neg_count(nb_cells, 0); // number of cells targeting each cell in table_neg
    for (size_t i = 0; i < nb_cells; i++)
        if (table_neg[i] != -1)
            neg_count[table_neg[i]]++;

    for (auto& cl : cells(msh)) {

        auto cl_id = offset(msh, cl);
        if(msh.cut_cell_data[cl_id].cut != cut_type::T_KO_NEG)
            continue;

        size_t offset1 = cl_id;

        // are there cells that try to agglomerate with cl ?
        if (pos_source[offset1] == -1)
            continue;

        // at this point a cell tries to agglomerate with cl
        int cl1_agglo = table_neg.at(offset1);
        if (cl1_agglo == -1)
            continue;

        // -> check that no one else tries to agglomerate with cl1_agglo
        bool agglo = (neg_count[cl1_agglo] > 1);

        auto cl1_cut_type = msh.cut_cell_data.at(cl1_agglo).cut;
        if (!agglo && cl1_cut_type == cut_type::T_KO_POS) {
            continue;
        }

        // at this point we risk chain agglomerations
        // -> remove the target of cl
        neg_count[cl1_agglo]--;
        table_neg.at(offset1) = -1;
    }
    if (output_pairing)
        output_agglo_lists(msh, table_neg, table_pos, "agglo_three.okc");
}

template<CutMesh Mesh>
void
make_polynomial_extension(Mesh& msh, bool output_pairing = false, bool largest_subcell = false) {

    // Reset any previous pairing
    for (auto& ccd : msh.cut_cell_data) {
        ccd.stabilizing_cell.reset();
        ccd.dependent_ill_cut_neg_cells.clear();
        ccd.dependent_ill_cut_pos_cells.clear();
    }

    std::vector<int> table_neg, table_pos;
    size_t nb_cells = msh.cells_size();
    make_pairing_tables(msh, table_neg, table_pos, output_pairing, largest_subcell);

    ///////////////////////////////////////////////////////////////////////////// STEP 4
    // ALL BAD CUT CELLS MUST POINT TOWARDS A CELL
    std::vector<int> first_pos(nb_cells, -1), first_neg(nb_cells, -1);
    for (size_t i = 0; i < nb_cells; i++) {
        if (table_pos[i] != -1 && first_pos[table_pos[i]] == -1)
            first_pos[table_pos[i]] = i;
        if (table_neg[i] != -1 && first_neg[table_neg[i]] == -1)
            first_neg[table_neg[i]] = i;
    }

    for (auto& cl : cells(msh)) {
        auto offset_cl = offset(msh, cl);
        auto cl_cut_type = msh.cut_cell_data[offset_cl].cut;
        auto TN = table_neg.at(offset_cl);
        auto TP = table_pos.at(offset_cl);
        if (TN == -1 && cl_cut_type == cut_type::T_KO_NEG && first_pos[offset_cl] != -1)
            table_neg.at(offset_cl) = first_pos[offset_cl];
        if (TP == -1 && cl_cut_type == cut_type::T_KO_POS && first_neg[offset_cl] != -1)
            table_pos.at(offset_cl) = first_neg[offset_cl];
    }

    // FILLING THE STRUCTURES paired_cells / dependent_cells_neg / dependent_cells_pos / paire(T,i)
    for (auto& cl : cells(msh)) {
        auto offset_cl = offset(msh, cl);
        auto& ccd = msh.cut_cell_data[offset_cl];
        auto cl_cut_type = ccd.cut;
        if (cl_cut_type != cut_type::T_KO_NEG && cl_cut_type != cut_type::T_KO_POS)
            continue;
        if (table_neg.at(offset_cl) != -1)
            ccd.stabilizing_cell = table_neg.at(offset_cl);
        if (table_pos.at(offset_cl) != -1)
            ccd.stabilizing_cell = table_pos.at(offset_cl);
        if (!ccd.stabilizing_cell.has_value())
            throw std::logic_error("make_polynomial_extension: ill-cut cell without stabilizing cell");
        auto& good = msh.cut_cell_data[ccd.stabilizing_cell.value()];
        if (cl_cut_type == cut_type::T_KO_NEG)
            good.dependent_ill_cut_neg_cells.push_back(offset_cl);
        else
            good.dependent_ill_cut_pos_cells.push_back(offset_cl);
    }

    // Display of the arrows
    std::vector<int> table(nb_cells, -1);
    if (output_pairing) {
        output_agglo_lists_step4(msh, table_neg, table, "agglo_four.okc");
        output_agglo_lists_step4(msh, table, table_pos, "agglo_five.okc");
    }
}

/* A mesh subcell (T,i) of the paper: the part of the cell "cell" lying in the
 * subdomain "side", together with the ill-cut cells S whose side "side" is
 * stabilized by (T,i), i.e. S in N_i^{-1}(T). */
struct cut_subcell {
    size_t              cell;
    location            side;
    std::vector<size_t> dependent_cells;
};

/* Cell agglomeration of [12] (ProtoN's make_agglomeration): every ill-cut
 * cell is merged with its target (see make_pairing_tables). The mesh is not
 * modified: an agglomerated cell is described by the list of its cells. */
struct cut_agglomeration {
    std::vector<std::vector<size_t>> groups;        // cells of each agglomerated cell
    std::vector<int>                 group_of_cell; // -1 if the cell is not agglomerated
};

template<CutMesh Mesh>
cut_agglomeration
make_cell_agglomeration(Mesh& msh, bool output_agglomeration = false, bool largest_subcell = false) {

    size_t nb_cells = msh.cells_size();

    std::vector<int> table_neg, table_pos;
    make_pairing_tables(msh, table_neg, table_pos, output_agglomeration, largest_subcell);

    std::vector<std::vector<size_t>> cell_faces(nb_cells);
    for (auto& cl : msh) {
        auto c = offset(msh, cl);
        for (auto& fc : faces(msh, cl))
            cell_faces[c].push_back(offset(msh, fc));
        std::sort(cell_faces[c].begin(), cell_faces[c].end());
    }
    /* do the interfaces of the cut cells a and b meet on a cut face? */
    auto connected = [&](size_t a, size_t b) {
        std::vector<size_t> common;
        std::set_intersection(cell_faces[a].begin(), cell_faces[a].end(), cell_faces[b].begin(),
                              cell_faces[b].end(), std::back_inserter(common));
        for (auto f : common)
            if (msh.cut_face_data[f].loc == location::ON_INTERFACE)
                return true;
        return false;
    };
    auto is_cut_cell = [&](size_t c) { return msh.cut_cell_data[c].loc == location::ON_INTERFACE; };

    std::vector<std::vector<size_t>> groups;
    std::vector<int> group_of_cell(nb_cells, -1);

    auto add_to_group = [&](size_t g, size_t c) {
        groups[g].push_back(c);
        group_of_cell[c] = g;
    };

    // another cell is needed when the interface of the agglomerated cell
    // would not be connected (square interface): a cut cell connecting the
    // interfaces of g and c
    auto corner_cell = [&](size_t g, size_t c) -> int {
        if (!is_cut_cell(c))
            return -1;
        bool group_cut = false;
        for (auto m : groups[g]) {
            if (!is_cut_cell(m))
                continue;
            group_cut = true;
            if (connected(m, c))
                return -1;
        }
        if (!group_cut)
            return -1;
        for (auto n : msh.cut_cell_data[c].f_neighbors) {
            if (!is_cut_cell(n) || !connected(n, c))
                continue;
            for (auto m : groups[g])
                if (is_cut_cell(m) && connected(n, m))
                    return int(n);
        }
        /* the interfaces only meet on an edge or a vertex (3D): the integrals
         * of the agglomerated cell are sums over its cells, which do not need
         * a connected interface */
        return -1;
    };

    for (size_t c = 0; c < nb_cells; c++) {
        int target = -1;
        if (table_pos.at(c) != -1) target = table_pos.at(c);
        if (table_neg.at(c) != -1) target = table_neg.at(c);
        if (target == -1)
            continue;

        int gc = group_of_cell[c], gt = group_of_cell[target];
        if (gc != -1 && gt != -1) {
            /* both cells already agglomerated (3D): the groups are merged */
            if (gc != gt) {
                for (auto m : groups[gt])
                    add_to_group(gc, m);
                groups[gt].clear();
            }
            continue;
        }
        if (gc == -1 && gt == -1) {
            groups.push_back({});
            add_to_group(groups.size()-1, c);
            if (is_cut_cell(c) && is_cut_cell(target) && !connected(c, target)) {
                int corner = corner_cell(groups.size()-1, target);
                if (corner != -1 && group_of_cell[corner] == -1)
                    add_to_group(groups.size()-1, corner);
            }
            add_to_group(groups.size()-1, target);
            continue;
        }

        // only one of the two cells is already agglomerated
        size_t g = (gc != -1) ? gc : gt;
        size_t c2 = (gc != -1) ? target : c;
        int corner = corner_cell(g, c2);
        if (corner != -1 && group_of_cell[corner] == -1)
            add_to_group(g, corner);
        add_to_group(g, c2);
    }

    cut_agglomeration ret;
    ret.group_of_cell.assign(nb_cells, -1);
    for (auto& g : groups) {
        if (g.empty())
            continue;
        for (auto m : g)
            ret.group_of_cell[m] = ret.groups.size();
        ret.groups.push_back(g);
    }

    return ret;
}

/* Build the partition P_h = P_h^OK U P_h^KO of the mesh subcells (eq. (22)).
 * Must be called after make_polynomial_extension(). */
template<CutMesh Mesh>
std::pair<std::vector<cut_subcell>, std::vector<cut_subcell>>
make_subcell_partition(const Mesh& msh) {

    std::vector<cut_subcell> P_OK, P_KO;

    for (auto& cl : msh) {
        auto cl_id = offset(msh, cl);
        const auto& ccd = msh.cut_cell_data[cl_id];
        const auto& dp_neg = ccd.dependent_ill_cut_neg_cells;
        const auto& dp_pos = ccd.dependent_ill_cut_pos_cells;

        if (ccd.loc != location::ON_INTERFACE) {
            if (ccd.loc == location::IN_NEGATIVE_SIDE)
                P_OK.push_back( {cl_id, location::IN_NEGATIVE_SIDE, dp_neg} );
            else
                P_OK.push_back( {cl_id, location::IN_POSITIVE_SIDE, dp_pos} );
        }
        else if (ccd.cut == cut_type::T_OK) {
            P_OK.push_back( {cl_id, location::IN_NEGATIVE_SIDE, dp_neg} );
            P_OK.push_back( {cl_id, location::IN_POSITIVE_SIDE, dp_pos} );
        }
        else if (ccd.cut == cut_type::T_KO_NEG) {
            P_KO.push_back( {cl_id, location::IN_NEGATIVE_SIDE, {}} );
            P_OK.push_back( {cl_id, location::IN_POSITIVE_SIDE, dp_pos} );
        }
        else if (ccd.cut == cut_type::T_KO_POS) {
            P_OK.push_back( {cl_id, location::IN_NEGATIVE_SIDE, dp_neg} );
            P_KO.push_back( {cl_id, location::IN_POSITIVE_SIDE, {}} );
        }
        else
            throw std::logic_error("make_subcell_partition: cut type undefined");
    }

    return std::make_pair(P_OK, P_KO);
}

} // namespace disk
