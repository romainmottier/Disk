
#pragma once

#include "mesh.hpp"
#include "level_set.hpp"
#include "mesh_storage.hpp"
#include "diskpp/quadratures/quadratures.hpp"

namespace disk {

template<typename Mesh, typename Element>
std::vector<typename Mesh::node_type>
nodes(const Mesh& msh, const Element& elem)
{
    auto ptids = elem.point_ids();

    auto storage = msh.backend_storage();

    auto points_begin = msh.points_begin();
    auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto {
        return *std::next(storage->nodes.begin(), pi);
    };

    std::vector<typename Mesh::node_type> nds(ptids.size());
    std::transform(ptids.begin(), ptids.end(), nds.begin(), ptid_to_point);

    return nds;
}

template<mesh_2D MeshType>
location
locate(const MeshType& msh, const typename MeshType::cell_type& cl) {
    
    auto storage = msh.backend_storage();
    auto cl_id = offset(msh, cl);
    auto& cl_loc = storage -> cut_cell_data[cl_id].loc;
    
    assert(cl_loc != location::UNDEF);

    return cl_loc;

}

template<mesh_2D MeshType>
location
locate(const MeshType& msh, const typename MeshType::face_type& fc) {

    auto storage = msh.backend_storage();
    auto fc_id = offset(msh, fc);
    auto& fc_loc = storage -> cut_face_data[fc_id].loc;
    
    assert(fc_loc != location::UNDEF);

    return fc_loc;

}

template<mesh_2D MeshType>
location
locate(const MeshType& msh, const typename MeshType::node_type& nd) {
    
    auto storage = msh.backend_storage();
    auto nd_id = nd.point_ids()[0];
    auto& nd_loc = storage -> cut_node_data[nd_id].loc;
    
    assert(nd_loc != location::UNDEF);

    return nd_loc;

}

template<mesh_2D MeshType>
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

template<mesh_2D MeshType>
std::vector<temp_tri<typename MeshType::coordinate_type >>
triangulate(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    assert(is_cut(msh, cl));

    auto tp = collect_triangulation_points(msh, cl, where);
    // auto bar = barycenter(tp);
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

template<mesh_2D MeshType>
std::vector<typename MeshType::point_type>
collect_triangulation_points(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    typedef typename MeshType::point_type point_type;
    typedef typename MeshType::node_type  node_type;

    assert(is_cut(msh, cl));
    auto ns = nodes(msh, cl);

    std::vector<point_type> ret;

    auto node2pt = [&](const MeshType& msh, const node_type& n) -> auto {
        return points(msh, n);
    };

    auto insert_interface = [&](void) -> void {
        auto storage = msh.backend_storage();
        auto cl_id = offset(msh,cl);
        auto& cl_interface = storage -> cut_cell_data[cl_id].interface;
        if (where == location::IN_NEGATIVE_SIDE)
            ret.insert(ret.end(), cl_interface.begin(), cl_interface.end());
        else if (where == location::IN_POSITIVE_SIDE)
            ret.insert(ret.end(), cl_interface.rbegin(), cl_interface.rend());
        else
            throw std::logic_error("If you've got here there is some issue...");
    };

    bool case1 = locate(msh, ns.front()) == where && locate(msh, ns.back()) != where;
    bool case2 = locate(msh, ns.front()) != where && locate(msh, ns.back()) == where;
    bool case3 = locate(msh, ns.front()) != where && locate(msh, ns.back()) != where;
    //bool case4 = locate(msh, ns.front()) == where && locate(msh, ns.back()) == where;

    if ( case1 || case2 || case3 ) {
        for (size_t i = 0; i < ns.size(); i++)
            if ( locate(msh, ns[i]) == where )
                ret.push_back(points(msh, ns[i])[0]);

        insert_interface();
    }
    else  {
        size_t i = 0;
        while ( i < ns.size() && locate(msh, ns.at(i)) == where )
            ret.push_back(points(msh, ns[i++])[0]);
        insert_interface();
        while ( i < ns.size() && locate(msh, ns.at(i)) != where )
            i++;
        while ( i < ns.size() && locate(msh, ns.at(i)) == where )
            ret.push_back(points(msh, ns[i++])[0]);
    }

    return ret;
}

template<mesh_2D MeshType>
typename MeshType::point_type
tesselation_center(const MeshType& msh, const typename MeshType::cell_type& cl, location where) {

    auto fcs = faces(msh, cl);
    auto pts = points(msh, cl);
    auto nds = nodes(msh, cl);

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


template<mesh_2D MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, 2>>
integrate(const MeshType& msh, const typename MeshType::cell_type& cl, size_t degree, const location& where) {

    auto storage = msh.backend_storage();

    auto cl_id = offset(msh, cl);
    auto& cl_inte_n = storage -> cut_cell_data[cl_id].integration_n;
    auto& cl_inte_p = storage -> cut_cell_data[cl_id].integration_p;

    if (cl_inte_n.size() != 0 && where == location::IN_NEGATIVE_SIDE)
        return cl_inte_n;

    if(cl_inte_p.size() != 0 && where == location::IN_POSITIVE_SIDE)
        return cl_inte_p;

    return make_integrate(msh, cl, degree, where);

}

template<mesh_2D MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, 2>>
integrate(const MeshType& msh, const typename MeshType::face_type& fc, size_t degree, const location& where) {

    using T = typename MeshType::coordinate_type;

    std::vector<disk::quadrature_point<T, 2>> ret;
    if (locate(msh, fc) != where && locate(msh, fc) != location::ON_INTERFACE)
        return ret;

    if (!is_cut(msh, fc)) /* Element is not cut, use std. integration */
        return integrate(msh, fc, degree);

    auto pts = points(msh, fc, where);

    auto scale = pts[1] - pts[0];
    auto meas = scale.to_vector().norm();

    auto qps = edge_quadrature<T>(degree);

    for (auto itor = qps.begin(); itor != qps.end(); itor++) {
        auto qp = *itor;
        //auto p = qp.first.x() * scale + pts[0];
        auto t = qp.first.x();
        auto p = 0.5*(1-t)*pts[0] + 0.5*(1+t)*pts[1];
        auto w = qp.second * meas * 0.5;

        ret.push_back( std::make_pair(p, w) );
    }

    return ret;
}

template<mesh_2D MeshType>
std::vector<disk::quadrature_point<typename MeshType::coordinate_type, 2>>
make_integrate(const MeshType& msh, const typename MeshType::cell_type& cl, size_t degree, location where) {

    using T = typename MeshType::coordinate_type;

    std::vector<disk::quadrature_point<T, 2>> ret;

    if (locate(msh, cl) != where && locate(msh, cl) != location::ON_INTERFACE )
        return ret;

    if (!is_cut(msh, cl)) // Element is not cut, use std. integration 
        return integrate(msh, cl, degree);

    auto tris = triangulate(msh, cl, where);
    for (auto& tri : tris) {
        auto qpts = disk::quadrature::triangle_gauss(degree,tri.pts[0], tri.pts[1], tri.pts[2]);
        ret.insert(ret.end(), qpts.begin(), qpts.end());
    }

    return ret;
}

template<mesh_2D MeshType, typename Function>
void 
detect_node_position(const MeshType& msh, const Function& level_set_function) {

    auto storage = msh.backend_storage();
    if (!storage) {
        std::cout << "STORAGE NOT VALID";
        return;
    }
    size_t nb_nodes = storage -> nodes.size();
    storage -> cut_node_data.resize(nb_nodes);
    for (size_t i=0; i < nb_nodes; i++) {
        auto& nd = storage->nodes[i];
        auto pt = points(msh, nd);
        auto& node_loc = storage -> cut_node_data[i].loc;
        if (level_set_function(pt[0]) < 0 )
            node_loc = location::IN_NEGATIVE_SIDE;
        else
            node_loc = location::IN_POSITIVE_SIDE;
    }
}

template<mesh_2D MeshType , typename Function>
void 
detect_cut_faces(const MeshType& msh, const Function& level_set_function) {

    auto storage = msh.backend_storage();
    size_t nb_faces = storage -> edges.size();
    storage -> cut_face_data.resize(nb_faces);

    for (auto& fc : faces(msh)) {
        auto fc_id = offset(msh, fc);
        auto pts = points(msh, fc);
        auto l0 = level_set_function(pts[0]);
        auto l1 = level_set_function(pts[1]);
        auto& face_loc = storage -> cut_face_data[fc_id].loc;
        auto& face_node_inside = storage -> cut_face_data[fc_id].node_inside;
        auto& face_intersection_point = storage -> cut_face_data[fc_id].intersection_point;
        if (l0 >= 0 && l1 >= 0) {
           face_loc  = location::IN_POSITIVE_SIDE;
            continue;
        }
        if (l0 < 0 && l1 < 0) {
            face_loc = location::IN_NEGATIVE_SIDE;
            continue;
        }

        auto threshold = diameter(msh, fc) / 1e20;
        auto pm = find_zero_crossing(pts[0], pts[1], level_set_function, threshold);
        face_node_inside = ( l0 < 0 ) ? 0 : 1;
        face_loc = location::ON_INTERFACE;
        face_intersection_point = pm;
    }
}

template<mesh_2D MeshType, typename Function>
void
detect_cut_cells(MeshType& msh, const Function& level_set_function) {

    using T = typename MeshType::coordinate_type;
    typedef typename MeshType::face_type  face_type;
    typedef typename MeshType::point_type point_type;

    auto storage = msh.backend_storage();
    auto& cut_cell_data = storage->cut_cell_data;
    size_t nb_cells = storage -> surfaces.size();
    cut_cell_data.resize(nb_cells);
    size_t cell_i = 0;
    for (auto& cl : cells(msh)) {

        auto cl_id = offset(msh, cl);
        auto fcs = faces(msh, cl);
        std::array<std::pair<size_t, point_type>, 2>  cut_faces;

        size_t k = 0;
        for (size_t i = 0; i < fcs.size(); i++) {
            bool face_is_cut_Q = is_cut(msh, fcs[i]);
            auto fc_id = offset(msh,fcs[i]);
            auto& cut_face_data = storage -> cut_face_data[fc_id];
            auto& face_intersection_pt = cut_face_data.intersection_point;
            if (face_is_cut_Q)
                cut_faces.at(k++) = std::make_pair(i, face_intersection_pt);
        }

        /* If a face is cut, the cells that own the face are cut. Is this
         * unconditionally true? It should...fortunately this isn't avionics
         * software */

        auto& cell_loc = storage -> cut_cell_data[cl_id].loc;
        auto& cell_p0 = storage -> cut_cell_data[cl_id].p0;
        auto& cell_p1 = storage -> cut_cell_data[cl_id].p1;
        auto& cell_interface = storage -> cut_cell_data[cl_id].interface;
        if (k == 0) {
            auto is_positive = [&](const point_type& pt) -> bool {
                return level_set_function(pt) > 0;
            };

            auto pts = points(msh, cl);
            if ( std::all_of(pts.begin(), pts.end(), is_positive) )
                cell_loc = location::IN_POSITIVE_SIDE;
            else
                cell_loc = location::IN_NEGATIVE_SIDE;
        }

        if (k == 2) {
            cell_loc = location::ON_INTERFACE;
            auto p0 = cut_faces[0].second;
            auto p1 = cut_faces[1].second;
            auto pt = p1 - p0;
            auto pt_t = point<T,2>(-pt.y(), pt.x());
            auto pn = p0 + pt_t;

            if (level_set_function(pn) >= 0) {
                cell_p0 = p1;
                cell_p1 = p0;
            }
            else {
                cell_p0 = p0;
                cell_p1 = p1;
            }

            cell_interface.push_back(cell_p0);
            cell_interface.push_back(cell_p1);
        }

        if ( k != 0 && k != 2 )
            throw std::logic_error("invalid number of cuts in cell");

        cell_i++;
    }
}

template<typename T, typename Function>
point<T, 2>
find_zero_crossing(const point<T,2>& p0, const point<T,2>& p1, const Function& level_set_function, const T& threshold) {

    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;

    /* A threshold of 1/10000 the diameter of the element is considered
     * acceptable. Since with 24 iterations we reduce the error by 16384
     * and the worst case is that the two points are at the opposite sides
     * of the element, we put 30 as limit. */
    size_t max_iter = 50;

    do {
        auto la = level_set_function(pa);
        auto lb = level_set_function(pb);
        auto lm = level_set_function(pm);

        if ( (lb >= 0 && lm >= 0) || (lb < 0 && lm < 0) ) {   /* intersection is between pa and pm */
            pm_prev = pm;
            pb = pm;
            pm = (pa+pb)/2.0;
        }
        else {   /* intersection is between pm and pb */
            pm_prev = pm;
            pa = pm;
            pm = (pa+pb)/2.0;
        }

        x_diff_sq = (pm_prev.x() - pm.x()) * (pm_prev.x() - pm.x());
        y_diff_sq = (pm_prev.y() - pm.y()) * (pm_prev.y() - pm.y());

    } while ( (sqrt(x_diff_sq + y_diff_sq) > threshold) && max_iter-- );

    return pm;

    /* Affine zero crossing was like that: */
    //auto t = l0/(l0-l1);
    //auto ip = (pts[1] - pts[0]) * t + pts[0];

}

template<mesh_2D MeshType , typename Function>
void
detect_cut_type(MeshType& msh, const Function& level_set_function) {

    auto storage = msh.backend_storage();
    const auto threshold = 0.3;
    const auto threshold_cells = 0.3;

    for (auto& cl : cells(msh)) {

        auto fcs = faces(msh, cl);
        auto pts = points(msh, cl);
        auto nds = nodes(msh, cl);

        auto cl_id = offset(msh, cl);

        if (fcs.size() != 4)
            throw std::invalid_argument("This works only on quads for now");

        auto& cell_cut = storage -> cut_cell_data[cl_id].cut;
        if (!is_cut(msh, cl)) {
            cell_cut = cut_type::T_OK;
            continue;
        }

        // another criterion on the area of the cell
        if (measure(msh, cl, location::IN_NEGATIVE_SIDE) < threshold_cells*measure(msh, cl)) {
            cell_cut = cut_type::T_KO_NEG;
            continue;
        }
        else if (measure(msh, cl, location::IN_POSITIVE_SIDE) < threshold_cells * measure(msh, cl)) {
            cell_cut = cut_type::T_KO_POS;
            continue;
        }

        /* If it is a quadrilateral we have 6 possible configurations of the
         * element-cut intersection. */

        auto agglo_set_single_node = [&](size_t f1, size_t f2, size_t n) -> void {

            auto fc_id1 = offset(msh, fcs[f1]);
            auto fc_id2 = offset(msh, fcs[f2]);
            auto& f1_intersection = storage -> cut_face_data[fc_id1].intersection_point;
            auto& f2_intersection = storage -> cut_face_data[fc_id2].intersection_point;

            auto ma = measure(msh, fcs[f1]);
            auto pa = (pts[n] - f1_intersection);
            auto da = pa.to_vector().norm() / ma;

            auto mb = measure(msh, fcs[f2]);
            auto pb = (pts[n] - f2_intersection);
            auto db = pb.to_vector().norm() / mb;

            assert(da >= 0 && da <= 1);
            assert(db >= 0 && db <= 1);

            auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
            if ( std::min(da, db) > threshold ) {
                cl_cut_type = cut_type::T_OK;
                return;
            }

            if (locate(msh, nds[n]) == location::IN_NEGATIVE_SIDE )
                cl_cut_type = cut_type::T_KO_NEG;
            else
                cl_cut_type = cut_type::T_KO_POS;

        };
        
        auto agglo_set_double_node = [&](size_t f1, size_t f2, size_t n1, size_t n2) -> void {

            assert ( (f1 == 0 && f2 == 3) || ( f1 == 1 && f2 == 2 ) );

            auto fc_id1 = offset(msh,fcs[f1]);
            auto& intersection_pt1 = storage -> cut_face_data[fc_id1].intersection_point;
            auto ma = measure(msh, fcs[f1]);
            auto pa = (pts[n1] - intersection_pt1);
            auto da = pa.to_vector().norm() / ma;

            auto fc_id2 = offset(msh,fcs[f2]);
            auto& intersection_pt2 = storage -> cut_face_data[fc_id2].intersection_point;
            auto mb = measure(msh, fcs[f2]);
            auto pb = (pts[n2] - intersection_pt2);
            auto db = pb.to_vector().norm() / mb;

            assert(da >= 0 && da <= 1);
            assert(db >= 0 && db <= 1);

            auto m1 = std::max(da, db);
            auto m2 = std::max(1-da, 1-db);

            auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
            if (std::min(m1, m2) > threshold) {
                cl_cut_type = cut_type::T_OK;
                return;
            }

            if (locate(msh, nds[n1]) == location::IN_NEGATIVE_SIDE)
                cl_cut_type = (m1 <= threshold) ? cut_type::T_KO_NEG : cut_type::T_KO_POS;
            else
                cl_cut_type = (m2 <= threshold) ? cut_type::T_KO_NEG : cut_type::T_KO_POS;
        
        };

        if (is_cut(msh, fcs[0]) && is_cut(msh, fcs[1]))
                agglo_set_single_node(0,1,0);
                
        if ((is_cut(msh, fcs[0]) && is_cut(msh, fcs[2])))
                agglo_set_single_node(0,2,1);
                
        if ((is_cut(msh, fcs[1]) && is_cut(msh, fcs[3])))
                agglo_set_single_node(1,3,2);
                
        if ((is_cut(msh, fcs[2]) && is_cut(msh, fcs[3])))
                agglo_set_single_node(2,3,3);

        if ( is_cut(msh, fcs[0]) && is_cut(msh, fcs[3]) )
            agglo_set_double_node(0,3,0,2);

        if ( is_cut(msh, fcs[1]) && is_cut(msh, fcs[2]) )
            agglo_set_double_node(1,2,0,1);
    }
}

template<typename Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::cell_type& cl) {
    return false;
}

template<mesh_2D Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::cell_type& cl) {

    auto storage = msh.backend_storage();
    auto cl_id = offset(msh, cl);
    auto& cl_loc = storage -> cut_cell_data[cl_id].loc;
    
    if (cl_loc == location::ON_INTERFACE)
        return true;
    else 
        return false;

}

template<typename Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::face_type& fc) {
    return false;

}

template<mesh_2D Mesh>
bool
is_cut(const Mesh& msh, const typename Mesh::face_type& fc) {

    auto storage = msh.backend_storage();
    auto fc_id = offset(msh, fc);
    auto& fc_loc = storage -> cut_face_data[fc_id].loc;
    
    if (fc_loc == location::ON_INTERFACE)
        return true;
    else 
        return false;

}

// //// version for cartesian meshes -> very quick
// /* this creates Delta(T) */
// // there are at least two row and two columns of cells
template<mesh_2D Mesh>
void
make_neighbors_info_cartesian(const Mesh& msh) {

    using T = typename Mesh::coordinate_type;
    auto storage = msh.backend_storage();

    // Face neighbors
    auto conn = connectivity_via_faces(msh);
    for (auto& cl : msh) {
        auto cl_id = offset(msh, cl);
        auto& cl_f_neighbors = storage -> cut_cell_data[cl_id].f_neighbors;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            auto cl_neighbour = conn.neighbour_via(msh, cl, fc);
            if (cl_neighbour.second) {
                auto cl_neighbour_id = offset(msh, cl_neighbour.first);
                cl_f_neighbors.insert(cl_neighbour_id);
            }
        }
    }

    // Diagonal neighbors 
    // TO ADD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

template<mesh_2D MeshType, typename Function>
void
refine_interface(MeshType& msh, const Function& level_set_function, size_t levels) {

    if (levels == 0)
        return;

    auto storage = msh.backend_storage();
    size_t interface_points = iexp_pow(2, levels);

    for (auto& cl : cells(msh)) {
    
        if ( !is_cut(msh, cl) )
            continue;

        auto cl_id = offset(msh, cl);
        auto& cl_interface = storage -> cut_cell_data[cl_id].interface;
        auto& cl_p0 = storage -> cut_cell_data[cl_id].p0;
        auto& cl_p1 = storage -> cut_cell_data[cl_id].p1;
        
        cl_interface.resize(interface_points+1);
        cl_interface.at(0)                = cl_p0;
        cl_interface.at(interface_points) = cl_p1;
        refine_interface(msh, cl, level_set_function, 0, interface_points);

    }
}


template<mesh_2D MeshType, typename Function>
void
refine_interface(MeshType& msh, const typename MeshType::cell_type& cl, const Function& level_set_function, size_t min, size_t max) {
    
    if ((max-min) < 2)
        return;

    typedef typename MeshType::point_type point_type;
    auto storage = msh.backend_storage();

    size_t mid = (max+min)/2;
    auto cl_id = offset(msh, cl);
    auto& cl_interface = storage -> cut_cell_data[cl_id].interface;
    auto p0 = cl_interface.at(min);
    auto p1 = cl_interface.at(max);
    auto pm = (p0+p1)/2.0;
    auto pt = p1 - p0;
    auto pn = point_type(-pt.y(), pt.x());
    auto ps1 = pm + pn;
    auto ps2 = pm - pn;

    auto lm = level_set_function(pm);
    auto ls1 = level_set_function(ps1);
    auto ls2 = level_set_function(ps2);

    point_type ip;

    if ( !((lm >= 0 && ls1 >= 0) || (lm < 0 && ls1 < 0)) ) {
        auto threshold = diameter(msh, cl) / 1e20;
        ip = find_zero_crossing(pm, ps1, level_set_function, threshold);
    }
    else if ( !((lm >= 0 && ls2 >= 0) || (lm < 0 && ls2 < 0)) ) {
        auto threshold = diameter(msh, cl) / 1e20;
        ip = find_zero_crossing(pm, ps2, level_set_function, threshold);
    }
    else
        throw std::logic_error("interface not found in search range");

    cl_interface.at(mid) = ip;

    refine_interface(msh, cl, level_set_function, min, mid);
    refine_interface(msh, cl, level_set_function, mid, max);
}

////// loc_agglo
// container for information about local agglomerations
// main_cell is defined when they are at least 3 cells
template<mesh_2D MeshType>
class loc_agglo {

public:

    MeshType::cell_type main_cell;
    std::vector<MeshType::cell_type> cells;

    MeshType::cell_type new_cell;

    loc_agglo(MeshType::cell_type cl1, MeshType::cell_type cl2, MeshType::cell_type n_cell) {
        
        // check that the two cells are neighbors
        auto pts1 = cl1.point_ids();
        auto pts2 = cl2.point_ids();
        std::vector<size_t> com_nodes;
        for(size_t i = 0; i < pts1.size(); i++) {
            for(size_t j = 0; j < pts2.size(); j++) {
                if(pts1[i] == pts2[j]) com_nodes.push_back(pts1[i]);
            }
        }

        assert( com_nodes.size() > 0 );

        // init
        cells.push_back( cl1 );
        cells.push_back( cl2 );

        new_cell = n_cell;
    }

    bool is_agglo_possible(MeshType msh, size_t offset_) {

        auto storage = msh.backend_storage();

        // when there are 2 cells, the main cell is not defined -> check all cells
        if( cells.size() == 2 ) {

            auto cl0_id = offset(msh, cells[0]);
            auto cl1_id = offset(msh, cells[1]);
            auto cl0_f_neigh = storage -> cut_cell_data[cl0_id].f_neighbors;
            auto cl0_d_neigh = storage -> cut_cell_data[cl0_id].d_neighbors;
            auto cl1_f_neigh = storage -> cut_cell_data[cl1_id].f_neighbors;
            auto cl1_d_neigh = storage -> cut_cell_data[cl1_id].d_neighbors;
            
            return (cl0_f_neigh.find(offset_) != cl0_f_neigh.end())
                || (cl1_f_neigh.find(offset_) != cl1_f_neigh.end())
                || (cl0_d_neigh.find(offset_) != cl0_d_neigh.end())
                || (cl1_d_neigh.find(offset_) != cl1_d_neigh.end());
        }

        // else the main cell is defined -> can agglomerate only if it is a neighbor of that cell
        auto main_cell_id = offset(msh, main_cell);
        auto main_cell_f_neigh = storage -> cut_cell_data[main_cell_id].f_neighbors;
        auto main_cell_d_neigh = storage -> cut_cell_data[main_cell_id].d_neighbors;
        return (main_cell_f_neigh.find(offset_) != main_cell_f_neigh.end())
            || (main_cell_d_neigh.find(offset_) != main_cell_d_neigh.end());
    }

    bool is_in(MeshType::cell_type cl) {
        bool ret = false;
        for(size_t i = 0; i < cells.size(); i++) {
            if( cells[i] == cl) {
                ret = true;
                break;
            }
        }
        return ret;
    }

    void add_cell(MeshType msh, MeshType::cell_type cl_added, size_t offset_, MeshType::cell_type new_new_cell) {

        assert(is_agglo_possible(offset_));
        
        auto storage = msh.backend_storage();
        auto cl0_id = offset(msh, cells[0]);
        auto cl1_id = offset(msh, cells[1]);
        auto cl0_f_neigh = storage -> cut_cell_data[cl0_id].f_neighbors;
        auto cl0_d_neigh = storage -> cut_cell_data[cl0_id].d_neighbors;
        auto cl1_f_neigh = storage -> cut_cell_data[cl1_id].f_neighbors;
        auto cl1_d_neigh = storage -> cut_cell_data[cl1_id].d_neighbors;

        // check that the cell is not already in the agglomeration
        if( is_in(cl_added) )
            throw std::logic_error("Cell already agglomerated !!");
            
        // if there are only two cells, we also have to define the main cell
        if (cells.size() == 2) {
            if (cl0_f_neigh.find(offset_) != cl0_f_neigh.end() || cl0_d_neigh.find(offset_) != cl0_d_neigh.end() )
                main_cell = cells[0];
            else
                main_cell = cells[1];
        }

        cells.push_back(cl_added);
        new_cell = new_new_cell;
    }
};

////// find_good_neighbor
template<typename MeshType>
typename MeshType::cell_type
find_good_neighbor(const MeshType& msh, const typename MeshType::cell_type cl, const location where) {

    auto storage = msh.backend_storage();
    typename MeshType::cell_type best_neigh = cl;

    location other_where;
    if (where == location::IN_NEGATIVE_SIDE)
        other_where = location::IN_POSITIVE_SIDE;
    else
        other_where = location::IN_NEGATIVE_SIDE;

    // look for a neighboring cell to merge with
    typename MeshType::coordinate_type area = 1000000;

    auto f_neigh = storage -> cut_cell_data[offset(msh,cl)].f_neighbors;
    
    for (std::set<size_t>::iterator it = f_neigh.begin(); it != f_neigh.end(); ++it) {
        auto cl_n = msh[*it];

        // if cl_n is on the wrong size -> do not consider it
        if (locate(msh, cl_n) != where && locate(msh, cl_n) != location::ON_INTERFACE)
            continue;

        // if cl_n is a small cut of the same size -> do not consider it
        auto& cl_n_cut_type = storage -> cut_cell_data[offset(msh, cl_n)].cut;
        if (where == location::IN_NEGATIVE_SIDE && cl_n_cut_type == cut_type::T_KO_NEG)
            continue;
        if (where == location::IN_POSITIVE_SIDE && cl_n_cut_type == cut_type::T_KO_POS)
            continue;

        // search for the "best" neighbor -> the one with the smallest volume in other_where
        if( area > measure(msh, cl_n, other_where) )
        {
            area = measure(msh, cl_n, other_where);
            best_neigh = cl_n;
        }
    }

    if(best_neigh == cl) // no possible face agglomerations
    {
        // look for a diagonal agglomeration
        auto d_neigh = storage -> cut_cell_data[offset(msh,cl)].d_neighbors;

        for (std::set<size_t>::iterator it = d_neigh.begin(); it != d_neigh.end(); ++it) {
            auto cl_n = msh[*it];

            // if cl_n is on the wrong size or cut -> do not consider it
            if (locate(msh, cl_n) != where)
                continue;

            // search for the "best" neighbor -> the one with the smallest volume
            if( area > measure(msh, cl_n, other_where) ) {
                area = measure(msh, cl_n, other_where);
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
        size_t TN = table_neg.at(cp);
        if (TN != -1) {
            auto bar_cl = barycenter(msh,cl);
            auto bar_neigh = barycenter(msh,msh[TN]);

            auto vect = bar_neigh - bar_cl;

            output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                   << vect[0] << "   " << vect[1] << std::endl;
        }


        size_t TP = table_pos.at(cp);
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
        size_t TN = table_neg.at(cp);
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
        size_t TP = table_pos.at(cp);
        if( TP != -1) {
          auto h = diameter(msh, cl);
          auto bar_cl = barycenter(msh,cl);
          auto bar_neigh = barycenter(msh, msh[TP]);
          auto vect = bar_neigh - bar_cl;
        //   if (vect[0] <= 1e-5) {
        //     bar_cl[0] = bar_cl[0];
        //   }
        //   else if (vect[1] <= 1e-5) {
        //     bar_cl[1] = bar_cl[1]; 
        //   }
          output << bar_cl[0] << "   " << bar_cl[1] << "   0.   "
                 << vect[0] << "   " << vect[1] << std::endl;
        }
        cp++;
    }
    output.close();
}

template<typename Mesh, typename Function>
void
make_polynomial_extension(Mesh& msh, const Function& level_set_function) {

    auto storage = msh.backend_storage();
    
    // Initiate lists to store the pairing infos
    std::vector<int> table_neg, table_pos;
    size_t nb_cells = msh.cells_size();
    table_neg.resize(nb_cells);
    table_pos.resize(nb_cells);
    for(size_t i=0; i < nb_cells; i++) {
        table_neg.at(i) = -1;
        table_pos.at(i) = -1;
    }

    ///////////////////////   LOOK FOR NEIGHBORS  ////////////////
    size_t nb_step1 = 0;
    size_t nb_step2 = 0;

    // start the process for domain 1, and then domain 2
    for(size_t domain=1; domain < 3; domain++) {

        // loop on the cells
        for (auto& cl : cells(msh)) {

            auto cl_id = offset(msh, cl);
            auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
            location where;
            if(domain == 1)
                where = location::IN_NEGATIVE_SIDE;
            else if(domain == 2)
                where = location::IN_POSITIVE_SIDE;
            else 
                throw std::logic_error("pb with domain");
                
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
            bool already_agglo = false;
            size_t offset_cl = offset(msh, cl);
            for (size_t i = 0; i < table_pos.size(); i++) {
                if (table_pos.at(i) == offset_cl || table_neg.at(i) == offset_cl ) {
                    already_agglo = true;
                    break;
                }
            }
            if( already_agglo )
                continue;


            typename Mesh::cell_type best_neigh = find_good_neighbor(msh, cl, where);

            auto f_neigh = storage -> cut_cell_data[offset(msh,cl)].f_neighbors;

            // prepare agglomeration of cells cl and best_neigh
            size_t offset1 = offset(msh,cl);
            size_t offset2 = offset(msh,best_neigh);

            if(where == location::IN_NEGATIVE_SIDE) {
                table_neg.at(offset1) = offset2;
                nb_step1++;
            }
            else {
                table_pos.at(offset1) = offset2;
                nb_step2++;
            }
        }

        if(domain == 1)
            output_agglo_lists(msh, table_neg, table_pos, "agglo_one.okc");
        if(domain == 2)
            output_agglo_lists(msh, table_neg, table_pos, "agglo_two.okc");
    }
    //////////////   CHANGE THE AGGLO FOR THE CELLS OF DOMAIN 1 THAT ARE TARGETTED ///////
    size_t nb_step3 = 0;
    for (auto& cl : cells(msh)) {

        auto cl_id = offset(msh, cl);
        auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
        if(cl_cut_type != cut_type::T_KO_NEG)
            continue;

        size_t offset1 = offset(msh,cl);

        // are there cells that try to agglomerate with cl ?
        bool agglo = false;
        size_t cl2_offset;
        for(size_t i = 0; i < table_pos.size(); i++) {
            if(table_pos.at(i) == offset1) {
                agglo = true;
                cl2_offset = i;
                break;
            }
        }
        if(!agglo)
            continue;

        // at this point cl2_offset tries to agglomerate with cl
        size_t cl1_agglo = table_neg.at(offset1);
        
        // -> check that no one tries to agglomerate with cl1_agglo
        agglo = false;
        for(size_t i = 0; i < table_neg.size(); i++) {
            if( i == offset1)
                continue;

            if(table_neg.at(i) == cl1_agglo) {
                agglo = true;
                break;
            }
        }

        auto& cl1_cut_type = storage -> cut_cell_data[cl1_agglo].cut;
        if (!agglo && cl1_cut_type == cut_type::T_KO_POS) {
            continue;
        }

        // at this point we risk chain agglomerations
        // -> remove the target of cl
        nb_step3++;
        table_neg.at(offset1) = -1;
    }
    output_agglo_lists(msh, table_neg, table_pos, "agglo_three.okc");
    

    ///////////////////////////////////////////////////////////////////////////// STEP 4 
    // ALL BAD CUT CELLS MUST POINT TOWARDS A CELL

    // loop on the cells
    for (auto& cl : cells(msh)) {
        auto cl_id = offset(msh, cl);
        auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
        // FIND THE CELL OFFSET
        auto offset_cl = offset(msh,cl);
        auto TN = table_neg.at(offset_cl);
        auto TP = table_pos.at(offset_cl);
        if (TN == -1) {    
            if (cl_cut_type == cut_type::T_KO_NEG) {  
                for(size_t i = 0; i < table_pos.size(); i++) {
                    if(table_pos.at(i) == offset_cl) {
                        table_neg.at(offset_cl) = i;
                    break;
                    }
                }
            }
        }
        if (TP == -1) {
            if (cl_cut_type == cut_type::T_KO_POS) { 
                for(size_t i = 0; i < table_neg.size(); i++) {
                    if(table_neg.at(i) == offset_cl) {
                        table_pos.at(offset_cl) = i;
                    break;
                    }
                }
            }
        }
    }

    // FILLING THE STRUCTURES paired_cells / dependent_cells_neg / dependent_cells_pos / paire(T,i)
    for (auto& cl : cells(msh)) {
        auto cl_id = offset(msh, cl);
        auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
        auto& cl_paired_cell = storage -> cut_cell_data[cl_id].stabilizing_cell;
        auto offset_cl = offset(msh,cl);
        if (cl_cut_type == cut_type::T_KO_NEG) {
            if (table_neg.at(offset_cl) != -1) 
                cl_paired_cell = table_neg.at(offset_cl);
            if (table_pos.at(offset_cl) != -1) 
                cl_paired_cell = table_pos.at(offset_cl);
            const auto& good_cl = msh[cl_paired_cell.value()];
            auto& good_dependent_cells_neg = storage -> cut_cell_data[offset(msh,good_cl)].dependent_ill_cut_neg_cells;
            good_dependent_cells_neg.push_back(offset_cl);
        }
        else if (cl_cut_type == cut_type::T_KO_POS) {
            if (table_neg.at(offset_cl) != -1) 
                cl_paired_cell = table_neg.at(offset_cl);
            if (table_pos.at(offset_cl) != -1) 
                cl_paired_cell = table_pos.at(offset_cl);
            const auto& good_cl = msh[cl_paired_cell.value()];
            auto& good_dependent_cells_pos = storage -> cut_cell_data[offset(msh,good_cl)].dependent_ill_cut_pos_cells;
            good_dependent_cells_pos.push_back(offset_cl);
        }
    }

    // Display of the arrows
    std::vector<int> table;
    table.resize(nb_cells);
    for(size_t i=0; i < nb_cells; i++) 
        table.at(i) = -1;
    

    output_agglo_lists_step4(msh, table_neg, table, "agglo_four.okc");
    output_agglo_lists_step4(msh, table, table_pos, "agglo_five.okc");

}

template<typename Mesh>
void print_polynomial_extension(Mesh& msh) {

    auto storage = msh.backend_storage();

    for (auto& cl : cells(msh)) {
        auto cl_id = offset(msh, cl);
        std::cout << std::endl << "CELL " << cl_id << " IS ";
        auto& cl_loc = storage -> cut_cell_data[cl_id].loc;
        auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
        auto& cl_paired_cell = storage -> cut_cell_data[cl_id].stabilizing_cell;
        auto& dependent_cells_neg = storage -> cut_cell_data[cl_id].dependent_ill_cut_neg_cells;
        auto& dependent_cells_pos = storage -> cut_cell_data[cl_id].dependent_ill_cut_pos_cells;
        if (cl_loc != location::ON_INTERFACE) {
            std::cout << "UNCUT: " << std::endl;
        }
        else {
            std::cout << "CUT: ";
            if (cl_cut_type == cut_type::T_OK) {
                std::cout << "T_OK" << std::endl;
            }
            else if (cl_cut_type == cut_type::T_KO_NEG) {
                std::cout << "T_KO_NEG" << std::endl;
            }
            else if (cl_cut_type == cut_type::T_KO_POS) {
                std::cout << "T_KO_POS" << std::endl;
            }
        }
        if (cl_paired_cell.has_value()) {
            std::cout << "STABILIZING CELL: " << cl_paired_cell.value() << std::endl;
        }
        else {
            std::cout << "NO STABILIZING CELL" << std::endl;
        }
        for (auto dp_cl_neg : dependent_cells_neg) {
            std::cout << "DEPENDENT NEGATIVE CELL: " << dp_cl_neg << std::endl;
        }
        for (auto dp_cl_pos : dependent_cells_pos) {
            std::cout << "DEPENDENT POSITIVE CELL: " << dp_cl_pos << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl << std::endl;
    }
}

template<typename Mesh> std::pair<std::vector< std::tuple<double,location,std::vector<double>>>, std::vector< std::tuple<double,location,std::vector<double>>>> 
make_pair_KO_pair_OK(Mesh& msh) {

    std::vector< std::tuple<double,location,std::vector<double>>> PairOK;
    std::vector< std::tuple<double,location,std::vector<double>>> PairKO;

    auto storage = msh.backend_storage();

    for (auto &cl : msh.cells) {
        auto offset_cl = offset(msh,cl);
        auto& cl_loc = storage -> cut_cell_data[offset_cl].loc;
        auto& cl_cut_type = storage -> cut_cell_data[offset_cl].cut;
        auto& cl_paired_cell = storage -> cut_cell_data[offset_cl].stabilizing_cell;
        auto& dependent_cells_neg = storage -> cut_cell_data[offset_cl].dependent_ill_cut_neg_cells;
        auto& dependent_cells_pos = storage -> cut_cell_data[offset_cl].dependent_ill_cut_pos_cells;
        if (cl_loc != location::ON_INTERFACE) { 
            std::vector<double> dp_cells;
            location loc;
            if (cl_loc == location::IN_NEGATIVE_SIDE) {
                loc = location::IN_NEGATIVE_SIDE;
                for (auto& dp_cl: dependent_cells_neg) 
                    dp_cells.push_back(dp_cl);
            }
            else {
                loc = location::IN_POSITIVE_SIDE;
                for (auto& dp_cl: dependent_cells_pos) 
                    dp_cells.push_back(dp_cl);
            }
            auto tuple = std::make_tuple(offset_cl, loc, dp_cells);
            PairOK.push_back(tuple);
        }
        else if (cl_cut_type == cut_type::T_OK) { 
            std::vector<double> dp_cells_neg;
            std::vector<double> dp_cells_pos;
            for (auto& dp_cl: dependent_cells_neg) 
                dp_cells_neg.push_back(dp_cl);
            for (auto& dp_cl: dependent_cells_pos) 
                dp_cells_pos.push_back(dp_cl);
            auto tuple_neg = std::make_tuple(offset_cl, location::IN_NEGATIVE_SIDE, dp_cells_neg);
            auto tuple_pos = std::make_tuple(offset_cl, location::IN_POSITIVE_SIDE, dp_cells_pos);
            PairOK.push_back(tuple_neg);
            PairOK.push_back(tuple_pos);
        }
        else if (cl_cut_type == cut_type::T_KO_NEG) { 
            std::vector<double> dp_cells_neg;
            std::vector<double> dp_cells_pos;
            for (auto& dp_cl: dependent_cells_pos) 
                dp_cells_pos.push_back(dp_cl);
            auto tuple_neg = std::make_tuple(offset_cl, location::IN_NEGATIVE_SIDE, dp_cells_neg);
            auto tuple_pos = std::make_tuple(offset_cl, location::IN_POSITIVE_SIDE, dp_cells_pos);
            PairKO.push_back(tuple_neg);
            PairOK.push_back(tuple_pos);
        }
        else if (cl_cut_type == cut_type::T_KO_POS) {
            std::vector<double> dp_cells_neg;
            std::vector<double> dp_cells_pos;
            for (auto& dp_cl: dependent_cells_neg) 
                dp_cells_neg.push_back(dp_cl);
            auto tuple_neg = std::make_tuple(offset_cl, location::IN_NEGATIVE_SIDE, dp_cells_neg);
            auto tuple_pos = std::make_tuple(offset_cl, location::IN_POSITIVE_SIDE, dp_cells_pos);
            PairOK.push_back(tuple_neg);
            PairKO.push_back(tuple_pos);
        }
    }  

    auto Pair_OK_KO = std::make_pair(PairOK,PairKO);
    return Pair_OK_KO;

}

void modify_dependent_cells(std::vector<std::tuple<double, location, std::vector<double>>>& PairOK, std::vector<std::tuple<double, location, std::vector<double>>>& PairKO, double cell_index_1, double cell_index_2, location loc) {

    // ADD CELL2 TO THE DEPENDENT CELLS OF CELL1
    for (auto& cell_tuple : PairOK) {
        if (std::get<0>(cell_tuple) == cell_index_1 && std::get<1>(cell_tuple) == loc) {
            // Ajouter cell_index_2 aux cellules dépendantes de cell_index_1
            std::get<2>(cell_tuple).push_back(cell_index_2);
            break;
        }
    }

    // REMOVE CELL2 FROM POK AND ADD IT TO PKO
    auto it = PairOK.begin();
    while (it != PairOK.end()) {
        if (std::get<0>(*it) == cell_index_2 && std::get<1>(*it) == loc) {
            PairKO.push_back(*it);
            it = PairOK.erase(it);
        } 
        else 
            ++it;
    }
}

//////////////  MERGE_CELLS_FACE
/// merge cl1 and cl2 through the common face fc
// output : new cell
template<typename cell_type, typename face_type>
cell_type
merge_cells_face(cell_type cl1, cell_type cl2, face_type com_f)
{
    cell_type ret;

    size_t f_pt1 = com_f.ptids[0];
    size_t f_pt2 = com_f.ptids[1];

    // list of points
    std::vector<size_t> pts1, pts2;
    // in order to be consistent with the cell representation, start with the smallest index
    if(cl1.ptids[0] < cl2.ptids[0])
    {
        pts1 = cl1.ptids;
        pts2 = cl2.ptids;
    }
    else
    {
        pts1 = cl2.ptids;
        pts2 = cl1.ptids;
    }

    // write points of pts1 until we reach the common face
    size_t ref_pt = pts1[0];
    size_t cp = 0;
    ret.ptids.push_back(ref_pt);

    bool on_face = false;
    if(ref_pt == f_pt1 || ref_pt == f_pt2) on_face = true;

    while(!on_face)
    {
        cp++;
        ref_pt = pts1[cp];
        ret.ptids.push_back(ref_pt);
        if(ref_pt == f_pt1 || ref_pt == f_pt2) on_face = true;
    }

    // look for the corresponding point in pts2
    for(size_t i = 0; i < pts2.size(); i++)
    {
        if(ref_pt == pts2[i])
        {
            cp = i;
            break;
        }
    }
    // write points of pts2 until we reach once more the common face
    on_face = false;
    while( !on_face )
    {
        cp = (cp + 1) % pts2.size();
        ref_pt = pts2[cp];
        ret.ptids.push_back(ref_pt);
        if(ref_pt == f_pt1 || ref_pt == f_pt2) on_face = true;
    }
    // look for the corresponding point in pts1
    for (size_t i=0; i < pts1.size(); i++) {
        if (ref_pt == pts1[i]) {
            cp = i;
            break;
        }
    }

    // finish to write the points of pts1
    cp++;
    while(cp < pts1.size()) {
        ret.ptids.push_back(pts1[cp]);
        cp++;
    }

    return ret;
}

//////////////  MERGE_CELLS_TWO_FACES
/// merge cl1 and cl2 through the list com_f of common faces
// output : new cell
template<typename cell_type, typename face_type>
cell_type
merge_cells_two_faces(cell_type cl1, cell_type cl2, typename std::vector<face_type> com_f)
{
    cell_type ret;

    assert( com_f.size() == 2 );

    size_t f1_pt1 = com_f[0].ptids[0];
    size_t f1_pt2 = com_f[0].ptids[1];
    size_t f2_pt1 = com_f[1].ptids[0];
    size_t f2_pt2 = com_f[1].ptids[1];

    std::set<size_t> com_pts;
    if(f1_pt1 == f2_pt1 || f1_pt2 == f2_pt1)
    {
        com_pts.insert(f1_pt1);
        com_pts.insert(f1_pt2);
        com_pts.insert(f2_pt2);
    }
    else if(f1_pt1 == f2_pt2 || f1_pt2 == f2_pt2)
    {
        com_pts.insert(f1_pt1);
        com_pts.insert(f1_pt2);
        com_pts.insert(f2_pt1);
    }
    else
        throw std::logic_error("com_f : no common points");

    // list of points
    std::vector<size_t> pts1, pts2;
    // in order to be consistent with the cell representation, start with the smallest index
    if(cl1.ptids[0] < cl2.ptids[0])
    {
        pts1 = cl1.ptids;
        pts2 = cl2.ptids;
    }
    else
    {
        pts1 = cl2.ptids;
        pts2 = cl1.ptids;
    }

    // write points of pts1 until we reach the common face
    size_t ref_pt = pts1[0];
    size_t next_pt= pts1[1];
    size_t cp = 0;
    ret.ptids.push_back(ref_pt);

    bool on_face = false;
    if( com_pts.find(ref_pt) != com_pts.end() && com_pts.find(next_pt) != com_pts.end() )
        on_face = true;

    while(!on_face)
    {
        cp++;
        ref_pt = next_pt;
        next_pt= pts1[(cp+1)%pts1.size()];
        ret.ptids.push_back(ref_pt);
        if( com_pts.find(ref_pt) != com_pts.end() && com_pts.find(next_pt) != com_pts.end() )
            on_face = true;
    }

    // look for the corresponding point in pts2
    for(size_t i = 0; i < pts2.size(); i++)
    {
        if(ref_pt == pts2[i])
        {
            cp = i;
            break;
        }
    }

    // write points of pts2 until we reach once more the common face
    on_face = false;
    while( !on_face )
    {
        cp = (cp + 1) % pts2.size();
        ref_pt = pts2[cp];
        next_pt= pts2[(cp+1) % pts2.size()];
        ret.ptids.push_back(ref_pt);
        if( com_pts.find(ref_pt) != com_pts.end() && com_pts.find(next_pt) != com_pts.end() )
            on_face = true;
    }
    // look for the corresponding point in pts1
    for(size_t i=0; i < pts1.size(); i++)
    {
        if(ref_pt == pts1[i])
        {
            cp = i;
            break;
        }
    }

    // finish to write the points of pts1
    cp++;
    while(cp < pts1.size())
    {
        ret.ptids.push_back(pts1[cp]);
        cp++;
    }

    return ret;
}

//////////////  MERGE_CELLS
/// merge cl1 and cl2
// output : the agglomerated cell + list of faces to withdraw
/////  For the moment we can merge only cells that have at least a common face
/////  This procedure currently cannot be iterated
template<typename Mesh>
std::pair<typename Mesh::cell_type, std::vector<typename Mesh::face_type> >
merge_cells(Mesh& msh, const typename Mesh::cell_type cl1,
            const typename Mesh::cell_type cl2)
{

    auto storage = msh.backend_storage();
    auto cl1_id = offset(msh, cl1);
    auto& cl1_loc = storage -> cut_cell_data[cl1_id].loc;
    auto cl2_id = offset(msh, cl2);
    auto& cl2_loc = storage -> cut_cell_data[cl1_id].loc;

    //////////////////  TESTS ON INPUTS  //////////////////
    // verify that the two cells are different
    if(cl1 == cl2)
        throw std::invalid_argument("Cannot merge a cell with itself.");

    // identify the common faces
    const auto fcs1 = faces(msh, cl1);
    const auto fcs2 = faces(msh, cl2);

    std::vector<typename Mesh::face_type> com_faces;
    for(size_t i = 0; i < fcs1.size(); i++)
    {
        const auto fc1 = fcs1[i];
        for(size_t j = 0; j < fcs2.size(); j++)
        {
            const auto fc2 = fcs2[j];
            if(fc1 == fc2) com_faces.push_back(fc1);
        }
    }

    // identify the common nodes
    std::set<size_t> com_nodes;
    size_t com_n;
    auto cl1_ptids = cl1.point_ids();
    auto cl2_ptids = cl2.point_ids();
    for(size_t i = 0; i < cl1_ptids.size(); i++) {
        for(size_t j = 0; j < cl2_ptids.size(); j++) {
            if(cl1_ptids[i] == cl2_ptids[j]) {
                com_nodes.insert(cl1_ptids[i]);
                com_n = cl1_ptids[i];
            }
        }
    }

    // choose the agglomeration technique
    typename Mesh::cell_type cl;
    if(com_faces.size() == 0) {
        std::cout << "com nodes nb = " << com_nodes.size() << std::endl;
        if( com_nodes.size() == 1 ) {
            cl = merge_cells_diag(msh, cl1, cl2, com_n);
        }
        else
            throw std::invalid_argument("The cells have no common faces.");
    }
    else if( com_faces.size() == 1 ) {
        assert(com_nodes.size() == 2); // only the nodes of the common face
        cl = merge_cells_face(cl1, cl2, com_faces[0]);
    }
    else if( com_faces.size() == 2 ) {
        std::cout << "com nodes nb = " << com_nodes.size() << std::endl;
        assert(com_nodes.size() == 3); // only the case with two adjascent common faces
        cl = merge_cells_two_faces(cl1, cl2, com_faces);
    }
    if(com_faces.size() > 2)
        throw std::invalid_argument("The cells have too many common faces.");

    //////////////    COMPLETE THE MERGED CELL   //////////////

    ///////////// build the cell_cuthho_info (using the sub_cells info)
    // location
    auto cl_id = offset(msh, cl);
    auto& cl_loc = storage -> cut_cell_data[cl1_id].loc;
    if(cl1_loc == location::ON_INTERFACE || cl2_loc == location::ON_INTERFACE)
        cl_loc = location::ON_INTERFACE;
    else if(cl1_loc == location::IN_NEGATIVE_SIDE)
        cl_loc = location::IN_NEGATIVE_SIDE;
    else if(cl1_loc == location::IN_POSITIVE_SIDE)
        cl_loc = location::IN_POSITIVE_SIDE;
    else
        throw std::logic_error("we shouldn't arrive here (cuthho_info) !!!");

    // agglo_set ---> NOT DONE : needed to iterate the agglomeration procedure
    auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
    cl_cut_type = cut_type::T_OK;

    // p0, p1 and interface
    auto& cl_p0 = storage -> cut_cell_data[cl_id].p0;
    auto& cl_p1 = storage -> cut_cell_data[cl_id].p1;
    auto& cl_interface = storage -> cut_cell_data[cl_id].interface;
    auto& cl_distorted = storage -> cut_cell_data[cl_id].distorted;
    auto& cl_highlight = storage -> cut_cell_data[cl_id].highlight;
    auto& cl_integration_n = storage -> cut_cell_data[cl_id].integration_n;
    auto& cl_integration_p = storage -> cut_cell_data[cl_id].integration_p;
    auto& cl1_p0 = storage -> cut_cell_data[cl1_id].p0;
    auto& cl1_p1 = storage -> cut_cell_data[cl1_id].p1;
    auto& cl1_interface = storage -> cut_cell_data[cl1_id].interface;
    auto& cl1_distorted = storage -> cut_cell_data[cl1_id].distorted;
    auto& cl1_highlight = storage -> cut_cell_data[cl1_id].highlight;
    auto& cl2_p0 = storage -> cut_cell_data[cl2_id].p0;
    auto& cl2_p1 = storage -> cut_cell_data[cl2_id].p1;
    auto& cl2_interface = storage -> cut_cell_data[cl2_id].interface;
    auto& cl2_distorted = storage -> cut_cell_data[cl2_id].distorted;
    auto& cl2_highlight = storage -> cut_cell_data[cl2_id].highlight;
    bool cut1 = (cl1_loc == location::ON_INTERFACE);
    bool cut2 = (cl2_loc == location::ON_INTERFACE);
    if(cut1 && !cut2) {
        cl_interface = cl1_interface;
        cl_p0 = cl1_p0;
        cl_p1 = cl1_p1;
    }
    else if(!cut1 && cut2) {
        cl_interface = cl2_interface;
        cl_p0 = cl2_p0;
        cl_p1 = cl2_p1;
    }
    else if(cut1 && cut2) { 
        // MODIFY THE LOOP (i=1 instead of i=0) TO DON'T HAVE DUPLICATED INTERFACE POINTS
        if(cl1_p0[0] == cl2_p1[0] && cl1_p0[1] == cl2_p1[1]) {
            cl_interface = cl2_interface;
            for(size_t i = 1; i < cl1_interface.size(); i++ ) {
                cl_interface.push_back(cl1_interface[i]);
            }
            cl_p0 = cl2_p0;
            cl_p1 = cl1_p1;
        }
        else if(cl2_p0[0] == cl1_p1[0] && cl2_p0[1] == cl1_p1[1]) {
            cl_interface = cl1_interface;
            for(size_t i = 1; i < cl2_interface.size(); i++ ) {
                cl_interface.push_back(cl2_interface[i]);
            }
            cl_p0 = cl1_p0;
            cl_p1 = cl2_p1;
        }
        else
            throw std::logic_error("we shouldn't arrive here (interface) !!!");
    }
    // distorted -> has to be updated for more general merges (if a node is withdrawn)
    if (cl1_distorted || cl2_distorted )
        cl_distorted = true;

    // for tests
    cl_highlight = true;

    // integration -> save composite quadrature
    size_t degree_max = 8; //////// VERY IMPORTANT !!!!!!! -> max deg for quadratures = 8
    auto integration1_n = integrate(msh, cl1, degree_max, location::IN_NEGATIVE_SIDE);
    auto integration1_p = integrate(msh, cl1, degree_max, location::IN_POSITIVE_SIDE);

    auto integration2_n = integrate(msh, cl2, degree_max, location::IN_NEGATIVE_SIDE);
    auto integration2_p = integrate(msh, cl2, degree_max, location::IN_POSITIVE_SIDE);

    cl_integration_n = integration1_n;
    cl_integration_p = integration1_p;
    for(size_t i = 0; i < integration2_n.size(); i++)
        cl_integration_n.push_back( integration2_n.at(i) );

    for(size_t i = 0; i < integration2_p.size(); i++)
        cl_integration_p.push_back( integration2_p.at(i) );

    // neighbors ---> NOT DONE : needed to iterate the agglomeration procedure

    return std::make_pair(cl, com_faces);
}

//////////////  MERGE_CELLS_DIAG
/// merge cl1 and cl2 through the common node com_n
// output : new cell
template<typename Mesh>
typename Mesh::cell_type
merge_cells_diag(Mesh& msh, const typename Mesh::cell_type cl1, const typename Mesh::cell_type cl2, size_t com_n)
{

    auto storage = msh.backend_storage();
    auto cl1_id = offset(msh, cl1);
    auto& cl1_loc = storage -> cut_cell_data[cl1_id].loc;
    auto cl2_id = offset(msh, cl2);
    auto& cl2_loc = storage -> cut_cell_data[cl1_id].loc;

    // abort the process if both cells are on the interface
    // (this case is not yet supported for the list of points on the interface)
    if ( cl1_loc == location::ON_INTERFACE && cl2_loc == location::ON_INTERFACE )
        throw std::invalid_argument("Cannot merge diagonally two cells on the interface.");

    typename Mesh::cell_type ret;

    // list of points
    auto pts1 = cl1.point_ids(); 
    auto pts2 = cl2.point_ids();
    // in order to be consistent with the cell representation, start with the smallest index
    if(pts1[0] < pts2[0]) {
        pts1 = cl1.point_ids();
        pts2 = cl2.point_ids();
    }
    else {
        pts1 = cl2.point_ids();
        pts2 = cl1.point_ids();
    }

    // write points of pts1 until we reach the common face
    size_t ref_pt = pts1[0];
    size_t cp = 0;
    ret.ptids.push_back(ref_pt);

    bool on_interface = false;
    if(ref_pt == com_n) on_interface = true;

    while(!on_interface)
    {
        cp++;
        ref_pt = pts1[cp];
        ret.ptids.push_back(ref_pt);
        if(ref_pt == com_n) on_interface = true;
    }

    // look for the corresponding point in pts2
    for(size_t i = 0; i < pts2.size(); i++)
    {
        if(ref_pt == pts2[i])
        {
            cp = i;
            break;
        }
    }
    // write points of pts2 until we reach once more the common node
    on_interface = false;
    while( !on_interface )
    {
        cp = (cp + 1) % pts2.size();
        ref_pt = pts2[cp];
        ret.ptids.push_back(ref_pt);
        if(ref_pt == com_n) on_interface = true;
    }

    // look for the corresponding point in pts1
    for(size_t i=0; i < pts1.size(); i++)
    {
        if(ref_pt == pts1[i])
        {
            cp = i;
            break;
        }
    }

    // finish to write the points of pts1
    cp++;
    while(cp < pts1.size())
    {
        ret.ptids.push_back(pts1[cp]);
        cp++;
    }

    return ret;
}

/////// check_corner
// do another agglomeration if the agglomerated cell is cut four times
// needed for the case of a square interface
// cl_agglo : a cell that is already agglomerated with other cells
// cl : a cell that is not yet agglomerated
// return: a boolean: whether another agglomeration is needed
//         a cell: the cell added to the agglomeration
template<typename Mesh>
std::pair<bool, typename Mesh::cell_type >
check_corner(Mesh& msh, typename Mesh::cell_type cl_agglo, typename Mesh::cell_type cl)
{
    auto storage = msh.backend_storage();
    typename Mesh::cell_type ret_cl;
    
    // check is another agglo is needed
    bool need_other_agglo = true;

    auto& cl_agglo_p0 = storage -> cut_cell_data[offset(msh,cl_agglo)].p0;
    auto& cl_agglo_p1 = storage -> cut_cell_data[offset(msh,cl_agglo)].p1;
    auto& cl_p0 = storage -> cut_cell_data[offset(msh,cl)].p0;
    auto& cl_p1 = storage -> cut_cell_data[offset(msh,cl)].p1;

    if(cl_agglo_p0.x() == cl_p1.x() && cl_agglo_p0.y() == cl_p1.y() )
        need_other_agglo = false;
    if(cl_agglo_p1.x() == cl_p0.x() && cl_agglo_p1.y() == cl_p0.y() )
        need_other_agglo = false;
    if(cl_agglo_p0.x() == cl_p0.x() && cl_agglo_p0.y() == cl_p0.y() )
        need_other_agglo = false;
    if(cl_agglo_p1.x() == cl_p1.x() && cl_agglo_p1.y() == cl_p1.y() )
        need_other_agglo = false;

    if(!need_other_agglo)
        return std::make_pair(false, ret_cl);

    std::cout << "agglomerate one more cell !!" << std::endl;

    // look for a cell that shares two intersection points with cl_agglo and cl
    // this cell is searched among the face neighbors
    auto f_neigh = storage -> cut_cell_data[offset(msh,cl)].f_neighbors;
    typename Mesh::cell_type added_cell;
    bool found_added_cell = false;
    size_t offset_added_cell;
    for (std::set<size_t>::iterator it = f_neigh.begin(); it != f_neigh.end(); ++it) {
        
        auto cl_f = msh[*it];
        auto& cl_f_p0 = storage -> cut_cell_data[offset(msh,cl_f)].p0;
        auto& cl_f_p1 = storage -> cut_cell_data[offset(msh,cl_f)].p1;

        if(locate(msh, cl_f) != location::ON_INTERFACE)
            continue;

        if((cl_f_p0.x() == cl_p0.x() && cl_f_p0.y() == cl_p0.y() )
            || (cl_f_p1.x() == cl_p0.x() && cl_f_p1.y() == cl_p0.y() )
            || (cl_f_p0.x() == cl_p1.x() && cl_f_p0.y() == cl_p1.y() )
            || (cl_f_p1.x() == cl_p1.x() && cl_f_p1.y() == cl_p1.y() ) ) { // this cell is connected with cl
            if((cl_f_p0.x() == cl_agglo_p0.x() && cl_f_p0.y() == cl_agglo_p0.y() )
                || (cl_f_p1.x() == cl_agglo_p0.x() && cl_f_p1.y() == cl_agglo_p0.y())
                || (cl_f_p0.x() == cl_agglo_p1.x() && cl_f_p0.y() == cl_agglo_p1.y())
                || (cl_f_p1.x() == cl_agglo_p1.x() && cl_f_p1.y() == cl_agglo_p1.y())) { // this cell is connected with cl_agglo
                added_cell = cl_f;
                offset_added_cell = *it;
                found_added_cell = true;
                break;
            }
        }
    }
    
    if(!found_added_cell)    
        throw std::logic_error("added_cell not found !!");

    return std::make_pair(true, added_cell);

}

// // MAIN AGGLOMERATION ROUTINE
// // Resulting mesh obtained must have convex cells by merging sub_cells with one face in common
// // For non-convex cells -> modify measure and integrate
// // For merging more general cells -> modify merge_cells
// // For diagonal cells -> modify make_neighbors_info
// template<mesh_2D MeshType, typename Function>
// void
// make_agglomeration(MeshType& msh, const Function& level_set_function) {

//     // initiate lists to store the agglomeration infos
//     auto storage = msh.backend_storage();
//     std::vector<int> agglo_table_neg, agglo_table_pos;
//     size_t nb_cells = msh.cells_size();
//     agglo_table_neg.resize(nb_cells);
//     agglo_table_pos.resize(nb_cells);
    
//     for(size_t i=0; i < nb_cells; i++) {
//         agglo_table_neg.at(i) = -1;
//         agglo_table_pos.at(i) = -1;
//     }

//     ///////////////////////   LOOK FOR NEIGHBORS  ////////////////
//     size_t nb_step1 = 0;
//     size_t nb_step2 = 0;
//     // start the process for domain 1, and then domain 2
//     for(size_t domain=1; domain < 3; domain++) {

//         for (auto& cl : cells(msh)) {
//             auto cl_id = offset(msh, cl);
//             auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
//             location where;
//             if (domain == 1)
//                 where = location::IN_NEGATIVE_SIDE;
//             else if (domain == 2)
//                 where = location::IN_POSITIVE_SIDE;
//             else 
//                 throw std::logic_error("pb with domain");
                
//             if (cl_cut_type == cut_type::T_OK)
//                 continue;
//             else if (cl_cut_type == cut_type::UNDEF)
//                 throw std::logic_error("UNDEF agglo_set");
//             else if (cl_cut_type == cut_type::T_KO_NEG
//                     && where != location::IN_NEGATIVE_SIDE)
//                 continue;
//             else if (cl_cut_type == cut_type::T_KO_POS
//                     && where != location::IN_POSITIVE_SIDE)
//                 continue;


//             // if cl is already agglomerated : no need for further agglomerations
//             bool already_agglo = false;
//             size_t offset_cl = offset(msh, cl);
//             for (size_t i = 0; i < agglo_table_pos.size(); i++)
//             {
//                 if( agglo_table_pos.at(i) == offset_cl || agglo_table_neg.at(i) == offset_cl )
//                 {
//                     already_agglo = true;
//                     break;
//                 }
//             }
//             if( already_agglo )
//                 continue;


//             typename MeshType::cell_type best_neigh = find_good_neighbor(msh, cl, where);

//             auto f_neigh = storage -> cut_cell_data[offset(msh,cl)].f_neighbors;

//             // prepare agglomeration of cells cl and best_neigh
//             size_t offset1 = offset(msh,cl);
//             size_t offset2 = offset(msh,best_neigh);

//             if(where == location::IN_NEGATIVE_SIDE)
//             {
//                 agglo_table_neg.at(offset1) = offset2;
//                 nb_step1++;
//             }
//             else
//             {
//                 agglo_table_pos.at(offset1) = offset2;
//                 nb_step2++;
//             }
//         }

//         if(domain == 1)
//             output_agglo_lists(msh, agglo_table_neg, agglo_table_pos, "agglo_one.okc");
//         if(domain == 2)
//             output_agglo_lists(msh, agglo_table_neg, agglo_table_pos, "agglo_two.okc");
//     }
//     //////////////   CHANGE THE AGGLO FOR THE CELLS OF DOMAIN 1 THAT ARE TARGETTED ///////
//     size_t nb_step3 = 0;
//     for (auto& cl : cells(msh)) {
        
//         auto cl_id = offset(msh, cl);
//         auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
//         if(cl_cut_type != cut_type::T_KO_NEG)
//             continue;

//         size_t offset1 = offset(msh,cl);

//         // are there cells that try to agglomerate with cl ?
//         bool agglo = false;
//         size_t cl2_offset;
//         for(size_t i = 0; i < agglo_table_pos.size(); i++)
//         {
//             if(agglo_table_pos.at(i) == offset1)
//             {
//                 agglo = true;
//                 cl2_offset = i;
//                 break;
//             }
//         }
//         if(!agglo)
//             continue;

//         // at this point cl2_offset tries to agglomerate with cl
//         size_t cl1_agglo = agglo_table_neg.at(offset1);
        
//         // -> check that no one tries to agglomerate with cl1_agglo
//         agglo = false;
//         for(size_t i = 0; i < agglo_table_neg.size(); i++)
//         {
//             if( i == offset1)
//                 continue;

//             if(agglo_table_neg.at(i) == cl1_agglo)
//             {
//                 agglo = true;
//                 break;
//             }
//         }

//         auto& cl1_cut_type = storage -> cut_cell_data[cl1_agglo].cut;
//         if (!agglo && cl1_cut_type == cut_type::T_KO_POS)
//             continue;

//         // at this point we risk chain agglomerations
//         // -> remove the target of cl
//         nb_step3++;
//         agglo_table_neg.at(offset1) = -1;
//     }

//     output_agglo_lists(msh, agglo_table_neg, agglo_table_pos, "agglo_three.okc");
    
//     ///////////////////  BUILD LOCAL AGGLOMERATIONS  //////////////////
//     std::vector<loc_agglo<MeshType>> loc_agglos;
//     std::vector<typename MeshType::face_type> removed_faces;
//     std::vector<typename MeshType::cell_type> removed_cells;
//     for (auto& cl : cells(msh)) {
//         auto offset_cl = offset(msh, cl);

//         bool agglo = false;
//         size_t offset_neigh;
//         if(agglo_table_pos.at(offset_cl) != -1)
//         {
//             offset_neigh = agglo_table_pos.at(offset_cl);
//             agglo = true;
//         }
//         if(agglo_table_neg.at(offset_cl) != -1)
//         {
//             offset_neigh = agglo_table_neg.at(offset_cl);
//             agglo = true;
//         }

//         if(!agglo)
//             continue;

//         typename MeshType::cell_type neigh = msh[offset_neigh];

        
//         // test if one of the two cells is already agglomerated
//         size_t agglo_offset_cl, agglo_offset_neigh;
//         bool already_agglo_cl = false;
//         for (size_t i = 0; i < loc_agglos.size(); i++)
//         {
//             if( loc_agglos.at(i).is_in(cl) )
//             {
//                 already_agglo_cl = true;
//                 agglo_offset_cl = i;
//                 break;
//             }
//         }
//         bool already_agglo_neigh = false;
//         for (size_t i = 0; i < loc_agglos.size(); i++)
//         {
//             if( loc_agglos.at(i).is_in(neigh) )
//             {
//                 already_agglo_neigh = true;
//                 agglo_offset_neigh = i;
//                 break;
//             }
//         }

//         // the two cells can not be both already agglomerated in different agglomeration sets
//         if(already_agglo_cl && already_agglo_neigh)
//         {
//             if(agglo_offset_neigh != agglo_offset_cl)
//             {
//                 throw std::logic_error("Both cells already agglomerated !!");
//             }
//             else
//                 std::cout << "agglo already done" << std::endl;
//             // else : the two cells are already agglomerated together : DO NOTHING
//         }
//         else if(!already_agglo_cl && !already_agglo_neigh)
//         {
//             // create a new local agglomeration
//             auto MC = merge_cells(msh, cl, neigh);

//             loc_agglos.push_back(loc_agglo<MeshType>(cl, neigh, MC.first) );

//             for(size_t i=0; i<MC.second.size(); i++)
//             {
//                 removed_faces.push_back( MC.second[i] );
//             }
//             removed_cells.push_back(cl);
//             removed_cells.push_back(neigh);
//         }
//         else // only one cell is already agglomerated
//         {
//             typename MeshType::cell_type cl1, cl2;
//             size_t offset_cl2, agglo_offset;
//             if(already_agglo_cl)
//             {
//                 removed_cells.push_back(neigh);
//                 cl2 = neigh;
//                 offset_cl2 = offset_neigh;
//                 agglo_offset = agglo_offset_cl;
//             }
//             else if(already_agglo_neigh)
//             {
//                 removed_cells.push_back(cl);
//                 cl2 = cl;
//                 offset_cl2 = offset_cl;
//                 agglo_offset = agglo_offset_neigh;
//             }

//             // get the agglomerated cell
//             cl1 = loc_agglos.at(agglo_offset).new_cell;
            
//             // check if we need one more agglomeration here
//             auto CC = check_corner(msh, cl1, cl2);
//             if(CC.first)
//             {
//                 auto offset_added_cell = offset(msh, CC.second);

//                 auto MC_bis = merge_cells(msh, CC.second, cl1);
//                 loc_agglos.at(agglo_offset).add_cell(msh, CC.second, offset_added_cell, MC_bis.first);
//                 for(size_t i=0; i<MC_bis.second.size(); i++)
//                 {
//                     removed_faces.push_back( MC_bis.second[i] );
//                 }
//                 removed_cells.push_back(CC.second);
//             }

//             // end the merge procedure
//             auto MC = merge_cells(msh, cl2, loc_agglos.at(agglo_offset).new_cell);

//             loc_agglos.at(agglo_offset).add_cell(msh, cl2, offset_cl2, MC.first);

//             for(size_t i=0; i<MC.second.size(); i++)
//             {
//                 removed_faces.push_back( MC.second[i] );
//             }
//         }
//     }
    
//     //////////////////////   UPDATE THE MESH   ////////////////////////
//     size_t nb_cells_before = msh.cells_size();
//     size_t nb_cells_ok = 0;
//     size_t nb_cells_ko1 = 0;
//     size_t nb_cells_ko2 = 0;
//     size_t nb_cut_before = 0;
//     for (auto& cl : cells(msh)) {

//         if(locate(msh, cl) == location::ON_INTERFACE )
//             nb_cut_before++;
//         else
//             continue;

//         auto cl_id = offset(msh, cl);
//         auto& cl_cut_type = storage -> cut_cell_data[cl_id].cut;
//         if (cl_cut_type == cut_type::T_OK )
//             nb_cells_ok++;
//         else if ( cl_cut_type == cut_type::T_KO_NEG )
//             nb_cells_ko1++;
//         else if ( cl_cut_type == cut_type::T_KO_POS )
//             nb_cells_ko2++;
//         else
//             throw std::logic_error("We should not arrive here !!");
//     }

//     auto& cells = storage -> surfaces;
//     // Remove the agglomerated cells
//     typename std::vector<typename MeshType::cell_type>::iterator it_RC;
//     for(it_RC = removed_cells.begin(); it_RC != removed_cells.end(); it_RC++) {
//         cells.erase(std::remove(begin(cells), end(cells), *it_RC ), end(cells));
//     }

//     // Add new cells
//     for (size_t i = 0; i < loc_agglos.size(); i++) {
//         cells.push_back(loc_agglos.at(i).new_cell);
//     }

//     // Sort the new list of cells
//     std::sort(cells.begin(), cells.end());
    
//     // Remove faces
//     auto& faces = storage -> edges;
//     typename std::vector<typename MeshType::face_type>::iterator it_RF;
//     for(it_RF = removed_faces.begin(); it_RF != removed_faces.end(); it_RF++) {
//         faces.erase(std::remove(begin(faces), end(faces), *it_RF ), end(faces));
//     }

//     // sort the new list of faces
//     std::sort(faces.begin(), faces.end());

//     size_t nb_cells_after = msh.cells_size();
//     size_t nb_cut_after = 0;
//     for (auto& cl : msh) {
//         if (locate(msh, cl) == location::ON_INTERFACE)
//             nb_cut_after++;
//     }

//     ////////////  output some info
//     std::ofstream output_cells("output_cells.txt", std::ios::out | std::ios::trunc);
//     output_cells << " NB_CELLS_BEFORE = " << nb_cells_before << std::endl;
//     output_cells << " NB_CELLS_AFTER = " << nb_cells_after << std::endl;
//     output_cells << " NB_CUT_CELLS_BEFORE = " << nb_cut_before << std::endl;
//     output_cells << " NB_CUT_CELLS_AFTER = " << nb_cut_after << std::endl;

//     output_cells << " NB_CELLS_OK = " << nb_cells_ok << std::endl;
//     output_cells << " NB_CELLS_KO1 = " << nb_cells_ko1 << std::endl;
//     output_cells << " NB_CELLS_KO2 = " << nb_cells_ko2 << std::endl;

//     output_cells << " NB_CELLS_STEP_1 = " << nb_step1 << std::endl;
//     output_cells << " NB_CELLS_STEP_2 = " << nb_step2 << std::endl;
//     output_cells << " NB_CELLS_STEP_3 = " << nb_step3 << std::endl;

//     output_cells.close();
// }

}