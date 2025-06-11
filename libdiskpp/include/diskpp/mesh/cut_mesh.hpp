
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
    //bool case4 = location(msh, ns.front()) == where && location(msh, ns.back()) == where;

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

    size_t N = sqrt(msh.cells_size());

    //////////////////  face neighbors  ///////////////////
    // first row of cells -> look left
    for (size_t i = 1; i < N; i++) {

        auto cl1 = msh[i];
        auto cl1_id = offset(msh, cl1);
        auto& cl1_f_neighbors = storage -> cut_cell_data[cl1_id].f_neighbors;
        cl1_f_neighbors.insert(i-1);

        auto cl2 = msh[i-1];
        auto cl2_id = offset(msh, cl2);
        auto& cl2_f_neighbors = storage -> cut_cell_data[cl2_id].f_neighbors;
        cl2_f_neighbors.insert(i);

    }

    // other rows of cells
    for (size_t j = 1; j < N; j++)
    {
        // first cell of the row -> look bottom
        auto cl1 = msh[j*N];
        auto cl1_id = offset(msh, cl1);
        auto& cl1_f_neighbors = storage -> cut_cell_data[cl1_id].f_neighbors;
        cl1_f_neighbors.insert((j-1)*N);

        auto cl2 = msh[((j-1)*N)];
        auto cl2_id = offset(msh, cl2);
        auto& cl2_f_neighbors = storage -> cut_cell_data[cl2_id].f_neighbors;
        cl2_f_neighbors.insert( j*N );

        // other cells -> look left and bottom
        for (size_t i = 1; i < N; i++) {

            auto cl_c = msh[j*N + i]; // current
            auto cl_l = msh[j*N + i - 1]; // left
            auto cl_b = msh[(j-1)*N + i]; // bottom

            auto  cl_c_id = offset(msh, cl_c);
            auto& cl_c_f_neighbors = storage -> cut_cell_data[cl_c_id].f_neighbors;
            auto  cl_l_id = offset(msh, cl_l);
            auto& cl_l_f_neighbors = storage -> cut_cell_data[cl_l_id].f_neighbors;
            auto  cl_b_id = offset(msh, cl_b);
            auto& cl_b_f_neighbors = storage -> cut_cell_data[cl_b_id].f_neighbors;

            cl_c_f_neighbors.insert( j*N + i - 1 );
            cl_c_f_neighbors.insert( (j-1)*N + i );

            cl_l_f_neighbors.insert( j*N + i );
            cl_b_f_neighbors.insert( j*N + i );
        }
    }

    //////////////////////  diagonal neighbors  //////////////////////
    // first row of cells -> look left top
    for (size_t i = 1; i < N; i++) {
        auto cl1 = msh[i];
        auto cl1_id = offset(msh, cl1);
        auto& cl1_d_neighbors = storage -> cut_cell_data[cl1_id].f_neighbors;
        cl1_d_neighbors.insert(N + i-1);

        auto cl2 = msh[N + i-1];
        auto cl2_id = offset(msh, cl2);
        auto& cl2_d_neighbors = storage -> cut_cell_data[cl2_id].f_neighbors;
        cl2_d_neighbors.insert(i);
    }

    // other rows of cells
    for (size_t j = 1; j < N-1; j++)
    {
        // first cell of the row -> nothing to do
        // other cells -> look left top and left bottom
        for (size_t i = 1; i < N; i++)
        {
            auto cl_c = msh[j*N + i]; // current
            auto cl_l = msh[(j+1)*N + i - 1]; // left
            auto cl_b = msh[(j-1)*N + i - 1]; // bottom

            auto  cl_c_id = offset(msh, cl_c);
            auto& cl_c_d_neighbors = storage -> cut_cell_data[cl_c_id].f_neighbors;
            auto  cl_l_id = offset(msh, cl_l);
            auto& cl_l_d_neighbors = storage -> cut_cell_data[cl_l_id].f_neighbors;
            auto  cl_b_id = offset(msh, cl_b);
            auto& cl_b_d_neighbors = storage -> cut_cell_data[cl_b_id].f_neighbors;

            cl_c_d_neighbors.insert((j+1)*N + i - 1);
            cl_c_d_neighbors.insert((j-1)*N + i - 1);
            cl_l_d_neighbors.insert( j*N + i );
            cl_b_d_neighbors.insert( j*N + i );

        }
    }

    // last row -> look left bottom
    for (size_t i = 1; i < N; i++)
    {
        auto cl1 = msh[(N-1)*N + i];
        auto cl2 = msh[(N-2)*N + i-1];

        auto cl1_id = offset(msh, cl1);
        auto& cl1_d_neighbors = storage -> cut_cell_data[cl1_id].f_neighbors;

        auto cl2_id = offset(msh, cl2);
        auto& cl2_d_neighbors = storage -> cut_cell_data[cl2_id].f_neighbors;

        cl1_d_neighbors.insert((N-2)*N + i-1);
        cl2_d_neighbors.insert((N-1)*N + i);
    }
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

}