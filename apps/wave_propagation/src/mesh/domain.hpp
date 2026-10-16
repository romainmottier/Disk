//
//  domain.hpp
//
//  Mesh and media of the elasto-acoustic problems (namespace coupling_common, shared by the two formulations,
//  driven by the Lua parameters, parameters/lua_parameters.hpp): common types, Cartesian mesh of the rectangle
//  mesh.x x mesh.y with optional local refinement, FVCA5, PolyMesher or tagged polygonal mesh files; elastic cells
//  for x < 0 (or y < 0) and acoustic cells for x > 0 (or y > 0) with the interface x = 0 (or y = 0)
//  (problem.medium = "coupled"), media and materials given by the cell tags (problem.interface = "tags"), or a
//  single medium; fine region of the LTS / LI integrators.
//

#pragma once
#ifndef mesh_domain_hpp
#define mesh_domain_hpp

#include "../parameters/lua_parameters.hpp"
#include "diskpp/loaders/loader.hpp"
#include <random>
#include <cstdint>
#include <cstring>

namespace coupling_common {

    using RealType = double;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>   mesh_type;
    typedef typename mesh_type::point_type                                      point_type;
    typedef disk::BoundaryConditions<mesh_type, false>                          e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true>                           a_boundary_type;
    typedef Matrix<RealType, Dynamic, 1>                                        vector_type;

    const size_t bc_elastic_id  = 0;
    const size_t bc_acoustic_id = 1;

    // Polygonal mesh with one integer tag per cell (mesh.type = "polygons"), text file:
    //     # comment lines (optional)
    //     n_points n_cells
    //     x y                          (n_points lines)
    //     nv i_1 ... i_nv tag          (n_cells lines: vertex indices from 0, counterclockwise, then the cell tag)
    // A hanging node is a vertex of the polygons on both sides (the cell is a polygon with an extra vertex). The
    // boundary faces are the edges of a single cell; tags() gives the tag of each cell of the mesh after
    // move_to_mesh_storage (cells sorted as in the diskpp storage).
    template<typename T>
    class tagged_polygon_mesh_reader : public fitted_geometry_builder<disk::generic_mesh<T, 2>> {
        typedef disk::generic_mesh<T, 2>         mesh_type;
        typedef typename mesh_type::point_type   point_type;
        typedef typename mesh_type::node_type    node_type;
        typedef typename mesh_type::surface_type surface_type;
        struct polygon {
            std::vector<size_t>             m_member_nodes;
            std::set<std::array<size_t, 2>> m_member_edges;
        };
        std::string                         m_file;
        std::vector<point_type>             m_points;
        std::vector<node_type>              m_nodes;
        std::vector<std::array<size_t, 2>>  m_facets, m_boundary;
        std::vector<polygon>                m_polygons;
        std::vector<int>                    m_file_tags, m_tags;
    public:
        explicit tagged_polygon_mesh_reader(const std::string& file) : m_file(file) {}

        bool build_mesh() override {
            std::ifstream in(m_file);
            if (!in.good())
                throw std::runtime_error("cannot read the polygonal mesh " + m_file);
            std::string line;
            auto next = [&](std::string& l) {
                while (std::getline(in, l))
                    if (!l.empty() && l[0] != '#') return true;
                return false;
            };
            size_t n_points = 0, n_cells = 0;
            if (!next(line) || !(std::istringstream(line) >> n_points >> n_cells))
                throw std::runtime_error(m_file + ": header \"n_points n_cells\" expected");
            for (size_t i = 0; i < n_points; i++) {
                T x, y;
                if (!next(line) || !(std::istringstream(line) >> x >> y))
                    throw std::runtime_error(m_file + ": point " + std::to_string(i) + " expected");
                m_points.push_back(point_type(x, y));
                m_nodes.push_back(node_type(disk::point_identifier<2>(i)));
            }
            std::map<std::array<size_t, 2>, size_t> count;
            for (size_t c = 0; c < n_cells; c++) {
                if (!next(line))
                    throw std::runtime_error(m_file + ": cell " + std::to_string(c) + " expected (nv i_1 ... i_nv tag)");
                std::istringstream is(line);
                size_t nv = 0;
                if (!(is >> nv) || nv < 3)
                    throw std::runtime_error(m_file + ": cell " + std::to_string(c) + " expected (nv i_1 ... i_nv tag)");
                polygon pg;
                for (size_t i = 0; i < nv; i++) {
                    size_t id;
                    if (!(is >> id) || id >= n_points)
                        throw std::runtime_error(m_file + ": bad vertex in cell " + std::to_string(c));
                    pg.m_member_nodes.push_back(id);
                }
                int tag = 0;
                is >> tag;
                for (size_t i = 0; i < nv; i++) {
                    std::array<size_t, 2> e = {pg.m_member_nodes[i], pg.m_member_nodes[(i + 1) % nv]};
                    if (e[0] > e[1]) std::swap(e[0], e[1]);
                    pg.m_member_edges.insert(e);
                    count[e]++;
                }
                m_polygons.push_back(pg);
                m_file_tags.push_back(tag);
            }
            for (auto& ec : count) {
                m_facets.push_back(ec.first);
                if (ec.second == 1) m_boundary.push_back(ec.first);
                else if (ec.second > 2)
                    throw std::runtime_error(m_file + ": an edge belongs to more than two cells");
            }
            return true;
        }

        void move_to_mesh_storage(mesh_type& msh) override {
            std::vector<surface_type> surfaces;
            if (!fitted_geometry_builder<mesh_type>::fill_mesh_storage(msh, m_points, m_nodes, m_facets, m_boundary, m_polygons, surfaces))
                throw std::runtime_error(m_file + ": inconsistent mesh");
            std::vector<size_t> order(surfaces.size());
            for (size_t i = 0; i < order.size(); i++) order[i] = i;
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return surfaces[a] < surfaces[b]; });
            std::vector<surface_type> sorted;
            sorted.reserve(surfaces.size());
            m_tags.clear();
            for (size_t i : order) {
                sorted.push_back(surfaces[i]);
                m_tags.push_back(m_file_tags[i]);
            }
            msh.backend_storage()->surfaces = std::move(sorted);
        }

        const std::vector<int>& tags() const { return m_tags; }
    };

    // Cartesian mesh of the rectangle with (nx 2^level) x (ny 2^level) cells. The cells whose barycenter is within one
    // cell size (in x and in y) of a point of mesh.refine_points are refined mesh.refine_levels times (hanging nodes).
    // FVCA5 or PolyMesher meshes: the file mesh.files[level] (Lua indexing from 1); PolyMesher files are read with
    // polygon_2d_mesh_reader, the reader of the polygonal meshes of Burman-Duran-Ern 2022.
    inline void build_mesh(mesh_type& msh, const wave_parameters::mesh_data& m, size_t level, std::vector<int>* tags = nullptr) {
        if (m.type == "polygons") {
            tagged_polygon_mesh_reader<RealType> reader(m.files.at(level - 1));
            reader.build_mesh();
            reader.move_to_mesh_storage(msh);
            if (tags) *tags = reader.tags();
            return;
        }
        if (m.type == "polymesher") {
            const std::string& file = m.files.at(level - 1);
            if (!std::ifstream(file).good())
                throw std::runtime_error("cannot read the PolyMesher mesh " + file);
            polygon_2d_mesh_reader<RealType> reader;
            reader.set_poly_mesh_file(file);
            reader.build_mesh();
            reader.move_to_mesh_storage(msh);
            return;
        }
        if (m.type == "fvca5") {
            const std::string& file = m.files.at(level - 1);
            if (!disk::load_mesh_fvca5_2d<RealType>(file.c_str(), msh))
                throw std::runtime_error("cannot read the FVCA5 mesh " + file);
            return;
        }
        const RealType lx = m.x[1] - m.x[0], ly = m.y[1] - m.y[0];
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, m.nx, m.ny);
        mesh_builder.refine_mesh(level);
        mesh_builder.set_translation_data(m.x[0], m.y[0]);
        mesh_builder.build_mesh();
        if (m.refine_levels > 0 && !m.refine_box.empty()) {
            // graded refinement of the box (as in the LTS stability prototypes, refine_with_protection)
            const std::vector<double> b = m.refine_box;
            mesh_builder.refine_with_protection([b](const point_type& pt) {
                return pt.x() > b[0] && pt.x() < b[1] && pt.y() > b[2] && pt.y() < b[3];
            }, int(m.refine_levels));
        }
        else if (m.refine_levels > 0) {
            const RealType hx = lx / (RealType(m.nx) * std::pow(2.0, RealType(level)));
            const RealType hy = ly / (RealType(m.ny) * std::pow(2.0, RealType(level)));
            std::vector<size_t> cells;
            size_t n_cells = mesh_builder.get_nx() * mesh_builder.get_ny();
            for (size_t i = 0; i < n_cells; i++) {
                auto bar = mesh_builder.cell_barycenter(i);
                for (auto& pt : m.refine_points)
                    if (std::abs(bar.x() - pt[0]) < hx && std::abs(bar.y() - pt[1]) < hy) {
                        cells.push_back(i);
                        break;
                    }
            }
            mesh_builder.refine_cells(cells, m.refine_levels);
        }
        mesh_builder.move_to_mesh_storage(msh);
    }

    // Materials, interface and Dirichlet boundary ids of the coupled domain
    struct domain_data {
        std::map<size_t, elastic_material_data<RealType>>  e_material;
        std::map<size_t, acoustic_material_data<RealType>> a_material;
        std::set<size_t> interface_face_indexes;
        std::map<size_t, std::pair<size_t, size_t>> interface_cell_pair_indexes;
        RealType h_max = 0.0;
        RealType h_min = 0.0;
        RealType c_max = 0.0;       // largest wave speed of the materials (CFL)
    };

    // Elastic cells (rho = 1, vp = sqrt(3), vs = 1: mu = lambda = 1) and acoustic cells (rho = 1, vp = 1), the values of
    // the manufactured solutions: x < 0 elastic and x > 0 acoustic for medium = "coupled", with the interface faces
    // (x = 0) and their (elastic, acoustic) cell pair; all the cells of the medium otherwise. Boundary faces tagged
    // bc_elastic_id / bc_acoustic_id.
    inline domain_data classify_domain(mesh_type& msh, const std::string& medium) {
        domain_data dd;
        const RealType eps = 1.0e-10;
        auto acoustic_cell = [&](const point_type& bar) {
            return medium == "acoustic" || (medium == "coupled" && bar.x() > 0);
        };
        if (medium == "coupled") {
            for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
                const auto face = *face_it;
                if (std::fabs(barycenter(msh, face).x()) < eps)
                    dd.interface_face_indexes.insert(msh.lookup(face));
            }
        }
        dd.h_max = 0.0;
        dd.h_min = 1.0e10;
        for (auto& cell : msh) {
            auto cell_ind = msh.lookup(cell);
            auto bar = barycenter(msh, cell);
            RealType h = diameter(msh, cell);
            dd.h_max = std::max(dd.h_max, h);
            dd.h_min = std::min(dd.h_min, h);
            if (acoustic_cell(bar))
                dd.a_material.insert(std::make_pair(cell_ind, acoustic_material_data<RealType>(1.0, 1.0)));
            else
                dd.e_material.insert(std::make_pair(cell_ind, elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0)));
            for (auto& face : faces(msh, cell)) {
                auto fc_id = msh.lookup(face);
                if (dd.interface_face_indexes.count(fc_id)) {
                    if (bar.x() > 0)
                        dd.interface_cell_pair_indexes[fc_id].second = cell_ind;
                    else
                        dd.interface_cell_pair_indexes[fc_id].first = cell_ind;
                }
            }
        }
        for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
            auto face = *face_it;
            auto fc_id = msh.lookup(face);
            size_t id = acoustic_cell(barycenter(msh, face)) ? bc_acoustic_id : bc_elastic_id;
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{id, true};
        }
        return dd;
    }

    // Media and materials from the parameters: problem.interface = "x" or "y" (elastic cells x < 0 or y < 0 with
    // problem.elastic = {rho, vp, vs}, acoustic cells x > 0 or y > 0 with problem.acoustic = {rho, vp}) for a coupled
    // problem, all cells of the medium otherwise; "tags": medium and material of each cell from its mesh tag
    // (problem.zones). Interface faces: faces between an elastic and an acoustic cell, with their (elastic, acoustic)
    // cell pair; boundary faces tagged bc_elastic_id / bc_acoustic_id from the medium of their cell. For
    // problem.interface = "x" and the default materials, the same data as classify_domain(msh, medium).
    inline domain_data classify_domain(mesh_type& msh, const wave_parameters& prm, const std::vector<int>& tags = {}) {
        const auto& pb = prm.problem;
        domain_data dd;
        const size_t n_cells = msh.cells_size();
        std::vector<char> acoustic(n_cells, 0);
        std::vector<RealType> rho(n_cells), vp(n_cells), vs(n_cells);
        for (auto& cell : msh) {
            auto ci = msh.lookup(cell);
            auto bar = barycenter(msh, cell);
            if (pb.interface == "tags") {
                if (ci >= tags.size() || !pb.zones.count(tags[ci]))
                    throw std::runtime_error("problem.zones: no zone for the tag " + (ci < tags.size() ? std::to_string(tags[ci]) : std::string("?"))
                                             + " of cell " + std::to_string(ci));
                const auto& z = pb.zones.at(tags[ci]);
                acoustic[ci] = z.medium == "acoustic";
                rho[ci] = z.rho; vp[ci] = z.vp; vs[ci] = z.vs;
                continue;
            }
            if (pb.medium == "acoustic") acoustic[ci] = 1;
            else if (pb.medium == "elastic") acoustic[ci] = 0;
            else acoustic[ci] = (pb.interface == "x" ? bar.x() : bar.y()) > 0;
            if (acoustic[ci]) { rho[ci] = pb.rho_acoustic; vp[ci] = pb.vp_acoustic; vs[ci] = 0.0; }
            else              { rho[ci] = pb.rho_elastic;  vp[ci] = pb.vp_elastic;  vs[ci] = pb.vs_elastic; }
        }
        dd.h_max = 0.0;
        dd.h_min = 1.0e10;
        for (auto& cell : msh) {
            auto ci = msh.lookup(cell);
            RealType h = diameter(msh, cell);
            dd.h_max = std::max(dd.h_max, h);
            dd.h_min = std::min(dd.h_min, h);
            dd.c_max = std::max(dd.c_max, vp[ci]);
            if (acoustic[ci])
                dd.a_material.insert(std::make_pair(ci, acoustic_material_data<RealType>(rho[ci], vp[ci])));
            else
                dd.e_material.insert(std::make_pair(ci, elastic_material_data<RealType>(rho[ci], vp[ci], vs[ci])));
        }
        // interface faces: an elastic and an acoustic cell on either side
        std::map<size_t, std::vector<size_t>> face_cells;
        for (auto& cell : msh)
            for (auto& face : faces(msh, cell))
                face_cells[msh.lookup(face)].push_back(msh.lookup(cell));
        for (auto& fc : face_cells) {
            if (fc.second.size() != 2 || acoustic[fc.second[0]] == acoustic[fc.second[1]])
                continue;
            size_t e = acoustic[fc.second[0]] ? fc.second[1] : fc.second[0];
            size_t a = acoustic[fc.second[0]] ? fc.second[0] : fc.second[1];
            dd.interface_face_indexes.insert(fc.first);
            dd.interface_cell_pair_indexes[fc.first] = std::make_pair(e, a);
        }
        for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
            auto fc_id = msh.lookup(*face_it);
            const auto& cells = face_cells[fc_id];
            size_t id = (!cells.empty() && acoustic[cells[0]]) ? bc_acoustic_id : bc_elastic_id;
            msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{id, true};
        }
        return dd;
    }

    // Multilevel LTS (time.scheme = "mlts", both formulations): level of each cell (indexed by mesh cell id) from
    // x_T = h_max / h_T (time.level_measure = "h") or x_T = max_K (h_K / c_K) / (h_T / c_T) (time.level_measure = "h/c",
    // c the P-wave speed of the cell: heterogeneous media), ceil(log_r x_T) for the level ratio r = time.level_ratio
    // (time.level_ratios p_0, p_1, ...: first l with p_0 ... p_{l-1} >= x_T, the last ratio repeated); then, from the
    // finest level down, each level annexes time.overlap rings of vertex neighbours from the coarser ones (the overlap
    // region of the preprint; time.overlaps: o_l rings for the level l >= 1, the last one repeated). Returns the number of
    // levels.
    inline size_t mlts_cell_levels(mesh_type& msh, const domain_data& dd, const wave_parameters& prm, std::vector<int>& level) {
        auto storage = msh.backend_storage();
        const size_t nc = msh.cells_size();
        level.assign(nc, 0);
        const RealType log_r = std::log(RealType(prm.time.level_ratio));
        const std::vector<size_t>& ratios = prm.time.level_ratios;
        // h_T / c_T (level_measure "h/c")
        const bool by_cfl = prm.time.level_measure == "h/c";
        // size of a cell: its diameter, or 4 |T| / |dT| (time.level_size = "area": the side of a square, smaller for a
        // thin cell, e.g. a piece of a cell cut along a fitted interface, whose stiffness the diameter misses)
        const bool by_area = prm.time.level_size == "area" || prm.time.level_size == "width";
        auto size_of = [&](size_t c) -> RealType {
            const auto& cl = storage->surfaces[c];
            if (!by_area) return diameter(msh, cl);
            if (prm.time.level_size == "width") {
                // smallest width: min over the edges of the largest distance of the vertices to the line of the edge
                auto pts = points(msh, cl);
                RealType best = std::numeric_limits<RealType>::max();
                for (size_t i = 0; i < pts.size(); i++) {
                    const auto a = pts[i], b = pts[(i + 1) % pts.size()];
                    const RealType tx = b.x() - a.x(), ty = b.y() - a.y(), L = std::sqrt(tx * tx + ty * ty);
                    if (L == 0) continue;
                    RealType far = 0;
                    for (auto& q : pts) far = std::max(far, std::abs(-(q.x() - a.x()) * ty + (q.y() - a.y()) * tx) / L);
                    best = std::min(best, far);
                }
                return best;
            }
            RealType per = 0;
            for (auto& fc : faces(msh, cl)) per += measure(msh, fc);
            return 4 * measure(msh, cl) / per;
        };
        std::vector<RealType> s(by_cfl || by_area ? nc : 0);
        RealType s_max = 0;
        if (by_cfl || by_area) {
            for (size_t c = 0; c < nc; c++) {
                RealType vp = 1;
                if (by_cfl) {
                    auto ia = dd.a_material.find(c);
                    if (ia != dd.a_material.end()) { auto m = ia->second; vp = m.vp(); }
                    else { auto m = dd.e_material.at(c); vp = m.vp(); }
                }
                s[c] = size_of(c) / vp;
                s_max = std::max(s_max, s[c]);
            }
        }
        const RealType tol = RealType(1) + RealType(prm.time.level_tolerance);
        auto ratio_of = [&](size_t c) {
            const RealType x = (by_cfl || by_area) ? s_max / s[c] : dd.h_max / diameter(msh, storage->surfaces[c]);
            return prm.time.level_tolerance > 0 ? std::max(RealType(1), x / tol) : x;
        };
        int finest = 0;
        for (size_t c = 0; c < nc; c++) {
            if (ratios.empty()) {
                const RealType x = std::log(ratio_of(c)) / log_r;
                level[c] = std::max(0, int(std::ceil(x - 1.0e-6)));
            }
            else {
                // first l with log(p_0 ... p_{l-1}) >= log(x_T), the last ratio repeated
                const RealType x = std::log(ratio_of(c)) - 1.0e-6 * std::log(2.0);
                RealType logR = 0;
                int l = 0;
                while (logR < x) {
                    logR += std::log(RealType(ratios[std::min(size_t(l), ratios.size() - 1)]));
                    l++;
                }
                level[c] = l;
            }
            finest = std::max(finest, level[c]);
        }
        std::vector<std::vector<size_t>> point_cells(msh.points_size());
        for (size_t c = 0; c < nc; c++)
            for (auto pid : storage->surfaces[c].point_ids()) point_cells[size_t(pid)].push_back(c);
        const std::vector<size_t>& ov = prm.time.overlaps;
        for (int l = finest; l >= 1; l--)
            for (size_t ring = 0; ring < (ov.empty() ? prm.time.overlap : ov[std::min(size_t(l - 1), ov.size() - 1)]); ring++) {
                std::vector<size_t> annexed;
                for (size_t c = 0; c < nc; c++)
                    if (level[c] >= l)
                        for (auto pid : storage->surfaces[c].point_ids())
                            for (size_t n : point_cells[size_t(pid)])
                                if (level[n] < l) annexed.push_back(n);
                for (size_t n : annexed) level[n] = l;
            }
        return size_t(finest) + 1;
    }

    // Finite value, from the exponent bits: std::isfinite may be folded to true under -ffast-math (OPT_AGGRESSIVE_FP), the
    // default of diskpp, while an unstable run of the CFL studies ends with inf or NaN
    inline bool finite_value(double v) {
        std::uint64_t b;
        std::memcpy(&b, &v, sizeof(b));
        return ((b >> 52) & 0x7ffu) != 0x7ffu;
    }

    // Ratios p_0 .. p_{L-2} of the L levels of the multilevel LTS: time.level_ratio, or time.level_ratios with its last
    // ratio repeated
    inline std::vector<size_t> mlts_ratios(const wave_parameters& prm, size_t L) {
        const std::vector<size_t>& r = prm.time.level_ratios;
        std::vector<size_t> p(L > 0 ? L - 1 : 0, prm.time.level_ratio);
        if (!r.empty())
            for (size_t l = 0; l < p.size(); l++) p[l] = r[std::min(l, r.size() - 1)];
        return p;
    }

    // Fine region of the LTS / LI integrators: round(h_max/h_min) substeps (or time.substeps), cells smaller than
    // 0.75 h_max (1.25 h_max on a uniform mesh, everything fine) or than time.fine_size
    inline size_t fine_substeps(const wave_parameters& prm, const domain_data& dd) {
        if (prm.time.substeps > 0) return prm.time.substeps;
        size_t p = static_cast<size_t>(std::round(dd.h_max / dd.h_min));
        return p < 1 ? 1 : p;
    }

    inline RealType fine_size(const wave_parameters& prm, const domain_data& dd) {
        if (prm.time.fine_size > 0.0) return prm.time.fine_size;
        size_t p = static_cast<size_t>(std::round(dd.h_max / dd.h_min));
        return p <= 1 ? 1.25 * dd.h_max : 0.75 * dd.h_max;
    }

}

#endif /* mesh_domain_hpp */
