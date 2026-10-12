//
//  coupling_common.hpp
//
//  Setup shared by the two elasto-acoustic prototypes (HHOFirstOrder, HHOSecondOrder), driven by the Lua parameters
//  (common/lua_parameters.hpp): Cartesian mesh of the rectangle mesh.x x mesh.y with optional local refinement
//  (fine region of the LTS / LI integrators), FVCA5, PolyMesher or tagged polygonal mesh files; elastic cells for x < 0
//  (or y < 0) and acoustic cells for x > 0 (or y > 0) with the interface x = 0 (or y = 0) (problem.medium = "coupled"),
//  media and materials given by the cell tags (problem.interface = "tags"), or a single medium; Dirichlet data and
//  sources from a manufactured solution; convergence tables.
//

#pragma once
#ifndef coupling_common_hpp
#define coupling_common_hpp

#include "../../common/lua_parameters.hpp"
#include "diskpp/loaders/loader.hpp"

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
        if (m.refine_levels > 0) {
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

    // Manufactured solution of problem.solution (wave_parameters::solution_names)
    inline scal_vec_analytic_functions::EFunctionType solution_type(const std::string& name) {
        typedef scal_vec_analytic_functions::EFunctionType F;
        static const std::map<std::string, F> types = {
            {"interface_coupling",          F::EFunctionInterfaceCoupling},
            {"non_polynomial",              F::EFunctionNonPolynomial},
            {"acoustic_quadratic_in_time",  F::reproduction_acoustic},
            {"acoustic_quadratic_in_space", F::EFunctionAcousticQuadraticInSpace},
            {"acoustic_non_polynomial",     F::EFunctionAcousticNonPolynomial},
            {"acoustic_quadratic_in_space_cos", F::EFunctionAcousticQuadraticInSpaceCos},
            {"acoustic_non_polynomial_cos", F::EFunctionAcousticNonPolynomialCos},
            {"elastic_quadratic_in_time",   F::reproduction_elastic},
            {"elastic_quadratic_in_space",  F::EFunctionElasticQuadraticInSpace},
            {"elastic_non_polynomial",      F::EFunctionElasticNonPolynomial},
            {"coupled_sinusoidal",          F::EFunctionNonPolynomial_paper},
            {"coupled_poly_in_space",       F::EFunctionQuadraticInSpace},
            {"coupled_poly_in_time",        F::EFunctionQuadraticInTime},
            {"ricker",                      F::EFunctionNonPolynomial}};       // data replaced by zero (HHOFirstOrder)
        return types.at(name);
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

    // One row of a convergence table: step sizes (h = largest, h_min = smallest cell diameter) and errors
    struct conv_row {
        RealType h, dt;
        std::vector<RealType> errors;
        RealType h_min = 0.0;
    };

    // Prints the errors and the observed orders log(e_{i-1}/e_i) / log(r_{i-1}/r_i), with r = h or dt
    inline void print_table(const std::vector<std::string>& names, const std::vector<conv_row>& rows,
                            bool rate_in_dt, std::ostream& os) {
        auto width = [&](size_t j) { return int(std::max<size_t>(16, names[j].size() + 2)); };
        os << std::setw(12) << "h" << std::setw(12) << "dt";
        for (size_t j = 0; j < names.size(); j++)
            os << std::setw(width(j)) << names[j] << std::setw(8) << "order";
        os << std::endl;
        for (size_t i = 0; i < rows.size(); i++) {
            os << std::setw(12) << std::setprecision(4) << rows[i].h << std::setw(12) << std::setprecision(4) << rows[i].dt;
            for (size_t j = 0; j < rows[i].errors.size(); j++) {
                os << std::setw(width(j)) << std::setprecision(6) << rows[i].errors[j];
                if (i > 0 && rows[i].errors[j] > 0 && rows[i-1].errors[j] > 0) {
                    RealType r = rate_in_dt ? rows[i-1].dt / rows[i].dt : rows[i-1].h / rows[i].h;
                    os << std::setw(8) << std::setprecision(3) << std::log(rows[i-1].errors[j] / rows[i].errors[j]) / std::log(r);
                }
                else
                    os << std::setw(8) << "--";
            }
            os << std::endl;
        }
        os.flush();
    }

    // Data table: "# h h_min dt name..." header, then one row per level (errors without orders), full precision
    inline void write_table(const std::string& file, const std::vector<std::string>& names, const std::vector<conv_row>& rows) {
        std::ofstream os(file);
        os << "# h h_min dt";
        for (auto& n : names)
            os << " " << n;
        os << std::endl << std::setprecision(16) << std::scientific;
        for (auto& r : rows) {
            os << r.h << " " << (r.h_min > 0.0 ? r.h_min : r.h) << " " << r.dt;
            for (auto e : r.errors)
                os << " " << e;
            os << std::endl;
        }
    }

    // Energy-type norm ||x||_M = sqrt(x^T M x) of a vector of unknowns for a mass matrix M
    inline RealType mass_norm(const SparseMatrix<RealType>& M, const vector_type& x) {
        return std::sqrt(std::abs(x.dot(M * x)));
    }

    // Levels and time steps of the study: "run" (mesh.level, n_steps), "space" (mesh.level .. last_level, n_steps
    // doubled with each level if study.refine_dt, fixed otherwise), "time" (mesh.level, see the prototypes)
    inline std::vector<size_t> study_levels(const wave_parameters& prm) {
        std::vector<size_t> levels;
        size_t last = prm.study.type == "space" || prm.study.type == "spectrum" ? prm.study.last_level : prm.mesh.level;
        for (size_t l = prm.mesh.level; l <= last; l++)
            levels.push_back(l);
        return levels;
    }

    // time.dt_h > 0: dt = dt_h h_min on each mesh (smallest cell diameter, as the h of Burman-Duran-Ern 2022 on
    // polygonal meshes), with as many steps as needed to reach the final time
    inline size_t study_steps(const wave_parameters& prm, size_t level, RealType h_min) {
        if (prm.time.dt_h > 0.0)
            return size_t(std::ceil(prm.time.final_time / (prm.time.dt_h * h_min) - 1.0e-9));
        if (!prm.study.refine_dt) return prm.time.n_steps;
        if (prm.time.steps_growth != 2.0)
            return std::max<size_t>(1, size_t(std::llround(RealType(prm.time.n_steps) *
                                                            std::pow(prm.time.steps_growth, RealType(level - prm.mesh.level)))));
        return prm.time.n_steps * (size_t(1) << (level - prm.mesh.level));
    }

    // Columns of the convergence tables present for the medium of the problem: the names of the elastic, then of the
    // acoustic errors, filtered by problem.medium; select picks the matching entries of a full error vector
    struct error_columns {
        std::vector<std::string> names;
        std::vector<size_t> index;
        std::vector<RealType> select(const std::vector<RealType>& all) const {
            std::vector<RealType> e;
            for (size_t i : index) e.push_back(all[i]);
            return e;
        }
    };

    inline error_columns columns(const wave_parameters& prm, const std::vector<std::string>& elastic,
                                 const std::vector<std::string>& acoustic) {
        error_columns c;
        for (size_t i = 0; i < elastic.size(); i++)
            if (prm.has_elastic()) { c.names.push_back(elastic[i]); c.index.push_back(i); }
        for (size_t i = 0; i < acoustic.size(); i++)
            if (prm.has_acoustic()) { c.names.push_back(acoustic[i]); c.index.push_back(elastic.size() + i); }
        return c;
    }

}

#endif /* coupling_common_hpp */
