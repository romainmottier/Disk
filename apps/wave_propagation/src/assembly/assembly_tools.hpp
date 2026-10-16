//
//  assembly_tools.hpp
//
//  Tools of the fast mode shared by the assemblers (elastoacoustic_four_fields_assembler,
//  elastoacoustic_two_fields_assembler): local computations of a loop over the cells done by the OpenMP threads and
//  scattered in the order of the serial loop (ordered_loop), and the data kept by set_cache_rhs_data for the
//  right-hand sides assembled at every stage: quadrature points, weights and basis values of the source terms,
//  assembly indices and Dirichlet faces (quadrature data, factorized face mass) of the cells with boundary faces.
//  The cached right-hand sides use the arithmetic of disk::make_rhs and of scatter_e_bc_data / scatter_a_bc_data
//  (elastoacoustic_assembler_base): the same vectors, bit for bit. Geometric classes of the cells (cell_classes):
//  translated copies of a cell with the same material have the same local operators up to round-off, computed once
//  per class (structured meshes: a few classes).
//

#pragma once
#ifndef assembly_tools_hpp
#define assembly_tools_hpp

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

#include "assembly_index.hpp"
#include "../solvers/parallel_sparse.hpp"

template<typename Mesh, typename T>
struct assembly_tools {

    typedef decltype(disk::make_vector_monomial_basis(std::declval<const Mesh&>(), std::declval<const typename Mesh::cell_type&>(), size_t(0))) e_cell_basis_type;
    typedef decltype(disk::make_scalar_monomial_basis(std::declval<const Mesh&>(), std::declval<const typename Mesh::cell_type&>(), size_t(0))) a_cell_basis_type;
    typedef decltype(disk::make_vector_monomial_basis(std::declval<const Mesh&>(), std::declval<const typename Mesh::face_type&>(), size_t(0))) e_face_basis_type;
    typedef decltype(disk::make_scalar_monomial_basis(std::declval<const Mesh&>(), std::declval<const typename Mesh::face_type&>(), size_t(0))) a_face_basis_type;
    template<typename Basis>
    using phi_type = decltype(std::declval<Basis>().eval_functions(std::declval<typename Mesh::point_type>()));

    // Quadrature points, weights and basis values on a cell or a face
    template<typename Phi>
    struct quadrature {
        std::vector<typename Mesh::point_type> points;
        std::vector<T> weights;
        std::vector<Phi> phi;
    };
    template<typename Phi>
    struct bc_face_data {
        size_t face_i, fc_id;
        quadrature<Phi> quad;
        Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> llt;
    };
    template<typename Phi>
    struct bc_cell_data {
        std::vector<assembly_index> asm_map;
        std::vector<bc_face_data<Phi>> faces;
    };
    typedef std::vector<quadrature<phi_type<e_cell_basis_type>>> e_rhs_data;
    typedef std::vector<quadrature<phi_type<a_cell_basis_type>>> a_rhs_data;
    typedef std::vector<bc_cell_data<phi_type<e_face_basis_type>>> e_bc_data;
    typedef std::vector<bc_cell_data<phi_type<a_face_basis_type>>> a_bc_data;

    // local(i) for i = 0..n-1, then scatter(i, local(i)) in increasing i: the serial loop; in parallel, the local
    // matrices and vectors of a chunk of cells are computed by the OpenMP threads, then scattered in the same order
    // (same triplets and same sums: the same global matrices and vectors, bit for bit)
    template<typename Local, typename Scatter>
    static void ordered_loop(bool parallel, size_t n, Local&& local, Scatter&& scatter) {
        if (!parallel) {
            for (size_t i = 0; i < n; i++) scatter(i, local(i));
            return;
        }
        typedef decltype(local(size_t(0))) result_type;
        const size_t chunk = 1024;
        std::vector<result_type> res;
        for (size_t c0 = 0; c0 < n; c0 += chunk) {
            const size_t nc = std::min(chunk, n - c0);
            res.assign(nc, result_type());
            #pragma omp parallel for schedule(dynamic, 4)
            for (long i = 0; i < (long)nc; i++) res[i] = local(c0 + size_t(i));
            for (size_t i = 0; i < nc; i++) scatter(c0 + i, res[i]);
        }
    }

    // Geometric classes of the cells cell_ids: diameter h_T (logarithm rounded to 1e-9: cells of different sizes are in
    // different classes), points of the cell (cell order) and of its faces (faces(msh, cell) order, face point order)
    // relative to its barycenter, rounded to 1e-9 h_T (relative to the cell: the same resolution for the smallest cells of
    // strongly graded meshes as for the largest ones), and the material material(i) (exact values). class_of[i]: class
    // of the i-th cell; representatives: first cell of each class.
    template<typename Material>
    static std::vector<size_t> cell_classes(const Mesh& msh, const std::vector<size_t>& cell_ids, Material&& material,
                                            std::vector<size_t>& representatives) {
        auto storage = msh.backend_storage();
        const size_t n = cell_ids.size();
        std::vector<std::vector<long long>> sig(n);
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)n; i++) {
            auto& cell = storage->surfaces[cell_ids[i]];
            const auto bc = barycenter(msh, cell);
            const T h_T = T(diameter(msh, cell));
            const T q = T(1e-9) * h_T;       // far above the round-off of the coordinates, far below the differences of shapes
            auto& s = sig[i];
            s.push_back(std::llround(std::log(h_T) * T(1e9)));
            auto put = [&](const typename Mesh::point_type& p) {
                for (size_t d = 0; d < Mesh::dimension; d++) s.push_back(std::llround((p[d] - bc[d]) / q));
            };
            auto pts = points(msh, cell);
            s.push_back((long long)pts.size());
            for (auto& p : pts) put(p);
            auto fcs = faces(msh, cell);
            s.push_back((long long)fcs.size());
            for (auto& fc : fcs)
                for (auto& p : points(msh, fc)) put(p);
            for (T m : material(size_t(i))) {
                long long bits;
                std::memcpy(&bits, &m, sizeof(bits));
                s.push_back(bits);
            }
        }
        std::map<std::vector<long long>, size_t> index;
        std::vector<size_t> class_of(n);
        representatives.clear();
        for (size_t i = 0; i < n; i++) {
            auto ins = index.emplace(std::move(sig[i]), representatives.size());
            if (ins.second) representatives.push_back(i);
            class_of[i] = ins.first->second;
        }
        return class_of;
    }

    // Classes of the elastic and acoustic cells of an assembler (m_e_material, m_a_material order)
    struct classes {
        std::vector<size_t> e_ids, e_class_of, e_rep, a_ids, a_class_of, a_rep;
        std::map<size_t, size_t> e_pos, a_pos;    // position of a cell id in e_ids, a_ids
    };
    template<typename EMaterial, typename AMaterial>
    static classes make_classes(const Mesh& msh, const EMaterial& e_material, const AMaterial& a_material) {
        classes c;
        std::vector<std::vector<T>> e_keys, a_keys;
        for (auto& e : e_material) {
            c.e_pos[e.first] = c.e_ids.size();
            c.e_ids.push_back(e.first);
            auto m = e.second;
            e_keys.push_back({m.rho(), m.vp(), m.vs()});
        }
        for (auto& a : a_material) {
            c.a_pos[a.first] = c.a_ids.size();
            c.a_ids.push_back(a.first);
            auto m = a.second;
            a_keys.push_back({m.rho(), m.vp()});
        }
        c.e_class_of = cell_classes(msh, c.e_ids, [&](size_t i) { return e_keys[i]; }, c.e_rep);
        c.a_class_of = cell_classes(msh, c.a_ids, [&](size_t i) { return a_keys[i]; }, c.a_rep);
        return c;
    }

    // op(i) for the representative i of every class (in parallel); ops[class]
    template<typename Op>
    static std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> class_operators(const std::vector<size_t>& rep, Op&& op) {
        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> ops(rep.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (long c = 0; c < (long)rep.size(); c++) ops[c] = op(rep[c]);
        return ops;
    }

    // parallel_sparse::set_from_triplets (the matrix of setFromTriplets, bit for bit, built by the OpenMP threads)
    static void set_from_triplets(Eigen::SparseMatrix<T>& A, const std::vector<Eigen::Triplet<T>>& trips) {
        parallel_sparse::set_from_triplets(A, trips);
    }

    // Quadrature data of degree 2 deg (disk::make_rhs) of the cells keys
    template<typename Quadrature, typename MakeBasis>
    static void cell_quadratures(const Mesh& msh, const std::vector<size_t>& keys, size_t deg, MakeBasis make_basis,
                                 std::vector<Quadrature>& out, bool parallel) {
        auto storage = msh.backend_storage();
        out.assign(keys.size(), {});
        #pragma omp parallel for schedule(dynamic, 16) if (parallel)
        for (long i = 0; i < (long)keys.size(); i++) {
            auto& cell = storage->surfaces[keys[i]];
            auto cell_basis = make_basis(cell);
            for (auto& qp : integrate(msh, cell, 2 * deg)) {
                out[i].points.push_back(qp.point());
                out[i].weights.push_back(qp.weight());
                out[i].phi.push_back(cell_basis.eval_functions(qp.point()));
            }
        }
    }

    // Number of doubles of the quadrature data of cell_quadratures (basis of n_components x basis_size values)
    static size_t quadrature_doubles(const Mesh& msh, const std::vector<size_t>& keys, size_t deg, size_t values) {
        auto storage = msh.backend_storage();
        size_t doubles = 0;
        for (size_t key : keys)
            doubles += integrate(msh, storage->surfaces[key], 2 * deg).size() * (3 + values);
        return doubles;
    }

    // disk::make_rhs on the quadrature data q
    template<typename Quadrature, typename Function>
    static Eigen::Matrix<T, Eigen::Dynamic, 1> rhs(const Quadrature& q, size_t basis_size, const Function& rhs_fun) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> ret = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(basis_size);
        for (size_t k = 0; k < q.points.size(); k++) {
            const auto qp_f = disk::priv::inner_product(q.weights[k], rhs_fun(q.points[k]));
            ret += disk::priv::outer_product(q.phi[k], qp_f);
        }
        return ret;
    }

    // Assembly indices and Dirichlet faces of the cells elements with boundary faces (as scatter_e_bc_data,
    // scatter_a_bc_data): cell unknowns at cell_offset + n_cbs cell_index(cell), face unknowns at face_offset +
    // n_fbs compress(face)
    template<typename Cells, typename Bnd, typename CellIndex, typename Compress, typename MakeFaceBasis>
    static void prepare_bc_cells(const Mesh& msh, Cells& cells, const std::vector<size_t>& elements, Bnd& bnd,
                                 const CellIndex& cell_index, const Compress& compress, size_t n_cbs, size_t n_fbs,
                                 size_t cell_offset, size_t face_offset, MakeFaceBasis make_face_basis, bool parallel) {
        auto storage = msh.backend_storage();
        cells.assign(elements.size(), {});
        #pragma omp parallel for schedule(dynamic, 8) if (parallel)
        for (long i = 0; i < (long)elements.size(); i++) {
            auto& cl = storage->surfaces[elements[i]];
            auto& c = cells[i];
            const size_t cell_LHS_offset = cell_index.at(elements[i]) * n_cbs + cell_offset;
            for (size_t l = 0; l < n_cbs; l++) c.asm_map.push_back(assembly_index(cell_LHS_offset + l, true));
            auto fcs = faces(msh, cl);
            for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
                auto fc = fcs[face_i];
                auto fc_id = msh.lookup(fc);
                const size_t face_LHS_offset = face_offset + compress.at(fc_id) * n_fbs;
                const bool dirichlet = bnd.is_dirichlet_face(fc_id);
                for (size_t l = 0; l < n_fbs; l++) c.asm_map.push_back(assembly_index(face_LHS_offset + l, !dirichlet));
                if (!dirichlet) continue;
                typename std::decay_t<decltype(c.faces)>::value_type f;
                f.face_i = face_i;
                f.fc_id = fc_id;
                auto fb = make_face_basis(fc);
                for (auto& qp : integrate(msh, fc, 2 * fb.degree())) {
                    f.quad.points.push_back(qp.point());
                    f.quad.weights.push_back(qp.weight());
                    f.quad.phi.push_back(fb.eval_functions(qp.point()));
                }
                f.llt.compute(make_mass_matrix(msh, fc, fb));
                c.faces.push_back(std::move(f));
            }
        }
    }

    // Dirichlet data of a cell: projections of the data on its Dirichlet faces (make_rhs and the mass, as
    // scatter_e_bc_data), zero elsewhere
    template<typename Cell, typename Bnd>
    static Eigen::Matrix<T, Eigen::Dynamic, 1> bc_dirichlet_data(const Cell& c, Bnd& bnd, size_t n_cbs, size_t n_fbs) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> data = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(c.asm_map.size());
        for (auto& f : c.faces) {
            auto dirichlet_fun = bnd.dirichlet_boundary_func(f.fc_id);
            Eigen::Matrix<T, Eigen::Dynamic, 1> rhs = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n_fbs);
            for (size_t q = 0; q < f.quad.points.size(); q++) {
                const auto qp_f = disk::priv::inner_product(f.quad.weights[q], dirichlet_fun(f.quad.points[q]));
                rhs += disk::priv::outer_product(f.quad.phi[q], qp_f);
            }
            data.block(n_cbs + f.face_i * n_fbs, 0, n_fbs, 1) = f.llt.solve(rhs);
        }
        return data;
    }

    // RHS -= lhs (Dirichlet columns) data, in the order of scatter_e_bc_data
    static void scatter_bc_cached(Eigen::Matrix<T, Eigen::Dynamic, 1>& RHS, const std::vector<assembly_index>& asm_map,
                                  const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& data) {
        for (size_t i = 0; i < size_t(lhs.rows()); i++) {
            if (!asm_map[i].assemble()) continue;
            for (size_t j = 0; j < size_t(lhs.cols()); j++)
                if (!asm_map[j].assemble()) RHS(asm_map[i].vidx()) -= lhs(i,j) * data(j);
        }
    }
};

#endif /* assembly_tools_hpp */
