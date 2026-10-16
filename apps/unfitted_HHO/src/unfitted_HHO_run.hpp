/*
 * Unfitted HHO methods stabilized by polynomial extension for elliptic
 * interface problems (E. Burman, A. Ern, R. Mottier): convergence test on a
 * mesh of any type of diskpp (2D or 3D), see run_<mesh>.cpp.
 */

#pragma once

#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#ifdef HAVE_INTEL_MKL
#include <Eigen/PardisoSupport>
#include <mkl_service.h>
#endif

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/output/cuthho_output.hpp"
#include "diskpp/methods/hho"

#include "common/preprocessor.hpp"
#include "common/test_cases.hpp"
#include "common/operators.hpp"
#include "common/cuthho_assemblers.hpp"
#include "common/postprocessor.hpp"
#include "common/static_condensation.hpp"
#include "common/bvc.hpp"
#include "common/bvc_extension.hpp"

/* Wall clock timer (the timecounter of diskpp can only be included in one
 * translation unit) */
class wall_timer {
    std::chrono::steady_clock::time_point m_start;
    double m_elapsed = 0;
public:
    void tic() { m_start = std::chrono::steady_clock::now(); }
    double toc() {
        m_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count();
        return m_elapsed;
    }
    double elapsed() const { return m_elapsed; }
};

inline std::ostream& operator<<(std::ostream& os, const wall_timer& t) { return os << t.elapsed(); }

/* solver: 0 SimplicialLDLT (ProtoN), 1 Jacobi-preconditioned CG, 2 Pardiso LLT
 * (supernodal, multithreaded, if Intel MKL is available; SimplicialLDLT
 * otherwise). The Intel threading layer of MKL fails with the OpenMP runtime
 * of LLVM: the GNU threading layer is linked (see CMakeLists.txt). Should the
 * factorization fail anyway, it is retried with one thread, and
 * SimplicialLDLT is the last resort. */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
solve_linear_system(const Eigen::SparseMatrix<T>& A, const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, size_t solver_type)
{
#ifdef HAVE_INTEL_MKL
    if (solver_type == 2) {
        for (int attempt = 0; attempt < 2; attempt++) {
            int previous = (attempt == 1) ? mkl_set_num_threads_local(1) : 0;
            Eigen::PardisoLLT<Eigen::SparseMatrix<T>> solver;
            solver.compute(A);
            bool ok = (solver.info() == Eigen::Success);
            Eigen::Matrix<T, Eigen::Dynamic, 1> x;
            if (ok)
                x = solver.solve(b);
            if (attempt == 1)
                mkl_set_num_threads_local(previous);
            if (ok)
                return x;
            std::cout << "         Pardiso failed" << (attempt == 0 ? ", retried with one thread" : ", SimplicialLDLT used")
                      << std::endl;
        }
        solver_type = 0;
    }
#endif
    if (solver_type == 1) {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower|Eigen::Upper> cg;
        cg.setTolerance(1e-14);
        cg.compute(A);
        return cg.solve(b);
    }

    /* Direct solver of ProtoN */
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("Factorization of the stiffness matrix failed");
    return solver.solve(b);
}

/* Euclidean condition numbers of a symmetric matrix A and of D^-1/2 A D^-1/2,
 * D = diag(A) (dense, coarse meshes only). The second one does not depend on
 * the scaling of the unknowns: it governs the accuracy of a Cholesky
 * factorization and the convergence of the Jacobi-preconditioned CG. */
template<typename T>
std::pair<T, T>
condition_numbers(const Eigen::SparseMatrix<T>& A)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    auto cond = [](const matrix_type& M) {
        Eigen::SelfAdjointEigenSolver<matrix_type> es(M, Eigen::EigenvaluesOnly);
        auto ev = es.eigenvalues().cwiseAbs();
        return ev.maxCoeff() / ev.minCoeff();
    };
    matrix_type Ad(A);
    Eigen::Matrix<T, Eigen::Dynamic, 1> d = Ad.diagonal().cwiseAbs().cwiseSqrt().cwiseInverse();
    matrix_type As = d.asDiagonal() * Ad * d.asDiagonal();
    return std::make_pair(cond(Ad), cond(As));
}

/* Convergence test for the elliptic interface problem, stabilization by
 * polynomial extension (Figs. 6, 7, 9, 10, 11, 12) or, if agglomeration is
 * true, by cell agglomeration as in [12] (dashed lines of Figs. 6, 7, 12). */
template<typename Mesh>
void Elliptic_Interface(const simulation_options& opts, bool agglomeration)
{
    using T = typename Mesh::coordinate_type;
    using pre = preprocessor<Mesh>;
    using cut_mesh_type = typename pre::cut_mesh_type;
    using post = postprocessor<cut_mesh_type>;

    wall_timer tc, tck, tcl;
    tc.tic();

    auto level_set_function = pre::make_level_set_function(opts);

    params<T> parms(1.0, std::pow(10.0, T(opts.contrast)));

    std::string error_file_txt = "solution_error_file.txt";
    std::ofstream error_file(error_file_txt);
    error_file << std::setprecision(14);
    post::write_conv_sol(error_file_txt);

    size_t k_min = std::min(opts.k_min, opts.degree);
    size_t l_min = std::min(opts.l_min, opts.l_divs);

    for (size_t k = k_min; k <= opts.degree; k++) {

        tck.tic();
        std::cout << std::endl << bold << red << "   Polynomial degree k : " << k << reset << std::endl;
        error_file << std::endl << "Polynomial degree k : " << k << std::endl;

        mesh_parameters<cut_mesh_type> prm(k);
        prm.eta_N = (opts.eta_N >= 0) ? opts.eta_N : -opts.eta_N*(k+1)*(k+1);
        prm.basis.kind = (opts.basis == 0) ? disk::cut_basis_kind::cell_monomial
                       : (opts.basis == 1) ? disk::cut_basis_kind::subcell_monomial
                                           : disk::cut_basis_kind::orthonormal;
        if (agglomeration)
            prm.eta_gamma = prm.hdi.cell_degree(); // as in ProtoN, s^Gamma is scaled like s°

        typename post::errors previous{0.0, 0.0, 0.0, 0.0, 0.0};

        for (size_t l = l_min; l <= opts.l_divs; l++) {

            tcl.tic();
            std::cout << bold << cyan << "      Space refinment level -l : " << l << reset << std::endl;
            error_file << "Space refinment level -l : " << l << std::endl;

            /* timings: mesh, assembly (of which the cut cells), solve, errors */
            wall_timer tph;
            double t_mesh = 0, t_asm = 0, t_asm_cut = 0, t_solve = 0, t_err = 0;
            tph.tic();

            // MESH
            auto msh = pre::msh_generation(level_set_function, l, opts, geometry_options(opts, k), !agglomeration);
            prm.basis.cache = std::make_shared<disk::cut_basis_cache<T, Mesh::dimension>>(); // frames of this mesh
            prm.basis.cache->degree = prm.hdi.cell_degree();
            prm.basis.cache->face_degree = prm.hdi.face_degree();
            if (opts.cache_uncut)
                prm.uncut_cache = std::make_shared<typename decltype(prm.uncut_cache)::element_type>();
            disk::cut_agglomeration agglo;
            if (agglomeration)
                agglo = disk::make_cell_agglomeration(msh, opts.debug != 0, opts.pairing != 0);
            if (opts.debug)
                pre::msh_debug(msh, *level_set_function, l, agglo);

            t_mesh = tph.toc();

            auto test_case = select_test_case<cut_mesh_type>(opts.test_case, level_set_function, parms);
            if (l == l_min && k == k_min)
                std::cout << "         Test case: " << test_case.name << std::endl;
            std::cout << "         Mesh: " << msh.cells_size() << " cells, " << msh.faces_size() << " faces" << std::endl;

            // ASSEMBLY, SOLVE AND ERRORS
            auto sparsity_file = "sparsity_k" + std::to_string(k) + "_l" + std::to_string(l) + ".txt";
            auto solve = [&](auto& assembler) {
                tph.tic();
                assembler.finalize();
                if (opts.sparsity) {
                    assembler.write_block_pattern(sparsity_file);
                    std::cout << "         Block sparsity pattern written to " << sparsity_file << std::endl;
                }
                if (opts.conditioning) {
                    auto cond = condition_numbers(assembler.LHS);
                    std::cout << bold << yellow << "         Condition number: " << cond.first
                              << " (diagonally scaled: " << cond.second << ")" << reset << std::endl;
                    error_file << "condition number: " << cond.first << std::endl;
                    error_file << "condition number (diagonally scaled): " << cond.second << std::endl;
                }
                Eigen::Matrix<T, Eigen::Dynamic, 1> x;
                if (opts.condensation == 2) {
                    /* algebraic elimination of the cell unknowns, solve on the faces */
                    wall_timer tsc;
                    tsc.tic();
                    algebraic_condensation<T> sc(assembler.LHS, assembler.RHS, assembler.num_cell_dofs());
                    auto t_condense = tsc.toc();
                    std::cout << "         Algebraic static condensation (" << t_condense << " s): " << sc.num_clusters()
                              << " clusters of cell unknowns (largest: " << sc.max_cluster << "), " << sc.S.rows()
                              << " face unknowns" << std::endl;
                    if (opts.conditioning) {
                        auto cond = condition_numbers(sc.S);
                        std::cout << bold << yellow << "         Condition number of the condensed system: " << cond.first
                                  << " (diagonally scaled: " << cond.second << ")" << reset << std::endl;
                        error_file << "condition number (condensed): " << cond.first << std::endl;
                    }
                    x = sc.recover(solve_linear_system(sc.S, sc.g, opts.iterative_solver));
                }
                else
                    x = solve_linear_system(assembler.LHS, assembler.RHS, opts.iterative_solver);
                t_solve = tph.toc();
                return x;
            };
            /* solve of the system condensed patch by patch (-c 1) */
            auto solve_local = [&](local_condensation<T>& lc) {
                tph.tic();
                lc.finalize();
                std::cout << "         Local static condensation: " << lc.num_patches() << " patches (largest: "
                          << lc.max_patch << " cell unknowns), " << lc.S.rows() << " face unknowns" << std::endl;
                if (opts.conditioning) {
                    auto cond = condition_numbers(lc.S);
                    std::cout << bold << yellow << "         Condition number of the condensed system: " << cond.first
                              << " (diagonally scaled: " << cond.second << ")" << reset << std::endl;
                    error_file << "condition number (condensed): " << cond.first << std::endl;
                }
                auto x = lc.recover(solve_linear_system(lc.S, lc.g, opts.iterative_solver));
                t_solve = tph.toc();
                return x;
            };
            wall_timer tel;

            typename post::errors err;
            if (!agglomeration) {
                auto partition = disk::make_subcell_partition(msh);
                const auto& P_OK = partition.first;
                const auto& P_KO = partition.second;

                unfitted_elliptic_interface_assembler<cut_mesh_type> assembler(msh, prm.hdi, test_case.bcs_fun, prm.basis);
                assembler.record_pattern(opts.sparsity != 0);
                Eigen::Matrix<T, Eigen::Dynamic, 1> x_dof;
                if (opts.condensation == 1) {
                    /* patches: connected components of the relation S ~ N(S) */
                    std::vector<size_t> parent(msh.cells_size());
                    std::iota(parent.begin(), parent.end(), 0);
                    std::function<size_t(size_t)> find = [&](size_t c) {
                        return parent[c] == c ? c : parent[c] = find(parent[c]);
                    };
                    for (auto& cl : msh) {
                        const auto& ccd = msh.cut_cell_data[offset(msh, cl)];
                        if (ccd.stabilizing_cell.has_value())
                            parent[find(offset(msh, cl))] = find(ccd.stabilizing_cell.value());
                    }
                    std::map<size_t, std::vector<std::pair<const disk::cut_subcell*, bool>>> patches;
                    for (auto& P : P_OK)
                        patches[find(P.cell)].push_back(std::make_pair(&P, true));
                    for (auto& P : P_KO)
                        patches[find(P.cell)].push_back(std::make_pair(&P, false));

                    local_condensation<T> lc(assembler.num_cell_dofs(), assembler.LHS.rows());
                    for (auto& pt : patches) {
                        tel.tic();
                        bool cut = false;
                        std::vector<typename local_condensation<T>::element> elems;
                        for (auto& Pk : pt.second) {
                            const auto& P = *Pk.first;
                            auto contrib = Pk.second ? make_contrib_POK(msh, P, test_case, prm)
                                                     : make_contrib_PKO(msh, P, test_case, prm);
                            typename local_condensation<T>::element e;
                            e.asm_map = assembler.init_asm_map(msh, P, e.dirichlet);
                            e.lhs = contrib.first;
                            e.rhs = contrib.second;
                            elems.push_back(std::move(e));
                            cut = cut || disk::is_cut(msh, msh[P.cell]) || !P.dependent_cells.empty();
                        }
                        lc.condense_patch(elems);
                        auto t = tel.toc();
                        t_asm += t;
                        if (cut)
                            t_asm_cut += t;
                    }
                    x_dof = solve_local(lc);
                }
                else {
                    for (auto& P : P_OK) {
                        tel.tic();
                        auto contrib = make_contrib_POK(msh, P, test_case, prm);
                        assembler.assemble(msh, P, contrib.first, contrib.second);
                        auto t = tel.toc();
                        t_asm += t;
                        if (disk::is_cut(msh, msh[P.cell]) || !P.dependent_cells.empty())
                            t_asm_cut += t;
                    }
                    for (auto& P : P_KO) {
                        tel.tic();
                        auto contrib = make_contrib_PKO(msh, P, test_case, prm);
                        assembler.assemble(msh, P, contrib.first, contrib.second);
                        auto t = tel.toc();
                        t_asm += t;
                        t_asm_cut += t;
                    }
                    x_dof = solve(assembler);
                }
                tph.tic();
                err = post::compute_errors(msh, P_OK, prm, assembler, x_dof, test_case);
                t_err = tph.toc();
            }
            else {
                auto macros = make_macro_cells(msh, agglo);
                unfitted_agglomeration_assembler<cut_mesh_type> assembler(msh, macros, prm.hdi, test_case.bcs_fun, prm.basis);
                assembler.record_pattern(opts.sparsity != 0);
                Eigen::Matrix<T, Eigen::Dynamic, 1> x_dof;
                /* -c 1: one patch per agglomerated cell (both sides if cut) */
                local_condensation<T> lc(assembler.num_cell_dofs(), assembler.LHS.rows());
                for (size_t q = 0; q < macros.size(); q++) {
                    tel.tic();
                    std::vector<typename local_condensation<T>::element> elems;
                    for (auto side : macro_sides(macros[q])) {
                        auto contrib = make_contrib_agglo(msh, macros[q], side, test_case, prm);
                        if (opts.condensation == 1) {
                            typename local_condensation<T>::element e;
                            e.asm_map = assembler.init_asm_map(msh, q, macros[q], side, e.dirichlet);
                            e.lhs = contrib.first;
                            e.rhs = contrib.second;
                            elems.push_back(std::move(e));
                        }
                        else
                            assembler.assemble(msh, q, macros[q], side, contrib.first, contrib.second);
                    }
                    if (opts.condensation == 1)
                        lc.condense_patch(elems);
                    auto t = tel.toc();
                    t_asm += t;
                    if (macros[q].cut)
                        t_asm_cut += t;
                }
                x_dof = (opts.condensation == 1) ? solve_local(lc) : solve(assembler);
                tph.tic();
                err = post::compute_errors_agglo(msh, macros, prm, assembler, x_dof, test_case);
                t_err = tph.toc();
            }
            if (opts.refined_mesh) {
                /* locally refined mesh: size of the background mesh */
                err.h = 0;
                for (auto& cl : msh)
                    err.h = std::max(err.h, diameter(msh, cl));
            }
            post::write_errors(err, previous, error_file);
            std::cout << bold << yellow << "         H1-Error: " << err.H1 << "   L2-Error: " << err.L2 << reset << std::endl;
            if (previous.h > 0)
                std::cout << bold << yellow << "         order H1: " << std::log(previous.H1/err.H1)/std::log(previous.h/err.h)
                          << "   order L2: " << std::log(previous.L2/err.L2)/std::log(previous.h/err.h) << reset << std::endl;
            previous = err;

            tcl.toc();
            std::cout << "         Times (s): mesh " << t_mesh << ", assembly " << t_asm << " (cut cells " << t_asm_cut
                      << "), solve " << t_solve << ", errors " << t_err << std::endl;
            std::cout << bold << yellow << "         Run l = " << l << " completed: " << tcl << " seconds" << reset << std::endl;
        }

        error_file << std::endl << std::endl;
        tck.toc();
        std::cout << bold << cyan << "      Run k = " << k << " completed: " << tck << " seconds" << reset << std::endl;
    }

    tc.toc();
    std::cout << std::endl << bold << red << "   Run completed: " << tc << " seconds" << reset << std::endl << std::endl;
}

/* The convergence tests on each type of mesh (run_<mesh>.cpp) */
void run_quadrangles(const simulation_options& opts, bool agglomeration);
void run_triangles(const simulation_options& opts, bool agglomeration);
void run_polygons(const simulation_options& opts, bool agglomeration);
void run_hexahedra(const simulation_options& opts, bool agglomeration);
void run_tetrahedra(const simulation_options& opts, bool agglomeration);
void run_polyhedra(const simulation_options& opts, bool agglomeration);

/* The convergence test on the meshes of type Mesh or, if the interface does
 * not cut once some of their cells and -x > 0, on their local refinements
 * (polytopal meshes, run by run_polygons and run_polyhedra); -p2, -p3:
 * unfitted HHO with boundary value correction, stabilized by cell
 * agglomeration (common/bvc.hpp) or by polynomial extension
 * (common/bvc_extension.hpp) */
template<typename Mesh>
void run_mesh(const simulation_options& opts, bool agglomeration)
{
    if (!opts.refined_mesh && opts.refinement > 0 && preprocessor<Mesh>::needs_refinement(opts)) {
        std::cout << "        The interface does not cut once some cells: local refinement of these cells"
                  << std::endl;
        auto ropts = opts;
        ropts.refined_mesh = true;
        if constexpr (Mesh::dimension == 2)
            run_polygons(ropts, agglomeration);
        else
            run_polyhedra(ropts, agglomeration);
        return;
    }
    if (opts.prototype == 2 || opts.prototype == 3)
        bvc::run<Mesh, wall_timer>(opts);
    else
        Elliptic_Interface<Mesh>(opts, agglomeration);
}
