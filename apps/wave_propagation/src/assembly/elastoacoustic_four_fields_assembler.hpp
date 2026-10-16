//
//  elastoacoustic_four_fields_hpp
//  acoustics
//
//  Created by Omar Durán on 9/7/20.
//  Contributor: Romain Mottier


#pragma once
#ifndef elastoacoustic_four_fields_assembler_hpp
#define elastoacoustic_four_fields_assembler_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "assembly_index.hpp"
#include "../physics/acoustic_material_data.hpp"
#include "../physics/elastic_material_data.hpp"
#include "elastoacoustic_assembler_base.hpp"
#include "assembly_tools.hpp"
#include <map>
#include <memory>

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class elastoacoustic_four_fields_assembler : public elastoacoustic_assembler_base<Mesh> {
    
    typedef elastoacoustic_assembler_base<Mesh>      base;
    typedef typename base::e_boundary_type           e_boundary_type;
    typedef typename base::a_boundary_type           a_boundary_type;
    using T = typename Mesh::coordinate_type;
    typedef assembly_tools<Mesh, T>                  tools;
    using base::m_e_compress_indexes;
    using base::m_e_expand_indexes;
    using base::m_a_compress_indexes;
    using base::m_a_expand_indexes;
    using base::m_hho_di;
    using base::m_e_bnd;
    using base::m_a_bnd;
    using base::m_triplets;
    using base::m_c_triplets;
    using base::m_mass_triplets;
    using base::m_e_material;
    using base::m_a_material;
    using base::m_e_cell_index;
    using base::m_a_cell_index;
    using base::m_e_elements_with_bc_eges;
    using base::m_a_elements_with_bc_eges;
    using base::m_interface_cell_indexes;
    using base::m_n_elastic_cells;
    using base::m_n_acoustic_cells;
    using base::n_e_edges;
    using base::n_a_edges;
    using base::m_n_edges;
    using base::m_n_essential_edges;
    using base::m_n_e_cbs;
    using base::m_n_a_cbs;
    using base::m_n_elastic_cell_dof;
    using base::m_n_acoustic_cell_dof;
    using base::m_n_elastic_face_dof;
    using base::m_n_acoustic_face_dof;
    using base::m_hho_stabilization_Q;

    std::vector< Triplet<T> >           m_triplets_stab;
    bool        m_scaled_stabilization_Q;
    T           m_e_eta = 1.50;     // stabilization weights eta (elastic, acoustic): MEKG values by default (MEKG uses
    T           m_a_eta = 0.80;     //   0.95 / 0.4 in mixed order, not used here); 1 in Burman-Duran-Ern(-Steins)
    int         m_stab_size = 0;    // length of the 1/h scaling, see set_stabilization_size
    T           m_graded_href = 0;  // > 0: O(1) weight times max(1, m_graded_href / h_T) (set_stabilization_graded)
    T           m_hybrid_hc = 0;    // > 0: weight in 1/h on the cells smaller than m_hybrid_hc, O(1) on the others
    bool        m_parallel_assembly = false;    // see set_parallel_assembly

    // assembly_tools::ordered_loop with m_parallel_assembly: the serial loop, or the local matrices and vectors computed by
    // the OpenMP threads and scattered in the same order (the same global matrices and vectors, bit for bit)
    template<typename Local, typename Scatter>
    void ordered_loop(size_t n, Local&& local, Scatter&& scatter) {
        tools::ordered_loop(m_parallel_assembly, n, std::forward<Local>(local), std::forward<Scatter>(scatter));
    }


public:

    using base::LHS;
    using base::RHS;
    using base::MASS;
    using base::COUPLING;
    using base::assemble_coupling_terms;
    using base::build_cells_maps;
    using base::classify_cells;
    using base::mark_fine_region;
    using base::e_interface_operator;
    using base::finalize_coupling;
    using base::finalize_mass;
    using base::gather_a_dof_data;
    using base::gather_e_dof_data;
    using base::get_a_bc_conditions;
    using base::get_a_cell_basis_data;
    using base::get_a_compress;
    using base::get_a_expand;
    using base::get_a_face_basis_data;
    using base::get_a_face_dof;
    using base::get_a_material_data;
    using base::get_a_n_cells_dof;
    using base::get_acoustic_cells;
    using base::get_acoustic_faces;
    using base::get_e_bc_conditions;
    using base::get_e_cell_basis_data;
    using base::get_e_compress;
    using base::get_e_expand;
    using base::get_e_face_basis_data;
    using base::get_e_face_dof;
    using base::get_e_material_data;
    using base::get_e_n_cells_dof;
    using base::get_elastic_cells;
    using base::get_elastic_faces;
    using base::get_interfaces;
    using base::get_n_face_dof;
    using base::project_over_faces;
    using base::scatter_a_bc_data;
    using base::scatter_a_cell_dof_data;
    using base::scatter_a_data;
    using base::scatter_a_face_dof_data;
    using base::scatter_a_mass_data;
    using base::scatter_a_rhs_data;
    using base::scatter_e_bc_data;
    using base::scatter_e_cell_dof_data;
    using base::scatter_e_data;
    using base::scatter_e_face_dof_data;
    using base::scatter_e_mass_data;
    using base::scatter_e_rhs_data;
    using base::scatter_ea_interface_data;
    using base::set_hdg_stabilization;
    using base::set_hho_stabilization;
    using base::set_least_squares_stabilization;
    using base::set_interface_cell_indexes;

    SparseMatrix<T>         LHS_STAB;
    SparseMatrix<T>         P_cell;
    SparseMatrix<T>         IminusP_cell;
    SparseMatrix<T>         Pfacefine;
    SparseMatrix<T>         Pfacecoarse;
    SparseMatrix<T>         Pfine;
    SparseMatrix<T>         Pcoarse;
    std::vector<bool> cell_is_fine_silo;

    size_t nnz_Sff_fine;

    // Cell unknowns: symmetric stress tensor + velocity (elastic), velocity in grad P^{k+1} + pressure (acoustic)
    static size_t e_cell_basis_size(const disk::hho_degree_info& hho_di) {
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        return n_ten_cbs + n_vec_cbs;
    }
    
    static size_t a_cell_basis_size(const disk::hho_degree_info& hho_di) {
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        return n_vel_scal_cbs + n_scal_cbs;
    }

    elastoacoustic_four_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const e_boundary_type& e_bnd, const a_boundary_type& a_bnd, std::map<size_t,elastic_material_data<T>> & e_material, std::map<size_t,acoustic_material_data<T>> & a_material) : base(msh, hho_di, e_bnd, a_bnd, e_material, a_material, e_cell_basis_size(hho_di), a_cell_basis_size(hho_di)), m_scaled_stabilization_Q(false) {
        
        size_t system_size = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_n_acoustic_face_dof;
        LHS_STAB = SparseMatrix<T>( system_size, system_size ); //OPTIONAL
        P_cell = SparseMatrix<T>(m_n_elastic_cell_dof + m_n_acoustic_cell_dof, m_n_elastic_cell_dof + m_n_acoustic_cell_dof);
        IminusP_cell = SparseMatrix<T>(m_n_elastic_cell_dof + m_n_acoustic_cell_dof, m_n_elastic_cell_dof + m_n_acoustic_cell_dof);
        Pfacefine   = SparseMatrix<T>(m_n_elastic_face_dof + m_n_acoustic_face_dof, m_n_elastic_face_dof + m_n_acoustic_face_dof);
        Pfacecoarse = SparseMatrix<T>(m_n_elastic_face_dof + m_n_acoustic_face_dof, m_n_elastic_face_dof + m_n_acoustic_face_dof);
        Pfine   = SparseMatrix<T>(system_size, system_size);
        Pcoarse = SparseMatrix<T>(system_size, system_size);
    }
    
    // Assembly of global mass matrix; local elastic and acoustic cell mass matrix computation; 
    // Scattering of elastic and acoustic mass contributions into global triplets; Finalization of mass assembly.
    void assemble_mass(const Mesh& msh, bool add_vector_mass_Q = true) {
        
        auto storage = msh.backend_storage();
        MASS.setZero();
        
        if (!m_parallel_assembly) {
            // elastic block
            for (auto e_chunk : m_e_material) {
                size_t e_cell_ind = m_e_cell_index[e_chunk.first];
                auto& cell = storage->surfaces[e_chunk.first];
                elastic_material_data<T> material = e_chunk.second;
                Matrix<T, Dynamic, Dynamic> mass_matrix = e_mass_operator(material, msh, cell, add_vector_mass_Q);
                scatter_e_mass_data(e_cell_ind,msh, cell, mass_matrix);
            }

            // acoustic block
            for (auto a_chunk : m_a_material) {
                size_t a_cell_ind = m_a_cell_index[a_chunk.first];
                auto& cell = storage->surfaces[a_chunk.first];
                acoustic_material_data<T> material = a_chunk.second;
                Matrix<T, Dynamic, Dynamic> mass_matrix = a_mass_operator(material, msh, cell);
                scatter_a_mass_data(a_cell_ind, msh, cell, mass_matrix);
            }
            finalize_mass();
            return;
        }

        if (m_class_reuse) {
            // mass matrices once per geometric class
            auto& cl = cell_classes(msh);
            auto e_ops = tools::class_operators(cl.e_rep, [&](size_t r) {
                elastic_material_data<T> material = m_e_material.at(cl.e_ids[r]);
                return Matrix<T, Dynamic, Dynamic>(e_mass_operator(material, msh, storage->surfaces[cl.e_ids[r]], add_vector_mass_Q));
            });
            for (size_t i = 0; i < cl.e_ids.size(); i++)
                scatter_e_mass_data(m_e_cell_index[cl.e_ids[i]], msh, storage->surfaces[cl.e_ids[i]], e_ops[cl.e_class_of[i]]);
            auto a_ops = tools::class_operators(cl.a_rep, [&](size_t r) {
                acoustic_material_data<T> material = m_a_material.at(cl.a_ids[r]);
                return Matrix<T, Dynamic, Dynamic>(a_mass_operator(material, msh, storage->surfaces[cl.a_ids[r]]));
            });
            for (size_t i = 0; i < cl.a_ids.size(); i++)
                scatter_a_mass_data(m_a_cell_index[cl.a_ids[i]], msh, storage->surfaces[cl.a_ids[i]], a_ops[cl.a_class_of[i]]);
            finalize_mass_parallel();
            return;
        }

        // elastic block
        std::vector<std::pair<size_t, elastic_material_data<T>>> e_cells(m_e_material.begin(), m_e_material.end());
        ordered_loop(e_cells.size(), [&](size_t i) {
            elastic_material_data<T> material = e_cells[i].second;
            return Matrix<T, Dynamic, Dynamic>(e_mass_operator(material, msh, storage->surfaces[e_cells[i].first], add_vector_mass_Q));
        }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
            scatter_e_mass_data(m_e_cell_index[e_cells[i].first], msh, storage->surfaces[e_cells[i].first], mass_matrix);
        });

        // acoustic block
        std::vector<std::pair<size_t, acoustic_material_data<T>>> a_cells(m_a_material.begin(), m_a_material.end());
        ordered_loop(a_cells.size(), [&](size_t i) {
            acoustic_material_data<T> material = a_cells[i].second;
            return Matrix<T, Dynamic, Dynamic>(a_mass_operator(material, msh, storage->surfaces[a_cells[i].first]));
        }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
            scatter_a_mass_data(m_a_cell_index[a_cells[i].first], msh, storage->surfaces[a_cells[i].first], mass_matrix);
        });
        
        finalize_mass();
    }
    
    // Application of boundary conditions; local elastic and acoustic cell operator computation for cells with boundary faces; 
    // Scattering of Dirichlet contributions into the global RHS;
    void apply_bc(const Mesh& msh, bool explicit_scheme){
        
        #ifdef HAVE_INTEL_TBB2
                size_t n_cells = m_elements_with_bc_eges.size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh] (size_t & i){
                        size_t cell_ind = m_elements_with_bc_eges[i];
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
                        scatter_bc_data(msh, cell, laplacian_operator_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
            if (m_rc) {
                if (m_rc->bc_ops_scheme != int(explicit_scheme) && m_class_reuse) {
                    // local operators of the classes of the cells
                    auto& cl = cell_classes(msh);
                    auto e_ops = e_class_operators(msh, explicit_scheme);
                    auto a_ops = a_class_operators(msh, explicit_scheme);
                    m_rc->e_bc_ops.clear();
                    m_rc->a_bc_ops.clear();
                    for (size_t cell_ind : m_e_elements_with_bc_eges) m_rc->e_bc_ops.push_back(e_ops[cl.e_class_of[cl.e_pos.at(cell_ind)]]);
                    for (size_t cell_ind : m_a_elements_with_bc_eges) m_rc->a_bc_ops.push_back(a_ops[cl.a_class_of[cl.a_pos.at(cell_ind)]]);
                    m_rc->bc_ops_scheme = int(explicit_scheme);
                    prepare_bc_cells(msh);
                }
                if (m_rc->bc_ops_scheme != int(explicit_scheme)) {
                    m_rc->e_bc_ops.assign(m_e_elements_with_bc_eges.size(), Matrix<T, Dynamic, Dynamic>());
                    m_rc->a_bc_ops.assign(m_a_elements_with_bc_eges.size(), Matrix<T, Dynamic, Dynamic>());
                    ordered_loop(m_e_elements_with_bc_eges.size(), [&](size_t i) {
                        const size_t cell_ind = m_e_elements_with_bc_eges[i];
                        elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                        return Matrix<T, Dynamic, Dynamic>(e_mixed_operator(e_mat, msh, storage->surfaces[cell_ind], explicit_scheme));
                    }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& op) { m_rc->e_bc_ops[i] = op; });
                    ordered_loop(m_a_elements_with_bc_eges.size(), [&](size_t i) {
                        const size_t cell_ind = m_a_elements_with_bc_eges[i];
                        acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                        return Matrix<T, Dynamic, Dynamic>(a_mixed_operator(a_mat, msh, storage->surfaces[cell_ind], explicit_scheme));
                    }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& op) { m_rc->a_bc_ops[i] = op; });
                    m_rc->bc_ops_scheme = int(explicit_scheme);
                    prepare_bc_cells(msh);
                }
                // Dirichlet data of each cell (in parallel), then its contributions in the order of the serial loop
                const size_t n_e_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
                const size_t n_a_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
                ordered_loop(m_rc->e_bc_cells.size(), [&](size_t i) {
                    return tools::bc_dirichlet_data(m_rc->e_bc_cells[i], m_e_bnd, get_e_cell_basis_data(), n_e_fbs);
                }, [&](size_t i, const Matrix<T, Dynamic, 1>& data) { tools::scatter_bc_cached(RHS, m_rc->e_bc_cells[i].asm_map, m_rc->e_bc_ops[i], data); });
                ordered_loop(m_rc->a_bc_cells.size(), [&](size_t i) {
                    return tools::bc_dirichlet_data(m_rc->a_bc_cells[i], m_a_bnd, get_a_cell_basis_data(), n_a_fbs);
                }, [&](size_t i, const Matrix<T, Dynamic, 1>& data) { tools::scatter_bc_cached(RHS, m_rc->a_bc_cells[i].asm_map, m_rc->a_bc_ops[i], data); });
                return;
            }
            if (!m_parallel_assembly) {
                for (auto& cell_ind : m_e_elements_with_bc_eges) {
                    auto& cell = storage->surfaces[cell_ind];
                    size_t e_cell_ind = m_e_cell_index[cell_ind];
                    elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_mat,msh,cell,explicit_scheme);
                    scatter_e_bc_data(e_cell_ind, msh, cell, mixed_operator_loc);
                }

                for (auto& cell_ind : m_a_elements_with_bc_eges) {
                    auto& cell = storage->surfaces[cell_ind];
                    size_t a_cell_ind = m_a_cell_index[cell_ind];
                    acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_mat, msh, cell, explicit_scheme);
                    scatter_a_bc_data(a_cell_ind, msh, cell, mixed_operator_loc);
                }
                return;
            }
            ordered_loop(m_e_elements_with_bc_eges.size(), [&](size_t i) {
                const size_t cell_ind = m_e_elements_with_bc_eges[i];
                elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                return Matrix<T, Dynamic, Dynamic>(e_mixed_operator(e_mat, msh, storage->surfaces[cell_ind], explicit_scheme));
            }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mixed_operator_loc) {
                const size_t cell_ind = m_e_elements_with_bc_eges[i];
                scatter_e_bc_data(m_e_cell_index[cell_ind], msh, storage->surfaces[cell_ind], mixed_operator_loc);
            });

            ordered_loop(m_a_elements_with_bc_eges.size(), [&](size_t i) {
                const size_t cell_ind = m_a_elements_with_bc_eges[i];
                acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                return Matrix<T, Dynamic, Dynamic>(a_mixed_operator(a_mat, msh, storage->surfaces[cell_ind], explicit_scheme));
            }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mixed_operator_loc) {
                const size_t cell_ind = m_a_elements_with_bc_eges[i];
                scatter_a_bc_data(m_a_cell_index[cell_ind], msh, storage->surfaces[cell_ind], mixed_operator_loc);
            });
        
        #endif
        
    }

    // Assembly of global RHS vector; Local elastic and acoustic cell RHS computation; 
    // Scattering of local RHS contributions into the global vector; Application of boundary condition corrections.
    void assemble_rhs(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme) {
        assemble_rhs_impl(msh, e_rhs_fun, a_rhs_fun, explicit_scheme, nullptr, false);
    }

    // EXPERIMENTAL (not used by any validated driver): same as assemble_rhs,
    // but the acoustic cell loop is restricted to `active_a_cells` (indexed
    // by ORIGINAL mesh cell index, i.e. the same indexing as m_a_material's
    // keys / a_chunk.first -- true = assemble this cell's RHS, false = skip
    // it, leaving that cell's RHS entries at the zero set by RHS.setZero()
    // above). Intended for MLTS non-terminal bands, whose own_w_i is masked
    // to Pell_c[i]/Pell_f[i] downstream anyway -- assembling RHS on cells
    // outside that mask is pure waste for those bands (this is NOT valid to
    // use for the terminal band's own eval_F calls without separately
    // checking apply_bc's boundary-cell path, not audited here).
    void assemble_rhs_restricted(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme, const std::vector<char>& active_a_cells) {
        assemble_rhs_impl(msh, e_rhs_fun, a_rhs_fun, explicit_scheme, &active_a_cells, false);
    }

    // EXPERIMENTAL (not used by any validated driver): same as assemble_rhs,
    // but the acoustic cell loop runs under OpenMP. Safe without locks or a
    // reduction: scatter_a_rhs_data writes into a cell-exclusive dof range
    // (a_cell_ind*n_cbs offset, no two cells share a row), so concurrent
    // writes never race, same argument as pmv()'s row-major disjointness
    // (see that function's comment). UNLIKE pmv (thousands of calls, a few
    // hundred flops each -- measured slower under OpenMP due to per-call
    // thread-team overhead), each iteration here runs a_mixed_rhs (a
    // quadrature loop + basis evaluation + outer products), several orders
    // of magnitude more work per iteration -- the untested "batch into
    // fewer, larger parallel regions" idea from pmv's comment, applied at
    // the one call site in this file already coarse-grained enough to
    // plausibly amortize the thread-team overhead. Falls back to serial if
    // compiled without -fopenmp.
    void assemble_rhs_parallel(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme) {
        assemble_rhs_impl(msh, e_rhs_fun, a_rhs_fun, explicit_scheme, nullptr, true);
    }

    // EXPERIMENTAL: restriction + OpenMP combined (see the two methods above
    // for the individual rationale of each).
    void assemble_rhs_restricted_parallel(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme, const std::vector<char>& active_a_cells) {
        assemble_rhs_impl(msh, e_rhs_fun, a_rhs_fun, explicit_scheme, &active_a_cells, true);
    }

    // Shared body of the assemble_rhs variants: the acoustic loop is optionally restricted to active_a_cells
    // (indexed by mesh cell) and run under OpenMP.
    void assemble_rhs_impl(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme, const std::vector<char>* active_a_cells, bool parallel) {

        RHS.setZero();
        auto storage = msh.backend_storage();
        if (!m_parallel_assembly && !m_rc) {
            for (auto e_chunk : m_e_material) {
                size_t e_cell_ind = m_e_cell_index[e_chunk.first];
                auto& cell = storage->surfaces[e_chunk.first];
                Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
                scatter_e_rhs_data(e_cell_ind, msh, cell, f_loc);
            }

            std::vector<size_t> a_keys;
            a_keys.reserve(m_a_material.size());
            for (auto a_chunk : m_a_material) if (!active_a_cells || (*active_a_cells)[a_chunk.first]) a_keys.push_back(a_chunk.first);
            const int n = (int)a_keys.size();
            #pragma omp parallel for schedule(static) if(parallel)
            for (int k = 0; k < n; ++k) {
                size_t cell_ind = a_keys[(size_t)k];
                size_t a_cell_ind = m_a_cell_index[cell_ind];
                auto& cell = storage->surfaces[cell_ind];
                Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
                scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
            }
            apply_bc(msh, explicit_scheme);
            return;
        }
        std::vector<size_t> e_keys;
        e_keys.reserve(m_e_material.size());
        for (auto& e_chunk : m_e_material) e_keys.push_back(e_chunk.first);
        const bool cached = m_rc && !active_a_cells && !parallel && prepare_rhs_cache(msh, e_keys);
        ordered_loop(e_keys.size(), [&](size_t i) {
            if (cached) return e_mixed_rhs_cached(i, msh, storage->surfaces[e_keys[i]], e_rhs_fun);
            return Matrix<T, Dynamic, 1>(e_mixed_rhs(msh, storage->surfaces[e_keys[i]], e_rhs_fun));
        }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
            scatter_e_rhs_data(m_e_cell_index[e_keys[i]], msh, storage->surfaces[e_keys[i]], f_loc);
        });

        std::vector<size_t> a_keys;
        a_keys.reserve(m_a_material.size());
        for (auto a_chunk : m_a_material) if (!active_a_cells || (*active_a_cells)[a_chunk.first]) a_keys.push_back(a_chunk.first);
        const int n = (int)a_keys.size();
        if (parallel) {
            #pragma omp parallel for schedule(static)
            for (int k = 0; k < n; ++k) {
                size_t cell_ind = a_keys[(size_t)k];
                size_t a_cell_ind = m_a_cell_index[cell_ind];
                auto& cell = storage->surfaces[cell_ind];
                Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
                scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
            }
        }
        else
            ordered_loop(a_keys.size(), [&](size_t i) {
                if (cached) return a_mixed_rhs_cached(i, msh, storage->surfaces[a_keys[i]], a_rhs_fun);
                return Matrix<T, Dynamic, 1>(a_mixed_rhs(msh, storage->surfaces[a_keys[i]], a_rhs_fun));
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
                scatter_a_rhs_data(m_a_cell_index[a_keys[i]], msh, storage->surfaces[a_keys[i]], f_loc);
            });
        apply_bc(msh, explicit_scheme);
    }
    
    // Assembly of acoustic and elastic part separately; 
    void assemble(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme) {
        
        auto storage = msh.backend_storage();
        LHS.setZero();
        RHS.setZero();
        
        if (!m_parallel_assembly) {
            // Elastic block
            for (auto e_chunk : m_e_material) {
                size_t e_cell_ind = m_e_cell_index[e_chunk.first];
                auto& cell = storage->surfaces[e_chunk.first];
                Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_chunk.second, msh, cell, explicit_scheme);
                Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
                scatter_e_data(e_cell_ind, msh, cell, mixed_operator_loc, f_loc);
            }

            // Acoustic block
            for (auto a_chunk : m_a_material) {
                size_t a_cell_ind = m_a_cell_index[a_chunk.first];
                auto& cell = storage->surfaces[a_chunk.first];
                Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_chunk.second, msh, cell, explicit_scheme);
                Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
                scatter_a_data(a_cell_ind, msh, cell, mixed_operator_loc, f_loc);
            }

            finalize();
            return;
        }

        if (m_class_reuse) {
            // local operators once per geometric class, local vectors per cell
            auto& cl = cell_classes(msh);
            auto e_ops = e_class_operators(msh, explicit_scheme);
            ordered_loop(cl.e_ids.size(), [&](size_t i) {
                return Matrix<T, Dynamic, 1>(e_mixed_rhs(msh, storage->surfaces[cl.e_ids[i]], e_rhs_fun));
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
                scatter_e_data(m_e_cell_index[cl.e_ids[i]], msh, storage->surfaces[cl.e_ids[i]], e_ops[cl.e_class_of[i]], f_loc);
            });
            auto a_ops = a_class_operators(msh, explicit_scheme);
            ordered_loop(cl.a_ids.size(), [&](size_t i) {
                return Matrix<T, Dynamic, 1>(a_mixed_rhs(msh, storage->surfaces[cl.a_ids[i]], a_rhs_fun));
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
                scatter_a_data(m_a_cell_index[cl.a_ids[i]], msh, storage->surfaces[cl.a_ids[i]], a_ops[cl.a_class_of[i]], f_loc);
            });
            finalize_parallel();
            return;
        }

        // Elastic block
        typedef std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, 1>> local_type;
        std::vector<std::pair<size_t, elastic_material_data<T>>> e_cells(m_e_material.begin(), m_e_material.end());
        ordered_loop(e_cells.size(), [&](size_t i) {
            auto& cell = storage->surfaces[e_cells[i].first];
            elastic_material_data<T> material = e_cells[i].second;
            return local_type(e_mixed_operator(material, msh, cell, explicit_scheme), e_mixed_rhs(msh, cell, e_rhs_fun));
        }, [&](size_t i, const local_type& loc) {
            scatter_e_data(m_e_cell_index[e_cells[i].first], msh, storage->surfaces[e_cells[i].first], loc.first, loc.second);
        });

        // Acoustic block
        std::vector<std::pair<size_t, acoustic_material_data<T>>> a_cells(m_a_material.begin(), m_a_material.end());
        ordered_loop(a_cells.size(), [&](size_t i) {
            auto& cell = storage->surfaces[a_cells[i].first];
            acoustic_material_data<T> material = a_cells[i].second;
            return local_type(a_mixed_operator(material, msh, cell, explicit_scheme), a_mixed_rhs(msh, cell, a_rhs_fun));
        }, [&](size_t i, const local_type& loc) {
            scatter_a_data(m_a_cell_index[a_cells[i].first], msh, storage->surfaces[a_cells[i].first], loc.first, loc.second);
        });

        finalize();

    }
    
    // Computation of local elastic mixed operator; Strain reconstruction & Stabilization operator; 
    Matrix<T, Dynamic, Dynamic> e_mixed_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool explicit_scheme) {
            
        T rho = material.rho();
        T vp = material.vp();
        T vs = material.vs();
        T mu = rho * vs * vs;
        T lambda = rho * vp * vp - 2*mu;
    
        auto reconstruction_operator = strain_tensor_reconstruction(msh, cell);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        auto n_rows = R_operator.rows();
        auto n_cols = R_operator.cols();
        Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);

        if (explicit_scheme) {
            auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, scaled_on(msh, cell));
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }
        else {
            if(m_hho_stabilization_Q) {
                auto rec_for_stab = make_vector_hho_symmetric_laplacian(msh, cell, m_hho_di);
                // weight in 1/h_F (scaled) or O(1): diskpp scales the HHO stabilization by 1/h_T (hF = false), times h_T
                Matrix<T, Dynamic, Dynamic> stabilization_operator = make_vector_hho_stabilization(msh, cell, rec_for_stab.first, m_hho_di,
                                                                                                   m_scaled_stabilization_Q && m_stab_size != 1);
                if (!m_scaled_stabilization_Q)
                    stabilization_operator *= diameter(msh, cell);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
            else {
                auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, scaled_on(msh, cell));
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }
        
        if (m_graded_href > 0) S_operator *= std::max(T(1), m_graded_href / diameter(msh, cell));
        // COEF SPECTRAL RADIUS
        return R_operator + m_e_eta*(rho*vs)*S_operator;
    }
            
    // Computation of local symmetric tensor mass matrix; 
    Matrix<T, Dynamic, Dynamic>
    symmetric_tensor_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cell) {
            
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(ten_bs, ten_bs);
        
        auto qps = integrate(msh, cell, 2 * gradeg);
        // number of tensor components
        size_t dec = 0;
         if (dim == 3) {
             dec = 6;
         }
         else if (dim == 2) {
             dec = 3;
         }
         else {
             std::logic_error("Expected 3 >= dim > 1");
         }

         for (auto& qp : qps) {
             auto phi = ten_b.eval_functions(qp.point());
             for (size_t j = 0; j < ten_bs; j++) {
                auto qp_phi_j = disk::priv::inner_product(qp.weight(), phi[j]);
                for (size_t i = j; i < ten_bs; i += dec) {
                    mass_matrix(i, j) += disk::priv::inner_product(phi[i], qp_phi_j);
                }
             }
         }

        for (size_t j = 0; j < ten_bs; j++) {
            for (size_t i = 0; i < j; i++) {
                 mass_matrix(i, j) = mass_matrix(j, i);
            }
        }
        
        return mass_matrix;
    }
    
    // Computation of local symmetric tensor trace mass matrix;
    Matrix<T, Dynamic, Dynamic>
    symmetric_tensor_trace_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cell) {
        
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(ten_bs, ten_bs);
        
        auto qps = integrate(msh, cell, 2 * gradeg);
        // number of tensor components
        size_t dec = 0;
        if (dim == 3) {
            dec = 6;
        }
        else if (dim == 2) {
            dec = 3;
        }
        else {
            std::logic_error("Expected 3 >= dim > 1");
        }
        
        for (auto& qp : qps) {
            auto phi = ten_b.eval_functions(qp.point());
            for (size_t j = 0; j < ten_bs; j++) {
                auto identity = phi[j];
                identity.setZero();
                for(size_t d = 0; d < dim; d++) {
                    identity(d,d) = 1.0;
                }
                auto trace = phi[j].trace();
                auto trace_phi_j = disk::priv::inner_product(phi[j].trace(), identity);
                auto qp_phi_j = disk::priv::inner_product(qp.weight(), trace_phi_j);
                for (size_t i = 0; i < ten_bs; i ++) {
                    mass_matrix(i, j) += disk::priv::inner_product(phi[i], qp_phi_j);
                }
            }
        }
        
        return mass_matrix;
    }
    
    // Computation of local elastic RHS vector;
    // Assembly indices and Dirichlet faces of the cells with boundary faces (as scatter_e_bc_data, scatter_a_bc_data)
    void prepare_bc_cells(const Mesh& msh) {
        const size_t fdeg = m_hho_di.face_degree();
        const size_t n_e_fbs = disk::vector_basis_size(fdeg, Mesh::dimension - 1, Mesh::dimension);
        const size_t n_a_fbs = disk::scalar_basis_size(fdeg, Mesh::dimension - 1);
        const size_t faces0 = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
        tools::prepare_bc_cells(msh, m_rc->e_bc_cells, m_e_elements_with_bc_eges, m_e_bnd, m_e_cell_index, m_e_compress_indexes,
            get_e_cell_basis_data(), n_e_fbs, 0, faces0,
            [&](const typename Mesh::face_type& fc) { return make_vector_monomial_basis(msh, fc, fdeg); }, m_parallel_assembly);
        tools::prepare_bc_cells(msh, m_rc->a_bc_cells, m_a_elements_with_bc_eges, m_a_bnd, m_a_cell_index, m_a_compress_indexes,
            get_a_cell_basis_data(), n_a_fbs, m_n_elastic_cell_dof, faces0 + m_n_elastic_face_dof,
            [&](const typename Mesh::face_type& fc) { return make_scalar_monomial_basis(msh, fc, fdeg); }, m_parallel_assembly);
    }


    // Quadrature points, weights and basis values of the source terms of every cell (in the order of m_e_material,
    // m_a_material), computed at the first call; false (and nothing kept) above m_rhs_cache_limit doubles
    bool prepare_rhs_cache(const Mesh& msh, const std::vector<size_t>& e_keys) {
        if (m_rc->rhs_off) return false;
        if (m_rc->rhs_ready) return true;
        const auto celdeg = m_hho_di.cell_degree();
        std::vector<size_t> a_keys;
        for (auto& a_chunk : m_a_material) a_keys.push_back(a_chunk.first);
        size_t doubles = tools::quadrature_doubles(msh, e_keys, celdeg, 2 * disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension))
                       + tools::quadrature_doubles(msh, a_keys, celdeg, disk::scalar_basis_size(celdeg, Mesh::dimension));
        if (doubles > m_rhs_cache_limit) { m_rc->rhs_off = true; return false; }
        tools::cell_quadratures(msh, e_keys, celdeg, [&](const typename Mesh::cell_type& cell) { return make_vector_monomial_basis(msh, cell, celdeg); },
                                m_rc->e_rhs, m_parallel_assembly);
        tools::cell_quadratures(msh, a_keys, celdeg, [&](const typename Mesh::cell_type& cell) { return make_scalar_monomial_basis(msh, cell, celdeg); },
                                m_rc->a_rhs, m_parallel_assembly);
        m_rc->rhs_ready = true;
        return true;
    }

    // e_mixed_rhs, a_mixed_rhs (di = 0) on the cached quadrature data of the i-th cell (same operations)
    Matrix<T, Dynamic, 1> e_mixed_rhs_cached(size_t i, const Mesh& msh, const typename Mesh::cell_type& cell,
                                             std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> & rhs_fun) {
        auto ten_bs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(m_hho_di.cell_degree(), Mesh::dimension, Mesh::dimension);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(ten_bs + vec_bs);
        ret.block(ten_bs, 0, vec_bs, 1) = tools::rhs(m_rc->e_rhs[i], vec_bs, rhs_fun);
        return ret;
    }

    Matrix<T, Dynamic, 1> a_mixed_rhs_cached(size_t i, const Mesh& msh, const typename Mesh::cell_type& cell,
                                             std::function<double(const typename Mesh::point_type& )> & rhs_fun) {
        const auto rbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension) - 1;
        const auto sbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(sbs + rbs);
        ret.block(rbs, 0, sbs, 1) = tools::rhs(m_rc->a_rhs[i], sbs, rhs_fun);
        return ret;
    }

    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    e_mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0) {

        auto recdeg = m_hho_di.grad_degree();
        auto celdeg = m_hho_di.cell_degree();
        auto facdeg = m_hho_di.face_degree();

        auto ten_bs = disk::sym_matrix_basis_size(recdeg, Mesh::dimension, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        size_t n_cbs = ten_bs + vec_bs;
        auto cell_basis = make_vector_monomial_basis(msh, cell, celdeg);
        using T = typename Mesh::coordinate_type;

        Matrix<T, Dynamic, 1> ret_loc = Matrix<T, Dynamic, 1>::Zero(cell_basis.size());
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(n_cbs);

        const auto qps = integrate(msh, cell, 2 * (celdeg + di));
        for (auto& qp : qps) {
            const auto phi  = cell_basis.eval_functions(qp.point());
            const auto qp_f = disk::priv::inner_product(qp.weight(), rhs_fun(qp.point()));
            ret_loc += disk::priv::outer_product(phi, qp_f);
        }
        ret.block(ten_bs,0,vec_bs,1) = ret_loc;
        return ret;

    }
            
    // Computation of local elastic mass operator;
    Matrix<T, Dynamic, Dynamic> 
    e_mass_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool add_vector_mass_Q = true) {
            
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs = n_ten_cbs + n_vec_cbs;
            
        T rho = material.rho();
        T vp = material.vp();
        T vs = material.vs();
        T mu = rho * vs * vs;
        T lambda = rho * vp * vp - 2*mu;
    
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(n_cbs, n_cbs);
        
        // Symmetric stress tensor mass block
        // Stress tensor
        Matrix<T, Dynamic, Dynamic> mass_matrix_sigma  = symmetric_tensor_mass_matrix(msh, cell);
        // Tensor trace
        Matrix<T, Dynamic, Dynamic> mass_matrix_trace_sigma  = symmetric_tensor_trace_mass_matrix(msh, cell);
        
        // Constitutive relationship inverse
        mass_matrix_trace_sigma *= (lambda/(2.0*mu+2.0*lambda));
        mass_matrix_sigma -= mass_matrix_trace_sigma;
        mass_matrix_sigma *= (1.0/(2.0*mu));
        mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs) = mass_matrix_sigma;
        
        if (add_vector_mass_Q) {
            // vector velocity mass mass block
            auto vec_basis = disk::make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, Dynamic> mass_matrix_v = disk::make_mass_matrix(msh, cell, vec_basis);
            mass_matrix_v *= rho;
            mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs) = mass_matrix_v;
        }

        return mass_matrix;
    }
    
    // Computation of local acoustic mass operator;
    Matrix<T, Dynamic, Dynamic> 
    a_mass_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool add_scalar_mass_Q = true){
            
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        Matrix<T, Dynamic, Dynamic> mass_matrix = Matrix<T, Dynamic, Dynamic>::Zero(n_cbs, n_cbs);
        
        auto recdeg = m_hho_di.reconstruction_degree();
        auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
        auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);

        Matrix<T, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
        Matrix<T, Dynamic, Dynamic> mass_matrix_q = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
        mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);
        mass_matrix_q *= material.rho();
        mass_matrix.block(0, 0, n_vec_cbs, n_vec_cbs) = mass_matrix_q;
        
        if (add_scalar_mass_Q) {
            auto scal_basis = disk::make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, Dynamic> mass_matrix_v = disk::make_mass_matrix(msh, cell, scal_basis);
            mass_matrix_v *= (1.0/(material.rho()*material.vp()*material.vp()));
            mass_matrix.block(n_vec_cbs, n_vec_cbs, n_scal_cbs, n_scal_cbs) = mass_matrix_v;
        }

        return mass_matrix;
    }
    
    // Computation of local strain tensor reconstruction;
    std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>, Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
    strain_tensor_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell) {

        using T = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;

        const size_t N = Mesh::dimension;

        auto graddeg = m_hho_di.grad_degree();
        auto celdeg  = m_hho_di.cell_degree();
        auto facdeg  = m_hho_di.face_degree();

        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, graddeg);
        auto vec_b = make_vector_monomial_basis(msh, cell, celdeg);

        auto ten_bs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
        auto vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        auto fbs = disk::vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

        auto num_faces = howmany_faces(msh, cell);

        matrix_type gr_lhs = matrix_type::Zero(ten_bs, ten_bs);
        matrix_type gr_rhs = matrix_type::Zero(ten_bs, vec_bs + num_faces * fbs);
        
        const auto qps = integrate(msh, cell, 2 * graddeg);
        size_t dec = 0;
         if (N == 3) {
             dec = 6;
         }
         else if (N == 2) {
             dec = 3;
         }
         else {
             std::logic_error("Expected 3 >= dim > 1");
         }

         for (auto& qp : qps) {
             const auto gphi = ten_b.eval_functions(qp.point());
             for (size_t j = 0; j < ten_bs; j++) {
                auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (size_t i = j; i < ten_bs; i += dec) {
                         gr_lhs(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
                }
             }
         }

         for (size_t j = 0; j < ten_bs; j++) {    
             for (size_t i = 0; i < j; i++) {
                 gr_lhs(i, j) = gr_lhs(j, i);
             }
         }
         
         if (celdeg > 0) {
             const auto qpc = integrate(msh, cell, graddeg + celdeg - 1);
             for (auto& qp : qpc) {
                 const auto gphi    = ten_b.eval_functions(qp.point());
                 const auto dphi    = vec_b.eval_sgradients(qp.point());
                 const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);
                 gr_rhs.block(0, 0, ten_bs, vec_bs) += disk::priv::outer_product(gphi, qp_dphi);
             }
         }

         const auto fcs = faces(msh, cell);
         for (size_t i = 0; i < fcs.size(); i++) {
             const auto fc = fcs[i];
             const auto n  = normal(msh, cell, fc);
             const auto fb = make_vector_monomial_basis(msh, fc, facdeg);
             const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
             for (auto& qp : qps_f) {
                 const auto gphi = ten_b.eval_functions(qp.point());
                 const auto cphi = vec_b.eval_functions(qp.point());
                 const auto fphi = fb.eval_functions(qp.point());
                 const auto qp_gphi_n = disk::priv::inner_product(gphi, disk::priv::inner_product(qp.weight(), n));
                 gr_rhs.block(0, vec_bs + i * fbs, ten_bs, fbs) += disk::priv::outer_product(qp_gphi_n, fphi);
                 gr_rhs.block(0, 0, ten_bs, vec_bs) -= disk::priv::outer_product(qp_gphi_n, cphi);
             }
         }
            
        auto n_rows = gr_rhs.cols() + ten_bs;
        auto n_cols = gr_rhs.cols() + ten_bs;

        // Shrinking data
        matrix_type data_mixed = matrix_type::Zero(n_rows,n_cols);
        data_mixed.block(0, (ten_bs), ten_bs, n_cols-(ten_bs)) = -gr_rhs;
        data_mixed.block((ten_bs), 0, n_rows-(ten_bs), ten_bs) = gr_rhs.transpose();

        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
    
    // Computation of local acoustic mixed operator; Gradient reconstruction & Stabilization operator; 
    Matrix<T, Dynamic, Dynamic> 
    a_mixed_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell, bool explicit_scheme){
        
        auto reconstruction_operator = mixed_scalar_reconstruction(msh, cell);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        auto n_rows = R_operator.rows();
        auto n_cols = R_operator.cols();

        Matrix<T, Dynamic, Dynamic> S_operator = Matrix<T, Dynamic, Dynamic>::Zero(n_rows, n_cols);
        if (explicit_scheme) {
            auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, scaled_on(msh, cell));
            auto n_s_rows = stabilization_operator.rows();
            auto n_s_cols = stabilization_operator.cols();
            S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
        }
        else {
            if (m_hho_stabilization_Q) {
                // weight in 1/h_F (scaled) or O(1): diskpp scales the HHO stabilization by 1/h_T (hF = false), times h_T
                // diskpp: 1/h_T with hho_degree_info, the cell degree information variant scales by 1/h_F
                Matrix<T, Dynamic, Dynamic> stabilization_operator = m_scaled_stabilization_Q && m_stab_size == 2
                    ? make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, disk::CellDegreeInfo<Mesh>(msh, cell, m_hho_di.cell_degree(), m_hho_di.face_degree()), true)
                    : make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di);
                if (!m_scaled_stabilization_Q)
                    stabilization_operator *= diameter(msh, cell);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            } 
            else {
                auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, scaled_on(msh, cell));
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }

        if (m_graded_href > 0) S_operator *= std::max(T(1), m_graded_href / diameter(msh, cell));
        // COEF SPECTRAL RADIUS
        T rho = material.rho();
        T vp  = material.vp();

        return R_operator + (m_a_eta/(vp*rho))*S_operator;
    }

    // Computation of local scalar (acoustic) mixed reconstruction;
    std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>, Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
    mixed_scalar_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cell) {

        using T = typename Mesh::coordinate_type;
        typedef Matrix<T, Dynamic, Dynamic> matrix_type;
        typedef Matrix<T, Dynamic, 1>       vector_type;

        const size_t DIM = Mesh::dimension;

        const auto recdeg = m_hho_di.reconstruction_degree();
        const auto celdeg = m_hho_di.cell_degree();
        const auto facdeg = m_hho_di.face_degree();

        auto cb = make_scalar_monomial_basis(msh, cell, recdeg);

        const auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
        const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension);
        const auto fbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        const auto num_faces = howmany_faces(msh, cell);

        const matrix_type stiff  = make_stiffness_matrix(msh, cell, cb);
        matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
        matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

        gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
        gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

        const auto fcs = faces(msh, cell);
        for (size_t i = 0; i < fcs.size(); i++) {

            const auto fc = fcs[i];
            const auto n  = normal(msh, cell, fc);
            auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

            auto qps_f = integrate(msh, fc, recdeg - 1 + std::max(facdeg,celdeg));
            for (auto& qp : qps_f) {
                vector_type c_phi_tmp = cb.eval_functions(qp.point());
                vector_type c_phi = c_phi_tmp.head(cbs);
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type f_phi = fb.eval_functions(qp.point());
                gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
                gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
            }
        }

        auto vec_cell_size = gr_lhs.cols();
        auto nrows = gr_rhs.cols()+vec_cell_size;
        auto ncols = gr_rhs.cols()+vec_cell_size;
        
        // Shrinking data
        matrix_type data_mixed = matrix_type::Zero(nrows,ncols);
        data_mixed.block(0, vec_cell_size, vec_cell_size, ncols-vec_cell_size) = -gr_rhs;
        data_mixed.block(vec_cell_size, 0, nrows-vec_cell_size, vec_cell_size) = gr_rhs.transpose();
        
        matrix_type oper = gr_lhs.llt().solve(gr_rhs);
        return std::make_pair(oper, data_mixed);
    }
    
    // Computation of local acoustic RHS vector;
    Matrix<typename Mesh::coordinate_type, Dynamic, 1>
    a_mixed_rhs(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<double(const typename Mesh::point_type& )> & rhs_fun, size_t di = 0) {

        const auto recdeg = m_hho_di.reconstruction_degree();
        const auto celdeg = m_hho_di.cell_degree();
        const auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension) - 1;
        const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension) + rbs;
        auto cell_basis = make_scalar_monomial_basis(msh, cell, celdeg);
        using T = typename Mesh::coordinate_type;

        Matrix<T, Dynamic, 1> ret_loc = Matrix<T, Dynamic, 1>::Zero(cell_basis.size());
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

        const auto qps = integrate(msh, cell, 2 * (celdeg + di));
        for (auto& qp : qps) {
            const auto phi  = cell_basis.eval_functions(qp.point());
            const auto qp_f = disk::priv::inner_product(qp.weight(), rhs_fun(qp.point()));
            ret_loc += disk::priv::outer_product(phi, qp_f);
        }
        ret.block(rbs,0,cell_basis.size(),1) = ret_loc;
        return ret;
    }
    
    // Function projection onto cell dofs 
    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> v_fun, std::function<disk::static_matrix<T,2,2>(const typename Mesh::point_type& )> flux_fun, std::function<T(const typename Mesh::point_type& )> v_s_fun, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> flux_s_fun) {
        
        auto storage = msh.backend_storage();
        size_t n_dof = MASS.rows();
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);

        // elastic block
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_e_cbs = n_ten_cbs + n_vec_cbs;
    
        if (m_parallel_assembly) {
            // fast mode: projections of the cells by the OpenMP threads (independent cells: the same vector)
            std::vector<size_t> e_ids, a_ids;
            for (auto& e_chunk : m_e_material) e_ids.push_back(e_chunk.first);
            for (auto& a_chunk : m_a_material) a_ids.push_back(a_chunk.first);
            tools::ordered_loop(true, e_ids.size(), [&](size_t i) {
                auto& cell = storage->surfaces[e_ids[i]];
                Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_e_cbs);
                x_proj_dof.block(0, 0, n_ten_cbs, 1)          = project_ten_function(msh, cell, flux_fun);
                x_proj_dof.block(n_ten_cbs, 0, n_vec_cbs, 1)  = project_function(msh, cell, m_hho_di.cell_degree(), v_fun);
                return x_proj_dof;
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& x_proj_dof) {
                scatter_e_cell_dof_data(m_e_cell_index[e_ids[i]], msh, storage->surfaces[e_ids[i]], x_glob, x_proj_dof);
            });
            size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
            size_t n_v_s_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
            tools::ordered_loop(true, a_ids.size(), [&](size_t i) {
                auto& cell = storage->surfaces[a_ids[i]];
                Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_v_s_cbs + n_scal_cbs);
                x_proj_dof.block(0, 0, n_v_s_cbs, 1) = project_vec_function(msh, cell, flux_s_fun);
                x_proj_dof.block(n_v_s_cbs, 0, n_scal_cbs, 1) = project_function(msh, cell, m_hho_di.cell_degree(), v_s_fun);
                return x_proj_dof;
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& x_proj_dof) {
                scatter_a_cell_dof_data(m_a_cell_index[a_ids[i]], msh, storage->surfaces[a_ids[i]], x_glob, x_proj_dof);
            });
            return;
        }

        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_ten_dof = project_ten_function(msh, cell, flux_fun);
            Matrix<T, Dynamic, 1> x_proj_vec_dof = project_function(msh, cell, m_hho_di.cell_degree(), v_fun);
            Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_e_cbs);
            x_proj_dof.block(0, 0, n_ten_cbs, 1)          = x_proj_ten_dof;
            x_proj_dof.block(n_ten_cbs, 0, n_vec_cbs, 1)  = x_proj_vec_dof;
            scatter_e_cell_dof_data(e_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
        
        // acoustic block
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_v_s_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_a_cbs = n_v_s_cbs + n_scal_cbs;
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_vec_dof = project_vec_function(msh, cell, flux_s_fun);
            Matrix<T, Dynamic, 1> x_proj_sca_dof = project_function(msh, cell, m_hho_di.cell_degree(), v_s_fun);
            Matrix<T, Dynamic, 1> x_proj_dof = Matrix<T, Dynamic, 1>::Zero(n_a_cbs);
            x_proj_dof.block(0, 0, n_v_s_cbs, 1) = x_proj_vec_dof;
            x_proj_dof.block(n_v_s_cbs, 0, n_scal_cbs, 1) = x_proj_sca_dof;
            scatter_a_cell_dof_data(a_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
    
    }
    
    // Vector function projection onto cell dofs.
    Matrix<T, Dynamic, 1> project_vec_function(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun) {
    
            auto recdeg = m_hho_di.reconstruction_degree();
            auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
            auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
            Matrix<T, Dynamic, Dynamic> mass_matrix_q = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
            mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);

            Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(rbs-1);
            Matrix<T, 1, 2> f_vec = Matrix<T, Dynamic, Dynamic>::Zero(1, 2);
            const auto qps = integrate(msh, cell, 2*recdeg);
            for (auto& qp : qps) {
                auto dphi = rec_basis.eval_gradients(qp.point());
                f_vec(0,0) = vec_fun(qp.point())[0];
                f_vec(0,1) = vec_fun(qp.point())[1];
                for (size_t i = 0; i < rbs-1; i++){
                    Matrix<T, 2, 1> phi_i = dphi.block(i+1, 0, 1, 2).transpose();
                    rhs(i,0) = rhs(i,0) + (qp.weight() * f_vec*phi_i)(0,0);
                }
            }
            Matrix<T, Dynamic, 1> x_dof = mass_matrix_q.llt().solve(rhs);
            return x_dof;
    }

    // Symmetric tensor function projection onto cell dofs.
    Matrix<T, Dynamic, 1> project_ten_function(const Mesh& msh, const typename Mesh::cell_type& cell, std::function<disk::static_matrix<T, 2,2>(const typename Mesh::point_type& )> ten_fun) {
    
        Matrix<T, Dynamic, Dynamic> mass_matrix = symmetric_tensor_mass_matrix(msh, cell);
        size_t dim = Mesh::dimension;
        auto gradeg = m_hho_di.grad_degree();
        auto ten_bs = disk::sym_matrix_basis_size(gradeg, dim, dim);
        auto ten_b = make_sym_matrix_monomial_basis(msh, cell, gradeg);
        Matrix<T, Dynamic, 1> rhs = Matrix<T, Dynamic, 1>::Zero(ten_bs);

        const auto qps = integrate(msh, cell, 2 * gradeg);
        for (auto& qp : qps) {
            auto phi = ten_b.eval_functions(qp.point());
            disk::static_matrix<T, 2,2> sigma = ten_fun(qp.point());
            for (size_t i = 0; i < ten_bs; i++){
                auto qp_phi_i = disk::priv::inner_product(qp.weight(), phi[i]);
                rhs(i,0) += disk::priv::inner_product(qp_phi_i,sigma);
            }
        }
        Matrix<T, Dynamic, 1> x_dof = mass_matrix.llt().solve(rhs);
        return x_dof;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// Sparse matrices triplet assembly finalization //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void finalize() {
        LHS.setFromTriplets( m_triplets.begin(), m_triplets.end() );
        std::vector<Triplet<T>>().swap(m_triplets);    // clear and release the memory
        LHS_STAB.setFromTriplets( m_triplets_stab.begin(), m_triplets_stab.end() );
        std::vector<Triplet<T>>().swap(m_triplets_stab);    // clear and release the memory
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Stabilization setting ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void set_coupling_stabilization() {
        m_hho_stabilization_Q = false;
        clear_rhs_cache();
    }   
    
    void set_scaled_stabilization() {
        m_scaled_stabilization_Q = true;
        clear_rhs_cache();
    }

    // Local matrices and vectors of assemble, assemble_mass, assemble_rhs computed by the OpenMP threads (same global
    // matrices and vectors as the serial loops, bit for bit)
    void set_parallel_assembly(bool parallel) {
        m_parallel_assembly = parallel;
    }

    // Fast mode, with set_parallel_assembly: local operators and mass matrices computed once per geometric class of cells
    // (translated copies with the same material, assembly_tools::cell_classes), equal to those of each cell up to
    // round-off
    void set_class_reuse(bool reuse) {
        m_class_reuse = reuse;
    }

    // assemble_rhs called many times (source or Dirichlet data depending on time): the local operators of the cells with
    // Dirichlet faces, which do not depend on time, are computed once, and the quadrature points, weights, basis values
    // and face masses of the source terms and Dirichlet data are kept (same arithmetic: the same RHS, bit for bit).
    // Cleared by the stabilization setters of this class (call set_cache_rhs_data again after those of the base class).
    void set_cache_rhs_data(bool cache) {
        m_rc = cache ? std::make_unique<rhs_cache_data>() : nullptr;
    }

    void clear_rhs_cache() {
        if (m_rc) m_rc = std::make_unique<rhs_cache_data>();
    }


    // Length h of the 1/h scaling: 0 as in diskpp (here: cell diameter h_T, except for the vector HHO stabilization:
    // face diameter h_F), 1 cell diameter h_T, 2 face diameter h_F
    // Scaled (1/h) least-squares stabilization on the cell: set_scaled_stabilization, or a cell smaller than the size of
    // set_stabilization_hybrid
    bool scaled_on(const Mesh& msh, const typename Mesh::cell_type& cell) const {
        return m_scaled_stabilization_Q || (m_hybrid_hc > 0 && diameter(msh, cell) < m_hybrid_hc);
    }

    // Hybrid weight (space.scaling = "hybrid"): 1/h on the cells smaller than hc (the fine region of the LI scheme, around
    // a singularity: consistency of the stabilization), O(1) on the others (hyperbolic CFL of the explicit part)
    void set_stabilization_hybrid(T hc) {
        m_hybrid_hc = hc;
    }

    // Graded O(1) weight (space.scaling = "graded"): the O(1) stabilization of a cell T multiplied by max(1, href / h_T):
    // O(1) on the cells of size href, 1/h_T on the smaller ones (meshes graded towards a singularity)
    void set_stabilization_graded(T href) {
        m_graded_href = href;
    }

    void set_stabilization_size(int size) {
        m_stab_size = size;
        clear_rhs_cache();
    }

    // Weights eta of the elastic and acoustic stabilizations (default 1.5 and 0.8)
    void set_stabilization_weights(T e_eta, T a_eta) {
        m_e_eta = e_eta;
        m_a_eta = a_eta;
        clear_rhs_cache();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////// Getter functions for DOFs, basis sizes, materials, boundaries, and interface data //////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    size_t get_cell_basis_data(){
        
        // Elastic cell basis
        size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_cbs1 = n_ten_cbs + n_vec_cbs;
        
        // Acoustic cell basis 
        size_t n_vel_scal_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_cbs2 = n_vel_scal_cbs + n_scal_cbs;
        
        size_t n_cbs = n_cbs1 + n_cbs2;
        
        return n_cbs;
    }
    
    size_t get_n_faces() {
        return m_n_edges - m_n_essential_edges;
    }
    
    size_t get_face_basis_data() {
        size_t n_fbs1 = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        size_t n_fbs2 = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        size_t n_fbs = n_fbs1 + n_fbs2;
        return n_fbs;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////// LTS data structure
    
// Multilevel LTS: level of every unknown from the level of every cell (cell_level, indexed by mesh cell id): the unknowns
// of a cell take its level, those of a face the largest level of its cells (as the fine faces of mark_fine_region);
// Dirichlet faces have no unknowns
std::vector<int> dof_levels(const Mesh& msh, const std::vector<int>& cell_level) {
    const size_t n_e_cbs = get_e_cell_basis_data(), n_a_cbs = get_a_cell_basis_data();
    const size_t n_e_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
    const size_t n_a_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
    const size_t e_face_shift = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
    const size_t a_face_shift = e_face_shift + m_n_elastic_face_dof;
    std::vector<int> level(e_face_shift + m_n_elastic_face_dof + m_n_acoustic_face_dof, 0);
    std::vector<int> e_face_level(n_e_edges, 0), a_face_level(n_a_edges, 0);
    auto storage = msh.backend_storage();
    for (auto& chunk : m_e_material) {
        const size_t ofs = m_e_cell_index[chunk.first] * n_e_cbs;
        const int l = cell_level.at(chunk.first);
        for (size_t i = 0; i < n_e_cbs; i++) level[ofs + i] = l;
        for (auto& fc : faces(msh, storage->surfaces[chunk.first])) {
            const auto fc_id = msh.lookup(fc);
            if (m_e_bnd.is_dirichlet_face(fc_id)) continue;
            int& fl = e_face_level[m_e_compress_indexes.at(fc_id)];
            fl = std::max(fl, l);
        }
    }
    for (auto& chunk : m_a_material) {
        const size_t ofs = m_n_elastic_cell_dof + m_a_cell_index[chunk.first] * n_a_cbs;
        const int l = cell_level.at(chunk.first);
        for (size_t i = 0; i < n_a_cbs; i++) level[ofs + i] = l;
        for (auto& fc : faces(msh, storage->surfaces[chunk.first])) {
            const auto fc_id = msh.lookup(fc);
            if (m_a_bnd.is_dirichlet_face(fc_id)) continue;
            int& fl = a_face_level[m_a_compress_indexes.at(fc_id)];
            fl = std::max(fl, l);
        }
    }
    for (size_t cf = 0; cf < n_e_edges; cf++)
        for (size_t i = 0; i < n_e_fbs; i++) level[e_face_shift + cf * n_e_fbs + i] = e_face_level[cf];
    for (size_t cf = 0; cf < n_a_edges; cf++)
        for (size_t i = 0; i < n_a_fbs; i++) level[a_face_shift + cf * n_a_fbs + i] = a_face_level[cf];
    return level;
}

void assemble_P(const Mesh& msh, T h_c, size_t nb_layer = 1) {
        
        Pfine.setZero();
        Pcoarse.setZero();
        
        // --------------------------------------------------
        // DOF SIZES
        // --------------------------------------------------
        const size_t n_e_cbs = get_e_cell_basis_data();
        const size_t n_a_cbs = get_a_cell_basis_data();
        
        const size_t n_e_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        const size_t n_a_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
        auto storage = msh.backend_storage();
        
        // --------------------------------------------------
        // STEPS 1-3 : fine cells (mesh size, nb_layer vertex layers) and fine faces
        // --------------------------------------------------
        std::vector<bool> e_cell_is_fine, a_cell_is_fine, e_face_is_fine, a_face_is_fine;
        mark_fine_region(msh, h_c, nb_layer, e_cell_is_fine, a_cell_is_fine, e_face_is_fine, a_face_is_fine);
        
        // --------------------------------------------------
        // STEP 4 : fill triplets for Pfine and Pcoarse
        // --------------------------------------------------
        std::vector<Triplet<T>> Pfine_trips;
        std::vector<Triplet<T>> Pcoarse_trips;
        Pfine_trips.reserve(LHS.rows());
        Pcoarse_trips.reserve(LHS.rows());
        
        // Elastic cells
        for (auto& chunk : m_e_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_e_cell_index[cell_id];
            const size_t glob_ofs = cell_ind * n_e_cbs;
            for (size_t i = 0; i < n_e_cbs; ++i) {
                const size_t dof = glob_ofs + i;
                if (e_cell_is_fine[cell_ind])
                    Pfine_trips.emplace_back(dof, dof, T(1));
                else
                    Pcoarse_trips.emplace_back(dof, dof, T(1));
            }
        }
        
        // Acoustic cells
        for (auto& chunk : m_a_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_a_cell_index[cell_id];
            const size_t glob_ofs = m_n_elastic_cell_dof + cell_ind * n_a_cbs;
            for (size_t i = 0; i < n_a_cbs; ++i) {
                const size_t dof = glob_ofs + i;
                if (a_cell_is_fine[cell_ind])
                    Pfine_trips.emplace_back(dof, dof, T(1));
                else
                    Pcoarse_trips.emplace_back(dof, dof, T(1));
            }
        }
        
        // Elastic faces
        const size_t e_face_shift = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
        for (size_t cf = 0; cf < n_e_edges; ++cf) {
            for (size_t i = 0; i < n_e_fbs; ++i) {
                const size_t dof = e_face_shift + cf * n_e_fbs + i;
                if (e_face_is_fine[cf])
                    Pfine_trips.emplace_back(dof, dof, T(1));
                else
                    Pcoarse_trips.emplace_back(dof, dof, T(1));
            }
        }
        
        // Acoustic faces
        const size_t a_face_shift = e_face_shift + m_n_elastic_face_dof;
        for (size_t cf = 0; cf < n_a_edges; ++cf) {
            for (size_t i = 0; i < n_a_fbs; ++i) {
                const size_t dof = a_face_shift + cf * n_a_fbs + i;
                if (a_face_is_fine[cf])
                    Pfine_trips.emplace_back(dof, dof, T(1));
                else
                    Pcoarse_trips.emplace_back(dof, dof, T(1));
            }
        }
        
        // --------------------------------------------------
        // STEP 5 : expose fine/coarse flag for Silo post-processing
        // Indexed by global cell_id (msh order), consistent with Pfine/Pcoarse
        // including protection layers.
        // --------------------------------------------------
        cell_is_fine_silo.assign(storage->surfaces.size(), false);
        for (auto& chunk : m_e_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_e_cell_index[cell_id];
            cell_is_fine_silo[cell_id] = e_cell_is_fine[cell_ind];
        }
        for (auto& chunk : m_a_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_a_cell_index[cell_id];
            cell_is_fine_silo[cell_id] = a_cell_is_fine[cell_ind];
        }

        // --------------------------------------------------
        // FINALIZE
        // --------------------------------------------------
        Pfine.setFromTriplets(Pfine_trips.begin(), Pfine_trips.end());
        Pcoarse.setFromTriplets(Pcoarse_trips.begin(), Pcoarse_trips.end());
        Pfine.makeCompressed();
        Pcoarse.makeCompressed();

        // Cell and face blocks of the projectors (two-level LTS: erk_coupling_hho_scheme::compute_wi, erk_euler_LTS)
        const size_t n_cell_dof = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
        const size_t n_face_dof = m_n_elastic_face_dof + m_n_acoustic_face_dof;
        P_cell       = Pfine.topLeftCorner(n_cell_dof, n_cell_dof);
        IminusP_cell = Pcoarse.topLeftCorner(n_cell_dof, n_cell_dof);
        Pfacefine    = Pfine.bottomRightCorner(n_face_dof, n_face_dof);
        Pfacecoarse  = Pcoarse.bottomRightCorner(n_face_dof, n_face_dof);
    }


private:
    // finalize, finalize_mass of the fast mode: the same matrices, bit for bit, built by the OpenMP threads
    void finalize_parallel() {
        tools::set_from_triplets(LHS, m_triplets);
        std::vector<Triplet<T>>().swap(m_triplets);
        LHS_STAB.setFromTriplets(m_triplets_stab.begin(), m_triplets_stab.end());
        std::vector<Triplet<T>>().swap(m_triplets_stab);
    }
    void finalize_mass_parallel() {
        tools::set_from_triplets(MASS, m_mass_triplets);
        std::vector<Triplet<T>>().swap(m_mass_triplets);
    }

    // geometric classes of the cells (set_class_reuse), computed at the first use
    const typename tools::classes& cell_classes(const Mesh& msh) {
        if (!m_cls) m_cls = std::make_unique<typename tools::classes>(tools::make_classes(msh, m_e_material, m_a_material));
        return *m_cls;
    }
    std::vector<Matrix<T, Dynamic, Dynamic>> e_class_operators(const Mesh& msh, bool explicit_scheme) {
        auto& cl = cell_classes(msh);
        auto storage = msh.backend_storage();
        return tools::class_operators(cl.e_rep, [&](size_t r) {
            elastic_material_data<T> material = m_e_material.at(cl.e_ids[r]);
            return Matrix<T, Dynamic, Dynamic>(e_mixed_operator(material, msh, storage->surfaces[cl.e_ids[r]], explicit_scheme));
        });
    }
    std::vector<Matrix<T, Dynamic, Dynamic>> a_class_operators(const Mesh& msh, bool explicit_scheme) {
        auto& cl = cell_classes(msh);
        auto storage = msh.backend_storage();
        return tools::class_operators(cl.a_rep, [&](size_t r) {
            acoustic_material_data<T> material = m_a_material.at(cl.a_ids[r]);
            return Matrix<T, Dynamic, Dynamic>(a_mixed_operator(material, msh, storage->surfaces[cl.a_ids[r]], explicit_scheme));
        });
    }

    // set_cache_rhs_data: local operators of the cells with Dirichlet faces (scheme flag of their computation, -1: none),
    // quadrature points, weights and basis values of the source terms (indexed as m_e_material, m_a_material; not
    // kept beyond m_rhs_cache_limit doubles), assembly indices and Dirichlet faces (quadrature data, factorized mass)
    // of the cells with boundary faces. One pointer, declared last: the layout of the other members is unchanged (the
    // results of the assembly were seen to depend on it, under -ffast-math)
    struct rhs_cache_data {
        std::vector<Matrix<T, Dynamic, Dynamic>> e_bc_ops, a_bc_ops;
        int bc_ops_scheme = -1;
        typename tools::e_rhs_data e_rhs;
        typename tools::a_rhs_data a_rhs;
        bool rhs_ready = false, rhs_off = false;
        typename tools::e_bc_data e_bc_cells;
        typename tools::a_bc_data a_bc_cells;
    };
    static constexpr size_t m_rhs_cache_limit = size_t(1) << 27;    // doubles (1 GB)
    std::unique_ptr<rhs_cache_data> m_rc;
    // fast mode (set_class_reuse); declared last
    bool m_class_reuse = false;
    std::unique_ptr<typename tools::classes> m_cls;
};

#endif /* elastoacoustic_four_fields_assembler_hpp */
