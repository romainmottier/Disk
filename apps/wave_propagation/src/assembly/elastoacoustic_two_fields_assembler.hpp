//
//  elastoacoustic_two_fields_assembler.hpp
//  acoustics
//
//  Created by Omar Durán on 6/3/20.
//

#pragma once
#ifndef elastoacoustic_two_fields_assembler_hpp
#define elastoacoustic_two_fields_assembler_hpp

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
class elastoacoustic_two_fields_assembler : public elastoacoustic_assembler_base<Mesh>
{
    
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
    using base::get_e_cell_index;
    using base::get_a_cell_index;
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
    using base::finalize;

    // Cell unknowns: displacement (elastic), scalar potential (acoustic)
    static size_t e_cell_basis_size(const disk::hho_degree_info& hho_di) {
        return disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
    }
    
    static size_t a_cell_basis_size(const disk::hho_degree_info& hho_di) {
        return disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
    }

    elastoacoustic_two_fields_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const e_boundary_type& e_bnd, const a_boundary_type& a_bnd, std::map<size_t,elastic_material_data<T>> & e_material, std::map<size_t,acoustic_material_data<T>> & a_material)
        : base(msh, hho_di, e_bnd, a_bnd, e_material, a_material, e_cell_basis_size(hho_di), a_cell_basis_size(hho_di))
    {
    }

    // Weight gamma of the stabilization in the stiffness: 2 mu (R + gamma S) + lambda D (elastic), (R + gamma S) / rho
    // (acoustic). gamma = 1 by default; the face splitting of Steins et al. (M2AN 2023) converges for gamma > gamma*.
    void set_stabilization_weight(T gamma) {
        m_stab_weight = gamma;
        clear_rhs_cache();
    }

    T stabilization_weight() const {
        return m_stab_weight;
    }

    // Least-squares stabilization weighted by 1/h (default) or O(1); the HHO stabilization is always weighted by 1/h
    void set_stabilization_scaling(bool scaled) {
        m_scaled_stab = scaled;
        clear_rhs_cache();
    }

    // Length h of the 1/h scaling: 0 as in diskpp (face diameter h_F, except for the scalar HHO stabilization: cell
    // diameter h_T), 1 cell diameter h_T, 2 face diameter h_F
    void set_stabilization_size(int size) {
        m_stab_size = size;
        clear_rhs_cache();
    }

    // Fast mode (solver.fast): local matrices and vectors of assemble, assemble_mass, assemble_rhs computed by the
    // OpenMP threads (same global matrices and vectors as the serial loops, bit for bit)
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
    // Cleared by the stabilization setters (call set_cache_rhs_data again after those of the base class).
    void set_cache_rhs_data(bool cache) {
        m_rc = cache ? std::make_unique<rhs_cache_data>() : nullptr;
    }

    void clear_rhs_cache() {
        if (m_rc) m_rc = std::make_unique<rhs_cache_data>();
    }

    // Face-face part sigma* of the stabilization of a cell with unit weight (Steins et al., M2AN 2023, eq. (31)):
    // c_T sum_{F in F_T} h_F^{-1} (v_F, w_F)_F (no h_F^{-1} with the O(1) scaling) on the face unknowns of the cell
    // (faces(msh, cell) order), with
    // c_T = 2 mu (elastic) or 1 / rho (acoustic), the scaling of S in e_laplacian_operator / a_laplacian_operator.
    // Block-diagonal; it is the whole face-face block of the stabilization for the HDG (mixed-order) stabilization.
    Matrix<T, Dynamic, Dynamic> e_face_stabilization_mass(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
        T mu = material.rho()*material.vs()*material.vs();
        const size_t vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto fcs = faces(msh, cell);
        Matrix<T, Dynamic, Dynamic> D = Matrix<T, Dynamic, Dynamic>::Zero(fcs.size() * vfbs, fcs.size() * vfbs);
        for (size_t i = 0; i < fcs.size(); i++) {
            auto fb = make_vector_monomial_basis(msh, fcs[i], m_hho_di.face_degree());
            T h = m_scaled_stab ? (m_stab_size == 1 ? diameter(msh, cell) : diameter(msh, fcs[i])) : T(1);
            D.block(i * vfbs, i * vfbs, vfbs, vfbs) = (2.0 * mu / h) * make_mass_matrix(msh, fcs[i], fb);
        }
        return D;
    }

    Matrix<T, Dynamic, Dynamic> a_face_stabilization_mass(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
        const size_t sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        auto fcs = faces(msh, cell);
        Matrix<T, Dynamic, Dynamic> D = Matrix<T, Dynamic, Dynamic>::Zero(fcs.size() * sfbs, fcs.size() * sfbs);
        for (size_t i = 0; i < fcs.size(); i++) {
            auto fb = make_scalar_monomial_basis(msh, fcs[i], m_hho_di.face_degree());
            T h = m_scaled_stab ? (m_stab_size == 1 ? diameter(msh, cell) : diameter(msh, fcs[i])) : T(1);
            D.block(i * sfbs, i * sfbs, sfbs, sfbs) = (1.0 / (material.rho() * h)) * make_mass_matrix(msh, fcs[i], fb);
        }
        return D;
    }

    // Global sigma* (unit weight) on the face unknowns, from the local blocks of e/a_face_stabilization_mass
    SparseMatrix<T> face_stabilization_mass(const Mesh& msh) {

        auto storage = msh.backend_storage();
        const size_t vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        const size_t sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        const size_t e_face_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
        const size_t a_face_offset = e_face_offset + m_n_elastic_face_dof;
        std::vector<Triplet<T>> triplets;

        auto scatter = [&](size_t offset, const Matrix<T, Dynamic, Dynamic>& D, size_t i_loc, size_t fbs) {
            for (size_t i = 0; i < fbs; i++)
                for (size_t j = 0; j < fbs; j++)
                    triplets.push_back(Triplet<T>(offset + i, offset + j, D(i_loc + i, i_loc + j)));
        };

        for (auto& e_chunk : m_e_material) {
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> D = e_face_stabilization_mass(e_chunk.second, msh, cell);
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++) {
                auto fc_id = msh.lookup(fcs[i]);
                if (!m_e_bnd.is_dirichlet_face(fc_id))
                    scatter(e_face_offset + m_e_compress_indexes.at(fc_id) * vfbs, D, i * vfbs, vfbs);
            }
        }

        for (auto& a_chunk : m_a_material) {
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, Dynamic> D = a_face_stabilization_mass(a_chunk.second, msh, cell);
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++) {
                auto fc_id = msh.lookup(fcs[i]);
                if (!m_a_bnd.is_dirichlet_face(fc_id))
                    scatter(a_face_offset + m_a_compress_indexes.at(fc_id) * sfbs, D, i * sfbs, sfbs);
            }
        }

        SparseMatrix<T> D(LHS.rows(), LHS.cols());
        D.setFromTriplets(triplets.begin(), triplets.end());
        return D;
    }
    
    void scatter_e_interface_neumann_data(const Mesh& msh, const typename Mesh::face_type& face,
             const Matrix<T, Dynamic, 1>& interface_bc_data)
    {
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

        std::vector<assembly_index> asm_map_e;
        asm_map_e.reserve(vfbs);
        
        auto fc_id = msh.lookup(face);
        
        auto e_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*vfbs;
        bool e_dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < vfbs; i++){
            asm_map_e.push_back( assembly_index(e_face_LHS_offset+i, !e_dirichlet) );
        }
                    
        assert( asm_map_e.size() == interface_bc_data.rows());

        for (size_t i = 0; i < interface_bc_data.rows(); i++)
        {
            if (!asm_map_e[i].assemble())
                 continue;
             RHS(int(asm_map_e[i])) += interface_bc_data(i,0);
        }

    }
    
    
    void scatter_a_interface_neumann_data(const Mesh& msh, const typename Mesh::face_type& face,
             const Matrix<T, Dynamic, 1>& interface_bc_data)
    {

        auto sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        std::vector<assembly_index> asm_map_a;
        asm_map_a.reserve(sfbs);
        
        auto fc_id = msh.lookup(face);
        
        auto a_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*sfbs;
        bool a_dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < sfbs; i++){
            asm_map_a.push_back( assembly_index(a_face_LHS_offset+i, !a_dirichlet) );
        }
        
        assert( asm_map_a.size() == interface_bc_data.rows() );

        for (size_t i = 0; i < interface_bc_data.rows(); i++)
        {
            if (!asm_map_a[i].assemble())
                 continue;
             RHS(int(asm_map_a[i])) += interface_bc_data(i,0);
        }

    }
    
    void assemble(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun){
        
        auto storage = msh.backend_storage();
        LHS.setZero();
        RHS.setZero();
        
        if (m_parallel_assembly && m_class_reuse) {
            // local operators once per geometric class, local vectors per cell
            auto& cl = cell_classes(msh);
            auto e_ops = e_class_operators(msh);
            tools::ordered_loop(true, cl.e_ids.size(), [&](size_t i) {
                auto& cell = storage->surfaces[cl.e_ids[i]];
                return Matrix<T, Dynamic, 1>(make_rhs(msh, cell, make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree()), e_rhs_fun));
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
                scatter_e_data(m_e_cell_index[cl.e_ids[i]], msh, storage->surfaces[cl.e_ids[i]], e_ops[cl.e_class_of[i]], f_loc);
            });
            auto a_ops = a_class_operators(msh);
            tools::ordered_loop(true, cl.a_ids.size(), [&](size_t i) {
                auto& cell = storage->surfaces[cl.a_ids[i]];
                return Matrix<T, Dynamic, 1>(make_rhs(msh, cell, make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree()), a_rhs_fun));
            }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
                scatter_a_data(m_a_cell_index[cl.a_ids[i]], msh, storage->surfaces[cl.a_ids[i]], a_ops[cl.a_class_of[i]], f_loc);
            });
            finalize_parallel();
            return;
        }
        if (m_parallel_assembly) {
            typedef std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, 1>> local_type;
            std::vector<std::pair<size_t, elastic_material_data<T>>> e_cells(m_e_material.begin(), m_e_material.end());
            tools::ordered_loop(true, e_cells.size(), [&](size_t i) {
                auto& cell = storage->surfaces[e_cells[i].first];
                elastic_material_data<T> material = e_cells[i].second;
                auto cell_basis = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
                return local_type(e_laplacian_operator(material, msh, cell), make_rhs(msh, cell, cell_basis, e_rhs_fun));
            }, [&](size_t i, const local_type& loc) {
                scatter_e_data(m_e_cell_index[e_cells[i].first], msh, storage->surfaces[e_cells[i].first], loc.first, loc.second);
            });
            std::vector<std::pair<size_t, acoustic_material_data<T>>> a_cells(m_a_material.begin(), m_a_material.end());
            tools::ordered_loop(true, a_cells.size(), [&](size_t i) {
                auto& cell = storage->surfaces[a_cells[i].first];
                acoustic_material_data<T> material = a_cells[i].second;
                auto cell_basis = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                return local_type(a_laplacian_operator(material, msh, cell), make_rhs(msh, cell, cell_basis, a_rhs_fun));
            }, [&](size_t i, const local_type& loc) {
                scatter_a_data(m_a_cell_index[a_cells[i].first], msh, storage->surfaces[a_cells[i].first], loc.first, loc.second);
            });
            finalize();
            return;
        }

        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> vectorial_laplacian_operator_loc = e_laplacian_operator(e_chunk.second,msh,cell);
            auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, e_rhs_fun);
            scatter_e_data(e_cell_ind, msh, cell, vectorial_laplacian_operator_loc, f_loc);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            
            Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = a_laplacian_operator(a_chunk.second, msh, cell);
            auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
            Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, a_rhs_fun);
            scatter_a_data(a_cell_ind, msh, cell, laplacian_operator_loc, f_loc);
        }
        finalize();
    }
    
    void apply_bc_conditions_on_interface(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_vel_fun, std::function<T(const typename Mesh::point_type& )> a_vel_fun){
        
        auto storage = msh.backend_storage();
        // Applying transmission conditions as neumann data for the case of uncoupled problems
        for (auto chunk : m_interface_cell_indexes) {
            auto& face = storage->edges[chunk.first];
            auto& e_cell = storage->surfaces[chunk.second.first];
            auto& a_cell = storage->surfaces[chunk.second.second];
            
            Matrix<T, Dynamic, 1> e_neumann_operator_loc = e_neumman_bc_operator(msh, face, e_cell, a_cell, a_vel_fun);
            scatter_e_interface_neumann_data(msh, face, e_neumann_operator_loc);
            
            Matrix<T, Dynamic, 1> a_neumann_operator_loc = a_neumman_bc_operator(msh, face, e_cell, a_cell, e_vel_fun);
            scatter_a_interface_neumann_data(msh, face, a_neumann_operator_loc);
            
        }
    }

    void assemble_mass(const Mesh& msh){
        
        auto storage = msh.backend_storage();
        MASS.setZero();
        
        if (m_parallel_assembly && m_class_reuse) {
            // mass matrices once per geometric class
            auto& cl = cell_classes(msh);
            auto e_ops = tools::class_operators(cl.e_rep, [&](size_t r) {
                elastic_material_data<T> material = m_e_material.at(cl.e_ids[r]);
                return e_mass_operator(material, msh, storage->surfaces[cl.e_ids[r]]);
            });
            for (size_t i = 0; i < cl.e_ids.size(); i++)
                scatter_e_mass_data(m_e_cell_index[cl.e_ids[i]], msh, storage->surfaces[cl.e_ids[i]], e_ops[cl.e_class_of[i]]);
            auto a_ops = tools::class_operators(cl.a_rep, [&](size_t r) {
                acoustic_material_data<T> material = m_a_material.at(cl.a_ids[r]);
                return a_mass_operator(material, msh, storage->surfaces[cl.a_ids[r]]);
            });
            for (size_t i = 0; i < cl.a_ids.size(); i++)
                scatter_a_mass_data(m_a_cell_index[cl.a_ids[i]], msh, storage->surfaces[cl.a_ids[i]], a_ops[cl.a_class_of[i]]);
            finalize_mass_parallel();
            return;
        }
        if (m_parallel_assembly) {
            std::vector<std::pair<size_t, elastic_material_data<T>>> e_cells(m_e_material.begin(), m_e_material.end());
            tools::ordered_loop(true, e_cells.size(), [&](size_t i) {
                elastic_material_data<T> material = e_cells[i].second;
                return e_mass_operator(material, msh, storage->surfaces[e_cells[i].first]);
            }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
                scatter_e_mass_data(m_e_cell_index[e_cells[i].first], msh, storage->surfaces[e_cells[i].first], mass_matrix);
            });
            std::vector<std::pair<size_t, acoustic_material_data<T>>> a_cells(m_a_material.begin(), m_a_material.end());
            tools::ordered_loop(true, a_cells.size(), [&](size_t i) {
                acoustic_material_data<T> material = a_cells[i].second;
                return a_mass_operator(material, msh, storage->surfaces[a_cells[i].first]);
            }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
                scatter_a_mass_data(m_a_cell_index[a_cells[i].first], msh, storage->surfaces[a_cells[i].first], mass_matrix);
            });
            finalize_mass();
            return;
        }

        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, Dynamic> mass_matrix = e_mass_operator(e_chunk.second,msh, cell);
            scatter_e_mass_data(e_cell_ind,msh, cell, mass_matrix);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, Dynamic> mass_matrix = a_mass_operator(a_chunk.second,msh, cell);
            scatter_a_mass_data(a_cell_ind,msh, cell, mass_matrix);
        }
        
        finalize_mass();
    }
    
    void assemble_rhs(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun){
        
        RHS.setZero();
        if (m_parallel_assembly || m_rc) {
            assemble_rhs_fast(msh, e_rhs_fun, a_rhs_fun);
            return;
        }
        #ifdef HAVE_INTEL_TBB2
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&msh,&rhs_fun] (size_t & cell_ind){
                        auto& cell = msh.backend_storage()->surfaces[cell_ind];
                        auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
                        Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
                        this->scatter_rhs_data(msh, cell, f_loc);
                }
            );
        #else
            auto storage = msh.backend_storage();
             for (auto e_chunk : m_e_material) {
                 size_t e_cell_ind = m_e_cell_index[e_chunk.first];
                 auto& cell = storage->surfaces[e_chunk.first];
                 auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
                 Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, e_rhs_fun);
                 scatter_e_rhs_data(e_cell_ind, msh, cell, f_loc);
             }
        
            for (auto a_chunk : m_a_material) {
                size_t a_cell_ind = m_a_cell_index[a_chunk.first];
                auto& cell = storage->surfaces[a_chunk.first];
                auto cell_basis   = make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
                Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, a_rhs_fun);
                scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
            }
        
        #endif
        apply_bc(msh);
    }
    
    void apply_bc(const Mesh& msh){
        
        if (m_rc) {
            apply_bc_cached(msh);
            return;
        }
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
            for (auto& cell_ind : m_e_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> vectorial_laplacian_operator_loc = e_laplacian_operator(e_mat,msh,cell);
                scatter_e_bc_data(m_e_cell_index[cell_ind], msh, cell, vectorial_laplacian_operator_loc);
            }
        
            for (auto& cell_ind : m_a_elements_with_bc_eges)
            {
                auto& cell = storage->surfaces[cell_ind];
                acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = a_laplacian_operator(a_mat, msh, cell);
                scatter_a_bc_data(m_a_cell_index[cell_ind], msh, cell, laplacian_operator_loc);
            }
        
        #endif
        
    }
    
    // assemble_rhs of the fast mode: local vectors in parallel (m_parallel_assembly) on the cached quadrature data (m_rc),
    // scattered in the order of the serial loop
    void assemble_rhs_fast(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )>& e_rhs_fun,
                           std::function<T(const typename Mesh::point_type& )>& a_rhs_fun) {
        auto storage = msh.backend_storage();
        const auto celdeg = m_hho_di.cell_degree();
        const size_t vec_bs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
        const size_t sca_bs = disk::scalar_basis_size(celdeg, Mesh::dimension);
        std::vector<size_t> e_keys, a_keys;
        for (auto& e_chunk : m_e_material) e_keys.push_back(e_chunk.first);
        for (auto& a_chunk : m_a_material) a_keys.push_back(a_chunk.first);
        bool cached = false;
        if (m_rc && !m_rc->rhs_off) {
            if (!m_rc->rhs_ready) {
                size_t doubles = tools::quadrature_doubles(msh, e_keys, celdeg, 2 * vec_bs) + tools::quadrature_doubles(msh, a_keys, celdeg, sca_bs);
                if (doubles > m_rhs_cache_limit)
                    m_rc->rhs_off = true;
                else {
                    tools::cell_quadratures(msh, e_keys, celdeg, [&](const typename Mesh::cell_type& cell) { return make_vector_monomial_basis(msh, cell, celdeg); },
                                            m_rc->e_rhs, m_parallel_assembly);
                    tools::cell_quadratures(msh, a_keys, celdeg, [&](const typename Mesh::cell_type& cell) { return make_scalar_monomial_basis(msh, cell, celdeg); },
                                            m_rc->a_rhs, m_parallel_assembly);
                    m_rc->rhs_ready = true;
                }
            }
            cached = m_rc->rhs_ready;
        }
        tools::ordered_loop(m_parallel_assembly, e_keys.size(), [&](size_t i) {
            if (cached) return tools::rhs(m_rc->e_rhs[i], vec_bs, e_rhs_fun);
            auto& cell = storage->surfaces[e_keys[i]];
            return Matrix<T, Dynamic, 1>(make_rhs(msh, cell, make_vector_monomial_basis(msh, cell, celdeg), e_rhs_fun));
        }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
            scatter_e_rhs_data(m_e_cell_index[e_keys[i]], msh, storage->surfaces[e_keys[i]], f_loc);
        });
        tools::ordered_loop(m_parallel_assembly, a_keys.size(), [&](size_t i) {
            if (cached) return tools::rhs(m_rc->a_rhs[i], sca_bs, a_rhs_fun);
            auto& cell = storage->surfaces[a_keys[i]];
            return Matrix<T, Dynamic, 1>(make_rhs(msh, cell, make_scalar_monomial_basis(msh, cell, celdeg), a_rhs_fun));
        }, [&](size_t i, const Matrix<T, Dynamic, 1>& f_loc) {
            scatter_a_rhs_data(m_a_cell_index[a_keys[i]], msh, storage->surfaces[a_keys[i]], f_loc);
        });
        apply_bc(msh);
    }

    // apply_bc with the cached local operators and Dirichlet faces of the cells with boundary faces
    void apply_bc_cached(const Mesh& msh) {
        auto storage = msh.backend_storage();
        const size_t fdeg = m_hho_di.face_degree();
        const size_t n_e_fbs = disk::vector_basis_size(fdeg, Mesh::dimension - 1, Mesh::dimension);
        const size_t n_a_fbs = disk::scalar_basis_size(fdeg, Mesh::dimension - 1);
        if (!m_rc->bc_ready) {
            if (m_class_reuse) {
                // local operators of the classes of the cells
                auto& cl = cell_classes(msh);
                auto e_ops = e_class_operators(msh);
                auto a_ops = a_class_operators(msh);
                m_rc->e_bc_ops.clear();
                m_rc->a_bc_ops.clear();
                for (size_t cell_ind : m_e_elements_with_bc_eges) m_rc->e_bc_ops.push_back(e_ops[cl.e_class_of[cl.e_pos.at(cell_ind)]]);
                for (size_t cell_ind : m_a_elements_with_bc_eges) m_rc->a_bc_ops.push_back(a_ops[cl.a_class_of[cl.a_pos.at(cell_ind)]]);
            }
            else {
                m_rc->e_bc_ops.assign(m_e_elements_with_bc_eges.size(), Matrix<T, Dynamic, Dynamic>());
                m_rc->a_bc_ops.assign(m_a_elements_with_bc_eges.size(), Matrix<T, Dynamic, Dynamic>());
                tools::ordered_loop(m_parallel_assembly, m_e_elements_with_bc_eges.size(), [&](size_t i) {
                    const size_t cell_ind = m_e_elements_with_bc_eges[i];
                    elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
                    return e_laplacian_operator(e_mat, msh, storage->surfaces[cell_ind]);
                }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& op) { m_rc->e_bc_ops[i] = op; });
                tools::ordered_loop(m_parallel_assembly, m_a_elements_with_bc_eges.size(), [&](size_t i) {
                    const size_t cell_ind = m_a_elements_with_bc_eges[i];
                    acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
                    return a_laplacian_operator(a_mat, msh, storage->surfaces[cell_ind]);
                }, [&](size_t i, const Matrix<T, Dynamic, Dynamic>& op) { m_rc->a_bc_ops[i] = op; });
            }
            const size_t faces0 = m_n_elastic_cell_dof + m_n_acoustic_cell_dof;
            tools::prepare_bc_cells(msh, m_rc->e_bc_cells, m_e_elements_with_bc_eges, m_e_bnd, m_e_cell_index, m_e_compress_indexes,
                get_e_cell_basis_data(), n_e_fbs, 0, faces0,
                [&](const typename Mesh::face_type& fc) { return make_vector_monomial_basis(msh, fc, fdeg); }, m_parallel_assembly);
            tools::prepare_bc_cells(msh, m_rc->a_bc_cells, m_a_elements_with_bc_eges, m_a_bnd, m_a_cell_index, m_a_compress_indexes,
                get_a_cell_basis_data(), n_a_fbs, m_n_elastic_cell_dof, faces0 + m_n_elastic_face_dof,
                [&](const typename Mesh::face_type& fc) { return make_scalar_monomial_basis(msh, fc, fdeg); }, m_parallel_assembly);
            m_rc->bc_ready = true;
        }
        tools::ordered_loop(m_parallel_assembly, m_rc->e_bc_cells.size(), [&](size_t i) {
            return tools::bc_dirichlet_data(m_rc->e_bc_cells[i], m_e_bnd, get_e_cell_basis_data(), n_e_fbs);
        }, [&](size_t i, const Matrix<T, Dynamic, 1>& data) { tools::scatter_bc_cached(RHS, m_rc->e_bc_cells[i].asm_map, m_rc->e_bc_ops[i], data); });
        tools::ordered_loop(m_parallel_assembly, m_rc->a_bc_cells.size(), [&](size_t i) {
            return tools::bc_dirichlet_data(m_rc->a_bc_cells[i], m_a_bnd, get_a_cell_basis_data(), n_a_fbs);
        }, [&](size_t i, const Matrix<T, Dynamic, 1>& data) { tools::scatter_bc_cached(RHS, m_rc->a_bc_cells[i].asm_map, m_rc->a_bc_ops[i], data); });
    }

    Matrix<T, Dynamic, Dynamic> e_mass_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
        auto vec_basis = disk::make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
        Matrix<T, Dynamic, Dynamic> mass_matrix = disk::make_mass_matrix(msh, cell, vec_basis);
        mass_matrix *= (material.rho());
        return mass_matrix;
    }
    
    Matrix<T, Dynamic, Dynamic> a_mass_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){

        auto scal_basis = disk::make_scalar_monomial_basis(msh, cell, m_hho_di.cell_degree());
        Matrix<T, Dynamic, Dynamic> mass_matrix = disk::make_mass_matrix(msh, cell, scal_basis);
        mass_matrix *= (1.0/(material.rho()*material.vp()*material.vp()));
        return mass_matrix;
    }
    
    
    Matrix<T, Dynamic, Dynamic> e_laplacian_operator(elastic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){
           
        T mu = material.rho()*material.vs()*material.vs();
        T lambda = material.rho()*material.vp()*material.vp() - 2.0*mu;
        auto reconstruction_operator   = make_matrix_hho_symmetric_gradrec(msh, cell, m_hho_di);
        auto rec_for_stab   = make_vector_hho_symmetric_laplacian(msh, cell, m_hho_di);
        auto divergence_operator = make_hho_divergence_reconstruction(msh, cell, m_hho_di);

        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        Matrix<T, Dynamic, Dynamic> S_operator;
        if(m_hho_stabilization_Q)
        {
            auto stabilization_operator    = make_vector_hho_stabilization(msh, cell, rec_for_stab.first, m_hho_di, m_stab_size != 1);
            S_operator = stabilization_operator;
        }else{
            auto stabilization_operator    = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_stab_size != 1, m_scaled_stab);
            S_operator = stabilization_operator;
        }
        return 2.0 * mu * (R_operator + m_stab_weight * S_operator) + lambda * divergence_operator.second;
    }
    
    Matrix<T, Dynamic, Dynamic> a_laplacian_operator(acoustic_material_data<T> & material, const Mesh& msh, const typename Mesh::cell_type& cell){

        auto reconstruction_operator   = make_scalar_hho_laplacian(msh, cell, m_hho_di);
        Matrix<T, Dynamic, Dynamic> R_operator = reconstruction_operator.second;
        
        Matrix<T, Dynamic, Dynamic> S_operator;
        if(m_hho_stabilization_Q)
        {
            // diskpp: 1/h_T with hho_degree_info, the cell degree information variant scales by 1/h_F
            Matrix<T, Dynamic, Dynamic> stabilization_operator = m_stab_size == 2
                ? make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, disk::CellDegreeInfo<Mesh>(msh, cell, m_hho_di.cell_degree(), m_hho_di.face_degree()), true)
                : make_scalar_hho_stabilization(msh, cell, reconstruction_operator.first, m_hho_di);
            S_operator = stabilization_operator;
        }else{
            auto stabilization_operator    = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_stab_size != 1, m_scaled_stab);
            S_operator = stabilization_operator;
        }
        return (1.0/material.rho())*(R_operator + m_stab_weight * S_operator);
    }
    
    
    Matrix<T, Dynamic, 1> e_neumman_bc_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell, std::function<T(const typename Mesh::point_type& )> a_vel_fun){

        Matrix<T, Dynamic, Dynamic> neumann_operator;
        auto facdeg = m_hho_di.face_degree();
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

        neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(vfbs, 1);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        const auto qps = integrate(msh, face, 2 * facdeg);
        const auto n = disk::normal(msh, e_cell, face);
        for (auto& qp : qps)
        {
            const auto v_f_phi = vfb.eval_functions(qp.point());
            const auto s_f_phi = sfb.eval_functions(qp.point());

            assert(v_f_phi.rows() == vfbs);
            auto q_n = disk::priv::inner_product(qp.weight(), n);
            const auto n_dot_v_f_phi = disk::priv::inner_product(v_f_phi,q_n);
            const auto result = (disk::priv::outer_product(n_dot_v_f_phi, a_vel_fun(qp.point()))).eval();
            neumann_operator += result;
            
        }
        return neumann_operator;
        }
    
    Matrix<T, Dynamic, Dynamic> a_neumman_bc_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_vel_fun){

        Matrix<T, Dynamic, Dynamic> neumann_operator;
        auto facdeg = m_hho_di.face_degree();
        auto sfbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(sfbs, 1);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        const auto qps = integrate(msh, face, 2 * facdeg);
        const auto n = disk::normal(msh, a_cell, face);
        for (auto& qp : qps)
        {
            const auto v_f_phi = vfb.eval_functions(qp.point());
            const auto s_f_phi = sfb.eval_functions(qp.point());

            assert(s_f_phi.rows() == sfbs);
            auto q_n = disk::priv::inner_product(qp.weight(), n);
            const auto n_dot_v_f = disk::priv::inner_product(e_vel_fun(qp.point()),q_n);
            const auto result = (disk::priv::inner_product(n_dot_v_f, s_f_phi)).eval();
            neumann_operator += result;
        }
        return neumann_operator;
    }
    
    void project_over_cells(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> vec_fun, std::function<T(const typename Mesh::point_type& )> scal_fun){
        
        auto storage = msh.backend_storage();
        size_t n_dof = MASS.rows();
        x_glob = Matrix<T, Dynamic, 1>::Zero(n_dof);

        // elastic block
        for (auto e_chunk : m_e_material) {
            size_t e_cell_ind = m_e_cell_index[e_chunk.first];
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cell, m_hho_di.cell_degree(), vec_fun);
            scatter_e_cell_dof_data(e_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            size_t a_cell_ind = m_a_cell_index[a_chunk.first];
            auto& cell = storage->surfaces[a_chunk.first];
            Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cell, m_hho_di.cell_degree(), scal_fun);
            scatter_a_cell_dof_data(a_cell_ind, msh, cell, x_glob, x_proj_dof);
        }
    
    }

private:

    T m_stab_weight = 1.0;
    bool m_scaled_stab = true;
    int  m_stab_size = 0;           // length of the 1/h scaling, see set_stabilization_size

    // finalize, finalize_mass of the fast mode: the same matrices, bit for bit, built by the OpenMP threads
    void finalize_parallel() {
        tools::set_from_triplets(LHS, m_triplets);
        std::vector<Triplet<T>>().swap(m_triplets);
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
    std::vector<Matrix<T, Dynamic, Dynamic>> e_class_operators(const Mesh& msh) {
        auto& cl = cell_classes(msh);
        auto storage = msh.backend_storage();
        return tools::class_operators(cl.e_rep, [&](size_t r) {
            elastic_material_data<T> material = m_e_material.at(cl.e_ids[r]);
            return e_laplacian_operator(material, msh, storage->surfaces[cl.e_ids[r]]);
        });
    }
    std::vector<Matrix<T, Dynamic, Dynamic>> a_class_operators(const Mesh& msh) {
        auto& cl = cell_classes(msh);
        auto storage = msh.backend_storage();
        return tools::class_operators(cl.a_rep, [&](size_t r) {
            acoustic_material_data<T> material = m_a_material.at(cl.a_ids[r]);
            return a_laplacian_operator(material, msh, storage->surfaces[cl.a_ids[r]]);
        });
    }

    // fast mode (set_parallel_assembly, set_cache_rhs_data): local operators of the cells with Dirichlet faces,
    // quadrature data of the source terms (indexed as m_e_material, m_a_material; not kept beyond m_rhs_cache_limit
    // doubles), assembly indices and Dirichlet faces of the cells with boundary faces
    bool m_parallel_assembly = false;
    struct rhs_cache_data {
        std::vector<Matrix<T, Dynamic, Dynamic>> e_bc_ops, a_bc_ops;
        bool bc_ready = false;
        typename tools::e_rhs_data e_rhs;
        typename tools::a_rhs_data a_rhs;
        bool rhs_ready = false, rhs_off = false;
        typename tools::e_bc_data e_bc_cells;
        typename tools::a_bc_data a_bc_cells;
    };
    static constexpr size_t m_rhs_cache_limit = size_t(1) << 27;    // doubles (1 GB)
    std::unique_ptr<rhs_cache_data> m_rc;
    bool m_class_reuse = false;
    std::unique_ptr<typename tools::classes> m_cls;
};

#endif /* elastoacoustic_two_fields_assembler_hpp */
