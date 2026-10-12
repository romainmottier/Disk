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
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/elastic_material_data.hpp"
#include "../common/elastoacoustic_assembler_base.hpp"
#include <map>

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
    }

    T stabilization_weight() const {
        return m_stab_weight;
    }

    // Least-squares stabilization weighted by 1/h (default) or O(1); the HHO stabilization is always weighted by 1/h
    void set_stabilization_scaling(bool scaled) {
        m_scaled_stab = scaled;
    }

    // Length h of the 1/h scaling: 0 as in diskpp (face diameter h_F, except for the scalar HHO stabilization: cell
    // diameter h_T), 1 cell diameter h_T, 2 face diameter h_F
    void set_stabilization_size(int size) {
        m_stab_size = size;
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
    
};

#endif /* elastoacoustic_two_fields_assembler_hpp */
