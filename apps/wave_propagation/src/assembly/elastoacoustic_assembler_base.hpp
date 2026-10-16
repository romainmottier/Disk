//
//  elastoacoustic_assembler_base.hpp
//  acoustics
//
//  Degrees of freedom bookkeeping shared by the elasto-acoustic assemblers of the first-order
//  (elastoacoustic_four_fields_assembler) and second-order (elastoacoustic_two_fields_assembler) formulations.
//  Global layout: elastic cells, acoustic cells, elastic faces, acoustic faces. Each elastic (resp. acoustic)
//  cell owns n_e_cbs (resp. n_a_cbs) unknowns, set by the formulation; elastic faces carry vector polynomials
//  and acoustic faces scalar ones.
//

#pragma once
#ifndef elastoacoustic_assembler_base_hpp
#define elastoacoustic_assembler_base_hpp

#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "assembly_index.hpp"
#include "../physics/acoustic_material_data.hpp"
#include "../physics/elastic_material_data.hpp"
#include <map>
#include <set>

template<typename Mesh>
class elastoacoustic_assembler_base {

protected:

    typedef disk::BoundaryConditions<Mesh, false>    e_boundary_type;
    typedef disk::BoundaryConditions<Mesh, true>     a_boundary_type;
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>                 m_e_compress_indexes;
    std::vector<size_t>                 m_e_expand_indexes;
    
    std::vector<size_t>                 m_a_compress_indexes;
    std::vector<size_t>                 m_a_expand_indexes;

    disk::hho_degree_info               m_hho_di;
    e_boundary_type                     m_e_bnd;
    a_boundary_type                     m_a_bnd;
    std::vector< Triplet<T> >           m_triplets;
    std::vector< Triplet<T> >           m_c_triplets;
    std::vector< Triplet<T> >           m_mass_triplets;
    std::map<size_t,elastic_material_data<T>> m_e_material;
    std::map<size_t,acoustic_material_data<T>> m_a_material;
    std::map<size_t,size_t> m_e_cell_index;
    std::map<size_t,size_t> m_a_cell_index;
    std::vector< size_t >               m_e_elements_with_bc_eges;
    std::vector< size_t >               m_a_elements_with_bc_eges;
    std::map<size_t,std::pair<size_t,size_t>>   m_interface_cell_indexes;

    size_t      m_n_elastic_cells;
    size_t      m_n_acoustic_cells;

    size_t      n_e_edges;
    size_t      n_a_edges;

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    size_t      m_n_e_cbs;                      // cell unknowns per elastic cell (set by the formulation)
    size_t      m_n_a_cbs;                      // cell unknowns per acoustic cell (set by the formulation)
    size_t      m_n_elastic_cell_dof;
    size_t      m_n_acoustic_cell_dof;
    size_t      m_n_elastic_face_dof;
    size_t      m_n_acoustic_face_dof;
    bool        m_hho_stabilization_Q;
    bool        m_legacy_interface_quadrature = false;  // coupling integrated with the degree-k rule of the papers' code

    // Identification of Dirichlet faces; Construction of compressed face index maps; Computation of dofs counts; 
    // Initialization of global matrices and RHS; Setup of cell classification and local-to-global mapping structures.
    elastoacoustic_assembler_base(const Mesh& msh, const disk::hho_degree_info& hho_di, const e_boundary_type& e_bnd, const a_boundary_type& a_bnd, std::map<size_t,elastic_material_data<T>> & e_material, std::map<size_t,acoustic_material_data<T>> & a_material, size_t n_e_cbs, size_t n_a_cbs) : m_hho_di(hho_di), m_e_bnd(e_bnd), m_a_bnd(a_bnd), m_e_material(e_material), m_a_material(a_material), m_n_e_cbs(n_e_cbs), m_n_a_cbs(n_a_cbs), m_hho_stabilization_Q(true) {
            
        auto storage = msh.backend_storage();
        auto is_e_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return e_bnd.is_dirichlet_face(fc_id);
        };
        
        auto is_a_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return a_bnd.is_dirichlet_face(fc_id);
        };

        size_t n_e_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_e_dirichlet);
        size_t n_a_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_a_dirichlet);
        
        std::set<size_t> e_egdes;
        for (auto &chunk : m_e_material) {
            size_t cell_i = chunk.first;
            auto& cell = storage->surfaces[cell_i];
            auto cell_faces = faces(msh,cell);
            for (auto &face : cell_faces) {
                if (!is_e_dirichlet(face)) {
                    auto fc_id = msh.lookup(face);
                    e_egdes.insert(fc_id);
                }
            }
        }
        n_e_edges = e_egdes.size();
        std::set<size_t> a_egdes;
        for (auto &chunk : m_a_material) {
            size_t cell_i = chunk.first;
            auto& cell = storage->surfaces[cell_i];
            auto cell_faces = faces(msh,cell);
            for (auto &face : cell_faces) {
                if (!is_a_dirichlet(face)) {
                    auto fc_id = msh.lookup(face);
                    a_egdes.insert(fc_id);
                }
            }
        }
        n_a_edges = a_egdes.size();
        
        m_n_edges = msh.faces_size();
        m_n_essential_edges = n_e_essential_edges + n_a_essential_edges;

        m_e_compress_indexes.resize( m_n_edges );
        m_e_expand_indexes.resize( m_n_edges - m_n_essential_edges );
        
        m_a_compress_indexes.resize( m_n_edges );
        m_a_expand_indexes.resize( m_n_edges - m_n_essential_edges );

        size_t e_compressed_offset = 0;
        for (auto face_id : e_egdes) {
            m_e_compress_indexes.at(face_id) = e_compressed_offset;
            m_e_expand_indexes.at(e_compressed_offset) = face_id;
            e_compressed_offset++;
        }
        size_t a_compressed_offset = 0;
        for (auto face_id : a_egdes) {
            m_a_compress_indexes.at(face_id) = a_compressed_offset;
            m_a_expand_indexes.at(a_compressed_offset) = face_id;
            a_compressed_offset++;
        }
    
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        size_t n_s_cbs = get_a_cell_basis_data();
        size_t n_s_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        m_n_elastic_cell_dof = (n_cbs * m_e_material.size());
        m_n_acoustic_cell_dof = (n_s_cbs * m_a_material.size());
        
        m_n_elastic_face_dof = (n_fbs * n_e_edges);
        m_n_acoustic_face_dof = (n_s_fbs * n_a_edges);
        size_t system_size = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_n_acoustic_face_dof;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
        MASS = SparseMatrix<T>( system_size, system_size );
        COUPLING = SparseMatrix<T>( system_size, system_size );

        classify_cells(msh);
        build_cells_maps();
    }

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;
    SparseMatrix<T>         MASS;
    SparseMatrix<T>         COUPLING;
    
    // Elastic cell– and face–dofs mapping; construction of assembly indices with Dirichlet filtering; Computation of global offsets for elastic unknowns; 
    // Insertion of local matrix entries into global triplets; Insertion of local RHS contributions into the global vector.
    void scatter_e_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {

        auto fcs = faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset = e_cell_ind * n_cbs;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }
            
        for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
            bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
            for (size_t i = 0; i < n_fbs; i++) {
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
        }
            
        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < lhs.cols(); j++) {
                if ( asm_map[j].assemble() ) {
                    if (m_skip_zero_entries && lhs(i,j) == T(0)) continue;
                    m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                }
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            RHS(asm_map[i].vidx()) += rhs(i);
        }
    }
    
    // Acoustic cell– and face–dofs mapping; Construction of assembly indices with Dirichlet filtering; Computation of global offsets for acoustic unknowns; 
    // Insertion of local matrix entries into global triplets; Insertion of local RHS contributions into the global vector.
    void scatter_a_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {

        auto fcs = faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }
        
        for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
            bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
            for (size_t i = 0; i < n_fbs; i++) {
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) {
                    if (m_skip_zero_entries && lhs(i,j) == T(0)) continue;
                    m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                }
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            RHS(asm_map[i].vidx()) += rhs(i);
        }
        
    }
    
    // Elastic cell–DoF mapping; Construction of assembly indices; Computation of global offsets for elastic cell unknowns; 
    // Insertion of local mass-matrix entries into global mass triplets.
    void scatter_e_mass_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
        
        size_t n_cbs = get_e_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);        
        auto cell_LHS_offset = e_cell_ind * n_cbs;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back(assembly_index(cell_LHS_offset+i, true));
        }

        assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );
        for (size_t i = 0; i < mass_matrix.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < mass_matrix.cols(); j++) {
                if (asm_map[j].assemble()) {
                    if (m_skip_zero_entries && mass_matrix(i,j) == T(0)) continue;
                    m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
                }
            }
        }
    }
    
    // Acoustic cell–DoF mapping; Construction of assembly indices; Computation of global offsets for acoustic cell unknowns; 
    // Insertion of local mass-matrix entries into global mass triplets.
    void scatter_a_mass_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& mass_matrix) {

        size_t n_cbs = get_a_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);
        auto cell_LHS_offset = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back(assembly_index(cell_LHS_offset+i, true));
        }

        assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );
        for (size_t i = 0; i < mass_matrix.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < mass_matrix.cols(); j++) {
                if (asm_map[j].assemble()) {
                    if (m_skip_zero_entries && mass_matrix(i,j) == T(0)) continue;
                    m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
                }
            }
        }
    }
    
    // Elastic cell–DoF mapping; Construction of assembly indices; Computation of global offsets for elastic cell unknowns; 
    // Insertion of local RHS contributions into the global vector.
    void scatter_e_rhs_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {
    
        size_t n_cbs = get_e_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = e_cell_ind * n_cbs;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }

        assert( asm_map.size() == rhs.rows() );
        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            RHS(asm_map[i].vidx()) += rhs(i);
        }
    }
    
    // Acoustic cell–DoF mapping; Construction of assembly indices; Computation of global offsets for acoustic cell unknowns; 
    // Insertion of local RHS contributions into the global vector.
    void scatter_a_rhs_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {

        size_t n_cbs = get_a_cell_basis_data();
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs);

        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }

        assert( asm_map.size() == rhs.rows());
        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            RHS(asm_map[i].vidx()) += rhs(i);
        }
    }
    
    // Elastic Dirichlet face identification; Construction of assembly indices for cell and face DoFs; 
    // Computation of Dirichlet face contributions via local mass matrices; Insertion of bc into the global RHS using local-to-global mapping.
    void scatter_e_bc_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs) {
    
        auto fcs = faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());

        auto cell_LHS_offset = e_cell_ind * n_cbs;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
            bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
            for (size_t i = 0; i < n_fbs; i++) {
                asm_map.push_back(assembly_index(face_LHS_offset+i, !dirichlet) );
            }
            if (dirichlet) {
                 auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                 auto dirichlet_fun  = m_e_bnd.dirichlet_boundary_func(fc_id);
                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                 Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                 dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
             }
        }

        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (!asm_map[j].assemble() ) {
                    RHS(asm_map[i].vidx()) -= lhs(i,j) * dirichlet_data(j);
                }
            }
        }
    }
    
    // Acoustic Dirichlet face identification; construction of assembly indices for cell and face DoFs; 
    // Computation of Dirichlet face contributions via local mass matrices; Insertion of bc into the global RHS using local-to-global mapping.
    void scatter_a_bc_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs) {

        auto fcs = faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        std::vector<assembly_index> asm_map;
        asm_map.reserve(n_cbs + n_fbs*fcs.size());
        
        auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        for (size_t i = 0; i < n_cbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        }
        
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
        for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
            auto fc = fcs[face_i];
            auto fc_id = msh.lookup(fc);
            auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;            
            bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
            for (size_t i = 0; i < n_fbs; i++) {
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
            if (dirichlet) {
                 auto fb = make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                 auto dirichlet_fun  = m_a_bnd.dirichlet_boundary_func(fc_id);
                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                 Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                 dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
             }
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble()) {
                continue;
            }
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (!asm_map[j].assemble()) {
                    RHS(asm_map[i].vidx()) -= lhs(i,j) * dirichlet_data(j);
                }
            }
        }
    }
    
    // Elasto–acoustic interface face identification; Construction of Interface matrix.
    void scatter_ea_interface_data(const Mesh& msh, const typename Mesh::face_type& face, const Matrix<T, Dynamic, Dynamic>& interface_matrix) {
        
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);

        std::vector<assembly_index> asm_map_e, asm_map_a;
        asm_map_e.reserve(vfbs);
        asm_map_a.reserve(sfbs);
        auto fc_id = msh.lookup(face);
        
        auto e_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*vfbs;
        bool e_dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < vfbs; i++){
            asm_map_e.push_back( assembly_index(e_face_LHS_offset+i, !e_dirichlet) );
        }
            
        auto a_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*sfbs;
        bool a_dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
        for (size_t i = 0; i < sfbs; i++){
            asm_map_a.push_back( assembly_index(a_face_LHS_offset+i, !a_dirichlet) );
        }
        
        assert( asm_map_e.size() == interface_matrix.rows() && asm_map_a.size() == interface_matrix.cols() );
        for (size_t i = 0; i < interface_matrix.rows(); i++) {
            for (size_t j = 0; j < interface_matrix.cols(); j++) {
              m_c_triplets.push_back( Triplet<T>(asm_map_e[i].vidx(), asm_map_a[j].vidx(),   interface_matrix(i,j)) );
              m_c_triplets.push_back( Triplet<T>(asm_map_a[j].vidx(), asm_map_e[i].vidx(), - interface_matrix(i,j)) );
            }
        }

    }
    
    // Assembly of global coupling matrix; Local elasto–acoustic interface operator computation for each interface face; 
    // Scattering of interface contributions into global coupling triplets; finalization of coupling assembly.
    void assemble_coupling_terms(const Mesh& msh){
        
        auto storage = msh.backend_storage();
        COUPLING.setZero();

        // coupling blocks
        for (auto chunk : m_interface_cell_indexes) {
            auto& face = storage->edges[chunk.first];
            auto& e_cell = storage->surfaces[chunk.second.first];
            auto& a_cell = storage->surfaces[chunk.second.second];
            Matrix<T, Dynamic, Dynamic> interface_operator_loc = e_interface_operator(msh, face, e_cell, a_cell);
            scatter_ea_interface_data(msh, face, interface_operator_loc);
        }

        finalize_coupling();
    }
    
    // Computation of local elasto–acoustic interface operator; Quadrature over the face to integrate normal component of elastic basis with acoustic basis; 
    Matrix<T, Dynamic, Dynamic> 
    e_interface_operator(const Mesh& msh, const typename Mesh::face_type& face, const typename Mesh::cell_type& e_cell, const typename Mesh::cell_type& a_cell) {

        Matrix<T, Dynamic, Dynamic> interface_operator;
        auto facdeg = m_hho_di.face_degree();
        auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto sfbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        interface_operator = Matrix<T, Dynamic, Dynamic>::Zero(vfbs, sfbs);
        
        auto vfb = make_vector_monomial_basis(msh, face, facdeg);
        auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
        // the integrand is a product of two face polynomials of degree k: a rule of degree k (k/2 + 1 Gauss points)
        // gave a rank-deficient operator for k >= 1
        const auto qps = integrate(msh, face, m_legacy_interface_quadrature ? facdeg : 2 * facdeg);
        const auto n = disk::normal(msh, e_cell, face);
        for (auto& qp : qps) {
            const auto v_f_phi = vfb.eval_functions(qp.point());
            const auto s_f_phi = sfb.eval_functions(qp.point());
            assert(v_f_phi.rows() == vfbs);
            assert(s_f_phi.rows() == sfbs);
            auto q_n = disk::priv::inner_product(qp.weight(), n);
            const auto n_dot_v_f_phi = disk::priv::inner_product(v_f_phi,q_n);
            const auto result = (-1.0*disk::priv::outer_product(n_dot_v_f_phi, s_f_phi)).eval();
            interface_operator += result;
        }
        return interface_operator;
    }
    
    // Fine region measured by the smallest width of the cells instead of their diameter (time.fine_measure = "width"):
    // the stable step of a cell scales with its width (thin cut cells, slivers), not with its diameter
    void set_fine_measure_width(bool w) { m_fine_by_width = w; }
    bool m_fine_by_width = false;
    // smallest width (altitude) of a convex polygonal cell: min over its edges of the largest distance of its vertices
    // to the line of the edge
    static T cell_width(const Mesh& msh, const typename Mesh::cell_type& cl) {
        auto pts = points(msh, cl);
        const size_t n = pts.size();
        T best = std::numeric_limits<T>::max();
        for (size_t i = 0; i < n; i++) {
            const auto a = pts[i], b = pts[(i + 1) % n];
            const T tx = b.x() - a.x(), ty = b.y() - a.y(), L = std::sqrt(tx * tx + ty * ty);
            if (L == 0) continue;
            T far = 0;
            for (auto& q : pts) far = std::max(far, std::abs(-(q.x() - a.x()) * ty + (q.y() - a.y()) * tx) / L);
            best = std::min(best, far);
        }
        return best;
    }
    T fine_measure(const Mesh& msh, const typename Mesh::cell_type& cl) const {
        return m_fine_by_width ? cell_width(msh, cl) : diameter(msh, cl);
    }

    // LTS fine region: cells whose diameter (rounded to 1e-3) is below h_c, extended by nb_layer layers of cells sharing
    // a vertex with a fine cell, and faces with a fine neighbour cell. Flags in the elastic / acoustic cell numbering
    // (m_e_cell_index, m_a_cell_index) and in the compressed face numbering of the unknowns (m_e/a_compress_indexes).
    void mark_fine_region(const Mesh& msh, T h_c, size_t nb_layer,
                          std::vector<bool>& e_cell_is_fine, std::vector<bool>& a_cell_is_fine,
                          std::vector<bool>& e_face_is_fine, std::vector<bool>& a_face_is_fine) {

        const T hc = std::round(h_c * T(1000)) / T(1000);
        auto storage = msh.backend_storage();

        // --------------------------------------------------
        // STEP 1 : mark cells as fine based on mesh size
        // --------------------------------------------------
        e_cell_is_fine.assign(m_e_cell_index.size(), false);
        a_cell_is_fine.assign(m_a_cell_index.size(), false);
        
        for (auto& chunk : m_e_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_e_cell_index[cell_id];
            auto& cell = storage->surfaces[cell_id];
            const T h = std::round(fine_measure(msh, cell) * T(1000)) / T(1000);
            if (h < hc) {
                e_cell_is_fine[cell_ind] = true;
            }
        }
        
        for (auto& chunk : m_a_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_a_cell_index[cell_id];
            auto& cell = storage->surfaces[cell_id];
            const T h = std::round(fine_measure(msh, cell) * T(1000)) / T(1000);
            if (h < hc) {
                a_cell_is_fine[cell_ind] = true;
            }
        }
        
        // --------------------------------------------------
        // STEP 2 : extend fine region by nb_layer layers
        // A coarse cell is marked fine if it shares at least one
        // point (vertex) with a fine cell — this captures diagonal
        // neighbours (corners) that share only a vertex, not a face.
        // The rings are geometric: a fine cell of either medium
        // annexes the cells of both media (the fine / coarse
        // partition does not depend on the elasto-acoustic
        // interface; same rings as the overlap of the multilevel
        // LTS, mesh/domain.hpp, mlts_cell_levels).
        // --------------------------------------------------
        for (size_t layer = 0; layer < nb_layer; ++layer)
        {
            // Collect all point ids belonging to currently fine elastic cells (flags indexed by point id)
            std::vector<char> e_fine_points(msh.points_size(), 0);
            for (auto& chunk : m_e_material)
            {
                const size_t cell_id  = chunk.first;
                const size_t cell_ind = m_e_cell_index[cell_id];
                if (!e_cell_is_fine[cell_ind]) continue;
                auto& cell = storage->surfaces[cell_id];
                for (auto pt_id : cell.point_ids())
                    e_fine_points[size_t(pt_id)] = 1;
            }

            // Collect all point ids belonging to currently fine acoustic cells
            std::vector<char> a_fine_points(msh.points_size(), 0);
            for (auto& chunk : m_a_material)
            {
                const size_t cell_id  = chunk.first;
                const size_t cell_ind = m_a_cell_index[cell_id];
                if (!a_cell_is_fine[cell_ind]) continue;
                auto& cell = storage->surfaces[cell_id];
                for (auto pt_id : cell.point_ids())
                    a_fine_points[size_t(pt_id)] = 1;
            }

            // Mark any elastic cell that shares a point with a fine cell (of either medium)
            for (auto& chunk : m_e_material)
            {
                const size_t cell_id  = chunk.first;
                const size_t cell_ind = m_e_cell_index[cell_id];
                auto& cell = storage->surfaces[cell_id];
                for (auto pt_id : cell.point_ids())
                {
                    if (e_fine_points[size_t(pt_id)] || a_fine_points[size_t(pt_id)])
                    {
                        e_cell_is_fine[cell_ind] = true;
                        break;
                    }
                }
            }

            // Mark any acoustic cell that shares a point with a fine cell (of either medium)
            for (auto& chunk : m_a_material)
            {
                const size_t cell_id  = chunk.first;
                const size_t cell_ind = m_a_cell_index[cell_id];
                auto& cell = storage->surfaces[cell_id];
                for (auto pt_id : cell.point_ids())
                {
                    if (a_fine_points[size_t(pt_id)] || e_fine_points[size_t(pt_id)])
                    {
                        a_cell_is_fine[cell_ind] = true;
                        break;
                    }
                }
            }
        }

        // --------------------------------------------------
        // STEP 3 : mark faces — fine if at least one neighbor cell is fine
        // --------------------------------------------------
        e_face_is_fine.assign(n_e_edges, false);
        a_face_is_fine.assign(n_a_edges, false);
        
        for (auto& chunk : m_e_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_e_cell_index[cell_id];
            auto& cell = storage->surfaces[cell_id];
            auto fcs = faces(msh, cell);
            for (auto& fc : fcs) {
                const auto fc_id = msh.lookup(fc);
                if (m_e_bnd.is_dirichlet_face(fc_id)) continue;
                const size_t cf = m_e_compress_indexes.at(fc_id);
                if (e_cell_is_fine[cell_ind]) e_face_is_fine[cf] = true;
            }
        }
        
        for (auto& chunk : m_a_material) {
            const size_t cell_id  = chunk.first;
            const size_t cell_ind = m_a_cell_index[cell_id];
            auto& cell = storage->surfaces[cell_id];
            auto fcs = faces(msh, cell);
            for (auto& fc : fcs) {
                const auto fc_id = msh.lookup(fc);
                if (m_a_bnd.is_dirichlet_face(fc_id)) continue;
                const size_t cf = m_a_compress_indexes.at(fc_id);
                if (a_cell_is_fine[cell_ind]) a_face_is_fine[cf] = true;
            }
        }
    }

    // Classification of mesh cells by boundary conditions: elastic (resp. acoustic) cells with at least one
    // elastic (resp. acoustic) Dirichlet face, in mesh order. Used by apply_bc.
    void classify_cells(const Mesh& msh) {

        auto has_dirichlet_face = [&](const typename Mesh::cell_type& cell, auto& bnd) {
            for (auto& fc : faces(msh, cell)) {
                if (bnd.is_dirichlet_face(msh.lookup(fc)))
                    return true;
            }
            return false;
        };

        m_e_elements_with_bc_eges.clear();
        m_a_elements_with_bc_eges.clear();
        for (auto& cell : msh) {
            auto cell_ind = msh.lookup(cell);
            if (m_e_material.count(cell_ind) && has_dirichlet_face(cell, m_e_bnd))
                m_e_elements_with_bc_eges.push_back(cell_ind);
        }
        for (auto& cell : msh) {
            auto cell_ind = msh.lookup(cell);
            if (m_a_material.count(cell_ind) && has_dirichlet_face(cell, m_a_bnd))
                m_a_elements_with_bc_eges.push_back(cell_ind);
        }
    }
    
    // Construction of cell-to-index maps; 
    void build_cells_maps() {
        
        // elastic data
        size_t e_cell_ind = 0;
        for (auto chunk : m_e_material) {
            m_e_cell_index.insert(std::make_pair(chunk.first,e_cell_ind));
            e_cell_ind++;
        }
        m_n_elastic_cells = e_cell_ind;

        // acoustic data
        size_t a_cell_ind = 0;
        for (auto chunk : m_a_material) {
            m_a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
            a_cell_ind++;
        }
        m_n_acoustic_cells = a_cell_ind;

    }
    
    // Elastic cell dofs update with projected values.
    void
    scatter_e_cell_dof_data(size_t e_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell, Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const {
        size_t n_cbs = get_e_cell_basis_data();
        auto cell_ofs = e_cell_ind * n_cbs;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    // Acoustic cell dofs update with projected values.
    void
    scatter_a_cell_dof_data(size_t a_cell_ind, const Mesh& msh, const typename Mesh::cell_type& cell, Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const {
        size_t n_cbs = get_a_cell_basis_data();
        auto cell_ofs = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
        x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
    }
    
    // Function projection onto face dofs 
    void project_over_faces(const Mesh& msh, Matrix<T, Dynamic, 1> & x_glob, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> vec_fun, std::function<T(const typename Mesh::point_type& )> scal_fun) {

        auto storage = msh.backend_storage();

        // elastic block
        for (auto e_chunk : m_e_material) {
            auto& cell = storage->surfaces[e_chunk.first];
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++) {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_e_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q) {
                    continue;
                }
                Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, face, m_hho_di.face_degree(), vec_fun);
                scatter_e_face_dof_data(msh, face, x_glob, x_proj_dof);
            }
        }
        
        // acoustic block
        for (auto a_chunk : m_a_material) {
            auto& cell = storage->surfaces[a_chunk.first];
            auto fcs = faces(msh, cell);
            for (size_t i = 0; i < fcs.size(); i++) {
                auto face = fcs[i];
                auto fc_id = msh.lookup(face);
                bool is_dirichlet_Q = m_a_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q) {
                    continue;
                }
                Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, face, m_hho_di.face_degree(), scal_fun);
                scatter_a_face_dof_data(msh, face, x_glob, x_proj_dof);
            }
        }
    }
    
    // Updating elastic face dofs projection
    void scatter_e_face_dof_data(const Mesh& msh, const typename Mesh::face_type& face, Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const {
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        auto fc_id = msh.lookup(face);
        auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
    }
    
    // Updating acoustic face dofs projection
    void scatter_a_face_dof_data( const Mesh& msh, const typename Mesh::face_type& face, Matrix<T, Dynamic, 1>& x_glob, Matrix<T, Dynamic, 1> x_proj_dof) const {
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        auto fc_id = msh.lookup(face);
        auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
        x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
    }
    
    // Gathering of elastic cell and face DoFs.
    Matrix<T, Dynamic, 1> gather_e_dof_data(size_t cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& x_glob) const {
                        
        auto e_cell_ind = m_e_cell_index.find(cell_ind)->second;
        auto num_faces = howmany_faces(msh, cl);
        size_t n_cbs = get_e_cell_basis_data();
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        
        Matrix<T, Dynamic, 1> x_el(n_cbs + num_faces * n_fbs );
        x_el.block(0, 0, n_cbs, 1) = x_glob.block(e_cell_ind * n_cbs, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++) {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) {
                throw std::invalid_argument("This is a bug: face not found");
            }
            const auto face_id = eid.second;

            if (m_e_bnd.is_dirichlet_face( face_id)) {
                auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
                auto dirichlet_fun = m_e_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
                x_el.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else {
                auto fc_id = msh.lookup(fc);
                auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_e_compress_indexes.at(fc_id)*n_fbs;
                x_el.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
            }
        }
        return x_el;
    }
    
    // Gathering of acoustic cell and face DoFs.
    Matrix<T, Dynamic, 1>
    gather_a_dof_data(size_t cell_ind, const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& x_glob) const {

        auto a_cell_ind = m_a_cell_index.find(cell_ind)->second;
        auto num_faces = howmany_faces(msh, cl);
        size_t n_cbs = get_a_cell_basis_data();
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
        Matrix<T, Dynamic, 1> x_local(n_cbs + num_faces * n_fbs );
        x_local.block(0, 0, n_cbs, 1) = x_glob.block(a_cell_ind * n_cbs + m_n_elastic_cell_dof, 0, n_cbs, 1);
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++) {
            auto fc = fcs[i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) {
                throw std::invalid_argument("This is a bug: face not found");
            }
            const auto face_id = eid.second;
            if (m_a_bnd.is_dirichlet_face( face_id)) {
                auto fb = disk::make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, m_hho_di.face_degree());
                auto velocity = m_a_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, velocity, m_hho_di.face_degree());
                x_local.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
            }
            else {
                auto fc_id = msh.lookup(fc);
                auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof + m_a_compress_indexes.at(fc_id)*n_fbs;
                x_local.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
            }
        }
        return x_local;
    }
    
    void finalize_mass() {
        MASS.setFromTriplets( m_mass_triplets.begin(), m_mass_triplets.end() );
        std::vector<Triplet<T>>().swap(m_mass_triplets);    // clear and release the memory
    }
    
    void finalize_coupling() {
        COUPLING.setFromTriplets( m_c_triplets.begin(), m_c_triplets.end() );
        std::vector<Triplet<T>>().swap(m_c_triplets);    // clear and release the memory
    }
    
    void set_hdg_stabilization() {
        if(m_hho_di.cell_degree() > m_hho_di.face_degree()) {
            m_hho_stabilization_Q = false;                                            
        }
        else {
            m_hho_stabilization_Q = true;
        }
    }
    
    // Interface coupling integrated with the degree-k rule of the original code (k/2 + 1 Gauss points, rank-deficient
    // for k >= 1) instead of the exact degree-2k rule: only to reproduce published results (space.interface_quadrature)
    void set_legacy_interface_quadrature(bool legacy_Q) {
        m_legacy_interface_quadrature = legacy_Q;
    }

    void set_interface_cell_indexes(std::map<size_t,std::pair<size_t,size_t>> & interface_cell_indexes) {
        m_interface_cell_indexes = interface_cell_indexes;
    }
    
    void set_hho_stabilization() {
        m_hho_stabilization_Q = true;
    }

    // Zero entries of the local matrices not stored in LHS and MASS (the matrix-free operator; the other integrators
    // keep the pattern of the local matrices, used by the orderings of the direct solvers)
    void set_skip_zero_entries(bool skip) {
        m_skip_zero_entries = skip;
    }

    // Least-squares stabilization in any order: plain least squares in equal order, Lehrenfeld-Schoberl in mixed order
    void set_least_squares_stabilization() {
        m_hho_stabilization_Q = false;
    }
    
    e_boundary_type & get_e_bc_conditions() {
        return m_e_bnd;
    }
    
    a_boundary_type & get_a_bc_conditions() {
        return m_a_bnd;
    }
    
    std::map<size_t,elastic_material_data<T>> & get_e_material_data() {
        return m_e_material;
    }
    
    std::map<size_t,acoustic_material_data<T>> & get_a_material_data() {
        return m_a_material;
    }
    
    size_t get_a_n_cells_dof() const {
        return m_n_acoustic_cell_dof;
    }
    
    size_t get_e_n_cells_dof() const {
        return m_n_elastic_cell_dof;
    }
    
    size_t get_n_face_dof() const {
        size_t n_face_dof = m_n_elastic_face_dof + m_n_acoustic_face_dof;
        return n_face_dof;
    }
    
    size_t get_e_face_dof() const {
        return m_n_elastic_face_dof;
    }
    
    size_t get_a_face_dof() const {
        return m_n_acoustic_face_dof;
    }
    
    size_t get_e_face_basis_data() {
        size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        return n_fbs;
    }
    
    size_t get_a_face_basis_data() {
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        return n_fbs;
    }
    
    size_t get_acoustic_cells() {
        return m_n_acoustic_cells;
    }
    
    size_t get_elastic_cells() {
        return m_n_elastic_cells;
    }
    
    size_t get_elastic_faces() {
        return n_e_edges;
    }
    
    size_t get_acoustic_faces() {
        return n_a_edges;
    }
    
    std::map<size_t,std::pair<size_t,size_t>>  get_interfaces() {
        return m_interface_cell_indexes;
    }
    
    std::vector<size_t> get_e_compress() {
        return m_e_compress_indexes;
    }

    // mesh cell id -> index of the elastic / acoustic cell (order of the cell unknowns)
    const std::map<size_t,size_t>& get_e_cell_index() const {
        return m_e_cell_index;
    }

    const std::map<size_t,size_t>& get_a_cell_index() const {
        return m_a_cell_index;
    }
    
    std::vector<size_t> get_a_compress() {
        return m_a_compress_indexes;
    }
    
    std::vector<size_t> get_e_expand() {
        return m_e_expand_indexes;
    }
    
    std::vector<size_t> get_a_expand() {
        return m_a_expand_indexes;
    }
    
    void finalize() {
        LHS.setFromTriplets( m_triplets.begin(), m_triplets.end() );
        std::vector<Triplet<T>>().swap(m_triplets);    // clear and release the memory
    }
    
    size_t get_e_cell_basis_data() const {
        return m_n_e_cbs;
    }
    
    size_t get_a_cell_basis_data() const {
        return m_n_a_cbs;
    }

protected:
    bool        m_skip_zero_entries = false;    // see set_skip_zero_entries

};

#endif /* elastoacoustic_assembler_base_hpp */
