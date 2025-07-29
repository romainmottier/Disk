
#include <iostream> 
#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "../common/assembly_index.hpp"
#include "../common/material_data.hpp"

using Tuple = std::tuple<double, disk::location, std::vector<double>>;

template<typename Mesh>
class unfitted_elliptic_interface_assembler {

    typedef disk::BoundaryConditions<Mesh, true> boundary_type;
    using T = typename Mesh::coordinate_type;

    std::vector<size_t>             m_compress_indexes;
    std::vector<size_t>             m_expand_indexes;
    disk::hho_degree_info           m_hho_di;
    const boundary_type&                  m_bnd;
    std::vector< Triplet<T> >       m_triplets;
    std::vector< Triplet<T> >       m_mass_triplets;
    std::vector< material_data<T> > m_material;
    std::vector< size_t >           m_elements_with_bc_eges;

    size_t      m_n_edges;
    size_t      m_n_essential_edges;
    bool        m_hho_stabilization_Q;

    //unfitted
    std::vector<size_t> m_face_table;
    std::vector<size_t> m_cell_table;
    disk::location loc_zone; 
    size_t m_num_cells;


public:

    SparseMatrix<T>       LHS;
    Matrix<T, Dynamic, 1> RHS;
    SparseMatrix<T>       MASS;


    unfitted_elliptic_interface_assembler(const Mesh& msh, const disk::hho_degree_info& hho_di, const boundary_type& bnd) : m_hho_di(hho_di), m_bnd(bnd), m_hho_stabilization_Q(true) {

        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return m_bnd.is_dirichlet_face(fc_id);
        };
        
        m_cell_table.clear();
        size_t m_num_cells = 0;
        for (auto& cl : msh) {
            m_cell_table.push_back(m_num_cells);
            if (disk::locate(msh, cl) == disk::location::ON_INTERFACE)
                m_num_cells += 2; 
            else
                m_num_cells += 1;
        }

        m_n_edges = 0;
        for (auto& fc : faces(msh)) {
            bool is_cut = (disk::locate(msh, fc) == disk::location::ON_INTERFACE);
            if (is_cut) 
                m_n_edges += 2;
            else
                m_n_edges += 1;
        }
            
        m_compress_indexes.resize(m_n_edges, static_cast<size_t>(-1));
        m_expand_indexes.resize(m_n_edges - m_n_essential_edges, static_cast<size_t>(-1));
        m_face_table.resize(msh.faces_size(), static_cast<size_t>(-1));
        
        size_t fc_cpt = 0;    
        size_t compressed_face_cpt = 0; 
        m_face_table.resize(msh.faces_size(), static_cast<size_t>(-1));
        m_compress_indexes.resize(m_n_edges, static_cast<size_t>(-1));
        m_expand_indexes.resize(m_n_edges - m_n_essential_edges, static_cast<size_t>(-1));
        for (size_t i = 0; i < msh.faces_size(); ++i) {
            auto& fc = *std::next(msh.faces_begin(), i);
            bool is_cut = (disk::locate(msh, fc) == disk::location::ON_INTERFACE);
            m_face_table[i] = fc_cpt;
            size_t nb_copies = is_cut ? 2 : 1;
            for (size_t c = 0; c < nb_copies; ++c) {
                if (!is_dirichlet(fc)) {
                    m_compress_indexes[fc_cpt] = compressed_face_cpt;
                    m_expand_indexes[compressed_face_cpt] = fc_cpt;
                    compressed_face_cpt++;
                }
                fc_cpt++;
            }
        }
        
        size_t n_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
        size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
        
        size_t system_size = n_cbs * m_num_cells + n_fbs * (m_n_edges - m_n_essential_edges);
        
        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = Matrix<T, Dynamic, 1>::Zero(system_size);
        MASS = SparseMatrix<T>(system_size, system_size);
        
        classify_cells(msh);
    }
    
    void classify_cells(const Mesh& msh){
        
        m_elements_with_bc_eges.clear();
        size_t cell_ind = 0;
        for (auto& cell : msh) {
            auto face_list = faces(msh, cell);
            for (size_t face_i = 0; face_i < face_list.size(); face_i++) {
                auto fc = face_list[face_i];
                auto fc_id = msh.lookup(fc);
                bool is_dirichlet_Q = m_bnd.is_dirichlet_face(fc_id);
                if (is_dirichlet_Q) {
                    m_elements_with_bc_eges.push_back(cell_ind);
                    break;
                }
            }
            cell_ind++;
        }
    }

    size_t
    face_SOL_offset(const Mesh& msh, const typename Mesh::face_type& fc) {
        auto cbs = scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);     
        auto fbs = scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);                
        auto face_offset = offset(msh, fc);
        return m_num_cells * cbs + m_face_table.at(face_offset) * fbs;
    }

    std::vector<assembly_index>
    init_asm_map_poly_ext(const Mesh& msh, Tuple P) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        
        // DOFS
        auto celdeg = m_hho_di.cell_degree();
        auto facdeg = m_hho_di.face_degree();
        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);                
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();      
        std::vector<assembly_index> asm_map;

        // CELL OFFSET
        size_t cell_offset = m_cell_table.at(offset(msh, cl)); 
        size_t cell_LHS_offset = cell_offset*cbs;
        auto cbs_cut = 2*cbs; // CELL DEGREES OF FREEDOM OF DEPENDENT CELLS 

        ///////////////////////////////////////////// ASSEMBLY OF THE DOFS OF THE CURRENT CELLS 
        // CELL DOFS
        auto cbs_loc = cbs; 
        if (is_cut(msh,cl))
            cbs_loc = 2*cbs;
        for (size_t i = 0; i < cbs_loc; i++)
            asm_map.push_back(assembly_index(cell_LHS_offset+i, true));
        // FACES DOFS
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = face_SOL_offset(msh, fc);
            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back(assembly_index(face_LHS_offset+i, !is_dirichlet(fc)));
        }
        // ASSEMBLY OF THE FACES IN THE POSITIVE SIDE IF THE CELL IS CUT 
        if(is_cut(msh,cl)) {
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto d = (disk::locate(msh, fc) == disk::location::ON_INTERFACE) ? fbs : 0;
                auto face_LHS_offset = face_SOL_offset(msh, fc) + d;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !is_dirichlet(fc)) );
            }
        }
        ///////////////////////////////////////////// ASSEMBLY OF THE DEPENDEND DOFS 
        // DEPENDENT CELLS = CELLS STABILIZED BY THE CURRENT CELL
        // LOOP OVER DEPENDENT CELLS
        auto dp_cells = std::get<2>(P);
        for (auto dp_cl : dp_cells) {
            // CELL DOFS
            auto dp_cell = msh.cells[dp_cl];
            cell_offset = m_cell_table.at(offset(msh, dp_cell)); 
            cell_LHS_offset = cell_offset*cbs;
            for (size_t i = 0; i < cbs_cut; i++)
                asm_map.push_back(assembly_index(cell_LHS_offset+i, true));
            // FACES DOFS
            fcs = faces(msh, dp_cell);
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto face_LHS_offset = face_SOL_offset(msh, fc);
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !is_dirichlet(fc)));
            }
            // ASSEMBLY OF THE FACES IN THE POSITIVE SIDE 
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto d = (disk::locate(msh, fc) == disk::location::ON_INTERFACE) ? fbs : 0;
                auto face_LHS_offset = face_SOL_offset(msh, fc) + d;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !is_dirichlet(fc)) );
            }
        }
         
        return asm_map;
    }
    
};

    // void scatter_data(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {

    //     auto fcs = faces(msh, cl);
    //     size_t cbs = disk::scalar_basis_size(m_hho_di.cell_degree(),Mesh::dimension);
    //     size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension-1);
    //     std::vector<assembly_index> asm_map;
        
    //     asm_map.reserve(n_cbs + n_fbs*fcs.size());

    //     auto cell_offset        = offset(msh, cl);
    //     auto cell_LHS_offset    = cell_offset * n_cbs;

    //     for (size_t i = 0; i < n_cbs; i++)
    //         asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
            
    //     for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    //     {
    //         auto fc = fcs[face_i];
    //         auto face_offset = offset(msh, fc);
    //         auto face_LHS_offset = n_cbs * msh.cells_size() + m_compress_indexes.at(face_offset)*n_fbs;

    //         auto fc_id = msh.lookup(fc);
    //         bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

    //         for (size_t i = 0; i < n_fbs; i++)
    //             asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
    //     }
            
    //     assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

    //     for (size_t i = 0; i < lhs.rows(); i++)
    //     {
    //         if (!asm_map[i].assemble())
    //             continue;

    //         for (size_t j = 0; j < lhs.cols(); j++)
    //         {
    //             if ( asm_map[j].assemble() )
    //                 m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
    //         }
    //     }

    //     for (size_t i = 0; i < rhs.rows(); i++)
    //     {
    //         if (!asm_map[i].assemble())
    //             continue;
    //         RHS(int(asm_map[i])) += rhs(i);
    //     }

    // }
            
    // void scatter_bc_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    //          const Matrix<T, Dynamic, Dynamic>& lhs)
    // {
    
    //     auto fcs = faces(msh, cl);
    //     size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
    //     size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
    //     std::vector<assembly_index> asm_map;
    //     asm_map.reserve(n_cbs + n_fbs*fcs.size());

    //     auto cell_offset        = offset(msh, cl);
    //     auto cell_LHS_offset    = cell_offset * n_cbs;

    //     for (size_t i = 0; i < n_cbs; i++)
    //         asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
    //     Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
    //     for (size_t face_i = 0; face_i < fcs.size(); face_i++)
    //     {
    //         auto fc = fcs[face_i];
    //         auto face_offset = offset(msh, fc);
    //         auto face_LHS_offset = n_cbs * msh.cells_size() + m_compress_indexes.at(face_offset)*n_fbs;

    //         auto fc_id = msh.lookup(fc);
    //         bool dirichlet = m_bnd.is_dirichlet_face(fc_id);

    //         for (size_t i = 0; i < n_fbs; i++)
    //             asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            
    //         if (dirichlet)
    //          {
    //              auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
    //              auto dirichlet_fun  = m_bnd.dirichlet_boundary_func(fc_id);

    //              Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
    //              Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
    //              dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
    //          }
            
    //     }

    //     assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

    //     for (size_t i = 0; i < lhs.rows(); i++)
    //     {
    //         if (!asm_map[i].assemble())
    //             continue;

    //         for (size_t j = 0; j < lhs.cols(); j++)
    //         {
    //             if ( !asm_map[j].assemble() )
    //                 RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
    //         }
    //     }

    // }
            
    // void scatter_rhs_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    //          const Matrix<T, Dynamic, 1>& rhs)
    // {
    
    //     size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
    //     std::vector<assembly_index> asm_map;
    //     asm_map.reserve(n_cbs);

    //     auto cell_offset        = offset(msh, cl);
    //     auto cell_LHS_offset    = cell_offset * n_cbs;

    //     for (size_t i = 0; i < n_cbs; i++)
    //         asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

    //     assert( asm_map.size() == rhs.rows() );

    //     for (size_t i = 0; i < rhs.rows(); i++)
    //     {
    //         if (!asm_map[i].assemble())
    //             continue;
    //         RHS(int(asm_map[i])) += rhs(i);
    //     }

    // }
            
    // void scatter_mass_data(const Mesh& msh, const typename Mesh::cell_type& cl,
    //          const Matrix<T, Dynamic, Dynamic>& mass_matrix)
    // {
    //     size_t n_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
    //     std::vector<assembly_index> asm_map;
    //     asm_map.reserve(n_cbs);

    //     auto cell_offset        = offset(msh, cl);
    //     auto cell_LHS_offset    = cell_offset * n_cbs;

    //     for (size_t i = 0; i < n_cbs; i++)
    //         asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

    //     assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );

    //     for (size_t i = 0; i < mass_matrix.rows(); i++)
    //     {
    //         if (!asm_map[i].assemble())
    //             continue;

    //         for (size_t j = 0; j < mass_matrix.cols(); j++)
    //         {
    //             if ( asm_map[j].assemble() )
    //                 m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
    //         }
    //     }

    // }
                      
    // size_t
    // face_SOL_offset(const Mesh& msh, const typename Mesh::face_type& fc) {
    //     auto facdeg = m_hho_di.face_degree();
    //     auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);                
    //     auto cbs = loc_cbs; 
    //     auto face_offset = offset(msh, fc);
    //     return num_cells * cbs + face_table.at(face_offset) * fbs;
    // }

    






//     size_t
//     n_dof(const Mesh& msh, const typename Mesh::cell_type& cl) {

//         bool double_unknowns = (location(msh, cl) == location::ON_INTERFACE && loc_zone == location::ON_INTERFACE );
//         auto facdeg = m_hho_di.face_degree();
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);                
//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
//         auto f_dofs = num_faces * fbs;
//         auto cbs = loc_cbs;
//         auto loc_size = cbs + f_dofs;
//         if( double_unknowns )
//             loc_size = 2 * loc_size;
//         return loc_size;
//     }

//     Matrix<T, Dynamic, 1>
//     get_dirichlet_data_poly_ext(const Mesh& msh, Tuple P) {

//         // CELL INFOS
//         auto cell_index = std::get<0>(P);
//         auto cl = msh.cells[cell_index];
//         auto dp_cells = std::get<2>(P);

//         // DOFS
//         auto celdeg = m_hho_di.cell_degree();
//         auto facdeg = m_hho_di.face_degree();
//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);                
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);  
//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
//         auto current_dofs = cbs + num_faces*fbs;
//         if (is_cut(msh,cl)) 
//             current_dofs = 2*current_dofs;
//         auto extended_dofs = 2*(cbs + num_faces*fbs);
//         auto nb_dp_cells = dp_cells.size();
//         auto local_dofs = current_dofs + nb_dp_cells*extended_dofs; 

//         Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(local_dofs);

//         // LOOP OVER FACES
//         for (size_t face_i = 0; face_i < num_faces; face_i++) {
//             auto fc = fcs[face_i];
//             auto face_LHS_offset = face_SOL_offset(msh, fc);
//             bool in_dom = true;
//             if (loc_zone != location::ON_INTERFACE) {
//                 location loc_fc = location(msh, fc);
//                 bool in_dom = (loc_fc == location::ON_INTERFACE || loc_fc == loc_zone);
//             }
//             bool dirichlet = fc.is_boundary && fc.bndtype == DirichletType::DIRICHLET && in_dom;
//             if (dirichlet && loc_zone == location::ON_INTERFACE ) {
//                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
//                 Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, dir_func);
//                 dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
//             }
//             if (dirichlet && loc_zone != location::ON_INTERFACE ) {
//                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg, loc_zone);
//                 Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, loc_zone, dir_func);
//                 dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
//             }
//         }
//         return dirichlet_data;
//     }
    
//     void
//     assemble_bis_poly_ext(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {

//         // CELL INFOS
//         auto cell_index = std::get<0>(P);
//         auto cl = msh.cells[cell_index];

//         if( !(location(msh, cl) == loc_zone
//               || location(msh, cl) == location::ON_INTERFACE
//               || loc_zone == location::ON_INTERFACE ) )
//             return;
        
//         auto asm_map = init_asm_map_ext(msh, P);
//         auto dirichlet_data = get_dirichlet_data_ext(msh, P);
//         assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());
//         assert(dirichlet_data.size() == lhs.cols());

//         // ASSEMBLY OF STIFFNESS MATRIX
//         for (size_t i = 0; i < lhs.rows(); i++) {
//             if (!asm_map[i].assemble())
//                 continue;
//             for (size_t j = 0; j < lhs.cols(); j++) {
//                 if (asm_map[j].assemble()) {
//                     triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
//                 }
//                 else {
// 	                int itmp=asm_map[i];
//                     RHS(itmp) -= lhs(i,j)*dirichlet_data(j);
// 		        }
//             }
//         }
        
//         // ASSEMBLY OF THE RHS
//         for (size_t i = 0; i < rhs.rows(); i++) {
//             if (!asm_map[i].assemble())
//                 continue;
//             RHS[asm_map[i]] += rhs(i);
//         }

//     }


//     Matrix<T, Dynamic, 1>
//     get_solF(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution) {
        
//         bool double_unknowns = ( location(msh, cl) == location::ON_INTERFACE && loc_zone == location::ON_INTERFACE );
//         auto facdeg = m_hho_di.face_degree();
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);  
//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
//         size_t f_dofs = num_faces*fbs;
//         if( double_unknowns )
//             f_dofs = 2 * f_dofs;
            
//         Matrix<T, Dynamic, 1> solF = Matrix<T, Dynamic, 1>::Zero( f_dofs );
//         for (size_t face_i = 0; face_i < num_faces; face_i++) {
//             auto fc = fcs[face_i];
//             if (loc_zone != location::ON_INTERFACE) {
//                 auto loc_fc = location(msh, fc);
//                 if (!(loc_fc == location::ON_INTERFACE || loc_fc == loc_zone) )
//                     continue;
//             }

//             auto face_LHS_offset = face_SOL_offset(msh, fc);
//             if (location(msh, fc) == location::ON_INTERFACE && loc_zone == location::ON_INTERFACE) {
//                 solF.block(face_i*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
//                 solF.block( (num_faces+face_i)*fbs, 0, fbs, 1)
//                     = solution.block(face_LHS_offset + fbs, 0, fbs, 1);
//                 continue;
//             }

//             bool dirichlet = fc.is_boundary && fc.bndtype == DirichletType::DIRICHLET;
//             if (dirichlet) {
//                 Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
//                 Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, facdeg, dir_func);
//                 solF.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
//                 continue;
//             }

//             if (location(msh, cl) == location::ON_INTERFACE && location(msh, fc) == location::IN_POSITIVE_SIDE && loc_zone == location::ON_INTERFACE) {
//                 solF.block((num_faces+face_i)*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
//                 continue;
//             }

//             solF.block(face_i*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
        
//         }
//         return solF;
//     }
         
//     void 
//     finalize(void) {
//         LHS.setFromTriplets( triplets.begin(), triplets.end() );
//         triplets.clear();
//         MASS.setFromTriplets( triplets_mass.begin(), triplets_mass.end() );
//         triplets_mass.clear();
//     }




// template<typename Mesh, typename Function>
// class virt_interface_assembler : public virt_scalar_assembler<Mesh, Function> {

//     using T = typename Mesh::coordinate_type;

// public:

//     virt_interface_assembler(const Mesh& msh, const Function& dirichlet_bf, hho_degree_info hdi) : virt_scalar_assembler<Mesh, Function>(msh, dirichlet_bf, hdi) {

//         this->loc_zone = location::ON_INTERFACE;
//         auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
//             return fc.is_boundary && fc.bndtype == DirichletType::DIRICHLET;
//         };
//         size_t loc_num_cells = 0; 
//         for (auto& cl : msh.cells) {
//             this->cell_table.push_back(loc_num_cells);
//             if (location(msh, cl) == location::ON_INTERFACE)
//                 loc_num_cells += 2;
//             else
//                 loc_num_cells += 1;
//         }
//         this->num_cells = loc_num_cells;
//         assert(this->cell_table.size() == msh.cells.size());

//         size_t num_all_faces = 0; 
//         for (auto& fc : msh.faces) {
//             if (location(msh, fc) == location::ON_INTERFACE)
//                 num_all_faces += 2;
//             else
//                 num_all_faces += 1;
//         }
            
//         size_t num_dirichlet_faces = 0; 
//         for (auto& fc : msh.faces) {
//             if(fc.is_boundary && fc.bndtype == DirichletType::DIRICHLET){
//                 if (location(msh, fc) == location::ON_INTERFACE)
//                         num_dirichlet_faces += 2;
//                     else
//                         num_dirichlet_faces += 1;
//             }
//         }
        
//         this->num_other_faces = num_all_faces - num_dirichlet_faces;
//         this->face_table.resize( msh.faces.size() );

//         size_t compressed_offset = 0;
//         for (size_t i = 0; i < msh.faces.size(); i++) {
//             auto fc = msh.faces.at(i);
//             if ( !is_dirichlet(fc) ) {
//                 this->face_table.at(i) = compressed_offset;
//                 if ( location(msh, fc) == location::ON_INTERFACE )
//                     compressed_offset += 2;
//                 else
//                     compressed_offset += 1;
//             }
//         }
//     }

// };


//     interface_assembler(const Mesh& msh, const Function& dirichlet_bf, hho_degree_info hdi) : virt_interface_assembler<Mesh, Function>(msh, dirichlet_bf, hdi) {

//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();
//         auto graddeg = this->di.grad_degree();
//         auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
//         auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);                
//         this->loc_cbs = cbs;
//         auto system_size = cbs * this->num_cells + fbs * this->num_other_faces;
//         this->LHS = SparseMatrix<T>(system_size, system_size);
//         this->RHS = Matrix<T, Dynamic, 1>::Zero(system_size);        
//         this->MASS = SparseMatrix<T>(system_size, system_size);

//     }

//     void
//     assemble_poly_ext(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs, 
//     const Matrix<T, Dynamic, 1>& rhs) {
//         this->assemble_bis_poly_ext(msh, P, lhs, rhs);
//     }

//     Matrix<T, Dynamic, 1>
//     take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution, location where) {

//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();
//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);
//         auto cell_offset = offset(msh, cl);
//         size_t cell_SOL_offset;
//         if (location(msh, cl) == location::ON_INTERFACE) {
//             if (where == location::IN_NEGATIVE_SIDE)
//                 cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
//             else if (where == location::IN_POSITIVE_SIDE)
//                 cell_SOL_offset = this->cell_table.at(cell_offset) * cbs + cbs;
//             else
//                 throw std::invalid_argument("Invalid location");
//         }
//         else 
//             cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;

//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();

//         Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + num_faces*fbs);
//         ret.block(0, 0, cbs, 1) = solution.block(cell_SOL_offset, 0, cbs, 1);
//         auto solF = this->get_solF(msh, cl, solution);
//         if (where == location::IN_NEGATIVE_SIDE)
//             ret.tail(num_faces * fbs) = solF.head(num_faces * fbs);
//         else
//             ret.tail(num_faces * fbs) = solF.tail(num_faces * fbs);

//         return ret;

//     }
    
//     Matrix<T, Dynamic, 1>
//     gather_proj(const Mesh& msh, Tuple P, hho_degree_info hdi, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
//         // CELL INFOS 
//         auto cell_index = std::get<0>(P);
//         auto loc = std::get<1>(P);
//         auto cl = msh.cells[cell_index];

//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();
//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);
//         auto fbs_cut = 2*fbs;

//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
//         auto uncut_dofs = cbs + num_faces*fbs;
//         auto current_dofs = uncut_dofs;
//         if (is_cut(msh,cl)) 
//             current_dofs = 2*current_dofs;
//         auto extended_dofs = 2*(cbs + num_faces*fbs);
//         auto dp_cells = std::get<2>(P);
//         auto local_dofs = current_dofs + dp_cells.size()*extended_dofs; 

//         Matrix<T, Dynamic, 1> dofs = Matrix<T, Dynamic, 1>::Zero(local_dofs);

//         if (!is_cut(msh, cl)) 
//             dofs.block(0, 0, cbs+num_faces*fbs, 1) = project_function(msh, cl, hdi, scal_fun);
//         else {
//             dofs.block(0, 0, cbs, 1) = project_function(msh, cl, hdi, location::IN_NEGATIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
//             dofs.block(cbs, 0, cbs, 1) = project_function(msh, cl, hdi, location::IN_POSITIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
//             dofs.block(2*cbs, 0, num_faces*fbs, 1) = project_function(msh, cl, hdi, location::IN_NEGATIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
//             dofs.block(2*cbs+num_faces*fbs, 0, num_faces*fbs, 1) = project_function(msh, cl, hdi, location::IN_POSITIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
//         }

//         // LOOP OVER DEPENDENT CELLS 
//         auto offset_dofs = current_dofs;  
//         for (auto &dp_cl : dp_cells) {
//             auto dp_cell = msh.cells[dp_cl];
//             dofs.block(offset_dofs, 0, cbs, 1) = project_function(msh, dp_cell, hdi, location::IN_NEGATIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
//             dofs.block(offset_dofs+cbs, 0, cbs, 1) = project_function(msh, dp_cell, hdi, location::IN_POSITIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
//             dofs.block(offset_dofs+2*cbs, 0, num_faces*fbs, 1) = project_function(msh, dp_cell, hdi, location::IN_NEGATIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
//             dofs.block(offset_dofs+2*cbs+num_faces*fbs, 0, num_faces*fbs, 1) = project_function(msh, dp_cell, hdi, location::IN_POSITIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
//             offset_dofs += extended_dofs;
//         }

//         return dofs;

//     }

//     Matrix<T, Dynamic, 1>
//     gather_cell_dof(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution, location where) {
        
//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();

//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

//         auto cell_offset = offset(msh, cl);
//         size_t cell_SOL_offset;
//         if (location(msh, cl) == location::ON_INTERFACE) {
//             if (where == location::IN_NEGATIVE_SIDE)
//                 cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
//             else if (where == location::IN_POSITIVE_SIDE)
//                 cell_SOL_offset = this->cell_table.at(cell_offset) * cbs + cbs;
//             else
//                 throw std::invalid_argument("Invalid location");
//         }
//         else {
//             cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
//         }
//         return solution.block(cell_SOL_offset, 0, cbs, 1);
//     }

//     void project_over_cells_and_faces(const Mesh& msh, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {

//         for (auto& cl : msh.cells) {
//             if( location(msh, cl) != location::ON_INTERFACE ) 
//                 project_over_uncutcells(msh, cl, hho_di, x_glob, scal_fun); 
//             else
//                 project_over_cutcells(msh, cl, hho_di, x_glob, scal_fun);
//         }
//     }

//     void project_over_uncutcells(const Mesh& msh, const typename Mesh::cell_type& cl, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
//         Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cl, hho_di, scal_fun);

//         // HHO DISCRETIZATION INFOS
//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();
//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);
//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
  
//         // CELL DOFS 
//         auto cell_offset = offset(msh, cl);
//         size_t cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
//         x_glob.block(cell_SOL_offset, 0, cbs, 1) = x_proj_dof.block(0, 0, cbs, 1);

//         // FACE DOFS 
//         for (size_t face_i = 0; face_i < num_faces; face_i++) {
//             auto fc = fcs[face_i];
//             auto face_LHS_offset = this->face_SOL_offset(msh, fc);
//             x_glob.block(face_LHS_offset, 0, fbs, 1) = x_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
//         }

//     }
            
//     void project_over_cutcells(const Mesh& msh, const typename Mesh::cell_type& cl, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
//         Matrix<T, Dynamic, 1> x_neg_proj_dof = project_function(msh, cl, hho_di, location::IN_NEGATIVE_SIDE, scal_fun);
//         Matrix<T, Dynamic, 1> x_pos_proj_dof = project_function(msh, cl, hho_di, location::IN_POSITIVE_SIDE, scal_fun);
            
        
        
//         auto celdeg = this->di.cell_degree();
//         auto facdeg = this->di.face_degree();
//         auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
//         auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1); 
//         auto fcs = faces(msh, cl);
//         auto num_faces = fcs.size();
//         auto cell_offset = offset(msh, cl);
//         size_t cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
//         x_glob.block(cell_SOL_offset, 0, cbs, 1) = x_neg_proj_dof.block(0, 0, cbs, 1);
//         x_glob.block(cell_SOL_offset + cbs, 0, cbs, 1) = x_pos_proj_dof.block(0, 0, cbs, 1);
        
//         // FACE DOFS 
//         for (size_t face_i = 0; face_i < num_faces; face_i++) {
//             auto fc = fcs[face_i];
//             auto face_LHS_offset = this->face_SOL_offset(msh, fc);
//             x_glob.block(face_LHS_offset, 0, fbs, 1) = x_neg_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
//             x_glob.block(face_LHS_offset + fbs, 0, fbs, 1) = x_pos_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
//         }

//     }
    
//     void
//     assemble_mass(const Mesh& msh, const typename Mesh::cell_type& cl, 
//                   const Matrix<T, Dynamic, Dynamic>& mass) {
        
//         this->assemble_bis_mass(msh, cl, mass);
    
//     }

//     void
//     assemble_rhs(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {
//         this->assemble_rhs_bis(msh, cl, rhs);
//     }

