/*
 *       /\        Matteo Cicuttin (C) 2017,2018; Guillaume Delay 2018,2019
 *      /__\       matteo.cicuttin@enpc.fr        guillaume.delay@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is ProtoN, a library for fast Prototyping of
 *  /_\/_\/_\/_\   Numerical methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

using Tuple = std::tuple<double,element_location,std::vector<double>>;

template<typename Mesh, typename Function>
class virt_scalar_assembler {

    using T = typename Mesh::coordinate_type;

protected:
    
    std::vector< Triplet<T>> triplets;
    std::vector< Triplet<T>> triplets_mass;
    std::vector<size_t> face_table;
    std::vector<size_t> cell_table;
    hho_degree_info di;
    Function dir_func;
    element_location loc_zone; 
    size_t num_cells, num_other_faces, loc_cbs, loc_gbs;

    // DEBUG SCHEME
    std::vector< Triplet<T>> triplets_GRAD;
    std::vector< Triplet<T>> triplets_GRAD_GRAD;
    std::vector< Triplet<T>> triplets_sparsity;
    std::vector< Triplet<T>> triplets_zip;

public:

    SparseMatrix<T> LHS;
    SparseMatrix<T> MASS;
    Matrix<T, Dynamic, 1> RHS;

    // DEBUG SCHEME
    Matrix<T, Dynamic, 1> GRAD;
    SparseMatrix<T> GLOBAL_GRAD_GRAD;
    SparseMatrix<T> SPARSITY;
    Matrix<T, Dynamic, 1> CONDITIONING;
    SparseMatrix<T> Kg_ZIP;

    auto get_cell_table() const { return cell_table; }

    virt_scalar_assembler(const Mesh& msh, const Function& dirichlet_bf, hho_degree_info hdi) : dir_func(dirichlet_bf), di(hdi) {
    }

    size_t
    face_SOL_offset(const Mesh& msh, const typename Mesh::face_type& fc) {
        auto facdeg = di.face_degree();
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto cbs = loc_cbs; // cbs = 0 if static condensation
        auto face_offset = offset(msh, fc);
        return num_cells * cbs + face_table.at(face_offset) * fbs;
    }

    std::vector<assembly_index>
    init_asm_map(const Mesh& msh, const typename Mesh::cell_type& cl) {
        bool double_unknowns = ( location(msh, cl) == element_location::ON_INTERFACE && loc_zone == element_location::ON_INTERFACE );
        std::vector<assembly_index> asm_map;
        auto facdeg = di.face_degree();
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto f_dofs = num_faces * fbs;
        auto cbs = loc_cbs;
        auto loc_size = cbs + f_dofs;
        if( double_unknowns )
            loc_size = 2 * loc_size;
        asm_map.reserve( loc_size );
        size_t cell_offset = cell_table.at( offset(msh, cl) );
        size_t cell_LHS_offset = cell_offset * cbs;

        if( double_unknowns )
            cbs = 2 * cbs;

        for (size_t i = 0; i < cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
        
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = face_SOL_offset(msh, fc);
            bool in_dom = true;
            if( loc_zone != element_location::ON_INTERFACE ) {
                element_location loc_fc = location(msh, fc);
                in_dom = (loc_fc == element_location::ON_INTERFACE || loc_fc == loc_zone);
            }
            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET && in_dom;
            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
        }

        if (double_unknowns) {
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto d = (location(msh, fc) == element_location::ON_INTERFACE) ? fbs : 0;
                auto face_LHS_offset = face_SOL_offset(msh, fc) + d;
                bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
                if ( dirichlet )
                    std::cout << "Dirichlet boundary on cut cell detected." << std::endl;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
        }

        return asm_map;
    }

    std::vector<assembly_index>
    init_asm_map_ext(const Mesh& msh, Tuple P) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        
        // DOFS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();      
        std::vector<assembly_index> asm_map;

        // CELL OFFSET
        size_t cell_offset = cell_table.at(offset(msh, cl)); 
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
            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back(assembly_index(face_LHS_offset+i, !dirichlet));
        }
        // ASSEMBLY OF THE FACES IN THE POSITIVE SIDE IF THE CELL IS CUT 
        if(is_cut(msh,cl)) {
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto d = (location(msh, fc) == element_location::ON_INTERFACE) ? fbs : 0;
                auto face_LHS_offset = face_SOL_offset(msh, fc) + d;
                bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
        }
        ///////////////////////////////////////////// ASSEMBLY OF THE DEPENDEND DOFS 
        // DEPENDENT CELLS = CELLS STABILIZED BY THE CURRENT CELL
        // LOOP OVER DEPENDENT CELLS
        auto dp_cells = std::get<2>(P);
        for (auto dp_cl : dp_cells) {
            // CELL DOFS
            auto dp_cell = msh.cells[dp_cl];
            cell_offset = cell_table.at(offset(msh, dp_cell)); 
            cell_LHS_offset = cell_offset*cbs;
            for (size_t i = 0; i < cbs_cut; i++)
                asm_map.push_back(assembly_index(cell_LHS_offset+i, true));
            // FACES DOFS
            fcs = faces(msh, dp_cell);
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto face_LHS_offset = face_SOL_offset(msh, fc);
                bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet));
            }
            // ASSEMBLY OF THE FACES IN THE POSITIVE SIDE 
            for (size_t face_i = 0; face_i < num_faces; face_i++) {
                auto fc = fcs[face_i];
                auto d = (location(msh, fc) == element_location::ON_INTERFACE) ? fbs : 0;
                auto face_LHS_offset = face_SOL_offset(msh, fc) + d;
                bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
                for (size_t i = 0; i < fbs; i++)
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
            }
        }
         
        return asm_map;
    }

    size_t
    n_dof(const Mesh& msh, const typename Mesh::cell_type& cl) {

        bool double_unknowns = ( location(msh, cl) == element_location::ON_INTERFACE && loc_zone == element_location::ON_INTERFACE );
        auto facdeg = di.face_degree();
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto f_dofs = num_faces * fbs;
        auto cbs = loc_cbs;
        auto loc_size = cbs + f_dofs;
        if( double_unknowns )
            loc_size = 2 * loc_size;

        return loc_size;

    }

    Matrix<T, Dynamic, 1>
    get_dirichlet_data(const Mesh& msh, const typename Mesh::cell_type& cl) {

        bool double_unknowns = ( location(msh, cl) == element_location::ON_INTERFACE && loc_zone == element_location::ON_INTERFACE );
        auto facdeg = di.face_degree();
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto f_dofs = num_faces * fbs;
        auto cbs = loc_cbs;
        auto loc_size = cbs + f_dofs;

        if( double_unknowns )
            loc_size *= 2;

        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero( loc_size );
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = face_SOL_offset(msh, fc);
            bool in_dom = true;
            if( loc_zone != element_location::ON_INTERFACE ); 
            {
                element_location loc_fc = location(msh, fc);
                bool in_dom = (loc_fc == element_location::ON_INTERFACE || loc_fc == loc_zone);
            }

            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET && in_dom;
            if( dirichlet && double_unknowns )
                std::cout << "Dirichlet boundary on cut cell detected." << std::endl;
            if (dirichlet && loc_zone == element_location::ON_INTERFACE) {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, dir_func);
                dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
            }
            if (dirichlet && loc_zone != element_location::ON_INTERFACE) {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg, loc_zone);
                Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, loc_zone, dir_func);
                dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
            }
        }

        return dirichlet_data;
    }

    Matrix<T, Dynamic, 1>
    get_dirichlet_data_ext(const Mesh& msh, Tuple P) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        auto dp_cells = std::get<2>(P);

        // DOFS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto current_dofs = cbs + num_faces*fbs;
        if (is_cut(msh,cl)) 
            current_dofs = 2*current_dofs;
        auto extended_dofs = 2*(cbs + num_faces*fbs);
        auto nb_dp_cells = dp_cells.size();
        auto local_dofs = current_dofs + nb_dp_cells*extended_dofs; 

        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(local_dofs);

        // LOOP OVER FACES
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = face_SOL_offset(msh, fc);
            bool in_dom = true;
            if (loc_zone != element_location::ON_INTERFACE) {
                element_location loc_fc = location(msh, fc);
                bool in_dom = (loc_fc == element_location::ON_INTERFACE || loc_fc == loc_zone);
            }
            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET && in_dom;
            if (dirichlet && loc_zone == element_location::ON_INTERFACE ) {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, dir_func);
                dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
            }
            if (dirichlet && loc_zone != element_location::ON_INTERFACE ) {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg, loc_zone);
                Matrix<T, Dynamic, 1> loc_rhs = make_rhs(msh, fc, facdeg, loc_zone, dir_func);
                dirichlet_data.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(loc_rhs);
            }
        }
        return dirichlet_data;
    }

    void
    assemble_bis(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {
        
        if (!(location(msh, cl) == loc_zone || location(msh, cl) == element_location::ON_INTERFACE || loc_zone == element_location::ON_INTERFACE ) )
            return;

        auto asm_map = init_asm_map(msh, cl);
        auto dirichlet_data = get_dirichlet_data(msh, cl);
        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        // LHS
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS[asm_map[i]] -= lhs(i,j)*dirichlet_data(j);
            }
        }

        // RHS
        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            RHS[asm_map[i]] += rhs(i);
        }
    }

    void
    assemble_bis_ext(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];

        if( !(location(msh, cl) == loc_zone
              || location(msh, cl) == element_location::ON_INTERFACE
              || loc_zone == element_location::ON_INTERFACE ) )
            return;
        
        auto asm_map = init_asm_map_ext(msh, P);
        auto dirichlet_data = get_dirichlet_data_ext(msh, P);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());
        assert(dirichlet_data.size() == lhs.cols());

        // ASSEMBLY OF STIFFNESS MATRIX
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) {
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                }
                else {
	                int itmp=asm_map[i];
                    RHS(itmp) -= lhs(i,j)*dirichlet_data(j);
		        }
            }
        }
        
        // ASSEMBLY OF THE RHS
        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            RHS[asm_map[i]] += rhs(i);
        }

    }

    void
    assemble_sparsity(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs) {

        // DOFS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto current_dofs = cbs + num_faces*fbs;
        if (is_cut(msh,cl)) 
            current_dofs = 2*current_dofs;
        
        auto asm_map = init_asm_map(msh, cl);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());

        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) {
                    triplets_sparsity.push_back(Triplet<T>(asm_map[i],asm_map[j],1));
                }
            }
        }
    }

    void
    assemble_sparsity(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        auto dp_cells = std::get<2>(P);

        // DOFS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto current_dofs = cbs + num_faces*fbs;
        if (is_cut(msh,cl)) 
            current_dofs = 2*current_dofs;
        auto extended_dofs = 2*(cbs + num_faces*fbs);
        auto nb_dp_cells = dp_cells.size();
        auto local_dofs = current_dofs + nb_dp_cells*extended_dofs; 
        
        auto asm_map = init_asm_map_ext(msh, P);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());

        // ASSEMBLY OF THE LOCAL CONTRIBUTIONS
        for (size_t i = 0; i < current_dofs; i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < current_dofs; j++) {
                if (asm_map[j].assemble()) 
                    triplets_sparsity.push_back(Triplet<T>(asm_map[i],asm_map[j],1));
            }
        }
        // ASSEMBLY OF THE EXTENDED CONTRIBUTIONS
        for (size_t i = current_dofs; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < current_dofs; j++) {
                if (asm_map[j].assemble()) {
                    triplets_sparsity.push_back( Triplet<T>(asm_map[i], asm_map[j], 200));
                    triplets_sparsity.push_back( Triplet<T>(asm_map[j], asm_map[i], 200));
                }
            }
        }
        for (size_t i = 0; i < current_dofs; i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = current_dofs; j < lhs.cols(); j++) {                
                if (asm_map[j].assemble()) 
                    triplets_sparsity.push_back( Triplet<T>(asm_map[i], asm_map[j], 200));
            }
        }
        for (size_t i = current_dofs; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = current_dofs; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) 
                    triplets_sparsity.push_back( Triplet<T>(asm_map[i], asm_map[j], 200));
            }
        }
    }

    template <typename T>
    SparseMatrix<T>
    condensed_Kg(const Mesh& msh, const SparseMatrix<T>& Kg) {
        
        Matrix<T, Dynamic, Dynamic> denseKg = Kg.toDense();  

        // DOFS INFOS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        
        // BLOCK INFOS
        auto CC_BLOCK_SIZE = this->num_cells*cbs;
        auto FF_BLOCK_SIZE = this->num_other_faces*fbs;
        std::cout << yellow << bold << "         " << "Sparsity profile: " << this->num_cells << "cells" << reset << std::endl;

        // COMPUTATION OF THE FROBENIUS NORM OF CELL-CELL CONTRIBS
        for (size_t cl_i = 0; cl_i < this->num_cells; cl_i++) {
            for (size_t cl_j = 0; cl_j < this->num_cells; cl_j++) {
                auto frobenius_norm = 0.0;
                for (size_t i = 0; i < cbs; i++) {
                    for (size_t j = 0; j < cbs; j++) {
                        frobenius_norm += denseKg(cl_i*cbs+i, cl_j*cbs+j)*denseKg(cl_i*cbs+i, cl_j*cbs+j);
                    }
                }
                // frobenius_norm = std::sqrt(frobenius_norm);
                if (frobenius_norm >= 1e-8)
                    triplets_zip.push_back(Triplet<T>(cl_i, cl_j, frobenius_norm));
            }  
        }
        
        // COMPUTATION OF THE FROBENIUS NORM OF CELL-FACE AND FACE-CELL CONTRIBS
        for (size_t cl_i = 0; cl_i < this->num_cells; cl_i++) {
            for (size_t cl_j = 0; cl_j < this->num_other_faces; cl_j++) {
                auto frobenius_norm = 0.0;
                for (size_t i = 0; i < cbs; i++) {
                    for (size_t j = 0; j < fbs; j++) 
                        frobenius_norm += denseKg(cl_i*cbs+i, CC_BLOCK_SIZE+cl_j*fbs+j)*denseKg(cl_i*cbs+i, CC_BLOCK_SIZE+cl_j*fbs+j);
                }
                // frobenius_norm = std::sqrt(frobenius_norm);
                if (frobenius_norm >= 1e-8) {
                    triplets_zip.push_back(Triplet<T>(cl_i, this->num_cells+cl_j, frobenius_norm));
                    triplets_zip.push_back(Triplet<T>(cl_j + this->num_cells, cl_i, frobenius_norm));
                }
            }  
        }

        
        // COMPUTATION OF THE FROBENIUS NORM OF FACE-FACE CONTRIBS
        for (size_t cl_i = 0; cl_i < this->num_other_faces; cl_i++) {
            for (size_t cl_j = 0; cl_j < this->num_other_faces; cl_j++) {
                auto frobenius_norm = 0.0;
                for (size_t i = 0; i < fbs; i++) {
                    for (size_t j = 0; j < fbs; j++) 
                        frobenius_norm += denseKg(CC_BLOCK_SIZE+cl_i*fbs+i, CC_BLOCK_SIZE+cl_j*fbs+j)*denseKg(CC_BLOCK_SIZE+cl_i*fbs+i, CC_BLOCK_SIZE+cl_j*fbs+j);
                }
                // frobenius_norm = std::sqrt(frobenius_norm);
                if (frobenius_norm >= 1e-8)
                    triplets_zip.push_back(Triplet<T>(this->num_cells+cl_i, this->num_cells+cl_j, frobenius_norm));
            }  
        }

        Kg_ZIP.setFromTriplets(triplets_zip.begin(), triplets_zip.end() );
        triplets_zip.clear();
        
        return Kg_ZIP;
                
    }
    
    void
    assemble_rhs_bis(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {

        if( !(location(msh, cl) == loc_zone || location(msh, cl) == element_location::ON_INTERFACE || loc_zone == element_location::ON_INTERFACE ) )
            return;

        auto asm_map = init_asm_map(msh, cl);
        for (size_t i = 0; i < rhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            RHS[asm_map[i]] += rhs(i);
        }
        
    }

    void
    assemble_bis_mass(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& mass) {

        if( !(location(msh, cl) == loc_zone || location(msh, cl) == element_location::ON_INTERFACE || loc_zone == element_location::ON_INTERFACE ) )
            return;

        auto asm_map = init_asm_map(msh, cl);
        auto dirichlet_data = get_dirichlet_data(msh, cl);
        assert( asm_map.size() == mass.rows() && asm_map.size() == mass.cols() );

        // MASS
        for (size_t i = 0; i < mass.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < mass.cols(); j++) {
                if ( asm_map[j].assemble() )
                    triplets_mass.push_back( Triplet<T>(asm_map[i], asm_map[j], mass(i,j)) );
            }
        }
    }
            
    void 
    assemble_grad_ext(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, 1>& contrib) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        auto loc = std::get<1>(P);

        // Cell offset
        size_t cell_offset = cell_table.at(offset(msh, cl)); 
        size_t cell_LHS_offset = cell_offset * contrib.rows();
        if (is_cut(msh, cl) && loc == element_location::IN_POSITIVE_SIDE)
            cell_LHS_offset += contrib.rows();
        GRAD.block(cell_LHS_offset, 0, contrib.rows(), 1) = contrib;

    }

    void 
    grad_contrib_assembly(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs) {
    
        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto loc = std::get<1>(P);
        auto cl = msh.cells[cell_index];

        // DOFS INFOS
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto graddeg = di.grad_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto gbs = vector_cell_basis<Mesh,T>::size(graddeg);

        // BLOCK INFOS
        auto CC_BLOCK_SIZE = this->num_cells*cbs;
        auto FF_BLOCK_SIZE = this->num_other_faces*fbs;

        auto asm_map = init_asm_map_ext(msh, P);

        // ASSEMBLY OF STIFFNESS MATRIX
        auto offset_i = cell_table.at(offset(msh, cl)); 
        if (is_cut(msh,cl) && loc==element_location::IN_POSITIVE_SIDE)
            offset_i++;
        for (size_t j = 0; j < lhs.cols(); j++) {
            if (!asm_map[j].assemble())
                continue;
            for (size_t i = 0; i < gbs; i++) {
                triplets_GRAD.push_back( Triplet<T>(offset_i*gbs+i, asm_map[j], lhs(i,j)) );
            }
        }
    }

    void 
    assemble_grad_grad_bis_extended(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];

        if( !(location(msh, cl) == loc_zone
              || location(msh, cl) == element_location::ON_INTERFACE
              || loc_zone == element_location::ON_INTERFACE ) )
            return;
        
        auto asm_map = init_asm_map_ext(msh, P);
        auto dirichlet_data = get_dirichlet_data_ext(msh, P);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());
        assert(dirichlet_data.size() == lhs.cols());

        // ASSEMBLY OF STIFFNESS MATRIX
        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) {
                    triplets_GRAD_GRAD.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                }
            }
        }
    }

    void assemble_conditioning(const Mesh& msh, Tuple P, double contrib) {

        // CELL INFOS
        auto cell_index = std::get<0>(P);
        auto cl = msh.cells[cell_index];
        auto loc = std::get<1>(P);
        
        // Cell offset
        size_t cell_offset = cell_table.at(offset(msh, cl)); 
        if (is_cut(msh, cl) && loc == element_location::IN_POSITIVE_SIDE)
            cell_offset += 1;
        CONDITIONING(cell_offset) = contrib;
        
    }

    Matrix<T, Dynamic, 1>
    get_solF(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution) {

        bool double_unknowns = ( location(msh, cl) == element_location::ON_INTERFACE && loc_zone == element_location::ON_INTERFACE );
        auto facdeg = di.face_degree();
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        size_t f_dofs = num_faces*fbs;
        if( double_unknowns )
            f_dofs = 2 * f_dofs;

        Matrix<T, Dynamic, 1> solF = Matrix<T, Dynamic, 1>::Zero( f_dofs );
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            if (loc_zone != element_location::ON_INTERFACE) {
                auto loc_fc = location(msh, fc);
                if (!(loc_fc == element_location::ON_INTERFACE || loc_fc == loc_zone) )
                    continue;
            }

            auto face_LHS_offset = face_SOL_offset(msh, fc);
            if (location(msh, fc) == element_location::ON_INTERFACE && loc_zone == element_location::ON_INTERFACE) {
                solF.block(face_i*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
                solF.block( (num_faces+face_i)*fbs, 0, fbs, 1)
                    = solution.block(face_LHS_offset + fbs, 0, fbs, 1);
                continue;
            }

            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
            if (dirichlet) {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, facdeg, dir_func);
                solF.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
                continue;
            }

            if( location(msh, cl) == element_location::ON_INTERFACE && location(msh, fc) == element_location::IN_POSITIVE_SIDE && loc_zone == element_location::ON_INTERFACE) {
                solF.block((num_faces+face_i)*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
                continue;
            }
            solF.block(face_i*fbs, 0, fbs, 1) = solution.block(face_LHS_offset, 0, fbs, 1);
        }

        return solF;

    }
         
    void 
    finalize(void) {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
        MASS.setFromTriplets( triplets_mass.begin(), triplets_mass.end() );
        triplets_mass.clear();
        Kg_ZIP.setFromTriplets( triplets_zip.begin(), triplets_zip.end() );
        triplets_zip.clear();
        SPARSITY.setFromTriplets( triplets_sparsity.begin(), triplets_sparsity.end() );
        triplets_sparsity.clear();
        GLOBAL_GRAD_GRAD.setFromTriplets( triplets_GRAD_GRAD.begin(), triplets_GRAD_GRAD.end() );
        triplets_GRAD_GRAD.clear();
    }

};

template<typename Mesh, typename Function>
class virt_interface_assembler : public virt_scalar_assembler<Mesh, Function> {

    using T = typename Mesh::coordinate_type;

public:

    virt_interface_assembler(const Mesh& msh, const Function& dirichlet_bf, hho_degree_info hdi) : virt_scalar_assembler<Mesh, Function>(msh, dirichlet_bf, hdi) {

        this->loc_zone = element_location::ON_INTERFACE;
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
        };
        size_t loc_num_cells = 0; 
        for (auto& cl : msh.cells) {
            this->cell_table.push_back( loc_num_cells );
            if (location(msh, cl) == element_location::ON_INTERFACE)
                loc_num_cells += 2;
            else
                loc_num_cells += 1;
        }
        this->num_cells = loc_num_cells;
        assert(this->cell_table.size() == msh.cells.size());

        size_t num_all_faces = 0; 
        for (auto& fc : msh.faces) {
            if (location(msh, fc) == element_location::ON_INTERFACE)
                num_all_faces += 2;
            else
                num_all_faces += 1;
        }
            
        size_t num_dirichlet_faces = 0; /* counts faces with dup. unknowns */
        for (auto& fc : msh.faces) {
            if(fc.is_boundary && fc.bndtype == boundary::DIRICHLET){
                if (location(msh, fc) == element_location::ON_INTERFACE)
                        num_dirichlet_faces += 2;
                    else
                        num_dirichlet_faces += 1;
            }
        }
        
        this->num_other_faces = num_all_faces - num_dirichlet_faces;
        this->face_table.resize( msh.faces.size() );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < msh.faces.size(); i++) {
            auto fc = msh.faces.at(i);
            if ( !is_dirichlet(fc) ) {
                this->face_table.at(i) = compressed_offset;
                if ( location(msh, fc) == element_location::ON_INTERFACE )
                    compressed_offset += 2;
                else
                    compressed_offset += 1;
            }
        }
    }
};

template<typename Mesh, typename Function>
class interface_assembler : public virt_interface_assembler<Mesh, Function> {

    using T = typename Mesh::coordinate_type;

public:

    interface_assembler(const Mesh& msh, const Function& dirichlet_bf, hho_degree_info hdi) : virt_interface_assembler<Mesh, Function>(msh, dirichlet_bf, hdi) {
        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();
        auto graddeg = this->di.grad_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto gbs = vector_cell_basis<Mesh,T>::size(graddeg);
        this->loc_cbs = cbs;
        auto system_size = cbs * this->num_cells + fbs * this->num_other_faces;
        this->LHS = SparseMatrix<T>( system_size, system_size );
        this->RHS = Matrix<T, Dynamic, 1>::Zero( system_size );        
        this->MASS = SparseMatrix<T>(system_size, system_size);
        // DEBUG
        this->loc_gbs = gbs;
        this->GRAD = Matrix<T, Dynamic, 1>::Zero(this->num_cells * gbs);
        this->GLOBAL_GRAD_GRAD = SparseMatrix<T>(system_size, system_size);
        this->SPARSITY = SparseMatrix<T>(system_size, system_size);
        this->Kg_ZIP = SparseMatrix<T>(this->num_cells + this->num_other_faces, this->num_cells + this->num_other_faces);

    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs) {
        this->assemble_bis(msh, cl, lhs, rhs);
    }

    void
    assemble_ext(const Mesh& msh, Tuple P, const Matrix<T, Dynamic, Dynamic>& lhs, 
    const Matrix<T, Dynamic, 1>& rhs) {

        this->assemble_bis_ext(msh, P, lhs, rhs);
    
    }
    
    Matrix<T, Dynamic, 1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution, element_location where) {

        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto cell_offset = offset(msh, cl);
        size_t cell_SOL_offset;
        if (location(msh, cl) == element_location::ON_INTERFACE) {
            if (where == element_location::IN_NEGATIVE_SIDE)
                cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
            else if (where == element_location::IN_POSITIVE_SIDE)
                cell_SOL_offset = this->cell_table.at(cell_offset) * cbs + cbs;
            else
                throw std::invalid_argument("Invalid location");
        }
        else 
            cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + num_faces*fbs);
        ret.block(0, 0, cbs, 1) = solution.block(cell_SOL_offset, 0, cbs, 1);
        auto solF = this->get_solF(msh, cl, solution);
        if (where == element_location::IN_NEGATIVE_SIDE)
            ret.tail(num_faces * fbs) = solF.head(num_faces * fbs);
        else
            ret.tail(num_faces * fbs) = solF.tail(num_faces * fbs);

        return ret;

    }
    
    Matrix<T, Dynamic, 1>
    gather_proj(const Mesh& msh, Tuple P, hho_degree_info hdi, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
        // CELL INFOS 
        auto cell_index = std::get<0>(P);
        auto loc = std::get<1>(P);
        auto cl = msh.cells[cell_index];

        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fbs_cut = 2*fbs;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto uncut_dofs = cbs + num_faces*fbs;
        auto current_dofs = uncut_dofs;
        if (is_cut(msh,cl)) 
            current_dofs = 2*current_dofs;
        auto extended_dofs = 2*(cbs + num_faces*fbs);
        auto dp_cells = std::get<2>(P);
        auto local_dofs = current_dofs + dp_cells.size()*extended_dofs; 

        Matrix<T, Dynamic, 1> dofs = Matrix<T, Dynamic, 1>::Zero(local_dofs);

        if (!is_cut(msh, cl)) 
            dofs.block(0, 0, cbs+num_faces*fbs, 1) = project_function(msh, cl, hdi, scal_fun);
        else {
            dofs.block(0, 0, cbs, 1) = project_function(msh, cl, hdi, element_location::IN_NEGATIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
            dofs.block(cbs, 0, cbs, 1) = project_function(msh, cl, hdi, element_location::IN_POSITIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
            dofs.block(2*cbs, 0, num_faces*fbs, 1) = project_function(msh, cl, hdi, element_location::IN_NEGATIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
            dofs.block(2*cbs+num_faces*fbs, 0, num_faces*fbs, 1) = project_function(msh, cl, hdi, element_location::IN_POSITIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
        }

        // LOOP OVER DEPENDENT CELLS 
        auto offset_dofs = current_dofs;  
        for (auto &dp_cl : dp_cells) {
            auto dp_cell = msh.cells[dp_cl];
            dofs.block(offset_dofs, 0, cbs, 1) = project_function(msh, dp_cell, hdi, element_location::IN_NEGATIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
            dofs.block(offset_dofs+cbs, 0, cbs, 1) = project_function(msh, dp_cell, hdi, element_location::IN_POSITIVE_SIDE, scal_fun).block(0, 0, cbs, 1);
            dofs.block(offset_dofs+2*cbs, 0, num_faces*fbs, 1) = project_function(msh, dp_cell, hdi, element_location::IN_NEGATIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
            dofs.block(offset_dofs+2*cbs+num_faces*fbs, 0, num_faces*fbs, 1) = project_function(msh, dp_cell, hdi, element_location::IN_POSITIVE_SIDE, scal_fun).block(cbs, 0, num_faces*fbs, 1);
            offset_dofs += extended_dofs;
        }

        return dofs;

    }

    Matrix<T, Dynamic, 1>
    gather_cell_dof(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& solution, element_location where) {
        
        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();

        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);

        auto cell_offset = offset(msh, cl);
        size_t cell_SOL_offset;
        if (location(msh, cl) == element_location::ON_INTERFACE) {
            if (where == element_location::IN_NEGATIVE_SIDE)
                cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
            else if (where == element_location::IN_POSITIVE_SIDE)
                cell_SOL_offset = this->cell_table.at(cell_offset) * cbs + cbs;
            else
                throw std::invalid_argument("Invalid location");
        }
        else {
            cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
        }
        return solution.block(cell_SOL_offset, 0, cbs, 1);
    }

    void project_over_cells_and_faces(const Mesh& msh, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {

        for (auto& cl : msh.cells) {
            if( location(msh, cl) != element_location::ON_INTERFACE ) 
                project_over_uncutcells(msh, cl, hho_di, x_glob, scal_fun); 
            else
                project_over_cutcells(msh, cl, hho_di, x_glob, scal_fun);
        }
    }

    void project_over_uncutcells(const Mesh& msh, const typename Mesh::cell_type& cl, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
        Matrix<T, Dynamic, 1> x_proj_dof = project_function(msh, cl, hho_di, scal_fun);

        // HHO DISCRETIZATION INFOS
        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
  
        // CELL DOFS 
        auto cell_offset = offset(msh, cl);
        size_t cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
        x_glob.block(cell_SOL_offset, 0, cbs, 1) = x_proj_dof.block(0, 0, cbs, 1);

        // FACE DOFS 
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = this->face_SOL_offset(msh, fc);
            x_glob.block(face_LHS_offset, 0, fbs, 1) = x_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
        }

    }
            
    void project_over_cutcells(const Mesh& msh, const typename Mesh::cell_type& cl, hho_degree_info hho_di, Matrix<T, Dynamic, 1> & x_glob, std::function<T(const typename Mesh::point_type& )> scal_fun) {
            
        Matrix<T, Dynamic, 1> x_neg_proj_dof = project_function(msh, cl, hho_di, element_location::IN_NEGATIVE_SIDE, scal_fun);
        Matrix<T, Dynamic, 1> x_pos_proj_dof = project_function(msh, cl, hho_di, element_location::IN_POSITIVE_SIDE, scal_fun);
            
        
        
        auto celdeg = this->di.cell_degree();
        auto facdeg = this->di.face_degree();
        auto cbs = cell_basis<Mesh,T>::size(celdeg);
        auto fbs = face_basis<Mesh,T>::size(facdeg); 
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto cell_offset = offset(msh, cl);
        size_t cell_SOL_offset = this->cell_table.at(cell_offset) * cbs;
        x_glob.block(cell_SOL_offset, 0, cbs, 1) = x_neg_proj_dof.block(0, 0, cbs, 1);
        x_glob.block(cell_SOL_offset + cbs, 0, cbs, 1) = x_pos_proj_dof.block(0, 0, cbs, 1);
        
        // FACE DOFS 
        for (size_t face_i = 0; face_i < num_faces; face_i++) {
            auto fc = fcs[face_i];
            auto face_LHS_offset = this->face_SOL_offset(msh, fc);
            x_glob.block(face_LHS_offset, 0, fbs, 1) = x_neg_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
            x_glob.block(face_LHS_offset + fbs, 0, fbs, 1) = x_pos_proj_dof.block(cbs + face_i*fbs, 0, fbs, 1);
        }

    }
    
    void
    assemble_mass(const Mesh& msh, const typename Mesh::cell_type& cl, 
                  const Matrix<T, Dynamic, Dynamic>& mass) {
        
        this->assemble_bis_mass(msh, cl, mass);
    
    }

    void
    assemble_rhs(const Mesh& msh, const typename Mesh::cell_type& cl, const Matrix<T, Dynamic, 1>& rhs) {
        this->assemble_rhs_bis(msh, cl, rhs);
    }

};
            
template<typename Mesh, typename Function>
auto make_interface_assembler(const Mesh& msh, Function dirichlet_bf, hho_degree_info hdi) {
    return interface_assembler<Mesh, Function>(msh, dirichlet_bf, hdi);
}
