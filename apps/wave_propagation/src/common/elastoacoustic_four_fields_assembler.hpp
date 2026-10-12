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
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/elastic_material_data.hpp"
#include "../common/elastoacoustic_assembler_base.hpp"
#include <map>

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class elastoacoustic_four_fields_assembler : public elastoacoustic_assembler_base<Mesh> {
    
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

    std::vector< Triplet<T> >           m_triplets_stab;
    bool        m_scaled_stabilization_Q;
    T           m_e_eta = 1.50;     // stabilization weights eta (elastic, acoustic): MEKG values by default (MEKG uses
    T           m_a_eta = 0.80;     //   0.95 / 0.4 in mixed order, not used here); 1 in Burman-Duran-Ern(-Steins)
    int         m_stab_size = 0;    // length of the 1/h scaling, see set_stabilization_size


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
    }
    
    // Assembly of acoustic and elastic part separately; 
    void assemble(const Mesh& msh, std::function<disk::static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun, std::function<T(const typename Mesh::point_type& )> a_rhs_fun, bool explicit_scheme) {
        
        auto storage = msh.backend_storage();
        LHS.setZero();
        RHS.setZero();
        
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
            auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, m_scaled_stabilization_Q);
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
                auto stabilization_operator = make_vector_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }
        
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
            auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, m_scaled_stabilization_Q);
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
                auto stabilization_operator = make_scalar_hdg_stabilization(msh, cell, m_hho_di, m_stab_size == 2, m_scaled_stabilization_Q);
                auto n_s_rows = stabilization_operator.rows();
                auto n_s_cols = stabilization_operator.cols();
                S_operator.block(n_rows-n_s_rows, n_cols-n_s_cols, n_s_rows, n_s_cols) = stabilization_operator;
            }
        }

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
    }   
    
    void set_scaled_stabilization() {
        m_scaled_stabilization_Q = true;
    }

    // Length h of the 1/h scaling: 0 as in diskpp (here: cell diameter h_T, except for the vector HHO stabilization:
    // face diameter h_F), 1 cell diameter h_T, 2 face diameter h_F
    void set_stabilization_size(int size) {
        m_stab_size = size;
    }

    // Weights eta of the elastic and acoustic stabilizations (default 1.5 and 0.8)
    void set_stabilization_weights(T e_eta, T a_eta) {
        m_e_eta = e_eta;
        m_a_eta = a_eta;
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

};

#endif /* elastoacoustic_four_fields_assembler_hpp */
