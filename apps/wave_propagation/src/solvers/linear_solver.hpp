//
//  linear_solver.hpp
//  acoustics
//
//  Created by Omar Durán on 5/2/20.
//

#pragma once
#ifndef linear_solver_hpp
#define linear_solver_hpp

#include "diskpp/solvers/solver.hpp"
#include "parallel_sparse.hpp"

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename T>
class linear_solver
{
    
    private:

    SparseMatrix<T> m_Kcc;
    SparseMatrix<T> m_Kcf;
    SparseMatrix<T> m_Kfc;
    SparseMatrix<T> m_Kff;
    Matrix<T, Dynamic, 1> m_Fc;
    Matrix<T, Dynamic, 1> m_Ff;
    
    SparseMatrix<T> m_Kcc_inv;
    SparseMatrix<T> m_K;
    Matrix<T, Dynamic, 1> m_F;
    
    #ifdef HAVE_INTEL_MKL
        PardisoLU<Eigen::SparseMatrix<T>>  m_analysis;
    #else
        SparseLU<Eigen::SparseMatrix<T>> m_analysis;
    #endif
    
    #ifdef HAVE_INTEL_MKL
        PardisoLDLT<Eigen::SparseMatrix<T>>  m_symm_analysis;
    #else
        SimplicialLDLT<Eigen::SparseMatrix<T>> m_symm_analysis;
    #endif
    
    ConjugateGradient<SparseMatrix<T>> m_analysis_cg;
    BiCGSTAB<SparseMatrix<T>,IncompleteLUT<double>> m_analysis_bi_cg;

    size_t m_n_c_dof;
    size_t m_n_f_dof;
    bool m_global_sc_Q;
    bool m_is_decomposed_Q; 
    std::pair<bool,bool> m_iterative_solver_Q;
    T m_tolerance;
    
    
    void DecomposeK(){
        if (!m_iterative_solver_Q.first) {
            if (m_iterative_solver_Q.second) {
                m_symm_analysis.analyzePattern(m_K);
                m_symm_analysis.factorize(m_K);
            }
            else{ 
                #ifdef HAVE_INTEL_MKL
                if (m_fast) {
                    m_analysis.pardisoParameterArray()[7] = 0;     // no iterative refinement (3x faster solves here)
                    m_analysis.pardisoParameterArray()[24] = 2;    // parallel forward / backward substitution
                }
                #endif
                m_analysis.analyzePattern(m_K);
                m_analysis.factorize(m_K);
            }
        }
    }
    
    Matrix<T, Dynamic, 1> solve_global(Matrix<T, Dynamic, 1> & Fg){
        if (m_iterative_solver_Q.first) {
            if (m_iterative_solver_Q.second) {
                Matrix<T, Dynamic, 1> x_dof = m_analysis_cg.solve(Fg);
                std::cout << "Number of iterations (CG): " << m_analysis_cg.iterations() << std::endl;
                std::cout << "Estimated error: " << m_analysis_cg.error() << std::endl;
                return x_dof;
            }
            else {
                Matrix<T, Dynamic, 1> x_dof = m_analysis_bi_cg.solve(Fg);
                std::cout << "Number of iterations (BiCG): " << m_analysis_bi_cg.iterations() << std::endl;
                std::cout << "Estimated error: " << m_analysis_bi_cg.error() << std::endl;
                return x_dof;
            }
        }
        else {
            if (m_iterative_solver_Q.second) {
                return m_symm_analysis.solve(Fg);
            }else{
                return m_analysis.solve(Fg);
            }
        }
    }
    
    void scatter_segments(Matrix<T, Dynamic, 1> & Fg){
        assert(m_n_c_dof + m_n_f_dof == Fg.rows());
        m_Fc = Fg.block(0, 0, m_n_c_dof, 1);
        m_Ff = Fg.block(m_n_c_dof, 0, m_n_f_dof, 1);
    }
    
    Matrix<T, Dynamic, 1> solve_sc(Matrix<T, Dynamic, 1> & Fg){
        
        // Split vector into cells and faces dofs
        scatter_segments(Fg);
        
        // Condense vector
        Matrix<T, Dynamic, 1> delta_c;
        if (m_fast) {
            delta_c = m_Kcc_inv_rm*m_Fc;
            m_F = m_Ff - m_Kfc_rm*delta_c;
        }
        else {
            delta_c = m_Kcc_inv*m_Fc;
            m_F = m_Ff - m_Kfc*delta_c;
        }
        
        // face linear solver step
        Matrix<T, Dynamic, 1> x_n_f_dof;
        if (m_iterative_solver_Q.first) {
            if (m_iterative_solver_Q.second) {
                x_n_f_dof = m_analysis_cg.solve(m_F);
                std::cout << "Number of iterations (CG): " << m_analysis_cg.iterations() << std::endl;
                std::cout << "Estimated error: " << m_analysis_cg.error() << std::endl;
            }
            else {
                // std::cout << "HERE" << std::endl;
                x_n_f_dof = m_analysis_bi_cg.solve(m_F);
                std::cout << bold << cyan << "      Iterative solver for stage s: " << reset << std::endl;
                std::cout << yellow << bold << "         Number of iterations (BiCG): " << m_analysis_bi_cg.iterations() << std::endl;
                std::cout << yellow << bold << "         Estimated error: " << m_analysis_bi_cg.error() << std::endl;
            }
        }
        else {
            if (m_iterative_solver_Q.second) {
                x_n_f_dof = m_symm_analysis.solve(m_F);
            }
            else {
                x_n_f_dof = m_analysis.solve(m_F);
            }
        }
        
        
        Matrix<T, Dynamic, 1> Kcf_x_f_dof = m_fast ? Matrix<T, Dynamic, 1>(m_Kcf_rm*x_n_f_dof) : Matrix<T, Dynamic, 1>(m_Kcf*x_n_f_dof);
        Matrix<T, Dynamic, 1> delta_f = m_fast ? Matrix<T, Dynamic, 1>(m_Kcc_inv_rm*Kcf_x_f_dof) : Matrix<T, Dynamic, 1>(m_Kcc_inv*Kcf_x_f_dof);
        Matrix<T, Dynamic, 1> x_n_c_dof = delta_c - delta_f;
        
        // Composing global solution
        Matrix<T, Dynamic, 1> x_dof = Matrix<T, Dynamic, 1>::Zero(m_n_c_dof+m_n_f_dof,1);
        x_dof.block(0, 0, m_n_c_dof, 1) = x_n_c_dof;
        x_dof.block(m_n_c_dof, 0, m_n_f_dof, 1) = x_n_f_dof;
        return x_dof;

    }
    
    void scatter_blocks(SparseMatrix<T> & Kg){
        
        m_n_c_dof = Kg.rows() - m_n_f_dof;
        
        // scattering matrix blocks
        m_Kcc = Kg.block(0, 0, m_n_c_dof, m_n_c_dof);
        m_Kcf = Kg.block(0, m_n_c_dof, m_n_c_dof, m_n_f_dof);
        m_Kfc = Kg.block(m_n_c_dof,0, m_n_f_dof, m_n_c_dof);
        m_Kff = Kg.block(m_n_c_dof,m_n_c_dof, m_n_f_dof, m_n_f_dof);
        m_is_decomposed_Q = false;
    }
    
    public:
    
    linear_solver() : m_global_sc_Q(false), m_is_decomposed_Q(false)  {
        m_iterative_solver_Q = std::make_pair(false,false);
    }
    
    linear_solver(SparseMatrix<T> & Kg) : m_K(Kg), m_global_sc_Q(false), m_is_decomposed_Q(false)  {
        m_iterative_solver_Q = std::make_pair(false,false);
    }
    
    linear_solver(SparseMatrix<T> & Kg, size_t n_f_dof) :
        m_n_f_dof(n_f_dof),
        m_global_sc_Q(true),
        m_is_decomposed_Q(false) {
        scatter_blocks(Kg);
    }
    
    void set_Kg(SparseMatrix<T> & Kg){
        m_global_sc_Q = false;
        m_K = Kg;
    }

    void set_Kg(SparseMatrix<T> & Kg, size_t n_f_dof){
        m_global_sc_Q = true;
        m_n_f_dof = n_f_dof;
        scatter_blocks(Kg);
    }
    
    // Fast mode (solver.matrix_free with sdirk): Pardiso without iterative refinement and with parallel substitutions,
    // products of the static condensation on row-major copies (parallel); results equal up to round-off
    void set_fast(bool fast) {
        m_fast = fast;
    }

    void set_direct_solver(bool symmetric_matrix_Q = false){
        m_iterative_solver_Q = std::make_pair(false,symmetric_matrix_Q);
    }
    
    void set_iterative_solver(bool symmetric_matrix_Q = false, T tolerance = 1.0e-11){
        m_iterative_solver_Q = std::make_pair(true, symmetric_matrix_Q);
        m_is_decomposed_Q = true;
        m_tolerance = tolerance;
    }
    
    SparseMatrix<T> & Kcc(){
        return m_Kcc;
    }
    
    SparseMatrix<T> & Kcf(){
        return m_Kcf;
    }
    
    SparseMatrix<T> & Kfc(){
        return m_Kfc;
    }

    Matrix<T, Dynamic, 1> & Fc(){
        return m_Fc;
    }
        
    void condense_equations(std::vector<std::pair<size_t,size_t>> vec_cell_basis_data){
        
        if (!m_global_sc_Q) {
            return;
        }
        if (m_fast) {
            condense_equations_fast(vec_cell_basis_data);
            return;
        }
        
        size_t nnz_cc = 0;
        for (auto chunk : vec_cell_basis_data) {
            nnz_cc += chunk.second*chunk.second*chunk.first;
        }
        std::vector< Triplet<T> > triplets_cc;
        triplets_cc.resize(nnz_cc);
        m_Kcc_inv = SparseMatrix<T>( m_n_c_dof, m_n_c_dof );
        #ifdef HAVE_INTEL_TBB
        size_t stride_n_block_eq = 0;
        size_t stride_n_block_l = 0;
        for (auto chunk : vec_cell_basis_data) {
                size_t n_cells = chunk.first;
                size_t n_cbs   = chunk.second;
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [this,&triplets_cc,&n_cbs,&stride_n_block_eq,&stride_n_block_l] (size_t & cell_ind){
                    
                    size_t stride_eq = cell_ind * n_cbs + stride_n_block_eq;
                    size_t stride_l = cell_ind * n_cbs * n_cbs + stride_n_block_l;

                    SparseMatrix<T> K_cc_loc = m_Kcc.block(stride_eq, stride_eq, n_cbs, n_cbs);
                    SparseLU<SparseMatrix<T>> analysis_cc;
                    analysis_cc.analyzePattern(K_cc_loc);
                    analysis_cc.factorize(K_cc_loc);
                    Matrix<T, Dynamic, Dynamic> K_cc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_cbs, n_cbs));
            
                    size_t l = 0;
                    for (size_t i = 0; i < K_cc_inv_loc.rows(); i++)
                    {
                        for (size_t j = 0; j < K_cc_inv_loc.cols(); j++)
                        {
                            triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, K_cc_inv_loc(i,j));
                            l++;
                        }
                    }
                }
            );
            stride_n_block_eq   += n_cells * n_cbs;
            stride_n_block_l    += n_cells * n_cbs * n_cbs;
        }
        #else

        size_t stride_n_block_eq = 0;
        size_t stride_n_block_l = 0;
        for (auto chunk : vec_cell_basis_data) {
            size_t n_cells = chunk.first;
            size_t n_cbs   = chunk.second;
            for (size_t cell_ind = 0; cell_ind < n_cells; cell_ind++) {

                size_t stride_eq = cell_ind * n_cbs + stride_n_block_eq;
                size_t stride_l = cell_ind * n_cbs * n_cbs + stride_n_block_l;
                
                SparseMatrix<T> K_cc_loc = m_Kcc.block(stride_eq, stride_eq, n_cbs, n_cbs);
                SparseLU<SparseMatrix<T>> analysis_cc;
                analysis_cc.analyzePattern(K_cc_loc);
                analysis_cc.factorize(K_cc_loc);
                Matrix<T, Dynamic, Dynamic> K_cc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(n_cbs, n_cbs));
        
                size_t l = 0;
                for (size_t i = 0; i < K_cc_inv_loc.rows(); i++) {
                    for (size_t j = 0; j < K_cc_inv_loc.cols(); j++) {
                        triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, K_cc_inv_loc(i,j));
                        l++;
                    }
                }

            }
            stride_n_block_eq   += n_cells * n_cbs;
            stride_n_block_l    += n_cells * n_cbs * n_cbs;
        }
        #endif
        
        m_Kcc_inv.setFromTriplets( triplets_cc.begin(), triplets_cc.end() );
        triplets_cc.clear();
        m_K = m_Kff - m_Kfc*m_Kcc_inv*m_Kcf;
        if (m_fast) {
            // row-major copies for the products of the stages (OpenMP in Eigen); the column-major ones released
            m_Kcc_inv_rm = m_Kcc_inv; m_Kfc_rm = m_Kfc; m_Kcf_rm = m_Kcf;
            m_Kcc_inv = SparseMatrix<T>(); m_Kfc = SparseMatrix<T>(); m_Kcf = SparseMatrix<T>();
        }
        m_is_decomposed_Q = false;
        return;

    }
    
    // After condense_equations: the cell block (inverted in m_Kcc_inv) and the face block (in the condensed m_K) are
    // not used by solve_sc
    void release_condensed_blocks(){
        m_Kcc = SparseMatrix<T>();
        m_Kff = SparseMatrix<T>();
    }

    void condense_equations(std::pair<size_t,size_t> cell_basis_data){
        std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(1);
        vec_cell_basis_data[0] = cell_basis_data;
        condense_equations(vec_cell_basis_data);
    }
        
    void factorize(){
        if (m_is_decomposed_Q) {
            if (m_iterative_solver_Q.first) {
                if (m_iterative_solver_Q.second) {
                    m_analysis_cg.compute(m_K);
                    m_analysis_cg.setTolerance(m_tolerance);
                    m_analysis_cg.setMaxIterations(m_K.rows());
                }else{
                    m_analysis_bi_cg.compute(m_K);
                    m_analysis_bi_cg.setTolerance(m_tolerance);
                    m_analysis_bi_cg.setMaxIterations(m_K.rows());
                }
            }
            return;
        }
        DecomposeK(); 
        m_is_decomposed_Q = true;
    }
    
    size_t n_equations(){
        size_t n_equations = m_K.rows();
        return n_equations;
    }

    Matrix<T, Dynamic, 1> solve(Matrix<T, Dynamic, 1> & Fg){
        if (m_global_sc_Q) 
            return solve_sc(Fg);
        else 
            return solve_global(Fg);
    }
    


    private:
    // condense_equations of the fast mode: dense inverses of the cell blocks and local Schur complements
    // K_fT K_TT^-1 K_Tf of the cells (on the faces of the cell) computed by the OpenMP threads, summed in the order of
    // the cells (equal to the products of condense_equations up to round-off); row-major copies for solve_sc
    void condense_equations_fast(const std::vector<std::pair<size_t,size_t>>& vec_cell_basis_data) {
        typedef Matrix<T, Dynamic, Dynamic> dense_type;
        std::vector<std::pair<size_t,size_t>> cells;    // (first unknown, size)
        size_t r = 0;
        for (auto chunk : vec_cell_basis_data)
            for (size_t i = 0; i < chunk.first; i++) { cells.push_back({r, chunk.second}); r += chunk.second; }
        Eigen::SparseMatrix<T, Eigen::RowMajor> Kcf_rm = m_Kcf;
        std::vector<std::vector<Triplet<T>>> inv_t(cells.size()), schur_t(cells.size());
        #pragma omp parallel for schedule(dynamic, 16)
        for (long c = 0; c < (long)cells.size(); c++) {
            const size_t r0 = cells[c].first, n = cells[c].second;
            const dense_type inv = dense_type(m_Kcc.block(r0, r0, n, n)).partialPivLu().inverse();
            // faces of the cell: columns of its rows of Kcf
            std::vector<int> J;
            for (size_t i = r0; i < r0 + n; i++)
                for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(Kcf_rm, Eigen::Index(i)); it; ++it)
                    J.push_back(int(it.col()));
            std::sort(J.begin(), J.end());
            J.erase(std::unique(J.begin(), J.end()), J.end());
            dense_type Bcf = dense_type::Zero(n, J.size()), Bfc = dense_type::Zero(J.size(), n);
            for (size_t i = r0; i < r0 + n; i++)
                for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(Kcf_rm, Eigen::Index(i)); it; ++it)
                    Bcf(i - r0, std::lower_bound(J.begin(), J.end(), int(it.col())) - J.begin()) = it.value();
            for (size_t j = r0; j < r0 + n; j++)
                for (typename SparseMatrix<T>::InnerIterator it(m_Kfc, Eigen::Index(j)); it; ++it) {
                    auto p = std::lower_bound(J.begin(), J.end(), int(it.row()));
                    if (p != J.end() && *p == int(it.row())) Bfc(p - J.begin(), j - r0) = it.value();
                }
            const dense_type S = Bfc * (inv * Bcf);
            inv_t[c].reserve(n * n);
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++) inv_t[c].emplace_back(r0 + i, r0 + j, inv(i, j));
            schur_t[c].reserve(J.size() * J.size());
            for (size_t a = 0; a < J.size(); a++)
                for (size_t b = 0; b < J.size(); b++) schur_t[c].emplace_back(J[a], J[b], -S(a, b));
        }
        std::vector<Triplet<T>> t;
        for (auto& v : inv_t) t.insert(t.end(), v.begin(), v.end());
        std::vector<std::vector<Triplet<T>>>().swap(inv_t);
        m_Kcc_inv = SparseMatrix<T>(m_n_c_dof, m_n_c_dof);
        parallel_sparse::set_from_triplets(m_Kcc_inv, t);
        t.clear();
        for (int k = 0; k < m_Kff.outerSize(); ++k)
            for (typename SparseMatrix<T>::InnerIterator it(m_Kff, k); it; ++it) t.emplace_back(it.row(), it.col(), it.value());
        for (auto& v : schur_t) t.insert(t.end(), v.begin(), v.end());
        std::vector<std::vector<Triplet<T>>>().swap(schur_t);
        m_K = SparseMatrix<T>(m_n_f_dof, m_n_f_dof);
        parallel_sparse::set_from_triplets(m_K, t);
        m_Kcc_inv_rm = m_Kcc_inv; m_Kfc_rm = m_Kfc; m_Kcf_rm = Kcf_rm;
        m_Kcc_inv = SparseMatrix<T>(); m_Kfc = SparseMatrix<T>(); m_Kcf = SparseMatrix<T>();
        m_is_decomposed_Q = false;
    }

    // fast mode (set_fast); declared last
    bool m_fast = false;
    Eigen::SparseMatrix<T, Eigen::RowMajor> m_Kcc_inv_rm, m_Kfc_rm, m_Kcf_rm;
};

#endif /* linear_solver_hpp */
