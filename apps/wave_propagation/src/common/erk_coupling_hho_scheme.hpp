#ifndef erk_coupling_hho_scheme_hpp
#define erk_coupling_hho_scheme_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "../common/assembly_index.hpp"
#include <iomanip>  // Pour std::setw


template<typename T>
class erk_coupling_hho_scheme {
    
    private:
    
    SparseMatrix<T> m_Mc;
    SparseMatrix<T> m_Kcc;
    SparseMatrix<T> m_Kcf;
    SparseMatrix<T> m_Kfc;
    SparseMatrix<T> m_Sff;
    SparseMatrix<T> m_Scc;
    SparseMatrix<T> m_Mc_inv;
    SparseMatrix<T> m_Sff_inv;
    SparseMatrix<T> m_inv_Sff;
    SparseMatrix<T> m_Cg;
    
    std::vector<int> m_coarse_active_dofs; 
    std::vector<int> m_coarse_active_c;      // active cell indices (local to [0, n_c_dof))
    std::vector<int> m_coarse_active_f;      // active face indices (local to [0, n_f_dof)
    std::vector<int> m_fine_active_dofs;   // active global indices of Pfine
    std::vector<int> m_fine_active_c;      // active cell indices (local to [0, n_c_dof))
    std::vector<int> m_fine_active_f;      // active face indices (local to [0, n_f_dof))
    SparseMatrix<T>  m_Kcc_fine;           // Kcc restricted to active cell dofs
    SparseMatrix<T>  m_Kcf_fine;           // Kcf restricted to active cell x active face dofs
    SparseMatrix<T>  m_Kfc_fine;           // Kfc restricted to active face x active cell dofs
    SparseMatrix<T>  m_Mc_inv_fine;        // Mc_inv restricted to active cell dofs
    SparseMatrix<T>  m_Sff_inv_fine;       // Sff_inv restricted to active face dofs
    SparseMatrix<T>  m_Kcc_coarse;           // Kcc restricted to active cell dofs
    SparseMatrix<T>  m_Kcf_coarse;           // Kcf restricted to active cell x active face dofs
    SparseMatrix<T>  m_Kfc_coarse;           // Kfc restricted to active face x active cell dofs
    SparseMatrix<T>  m_Mc_inv_coarse;        // Mc_inv restricted to active cell dofs
    SparseMatrix<T>  m_Sff_inv_coarse;       // Sff_inv restricted to active face dofs

    Matrix<T, Dynamic, 1> m_Fc;
    
    #ifdef HAVE_INTEL_MKL
    PardisoLDLT<SparseMatrix<T>>  m_analysis_f;
    #else
    SimplicialLDLT<SparseMatrix<T>> m_analysis_f;
    #endif
    
    ConjugateGradient<SparseMatrix<T>> m_analysis_cg;
    
    size_t m_n_ec_dof;
    size_t m_n_ac_dof;
    size_t m_n_c_dof;
    size_t m_n_ef_dof;
    size_t m_n_af_dof;
    size_t m_n_f_dof;
    
    bool m_sff_is_block_diagonal_Q;
    bool m_iterative_solver_Q;
    
    public:
    
    erk_coupling_hho_scheme(SparseMatrix<T> & Kg, Matrix<T, Dynamic, 1> & Fg, SparseMatrix<T> & Mg, SparseMatrix<T> & Cg, size_t elastic_cell_dofs, size_t acoustic_cell_dofs, size_t e_face_dofs, size_t a_face_dofs) {
        
        m_n_ec_dof = elastic_cell_dofs;
        m_n_ac_dof = acoustic_cell_dofs;
        // m_n_c_dof  = m_n_ec_dof + m_n_ac_dof;
        m_n_ef_dof = e_face_dofs;
        m_n_af_dof = a_face_dofs;
        m_n_f_dof  = e_face_dofs + a_face_dofs;
        m_n_c_dof  = Kg.rows() - m_n_f_dof;
        
        m_Mc  = Mg.block(        0,         0, m_n_c_dof, m_n_c_dof);
        m_Kcc = Kg.block(        0,         0, m_n_c_dof, m_n_c_dof);
        m_Kcf = Kg.block(        0, m_n_c_dof, m_n_c_dof, m_n_f_dof);
        m_Kfc = Kg.block(m_n_c_dof,         0, m_n_f_dof, m_n_c_dof);
        m_Sff = Kg.block(m_n_c_dof, m_n_c_dof, m_n_f_dof, m_n_f_dof);
        
        m_Cg  = Cg.block(m_n_c_dof, m_n_c_dof, m_n_f_dof, m_n_f_dof);
        
        m_Fc  = Fg.block(0, 0, m_n_c_dof, 1);
        
        m_sff_is_block_diagonal_Q   = true;
        m_iterative_solver_Q        = false;
        
    }
    
    void Mcc_inverse(size_t e_cells, size_t a_cells, size_t e_cbs, size_t a_cbs) {
        
        size_t nnz_cc  = e_cbs*e_cbs*e_cells + a_cbs*a_cbs*a_cells;
        std::vector< Triplet<T> > triplets_cc;
        triplets_cc.resize(nnz_cc);
        m_Mc_inv = SparseMatrix<T>(m_n_c_dof, m_n_c_dof);
        
        #ifdef HAVE_INTEL_TBB
        
        tbb::parallel_for(size_t(0), size_t(e_cells), size_t(1), [this,&triplets_cc,&e_cbs] (size_t & cell_ind) {
            
            size_t stride_eq = cell_ind * e_cbs;
            size_t stride_l  = cell_ind * e_cbs * e_cbs;
            
            SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, e_cbs, e_cbs);
            SparseLU<SparseMatrix<T>> analysis_cc;
            analysis_cc.analyzePattern(m_Mc_loc);
            analysis_cc.factorize(m_Mc_loc);
            Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_cbs, e_cbs));
            
            size_t l = 0;
            for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                }
            }
        });
        tbb::parallel_for(size_t(0), size_t(a_cells), size_t(1), [this,&triplets_cc,&a_cbs,&e_cbs,&e_cells] (size_t & cell_ind) {
            
            size_t stride_eq = cell_ind*a_cbs       + e_cells*e_cbs;
            size_t stride_l  = cell_ind*a_cbs*a_cbs + e_cells*e_cbs*e_cbs;
            
            SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, a_cbs, a_cbs);
            SparseLU<SparseMatrix<T>> analysis_cc;
            analysis_cc.analyzePattern(m_Mc_loc);
            analysis_cc.factorize(m_Mc_loc);
            Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_cbs, a_cbs));
            
            size_t l = 0;
            for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                }
            }
        });
        
        #else
        
        for (size_t cell_ind = 0; cell_ind < e_cells; cell_ind++) {
            size_t stride_eq = cell_ind * e_cbs;
            size_t stride_l  = cell_ind * e_cbs * e_cbs;
            
            SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, e_cbs, e_cbs);
            SparseLU<SparseMatrix<T>> analysis_cc;
            analysis_cc.analyzePattern(m_Mc_loc);
            analysis_cc.factorize(m_Mc_loc);
            Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_cbs, e_cbs));
            
            size_t l = 0;
            for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                }
            }
        }
        for (size_t cell_ind = 0; cell_ind < a_cells; cell_ind++) {
            size_t stride_eq = cell_ind*a_cbs       + e_cells*e_cbs;
            size_t stride_l  = cell_ind*a_cbs*a_cbs + e_cells*e_cbs*e_cbs;
            
            SparseMatrix<T> m_Mc_loc = m_Mc.block(stride_eq, stride_eq, a_cbs, a_cbs);
            SparseLU<SparseMatrix<T>> analysis_cc;
            analysis_cc.analyzePattern(m_Mc_loc);
            analysis_cc.factorize(m_Mc_loc);
            Matrix<T, Dynamic, Dynamic> m_Mc_inv_loc = analysis_cc.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_cbs, a_cbs));
            
            size_t l = 0;
            for (size_t i = 0; i < m_Mc_inv_loc.rows(); i++) {
                for (size_t j = 0; j < m_Mc_inv_loc.cols(); j++) {
                    triplets_cc[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, m_Mc_inv_loc(i,j));
                    l++;
                }
            }
        }
        
        #endif
        
        m_Mc_inv.setFromTriplets(triplets_cc.begin(), triplets_cc.end());
        triplets_cc.clear();
        return;
        
    } 
    
    
    void Sff_inverse(size_t e_faces, size_t a_faces, size_t e_fbs, size_t a_fbs, std::vector<size_t> e_compress, std::vector<size_t> a_compress, std::set<size_t> elastic_internal_faces, std::set<size_t> acoustic_internal_faces, std::set<size_t> interfaces_index) {
        
        size_t n_interfaces = interfaces_index.size();                                          // Number of interfaces
        size_t nnz_ff = e_fbs*e_fbs*e_faces + a_fbs*a_fbs*a_faces + 2*e_fbs*a_fbs*n_interfaces; // Number of nonzeros
        std::vector< Triplet<T> > triplets_ff;
        triplets_ff.resize(nnz_ff);
        m_Sff_inv = SparseMatrix<T>(m_n_f_dof, m_n_f_dof);                                      // size: number of faces x number of faces
        
        // Inversion of elastic stabilization 
        for (size_t face_ind = 0; face_ind < e_faces; face_ind++) {
            // std::cout << "Elastic face: " << face_ind << std::endl << std::endl;
            size_t stride_eq = face_ind * e_fbs;
            size_t stride_l  = face_ind * e_fbs * e_fbs;
            SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, e_fbs, e_fbs);
            SparseLU<SparseMatrix<T>> analysis_ff;
            analysis_ff.analyzePattern(S_ff_loc);
            analysis_ff.factorize(S_ff_loc);
            Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(e_fbs, e_fbs));
            size_t l = 0;
            for (size_t i = 0; i < S_ff_inv_loc.rows(); i++) {
                for (size_t j = 0; j < S_ff_inv_loc.cols(); j++) {
                    triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
                    l++;
                }
            }
        }
        
        // Inversion of acoustic stabilization 
        for (size_t face_ind = 0; face_ind < a_faces; face_ind++) {
            // std::cout << "Acoutic face: " << e_faces + face_ind  << std::endl << std::endl; 
            size_t stride_eq = e_faces*e_fbs       + face_ind*a_fbs ;
            size_t stride_l  = e_faces*e_fbs*e_fbs + face_ind*a_fbs*a_fbs;   
            SparseMatrix<T> S_ff_loc = m_Sff.block(stride_eq, stride_eq, a_fbs, a_fbs);
            SparseLU<SparseMatrix<T>> analysis_ff;
            analysis_ff.analyzePattern(S_ff_loc);
            analysis_ff.factorize(S_ff_loc);
            Matrix<T, Dynamic, Dynamic> S_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(a_fbs, a_fbs));  
            size_t l = 0;
            for (size_t i = 0; i < S_ff_inv_loc.rows(); i++) {
                for (size_t j = 0; j < S_ff_inv_loc.cols(); j++) {
                    triplets_ff[stride_l+l] = Triplet<T>(stride_eq+i, stride_eq+j, S_ff_inv_loc(i,j));
                    l++;
                }
            }
        }
        
        // Inversion of coupling terms 
        size_t cpt = 0;
        for (auto face : interfaces_index) {                                     // Parcours des interfaces
            size_t e_face_LHS_offset = e_compress.at(face)*e_fbs;                  // Indice de la face elastique
            size_t a_face_LHS_offset = e_faces*e_fbs + a_compress.at(face)*a_fbs;  // Indice de la face acoustique
            // std::cout << "Interface: " << face  << std::endl;
            // std::cout << "Elastic interface: "  << e_compress.at(face) << std::endl;
            // std::cout << "Acoustic interface: " << e_faces + a_compress.at(face) << std::endl << std::endl; 
            size_t fbs = e_fbs + a_fbs;
            size_t e_stride_l = e_compress.at(face)*e_fbs*e_fbs;
            size_t a_stride_l = e_faces*e_fbs*e_fbs + a_compress.at(face)*a_fbs*a_fbs;
            size_t i_stride_l = e_faces*e_fbs*e_fbs + a_faces*a_fbs*a_fbs + 2*cpt*a_fbs*e_fbs;
            
            // Extraction du bloc stabilisation local
            Matrix<T, Dynamic, Dynamic> dense_SC_ff(fbs, fbs);
            SparseMatrix<T> elastic_stab  = m_Sff.block(e_face_LHS_offset, e_face_LHS_offset, e_fbs, e_fbs);
            SparseMatrix<T> acoustic_stab = m_Sff.block(a_face_LHS_offset, a_face_LHS_offset, a_fbs, a_fbs);
            dense_SC_ff.block(0, 0, e_fbs, e_fbs)         = elastic_stab;
            dense_SC_ff.block(e_fbs, e_fbs, a_fbs, a_fbs) = acoustic_stab;
            
            // Extraction du bloc coupling
            SparseMatrix<T> coupling_ela  = m_Cg.block(e_face_LHS_offset, a_face_LHS_offset, e_fbs, a_fbs);
            SparseMatrix<T> coupling_acou = m_Cg.block(a_face_LHS_offset, e_face_LHS_offset, a_fbs, e_fbs);
            dense_SC_ff.block(0, e_fbs, e_fbs, a_fbs) = coupling_ela;
            dense_SC_ff.block(e_fbs, 0, a_fbs, e_fbs) = coupling_acou;
            
            // Inversion
            SparseMatrix<T> SC_ff_loc = dense_SC_ff.sparseView();
            SparseLU<SparseMatrix<T>> analysis_ff;
            analysis_ff.analyzePattern(SC_ff_loc);
            analysis_ff.factorize(SC_ff_loc);
            Matrix<T, Dynamic, Dynamic> SC_ff_inv_loc = analysis_ff.solve(Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs));
            
            size_t l = 0;
            for (size_t i = 0; i < e_fbs; i++) {
                for (size_t j = 0; j < e_fbs; j++) {
                    triplets_ff[e_stride_l+l] = Triplet<T>(e_face_LHS_offset+i, e_face_LHS_offset+j, SC_ff_inv_loc(i,j));
                    l++;
                }
            } 
            l = 0;
            for (size_t i = 0; i < a_fbs; i++) {
                for (size_t j = 0; j < a_fbs; j++) {
                    triplets_ff[a_stride_l+l] = Triplet<T>(a_face_LHS_offset+i, a_face_LHS_offset+j, SC_ff_inv_loc(e_fbs+i,e_fbs+j));
                    l++;
                }
            }
            l = 0;
            // Upper right bloc
            for (size_t i = 0; i < e_fbs; i++) {
                for (size_t j = 0; j < a_fbs; j++) {
                    triplets_ff[i_stride_l+l] = Triplet<T>(e_face_LHS_offset+i, a_face_LHS_offset+j, SC_ff_inv_loc(i,e_fbs+j));
                    l++;
                }
            }
            // Lower left bloc
            for (size_t i = 0; i < a_fbs; i++) {
                for (size_t j = 0; j < e_fbs; j++) {
                    triplets_ff[i_stride_l+l] = Triplet<T>(a_face_LHS_offset+i, e_face_LHS_offset+j, SC_ff_inv_loc(e_fbs+i,j));
                    l++;
                }
            }
            cpt++;
        }        
        
        m_Sff_inv.setFromTriplets(triplets_ff.begin(), triplets_ff.end());
        triplets_ff.clear();
        return;
        
    }
    
    void inverse_Sff() {
        
        // Bi CG
        BiCGSTAB<SparseMatrix<double>> solverBiCG;
        solverBiCG.compute(m_Sff);
        if (solverBiCG.info() != Success) {
            std::cout << "Error: Matrix decomposition failed, the matrix may not be invertible";
        }
        SparseMatrix<double> identity(m_Sff.rows(), m_Sff.cols());
        identity.setIdentity();
        m_inv_Sff = solverBiCG.solve(identity);
        if (solverBiCG.info() != Success) {
            std::cout << "Error: Solving the system failed, the matrix may not be invertible";
        }
        
        return;
        
    }
    
    void setIterativeSolver(T tolerance = 1.0e-11){
        m_iterative_solver_Q = true;
        m_analysis_cg.setTolerance(tolerance);
    }
    
    void DecomposeFaceTerm(){
        
        if (m_iterative_solver_Q) {
            m_analysis_cg.compute(m_Sff);
            m_analysis_cg.setMaxIterations(m_Sff.rows());
        }
        
        else {
            m_analysis_f.analyzePattern(m_Sff);
            m_analysis_f.factorize(m_Sff);
        }
        m_sff_is_block_diagonal_Q = false;
    }
    
    void refresh_faces_unknowns(Matrix<T, Dynamic, 1> & x) {
        
        Matrix<T, Dynamic, 1> x_c_dof = x.block(0, 0, m_n_c_dof, 1);
        
        // Faces update from cells data
        Matrix<T, Dynamic, 1> RHSf = Kfc()*x_c_dof;
        if (m_sff_is_block_diagonal_Q) {
            x.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf;
        }
        else { 
            inverse_Sff();
            x.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_inv_Sff * RHSf;
        }
    }
       
//     void erk_weight_LTS_coarse(const Matrix<T, Dynamic, 1> &y, const Eigen::SparseMatrix<double> &Pcoarse, std::vector<Matrix<T, Dynamic, 1>> &w, const Matrix<T, Dynamic, 1> &Fn, const Matrix<T, Dynamic, 1> &Fn12, const Matrix<T, Dynamic, 1> &Fn1, const T dt) {
        
//         // Precompute active indices of the boolean diagonal projector Pcoarse
//         // Replaces all Pcoarse * v matvecs with O(nnz) index masking
//         m_coarse_active_dofs.clear();
//         for (int i = 0; i < Pcoarse.rows(); ++i)
//         if (Pcoarse.coeff(i, i) > 0.5)
//         m_coarse_active_dofs.push_back(i);
        
//         // Quadratic Lagrange interpolation coefficients
//         Matrix<T, Dynamic, 1> F0 =  Fn;
//         Matrix<T, Dynamic, 1> F1 = (-3*Fn + 4*Fn12 - Fn1) / dt;
//         Matrix<T, Dynamic, 1> F2 = ( 4*Fn - 8*Fn12 + 4*Fn1) / (dt*dt);
        
//         // Pcoarse is a boolean diagonal projector with static sparsity —
//         // application reduces to index masking, no matvec needed
//         auto apply_Pcoarse = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
//             Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(v.rows());
//             for (int i : m_coarse_active_dofs)
//             out(i) = v(i);
//             return out;
//         };
        
//         // Coarse projections — now O(nnz) index copies instead of O(n) sparse matvec
//         Matrix<T, Dynamic, 1> IPF0 = apply_Pcoarse(F0);
//         Matrix<T, Dynamic, 1> IPF1 = apply_Pcoarse(F1);
//         Matrix<T, Dynamic, 1> IPF2 = apply_Pcoarse(F2);
        
//         // Apply weight operator with y=0 : reduces to Mc_inv * Fc then Sff_inv * Kfc
//         // Avoids Kcc*0, Kcf*0, SetFg, ZeroFc entirely
//         auto weight_zero_y = [&](const Matrix<T, Dynamic, 1>& Fc_in, Matrix<T, Dynamic, 1>& out) {
//             out.resize(y.rows());
//             Matrix<T, Dynamic, 1> kc = m_Mc_inv * Fc_in.head(m_n_c_dof);
//             out.head(m_n_c_dof) = kc;
//             Matrix<T, Dynamic, 1> RHSf = Kfc() * kc;
//             out.tail(m_n_f_dof) = m_sff_is_block_diagonal_Q
//             ? -m_Sff_inv * RHSf
//             : -m_inv_Sff * RHSf;
//         };
        
//         Matrix<T, Dynamic, 1> F0_full,   F1_full,   F2_full;
//         Matrix<T, Dynamic, 1> IPF0_full, IPF1_full, IPF2_full;
        
//         weight_zero_y(F0,   F0_full);
//         weight_zero_y(F1,   F1_full);
//         weight_zero_y(F2,   F2_full);
//         weight_zero_y(IPF0, IPF0_full);
//         weight_zero_y(IPF1, IPF1_full);
//         weight_zero_y(IPF2, IPF2_full);
        
//         // Coarse block of Pcoarse restricted to coarse dofs
//         // Pc * v reduces to the same index masking, restricted to [0, m_n_c_dof)
//         auto apply_Pc = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
//             Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(m_n_c_dof);
//             for (int i : m_coarse_active_dofs)
//             if (i < m_n_c_dof)
//             out(i) = v(i);
//             return out;
//         };
        
//         Matrix<T, Dynamic, 1> MinvF0 = IPF0_full.head(m_n_c_dof);
//         Matrix<T, Dynamic, 1> MinvF1 = IPF1_full.head(m_n_c_dof);
//         Matrix<T, Dynamic, 1> MinvF2 = IPF2_full.head(m_n_c_dof);
        
//         // B^i y chains (homogeneous, Fc = 0)
//         Matrix<T, Dynamic, 1> B0y = y;
//         Matrix<T, Dynamic, 1> B1y, B2y, B3y;
//         erk_weight(B0y, B1y);
//         erk_weight(B1y, B2y);
//         erk_weight(B2y, B3y);
        
//         // B^i F chains — full-space F as argument
//         Matrix<T, Dynamic, 1> BF0, B2F0, BF1;
//         erk_weight(F0_full, BF0);
//         erk_weight(BF0,     B2F0);
//         erk_weight(F1_full, BF1);
        
//         // Taylor expansion arguments for each w coefficient
//         Matrix<T, Dynamic, 1> arg0 = B0y;
//         Matrix<T, Dynamic, 1> arg1 = B1y  + F0_full;
//         Matrix<T, Dynamic, 1> arg2 = B2y  + BF0    + F1_full;
//         Matrix<T, Dynamic, 1> arg3 = B3y  + B2F0   + BF1 + F2_full;
        
//         // Compute one w coefficient from its Taylor argument and optional coarse force term
//         auto compute_one_w = [&](const Matrix<T, Dynamic, 1> &arg, const Matrix<T, Dynamic, 1> *MinvFext, Matrix<T, Dynamic, 1> &wi) {
            
//             // Pcoarse * arg — index masking instead of sparse matvec
//             Matrix<T, Dynamic, 1> Ptmp   = apply_Pcoarse(arg);
//             Matrix<T, Dynamic, 1> Ptmp_c = Ptmp.head(m_n_c_dof);
//             Matrix<T, Dynamic, 1> Ptmp_f = Ptmp.tail(m_n_f_dof);
            
//             Matrix<T, Dynamic, 1> wi_c = m_Mc_inv * (-Kcc()*Ptmp_c - Kcf()*Ptmp_f);
            
//             if (MinvFext)
//             wi_c += apply_Pc(*MinvFext);
            
//             wi = Ptmp;
//             wi.head(m_n_c_dof) = wi_c;
            
//             Matrix<T, Dynamic, 1> RHSf = Kfc() * wi_c;
//             wi.tail(m_n_f_dof) = m_sff_is_block_diagonal_Q
//             ? -m_Sff_inv * RHSf
//             : -m_inv_Sff * RHSf;
//         };
        
//         compute_one_w(arg0, &MinvF0, w[0]);
//         compute_one_w(arg1, &MinvF1, w[1]);
//         compute_one_w(arg2, &MinvF2, w[2]);
//         compute_one_w(arg3, nullptr,  w[3]);
//     }
    
    void erk_weight_LTS_coarse_old(const Matrix<T, Dynamic, 1> &y, const Eigen::SparseMatrix<double> &Pcoarse, std::vector<Matrix<T, Dynamic, 1>> &w) {
        
        Matrix<T, Dynamic, 1> k = y;                 
        Matrix<T, Dynamic, 1> Biy = y; // i= 0   
        for(int i=0; i<4; i++) {
            // COMPUTATION OF B^iy
            if (i != 0) {
                erk_weight(Biy, k);
                Biy = k;
            }
            // COMPUTATION OF w_i
            k = Pcoarse * Biy;
            erk_weight(k, w[i]);
        }
    }
    
    void erk_weight(Matrix<T, Dynamic, 1> & y, Matrix<T, Dynamic, 1> & k) {
        
        k=y;
        Matrix<T, Dynamic, 1> y_c_dof = y.block(0, 0, m_n_c_dof, 1);
        Matrix<T, Dynamic, 1> y_f_dof = y.block(m_n_c_dof, 0, m_n_f_dof, 1);
        
        ////////// CELLS UPDATE
        Matrix<T, Dynamic, 1> RHSc = Fc() - Kcc()*y_c_dof - Kcf()*y_f_dof;
        Matrix<T, Dynamic, 1> k_c_dof = m_Mc_inv * RHSc;
        k.block(0, 0, m_n_c_dof, 1) = k_c_dof;
        
        // FACES UPDATE 
        Matrix<T, Dynamic, 1> RHSf = Kfc()*k_c_dof ;
        if (m_sff_is_block_diagonal_Q) {
            k.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_Sff_inv * RHSf; 
        }
        else {
            k.block(m_n_c_dof, 0, m_n_f_dof, 1) = - m_inv_Sff * RHSf; 
        }
        
    }
    
void erk_weight_LTS_coarse(const Matrix<T, Dynamic, 1> &y,
                            const Eigen::SparseMatrix<double> &Pcoarse,
                            std::vector<Matrix<T, Dynamic, 1>> &w,
                            const Matrix<T, Dynamic, 1> &Fn,
                            const Matrix<T, Dynamic, 1> &Fn12,
                            const Matrix<T, Dynamic, 1> &Fn1,
                            const T dt) {

    // Pcoarse = (I-P) in the algorithm — projects onto the coarse subspace
    // B = erk_weight

    // Build active index list once from the boolean diagonal projector
    if (m_coarse_active_dofs.empty()) {
        for (int i = 0; i < Pcoarse.rows(); ++i) {
            if (Pcoarse.coeff(i, i) > 0.5) {
                m_coarse_active_dofs.push_back(i);
            }
        }
    }

    // (I-P) * v — index masking, no sparse matvec
    auto IP = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(v.rows());
        for (int i : m_coarse_active_dofs) {
            out(i) = v(i);
        }
        return out;
    };

    // (I-P) restricted to cell dofs — for the Pc * MinvF term
    auto IP_c = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(m_n_c_dof);
        for (int i : m_coarse_active_dofs) {
            if (i < m_n_c_dof) {
                out(i) = v(i);
            }
        }
        return out;
    };

    // B applied with y=0 — reduces to Mc_inv * Fc then Sff_inv * Kfc * kc
    // Used to compute B * (I-P) * Fi with zero displacement argument
    auto B_zero_y = [&](const Matrix<T, Dynamic, 1>& Fc_in,
                         Matrix<T, Dynamic, 1>& out) {
        out.resize(y.rows());
        Matrix<T, Dynamic, 1> kc = m_Mc_inv * Fc_in.head(m_n_c_dof);
        out.head(m_n_c_dof) = kc;
        Matrix<T, Dynamic, 1> RHSf = Kfc() * kc;
        out.tail(m_n_f_dof) = m_sff_is_block_diagonal_Q
                               ? -m_Sff_inv * RHSf
                               : -m_inv_Sff * RHSf;
    };

    // Quadratic Lagrange interpolation of F over [tn, tn+1]
    // F0 = Fn, F1 = dF/dt, F2 = d²F/dt²  (coefficients, not values)
    Matrix<T, Dynamic, 1> F0 =  Fn;
    Matrix<T, Dynamic, 1> F1 = (-3*Fn + 4*Fn12 - Fn1) / dt;
    Matrix<T, Dynamic, 1> F2 = ( 4*Fn - 8*Fn12 + 4*Fn1) / (dt*dt);
    Matrix<T, Dynamic, 1> IPF0 = IP(F0);
    Matrix<T, Dynamic, 1> IPF1 = IP(F1);
    Matrix<T, Dynamic, 1> IPF2 = IP(F2);
    Matrix<T, Dynamic, 1> BIPFn_0, BIPFn_1, BIPFn_2;
    Matrix<T, Dynamic, 1> BFn_0,   BFn_1,   BFn_2;
    B_zero_y(IPF0, BIPFn_0);   // B(I-P)F0
    B_zero_y(IPF1, BIPFn_1);   // B(I-P)F1
    B_zero_y(IPF2, BIPFn_2);   // B(I-P)F2
    B_zero_y(F0,   BFn_0);     // B*F0  (full, for B^i F chains)
    B_zero_y(F1,   BFn_1);     // B*F1
    B_zero_y(F2,   BFn_2);     // B*F2

    // MinvF = (I-P) * B * (I-P) * Fi restricted to cell dofs
    // = IP_c applied to BIPFn_i  (used in w_{n,i} = B(I-P)(arg) + (I-P)(Fi) term)
    Matrix<T, Dynamic, 1> MinvF0 = BIPFn_0.head(m_n_c_dof);
    Matrix<T, Dynamic, 1> MinvF1 = BIPFn_1.head(m_n_c_dof);
    Matrix<T, Dynamic, 1> MinvF2 = BIPFn_2.head(m_n_c_dof);

    // B^i * yn chains — Algorithm step 2: B^0 yn = yn, B^1 yn, B^2 yn, B^3 yn
    Matrix<T, Dynamic, 1> B0yn = y;
    Matrix<T, Dynamic, 1> B1yn, B2yn, B3yn;
    erk_weight(B0yn, B1yn);
    erk_weight(B1yn, B2yn);
    erk_weight(B2yn, B3yn);

    // B^i * F0 chains — for the BF and B²F terms in w_{n,2} and w_{n,3}
    Matrix<T, Dynamic, 1> BF0, B2F0, BF1;
    erk_weight(BFn_0, BF0);
    erk_weight(BF0,   B2F0);
    erk_weight(BFn_1, BF1);

    // arg for w_{n,0} = B^0 yn
    // arg for w_{n,1} = B^1 yn + F0
    // arg for w_{n,2} = B^2 yn + B*F0 + F1
    // arg for w_{n,3} = B^3 yn + B²*F0 + B*F1 + F2
    Matrix<T, Dynamic, 1> arg0 = B0yn;
    Matrix<T, Dynamic, 1> arg1 = B1yn + BFn_0;
    Matrix<T, Dynamic, 1> arg2 = B2yn + BF0   + BFn_1;
    Matrix<T, Dynamic, 1> arg3 = B3yn + B2F0  + BF1   + BFn_2;

    // Compute w_{n,i} = B(I-P)(arg_i) + (I-P)(F_i)
    // = erk_weight applied to (I-P)*arg, plus coarse force correction
    auto compute_w = [&](const Matrix<T, Dynamic, 1>& arg,
                          const Matrix<T, Dynamic, 1>* MinvFext,
                          Matrix<T, Dynamic, 1>& wi) {

        // (I-P) * arg
        Matrix<T, Dynamic, 1> IParg   = IP(arg);
        Matrix<T, Dynamic, 1> IParg_c = IParg.head(m_n_c_dof);
        Matrix<T, Dynamic, 1> IParg_f = IParg.tail(m_n_f_dof);

        // B * (I-P) * arg
        Matrix<T, Dynamic, 1> wi_c = m_Mc_inv * (-Kcc()*IParg_c - Kcf()*IParg_f);

        // Add (I-P) * Fi term (coarse force correction)
        if (MinvFext) {
            wi_c += IP_c(*MinvFext);
        }

        // Assemble wi — starts from IParg (coarse mask), cell part overwritten by solve
        wi = IParg;
        wi.head(m_n_c_dof) = wi_c;
        Matrix<T, Dynamic, 1> RHSf = Kfc() * wi_c;
        wi.tail(m_n_f_dof) = m_sff_is_block_diagonal_Q
                              ? -m_Sff_inv * RHSf
                              : -m_inv_Sff * RHSf;
    };

    // Algorithm 3 step 2 — compute w_{n,0} ... w_{n,3}
    compute_w(arg0, &MinvF0, w[0]);
    compute_w(arg1, &MinvF1, w[1]);
    compute_w(arg2, &MinvF2, w[2]);
    compute_w(arg3, nullptr,  w[3]);
}


void erk_weight_LTS_fine(Matrix<T, Dynamic, 1> &x_dof_n,
                          const Eigen::SparseMatrix<T> &Pfine,
                          const std::vector<Matrix<T, Dynamic, 1>> &w,
                          const Matrix<T, Dynamic, 1> &Fm,
                          const Matrix<T, Dynamic, 1> &Fmh,
                          const Matrix<T, Dynamic, 1> &Fm1,
                          const T tm,
                          const T dtau) {

    // Pfine = P in the algorithm — projects onto the fine subspace
    // B = erk_weight
    // Taylor_w(tau) = w_{n,0} + tau*w_{n,1} + (tau²/2)*w_{n,2} + (tau³/6)*w_{n,3}
    //               = coarse contribution evaluated at local time tau = m*dtau

    // Taylor expansion of w at local time tau — Algorithm 3 step 3
    auto Taylor_w = [&](T tau) -> Matrix<T, Dynamic, 1> {
        T tau2 = tau*tau, tau3 = tau*tau2;
        return w[0] + tau*w[1] + (tau2/2.0)*w[2] + (tau3/6.0)*w[3];
    };

    T tmh = tm + 0.5*dtau;
    T tm1 = tm + dtau;

    // One RK4 stage — matches Algorithm 3 step 3:
    // ki = Taylor_w(tau) + B*P*ỹ_stage + P*F_{n,tau}
    //    = coarse Taylor term + fine erk_weight(P*y) with P*F as force
    auto fine_stage = [&](const Matrix<T, Dynamic, 1>& y_stage,
                           T tau,
                           const Matrix<T, Dynamic, 1>& F_tau) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> Py = Pfine * y_stage;   // P * ỹ_stage
        Matrix<T, Dynamic, 1> PF = Pfine * F_tau;     // P * F_{n,tau}
        SetFg(PF);                                     // set fine force
        Matrix<T, Dynamic, 1> k;
        erk_weight(Py, k);                             // B * P * ỹ_stage
        ZeroFc();
        k += Taylor_w(tau);                            // + coarse Taylor term
        return k;
        
    };

    // Algorithm 3 step 3 — RK4 loop for m-th fine sub-step
    Matrix<T, Dynamic, 1> k1 = fine_stage(x_dof_n,                tm,  Fm);
    Matrix<T, Dynamic, 1> k2 = fine_stage(x_dof_n + 0.5*dtau*k1, tmh, Fmh);
    Matrix<T, Dynamic, 1> k3 = fine_stage(x_dof_n + 0.5*dtau*k2, tmh, Fmh);
    Matrix<T, Dynamic, 1> k4 = fine_stage(x_dof_n +     dtau*k3, tm1, Fm1);

    // Algorithm 3 step 3 — update ỹ_{m+1/p}
    x_dof_n += dtau * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
}





    #ifdef HAVE_INTEL_MKL
    PardisoLDLT<SparseMatrix<T>> & FacesAnalysis(){
        return m_analysis_f;
    }
    #else
    SimplicialLDLT<SparseMatrix<T>> & FacesAnalysis(){
        return m_analysis_f;
    }
    #endif
    
    SparseMatrix<T> & Mc(){
        return m_Mc;
    }
    
    SparseMatrix<T> & invMc(){
        return m_Mc_inv;
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
    
    SparseMatrix<T> & Sff(){
        return m_Sff;
    }
    
    Matrix<T, Dynamic, 1> & Fc(){
        return m_Fc;
    }
    
    void SetFg(Matrix<T, Dynamic, 1> & Fg){
        m_Fc = Fg.block(0, 0, m_n_c_dof, 1);
    }
    
    void ZeroFc() {
        m_Fc.setZero();
    }
    
    SparseMatrix<T> & SffInv(){
        return m_Sff_inv;
    }

    size_t n_c_dof() const { return m_n_c_dof; }
    size_t n_f_dof() const { return m_n_f_dof; }
    
    void compute_eigenvalues(std::ostream & simulation_log = std::cout){
        
        SparseMatrix<T> A_SCHUR = m_Kcc - m_Kcf*m_Sff_inv*m_Kfc;
        SparseMatrix<T> A = A_SCHUR.transpose()*m_Mc_inv*A_SCHUR;
        Spectra::SparseSymMatProd<double> opA(A);
        Spectra::SparseCholesky<double>   opB(m_Mc);
        Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> eigs(opA, opB, 1, 10);
        // Initialize and compute
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        std::cout << std::endl << bold << red << "   Computation of the eigenvalues: " << reset;
        std::cout << std::endl << bold << cyan << "      State of the computation: " << reset;
        bool debug = true;
        if (debug) {
            if(eigs.info() == Spectra::CompInfo::Successful)
            std::cout << cyan << bold << "Successful\n";
            if(eigs.info() == Spectra::CompInfo::NotComputed)
            std::cout << cyan << bold << "NotComputed\n";
            if(eigs.info() == Spectra::CompInfo::NotConverging)
            std::cout << cyan << bold << "NotConverging\n";
            if(eigs.info() == Spectra::CompInfo::NumericalIssue)
            std::cout << cyan << bold << "NumericalIssue\n";
        }
        eigs.eigenvalues();
        std::cout << bold << cyan << "      Eigenvalue found: " << reset << cyan << eigs.eigenvalues() << std::endl << std::endl; 
        simulation_log << "Eigenvalue found: " << eigs.eigenvalues() << std::endl;
        
    }
    
    void compute_eigenvalues_bis(SparseMatrix<T> LHS_STAB, std::pair<size_t,size_t> block_dimension, std::ostream & simulation_log = std::cout){
        
        auto ten_bs = block_dimension.first;         // Elastic block
        auto vec_cell_size = block_dimension.second; // Acoustic block
        SparseMatrix<T> m_Scc = LHS_STAB.block(0,0, m_n_c_dof, m_n_c_dof); // Stabilisation
        
        SparseMatrix<T> Delta = m_Kcf*m_Sff_inv*m_Kfc;
        SparseMatrix<T> A_SCHUR = m_Kcc - Delta;
        SparseMatrix<T> S_SCHUR = m_Scc - Delta;
        SparseMatrix<T> A = A_SCHUR.transpose()*m_Mc_inv*A_SCHUR - S_SCHUR.transpose()*m_Mc_inv*S_SCHUR;
        Spectra::SparseSymMatProd<double> opA(A);
        Spectra::SparseCholesky<double>   opB(m_Mc);
        Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> eigs(opA, opB, 1, 10);
        // Initialize and compute
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        bool debug = true;
        if (debug) {
            if(eigs.info() == Spectra::CompInfo::Successful)
            std::cout << "Successful\n";
            if(eigs.info() == Spectra::CompInfo::NotComputed)
            std::cout << "NotComputed\n";
            if(eigs.info() == Spectra::CompInfo::NotConverging)
            std::cout << "NotConverging\n";
            if(eigs.info() == Spectra::CompInfo::NumericalIssue)
            std::cout << "NumericalIssue\n";
        }
        eigs.eigenvalues();
        std::cout << std::endl;
        std::cout << bold << red << "   Eigenvalue found: " << reset << eigs.eigenvalues();
        simulation_log << "Eigenvalue found: " << eigs.eigenvalues() << std::endl;
        
    }
    

void build_LTS_subspaces(const Eigen::SparseMatrix<T>& Pcoarse,
                          const Eigen::SparseMatrix<T>& Pfine) {

    // -------------------------------------------------------
    // STEP 1 : collecter les indices actifs via InnerIterator
    // -------------------------------------------------------
    m_coarse_active_dofs.clear(); m_coarse_active_c.clear(); m_coarse_active_f.clear();
    m_fine_active_dofs.clear();   m_fine_active_c.clear();   m_fine_active_f.clear();

    for (int k = 0; k < Pcoarse.outerSize(); ++k)
        for (typename SparseMatrix<T>::InnerIterator it(Pcoarse, k); it; ++it) {
            int i = (int)it.row();
            m_coarse_active_dofs.push_back(i);
            if (i < (int)m_n_c_dof) m_coarse_active_c.push_back(i);
            else                     m_coarse_active_f.push_back(i - (int)m_n_c_dof);
        }

    for (int k = 0; k < Pfine.outerSize(); ++k)
        for (typename SparseMatrix<T>::InnerIterator it(Pfine, k); it; ++it) {
            int i = (int)it.row();
            m_fine_active_dofs.push_back(i);
            if (i < (int)m_n_c_dof) m_fine_active_c.push_back(i);
            else                     m_fine_active_f.push_back(i - (int)m_n_c_dof);
        }

    // std::cout << "LTS subspaces — indices collected:" << std::endl;
    // std::cout << "  m_n_c_dof = " << m_n_c_dof << std::endl;
    // std::cout << "  m_n_f_dof = " << m_n_f_dof << std::endl;
    // std::cout << "  coarse: " << m_coarse_active_c.size() << " cell dofs, "
    //                            << m_coarse_active_f.size() << " face dofs" << std::endl;
    // std::cout << "  fine:   " << m_fine_active_c.size()   << " cell dofs, "
    //                            << m_fine_active_f.size()   << " face dofs" << std::endl;

    // Verification
for (int i : m_fine_active_c)   assert(i >= 0 && i < (int)m_n_c_dof);
for (int i : m_fine_active_f)   assert(i >= 0 && i < (int)m_n_f_dof);
for (int i : m_coarse_active_c) assert(i >= 0 && i < (int)m_n_c_dof);
for (int i : m_coarse_active_f) assert(i >= 0 && i < (int)m_n_f_dof);

    // -------------------------------------------------------
    // STEP 2 : lambda d'extraction — itère sur tous les nnz
    // correct indépendamment du storage order ColMajor/RowMajor
    // -------------------------------------------------------
    auto extract = [](const SparseMatrix<T>& M,
                      const std::vector<int>& rows,
                      const std::vector<int>& cols) -> SparseMatrix<T> {

        // tables de lookup global -> local
        std::unordered_map<int,int> row_map;
        row_map.reserve(rows.size());
        for (int i = 0; i < (int)rows.size(); ++i)
            row_map[rows[i]] = i;

        std::unordered_map<int,int> col_map;
        col_map.reserve(cols.size());
        for (int j = 0; j < (int)cols.size(); ++j)
            col_map[cols[j]] = j;

        std::vector<Triplet<T>> trips;
        trips.reserve(rows.size() * 8);

        // Itérer sur tous les non-zeros — it.row() et it.col() toujours corrects
        for (int k = 0; k < M.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(M, k); it; ++it) {
                auto row_it = row_map.find((int)it.row());
                if (row_it == row_map.end()) continue;
                auto col_it = col_map.find((int)it.col());
                if (col_it == col_map.end()) continue;
                trips.emplace_back(row_it->second, col_it->second, it.value());
            }
        }

        SparseMatrix<T> out((int)rows.size(), (int)cols.size());
        out.setFromTriplets(trips.begin(), trips.end());
        return out;
    };

    // -------------------------------------------------------
    // STEP 3 : sous-blocs fine
    // -------------------------------------------------------
    // std::cout << "  extracting fine blocks..." << std::endl;
    m_Kcc_fine     = extract(m_Kcc,     m_fine_active_c, m_fine_active_c);
    m_Kcf_fine     = extract(m_Kcf,     m_fine_active_c, m_fine_active_f);
    m_Kfc_fine     = extract(m_Kfc,     m_fine_active_f, m_fine_active_c);
    m_Mc_inv_fine  = extract(m_Mc_inv,  m_fine_active_c, m_fine_active_c);
    m_Sff_inv_fine = extract(m_Sff_inv, m_fine_active_f, m_fine_active_f);

    // std::cout << "    Kcc_fine:     " << m_Kcc_fine.rows()     << " x " << m_Kcc_fine.cols()     << std::endl;
    // std::cout << "    Kcf_fine:     " << m_Kcf_fine.rows()     << " x " << m_Kcf_fine.cols()     << std::endl;
    // std::cout << "    Kfc_fine:     " << m_Kfc_fine.rows()     << " x " << m_Kfc_fine.cols()     << std::endl;
    // std::cout << "    Mc_inv_fine:  " << m_Mc_inv_fine.rows()  << " x " << m_Mc_inv_fine.cols()  << std::endl;
    // std::cout << "    Sff_inv_fine: " << m_Sff_inv_fine.rows() << " x " << m_Sff_inv_fine.cols() << std::endl;

    // Vérifier que les sous-blocs ne sont pas vides si la zone est non vide
    if (!m_fine_active_c.empty()) {
        // std::cout << "    Kcc_fine nnz:     " << m_Kcc_fine.nonZeros()     << std::endl;
        // std::cout << "    Mc_inv_fine nnz:  " << m_Mc_inv_fine.nonZeros()  << std::endl;
        // std::cout << "    Sff_inv_fine nnz: " << m_Sff_inv_fine.nonZeros() << std::endl;
    }

    // -------------------------------------------------------
    // STEP 4 : sous-blocs coarse
    // -------------------------------------------------------
    // std::cout << "  extracting coarse blocks..." << std::endl;
    m_Kcc_coarse     = extract(m_Kcc,     m_coarse_active_c, m_coarse_active_c);
    m_Kcf_coarse     = extract(m_Kcf,     m_coarse_active_c, m_coarse_active_f);
    m_Kfc_coarse     = extract(m_Kfc,     m_coarse_active_f, m_coarse_active_c);
    m_Mc_inv_coarse  = extract(m_Mc_inv,  m_coarse_active_c, m_coarse_active_c);
    m_Sff_inv_coarse = extract(m_Sff_inv, m_coarse_active_f, m_coarse_active_f);

    // std::cout << "    Kcc_coarse:     " << m_Kcc_coarse.rows()     << " x " << m_Kcc_coarse.cols()     << std::endl;
    // std::cout << "    Kcf_coarse:     " << m_Kcf_coarse.rows()     << " x " << m_Kcf_coarse.cols()     << std::endl;
    // std::cout << "    Kfc_coarse:     " << m_Kfc_coarse.rows()     << " x " << m_Kfc_coarse.cols()     << std::endl;
    // std::cout << "    Mc_inv_coarse:  " << m_Mc_inv_coarse.rows()  << " x " << m_Mc_inv_coarse.cols()  << std::endl;
    // std::cout << "    Sff_inv_coarse: " << m_Sff_inv_coarse.rows() << " x " << m_Sff_inv_coarse.cols() << std::endl;

    // if (!m_coarse_active_c.empty()) {
    //     std::cout << "    Kcc_coarse nnz:     " << m_Kcc_coarse.nonZeros()     << std::endl;
    //     std::cout << "    Mc_inv_coarse nnz:  " << m_Mc_inv_coarse.nonZeros()  << std::endl;
    //     std::cout << "    Sff_inv_coarse nnz: " << m_Sff_inv_coarse.nonZeros() << std::endl;
    // }

    // std::cout << "LTS subspaces built successfully." << std::endl;
}
void erk_weight_restricted(const Matrix<T,Dynamic,1>& y,
                            Matrix<T,Dynamic,1>& k,
                            const std::vector<int>& active_c,
                            const std::vector<int>& active_f,
                            const SparseMatrix<T>& Kcc_loc,
                            const SparseMatrix<T>& Kcf_loc,
                            const SparseMatrix<T>& Kfc_loc,
                            const SparseMatrix<T>& Mc_inv_loc,
                            const SparseMatrix<T>& Sff_inv_loc) {

    const int nc = (int)active_c.size();
    const int nf = (int)active_f.size();

    Matrix<T,Dynamic,1> yc(nc), Fc_loc(nc);
    for (int i = 0; i < nc; ++i) {
        yc(i)     = y(active_c[i]);
        Fc_loc(i) = m_Fc(active_c[i]);
    }

    Matrix<T,Dynamic,1> yf(nf);
    for (int i = 0; i < nf; ++i)
        yf(i) = y(m_n_c_dof + active_f[i]);

    Matrix<T,Dynamic,1> kc = Mc_inv_loc * (Fc_loc - Kcc_loc*yc - Kcf_loc*yf);
    Matrix<T,Dynamic,1> kf = -Sff_inv_loc * (Kfc_loc * kc);

    k = Matrix<T,Dynamic,1>::Zero(y.rows());
    for (int i = 0; i < nc; ++i) k(active_c[i])             = kc(i);
    for (int i = 0; i < nf; ++i) k(m_n_c_dof + active_f[i]) = kf(i);

}


void erk_weight_LTS_fine_restricted(Matrix<T, Dynamic, 1> &x_dof_n,
                                     const Eigen::SparseMatrix<T> &Pfine,
                                     const std::vector<Matrix<T, Dynamic, 1>> &w,
                                     const Matrix<T, Dynamic, 1> &Fm,
                                     const Matrix<T, Dynamic, 1> &Fmh,
                                     const Matrix<T, Dynamic, 1> &Fm1,
                                     const T tm,
                                     const T dtau) {

    auto Taylor_w = [&](T tau) -> Matrix<T, Dynamic, 1> {
        T tau2 = tau*tau, tau3 = tau*tau2;
        return w[0] + tau*w[1] + (tau2/2.0)*w[2] + (tau3/6.0)*w[3];
    };

    T tmh = tm + 0.5*dtau;
    T tm1 = tm +     dtau;

    auto fine_stage = [&](const Matrix<T, Dynamic, 1>& y_stage,
                           T tau,
                           const Matrix<T, Dynamic, 1>& F_tau) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> Py = Pfine * y_stage;
        Matrix<T, Dynamic, 1> PF = Pfine * F_tau;
        SetFg(PF);
        Matrix<T, Dynamic, 1> k;
        erk_weight_restricted(Py, k,
            m_fine_active_c, m_fine_active_f,
            m_Kcc_fine, m_Kcf_fine, m_Kfc_fine,
            m_Mc_inv_fine, m_Sff_inv_fine);
        ZeroFc();
        k += Taylor_w(tau);
        return k;
    };

    Matrix<T, Dynamic, 1> k1 = fine_stage(x_dof_n,                tm,  Fm);
    Matrix<T, Dynamic, 1> k2 = fine_stage(x_dof_n + 0.5*dtau*k1, tmh, Fmh);
    Matrix<T, Dynamic, 1> k3 = fine_stage(x_dof_n + 0.5*dtau*k2, tmh, Fmh);
    Matrix<T, Dynamic, 1> k4 = fine_stage(x_dof_n +     dtau*k3, tm1, Fm1);

    // *** CORRECTION : mettre à jour seulement les dofs fins ***
    // Les dofs coarses sont gérés par erk_weight_LTS_coarse_restricted
    Matrix<T, Dynamic, 1> dx = dtau * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    for (int i : m_fine_active_c)             x_dof_n(i)             += dx(i);
    for (int i : m_fine_active_f) x_dof_n(m_n_c_dof + i) += dx(m_n_c_dof + i);
}


void erk_weight_LTS_coarse_restricted(const Matrix<T, Dynamic, 1> &y,
                                       const Eigen::SparseMatrix<double> &Pcoarse,
                                       std::vector<Matrix<T, Dynamic, 1>> &w,
                                       const Matrix<T, Dynamic, 1> &Fn,
                                       const Matrix<T, Dynamic, 1> &Fn12,
                                       const Matrix<T, Dynamic, 1> &Fn1,
                                       const T dt) {

    auto IP = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(v.rows());
        for (int i : m_coarse_active_dofs) out(i) = v(i);
        return out;
    };

    auto IP_c = [&](const Matrix<T, Dynamic, 1>& v) -> Matrix<T, Dynamic, 1> {
        Matrix<T, Dynamic, 1> out = Matrix<T, Dynamic, 1>::Zero(m_n_c_dof);
        for (int i : m_coarse_active_dofs)
            if (i < (int)m_n_c_dof) out(i) = v(i);
        return out;
    };

    auto B_zero_y = [&](const Matrix<T, Dynamic, 1>& Fc_in, Matrix<T, Dynamic, 1>& out) {
        out.resize(y.rows());
        Matrix<T, Dynamic, 1> kc = m_Mc_inv * Fc_in.head(m_n_c_dof);
        out.head(m_n_c_dof) = kc;
        Matrix<T, Dynamic, 1> RHSf = Kfc() * kc;
        out.tail(m_n_f_dof) = m_sff_is_block_diagonal_Q ? -m_Sff_inv * RHSf : -m_inv_Sff * RHSf;
    };

    Matrix<T, Dynamic, 1> F0 =  Fn;
    Matrix<T, Dynamic, 1> F1 = (-3*Fn + 4*Fn12 - Fn1) / dt;
    Matrix<T, Dynamic, 1> F2 = ( 4*Fn - 8*Fn12 + 4*Fn1) / (dt*dt);

    Matrix<T, Dynamic, 1> IPF0 = IP(F0), IPF1 = IP(F1), IPF2 = IP(F2);

    Matrix<T, Dynamic, 1> BIPFn_0, BIPFn_1, BIPFn_2, BFn_0, BFn_1, BFn_2;
    B_zero_y(IPF0, BIPFn_0); B_zero_y(IPF1, BIPFn_1); B_zero_y(IPF2, BIPFn_2);
    B_zero_y(F0,   BFn_0);   B_zero_y(F1,   BFn_1);   B_zero_y(F2,   BFn_2);

    Matrix<T, Dynamic, 1> MinvF0 = BIPFn_0.head(m_n_c_dof);
    Matrix<T, Dynamic, 1> MinvF1 = BIPFn_1.head(m_n_c_dof);
    Matrix<T, Dynamic, 1> MinvF2 = BIPFn_2.head(m_n_c_dof);

    Matrix<T, Dynamic, 1> B0yn = y, B1yn, B2yn, B3yn;
    erk_weight(B0yn, B1yn); erk_weight(B1yn, B2yn); erk_weight(B2yn, B3yn);

    Matrix<T, Dynamic, 1> BF0, B2F0, BF1;
    erk_weight(BFn_0, BF0); erk_weight(BF0, B2F0); erk_weight(BFn_1, BF1);

    Matrix<T, Dynamic, 1> arg0 = B0yn;
    Matrix<T, Dynamic, 1> arg1 = B1yn + BFn_0;
    Matrix<T, Dynamic, 1> arg2 = B2yn + BF0  + BFn_1;
    Matrix<T, Dynamic, 1> arg3 = B3yn + B2F0 + BF1 + BFn_2;

    // *** SEUL CHANGEMENT : compute_w utilise les sous-blocs coarse ***
    const int nc = (int)m_coarse_active_c.size();
    const int nf = (int)m_coarse_active_f.size();

    auto compute_w_restricted = [&](const Matrix<T, Dynamic, 1>& arg,
                                     const Matrix<T, Dynamic, 1>* MinvFext,
                                     Matrix<T, Dynamic, 1>& wi) {

        Matrix<T, Dynamic, 1> IParg = IP(arg);

        // Gather
        Matrix<T, Dynamic, 1> IParg_c_loc(nc), IParg_f_loc(nf);
        for (int i = 0; i < nc; ++i) IParg_c_loc(i) = IParg(m_coarse_active_c[i]);
        for (int i = 0; i < nf; ++i) IParg_f_loc(i) = IParg(m_n_c_dof + m_coarse_active_f[i]);

        // Produit restreint
        Matrix<T, Dynamic, 1> wi_c_loc = m_Mc_inv_coarse *
            (-m_Kcc_coarse*IParg_c_loc - m_Kcf_coarse*IParg_f_loc);

        // Correction force coarse
        if (MinvFext) {
            Matrix<T, Dynamic, 1> tmp = IP_c(*MinvFext);
            Matrix<T, Dynamic, 1> IPF_loc(nc);
            for (int i = 0; i < nc; ++i) IPF_loc(i) = tmp(m_coarse_active_c[i]);
            wi_c_loc += IPF_loc;
        }

        // Faces restreintes
        Matrix<T, Dynamic, 1> kf_loc = -m_Sff_inv_coarse * (m_Kfc_coarse * wi_c_loc);

        // Scatter
        wi = IParg;
        for (int i = 0; i < nc; ++i) wi(m_coarse_active_c[i])             = wi_c_loc(i);
        for (int i = 0; i < nf; ++i) wi(m_n_c_dof + m_coarse_active_f[i]) = kf_loc(i);
    };

    compute_w_restricted(arg0, &MinvF0, w[0]);
    compute_w_restricted(arg1, &MinvF1, w[1]);
    compute_w_restricted(arg2, &MinvF2, w[2]);
    compute_w_restricted(arg3, nullptr,  w[3]);
}


// =====================================================================
// Multi-level LTS-RK4 building blocks.
//
// PURELY ADDITIVE: new struct + new methods only, none of the existing
// methods above are touched, so the 11 existing 2-level LTS callers
// (ERK4_LTS.hpp, ERK4_LTS_conv_test.hpp, ERK4_LTS_stab*.hpp,
// ERK4_LTS_Lshape_conv_test.hpp, ERK4_LTS_Lshape_MMS_conv_test.hpp,
// ERK4_LTS_optimised.hpp, ERK4_LTS_SSTAB.hpp, ...) and assemble_P in
// elastoacoustic_four_fields_assembler.hpp are completely unaffected.
//
// Unlike build_LTS_subspaces() (single hardcoded coarse/fine member-slot
// pair, only one alive at a time), build_LTS_subblock() returns a
// self-contained, independently-alive LTS_subblock_set per call, so a
// std::vector<LTS_subblock_set> (one per multi-level band) can coexist.
//
// Unlike erk_weight_LTS_coarse_restricted/erk_weight_LTS_fine_restricted
// (which only restrict the FINAL solve, still doing the B-chain and the
// scatter over the WHOLE mesh / only "coarse-classified" faces
// respectively -- see session notes on Bug A/B), erk_weight_LTS_coarse_v2
// restricts EVERY step (including the B^i-chain) to a band's own dofs
// plus a halo wide enough to stay exact (Kcc is exactly block-diagonal
// per cell in HHO -- cells only couple through faces -- so a halo in the
// face direction, expanded a few rings via the Kfc/Kcf sparsity graph,
// suffices), and scatters its output over the FULL active set (own +
// halo), fixing the "interface face dropped" bug. erk_weight_LTS_fine_v2
// only needs to be exact at a band's own dofs (the Taylor polynomial
// additive term already carries every other band's influence), so no
// halo is required there, and it scatters exactly at the band's own
// positions (fixing the "coarse dofs never written" bug).
// =====================================================================

struct LTS_subblock_set {
    std::vector<int> active_c, active_f;   // band's OWN dofs first, halo dofs after
    size_t n_own_c = 0, n_own_f = 0;       // how many of the above are this band's own (rest is halo context)
    SparseMatrix<T> Kcc, Kcf, Kfc, Mc_inv, Sff_inv;             // sized to active_c/active_f (own+halo) -- used by the "coarse role" (erk_weight_LTS_coarse_v2), which needs the halo for its B-chain to stay exact
    SparseMatrix<T> Kcc_own, Kcf_own, Kfc_own, Mc_inv_own, Sff_inv_own; // sized to just own_c/own_f (no halo) -- used by the "fine role" (erk_weight_LTS_fine_v2), which needs no halo (the ancestor/descendant influence comes in purely additively via the Taylor term, never through Kcf/Kfc directly)
};

private:

static SparseMatrix<T> extract_block(const SparseMatrix<T>& M,
                                      const std::vector<int>& rows,
                                      const std::vector<int>& cols) {
    std::unordered_map<int,int> row_map; row_map.reserve(rows.size());
    for (int i = 0; i < (int)rows.size(); ++i) row_map[rows[i]] = i;
    std::unordered_map<int,int> col_map; col_map.reserve(cols.size());
    for (int j = 0; j < (int)cols.size(); ++j) col_map[cols[j]] = j;

    std::vector<Triplet<T>> trips;
    trips.reserve(rows.size() * 8);
    for (int k = 0; k < M.outerSize(); ++k) {
        for (typename SparseMatrix<T>::InnerIterator it(M, k); it; ++it) {
            auto row_it = row_map.find((int)it.row());
            if (row_it == row_map.end()) continue;
            auto col_it = col_map.find((int)it.col());
            if (col_it == col_map.end()) continue;
            trips.emplace_back(row_it->second, col_it->second, it.value());
        }
    }
    SparseMatrix<T> out((int)rows.size(), (int)cols.size());
    out.setFromTriplets(trips.begin(), trips.end());
    return out;
}

public:

// Builds a self-contained restricted sub-block set for one multi-level
// band: `own_c`/`own_f` are the band's own (local, global-indexed) cell
// and face dof indices; `halo_rings` expands the set by that many
// cell<->face BFS hops through the Kfc sparsity graph (a face is added
// if it touches an already-included cell, a cell is added if it touches
// an already-included face), so that every B-application used inside
// erk_weight_LTS_coarse_v2 (up to a 3-deep chain, plus the force chains)
// stays exact within the band's own dofs. Generous by default (safety
// over squeezing out the last bit of speed) -- verified numerically by
// the L=2 exact-reproduction gate before being trusted for L>2.
LTS_subblock_set build_LTS_subblock(const std::vector<int>& own_c,
                                     const std::vector<int>& own_f,
                                     int halo_rings = 6) const {

    std::set<int> set_c(own_c.begin(), own_c.end());
    std::set<int> set_f(own_f.begin(), own_f.end());

    for (int ring = 0; ring < halo_rings; ++ring) {
        std::set<int> new_c, new_f;
        for (int k = 0; k < m_Kfc.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Kfc, k); it; ++it) {
                int face = (int)it.row(), cell = (int)it.col();
                bool cell_in = set_c.count(cell) > 0;
                bool face_in = set_f.count(face) > 0;
                if (cell_in && !face_in) new_f.insert(face);
                if (face_in && !cell_in) new_c.insert(cell);
            }
        }
        if (new_c.empty() && new_f.empty()) break;
        set_c.insert(new_c.begin(), new_c.end());
        set_f.insert(new_f.begin(), new_f.end());
    }

    std::set<int> own_c_set(own_c.begin(), own_c.end());
    std::set<int> own_f_set(own_f.begin(), own_f.end());
    std::vector<int> halo_c, halo_f;
    for (int c : set_c) if (!own_c_set.count(c)) halo_c.push_back(c);
    for (int f : set_f) if (!own_f_set.count(f)) halo_f.push_back(f);

    LTS_subblock_set out;
    out.n_own_c = own_c.size();
    out.n_own_f = own_f.size();
    out.active_c = own_c; out.active_c.insert(out.active_c.end(), halo_c.begin(), halo_c.end());
    out.active_f = own_f; out.active_f.insert(out.active_f.end(), halo_f.begin(), halo_f.end());

    out.Kcc     = extract_block(m_Kcc,     out.active_c, out.active_c);
    out.Kcf     = extract_block(m_Kcf,     out.active_c, out.active_f);
    out.Kfc     = extract_block(m_Kfc,     out.active_f, out.active_c);
    out.Mc_inv  = extract_block(m_Mc_inv,  out.active_c, out.active_c);
    out.Sff_inv = extract_block(m_Sff_inv, out.active_f, out.active_f);

    out.Kcc_own     = extract_block(m_Kcc,     own_c, own_c);
    out.Kcf_own     = extract_block(m_Kcf,     own_c, own_f);
    out.Kfc_own     = extract_block(m_Kfc,     own_f, own_c);
    out.Mc_inv_own  = extract_block(m_Mc_inv,  own_c, own_c);
    out.Sff_inv_own = extract_block(m_Sff_inv, own_f, own_f);
    return out;
}

// Band i's "coarse role": produces a cubic-in-time Taylor polynomial
// w[0..3], valid over LOCAL time [0,dt] (dt = the length of the CURRENT
// recursion interval, tau=0 <-> the state y at this interval's start),
// approximating band i's own dofs' evolution PLUS its leak into
// neighbouring bands' interface dofs (via the halo in `blocks.active_f`)
// -- exactly what a descendant (finer) level needs to know about band i
// without ever calling band i's own B-operator again. y is the FULL,
// global-size current state (same convention as erk_weight_LTS_coarse).
// Output w[i] are full-length (m_n_c_dof+m_n_f_dof) vectors, zero
// outside blocks.active_c/active_f.
void erk_weight_LTS_coarse_v2(const Matrix<T, Dynamic, 1> &y,
                               const LTS_subblock_set &blocks,
                               std::vector<Matrix<T, Dynamic, 1>> &w,
                               const Matrix<T, Dynamic, 1> &Fn,
                               const Matrix<T, Dynamic, 1> &Fn12,
                               const Matrix<T, Dynamic, 1> &Fn1,
                               const T dt) const {

    const int nc = (int)blocks.active_c.size();
    const int nf = (int)blocks.active_f.size();
    const size_t N = y.rows();

    auto gather_c = [&](const Matrix<T,Dynamic,1>& v) -> Matrix<T,Dynamic,1> {
        Matrix<T,Dynamic,1> out(nc);
        for (int i = 0; i < nc; ++i) out(i) = v(blocks.active_c[i]);
        return out;
    };
    auto gather_f = [&](const Matrix<T,Dynamic,1>& v) -> Matrix<T,Dynamic,1> {
        Matrix<T,Dynamic,1> out(nf);
        for (int i = 0; i < nf; ++i) out(i) = v(m_n_c_dof + blocks.active_f[i]);
        return out;
    };
    auto scatter = [&](const Matrix<T,Dynamic,1>& kc_loc, const Matrix<T,Dynamic,1>& kf_loc) -> Matrix<T,Dynamic,1> {
        Matrix<T,Dynamic,1> out = Matrix<T,Dynamic,1>::Zero(N);
        for (int i = 0; i < nc; ++i) out(blocks.active_c[i])             = kc_loc(i);
        for (int i = 0; i < nf; ++i) out(m_n_c_dof + blocks.active_f[i]) = kf_loc(i);
        return out;
    };

    // B restricted to this band's blocks, y=0 (pure force term):
    // returns (kc_loc, kf_loc) local vectors, NOT scattered.
    auto B_zero_y_loc = [&](const Matrix<T,Dynamic,1>& Fc_loc,
                             Matrix<T,Dynamic,1>& kc_loc, Matrix<T,Dynamic,1>& kf_loc) {
        kc_loc = blocks.Mc_inv * Fc_loc;
        kf_loc = -blocks.Sff_inv * (blocks.Kfc * kc_loc);
    };
    // B restricted to this band's blocks, general y (local vectors in/out).
    auto B_loc = [&](const Matrix<T,Dynamic,1>& yc_loc, const Matrix<T,Dynamic,1>& yf_loc, const Matrix<T,Dynamic,1>& Fc_loc,
                      Matrix<T,Dynamic,1>& kc_loc, Matrix<T,Dynamic,1>& kf_loc) {
        kc_loc = blocks.Mc_inv * (Fc_loc - blocks.Kcc*yc_loc - blocks.Kcf*yf_loc);
        kf_loc = -blocks.Sff_inv * (blocks.Kfc * kc_loc);
    };

    Matrix<T,Dynamic,1> zeroFc = Matrix<T,Dynamic,1>::Zero(nc);

    // Zero out entries beyond this band's own dofs (own entries are always
    // first in the gathered/local vectors, by construction of active_c/f).
    auto mask_own_c = [&](const Matrix<T,Dynamic,1>& v) {
        Matrix<T,Dynamic,1> out = v;
        for (size_t i = blocks.n_own_c; i < (size_t)out.rows(); ++i) out(i) = 0;
        return out;
    };
    auto mask_own_f = [&](const Matrix<T,Dynamic,1>& v) {
        Matrix<T,Dynamic,1> out = v;
        for (size_t i = blocks.n_own_f; i < (size_t)out.rows(); ++i) out(i) = 0;
        return out;
    };

    // Quadratic Lagrange interpolation of F over [0,dt] (local band force,
    // RAW/unmasked -- matches erk_weight_LTS_coarse's F0/F1/F2, used as-is
    // in the B^i*yn / B^i*F "arg" chains below).
    Matrix<T,Dynamic,1> Fn_c   = gather_c(Fn),   Fn12_c = gather_c(Fn12),  Fn1_c = gather_c(Fn1);
    Matrix<T,Dynamic,1> F0 =  Fn_c;
    Matrix<T,Dynamic,1> F1 = (-3*Fn_c + 4*Fn12_c - Fn1_c) / dt;
    Matrix<T,Dynamic,1> F2 = ( 4*Fn_c - 8*Fn12_c + 4*Fn1_c) / (dt*dt);

    // BFn_i = B_zero_y(Fi) -- RAW F, feeds the arg1/arg2/arg3 chains
    Matrix<T,Dynamic,1> BFn_0_c, BFn_0_f, BFn_1_c, BFn_1_f, BFn_2_c, BFn_2_f;
    B_zero_y_loc(F0, BFn_0_c, BFn_0_f);
    B_zero_y_loc(F1, BFn_1_c, BFn_1_f);
    B_zero_y_loc(F2, BFn_2_c, BFn_2_f);

    // IPFi = Fi masked to this band's OWN dofs only (matches erk_weight_LTS_coarse's
    // IP(Fi)); BIPFn_i = B_zero_y(IPFi); MinvFi = its cell part -- the SEPARATE
    // force-correction term added directly inside compute_w (distinct from the
    // RAW-F-based BFn_i chain above).
    Matrix<T,Dynamic,1> IPF0 = mask_own_c(F0), IPF1 = mask_own_c(F1), IPF2 = mask_own_c(F2);
    Matrix<T,Dynamic,1> BIPFn_0_c, BIPFn_0_f, BIPFn_1_c, BIPFn_1_f, BIPFn_2_c, BIPFn_2_f;
    B_zero_y_loc(IPF0, BIPFn_0_c, BIPFn_0_f);
    B_zero_y_loc(IPF1, BIPFn_1_c, BIPFn_1_f);
    B_zero_y_loc(IPF2, BIPFn_2_c, BIPFn_2_f);
    const Matrix<T,Dynamic,1>& MinvF0 = BIPFn_0_c;
    const Matrix<T,Dynamic,1>& MinvF1 = BIPFn_1_c;
    const Matrix<T,Dynamic,1>& MinvF2 = BIPFn_2_c;

    // B^i * yn chains, restricted to this band (UNMASKED state, matches
    // erk_weight_LTS_coarse's B0yn=y, B1yn=B(y), ...)
    Matrix<T,Dynamic,1> y0c = gather_c(y), y0f = gather_f(y);
    Matrix<T,Dynamic,1> B1c, B1f, B2c, B2f, B3c, B3f;
    B_loc(y0c, y0f, zeroFc, B1c, B1f);
    B_loc(B1c, B1f, zeroFc, B2c, B2f);
    B_loc(B2c, B2f, zeroFc, B3c, B3f);

    // B^i * F0 chains: BF0=B(BFn_0), B2F0=B(BF0), BF1=B(BFn_1)
    Matrix<T,Dynamic,1> BF0_c, BF0_f, B2F0_c, B2F0_f, BF1_c, BF1_f;
    B_loc(BFn_0_c, BFn_0_f, zeroFc, BF0_c, BF0_f);
    B_loc(BF0_c,   BF0_f,   zeroFc, B2F0_c, B2F0_f);
    B_loc(BFn_1_c, BFn_1_f, zeroFc, BF1_c, BF1_f);

    // Taylor arguments (local, cell+face parts) -- matches erk_weight_LTS_coarse:
    //   arg0 = B0yn
    //   arg1 = B1yn + BFn_0
    //   arg2 = B2yn + BF0  + BFn_1
    //   arg3 = B3yn + B2F0 + BF1  + BFn_2
    Matrix<T,Dynamic,1> arg0_c = y0c,                          arg0_f = y0f;
    Matrix<T,Dynamic,1> arg1_c = B1c + BFn_0_c,                arg1_f = B1f + BFn_0_f;
    Matrix<T,Dynamic,1> arg2_c = B2c + BF0_c  + BFn_1_c,       arg2_f = B2f + BF0_f  + BFn_1_f;
    Matrix<T,Dynamic,1> arg3_c = B3c + B2F0_c + BF1_c + BFn_2_c, arg3_f = B3f + B2F0_f + BF1_f + BFn_2_f;

    // wi = B((I-P)*arg) + (I-P)*Fi, matching erk_weight_LTS_coarse's compute_w:
    // mask arg down to this band's OWN dofs (zero halo) before the final B
    // application -- the halo only exists to make the ABOVE chains exact;
    // the halo cells' own dynamics are someone else's (a different band's)
    // responsibility. wi_c is then trivially zero beyond own_c (Kcc is
    // diagonal-per-cell and both arg_c/arg_f/MinvFext are zero there), so
    // scattering wi_f over the FULL active_f (own+halo) correctly captures
    // this band's leak into a neighbouring band's interface face (halo
    // faces adjacent to an own cell pick up a nonzero Kfc*wi_c row) without
    // any extra bookkeeping.
    auto compute_w = [&](const Matrix<T,Dynamic,1>& arg_c, const Matrix<T,Dynamic,1>& arg_f,
                          const Matrix<T,Dynamic,1>* MinvFext, Matrix<T,Dynamic,1>& wi) {
        Matrix<T,Dynamic,1> IParg_c = mask_own_c(arg_c);
        Matrix<T,Dynamic,1> IParg_f = mask_own_f(arg_f);
        Matrix<T,Dynamic,1> wi_c = blocks.Mc_inv * (-blocks.Kcc*IParg_c - blocks.Kcf*IParg_f);
        if (MinvFext) wi_c += *MinvFext;
        Matrix<T,Dynamic,1> wi_f = -blocks.Sff_inv * (blocks.Kfc * wi_c);
        wi = scatter(wi_c, wi_f);
    };

    compute_w(arg0_c, arg0_f, &MinvF0, w[0]);
    compute_w(arg1_c, arg1_f, &MinvF1, w[1]);
    compute_w(arg2_c, arg2_f, &MinvF2, w[2]);
    compute_w(arg3_c, arg3_f, nullptr,  w[3]);
}

// Band i's "fine role": genuine RK4 sub-step of size dtau, mirroring
// erk_weight_LTS_fine exactly (just restricted): each stage is
// B(P_own * y_stage) + Taylor_w(tau), where P_own masks the INPUT down to
// this band's own dofs (matching Pfine*y_stage in the original) before
// applying B with the HALO-INCLUSIVE operators -- crucially, B(P_own*y)
// is generally NONZERO not just at this band's own rows but also at
// halo (neighbouring band) rows that couple to an own face (a coarse
// cell bordering the interface picks up a real contribution from the
// fine side's Kcf coupling to that shared face, exactly as the original,
// unmasked erk_weight_LTS_fine's output naturally is over the WHOLE
// vector). So the RK4 update is accumulated over active_c/active_f
// (own + halo), not just own -- this is what fixes Bug A properly (an
// earlier attempt that only wrote "own" positions, compensating with a
// separate closed-form integral for the other band, was WRONG: it missed
// exactly this halo leak, confirmed by a direct side-by-side numerical
// mismatch against the trusted unrestricted algorithm).
void erk_weight_LTS_fine_v2(Matrix<T, Dynamic, 1> &x_dof_n,
                             const LTS_subblock_set &blocks,
                             const std::vector<Matrix<T, Dynamic, 1>> &w,
                             const Matrix<T, Dynamic, 1> &Fm,
                             const Matrix<T, Dynamic, 1> &Fmh,
                             const Matrix<T, Dynamic, 1> &Fm1,
                             const T tm,
                             const T dtau) const {

    const int nc = (int)blocks.active_c.size();
    const int nf = (int)blocks.active_f.size();

    auto Taylor_w = [&](T tau) -> Matrix<T,Dynamic,1> {
        T tau2 = tau*tau, tau3 = tau*tau2;
        Matrix<T,Dynamic,1> out(nc + nf);
        for (int i = 0; i < nc; ++i) {
            int g = blocks.active_c[i];
            out(i) = w[0](g) + tau*w[1](g) + (tau2/2.0)*w[2](g) + (tau3/6.0)*w[3](g);
        }
        for (int i = 0; i < nf; ++i) {
            size_t g = m_n_c_dof + blocks.active_f[i];
            out(nc+i) = w[0](g) + tau*w[1](g) + (tau2/2.0)*w[2](g) + (tau3/6.0)*w[3](g);
        }
        return out;
    };

    T tmh = tm + 0.5*dtau;
    T tm1 = tm + dtau;

    // gather own+halo local state from x_dof_n (own+halo -- matches B's
    // full active_c/active_f domain; masking to own-only happens INSIDE
    // fine_stage, mirroring Pfine*y_stage in the original)
    auto gather = [&](const Matrix<T,Dynamic,1>& full) -> Matrix<T,Dynamic,1> {
        Matrix<T,Dynamic,1> out(nc + nf);
        for (int i = 0; i < nc; ++i) out(i)      = full(blocks.active_c[i]);
        for (int i = 0; i < nf; ++i) out(nc + i) = full(m_n_c_dof + blocks.active_f[i]);
        return out;
    };
    auto mask_own = [&](const Matrix<T,Dynamic,1>& v) {
        Matrix<T,Dynamic,1> out = v;
        for (size_t i = blocks.n_own_c; i < (size_t)nc; ++i) out(i) = 0;
        for (size_t i = blocks.n_own_f; i < (size_t)nf; ++i) out(nc + i) = 0;
        return out;
    };

    auto fine_stage = [&](const Matrix<T,Dynamic,1>& y_stage, T tau, const Matrix<T,Dynamic,1>& F_tau) -> Matrix<T,Dynamic,1> {
        Matrix<T,Dynamic,1> Py = mask_own(y_stage);
        Matrix<T,Dynamic,1> PF = mask_own(gather(F_tau));
        Matrix<T,Dynamic,1> yc = Py.head(nc), yf = Py.tail(nf), Fc_loc = PF.head(nc);
        Matrix<T,Dynamic,1> kc = blocks.Mc_inv * (Fc_loc - blocks.Kcc*yc - blocks.Kcf*yf);
        Matrix<T,Dynamic,1> kf = -blocks.Sff_inv * (blocks.Kfc * kc);
        Matrix<T,Dynamic,1> k(nc+nf);
        k.head(nc) = kc; k.tail(nf) = kf;
        return k + Taylor_w(tau);
    };

    Matrix<T,Dynamic,1> y0 = gather(x_dof_n);
    Matrix<T,Dynamic,1> k1 = fine_stage(y0,                tm,  Fm);
    Matrix<T,Dynamic,1> k2 = fine_stage(y0 + 0.5*dtau*k1, tmh, Fmh);
    Matrix<T,Dynamic,1> k3 = fine_stage(y0 + 0.5*dtau*k2, tmh, Fmh);
    Matrix<T,Dynamic,1> k4 = fine_stage(y0 +     dtau*k3, tm1, Fm1);

    Matrix<T,Dynamic,1> dy = dtau * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    for (int i = 0; i < nc; ++i) x_dof_n(blocks.active_c[i])             += dy(i);
    for (int i = 0; i < nf; ++i) x_dof_n(m_n_c_dof + blocks.active_f[i]) += dy(nc + i);
}

// Companion to erk_weight_LTS_coarse_v2 + erk_weight_LTS_fine_v2: advances
// `blocks_coarse`'s own dofs that are NOT reached by the descendant (fine)
// band's halo (`descendant_active_c`/`descendant_active_f`, i.e. its
// LTS_subblock_set::active_c/active_f, own+halo) via the closed-form
// integral of wk over [0,dt]: wk[0]*dt + wk[1]*dt^2/2 + wk[2]*dt^3/6 +
// wk[3]*dt^4/24. This is EXACT, not an approximation of the exact
// algorithm: for a cell with every face coarse-classified (no coupling to
// any fine-region face), B(Pfine*y)=0 identically, so its whole-macro-step
// update is exactly this integral of its own cubic Taylor polynomial
// (Simpson's rule -- which is what the RK4 update at such positions
// reduces to -- is exact for a cubic).
//
// Dofs the descendant's halo DOES reach must NOT also get this update:
// erk_weight_LTS_fine_v2's own scatter already includes both the Taylor
// term (via `w` argument, i.e. this same wk) AND the genuine B(Pfine*y)
// coupling correction for them -- adding this integral there too would
// double-count the Taylor part. This split (halo-covered dofs handled by
// the descendant's own scatter; only genuinely-uncovered dofs need this
// separate step) was found necessary and sufficient by a direct,
// validated side-by-side comparison against the original 2-level scheme
// (see ERK4_LTS_v2_L2_validation_test.hpp) -- an earlier attempt applying
// this integral to ALL of a band's own dofs unconditionally, and a
// separate attempt relying solely on the descendant's scatter with no
// fallback at all, were both numerically wrong.
void erk_weight_LTS_coarse_advance_uncovered(Matrix<T, Dynamic, 1> &x_dof_n,
                                              const LTS_subblock_set &blocks_coarse,
                                              const std::vector<Matrix<T, Dynamic, 1>> &wk,
                                              const std::vector<int> &descendant_active_c,
                                              const std::vector<int> &descendant_active_f,
                                              const T dt) const {
    std::set<int> cov_c(descendant_active_c.begin(), descendant_active_c.end());
    std::set<int> cov_f(descendant_active_f.begin(), descendant_active_f.end());
    T dt2 = dt*dt, dt3 = dt2*dt, dt4 = dt3*dt;
    for (size_t i = 0; i < blocks_coarse.n_own_c; ++i) {
        int g = blocks_coarse.active_c[i];
        if (cov_c.count(g)) continue;
        x_dof_n(g) += wk[0](g)*dt + wk[1](g)*dt2/2.0 + wk[2](g)*dt3/6.0 + wk[3](g)*dt4/24.0;
    }
    for (size_t i = 0; i < blocks_coarse.n_own_f; ++i) {
        int local = blocks_coarse.active_f[i];
        if (cov_f.count(local)) continue;
        int g = (int)m_n_c_dof + local;
        x_dof_n(g) += wk[0](g)*dt + wk[1](g)*dt2/2.0 + wk[2](g)*dt3/6.0 + wk[3](g)*dt4/24.0;
    }
}


};




#endif /* erk_hho_scheme_hpp */
