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
    
    SparseMatrix<T> m_Mcc_coarse;
    SparseMatrix<T> m_Kcc_coarse;
    SparseMatrix<T> m_Kcf_coarse;
    SparseMatrix<T> m_Mcc_fine;
    SparseMatrix<T> m_Kcc_fine;
    SparseMatrix<T> m_Kcf_fine;
    SparseMatrix<T> m_Kfc_fine;
    SparseMatrix<T> m_Sff_fine;
    SparseMatrix<T> m_Sff_inv_fine;  
    SparseMatrix<T> m_Mc_inv_fine;  
    
    std::vector<size_t> m_fine_c_indices;
    std::vector<size_t> m_fine_f_indices;
    std::vector<size_t> m_coarse_c_indices;
    std::vector<size_t> m_coarse_f_indices;
    


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
    
    void erk_weight_LTS_coarse(const Matrix<T, Dynamic, 1> &y, const Eigen::SparseMatrix<double> &Pcoarse, std::vector<Matrix<T, Dynamic, 1>> &w, const Matrix<T, Dynamic, 1> &Fn, const Matrix<T, Dynamic, 1> &Fn12, const Matrix<T, Dynamic, 1> &Fn1, const T dt) {
        
        // Coefficients Lagrange quadratique
        Matrix<T, Dynamic, 1> F0 =  Fn;
        Matrix<T, Dynamic, 1> F1 = (-3*Fn + 4*Fn12 - Fn1) / dt;
        Matrix<T, Dynamic, 1> F2 = ( 4*Fn - 8*Fn12 + 4*Fn1) / (dt*dt);
        
        // (I-P)Fi pour les termes extérieurs des w
        Matrix<T, Dynamic, 1> IPF0 = Pcoarse * F0;
        Matrix<T, Dynamic, 1> IPF1 = Pcoarse * F1;
        Matrix<T, Dynamic, 1> IPF2 = Pcoarse * F2;
        
        // F complet pour les chaînes B^i F
        // (I-P)F pour les termes extérieurs des w
        Matrix<T, Dynamic, 1> zero = Matrix<T, Dynamic, 1>::Zero(y.rows());
        Matrix<T, Dynamic, 1> F0_full, F1_full, F2_full;
        Matrix<T, Dynamic, 1> IPF0_full, IPF1_full, IPF2_full;
        
        SetFg(F0);   erk_weight(zero, F0_full);   ZeroFc();
        SetFg(F1);   erk_weight(zero, F1_full);   ZeroFc();
        SetFg(F2);   erk_weight(zero, F2_full);   ZeroFc();
        SetFg(IPF0); erk_weight(zero, IPF0_full); ZeroFc();
        SetFg(IPF1); erk_weight(zero, IPF1_full); ZeroFc();
        SetFg(IPF2); erk_weight(zero, IPF2_full); ZeroFc();
        
        auto Pc = Pcoarse.block(0, 0, m_n_c_dof, m_n_c_dof);
        Matrix<T, Dynamic, 1> MinvF0 = IPF0_full.block(0, 0, m_n_c_dof, 1);
        Matrix<T, Dynamic, 1> MinvF1 = IPF1_full.block(0, 0, m_n_c_dof, 1);
        Matrix<T, Dynamic, 1> MinvF2 = IPF2_full.block(0, 0, m_n_c_dof, 1);
        
        // Chaînes B^i y (homogène, Fc = 0)
        Matrix<T, Dynamic, 1> B0y = y;
        Matrix<T, Dynamic, 1> B1y, B2y, B3y;
        erk_weight(B0y, B1y);
        erk_weight(B1y, B2y);
        erk_weight(B2y, B3y);
        
        // Chaînes B^i F — F complet dans l'argument
        Matrix<T, Dynamic, 1> BF0, B2F0, BF1;
        erk_weight(F0_full, BF0);
        erk_weight(BF0,     B2F0);
        erk_weight(F1_full, BF1);
        
        // Arguments de B(I-P)(·)
        Matrix<T, Dynamic, 1> arg0 = B0y;
        Matrix<T, Dynamic, 1> arg1 = B1y + F0_full;
        Matrix<T, Dynamic, 1> arg2 = B2y + BF0  + F1_full;
        Matrix<T, Dynamic, 1> arg3 = B3y + B2F0 + BF1 + F2_full;
        
        auto compute_one_w = [&](const Matrix<T, Dynamic, 1> &arg, const Matrix<T, Dynamic, 1> *MinvFext, Matrix<T, Dynamic, 1> &wi) {
            
            Matrix<T, Dynamic, 1> Ptmp   = Pcoarse * arg;
            Matrix<T, Dynamic, 1> Ptmp_c = Ptmp.block(0, 0, m_n_c_dof, 1);
            Matrix<T, Dynamic, 1> Ptmp_f = Ptmp.block(m_n_c_dof, 0, m_n_f_dof, 1);
            
            Matrix<T, Dynamic, 1> wi_c = m_Mc_inv * (-Kcc()*Ptmp_c - Kcf()*Ptmp_f);
            
            if (MinvFext)
            wi_c += Pc * (*MinvFext);
            
            wi = Ptmp;
            wi.block(0, 0, m_n_c_dof, 1) = wi_c;
            Matrix<T, Dynamic, 1> RHSf = Kfc() * wi_c;
            if (m_sff_is_block_diagonal_Q)
            wi.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_Sff_inv * RHSf;
            else
            wi.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_inv_Sff * RHSf;
        };
        
        compute_one_w(arg0, &MinvF0, w[0]);
        compute_one_w(arg1, &MinvF1, w[1]);
        compute_one_w(arg2, &MinvF2, w[2]);
        compute_one_w(arg3, nullptr,  w[3]);
    }
    
    void erk_weight_LTS_fine(Matrix<T, Dynamic, 1> &x_dof_n, const Eigen::SparseMatrix<T> &Pfine, const std::vector<Matrix<T, Dynamic, 1>> &w, const Matrix<T, Dynamic, 1> &Fm, const Matrix<T, Dynamic, 1> &Fmh, const Matrix<T, Dynamic, 1> &Fm1, const T tm, const T dtau) {
        
        auto Taylor_w = [&](T tau) -> Matrix<T, Dynamic, 1> {
            T tau2 = tau*tau, tau3 = tau*tau2;
            return w[0] + tau*w[1] + (tau2/2)*w[2] + (tau3/6)*w[3];
        };
        
        auto fine_stage = [&](const Matrix<T, Dynamic, 1> &y_stage, T tau, const Matrix<T, Dynamic, 1> &F_tau) -> Matrix<T, Dynamic, 1> {
            Matrix<T, Dynamic, 1> Py = Pfine * y_stage;
            Matrix<T, Dynamic, 1> PF = Pfine * F_tau;
            SetFg(PF);
            Matrix<T, Dynamic, 1> k;
            erk_weight(Py, k);
            ZeroFc();
            k += Taylor_w(tau);
            return k;
        };
        
        T tmh = tm + 0.5*dtau;
        T tm1 = tm +     dtau;
        
        Matrix<T, Dynamic, 1> k0 = fine_stage(x_dof_n,                tm,  Fm);
        Matrix<T, Dynamic, 1> k1 = fine_stage(x_dof_n + 0.5*dtau*k0, tmh, Fmh);
        Matrix<T, Dynamic, 1> k2 = fine_stage(x_dof_n + 0.5*dtau*k1, tmh, Fmh);
        Matrix<T, Dynamic, 1> k3 = fine_stage(x_dof_n +     dtau*k2, tm1, Fm1);
        
        x_dof_n += dtau * (k0 + 2*k1 + 2*k2 + k3) / 6;
    }
    
    void build_fine_submatrices() {
        
        size_t nfc = m_fine_c_indices.size();
        size_t nff = m_fine_f_indices.size();
        
        // Sous-matrices fines extraites depuis Kcc, Kcf, Kfc, Sff_inv, Mc_inv
        // Kcc_ff : lignes et colonnes fines de Kcc
        // Kcf_ff : lignes fines de Kcc, colonnes fines de Kcf
        // etc.
        
        // Construction par triplets
        std::vector<Triplet<T>> trips_Kcc, trips_Kcf, trips_Kfc, trips_Sff, trips_Mc;
        
        // Map global → local pour les indices fins
        std::vector<int> c_global_to_local(m_n_c_dof, -1);
        std::vector<int> f_global_to_local(m_n_f_dof, -1);
        for (size_t li = 0; li < nfc; ++li) c_global_to_local[m_fine_c_indices[li]] = li;
        for (size_t li = 0; li < nff; ++li) f_global_to_local[m_fine_f_indices[li] - m_n_c_dof] = li;
        
        // Kcc_fine : nfc × nfc
        for (int k = 0; k < m_Kcc.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Kcc, k); it; ++it) {
                int r = c_global_to_local[it.row()];
                int c = c_global_to_local[it.col()];
                if (r >= 0 && c >= 0)
                trips_Kcc.emplace_back(r, c, it.value());
            }
        }
        m_Kcc_fine.resize(nfc, nfc);
        m_Kcc_fine.setFromTriplets(trips_Kcc.begin(), trips_Kcc.end());
        
        // Kcf_fine : nfc × nff
        for (int k = 0; k < m_Kcf.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Kcf, k); it; ++it) {
                int r = c_global_to_local[it.row()];
                int c = f_global_to_local[it.col()];
                if (r >= 0 && c >= 0)
                trips_Kcf.emplace_back(r, c, it.value());
            }
        }
        m_Kcf_fine.resize(nfc, nff);
        m_Kcf_fine.setFromTriplets(trips_Kcf.begin(), trips_Kcf.end());
        
        // Kfc_fine : nff × nfc
        for (int k = 0; k < m_Kfc.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Kfc, k); it; ++it) {
                int r = f_global_to_local[it.row()];
                int c = c_global_to_local[it.col()];
                if (r >= 0 && c >= 0)
                trips_Kfc.emplace_back(r, c, it.value());
            }
        }
        m_Kfc_fine.resize(nff, nfc);
        m_Kfc_fine.setFromTriplets(trips_Kfc.begin(), trips_Kfc.end());
        
        // Sff_inv_fine : nff × nff
        for (int k = 0; k < m_Sff_inv.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Sff_inv, k); it; ++it) {
                int r = f_global_to_local[it.row()];
                int c = f_global_to_local[it.col()];
                if (r >= 0 && c >= 0)
                trips_Sff.emplace_back(r, c, it.value());
            }
        }
        m_Sff_inv_fine.resize(nff, nff);
        m_Sff_inv_fine.setFromTriplets(trips_Sff.begin(), trips_Sff.end());
        
        // Mc_inv_fine : nfc × nfc
        for (int k = 0; k < m_Mc_inv.outerSize(); ++k) {
            for (typename SparseMatrix<T>::InnerIterator it(m_Mc_inv, k); it; ++it) {
                int r = c_global_to_local[it.row()];
                int c = c_global_to_local[it.col()];
                if (r >= 0 && c >= 0)
                trips_Mc.emplace_back(r, c, it.value());
            }
        }
        m_Mc_inv_fine.resize(nfc, nfc);
        m_Mc_inv_fine.setFromTriplets(trips_Mc.begin(), trips_Mc.end());
    }
    
    void erk_weight_LTS_coarse_optimised(
        const Matrix<T, Dynamic, 1> &y,
        const Eigen::SparseMatrix<double> &Pcoarse,
        std::vector<Matrix<T, Dynamic, 1>> &w,
        const Matrix<T, Dynamic, 1> &Fn,
        const Matrix<T, Dynamic, 1> &Fn12,
        const Matrix<T, Dynamic, 1> &Fn1,
        const T dt) {

    Matrix<T, Dynamic, 1> F0 =  Fn;
    Matrix<T, Dynamic, 1> F1 = (-3*Fn + 4*Fn12 - Fn1) / dt;
    Matrix<T, Dynamic, 1> F2 = ( 4*Fn - 8*Fn12 + 4*Fn1) / (dt*dt);

    // apply_B_source : B(0) avec Fc=F — évite Kcc*0 et Kcf*0
    auto apply_B_source = [&](const Matrix<T, Dynamic, 1> &F,
                               Matrix<T, Dynamic, 1> &out) {
        out.resize(y.rows());
        out.setZero();
        Matrix<T, Dynamic, 1> k_c = m_Mc_inv * F.block(0, 0, m_n_c_dof, 1);
        out.block(0, 0, m_n_c_dof, 1) = k_c;
        Matrix<T, Dynamic, 1> RHSf = Kfc() * k_c;
        if (m_sff_is_block_diagonal_Q)
            out.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_Sff_inv * RHSf;
        else
            out.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_inv_Sff * RHSf;
    };

    Matrix<T, Dynamic, 1> F0_full, F1_full, F2_full;
    apply_B_source(F0, F0_full);
    apply_B_source(F1, F1_full);
    apply_B_source(F2, F2_full);

    // (I-P)F via produit matriciel — nécessaire car Kcc couple fins et grossiers
    Matrix<T, Dynamic, 1> MinvF0 = m_Mc_inv * (Pcoarse * F0).block(0, 0, m_n_c_dof, 1);
    Matrix<T, Dynamic, 1> MinvF1 = m_Mc_inv * (Pcoarse * F1).block(0, 0, m_n_c_dof, 1);
    Matrix<T, Dynamic, 1> MinvF2 = m_Mc_inv * (Pcoarse * F2).block(0, 0, m_n_c_dof, 1);

    Matrix<T, Dynamic, 1> B0y = y;
    Matrix<T, Dynamic, 1> B1y, B2y, B3y;
    erk_weight(B0y, B1y);
    erk_weight(B1y, B2y);
    erk_weight(B2y, B3y);

    Matrix<T, Dynamic, 1> BF0, B2F0, BF1;
    erk_weight(F0_full, BF0);
    erk_weight(BF0,     B2F0);
    erk_weight(F1_full, BF1);

    Matrix<T, Dynamic, 1> arg0 = B0y;
    Matrix<T, Dynamic, 1> arg1 = B1y + F0_full;
    Matrix<T, Dynamic, 1> arg2 = B2y + BF0  + F1_full;
    Matrix<T, Dynamic, 1> arg3 = B3y + B2F0 + BF1 + F2_full;

    auto compute_one_w = [&](const Matrix<T, Dynamic, 1> &arg,
                              const Matrix<T, Dynamic, 1> *MinvFext,
                              Matrix<T, Dynamic, 1> &wi) {

        // Pcoarse * arg — produit matriciel obligatoire
        Matrix<T, Dynamic, 1> Ptmp   = Pcoarse * arg;
        Matrix<T, Dynamic, 1> Ptmp_c = Ptmp.block(0, 0, m_n_c_dof, 1);
        Matrix<T, Dynamic, 1> Ptmp_f = Ptmp.block(m_n_c_dof, 0, m_n_f_dof, 1);

        Matrix<T, Dynamic, 1> wi_c = m_Mc_inv * (-Kcc()*Ptmp_c - Kcf()*Ptmp_f);

        if (MinvFext)
            wi_c += *MinvFext;

        wi = Ptmp;
        wi.block(0, 0, m_n_c_dof, 1) = wi_c;
        Matrix<T, Dynamic, 1> RHSf = Kfc() * wi_c;
        if (m_sff_is_block_diagonal_Q)
            wi.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_Sff_inv * RHSf;
        else
            wi.block(m_n_c_dof, 0, m_n_f_dof, 1) = -m_inv_Sff * RHSf;
    };

    compute_one_w(arg0, &MinvF0, w[0]);
    compute_one_w(arg1, &MinvF1, w[1]);
    compute_one_w(arg2, &MinvF2, w[2]);
    compute_one_w(arg3, nullptr,  w[3]);
}


void erk_weight_LTS_fine_optimised(
        Matrix<T, Dynamic, 1> &x_dof_n,
        const std::vector<Matrix<T, Dynamic, 1>> &w,
        const Matrix<T, Dynamic, 1> &Fm,
        const Matrix<T, Dynamic, 1> &Fmh,
        const Matrix<T, Dynamic, 1> &Fm1,
        const T tm,
        const T dtau) {

    size_t nfc = m_fine_c_indices.size();
    size_t nff = m_fine_f_indices.size();

    // Extraction des composantes fines d'un vecteur global
    auto extract_fine = [&](const Matrix<T, Dynamic, 1> &v)
                         -> std::pair<Matrix<T,Dynamic,1>, Matrix<T,Dynamic,1>> {
        Matrix<T, Dynamic, 1> vc(nfc), vf(nff);
        for (size_t i = 0; i < nfc; ++i) vc(i) = v(m_fine_c_indices[i]);
        for (size_t i = 0; i < nff; ++i) vf(i) = v(m_fine_f_indices[i]);
        return {vc, vf};
    };

    // Taylor_w sur les indices fins seulement
    auto Taylor_w_fine = [&](T tau) -> std::pair<Matrix<T,Dynamic,1>, Matrix<T,Dynamic,1>> {
        T tau2 = tau*tau, tau3 = tau*tau2;
        Matrix<T, Dynamic, 1> tw = w[0] + tau*w[1] + (tau2/2)*w[2] + (tau3/6)*w[3];
        return extract_fine(tw);
    };

    // Un stage fin — tout en taille réduite nfc/nff
    auto fine_stage = [&](const Matrix<T, Dynamic, 1> &y_global,
                           T tau,
                           const Matrix<T, Dynamic, 1> &F_tau)
                       -> Matrix<T, Dynamic, 1> {

        auto [yc, yf] = extract_fine(y_global);
        auto [Fc, Ff] = extract_fine(F_tau);
        auto [twc, twf] = Taylor_w_fine(tau);

        // k_c = Mc_inv_fine * (Fc - Kcc_fine*yc - Kcf_fine*yf)
        Matrix<T, Dynamic, 1> k_c = m_Mc_inv_fine * (Fc - m_Kcc_fine*yc - m_Kcf_fine*yf);

        // k_f = -Sff_inv_fine * Kfc_fine * k_c
        Matrix<T, Dynamic, 1> k_f = -m_Sff_inv_fine * (m_Kfc_fine * k_c);

        // Taylor_w ajouté sur les composantes fines
        k_c += twc;
        k_f += twf;

        // Reconstruction du vecteur global — DDL grossiers inchangés
        Matrix<T, Dynamic, 1> k = Matrix<T, Dynamic, 1>::Zero(y_global.rows());
        for (size_t i = 0; i < nfc; ++i) k(m_fine_c_indices[i]) = k_c(i);
        for (size_t i = 0; i < nff; ++i) k(m_fine_f_indices[i]) = k_f(i);
        return k;
    };

    T tmh = tm + 0.5*dtau;
    T tm1 = tm +     dtau;

    Matrix<T, Dynamic, 1> k0 = fine_stage(x_dof_n,                tm,  Fm);
    Matrix<T, Dynamic, 1> k1 = fine_stage(x_dof_n + 0.5*dtau*k0, tmh, Fmh);
    Matrix<T, Dynamic, 1> k2 = fine_stage(x_dof_n + 0.5*dtau*k1, tmh, Fmh);
    Matrix<T, Dynamic, 1> k3 = fine_stage(x_dof_n +     dtau*k2, tm1, Fm1);

    x_dof_n += dtau * (k0 + 2*k1 + 2*k2 + k3) / 6;
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
    
    
};




#endif /* erk_hho_scheme_hpp */
