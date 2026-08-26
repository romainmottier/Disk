#include <iostream>
#include <libgen.h>
#include <libgen.h>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <list>
#include <getopt.h>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
using namespace Eigen;

#include "diskpp/common/timecounter.hpp"
#include "diskpp/methods/hho"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/output/silo.hpp"

// application common sources
#include "common/display_settings.hpp"
#include "common/fitted_geometry_builders.hpp"
#include "common/linear_solver.hpp"
#include "common/acoustic_material_data.hpp"
#include "common/elastic_material_data.hpp"
#include "common/scal_vec_analytic_functions.hpp"
#include "common/preprocessor.hpp"
#include "common/postprocessor.hpp"

// RK schemes
#include "common/dirk_hho_scheme.hpp"
#include "common/dirk_butcher_tableau.hpp"
#include "common/erk_butcher_tableau.hpp"
#include "common/erk_hho_scheme.hpp"
#include "common/erk_coupling_hho_scheme.hpp"

// PROTOTYPES:

   // Computation of an empirical CFL criteria                    
   #include "prototypes/acoustic/EAcoustic_CFL.hpp"                   // CFl - Acoustic                      
   #include "prototypes/elastic/EElasticity_CFL.hpp"                  // CFl - Linear Elasticity  
   #include "prototypes/coupling/CFL/EHHOFirstOrderCFL.hpp"           // CFl - Elasto-Acoustic Coupling 

   // Stability study & Spectral radius computation:
   #include "prototypes/acoustic/EAcoustic_stability.hpp"             // Acoustic
   #include "prototypes/elastic/EElastic_stability.hpp"               // Linear Elasticity
   #include "prototypes/coupling/EHHOFirstOrder_stability.hpp"        // Elasto-Acoustic Coupling                   
   
   // Convergence test on sinusoidal analytical solution 
   #include "prototypes/acoustic/EAcoustic_conv_test.hpp"             // Explicit Acoustic               
   #include "prototypes/acoustic/IAcoustic_conv_test.hpp"             // Implicit Acoustic               
   #include "prototypes/elastic/IElastic_conv_test.hpp"               // Implicit Elastic 
   #include "prototypes/coupling/Conv_Tests/IHHOFirstOrder.hpp"                // Implicit Coupling                         
   #include "prototypes/coupling/Conv_Tests/IHHOFirstOrder_conv_tests.hpp"     // Explicit Coupling    
   #include "prototypes/coupling/Conv_Tests/EHHOFirstOrder.hpp"                // Explicit Coupling    
   #include "prototypes/coupling/Conv_Tests/EHHOFirstOrder_conv_tests.hpp"     // Explicit Coupling    

   // Pulses for comparison with Gar6more  
   #include "prototypes/coupling/Pulse/HeterogeneousIHHOFirstOrder.hpp"        // Implicit Pulse (adimensional)
   #include "prototypes/coupling/Pulse/HeterogeneousEHHOFirstOrder.hpp"        // Explicit Pulse (adimensional)
   #include "prototypes/coupling/Pulse/ConicWavesIHHOFirstOrder.hpp"           // Implicit Pulse (geophysic) 
   #include "prototypes/coupling/Pulse/review_CMAME.hpp"           // Implicit Pulse (geophysic) 
   #include "prototypes/coupling/Pulse/ConicWavesEHHOFirstOrder.hpp"           // Implicit Pulse (geophysic) 

   // Sedimentary Basin
   #include "prototypes/coupling/Basin/BassinIHHOFirstOrder.hpp"               // Implicit Sedimentary Basin

   // Segmented Brain (MRI)
   #include "prototypes/coupling/Basin/BrainIHHOFirstOrder.hpp"                // Implicit Brain (MRI segmentation)
   #include "prototypes/coupling/Basin/BrainEHHOFirstOrder.hpp"                // Explicit Brain (MRI segmentation)

   // LTS
      // CONV TEST 
      #include "prototypes/LTS/ERK4_LTS.hpp"     
      #include "prototypes/LTS/ERK4_LTS_stab.hpp"     
      #include "prototypes/LTS/ERK4_LTS_stab_acou.hpp"
      #include "prototypes/LTS/LHS_spectrum_acou.hpp"
      #include "prototypes/LTS/ERK4_LTS_Lshape_conv_test.hpp"
      #include "prototypes/LTS/EllipticLshape_conv_test.hpp"
      #include "prototypes/LTS/ERK4_LTS_Lshape_MMS_conv_test.hpp"
      #include "prototypes/LTS/ERK4_LTS_v2_L2_validation_test.hpp"
      #include "prototypes/LTS/ERK4_LTS_v2_timing_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_timing_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_L2_validation_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_L3_validation_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_L5small_validation_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_L3_v2design_validation_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_L5small_v2design_validation_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_Lsweep_small_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_GlobalGD_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_RestrictedGD_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_DiagCompare_test.hpp"
      #include "prototypes/LTS/ERK4_MLTS_Lshape_conv_test.hpp"
      // -- multilevel LTS-RK4 investigation (2026), isolated below --
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_GlobalExact_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_RestrictedExact_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_RestrictedLeveled_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_MehlinExact_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_MehlinRestricted_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_FivePoint_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_HybridTerminal_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/lshape/ERK4_MLTS_Lshape_DiagN2_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_Square_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_Square_Pell_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_Square_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CornerSquare_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CornerSquareFixed_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CornerSquareRamp_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CornerSquareCoarse_PellSparse_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CenterSquareCoarse_PellSparse_conv_test.hpp"
      // -- PLOT ARTICLES: prototypes behind the published figures, see main()'s "PLOT ARTICLES" section --
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CenterSquare_TwoLevel_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CenterSquare_ClassicalTwoLevel_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CenterSquare_FixedLevel_conv_test.hpp"
      #include "prototypes/LTS/mlts_2026/square/ERK4_MLTS_CenterSquare_Adaptive_conv_test.hpp"
      #include "prototypes/LTS/ERK4_LTS_SSTAB.hpp"
      #include "prototypes/LTS/ERK4_LTS_conv_test.hpp"    
      #include "prototypes/LTS/ERK4_LTS_optimised.hpp"    
      #include "prototypes/LTS/ERK_LTS_stab.hpp" 
      // PULSE
      #include "prototypes/LTS/AcousticHeterogeneousPulse.hpp"           // ACOUSTIC ERK
      #include "prototypes/LTS/AcousticLTSEulerHeterogeneousPulse.hpp"   // ACOUSTIC EULER-LTS
      #include "prototypes/LTS/AcousticHeterogeneousPulse_LTS_RK4.hpp"   // ACOUSTIC ERK4-LTS
      #include "prototypes/LTS/HeterogeneousERK4_LTS_HHO_FirstOrder.hpp" // COUPLING ERK4-LTS
      #include "prototypes/LTS/HeterogeneousERK4_LTS_HHO_FirstOrder_stab.hpp" // COUPLING ERK4-LTS
      
int main(int argc, char **argv){

    DBSetDeprecateWarnings(0);
    
    // REGRESSION TESTS
    if (basename(argv[0]) == std::string("name1") ) {
        std::cout << "called with name1" << std::endl;
        return 0;
    }

    if (basename(argv[0]) == std::string("name2") ) {
        std::cout << "called with name2" << std::endl;
        return 0;
    }
    
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// HHO ERK ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// CFL TABLE:
   // EAcoustic_CFL(argc, argv); 
   // EElasticity_CFL(argc, argv);
   // EHHOFirstOrderCFL(argc, argv); 
 
// STABILITY STUDY & SPECTRAL RADIUS COMPUTATION:
   // EAcoustic_stability(argc, argv);
   // EElastic_stability(argc, argv);
   // EHHOFirstOrder_stability(argc, argv); 

// CV TESTS:
   // EAcousticFirstOrder(argc, argv);
   // IAcoustic_conv_test(argc, argv);
   // IElastic_conv_test(argc, argv);
   // IHHOFirstOrder(argc, argv);
   // IHHOFirstOrder_conv_tests(argc, argv);
   // EHHOFirstOrder(argc, argv);
   // EHHOFirstOrder_conv_tests(argc, argv);

// PULSE: 
   // HeterogeneousIHHOFirstOrder(argc, argv); 
   // HeterogeneousEHHOFirstOrder(argc, argv); 
   // ConicWavesIHHOFirstOrder(argc, argv);
   // ConicWavesIHHOFirstOrder_review(argc, argv);
   // ConicWavesEHHOFirstOrder(argc, argv);
   // ConicWavesEHHOFirstOrder_review(argc, argv);

// SEDIMENTARY BASIN:
   // BassinIHHOFirstOrder(argc, argv);
   // Test(argc, argv);
   // BassinEHHOFirstOrder(argc, argv); Not working

// SEGMENTED BRAIN (MRI):
   // BrainIHHOFirstOrder(argc, argv);
   // BrainEHHOFirstOrder(argc, argv);

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////// PLOT ARTICLES ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

   //   1 -
            // ERK4_MLTS_CenterSquare_ClassicalTwoLevel_conv_test(argc, argv);
   
   //   2 - Multilevel (Diaz-Grote-style Pell recursion), three designs:
   //    (A) p fixed & large enough for genuine L>=3, L left NATURAL
   //        (unforced): MLTS_PFAM=6 (p=64->L=3), 8 (p=256->L=4), 10 (p=1024->L=5)
   //    (B) p FIXED at 1024 (MLTS_PFAM=10), L FORCED explicitly via
   //        MLTS_FORCE_L=1,2,3,4,5 -- isolates the effect of L alone
   //        at constant mesh ratio
             ERK4_MLTS_CenterSquare_FixedLevel_conv_test(argc, argv);
   //
   //    (C) p variable (grows as 4^N via the centersquareadaptive_graded
   //        mesh family), L left NATURAL/unforced -- an intelligent
   //        p->L relation via the same threshold rule as (A)/(B), giving
   //        L=1,1,2,3,4 across N=0..4:
             // ERK4_MLTS_CenterSquare_Adaptive_conv_test(argc, argv);

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////// PLOT ARTICLES ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


// LOCAL TIME STEPPING
   // HeterogeneousEULER_LTS_HHO_FirstOrder(argc, argv);

   // TEST LTS PULSE:
   // AcousticHeterogeneousPulse(argc, argv);
   // AcousticLTSEulerHeterogeneousPulse(argc, argv);
   // AcousticHeterogeneousPulse_LTS_RK4(argc, argv); 

   // EHHOFirstOrder(argc, argv);
   // ERK4_LTS(argc, argv);
   // ERK4_LTS_conv_test(argc, argv);
   // ERK4_LTS_stab(argc, argv);
   // ERK4_LTS_stab_acou(argc, argv);
   // LHS_spectrum_acou(argc, argv);
   // ERK4_LTS_Lshape_conv_test(argc, argv);
   // EllipticLshape_conv_test(argc, argv);
   // ERK4_LTS_Lshape_MMS_conv_test(argc, argv);
   // ERK4_LTS_v2_L2_validation_test(argc, argv);
   // ERK4_LTS_v2_timing_test(argc, argv);
   // ERK4_MLTS_timing_test(argc, argv);
   // ERK4_MLTS_L2_validation_test(argc, argv);
   // ERK4_MLTS_L3_validation_test(argc, argv);
   // ERK4_MLTS_L5small_validation_test(argc, argv);
   // ERK4_MLTS_L3_v2design_validation_test(argc, argv);
   // ERK4_MLTS_timing_test(argc, argv);
   // ERK4_MLTS_Lshape_conv_test(argc, argv);
   //
   // -- multilevel LTS-RK4 investigation (2026), see prototypes/LTS/mlts_2026/ --
   // L-shape (mlts_2026/lshape/):
   // ERK4_MLTS_Lshape_GlobalExact_conv_test(argc, argv);       // proven-exact reference, O(n_dof)/call, no truncation
   // ERK4_MLTS_Lshape_RestrictedExact_conv_test(argc, argv);   // own+halo submatrix restriction (~5-9x residual)
   // ERK4_MLTS_Lshape_RestrictedLeveled_conv_test(argc, argv); // + P_ell masking on top of the halo (refuted: no better)
   // ERK4_MLTS_Lshape_MehlinExact_conv_test(argc, argv);       // literal P_ell+overlap, GLOBAL matrices (matches GlobalExact <0.5%)
   // ERK4_MLTS_Lshape_MehlinRestricted_conv_test(argc, argv);  // P_ell+overlap, TRUNCATED submatrix (still ~5x residual)
   // ERK4_MLTS_Lshape_PellSparse_conv_test(argc, argv);        // MehlinExact's math, sparse-matvec implementation (exact + faster)
   // ERK4_MLTS_Lshape_HybridTerminal_conv_test(argc, argv);    // global non-terminal + restricted terminal only
   // ERK4_MLTS_Lshape_FivePoint_conv_test(argc, argv);         // superseded 5-point mesh family attempt
   // ERK4_MLTS_Lshape_DiagN2_test(argc, argv);
   //
   // Square (mlts_2026/square/), smooth sinusoidal solution, no singularity:
   // ERK4_MLTS_Square_conv_test(argc, argv);                       // hybrid architecture, broken at deep L
   // ERK4_MLTS_Square_Pell_conv_test(argc, argv);                  // truncated-submatrix P_ell (broken, same as Lshape RestrictedLeveled)
   // ERK4_MLTS_Square_PellSparse_conv_test(argc, argv);            // exact sparse P_ell, CENTER-graded (anomalous with abrupt grading)
   // ERK4_MLTS_CornerSquare_PellSparse_conv_test(argc, argv);      // exact sparse P_ell, CORNER-graded, abrupt corner_boost grading
   // ERK4_MLTS_CornerSquareFixed_PellSparse_conv_test(argc, argv); // CORNER-graded, fixed local depth (p_global constant)
   // ERK4_MLTS_CornerSquareRamp_PellSparse_conv_test(argc, argv);  // CORNER-graded, local depth ramping +1 per N
   // ERK4_MLTS_CornerSquareCoarse_PellSparse_conv_test(argc, argv);// CORNER-graded, coarse start (n0=2) + progressive depth -- clean order 4
   // ERK4_MLTS_CenterSquareCoarse_PellSparse_conv_test(argc, argv);// CENTER-graded, same progressive construction -- ALSO clean order 4; also the source of the "p=2^N growing" curve in PLOT ARTICLES above
   // ERK4_MLTS_CenterSquare_TwoLevel_conv_test / ERK4_MLTS_CenterSquare_FixedLevel_conv_test -- moved to "PLOT ARTICLES" above
   //
   // ERK4_MLTS_RestrictedGD_test(argc, argv);
   // ERK4_MLTS_L5small_v2design_validation_test(argc, argv);
   // ERK4_MLTS_Lsweep_small_test(argc, argv);
   // ERK4_MLTS_GlobalGD_test(argc, argv);
   // ERK4_MLTS_RestrictedGD_test(argc, argv);
   // ERK4_MLTS_DiagCompare_test(argc, argv);
   // ERK_LTS_stab(argc, argv);
   // ERK4_LTS_SSTAB(argc, argv);
   // ERK4_LTS_optimised(argc, argv);

   // HeterogeneousEHHOFirstOrder(argc, argv); 
   // HeterogeneousERK4_LTS_HHO_FirstOrder(argc, argv);
   // HeterogeneousERK4_LTS_HHO_FirstOrder_stab(argc, argv);

   // BrainIHHOFirstOrder(argc, argv);
   // BrainEHHOFirstOrder(argc, argv);

}




