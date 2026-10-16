//
//  manufactured_solutions.hpp
//
//  Manufactured solution (scal_vec_analytic_functions) of the name problem.solution of the Lua parameters.
//

#pragma once
#ifndef physics_manufactured_solutions_hpp
#define physics_manufactured_solutions_hpp

#include "../mesh/domain.hpp"

namespace coupling_common {

    // Manufactured solution of problem.solution (wave_parameters::solution_names)
    inline scal_vec_analytic_functions::EFunctionType solution_type(const std::string& name) {
        typedef scal_vec_analytic_functions::EFunctionType F;
        static const std::map<std::string, F> types = {
            {"interface_coupling",          F::EFunctionInterfaceCoupling},
            {"non_polynomial",              F::EFunctionNonPolynomial},
            {"acoustic_quadratic_in_time",  F::reproduction_acoustic},
            {"acoustic_quadratic_in_space", F::EFunctionAcousticQuadraticInSpace},
            {"acoustic_non_polynomial",     F::EFunctionAcousticNonPolynomial},
            {"acoustic_quadratic_in_space_cos", F::EFunctionAcousticQuadraticInSpaceCos},
            {"acoustic_non_polynomial_cos", F::EFunctionAcousticNonPolynomialCos},
            {"lshape_corner",               F::EFunctionAcousticLShapeCorner},
            {"lshape_corner_sin",           F::EFunctionAcousticLShapeCornerSin},
            {"elastic_quadratic_in_time",   F::reproduction_elastic},
            {"elastic_quadratic_in_space",  F::EFunctionElasticQuadraticInSpace},
            {"elastic_non_polynomial",      F::EFunctionElasticNonPolynomial},
            {"coupled_sinusoidal",          F::EFunctionNonPolynomial_paper},
            {"coupled_poly_in_space",       F::EFunctionQuadraticInSpace},
            {"coupled_poly_in_time",        F::EFunctionQuadraticInTime},
            {"ricker",                      F::EFunctionNonPolynomial}};       // data replaced by zero (HHOFirstOrder)
        return types.at(name);
    }

}

#endif /* physics_manufactured_solutions_hpp */
