/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <vector>

namespace disk {

namespace mechanics {

/**
 * @brief Aitken acceleration
 *
 */
template < typename T >
class ConvergenceAcceleration {

    typedef dynamic_vector< T > vector_type;

    int n_iter;

    vector_type va_km, va_k;
    vector_type v_k, v_km;

  public:
    ConvergenceAcceleration() : n_iter( 0 ) {}

    vector_type aitken( const vector_type &v_kp ) {
        if ( n_iter == 0 ) {
            n_iter++;
            v_km = v_kp;
            return v_km;
        } else if ( n_iter == 1 ) {
            n_iter++;
            v_k = v_kp;
            return v_k;
        } else {
            n_iter++;
            const vector_type vt = v_kp - 2 * v_k + v_km;
            const vector_type dv = v_kp - v_k;

            // compute relaxation
            const T wr = dv.dot( vt ) / vt.squaredNorm();

            // compute accelerted solution
            const vector_type va_kp = wr * v_k + ( 1 - wr ) * v_kp;

            // update
            v_km = v_k;
            v_k = v_kp;

            return va_kp;
        }
    }

    vector_type relaxation( const vector_type &v_kp, const T omega = 0.5 ) {
        if ( n_iter == 0 ) {
            n_iter++;
            va_k = v_kp;
            return va_k;
        } else {
            n_iter++;

            // compute accelerted solution
            const vector_type va_kp = ( 1.0 - omega ) * va_k + omega * v_kp;

            // update
            va_k = va_kp;

            return va_kp;
        }
    }

    template < typename Func >
    void secant( const Func &func, const double ALF = 0.1, const int MAXIT = 10 ) {

        const T TOLX = std::numeric_limits< T >::epsilon();
        T p0, p1, f1, f0, rho_0, rho_1;
        T rho, rho_neg, rho_pos, rho_opt, rho_new, rho_cur;
        T f, f_opt, f_cur;
        bool b_pos;

        // fixed paramters - from code_aster
        const T rho_min = 1e-2, rho_max = 10., rho_excl = 0.9e-2;
        const T parmul = 3.0;

        // Doc:
        // https://codeaster.pages.pleiade.edf.fr/doc/docaster/manuals/man_r/r5/r5.03.01/Recherche_lin_aire.html
        // METHODE="SECANT"

        // Compute residual.dot(increment) (and update solution)
        // auto _f = [&fvec, &func, &dx, &xold, &n, &x]( const T &rho ) {
        //     for ( int j = 0; j < n; j++ )
        //         x[j] = xold[j] + rho * dx[j];
        //     const auto norm = func( x );

        //     T f = 0.0;
        //     for ( int j = 0; j < n; j++ )
        //         f += fvec[j] * dx[j];
        //     return f;
        // };

        // project bound on admissible interval
        auto _proj = [rho_min, rho_max, rho_excl]( T &rho ) {
            const T rho_tmp = rho;
            if ( rho_tmp < rho_min ) {
                rho = rho_min;
            }
            if ( rho_tmp > rho_max ) {
                rho = rho_max;
            }
            if ( rho_tmp < 0.0 && rho_tmp >= -rho_excl ) {
                rho = -rho_excl;
            }
            if ( rho_tmp >= 0 && rho_tmp <= rho_excl ) {
                rho = rho_excl;
            }
        };

        // initial values
        const T f_old = func( 0.0 );
        const T f_cvg = ALF * std::abs( f_old );
        const T sens = ( f_old <= 0.0 ) ? 1.0 : -1.0;

        rho_opt = 1.0, rho = sens * 1.0;
        rho_neg = 0.0, rho_pos = std::numeric_limits< T >::signaling_NaN();
        f_opt = 10e100;

        rho_0 = 0.0, rho_1 = rho_0;
        f0 = sens * f_old, f1 = f0;

        b_pos = false;

        for ( int its = 0; its < MAXIT; its++ ) {
            // Compute new residual
            try {
                f = func( rho );
            } catch ( ... ) {
                break;
            }

            rho_cur = sens * rho;
            f_cur = sens * f;

            // Store value
            rho_0 = rho_1, f0 = f1;
            rho_1 = rho_cur, f1 = f_cur;

            // Update bounds
            if ( f_cur < 0.0 ) {
                rho_neg = rho_cur;
            } else {
                b_pos = true;
                rho_pos = rho_cur;
            }

            // Optimal solution until now ?
            if ( std::abs( f_cur ) < std::abs( f_opt ) ) {
                rho_opt = rho_cur;
                f_opt = f_cur;
                _proj( rho_opt );
            }

            // Search maximal bound
            if ( b_pos ) {
                if ( std::abs( f1 ) >= std::abs( f0 ) ) {
                    // f is not decreased - use dichotomie
                    rho_new = 0.5 * ( rho_neg + rho_pos );
                } else {
                    // linear interpolation
                    if ( std::abs( rho_1 - rho_0 ) > TOLX ) {
                        p1 = ( f1 - f0 ) / ( rho_1 - rho_0 );
                        p0 = f0 - p1 * rho_0;

                        if ( std::abs( p1 ) <= std::abs( f0 ) / ( rho_pos + rho_0 ) ) {
                            rho_new = 0.5 * ( rho_neg + rho_pos );
                        } else {
                            rho_new = -p0 / p1;
                        }
                    } else {
                        // failed
                        break;
                    }
                }
            } else {
                rho_new = parmul * rho_cur;
            }

            // minimal bound
            if ( rho_new < rho_neg ) {
                if ( b_pos ) {
                    rho_new = 0.5 * ( rho_neg + rho_pos );
                } else {
                    // failed
                    break;
                }
            }

            // maximal bound
            if ( b_pos && rho_new > rho_pos ) {
                rho_new = 0.5 * ( rho_neg + rho_pos );
            }

            // project bound
            _proj( rho_new );

            // update
            rho = sens * rho_new;

            // Test convergence ?
            if ( std::abs( f_opt ) <= f_cvg ) {
                break;
            }
        }

        /* Return optimal value */
        // std::cout << "rho_opt: " << rho_opt << std::endl;
        f = func( rho_opt, false );
    }

    void aitken2() {}
};
} // namespace mechanics
} // namespace disk