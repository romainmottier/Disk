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
 * @brief LineSearch
 *
 */

// Référence of Aitken and Anderson algorithm.
// Isabelle Ramière, Thomas Helfer. Iterative residual-based vector methods to accelerate fixed
// point iterations.
// Computers & Mathematics with Applications, 2015, 70, pp.2210 - 2226.
// 10.1016/j.camwa.2015.08.025. cea-01403292

template < typename T >
class ConvergenceAcceleration {

    typedef dynamic_vector< T > vector_type;

    int n_iter;

    vector_type Xa_km, Xa_k;
    vector_type X_k, X_km;

  public:
    ConvergenceAcceleration() : n_iter( 0 ) {}

    vector_type aitken( const vector_type &X_kp ) {
        // Also called crossed secand method - eq 44.
        if ( n_iter == 0 ) {
            n_iter++;
            X_km = X_kp;
            Xa_km = X_km;
            return Xa_km;
        } else if ( n_iter == 1 ) {
            n_iter++;
            X_k = X_kp;
            Xa_k = X_k;
            return Xa_k;
        } else {
            n_iter++;
            const auto G_k = X_kp;
            const auto G_km = X_k;
            const auto dG_k = G_k - G_km;

            const auto dX_k = G_k - Xa_k;
            const auto dX_km = G_km - Xa_km;
            const auto ddX = dX_k - dX_km;

            // compute acceleration
            const T wr = ddX.dot( dG_k ) / ddX.squaredNorm();

            // compute accelerted solution
            const vector_type Xa_kp = G_k - wr * dX_k;

            // update
            X_km = X_k;
            X_k = X_kp;

            Xa_km = Xa_k;
            Xa_k = Xa_kp;

            return Xa_kp;
        }
    }

    vector_type relaxation( const vector_type &X_kp, const T omega = 0.5 ) {
        if ( n_iter == 0 ) {
            n_iter++;
            Xa_k = X_kp;
            return Xa_k;
        } else {
            n_iter++;

            // compute accelerted solution
            const vector_type Xa_kp = ( 1.0 - omega ) * Xa_k + omega * X_kp;

            // update
            Xa_k = Xa_kp;

            return Xa_kp;
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

    vector_type anderson( const vector_type &X_kp ) {
        // also called alternate decant method - eq.45
        if ( n_iter == 0 ) {
            n_iter++;
            X_km = X_kp;
            Xa_km = X_km;
            return Xa_km;
        } else if ( n_iter == 1 ) {
            n_iter++;
            X_k = X_kp;
            Xa_k = X_k;
            return Xa_k;
        } else {
            n_iter++;
            const auto G_k = X_kp;
            const auto G_km = X_k;

            const auto dX_k = G_k - Xa_k;
            const auto dX_km = G_km - Xa_km;
            const auto ddX = dX_k - dX_km;

            // compute acceleration
            const T wr = ddX.dot( dX_k ) / ddX.squaredNorm();

            // compute accelerted solution
            const vector_type Xa_kp = ( 1.0 - wr ) * G_k + wr * G_km;

            // update
            X_km = X_k;
            X_k = X_kp;

            Xa_km = Xa_k;
            Xa_k = Xa_kp;

            return Xa_kp;
        }
    }
};
} // namespace mechanics
} // namespace disk