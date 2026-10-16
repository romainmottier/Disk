//
//  studies.hpp
//
//  First-order formulation (HHOFirstOrder.hpp): integration with the outputs of a run (energy history, sensors,
//  snapshots), CFL study, spectral study, space and time convergence studies.
//

#pragma once
#ifndef first_order_studies_hpp
#define first_order_studies_hpp

#include "time_schemes.hpp"

namespace hho_first_order {

    // Base name of the outputs: the data table without its extension
    inline std::string output_base(const wave_parameters& prm) {
        std::string f = prm.table_file();
        return f.substr(0, f.size() - 4);
    }

    // Integration with the outputs of the run: history of the discrete energy 1/2 x^T M x (file prm.energy_file():
    // "# t E" header with max |E^n - E^0| / E^0, one line per step; acoustic and elastic parts E_acoustic E_elastic in
    // two more columns for a coupled domain), sensors and snapshots
    template<typename Problem>
    inline vector_type integrate_observed(Problem& pb, const wave_parameters& prm, vector_type x, size_t nt, size_t level = 0) {
        const RealType dt = prm.time.final_time / RealType(nt);
        const std::string base = output_base(prm);
        const bool energy_Q = prm.output.energy;
        std::ostringstream hist;
        hist << std::setprecision(16) << std::scientific;
        const RealType E0 = pb.energy(x);
        RealType var = 0.0;
        // acoustic and elastic parts when both media are present
        bool parts_Q = false;
        if constexpr (std::is_same_v<Problem, problem>)
            parts_Q = !pb.dd.e_material.empty() && !pb.dd.a_material.empty();
        auto energy_line = [&](RealType t, const vector_type& y) {
            std::array<RealType, 2> parts{};
            bool parts_ready = false;
            if constexpr (std::is_same_v<Problem, problem>) {
                if (parts_Q && pb.energy_from_parts()) {
                    parts = pb.energy_parts(y);
                    parts_ready = true;
                }
            }
            RealType E = parts_ready ? parts[0] + parts[1] : pb.energy(y);
            var = std::max(var, std::abs(E - E0) / E0);
            hist << t << " " << E;
            if constexpr (std::is_same_v<Problem, problem>) {
                if (parts_Q) {
                    if (!parts_ready) parts = pb.energy_parts(y);
                    hist << " " << parts[0] << " " << parts[1];
                }
            }
            hist << "\n";
        };
        std::unique_ptr<sensors> sens;
        std::unique_ptr<snapshots> snaps;
        std::unique_ptr<peak_map> peaks;
        if constexpr (std::is_same_v<Problem, problem>) {
            if (!prm.output.sensors.empty())   sens  = std::make_unique<sensors>(pb, prm, base);
            if (!prm.output.snapshots.empty()) snaps = std::make_unique<snapshots>(pb, prm, base, dt);
            if (prm.output.peak_samples > 0)   peaks = std::make_unique<peak_map>(pb, prm, base, x);
            if (sens)  sens->record(pb, 0, 0.0, x);
            if (snaps) snaps->maybe_write(pb, 0, 0.0, x);
        }
        if (energy_Q) energy_line(0.0, x);
        size_t step = 0;
        x = integrate(pb, prm, x, nt, [&](RealType t, const vector_type& y) {
            step++;
            if (energy_Q) energy_line(t, y);
            if constexpr (std::is_same_v<Problem, problem>) {
                if (sens)  sens->record(pb, step, t, y);
                if (snaps) snaps->maybe_write(pb, step, t, y);
                if (peaks) peaks->record(step, y);
            }
        }, level);
        if constexpr (std::is_same_v<Problem, problem>) {
            if (peaks) peaks->write();
        }
        if (energy_Q) {
            std::ofstream os(prm.energy_file());
            os << "# t E" << (parts_Q ? " E_acoustic E_elastic" : "")
               << " (discrete energy 1/2 x^T M x), max |E - E0| / E0 = " << var << std::endl << hist.str();
            std::cout << bold << cyan << "      energy: max |E^n - E^0| / E^0 = " << var << reset << std::endl;
        }
        return x;
    }

    // CFL study (Mottier-Ern-Guillot, Algorithm 1): ERK(s) with the initial data of the solution, homogeneous Dirichlet
    // data and no source; starting from time.n_steps (doubled until stable, then decreased by 10 % while stable), the
    // number of steps N is decreased by cfl_decrement N + 1 after each stable run, until a run where the discrete
    // energy grows by more than
    // cfl_threshold of the initial or of the previous energy (or is not finite); cfl_absolute: variation in absolute
    // value instead of growth. Table: last stable and first unstable N, CFL* = c_max dt / h with dt = T / (last stable
    // N), h = h_min and h_max (cell diameters).
    inline void cfl_study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        problem pb(prm, prm.mesh.level, true);
        const vector_type x0 = prm.study.cfl_random ? pb.random_state() : pb.initial_state();
        const RealType E0 = pb.energy(x0), tf = prm.time.final_time, thr = prm.study.cfl_threshold;
        const size_t s = prm.stages();
        Matrix<RealType, Dynamic, Dynamic> a; vector_type b, c;
        erk_butcher_tableau::erk_tables(int(s), a, b, c);
        const vector_type F = pb.eval_F(0.0);
        auto stable_run = [&](size_t N) {
            vector_type x = x0;
            RealType E_prev = E0;
            bool ok = true;
            // after each step: growth of the energy within the threshold, else the run stops
            auto check = [&](RealType, const vector_type& y) {
                RealType E = pb.energy(y);
                RealType d_prev = (E - E_prev) / E_prev, d_0 = (E - E0) / E0;
                if (prm.study.cfl_absolute) { d_prev = std::abs(d_prev); d_0 = std::abs(d_0); }
                ok = finite_value(E) && finite_value(d_prev) && finite_value(d_0) && d_prev <= thr && d_0 <= thr;
                E_prev = E;
                return ok;
            };
            if (prm.time.scheme != "erk")       // lts, li, mlts: the whole run, energy checked after each step
                integrate(pb, prm, x, N, [&](RealType t, const vector_type& y) { if (ok) check(t, y); });
            else if (prm.solver.matrix_free)
                pb.mf().erk(x, N, tf, a, b, c, [&](RealType t) { return pb.eval_F(t); }, zero_force(pb), check);
            else
                runge_kutta::erk(x, N, tf, a, b, c, [&](RealType, vector_type& y, vector_type& k) {
                    vector_type Fi = F;
                    pb.erk().SetFg(Fi);
                    pb.erk().erk_weight(y, k);
                }, check);
            return ok;
        };
        size_t N = prm.time.n_steps, stable = 0, unstable = 0;
        for (size_t tries = 0; tries < 12 && !stable_run(N); tries++)
            N *= 2;
        // coarse approach of the limit (steps of 10 %), then the decrements of Algorithm 1 from the last stable N
        for (size_t M = N - (N / 10 + 1); M > 0 && stable_run(M); M -= M / 10 + 1)
            N = M;
        for (size_t run = 0; run < prm.study.cfl_runs && N > 0; run++) {
            bool ok = stable_run(N);
            std::cout << bold << cyan << "      " << prm.scheme_label() << " N = " << N << ": " << (ok ? "stable" : "unstable") << reset << std::endl;
            if (!ok) { unstable = N; break; }
            stable = N;
            N -= size_t(RealType(N) * prm.study.cfl_decrement) + 1;
        }
        if (stable == 0)
            throw std::runtime_error("cfl study: no stable number of steps found from time.n_steps");
        const RealType dt = tf / RealType(stable);
        names = {"stable_steps", "unstable_steps", "cfl_h_min", "cfl_h_max"};
        rows.push_back({pb.dd.h_max, dt, {RealType(stable), RealType(unstable), pb.dd.c_max * dt / pb.dd.h_min,
                                          pb.dd.c_max * dt / pb.dd.h_max}, pb.dd.h_min});
        std::cout << bold << cyan << "      c_max = " << pb.dd.c_max << "  h_min = " << pb.dd.h_min << "  h_max = " << pb.dd.h_max << reset << std::endl;
        print_table(names, rows, false, std::cout);
        print_table(names, rows, false, log);
    }

    // Spectral radius of one step (study.type = "rho"): one run of time.n_steps steps of the scheme with random initial
    // cell data (all the modes, consistent faces) and no source; the energy grows or decays asymptotically like
    // rho(C)^(2n), C the amplification operator of a step, and rho is estimated from the second half of the run,
    // rho = (E_N / E_{N/2})^(1/N) (from the last finite energies if the run overflows). Table: rho and the steps used.
    inline void rho_study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        problem pb(prm, prm.mesh.level, true);
        vector_type x = pb.random_state();
        const size_t N = prm.time.n_steps;
        const RealType dt = prm.time.final_time / RealType(N);
        std::vector<RealType> logE(1, std::log(pb.energy(x)));
        bool finite = true;
        integrate(pb, prm, x, N, [&](RealType, const vector_type& y) {
            if (!finite) return;
            const RealType E = pb.energy(y);
            if (!finite_value(E) || E <= 0) { finite = false; return; }
            logE.push_back(std::log(E));
        });
        const size_t m = logE.size() - 1;
        RealType rho = 0;
        if (m >= 2) rho = std::exp((logE[m] - logE[m / 2]) / (2 * RealType(m - m / 2)));
        names = {"rho", "steps"};
        rows.push_back({pb.dd.h_max, dt, {rho, RealType(m)}, pb.dd.h_min});
        std::cout << bold << cyan << "      rho = " << rho << " (" << m << " steps)" << reset << std::endl;
        print_table(names, rows, false, std::cout);
        print_table(names, rows, false, log);
    }

    // Spectral study (Mottier-Ern-Khot-Guillot, Section 6.1): largest eigenvalue gamma of
    // K_schur^T M_T^-1 K_schur X = gamma M_T X, K_schur = K_TT - K_TF K_FF^-1 K_FT the Schur complement of the face block
    // of the stiffness (coupling included), M_T the cell mass matrix; spectral radius rho = sqrt(gamma), normalized by
    // sqrt(number of cells). Dense computation, on each mesh of the study.
    inline void spectrum_study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        typedef Matrix<RealType, Dynamic, Dynamic> dense_type;
        names = {"rho", "rho_normalized", "cells"};
        for (size_t l : study_levels(prm)) {
            problem pb(prm, l, true);
            auto& as = pb.assembler();
            const long n_c = long(as.get_e_n_cells_dof() + as.get_a_n_cells_dof()), n_f = long(as.get_n_face_dof());
            dense_type K = dense_type(as.LHS), M = dense_type(as.MASS);
            dense_type S = K.topLeftCorner(n_c, n_c);
            if (n_f > 0)
                S -= K.topRightCorner(n_c, n_f) * K.bottomRightCorner(n_f, n_f).partialPivLu().solve(dense_type(K.bottomLeftCorner(n_f, n_c)));
            dense_type Mc = M.topLeftCorner(n_c, n_c);
            dense_type A = S.transpose() * Mc.llt().solve(S);
            A = 0.5 * (A + A.transpose()).eval();
            Eigen::GeneralizedSelfAdjointEigenSolver<dense_type> es(A, Mc, Eigen::EigenvaluesOnly);
            const RealType gamma = es.eigenvalues().maxCoeff();
            const RealType cells = RealType(pb.msh.cells_size());
            rows.push_back({pb.dd.h_max, 0.0, {std::sqrt(gamma), std::sqrt(gamma) / std::sqrt(cells), cells}, pb.dd.h_min});
            std::cout << bold << cyan << "      l = " << l << "  cells = " << cells << "  cell dofs = " << n_c << "  face dofs = " << n_f << reset << std::endl;
            print_table(names, rows, false, std::cout);
        }
        print_table(names, rows, false, log);
    }


    // Space (or single run) and time studies of the first-order formulation for a problem type (problem: coupled
    // assembler; problem_three_fields: three-field formulation of the elastic medium): convergence rows and names
    template<typename Problem>
    inline void study(const wave_parameters& prm, std::ofstream& log, std::vector<conv_row>& rows, std::vector<std::string>& names) {
        error_columns all = columns(prm, {"elastic_v", "elastic_sigma"}, {"acoustic_p", "acoustic_q"});
        auto all_errors = [](const std::array<RealType, 5>& e) { return std::vector<RealType>{e[1], e[2], e[3], e[4]}; };
        const RealType tf = prm.time.final_time;
        // no error for the Ricker wavelet; energy history, sensors and snapshots of the last (or only) mesh
        const bool ricker = prm.problem.solution == "ricker";
        const bool observed = prm.output.energy || !prm.output.sensors.empty() || !prm.output.snapshots.empty()
                              || prm.output.peak_samples > 0;
        // components of the errors in the energy norms (output.errors = "energy", problem::energy_error_names)
        const bool energy_norms = prm.output.errors == "energy";
        auto final_errors = [&](Problem& pb, vector_type& x) -> std::vector<RealType> {
            if constexpr (std::is_same_v<Problem, problem>) {
                if (energy_norms)
                    return pb.energy_errors(x, tf);
            }
            return all.select(all_errors(pb.errors(x, tf)));
        };
        if (prm.study.type != "time") {
            names = ricker ? std::vector<std::string>{} : energy_norms ? problem::energy_error_names() : all.names;
            if (prm.output.cpu) { names.push_back("cpu"); names.push_back("dofs"); names.push_back("face_dofs"); }
            std::vector<size_t> levels = study_levels(prm);
            for (size_t l : levels) {
                auto t0 = std::chrono::steady_clock::now();
                Problem pb(prm, l);
                size_t nt = study_steps(prm, l, pb.dd.h_min);
                // output.profile: wall time of the time integration alone (operator built before)
                if constexpr (std::is_same_v<Problem, problem>) {
                    if (prm.output.profile && pb.m_matrix_free) pb.mf();
                }
                vector_type x0 = pb.initial_state();
                auto t1 = std::chrono::steady_clock::now();
                vector_type x = observed && l == levels.back() ? integrate_observed(pb, prm, std::move(x0), nt, l)
                                                               : integrate(pb, prm, std::move(x0), nt, {}, l);
                if (prm.output.profile)
                    std::cout << "      integration wall time (s): "
                              << std::chrono::duration<RealType>(std::chrono::steady_clock::now() - t1).count() << std::endl;
                RealType cpu = std::chrono::duration<RealType>(std::chrono::steady_clock::now() - t0).count();
                std::vector<RealType> e = ricker ? std::vector<RealType>{} : final_errors(pb, x);
                if (prm.output.cpu) {
                    e.push_back(cpu);
                    e.push_back(RealType(pb.assembler().LHS.rows()));     // unknowns before / after the static condensation
                    e.push_back(RealType(pb.n_face_dof()));               //   of the cell unknowns
                }
                rows.push_back({pb.dd.h_max, tf / nt, e, pb.dd.h_min});
                std::cout << bold << cyan << "      l = " << l << "  h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min
                          << "  nt = " << nt << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << reset << std::endl;
                print_table(names, rows, false, std::cout);
            }
            print_table(names, rows, false, log);
        }
        else {
            // against the exact solution only, or against a reference solution of the same integrator (mass norm of the
            // difference) with the exact errors of the velocity / scalar unknowns
            bool reference = prm.study.reference == "numerical";
            error_columns cols = reference ? columns(prm, {"elastic_v"}, {"acoustic_p"}) : all;
            names = energy_norms ? problem::energy_error_names() : cols.names;
            if (reference) names.insert(names.begin(), "mass_norm");
            Problem pb(prm, prm.mesh.level);
            size_t n0 = prm.time.n_steps;
            size_t n_levels = prm.study.time_levels;
            vector_type x_ref;
            if (reference) {
                size_t nt_ref = prm.study.reference_factor * n0 * (size_t(1) << (n_levels - 1));
                x_ref = integrate(pb, prm, pb.initial_state(), nt_ref);
                std::cout << bold << cyan << "      reference nt = " << nt_ref << reset << std::endl;
            }
            std::cout << bold << cyan << "      h_max = " << pb.dd.h_max << "  h_min = " << pb.dd.h_min
                      << "  fine substeps = " << pb.p << "  fine cells = " << pb.n_fine_cells << "/" << pb.msh.cells_size() << reset << std::endl;
            for (size_t i = 0; i < n_levels; i++) {
                size_t nt = n0 * (size_t(1) << i);
                vector_type x = integrate(pb, prm, pb.initial_state(), nt);
                auto e = pb.errors(x, tf);
                std::vector<RealType> errors = energy_norms ? final_errors(pb, x) : reference ? cols.select({e[1], e[3]}) : cols.select(all_errors(e));
                if (reference) errors.insert(errors.begin(), mass_norm(pb.assembler().MASS, x - x_ref));
                rows.push_back({pb.dd.h_max, tf / nt, errors, pb.dd.h_min});
                print_table(names, rows, true, std::cout);
            }
            print_table(names, rows, true, log);
        }
    }

}

#endif /* first_order_studies_hpp */
