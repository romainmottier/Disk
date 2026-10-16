//
//  HHOSecondOrder.hpp
//
//  Elasto-acoustic coupling, second-order formulation (elastoacoustic_two_fields_assembler: displacement in the
//  solid, velocity potential in the fluid), M u'' + C u' + K u = F(t). Every setting comes from the Lua input file
//  (parameters/lua_parameters.hpp, template share/second_order.lua). The Runge-Kutta integrators solve a first-order
//  system: the auxiliary variable is the cell velocity v = u_T',
//      u_T' = v,     M v' = F_T - K_TT u_T - K_TF u_F,     C_FF u_F' = F_F - K_FT u_T - K_FF u_F,
//  M being the cell mass. The coupling C = [[0, B], [-B^T, 0]] only acts on the interface faces: a QR factorization of
//  B on each interface face splits the elastic face unknowns into a normal part (range of B) and a tangential part.
//  The normal part and the acoustic interface unknowns, d, are dynamic (C_dd invertible); the other face unknowns, a
//  (all the other faces, tangential part on the interface), are algebraic, K_aa a = F_a - K_aT u_T - K_ad d.
//  State y = (u_T, v, d), y' = f(t, y) = B y + g(t).
//    space : face degree k, equal (k,k) or mixed (k+1,k) order, least-squares (plain least squares / Lehrenfeld-
//            Schoberl) or HHO stabilization, weight in 1/h or O(1), weight gamma (number, or "auto" = weight_factor
//            gamma*, gamma* the critical weight of the splitting computed cell by cell, Steins et al. Lemma 4.3)
//    faces : explicit integrators (erk, lts) and li: algebraic faces from the splitting of Steins et al. (M2AN 57,
//            2023), D a^{m+1} = r - (K_aa - D) a^m with D = gamma sigma*, the block-diagonal face part of the
//            stabilization (fully explicit, contraction gamma*/gamma, relative increment faces.tolerance), or from a
//            sparse factorization of K_aa ("direct"); sdirk solves them within its stage systems
//    time  : ERK(s), SDIRK(s,s+1), two-level LTS-RK4, LI-SDIRK(4,3), final time, number of time steps
//    mesh, problem (manufactured solution), study (one run, space-time convergence, time convergence), output
//  The SDIRK stages and the LI fine solves, (I - a_ii dt B P) Y = R, are solved with a, u_T eliminated, through the
//  sparse system (K + M / (a_ii dt)^2 + C / (a_ii dt)) restricted by the projector P (P = I for SDIRK). The explicit
//  integrators print omega_max (power iteration) and warn when dt is above their stability limit.
//  Steins' optimal explicit setting: mixed order, least squares, 1/h, weight "auto" with weight_factor in [3, 5],
//  splitting; implicit: mixed order, least squares (Lehrenfeld-Schoberl), 1/h, weight 1. With equal order and the
//  HHO stabilization, alpha = lambda_max(S_aa, sigma*_aa) - 1 >= 1 for k = 1, 3 in the solid: the splitting cannot
//  converge (use faces.solver = "direct").
//

#pragma once
#ifndef HHOSecondOrder_hpp
#define HHOSecondOrder_hpp

// system.hpp (first-order system), problem.hpp (problem), time_schemes.hpp (integrators), studies.hpp (summaries)
#include "studies.hpp"

void HHOSecondOrder(const wave_parameters& prm);

void HHOSecondOrder(const wave_parameters& prm) {

    using namespace hho_second_order;

    // algebraic faces of the explicit (erk, lts) and locally implicit (li) integrators: splitting or direct; sdirk
    // solves them within its stage systems
    bool splitting = prm.explicit_scheme() && prm.faces.solver == "splitting";
    std::cout << std::endl << bold << red << "   SECOND-ORDER ELASTO-ACOUSTIC COUPLING, " << prm.scheme_label()
              << ", " << prm.space_label() << (splitting ? ", faces: Steins splitting" : ", faces: direct") << reset << std::endl;
    std::string fname = prm.output_file();
    std::ofstream log(fname);
    log << "second-order formulation" << std::endl;
    prm.print(log);

    // stabilization weight: space.weight, or weight_factor gamma* (gamma* computed cell by cell on the first mesh of
    // the study, it only depends on the cell shapes)
    RealType gamma = prm.space.weight;
    if (prm.space.weight_auto) {
        problem probe(prm, prm.mesh.level, false, 1.0);
        RealType gamma_star = probe.critical_weight();
        gamma = prm.space.weight_factor * gamma_star;
        std::ostringstream os;
        os << "      gamma* = " << gamma_star << " (cells of the mesh l = " << prm.mesh.level << "), gamma = "
           << prm.space.weight_factor << " gamma* = " << gamma << std::endl;
        std::cout << os.str();
        log << os.str();
    }

    // errors at the final time: elastic displacement, stress (reconstruction) and velocity, acoustic potential, its
    // reconstructed gradient and its time derivative, L2 norms
    error_columns all = columns(prm, {"elastic_u", "elastic_sigma", "elastic_v"}, {"acoustic_phi", "acoustic_grad_phi", "acoustic_phi_t"});
    auto all_errors = [](const std::array<RealType, 7>& e) { return std::vector<RealType>{e[1], e[2], e[5], e[3], e[4], e[6]}; };
    const bool explicit_scheme = prm.explicit_scheme();

    // Newmark energy history of the last mesh of a space study / of the finest run of a time study
    std::ofstream energy_os;
    if (prm.output.energy && prm.time.scheme == "newmark")
        energy_os.open(prm.energy_file());
    auto run = [&](problem& pb, size_t nt, bool record) {
        if (prm.time.scheme == "newmark" && record && energy_os.is_open()) {
            energy_os.seekp(0);
            std::ostringstream hist;
            vector_type y = integrate_newmark(pb, prm, nt, &hist);
            // relative energy variation max_n |E^n - E^0| / E^0
            std::istringstream is(hist.str());
            RealType t, E, E0 = -1.0, var = 0.0;
            while (is >> t >> E) {
                if (E0 < 0.0) E0 = E;
                var = std::max(var, std::abs(E - E0) / E0);
            }
            energy_os << "# t E (discrete energy (69) of JSC21, Dirichlet faces = data), max |E - E0| / E0 = " << var << std::endl << hist.str();
            std::cout << bold << cyan << "      Newmark energy: max |E^n - E^0| / E^0 = " << var << reset << std::endl;
            return y;
        }
        return integrate(pb, prm, pb.initial_state(), nt);
    };

    const RealType tf = prm.time.final_time;
    std::vector<conv_row> rows;
    std::vector<std::string> names;
    if (prm.study.type == "cfl") {
        // CFL study (as in the first-order formulation, first_order/studies.hpp, but the energy is checked at the final
        // time only): initial data of the solution (problem.source = "zero", problem.dirichlet = "zero" for a pure
        // stability test), energy E = v^T M v + x^T K x of the state (energy_distance to 0); starting from
        // time.n_steps (doubled until stable, then decreased by 10 % while stable), N decreased by cfl_decrement N + 1
        // after each stable run until E(t_f) > (1 + cfl_threshold) E(0) or is not finite. Table: last stable and first
        // unstable N, CFL* = c_max dt / h (h = h_min, h_max).
        problem pb(prm, prm.mesh.level, splitting, gamma);
        std::string summary = problem_summary(pb, explicit_scheme);
        std::cout << summary;
        log << summary;
        const vector_type y0 = prm.study.cfl_random ? pb.random_state() : pb.initial_state(), z = vector_type::Zero(y0.size());
        auto energy = [&](const vector_type& y) { RealType e = pb.energy_distance(y, z, tf); return e * e; };
        const RealType E0 = energy(y0), thr = prm.study.cfl_threshold;
        auto stable_run = [&](size_t N) {
            RealType E = energy(integrate(pb, prm, y0, N));
            RealType d = (E - E0) / E0;
            if (prm.study.cfl_absolute) d = std::abs(d);
            return finite_value(E) && finite_value(d) && d <= thr;
        };
        size_t N = prm.time.n_steps, stable = 0, unstable = 0;
        for (size_t tries = 0; tries < 12 && !stable_run(N); tries++)
            N *= 2;
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
        RealType c_max = 0.0;
        for (auto& m : pb.dd.a_material) { auto a = m.second; c_max = std::max(c_max, a.vp()); }
        for (auto& m : pb.dd.e_material) { auto e = m.second; c_max = std::max(c_max, e.vp()); }
        const RealType dt = tf / RealType(stable);
        names = {"stable_steps", "unstable_steps", "cfl_h_min", "cfl_h_max"};
        rows.push_back({pb.dd.h_max, dt, {RealType(stable), RealType(unstable), c_max * dt / pb.dd.h_min, c_max * dt / pb.dd.h_max},
                        pb.dd.h_min});
        print_table(names, rows, false, std::cout);
        print_table(names, rows, false, log);
    }
    else if (prm.study.type != "time") {
        names = all.names;
        std::vector<size_t> levels = study_levels(prm);
        for (size_t l : levels) {
            problem pb(prm, l, splitting, gamma);
            size_t nt = study_steps(prm, l, pb.dd.h_min);
            std::cout << bold << cyan << "      l = " << l << "  nt = " << nt << reset << std::endl;
            std::string summary = problem_summary(pb, explicit_scheme);
            std::cout << summary;
            log << summary;
            if (explicit_scheme) check_time_step(pb, prm, tf / nt);
            vector_type y = run(pb, nt, l == levels.back());
            print_splitting_stats(pb.system(), std::cout);
            rows.push_back({pb.dd.h_max, tf / nt, all.select(all_errors(pb.errors(y, tf))), pb.dd.h_min});
            print_table(names, rows, false, std::cout);
        }
        print_table(names, rows, false, log);
    }
    else {
        // against the exact solution only, or against a reference solution of the same integrator (mass and energy
        // norms of the difference) with the exact errors of the displacement / potential
        bool reference = prm.study.reference == "numerical";
        error_columns cols = reference ? columns(prm, {"elastic_u"}, {"acoustic_phi"}) : all;
        names = cols.names;
        if (reference) names.insert(names.begin(), {"mass_norm", "energy_norm"});
        problem pb(prm, prm.mesh.level, splitting, gamma);
        std::string summary = problem_summary(pb, explicit_scheme);
        std::cout << summary;
        log << summary;
        size_t n0 = prm.time.n_steps;
        size_t n_levels = prm.study.time_levels;
        if (explicit_scheme) check_time_step(pb, prm, tf / n0);
        vector_type y_ref;
        if (reference) {
            size_t nt_ref = prm.study.reference_factor * n0 * (size_t(1) << (n_levels - 1));
            y_ref = run(pb, nt_ref, false);
            std::cout << bold << cyan << "      reference nt = " << nt_ref << reset << std::endl;
            print_splitting_stats(pb.system(), std::cout);
        }
        for (size_t i = 0; i < n_levels; i++) {
            size_t nt = n0 * (size_t(1) << i);
            vector_type y = run(pb, nt, i + 1 == n_levels);
            auto e = pb.errors(y, tf);
            std::vector<RealType> errors = reference ? cols.select({e[1], e[3]}) : cols.select(all_errors(e));
            if (reference) errors.insert(errors.begin(), {pb.mass_distance(y, y_ref), pb.energy_distance(y, y_ref, tf)});
            rows.push_back({pb.dd.h_max, tf / nt, errors, pb.dd.h_min});
            print_table(names, rows, true, std::cout);
        }
        print_splitting_stats(pb.system(), std::cout);
        print_table(names, rows, true, log);
    }
    write_table(prm.table_file(), names, rows);
    std::cout << bold << red << "   --> " << fname << ", " << prm.table_file() << reset << std::endl << std::endl;
}

#endif /* HHOSecondOrder_hpp */
