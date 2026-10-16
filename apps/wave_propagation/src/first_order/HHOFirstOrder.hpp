//
//  HHOFirstOrder.hpp
//
//  Elasto-acoustic coupling, first-order formulation (velocity-stress in the solid, velocity-pressure in the fluid,
//  elastoacoustic_four_fields_assembler): M y' + K y = F(t). Every setting comes from the Lua input file
//  (parameters/lua_parameters.hpp, template share/first_order.lua):
//    space : face degree k, equal (k,k) or mixed (k+1,k) order, least-squares stabilization (plain least squares in
//            equal order, Lehrenfeld-Schoberl in mixed order) or HHO stabilization, weight O(1) or in 1/h
//    time  : ERK(s), SDIRK(s,s+1), two-level LTS-RK4 (substeps in the locally refined region), LI-SDIRK(4,3) (coarse
//            ERK4 Taylor polynomial, SDIRK(3,4) in the fine region), final time, number of time steps
//    mesh, problem (manufactured solution), study (one run, space-time convergence, time convergence), output
//  The explicit integrators (erk, lts) and li eliminate the faces with the block-diagonal face block of the
//  least-squares stabilization (erk_coupling_hho_scheme::Sff_inverse); sdirk accepts both stabilizations.
//  solver.matrix_free (fast mode, same results up to round-off): erk, lts, li and the cfl study on the matrix-free
//  operator of time_schemes/mf_erk_operator.hpp (local blocks shared between cells, OpenMP); sdirk with parallel products and
//  Pardiso without iterative refinement (dirk_hho_scheme::set_fast); parallel assembly without the zeros of the local
//  matrices, time-independent data of F(t) kept. Without it, the same operations as before.
//  Papers' settings (MEG, MEKG): explicit = equal order, plain least squares, O(1); implicit = mixed order,
//  Lehrenfeld-Schoberl, 1/h.
//  space.elastic_formulation = "three_fields" (elastic medium, sdirk): the formulation and the assembler of
//  Burman-Duran-Ern-Steins 2021 (elastodynamic_three_fields_assembler: stress 2 mu eps, hydrostatic part lambda div u
//  and velocity), instead of the velocity-stress block of the coupled assembler.
//

#pragma once
#ifndef HHOFirstOrder_hpp
#define HHOFirstOrder_hpp

// problem.hpp (problems), observers.hpp (sensors, snapshots), time_schemes.hpp (integrators), studies.hpp (studies)
#include "studies.hpp"

void HHOFirstOrder(const wave_parameters& prm);

void HHOFirstOrder(const wave_parameters& prm) {

    using namespace hho_first_order;

    std::cout << std::endl << bold << red << "   FIRST-ORDER ELASTO-ACOUSTIC COUPLING, " << prm.scheme_label()
              << ", " << prm.space_label() << reset << std::endl;
    std::string fname = prm.output_file();
    std::ofstream log(fname);
    log << "first-order formulation" << std::endl;
    prm.print(log);

    // errors at the final time: elastic velocity and stress, acoustic scalar (time derivative of the potential) and
    // flux (gradient of the potential), L2 norms on the cell unknowns
    std::vector<conv_row> rows;
    std::vector<std::string> names;
    if (prm.study.type == "cfl")
        cfl_study(prm, log, rows, names);
    else if (prm.study.type == "rho")
        rho_study(prm, log, rows, names);
    else if (prm.study.type == "spectrum")
        spectrum_study(prm, log, rows, names);
    else if (prm.space.elastic_formulation == "three_fields")
        study<problem_three_fields>(prm, log, rows, names);
    else
        study<problem>(prm, log, rows, names);
    write_table(prm.table_file(), names, rows);
    std::cout << bold << red << "   --> " << fname << ", " << prm.table_file() << reset << std::endl << std::endl;
}

#endif /* HHOFirstOrder_hpp */
