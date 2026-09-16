//
//  LI_SDIRK_LTS_MeshConv_test.hpp
//
//  Rigorous convergence test for LI-SDIRK(4,3), mirroring
//  prototypes/LTS/ERK4_LTS_conv_test.hpp (the test used to produce the
//  paper's Section 6 LTS-ERK(4,p) convergence results) as closely as
//  possible, so the numbers are directly comparable:
//
//    - SAME mesh construction: cartesian_2d_mesh_builder.refine_mesh(l)
//      base mesh + the SAME small local refinement around pt1(-0.5,0.5)/
//      pt2(0.5,0.5) via find_cells + refine_cells(vec, substeps_Q) --
//      most of the domain stays coarse, only cells near these two points
//      get the extra local refinement (this is what makes "most of the
//      domain coarse, only a few cells fine" true here).
//    - SAME nt scaling law: nt = sim_data.m_nt_divs * 2^(l-1), so dt
//      shrinks proportionally to the coarse mesh size h_max as l grows,
//      keeping the coarse-region explicit operator's CFL number roughly
//      CONSTANT across refinement levels -- this is why the original
//      ERK4_LTS_conv_test.hpp does not need manual CFL bisection per
//      level the way the earlier fixed-mesh, Delta-t-only sweep
//      (LI_SDIRK_LTS_conv_test.hpp) did.
//    - SAME error routines: postprocessor::compute_errors_four_fields_
//      elastoacoustic[_energy_norm], SAME log file content format
//      (so the existing generate_LTS_convergence_script-style parsing
//      could be reused directly if desired).
//    - Loops over nb_layer in {1,2} (the assemble_P "protection layer"
//      overlap width), reporting both -- ERK4_LTS_conv_test.hpp itself
//      only ever used nb_layer=1 (assemble_P(msh,h_c) with the default).
//
//  ONLY the time integrator differs from ERK4_LTS_conv_test.hpp: the
//  p-substep explicit RK4 fine loop is replaced by ONE call to the new
//  dirk_weight_LTS_fine (implicit SDIRK(3,4), single global dt), and
//  erk_weight_LTS_coarse is called UNCHANGED exactly as in the original
//  (safe here: a fresh erk_coupling_hho_scheme is built per (l,nb_layer)
//  run, so its one-shot active-dof cache is never reused across a
//  different Pcoarse).
//
//  Analytical solution: EFunctionNonPolynomial (same as the original
//  test and as LI_SDIRK_LTS_conv_test.hpp) -- see the write-up for why
//  EFunctionPlaneWaveAcoustic was NOT used: its elastic-side methods
//  (Evaluate_u/v/a/f/sigma) have NO case for it in
//  scal_vec_analytic_functions.hpp and silently fall through to
//  "Function not implemented" returning an UNINITIALIZED
//  disk::static_vector (not even zero) -- using it as-is would corrupt
//  the elastic half of this coupled test. EFunctionNonPolynomial is
//  already sinusoidal in x,y AND t (e.g. x^2*sin(pi*y)*cos(pi*x/2)*
//  cos(sqrt(2)*pi*t) on the elastic side), just with a polynomial
//  envelope rather than a pure plane wave.
//
//  Created by Claude (Sonnet 5), for Romain Mottier -- 2026-09-14.
//

#include "dirk_coupling_hho_scheme.hpp"
#include <regex>

void LI_SDIRK_LTS_MeshConv_test(int argc, char **argv);

namespace li_sdirk_meshconv_detail {
    inline double combined_L2_error(const std::string & log_text) {
        std::regex re("L2-norm error = ([0-9eE+\\-\\.naninfNANIF]+)");
        auto begin = std::sregex_iterator(log_text.begin(), log_text.end(), re);
        auto end   = std::sregex_iterator();
        std::vector<double> vals;
        for (auto it = begin; it != end; ++it) vals.push_back(std::stod((*it)[1].str()));
        if (vals.size() < 2) return -1.0;
        return std::sqrt(vals[0]*vals[0] + vals[1]*vals[1]);
    }
    inline double combined_energy_error(const std::string & log_text) {
        // compute_errors_four_fields_elastoacoustic_energy_norm writes
        // "L2 errors of dG unknowns:" then "Acoustic region : <val>" and
        // "Elastic  region : <val>" (see postprocessor.hpp:106-108) --
        // NOT "dG-norm error = ". Match both "region : <value>" lines.
        std::regex re("region\\s*:\\s*([0-9eE+\\-\\.naninfNANIF]+)");
        auto begin = std::sregex_iterator(log_text.begin(), log_text.end(), re);
        auto end   = std::sregex_iterator();
        std::vector<double> vals;
        for (auto it = begin; it != end; ++it) vals.push_back(std::stod((*it)[1].str()));
        if (vals.size() < 2) return -1.0;
        return std::sqrt(vals[0]*vals[0] + vals[1]*vals[1]);
    }
}

void LI_SDIRK_LTS_MeshConv_test(int argc, char **argv){

    std::cout << std::endl << bold << red << "   LI-SDIRK(4,3) LTS -- MESH-FAMILY CONVERGENCE TEST (mirrors ERK4_LTS_conv_test.hpp)" << reset << std::endl;

    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;

    Matrix<double, Dynamic, Dynamic> a_tab;
    Matrix<double, Dynamic, 1> b_tab;
    Matrix<double, Dynamic, 1> c_tab;
    dirk_butcher_tableau::sdirk_tables(3, a_tab, b_tab, c_tab);
    RealType gamma = a_tab(0,0);

    for (size_t k = 3; k <= 3; k++) {

        std::cout << std::endl << bold << red << "   Polynomial degree k : " << k << reset << std::endl;
        size_t cell_k_degree = k;
        if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
        disk::hho_degree_info hho_di(cell_k_degree, k);

        for (size_t nb_layer : {size_t(1), size_t(2)}) {

            std::cout << std::endl << bold << red << "   ================================================================" << reset << std::endl;
            std::cout << bold << red << "   nb_layer (assemble_P protection layers) = " << nb_layer << reset << std::endl;
            std::cout << bold << red << "   ================================================================" << reset << std::endl;

            struct Row { size_t l; RealType h_max, h_min; size_t nt; RealType dt; double L2err, Eerr; };
            std::vector<Row> rows;

            for (size_t l = 1; l <= sim_data.m_n_divs; l++) {

                std::cout << bold << cyan << "      Space refinement level -l : " << l << reset << std::endl;

                // ##################################################  Mesh (SAME as ERK4_LTS_conv_test.hpp)
                mesh_type msh;
                {
                    RealType lx = 2.0, ly = 1.0;
                    size_t nx = 4, ny = 2;
                    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
                    mesh_builder.refine_mesh(l);
                    mesh_builder.set_translation_data(-1.0, 0.0);
                    mesh_builder.build_mesh();
                    mesh_builder.move_to_mesh_storage(msh);

                    cartesian_2d_mesh_builder<RealType> mesh_builder2(lx,ly,nx,ny);
                    mesh_builder2.refine_mesh(l);
                    mesh_builder2.set_translation_data(-1.0, 0.0);
                    mesh_builder2.build_mesh();
                    typename mesh_type::point_type pt1(-0.5, 0.5);
                    typename mesh_type::point_type pt2( 0.5, 0.5);
                    std::set<size_t> cell_indexes1 = postprocessor<mesh_type>::find_cells(pt1, msh, true);
                    std::set<size_t> cell_indexes2 = postprocessor<mesh_type>::find_cells(pt2, msh, true);
                    std::vector<size_t> vec;
                    vec.insert(vec.end(), cell_indexes1.begin(), cell_indexes1.end());
                    vec.insert(vec.end(), cell_indexes2.begin(), cell_indexes2.end());
                    mesh_builder2.refine_cells(vec, sim_data.m_substeps_Q);
                    mesh_builder2.move_to_mesh_storage(msh);
                }

                RealType h_max = 1e-5, h_min = 10;
                for (auto & cell : msh) {
                    RealType h_l = diameter(msh, cell);
                    if (h_l < h_min) h_min = h_l;
                    else if (h_l > h_max) h_max = h_l;
                }
                int p_ratio = static_cast<int>(std::round(h_max / h_min));
                auto h_c = (p_ratio == 1) ? 1.25*h_max : 0.75*h_max;
                std::cout << bold << yellow << "            h_max/h_min = " << p_ratio << reset << std::endl;

                // ##################################################  Time discretization (SAME scaling law)
                size_t nt = sim_data.m_nt_divs * (size_t)std::pow(2, (double)l - 1.0);
                RealType ti = 0.0, tf = 1.0, t = ti;
                RealType dt = (tf - ti) / RealType(nt);

                // ##################################################  Manufactured solution
                scal_vec_analytic_functions functions;
                functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
                auto u_fun = functions.Evaluate_u(t);
                auto v_fun = functions.Evaluate_v(t);
                auto f_fun = functions.Evaluate_f(t);
                auto flux_fun = functions.Evaluate_sigma(t);
                auto s_u_fun = functions.Evaluate_s_u(t);
                auto s_v_fun = functions.Evaluate_s_v(t);
                auto s_f_fun = functions.Evaluate_s_f(t);
                auto s_flux_fun = functions.Evaluate_s_q(t);

                // ##################################################  Materials / interface / BCs (SAME)
                auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
                    return elastic_material_data<RealType>(1.0, std::sqrt(3.0), 1.0);
                };
                auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
                    return acoustic_material_data<RealType>(1.0, 1.0);
                };
                std::map<size_t,elastic_material_data<RealType>> e_material;
                std::map<size_t,acoustic_material_data<RealType>> a_material;
                std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
                std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
                RealType eps = 1.0e-10;
                for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
                    mesh_type::point_type bar = barycenter(msh, *face_it);
                    auto fc_id = msh.lookup(*face_it);
                    if (std::fabs(bar.x()) < eps) interface_face_indexes.insert(fc_id);
                }
                for (auto & cell : msh) {
                    auto cell_ind = msh.lookup(cell);
                    mesh_type::point_type bar = barycenter(msh, cell);
                    if (bar.x() > 0) a_material.insert(std::make_pair(cell_ind, acoustic_mat_fun(bar)));
                    else e_material.insert(std::make_pair(cell_ind, elastic_mat_fun(bar)));
                    for (auto face : faces(msh, cell)) {
                        auto fc_id = msh.lookup(face);
                        if (interface_face_indexes.count(fc_id)) {
                            if (bar.x() > 0) interface_cell_pair_indexes[fc_id].second = cell_ind;
                            else interface_cell_pair_indexes[fc_id].first = cell_ind;
                        }
                    }
                }
                std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
                for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
                    mesh_type::point_type bar = barycenter(msh, *face_it);
                    auto fc_id = msh.lookup(*face_it);
                    if (interface_face_indexes.count(fc_id)) {
                        if (bar.y() > 0) acoustic_internal_faces.insert(fc_id);
                        else elastic_internal_faces.insert(fc_id);
                    }
                }
                size_t bc_elastic_id = 0, bc_acoustic_id = 1;
                for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
                    mesh_type::point_type bar = barycenter(msh, *face_it);
                    auto fc_id = msh.lookup(*face_it);
                    if (bar.x() > 0) {
                        msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_acoustic_id, true};
                        acoustic_bc_face_indexes.insert(fc_id);
                    } else {
                        msh.backend_storage()->boundary_info.at(fc_id) = disk::boundary_descriptor{bc_elastic_id, true};
                        elastic_bc_face_indexes.insert(fc_id);
                    }
                }
                e_boundary_type e_bnd(msh);
                a_boundary_type a_bnd(msh);
                e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
                a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);

                // ##################################################  Assembly
                auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
                assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
                assembler.set_coupling_stabilization();
                if (sim_data.m_scaled_stabilization_Q) assembler.set_scaled_stabilization();
                assembler.assemble_mass(msh);
                assembler.assemble_coupling_terms(msh);

                Matrix<RealType, Dynamic, 1> x_dof;
                assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
                assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);

                assembler.assemble(msh, f_fun, s_f_fun, true);
                assembler.LHS += assembler.COUPLING;

                size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
                size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
                size_t e_face_dofs        = assembler.get_e_face_dof();
                size_t a_face_dofs        = assembler.get_a_face_dof();

                erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
                                                         elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
                erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
                                   assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
                erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
                                   assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
                                   assembler.get_e_compress(), assembler.get_a_compress(),
                                   elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
                erk_an.refresh_faces_unknowns(x_dof);

                // ##################################################  Fine-region projector + implicit operator
                assembler.assemble_P(msh, h_c, nb_layer);
                dirk_coupling_hho_scheme<RealType> dirk_op;
                dirk_op.compose(erk_an, assembler.Pfine, gamma, dt);

                // ##################################################  Log
                std::ostringstream filename;
                filename << "li_sdirk_lts_l_" << l << "_n_" << sim_data.m_nt_divs
                         << "_k_" << k << "_s_sdirk34_nblayer_" << nb_layer
                         << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
                std::ofstream simulation_log(filename.str());
                sim_data.write_simulation_data(simulation_log);
                simulation_log << "Number of time steps = " << nt << std::endl;
                simulation_log << "Step size =  " << dt << std::endl;
                simulation_log << "nb_layer = " << nb_layer << std::endl;
                simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
                simulation_log << "Space step = " << h_max << std::endl;
                simulation_log.flush();

                // ##################################################  Source term lambda
                auto eval_F = [&](RealType t_abs) -> Matrix<RealType, Dynamic, 1> {
                    t = t_abs;
                    auto v_fun_t   = functions.Evaluate_v(t);
                    auto f_fun_t   = functions.Evaluate_f(t);
                    auto s_v_fun_t = functions.Evaluate_s_v(t);
                    auto s_f_fun_t = functions.Evaluate_s_f(t);
                    assembler.get_e_bc_conditions().updateDirichletFunction(v_fun_t, 0);
                    assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
                    assembler.assemble_rhs(msh, f_fun_t, s_f_fun_t, true);
                    return assembler.RHS;
                };
                std::function<Matrix<RealType, Dynamic, 1>(RealType)> eval_F_std = eval_F;

                // ##################################################  Time marching: LI-SDIRK(4,3)
                for (size_t it = 1; it <= nt; it++) {

                    RealType tn = dt*(it-1) + ti;
                    auto x_dof_n = x_dof;

                    std::vector<Matrix<RealType, Dynamic, 1>> w(4);
                    for (int i = 0; i < 4; ++i) { w[i].resize(x_dof.rows()); w[i].setZero(); }

                    erk_an.ZeroFc();
                    RealType tn12 = tn + 0.5*dt, tn1 = tn + dt;
                    Matrix<RealType, Dynamic, 1> Fn   = eval_F(tn);
                    Matrix<RealType, Dynamic, 1> Fn12 = eval_F(tn12);
                    Matrix<RealType, Dynamic, 1> Fn1  = eval_F(tn1);
                    erk_an.ZeroFc();
                    // Exactly as ERK4_LTS.hpp / ERK4_LTS_conv_test.hpp: erk_weight_LTS_coarse
                    // UNCHANGED, one fixed Pcoarse for this (l,nb_layer) run's whole lifetime.
                    erk_an.erk_weight_LTS_coarse(x_dof_n, assembler.Pcoarse, w, Fn, Fn12, Fn1, dt);

                    dirk_weight_LTS_fine<RealType>(erk_an, dirk_op, x_dof_n, assembler.Pfine, w, tn, dt, eval_F_std, a_tab, b_tab, c_tab, gamma);

                    x_dof = x_dof_n;
                    t = tn + dt;
                }

                // ##################################################  Errors
                t = tf;
                auto v_fun_T      = functions.Evaluate_v(t);
                auto flux_fun_T   = functions.Evaluate_sigma(t);
                auto s_v_fun_T    = functions.Evaluate_s_v(t);
                auto s_flux_fun_T = functions.Evaluate_s_q(t);

                std::ostringstream l2_log, e_log;
                postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler, x_dof, v_fun_T, flux_fun_T, s_v_fun_T, s_flux_fun_T, l2_log);
                postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic_energy_norm(msh, hho_di, assembler, x_dof, v_fun_T, flux_fun_T, s_v_fun_T, s_flux_fun_T, e_log);
                simulation_log << l2_log.str() << e_log.str();

                double L2err = li_sdirk_meshconv_detail::combined_L2_error(l2_log.str());
                double Eerr  = li_sdirk_meshconv_detail::combined_energy_error(e_log.str());

                std::cout << bold << cyan << "      l=" << l << "  h_max=" << h_max << "  h_min=" << h_min
                          << "  nt=" << nt << "  dt=" << dt
                          << "  L2err=" << std::setprecision(8) << L2err
                          << "  Eerr=" << std::setprecision(8) << Eerr << reset << std::endl;

                rows.push_back({l, h_max, h_min, nt, dt, L2err, Eerr});
            }

            std::cout << std::endl << bold << yellow << "   ---- Mesh-family convergence table (nb_layer=" << nb_layer << ") ----" << reset << std::endl;
            std::cout << std::setw(4) << "l" << std::setw(12) << "h_max" << std::setw(12) << "h_min"
                       << std::setw(8) << "nt" << std::setw(14) << "dt"
                       << std::setw(16) << "L2 error" << std::setw(10) << "order"
                       << std::setw(16) << "Energy error" << std::setw(10) << "order" << std::endl;
            for (size_t i = 0; i < rows.size(); i++) {
                std::cout << std::setw(4) << rows[i].l
                          << std::setw(12) << std::setprecision(6) << rows[i].h_max
                          << std::setw(12) << std::setprecision(6) << rows[i].h_min
                          << std::setw(8)  << rows[i].nt
                          << std::setw(14) << std::setprecision(6) << rows[i].dt
                          << std::setw(16) << std::setprecision(8) << rows[i].L2err;
                if (i > 0 && rows[i].L2err > 0 && rows[i-1].L2err > 0) {
                    double order = std::log(rows[i-1].L2err/rows[i].L2err) / std::log(rows[i-1].h_max/rows[i].h_max);
                    std::cout << std::setw(10) << std::setprecision(4) << order;
                } else std::cout << std::setw(10) << "--";
                std::cout << std::setw(16) << std::setprecision(8) << rows[i].Eerr;
                if (i > 0 && rows[i].Eerr > 0 && rows[i-1].Eerr > 0) {
                    double order = std::log(rows[i-1].Eerr/rows[i].Eerr) / std::log(rows[i-1].h_max/rows[i].h_max);
                    std::cout << std::setw(10) << std::setprecision(4) << order;
                } else std::cout << std::setw(10) << "--";
                std::cout << std::endl;
            }
        }
    }

    std::cout << std::endl << bold << red << "   LI-SDIRK(4,3) MESH-FAMILY CONVERGENCE TEST -- DONE" << reset << std::endl << std::endl;
}
