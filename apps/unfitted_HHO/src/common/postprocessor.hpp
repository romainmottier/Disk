/*
 * Postprocessing of the unfitted HHO method stabilized by polynomial
 * extension (ported from ProtoN). The text output format is the one of
 * ProtoN, so that the generated python scripts can plot the results.
 */

#pragma once
#ifndef postprocessor_hpp
#define postprocessor_hpp

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/bases/cut_bases.hpp"
#include "operators.hpp"
#include "agglo_operators.hpp"
#include "quadrature_matrices.hpp"

template<disk::CutMesh Mesh>
class postprocessor {

    using T = typename Mesh::coordinate_type;

public:

    struct errors {
        T h, H1, L2;
        T H1_own;   // energy error of Theorem 4.6, with the own unknowns of the ill-cut subcells
        T L2_own;   // L2 error with the own unknowns of the ill-cut subcells
    };

    /* kappa-weighted squared H1 and L2 errors of the discrete function of
     * dofs "dofs" in the basis cb, on the quadrature qps */
    template<typename Basis, typename QPs, typename TestCase>
    static std::pair<T, T>
    squared_errors(const Basis& cb, const Eigen::Matrix<T, Eigen::Dynamic, 1>& dofs, const QPs& qps,
                   const TestCase& test_case, T kappa)
    {
        static const size_t DIM = Mesh::dimension;
        if (qps.empty())
            return std::make_pair(T(0), T(0));
        auto V  = basis_values(cb, qps);
        auto dV = basis_gradients<DIM>(cb, qps);
        Eigen::Matrix<T, Eigen::Dynamic, 1> u = V * dofs;
        std::array<Eigen::Matrix<T, Eigen::Dynamic, 1>, DIM> du;
        for (size_t d = 0; d < DIM; d++)
            du[d] = dV[d] * dofs;
        T H1 = 0, L2 = 0;
        for (size_t q = 0; q < qps.size(); q++) {
            auto g = test_case.sol_grad(qps[q].point());
            T e2 = 0;
            for (size_t d = 0; d < DIM; d++)
                e2 += (g(d) - du[d](q)) * (g(d) - du[d](q));
            H1 += kappa * qps[q].weight() * e2;
            T dv = test_case.sol_fun(qps[q].point()) - u(q);
            L2 += kappa * qps[q].weight() * dv * dv;
        }
        return std::make_pair(H1, L2);
    }

    /* Errors of the cell unknowns, weighted by kappa (Theorem 4.6). As in
     * ProtoN, on the ill-cut subcells S^i the polynomial of the paired subcell
     * T^i is extended to S^i; H1_own and L2_own use instead the own unknowns of S^i. */
    template<typename Assembler, typename TestCase>
    static errors
    compute_errors(const Mesh& msh, const std::vector<disk::cut_subcell>& P_OK,
                   const mesh_parameters<Mesh>& prm, const Assembler& assembler,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& x_dof, const TestCase& test_case)
    {
        T H1_error = 0.0;
        T L2_error = 0.0;
        T H1_own   = 0.0;
        T L2_own   = 0.0;
        T h = std::numeric_limits<T>::max();

        auto celdeg = prm.hdi.cell_degree();

        for (auto& P : P_OK) {
            auto cl  = msh[P.cell];
            auto loc = P.side;
            T kappa = (loc == disk::location::IN_NEGATIVE_SIDE) ? test_case.parms.kappa_1
                                                                : test_case.parms.kappa_2;
            h = std::min(h, diameter(msh, cl));

            disk::cut_scalar_basis<Mesh> cb(msh, cl, celdeg, loc, prm.basis);
            Eigen::Matrix<T, Eigen::Dynamic, 1> cell_dofs = assembler.take_cell_dofs(msh, cl, x_dof, loc);

            auto e = squared_errors(cb, cell_dofs, integrate(msh, cl, 2*celdeg, loc), test_case, kappa);
            H1_error += e.first;
            L2_error += e.second;
            H1_own   += e.first;
            L2_own   += e.second;

            for (auto dp : P.dependent_cells) {
                auto dp_cl = msh[dp];
                auto qps = integrate(msh, dp_cl, 2*celdeg, loc);
                auto ee = squared_errors(cb, cell_dofs, qps, test_case, kappa);
                H1_error += ee.first;
                L2_error += ee.second;

                disk::cut_scalar_basis<Mesh> cb_dp(msh, dp_cl, celdeg, loc, prm.basis);
                Eigen::Matrix<T, Eigen::Dynamic, 1> dp_dofs = assembler.take_cell_dofs(msh, dp_cl, x_dof, loc);
                auto eo = squared_errors(cb_dp, dp_dofs, qps, test_case, kappa);
                H1_own += eo.first;
                L2_own += eo.second;
            }
        }

        return errors{h, std::sqrt(H1_error), std::sqrt(L2_error), std::sqrt(H1_own), std::sqrt(L2_own)};
    }

    /* Errors of the method stabilized by cell agglomeration (same norms) */
    template<typename Assembler, typename TestCase>
    static errors
    compute_errors_agglo(const Mesh& msh, const std::vector<macro_cell<Mesh>>& macros,
                         const mesh_parameters<Mesh>& prm, const Assembler& assembler,
                         const Eigen::Matrix<T, Eigen::Dynamic, 1>& x_dof, const TestCase& test_case)
    {
        T H1_error = 0.0;
        T L2_error = 0.0;
        T h = std::numeric_limits<T>::max();

        auto celdeg = prm.hdi.cell_degree();

        for (auto& cl : msh)
            h = std::min(h, diameter(msh, cl));

        for (size_t q = 0; q < macros.size(); q++) {
            const auto& Q = macros[q];
            for (auto side : macro_sides(Q)) {
                T kappa = (side == disk::location::IN_NEGATIVE_SIDE) ? test_case.parms.kappa_1
                                                                     : test_case.parms.kappa_2;
                macro_bases<Mesh> mb(msh, Q, side, prm);
                Eigen::Matrix<T, Eigen::Dynamic, 1> cell_dofs = assembler.take_cell_dofs(q, Q, x_dof, side);
                for (auto m : Q.members) {
                    auto e = squared_errors(mb.cb, cell_dofs, integrate(msh, msh[m], 2*celdeg, side), test_case, kappa);
                    H1_error += e.first;
                    L2_error += e.second;
                }
            }
        }

        return errors{h, std::sqrt(H1_error), std::sqrt(L2_error), std::sqrt(H1_error), std::sqrt(L2_error)};
    }

    static void
    write_errors(const errors& err, const errors& previous, std::ostream& error_file)
    {
        T orderH = std::log(previous.H1 / err.H1) / std::log(previous.h / err.h);
        T orderL = std::log(previous.L2 / err.L2) / std::log(previous.h / err.h);
        error_file << "Characteristic h size = " << err.h << std::endl;
        error_file << "L2-norm error = " << err.L2 << std::endl;
        error_file << "H1-norm error = " << err.H1 << std::endl;
        error_file << "order L2 = " << orderL << std::endl;
        error_file << "order H1 = " << orderH << std::endl;
        error_file << "H1-norm error with own unknowns on ill-cut subcells = " << err.H1_own << std::endl;
        error_file << "L2-norm error with own unknowns on ill-cut subcells = " << err.L2_own << std::endl << std::endl;
    }

    static void 
    write_conv_sol(const std::string& txtFilename) {

        if (txtFilename.size() < 4 || txtFilename.substr(txtFilename.size() - 4) != ".txt") {
            std::cerr << "Erreur : Le fichier d'entrée n'a pas une extension .txt valide." << std::endl;
            return;
        }
        
        std::string pyFilename = txtFilename.substr(0, txtFilename.size() - 4) + ".py";
        
        std::ofstream pyFile(pyFilename);
        if (!pyFile.is_open()) {
            std::cerr << "Erreur : Impossible de créer le fichier Python." << std::endl;
            return;
        }
        
      pyFile << "import matplotlib.pyplot as plt\n";
      pyFile << "import numpy as np\n";
      pyFile << "import sys\n\n";
      
      pyFile << "error_type = \"L2\"\n";
      pyFile << "if len(sys.argv) > 1:\n";
      pyFile << "    if sys.argv[1] in [\"L2\", \"H1\"]:\n";
      pyFile << "        error_type = sys.argv[1]\n";
      pyFile << "    else:\n";
      pyFile << "        print(\"Argument invalide. Utilisez 'L2' ou 'H1'.\")\n";
      pyFile << "        sys.exit(1)\n\n";
      
      pyFile << "filename = \"" << txtFilename << "\"\n";
      pyFile << "with open(filename, \"r\") as file:\n";
      pyFile << "    lines = file.readlines()\n\n";
      
      pyFile << "results_L2 = {}\n";
      pyFile << "results_H1 = {}\n";
      pyFile << "h_values = {}\n\n";
      
      pyFile << "current_degree = None\n";
      pyFile << "for line in lines:\n";
      pyFile << "    line = line.strip()\n";
      pyFile << "    if line.startswith(\"Polynomial degree k :\"):\n";
      pyFile << "        current_degree = int(line.split(\":\")[1].strip())\n";
      pyFile << "        if current_degree not in results_L2:\n";
      pyFile << "            results_L2[current_degree] = []\n";
      pyFile << "            results_H1[current_degree] = []\n";
      pyFile << "            h_values[current_degree] = []\n";
      pyFile << "    elif line.startswith(\"Characteristic h size =\"):\n";
      pyFile << "        h = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        h_values[current_degree].append(h)\n";
      pyFile << "    elif line.startswith(\"L2-norm error =\"):\n";
      pyFile << "        L2_error = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        results_L2[current_degree].append(L2_error)\n";
      pyFile << "    elif line.startswith(\"H1-norm error =\"):\n";
      pyFile << "        H1_error = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        results_H1[current_degree].append(H1_error)\n\n";
      
      pyFile << "plt.figure(figsize=(10, 6))\n\n";
      
      pyFile << "if error_type == \"L2\":\n";
      pyFile << "    results = results_L2\n";
      pyFile << "    ylabel = \"L2-norm Error\"\n";
      pyFile << "    title = \"L2-norm Error vs h\"\n";
      pyFile << "elif error_type == \"H1\":\n";
      pyFile << "    results = results_H1\n";
      pyFile << "    ylabel = \"H1-norm Error\"\n";
      pyFile << "    title = \"H1-norm Error vs h\"\n\n";
      
      pyFile << "for degree in sorted(results.keys()):\n";
      pyFile << "    h_values_np = np.array(h_values[degree])\n";
      pyFile << "    results_np = np.array(results[degree])\n";
      pyFile << "    log_h = np.log(h_values_np)\n";
      pyFile << "    log_error = np.log(results_np)\n";
      pyFile << "    slope, _ = np.polyfit(log_h, log_error, 1)\n";
      pyFile << "    plt.loglog(h_values_np, results_np, marker='o', label=f\"k={degree} {ylabel} (rate={slope:.2f})\")\n\n";
      
      pyFile << "plt.xlabel(\"h\", fontsize=12)\n";
      pyFile << "plt.ylabel(ylabel, fontsize=12)\n";
      pyFile << "plt.title(title, fontsize=14)\n";
      pyFile << "plt.legend()\n";
      pyFile << "plt.grid(which=\"both\", linestyle=\"--\", linewidth=0.5)\n";
      pyFile << "plt.tight_layout()\n\n";
      
      pyFile << "plt.show()\n";
      pyFile.close();
      
    }

};

#endif // postprocessor_hpp
