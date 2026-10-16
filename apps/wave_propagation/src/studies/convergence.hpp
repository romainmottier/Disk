//
//  convergence.hpp
//
//  Convergence studies of the two formulations: rows and printing of the tables (observed orders), data tables,
//  mesh levels and numbers of time steps of a study, selection of the error columns.
//

#pragma once
#ifndef studies_convergence_hpp
#define studies_convergence_hpp

#include "../mesh/domain.hpp"

namespace coupling_common {

    // One row of a convergence table: step sizes (h = largest, h_min = smallest cell diameter) and errors
    struct conv_row {
        RealType h, dt;
        std::vector<RealType> errors;
        RealType h_min = 0.0;
    };

    // Prints the errors and the observed orders log(e_{i-1}/e_i) / log(r_{i-1}/r_i), with r = h or dt
    inline void print_table(const std::vector<std::string>& names, const std::vector<conv_row>& rows,
                            bool rate_in_dt, std::ostream& os) {
        auto width = [&](size_t j) { return int(std::max<size_t>(16, names[j].size() + 2)); };
        os << std::setw(12) << "h" << std::setw(12) << "dt";
        for (size_t j = 0; j < names.size(); j++)
            os << std::setw(width(j)) << names[j] << std::setw(8) << "order";
        os << std::endl;
        for (size_t i = 0; i < rows.size(); i++) {
            os << std::setw(12) << std::setprecision(4) << rows[i].h << std::setw(12) << std::setprecision(4) << rows[i].dt;
            for (size_t j = 0; j < rows[i].errors.size(); j++) {
                os << std::setw(width(j)) << std::setprecision(6) << rows[i].errors[j];
                if (i > 0 && rows[i].errors[j] > 0 && rows[i-1].errors[j] > 0) {
                    RealType r = rate_in_dt ? rows[i-1].dt / rows[i].dt : rows[i-1].h / rows[i].h;
                    os << std::setw(8) << std::setprecision(3) << std::log(rows[i-1].errors[j] / rows[i].errors[j]) / std::log(r);
                }
                else
                    os << std::setw(8) << "--";
            }
            os << std::endl;
        }
        os.flush();
    }

    // Data table: "# h h_min dt name..." header, then one row per level (errors without orders), full precision
    inline void write_table(const std::string& file, const std::vector<std::string>& names, const std::vector<conv_row>& rows) {
        std::ofstream os(file);
        os << "# h h_min dt";
        for (auto& n : names)
            os << " " << n;
        os << std::endl << std::setprecision(16) << std::scientific;
        for (auto& r : rows) {
            os << r.h << " " << (r.h_min > 0.0 ? r.h_min : r.h) << " " << r.dt;
            for (auto e : r.errors)
                os << " " << e;
            os << std::endl;
        }
    }

    // Energy-type norm ||x||_M = sqrt(x^T M x) of a vector of unknowns for a mass matrix M
    inline RealType mass_norm(const SparseMatrix<RealType>& M, const vector_type& x) {
        return std::sqrt(std::abs(x.dot(M * x)));
    }

    // Levels and time steps of the study: "run" (mesh.level, n_steps), "space" (mesh.level .. last_level, n_steps
    // doubled with each level if study.refine_dt, fixed otherwise), "time" (mesh.level, see the prototypes)
    inline std::vector<size_t> study_levels(const wave_parameters& prm) {
        std::vector<size_t> levels;
        size_t last = prm.study.type == "space" || prm.study.type == "spectrum" ? prm.study.last_level : prm.mesh.level;
        for (size_t l = prm.mesh.level; l <= last; l++)
            levels.push_back(l);
        return levels;
    }

    // time.dt_h > 0: dt = dt_h h_min on each mesh (smallest cell diameter, as the h of Burman-Duran-Ern 2022 on
    // polygonal meshes), with as many steps as needed to reach the final time
    inline size_t study_steps(const wave_parameters& prm, size_t level, RealType h_min) {
        if (prm.time.dt_h > 0.0)
            return size_t(std::ceil(prm.time.final_time / (prm.time.dt_h * h_min) - 1.0e-9));
        if (!prm.study.refine_dt) return prm.time.n_steps;
        if (prm.time.steps_growth != 2.0)
            return std::max<size_t>(1, size_t(std::llround(RealType(prm.time.n_steps) *
                                                            std::pow(prm.time.steps_growth, RealType(level - prm.mesh.level)))));
        return prm.time.n_steps * (size_t(1) << (level - prm.mesh.level));
    }

    // Columns of the convergence tables present for the medium of the problem: the names of the elastic, then of the
    // acoustic errors, filtered by problem.medium; select picks the matching entries of a full error vector
    struct error_columns {
        std::vector<std::string> names;
        std::vector<size_t> index;
        std::vector<RealType> select(const std::vector<RealType>& all) const {
            std::vector<RealType> e;
            for (size_t i : index) e.push_back(all[i]);
            return e;
        }
    };

    inline error_columns columns(const wave_parameters& prm, const std::vector<std::string>& elastic,
                                 const std::vector<std::string>& acoustic) {
        error_columns c;
        for (size_t i = 0; i < elastic.size(); i++)
            if (prm.has_elastic()) { c.names.push_back(elastic[i]); c.index.push_back(i); }
        for (size_t i = 0; i < acoustic.size(); i++)
            if (prm.has_acoustic()) { c.names.push_back(acoustic[i]); c.index.push_back(elastic.size() + i); }
        return c;
    }

}

#endif /* studies_convergence_hpp */
