//
//  observers.hpp
//
//  First-order formulation (HHOFirstOrder.hpp): outputs during the integration, sensors (values at points of the
//  cells and of the interface) and snapshots (cell values at given times).
//

#pragma once
#ifndef first_order_observers_hpp
#define first_order_observers_hpp

#include "problem.hpp"

namespace hho_first_order {

    // Point and segment geometry of the sensors: is pt in the (closed) polygon of a cell, within a tolerance eps?
    inline bool cell_contains(const mesh_type& msh, const typename mesh_type::cell_type& cell, const point_type& pt, RealType eps) {
        auto pts = points(msh, cell);
        const size_t n = pts.size();
        bool inside = false;
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            const auto& a = pts[j];
            const auto& b = pts[i];
            // on the edge [a, b]
            RealType ex = b.x() - a.x(), ey = b.y() - a.y();
            RealType len2 = ex * ex + ey * ey;
            RealType s = len2 > 0.0 ? ((pt.x() - a.x()) * ex + (pt.y() - a.y()) * ey) / len2 : 0.0;
            s = std::min(1.0, std::max(0.0, s));
            RealType dx = a.x() + s * ex - pt.x(), dy = a.y() + s * ey - pt.y();
            if (dx * dx + dy * dy <= eps * eps)
                return true;
            if (((a.y() > pt.y()) != (b.y() > pt.y())) &&
                (pt.x() < (b.x() - a.x()) * (pt.y() - a.y()) / (b.y() - a.y()) + a.x()))
                inside = !inside;
        }
        return inside;
    }

    inline bool face_contains(const mesh_type& msh, const typename mesh_type::face_type& face, const point_type& pt, RealType eps) {
        auto pts = points(msh, face);
        const auto& a = pts[0];
        const auto& b = pts[1];
        RealType ex = b.x() - a.x(), ey = b.y() - a.y();
        RealType len2 = ex * ex + ey * ey;
        RealType s = len2 > 0.0 ? ((pt.x() - a.x()) * ex + (pt.y() - a.y()) * ey) / len2 : 0.0;
        if (s < -1.0e-12 || s > 1.0 + 1.0e-12) return false;
        RealType dx = a.x() + s * ex - pt.x(), dy = a.y() + s * ey - pt.y();
        return dx * dx + dy * dy <= eps * eps;
    }

    // Point sensors (output.sensors): one CSV file <table>_<name>.csv per sensor, a line "t,values" every
    // output.sensor_every steps. Fields: cell values at the point (acoustic_pressure p_T, acoustic_flux m_T,
    // elastic_velocity v_T, elastic_stress s_T: xx, xy, yy) in a cell of the corresponding medium containing the point,
    // and at a point of the interface the face unknowns interface_pressure p_F (acoustic face) and interface_velocity
    // v_F (elastic face) of the interface face containing it.
    class sensors {
        struct item {
            std::string name;
            point_type pt;
            std::vector<std::string> fields;
            long e_cell = -1, a_cell = -1;      // cells of the media containing the point
            long e_face = -1, a_face = -1;      // position of the interface face among the faces of the cell
            std::ofstream os;
        };
        std::vector<std::unique_ptr<item>> m_items;
        size_t m_every = 1;

    public:
        sensors(problem& pb, const wave_parameters& prm, const std::string& base) : m_every(prm.output.sensor_every) {
            auto& msh = pb.msh;
            const RealType eps = 1.0e-8 * std::max(RealType(1.0), pb.dd.h_max);
            for (auto& sn : prm.output.sensors) {
                auto it = std::make_unique<item>();
                it->name = sn.name;
                it->pt = point_type(sn.point[0], sn.point[1]);
                it->fields = sn.fields;
                for (auto& cell : msh) {
                    auto ci = msh.lookup(cell);
                    if (!cell_contains(msh, cell, it->pt, eps)) continue;
                    if (pb.dd.e_material.count(ci) && it->e_cell < 0) it->e_cell = long(ci);
                    if (pb.dd.a_material.count(ci) && it->a_cell < 0) it->a_cell = long(ci);
                }
                auto interface_face = [&](long ci) -> long {
                    auto& cell = msh.backend_storage()->surfaces[size_t(ci)];
                    auto fcs = faces(msh, cell);
                    for (size_t j = 0; j < fcs.size(); j++)
                        if (pb.dd.interface_face_indexes.count(msh.lookup(fcs[j])) && face_contains(msh, fcs[j], it->pt, eps))
                            return long(j);
                    return -1;
                };
                std::string header = "t";
                for (auto& f : it->fields) {
                    bool acoustic = f == "acoustic_pressure" || f == "acoustic_flux" || f == "interface_pressure";
                    long ci = acoustic ? it->a_cell : it->e_cell;
                    if (ci < 0)
                        throw std::runtime_error("sensor " + it->name + ": no " + std::string(acoustic ? "acoustic" : "elastic")
                                                 + " cell contains the point for the field " + f);
                    if (f == "interface_pressure") {
                        it->a_face = interface_face(it->a_cell);
                        if (it->a_face < 0) throw std::runtime_error("sensor " + it->name + ": the point is not on the interface");
                    }
                    if (f == "interface_velocity") {
                        it->e_face = interface_face(it->e_cell);
                        if (it->e_face < 0) throw std::runtime_error("sensor " + it->name + ": the point is not on the interface");
                    }
                    if (f == "acoustic_pressure")  header += ",p";
                    if (f == "acoustic_flux")      header += ",mx,my";
                    if (f == "elastic_velocity")   header += ",vx,vy";
                    if (f == "elastic_stress")     header += ",sxx,sxy,syy";
                    if (f == "interface_pressure") header += ",pF";
                    if (f == "interface_velocity") header += ",vFx,vFy";
                }
                it->os.open(base + "_" + it->name + ".csv");
                it->os << "# sensor " << it->name << " (" << sn.point[0] << ", " << sn.point[1] << ")" << std::endl
                       << header << std::endl << std::setprecision(12) << std::scientific;
                m_items.push_back(std::move(it));
            }
        }

        bool empty() const { return m_items.empty(); }

        void record(problem& pb, size_t step, RealType t, const vector_type& x) {
            if (step % m_every != 0) return;
            auto& msh = pb.msh;
            auto& as = pb.assembler();
            const auto& di = pb.hho_di;
            const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
            const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
            const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
            const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
            for (auto& it : m_items) {
                it->os << t;
                vector_type xe, xa;
                if (it->e_cell >= 0) xe = as.gather_e_dof_data(size_t(it->e_cell), msh, msh.backend_storage()->surfaces[size_t(it->e_cell)], x);
                if (it->a_cell >= 0) xa = as.gather_a_dof_data(size_t(it->a_cell), msh, msh.backend_storage()->surfaces[size_t(it->a_cell)], x);
                for (auto& f : it->fields) {
                    if (f == "acoustic_pressure") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                        it->os << "," << xa.segment(n_flux, n_scal).dot(cb.eval_functions(it->pt));
                    }
                    else if (f == "acoustic_flux") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto rb = disk::make_scalar_monomial_basis(msh, cell, di.reconstruction_degree());
                        auto dphi = rb.eval_gradients(it->pt);
                        RealType mx = 0.0, my = 0.0;
                        for (long i = 1; i < dphi.rows(); i++) { mx += xa(i - 1) * dphi(i, 0); my += xa(i - 1) * dphi(i, 1); }
                        it->os << "," << mx << "," << my;
                    }
                    else if (f == "elastic_velocity") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                        auto v = disk::eval(vector_type(xe.segment(n_ten, n_vec)), vb.eval_functions(it->pt));
                        it->os << "," << v(0) << "," << v(1);
                    }
                    else if (f == "elastic_stress") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto tb = disk::make_sym_matrix_monomial_basis(msh, cell, di.grad_degree());
                        auto st = disk::eval(vector_type(xe.head(n_ten)), tb.eval_functions(it->pt));
                        it->os << "," << st(0, 0) << "," << st(0, 1) << "," << st(1, 1);
                    }
                    else if (f == "interface_pressure") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->a_cell)];
                        auto fc = faces(msh, cell)[size_t(it->a_face)];
                        auto fb = disk::make_scalar_monomial_basis(msh, fc, di.face_degree());
                        const size_t n_fbs = disk::scalar_basis_size(di.face_degree(), 1);
                        const size_t a_cbs = as.get_a_cell_basis_data();
                        it->os << "," << xa.segment(a_cbs + size_t(it->a_face) * n_fbs, n_fbs).dot(fb.eval_functions(it->pt));
                    }
                    else if (f == "interface_velocity") {
                        auto& cell = msh.backend_storage()->surfaces[size_t(it->e_cell)];
                        auto fc = faces(msh, cell)[size_t(it->e_face)];
                        auto fb = disk::make_vector_monomial_basis(msh, fc, di.face_degree());
                        const size_t n_fbs = disk::vector_basis_size(di.face_degree(), 1, 2);
                        const size_t e_cbs = as.get_e_cell_basis_data();
                        auto v = disk::eval(vector_type(xe.segment(e_cbs + size_t(it->e_face) * n_fbs, n_fbs)), fb.eval_functions(it->pt));
                        it->os << "," << v(0) << "," << v(1);
                    }
                }
                it->os << "\n";
            }
        }
    };

    // Snapshots (output.snapshots): at the time step nearest to each requested time, the file
    // <table>_snapshot_<i>.dat with one line per cell, the acoustic pressure (acoustic cells) or the norm of the elastic
    // velocity (elastic cells) of the cell polynomial at the barycenter, then at the vertices of the cell; the cells
    // (medium, number of vertices, coordinates) once in <table>_cells.dat.
    class snapshots {
        std::vector<size_t> m_steps;
        std::vector<RealType> m_times;
        std::string m_base;
        size_t m_samples = 0;
        bool m_stress = false;      // output.snapshot_elastic = "stress"
    public:
        snapshots(problem& pb, const wave_parameters& prm, const std::string& base, RealType dt) : m_base(base) {
            if (prm.output.snapshots.empty()) return;
            m_samples = prm.output.snapshot_samples;
            m_stress = prm.output.snapshot_elastic == "stress";
            for (auto t : prm.output.snapshots) {
                m_steps.push_back(size_t(std::llround(t / dt)));
                m_times.push_back(t);
            }
            std::ofstream os(base + "_cells.dat");
            os << "# cells: medium (0 elastic, 1 acoustic), number of vertices, x_1 y_1 ... x_n y_n" << std::endl;
            os << std::setprecision(12);
            for (auto& cell : pb.msh) {
                auto ci = pb.msh.lookup(cell);
                auto pts = points(pb.msh, cell);
                os << (pb.dd.a_material.count(ci) ? 1 : 0) << " " << pts.size();
                for (auto& q : pts) os << " " << q.x() << " " << q.y();
                os << "\n";
            }
        }

        bool empty() const { return m_steps.empty(); }

        void maybe_write(problem& pb, size_t step, RealType t, const vector_type& x) {
            m_state = &x;
            for (size_t i = 0; i < m_steps.size(); i++) {
                if (m_steps[i] != step) continue;
                auto& msh = pb.msh;
                auto& as = pb.assembler();
                const auto& di = pb.hho_di;
                const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
                const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
                const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
                const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
                std::ofstream os(m_base + "_snapshot_" + std::to_string(i) + ".dat");
                os << "# t = " << std::setprecision(12) << t << " (requested " << m_times[i] << "): acoustic pressure "
                   << "(acoustic cells), |elastic velocity| (elastic cells) at the barycenter, then at the vertices" << std::endl;
                os << std::setprecision(8) << std::scientific;
                for (auto& cell : msh) {
                    auto ci = msh.lookup(cell);
                    auto pts = points(msh, cell);
                    std::vector<point_type> at(1, barycenter(msh, cell));
                    at.insert(at.end(), pts.begin(), pts.end());
                    if (pb.dd.a_material.count(ci)) {
                        vector_type xa = as.gather_a_dof_data(ci, msh, cell, x);
                        vector_type xp = xa.segment(n_flux, n_scal);
                        auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                        for (size_t j = 0; j < at.size(); j++)
                            os << (j ? " " : "") << xp.dot(cb.eval_functions(at[j]));
                    }
                    else {
                        vector_type xe = as.gather_e_dof_data(ci, msh, cell, x);
                        vector_type xv = xe.segment(n_ten, n_vec);
                        auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                        for (size_t j = 0; j < at.size(); j++) {
                            auto v = disk::eval(xv, vb.eval_functions(at[j]));
                            os << (j ? " " : "") << std::sqrt(v(0) * v(0) + v(1) * v(1));
                        }
                    }
                    os << "\n";
                }
                if (m_samples > 0) write_samples(pb, i);
            }
        }

        // the same fields on an m x m grid of the bounding box of each cell (centres of the sub-squares, rows from bottom
        // to top, each row from left to right): one line per cell, xmin xmax ymin ymax then the m^2 values
        void write_samples(problem& pb, size_t i) {
            auto& msh = pb.msh;
            auto& as = pb.assembler();
            const auto& di = pb.hho_di;
            const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
            const size_t n_vec  = disk::vector_basis_size(di.cell_degree(), 2, 2);
            const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
            const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
            const size_t m = m_samples;
            std::ofstream os(m_base + "_snapshot_" + std::to_string(i) + "_samples.dat");
            os << "# t requested " << m_times[i] << ": per cell xmin xmax ymin ymax, then the " << m << " x " << m
               << " values (acoustic pressure or " << (m_stress ? "mean stress -tr(sigma)/2" : "|elastic velocity|")
               << ") at the centres of the sub-squares, rows bottom to top" << std::endl;
            os << std::setprecision(7) << std::scientific;
            // stress mode: shear stress sqrt(((sxx - syy)/2)^2 + sxy^2) of the elastic cells (0 in the acoustic ones)
            std::ofstream osh;
            if (m_stress) {
                osh.open(m_base + "_snapshot_" + std::to_string(i) + "_shear.dat");
                osh << "# t requested " << m_times[i] << ": shear stress sqrt(((sxx - syy)/2)^2 + sxy^2) (elastic cells, 0 in "
                    << "the acoustic cells), format of _samples.dat" << std::endl << std::setprecision(7) << std::scientific;
            }
            for (auto& cell : msh) {
                auto ci = msh.lookup(cell);
                auto pts = points(msh, cell);
                RealType x0 = pts[0].x(), x1 = x0, y0 = pts[0].y(), y1 = y0;
                for (auto& q : pts) { x0 = std::min(x0, q.x()); x1 = std::max(x1, q.x()); y0 = std::min(y0, q.y()); y1 = std::max(y1, q.y()); }
                os << x0 << " " << x1 << " " << y0 << " " << y1;
                if (m_stress) osh << x0 << " " << x1 << " " << y0 << " " << y1;
                const bool acoustic = pb.dd.a_material.count(ci) > 0;
                vector_type xp, xv, xs;
                if (acoustic) xp = as.gather_a_dof_data(ci, msh, cell, *m_state).segment(n_flux, n_scal);
                else {
                    vector_type xe = as.gather_e_dof_data(ci, msh, cell, *m_state);
                    xv = xe.segment(n_ten, n_vec);
                    xs = xe.head(n_ten);
                }
                auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                auto vb = disk::make_vector_monomial_basis(msh, cell, di.cell_degree());
                auto tb = disk::make_sym_matrix_monomial_basis(msh, cell, di.grad_degree());
                for (size_t r = 0; r < m; r++)
                    for (size_t c = 0; c < m; c++) {
                        const point_type at((x0 * (2 * m - 2 * c - 1) + x1 * (2 * c + 1)) / (2 * m),
                                            (y0 * (2 * m - 2 * r - 1) + y1 * (2 * r + 1)) / (2 * m));
                        if (acoustic) {
                            os << " " << xp.dot(cb.eval_functions(at));
                            if (m_stress) osh << " 0";
                        }
                        else if (m_stress) {
                            auto st = disk::eval(xs, tb.eval_functions(at));
                            os << " " << -0.5 * (st(0, 0) + st(1, 1));
                            const RealType d = 0.5 * (st(0, 0) - st(1, 1));
                            osh << " " << std::sqrt(d * d + st(0, 1) * st(0, 1));
                        }
                        else {
                            auto v = disk::eval(xv, vb.eval_functions(at));
                            os << " " << std::sqrt(v(0) * v(0) + v(1) * v(1));
                        }
                    }
                os << "\n";
                if (m_stress) osh << "\n";
            }
        }
        const vector_type* m_state = nullptr;
    };

    // Peak map (output.peak_samples = m > 0): maximum over the time steps (every output.sensor_every steps, and the
    // initial state) of |p| in the acoustic cells and of the mean stress |tr(sigma)|/2 in the elastic cells, at the
    // centres of an m x m grid of the bounding box of each cell; <table>_peak_samples.dat, format of the snapshot
    // samples. The values at the sample points are linear in the cell unknowns: one small matrix per cell, built once.
    class peak_map {
        struct cell_data { long off; int n; Matrix<RealType, Dynamic, Dynamic> phi; RealType x0, x1, y0, y1; };
        std::vector<cell_data> m_cells;
        std::vector<RealType> m_peak;
        size_t m = 0, m_every = 1;
        std::string m_base;
    public:
        peak_map(problem& pb, const wave_parameters& prm, const std::string& base, const vector_type& x)
            : m(prm.output.peak_samples), m_every(std::max<size_t>(1, prm.output.sensor_every)), m_base(base) {
            auto& msh = pb.msh;
            auto& as = pb.assembler();
            const auto& di = pb.hho_di;
            const size_t n_ten  = disk::sym_matrix_basis_size(di.grad_degree(), 2, 2);
            const size_t n_flux = disk::scalar_basis_size(di.reconstruction_degree(), 2) - 1;
            const size_t n_scal = disk::scalar_basis_size(di.cell_degree(), 2);
            // global position of the first cell unknown: gathered from the vector of the indices
            vector_type idx(x.size());
            for (Eigen::Index k = 0; k < idx.size(); k++) idx[k] = RealType(k);
            m_cells.resize(msh.cells_size());
            for (auto& cell : msh) {
                auto ci = msh.lookup(cell);
                cell_data& cd = m_cells[ci];
                auto pts = points(msh, cell);
                cd.x0 = cd.x1 = pts[0].x(); cd.y0 = cd.y1 = pts[0].y();
                for (auto& q : pts) { cd.x0 = std::min(cd.x0, q.x()); cd.x1 = std::max(cd.x1, q.x()); cd.y0 = std::min(cd.y0, q.y()); cd.y1 = std::max(cd.y1, q.y()); }
                const bool acoustic = pb.dd.a_material.count(ci) > 0;
                const long first = acoustic ? long(as.gather_a_dof_data(ci, msh, cell, idx)(0))
                                            : long(as.gather_e_dof_data(ci, msh, cell, idx)(0));
                cd.off = acoustic ? first + long(n_flux) : first;
                cd.n = int(acoustic ? n_scal : n_ten);
                cd.phi.resize(m * m, cd.n);
                auto cb = disk::make_scalar_monomial_basis(msh, cell, di.cell_degree());
                auto tb = disk::make_sym_matrix_monomial_basis(msh, cell, di.grad_degree());
                for (size_t r = 0; r < m; r++)
                    for (size_t c = 0; c < m; c++) {
                        const point_type at((cd.x0 * (2 * m - 2 * c - 1) + cd.x1 * (2 * c + 1)) / (2 * m),
                                            (cd.y0 * (2 * m - 2 * r - 1) + cd.y1 * (2 * r + 1)) / (2 * m));
                        if (acoustic) cd.phi.row(r * m + c) = cb.eval_functions(at).transpose();
                        else {
                            auto T = tb.eval_functions(at);
                            for (int j = 0; j < cd.n; j++) cd.phi(r * m + c, j) = -0.5 * (T[j](0, 0) + T[j](1, 1));
                        }
                    }
            }
            m_peak.assign(m_cells.size() * m * m, RealType(0));
            record(0, x);
        }

        void record(size_t step, const vector_type& x) {
            if (step % m_every) return;
            const long nc = long(m_cells.size());
            #pragma omp parallel for schedule(static)
            for (long k = 0; k < nc; k++) {
                const cell_data& cd = m_cells[k];
                const vector_type v = cd.phi * x.segment(cd.off, cd.n);
                RealType* pk = m_peak.data() + size_t(k) * m * m;
                for (size_t q = 0; q < m * m; q++) pk[q] = std::max(pk[q], std::abs(v[q]));
            }
        }

        void write() const {
            std::ofstream os(m_base + "_peak_samples.dat");
            os << "# max over time of |p| (acoustic cells) and |tr(sigma)|/2 (elastic cells): per cell xmin xmax ymin ymax, "
               << "then the " << m << " x " << m << " values at the centres of the sub-squares, rows bottom to top" << std::endl;
            os << std::setprecision(7) << std::scientific;
            for (size_t k = 0; k < m_cells.size(); k++) {
                const cell_data& cd = m_cells[k];
                os << cd.x0 << " " << cd.x1 << " " << cd.y0 << " " << cd.y1;
                for (size_t q = 0; q < m * m; q++) os << " " << m_peak[k * m * m + q];
                os << "\n";
            }
        }
    };

}

#endif /* first_order_observers_hpp */
