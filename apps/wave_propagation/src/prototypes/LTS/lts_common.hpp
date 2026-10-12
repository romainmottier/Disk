//  Created by Romain Mottier
//
// Scaffolding shared by the acoustic LTS drivers (prototypes/LTS,
// prototypes/LTS/mlts_2026, prototypes/LTS_ERK_SDIRK). Every driver solves
// the first-order acoustic system (rho = vp = 1) on a poly2d mesh, with
// Dirichlet data on every boundary face and the coupled four-fields HHO
// assembler used acoustic-only; only the time integrator changes from one
// driver to the next (classical two-level, multilevel recursion,
// LI-SDIRK, ...).
//
// An exact solution (analytic_solution, slow_sine_solution,
// lshape_corner_solution) provides the Dirichlet/initial data, the
// right-hand side at any time and the fields used for the L2 errors.
#pragma once

#include <filesystem>
#include <functional>
#include <memory>
#include <sys/resource.h>

namespace lts {

using RealType       = double;
using mesh_type      = disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>;
using point_type     = typename mesh_type::point_type;
using vector_type    = Matrix<RealType, Dynamic, 1>;
using assembler_type = elastoacoustic_four_fields_assembler<mesh_type>;
using erk_type       = erk_coupling_hho_scheme<RealType>;

inline const std::string mesh_root = "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/";

inline const auto null_fun = [](const point_type&) -> disk::static_vector<double, 2> {
    return disk::static_vector<double, 2>{0, 0};
};
inline const auto null_flux_fun = [](const point_type&) -> disk::static_matrix<double, 2, 2> {
    return disk::static_matrix<double, 2, 2>::Zero(2, 2);
};
inline const auto null_s_fun = [](const point_type&) -> RealType { return 0.0; };

// Process CPU time (user + system) via getrusage: immune to wall-clock
// contention from unrelated processes.
inline double cpu_seconds() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec
         + (double)ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
}

inline RealType env_real(const char* name, RealType value) {
    if (const char* env = std::getenv(name)) value = std::atof(env);
    return value;
}

inline int env_int(const char* name, int value) {
    if (const char* env = std::getenv(name)) value = std::atoi(env);
    return value;
}

// Mesh levels N = n_min..n_max; MLTS_N_MAX (and MLTS_N_MIN when
// read_n_min) override the bounds.
inline std::vector<int> levels(int n_min, int n_max, bool read_n_min = false) {
    n_max = env_int("MLTS_N_MAX", n_max);
    if (read_n_min) n_min = env_int("MLTS_N_MIN", n_min);
    std::vector<int> level_N;
    for (int n = n_min; n <= n_max; ++n) level_N.push_back(n);
    return level_N;
}

inline disk::hho_degree_info degree_info(const simulation_data& sim_data) {
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) cell_k_degree++;
    return disk::hho_degree_info(cell_k_degree, sim_data.m_k_degree);
}

inline void load_mesh(mesh_type& msh, const std::string& mesh_fname) {
    polygon_2d_mesh_reader<RealType> mesh_builder;
    mesh_builder.set_poly_mesh_file(mesh_fname);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
}

inline void h_range(mesh_type& msh, RealType& h_max, RealType& h_min) {
    h_max = 0.0; h_min = 1.0e10;
    for (auto & cell : msh) {
        RealType h_l = diameter(msh, cell);
        if (h_l > h_max) h_max = h_l;
        if (h_l < h_min) h_min = h_l;
    }
}

// L2 error of the pressure: mass-matrix projection of the exact field
// against the cell pressure dofs (degree = cell_degree).
template<typename Function>
RealType pressure_l2_error(mesh_type& msh, const disk::hho_degree_info& hho_di, const vector_type& x_dof, const Function& p_exact) {
    size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
    size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
    size_t cell_dof   = n_scal_dof + n_vec_dof;
    RealType l2_error_sq = 0.0;
    size_t cell_i = 0;
    for (auto & cell : msh) {
        Matrix<RealType, Dynamic, 1> p_dof = x_dof.block(cell_i * cell_dof + n_vec_dof, 0, n_scal_dof, 1);
        auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
        Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
        Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, p_exact);
        Matrix<RealType, Dynamic, 1> proj = mass.llt().solve(rhs);
        Matrix<RealType, Dynamic, 1> diff = proj - p_dof;
        l2_error_sq += diff.dot(mass * diff);
        cell_i++;
    }
    return std::sqrt(l2_error_sq);
}

// L2 error of the velocity (the dual variable, at the start of each cell's
// dof block): its n_vec_dof coefficients are the gradients of a scalar
// basis of degree reconstruction_degree (the constant is skipped), so it
// is integrated by quadrature as in postprocessor.hpp's four-fields error
// routines, not projected.
template<typename Function>
RealType velocity_l2_error(mesh_type& msh, const disk::hho_degree_info& hho_di, const vector_type& x_dof, const Function& q_exact) {
    size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
    size_t n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
    size_t cell_dof   = n_scal_dof + n_vec_dof;
    RealType v_l2_error_sq = 0.0;
    size_t cell_i = 0;
    for (auto & cell : msh) {
        Matrix<RealType, Dynamic, 1> v_dof = x_dof.block(cell_i * cell_dof, 0, n_vec_dof, 1);
        auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
        auto int_rule = integrate(msh, cell, 2 * (hho_di.cell_degree() + 1));
        for (auto & qp : int_rule) {
            auto t_dphi = rec_basis.eval_gradients(qp.point());
            Matrix<RealType, 1, 2> vh = Matrix<RealType, 1, 2>::Zero();
            for (size_t i = 1; i < t_dphi.rows(); i++)
                vh = vh + v_dof(i - 1) * t_dphi.block(i, 0, 1, 2);
            disk::static_vector<RealType, 2> v_exact = q_exact(qp.point());
            Matrix<RealType, 1, 2> v_exact_row; v_exact_row << v_exact(0), v_exact(1);
            v_l2_error_sq += qp.weight() * (v_exact_row - vh).squaredNorm();
        }
        cell_i++;
    }
    return std::sqrt(v_l2_error_sq);
}

// ---------------------------------------------------------------------------
// Exact solutions
// ---------------------------------------------------------------------------

// Manufactured solution of scal_vec_analytic_functions (acoustic part):
// Dirichlet data and source are refreshed at every evaluation time. The
// evaluators capture `t` by reference, so the object must not move once
// the problem has been set up.
class analytic_solution {
    scal_vec_analytic_functions m_functions;
public:
    RealType t = 0.0;

    explicit analytic_solution(scal_vec_analytic_functions::EFunctionType type) {
        m_functions.set_function_type(type);
    }

    auto dirichlet0() { return m_functions.Evaluate_s_u(t); }
    auto p0()         { return m_functions.Evaluate_s_v(t); }
    auto q0()         { return m_functions.Evaluate_s_q(t); }
    void prepare(assembler_type&, mesh_type&) {}

    // Sets the pressure Dirichlet data at t_abs and returns the source.
    auto set_time(assembler_type& assembler, RealType t_abs) {
        t = t_abs;
        auto s_v_fun_t = m_functions.Evaluate_s_v(t);
        auto s_f_fun_t = m_functions.Evaluate_s_f(t);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
        return s_f_fun_t;
    }

    vector_type rhs(assembler_type& assembler, mesh_type& msh, RealType t_abs) {
        auto s_f_fun_t = set_time(assembler, t_abs);
        assembler.assemble_rhs(msh, null_fun, s_f_fun_t, true);
        return assembler.RHS;
    }

    auto p_exact(RealType tf) { t = tf; return m_functions.Evaluate_s_v(t); }
    auto q_exact(RealType tf) { t = tf; return m_functions.Evaluate_s_q(t); }

    scal_vec_analytic_functions& functions() { return m_functions; }
};

// Slow sine: s_u = x^2*sin(pi*x)*sin(pi*y)*sin(omega*t) -- SAME spatial
// ansatz as scal_vec_analytic_functions' EFunctionNonPolynomial acoustic
// side, generalized from the fixed omega=sqrt(2)*pi to a SLOWER omega=pi/2
// (~4x slower), for the acoustic wave equation u_tt = Delta(u) + f (vp=1).
// Derived by direct differentiation, cross-validated two ways: (1) sympy
// symbolic + numeric finite-difference PDE-residual check (residual ~3.5e-8
// at h=1e-5, consistent with pure FD truncation, i.e. exact); (2) at
// omega=sqrt(2)*pi this formula EXACTLY reproduces EFunctionNonPolynomial's
// own s_v/s_q formulas term-for-term:
//   s_v = d(s_u)/dt = omega*x^2*sin(pi x)*sin(pi y)*cos(omega t)
//   s_q = grad(s_u) = ( x*(pi*x*cos(pi x)+2*sin(pi x))*sin(pi y)*sin(omega t),
//                       pi*x^2*sin(pi x)*cos(pi y)*sin(omega t) )
//   s_f = u_tt - Delta(u)
//       = sin(omega t)*sin(pi y)*[ (2*pi^2-omega^2)*x^2*sin(pi x)
//                                   - 2*sin(pi x) - 4*pi*x*cos(pi x) ]
// The evaluators capture the time by value.
class slow_sine_solution {
    static constexpr RealType omega = M_PI/2.0;

    static auto s_u(RealType tt) {
        return [tt](const point_type& pt) -> RealType {
            RealType x = pt.x(), y = pt.y();
            return x*x*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(omega*tt);
        };
    }
    static auto s_v(RealType tt) {
        return [tt](const point_type& pt) -> RealType {
            RealType x = pt.x(), y = pt.y();
            return omega*x*x*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(omega*tt);
        };
    }
    static auto s_f(RealType tt) {
        return [tt](const point_type& pt) -> RealType {
            RealType x = pt.x(), y = pt.y();
            return std::sin(omega*tt)*std::sin(M_PI*y)*(
                (2.0*M_PI*M_PI - omega*omega)*x*x*std::sin(M_PI*x)
                - 2.0*std::sin(M_PI*x) - 4.0*M_PI*x*std::cos(M_PI*x));
        };
    }
    static auto s_q(RealType tt) {
        return [tt](const point_type& pt) -> disk::static_vector<RealType,2> {
            RealType x = pt.x(), y = pt.y();
            RealType qx = x*(M_PI*x*std::cos(M_PI*x) + 2.0*std::sin(M_PI*x))*std::sin(M_PI*y)*std::sin(omega*tt);
            RealType qy = M_PI*x*x*std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(omega*tt);
            return disk::static_vector<RealType,2>{qx, qy};
        };
    }

public:
    RealType t = 0.0;

    auto dirichlet0() { return s_u(t); }
    auto p0()         { return s_v(t); }
    auto q0()         { return s_q(t); }
    void prepare(assembler_type&, mesh_type&) {}

    vector_type rhs(assembler_type& assembler, mesh_type& msh, RealType t_abs) {
        t = t_abs;
        auto s_v_fun_t = s_v(t);
        auto s_f_fun_t = s_f(t);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun_t, 0);
        assembler.assemble_rhs(msh, null_fun, s_f_fun_t, true);
        return assembler.RHS;
    }

    auto p_exact(RealType tf) { t = tf; return s_v(t); }
    auto q_exact(RealType tf) { t = tf; return s_q(t); }
};

// L-shape re-entrant corner benchmark: p = -w sin(w t) phi, v = cos(w t)
// grad(phi), phi = r^nu sin(nu theta) around the corner (0.5, 0.5),
// nu = 2/3. Homogeneous Dirichlet and zero initial pressure; the
// right-hand side is a time-modulated combination of two unit vectors
// (boundary lift of phi, source phi) assembled once in prepare().
class lshape_corner_solution {
    static constexpr RealType corner_omega = 1.5 * M_PI;
    static constexpr RealType nu           = M_PI / corner_omega;
    static constexpr RealType w            = 2.0;

    vector_type m_rhs_bc_unit, m_rhs_src_unit;

public:
    static std::pair<RealType, RealType> polar(const point_type& pt) {
        RealType dx = pt.x() - 0.5, dy = pt.y() - 0.5;
        RealType r = std::sqrt(dx * dx + dy * dy);
        RealType theta = std::atan2(dy, dx) - 0.5 * M_PI;
        if (theta < 0.0) theta += 2.0 * M_PI;
        return {r, theta};
    }
    static RealType phi(RealType r, RealType theta) {
        if (r < 1.0e-12) return 0.0;
        return std::pow(r, nu) * std::sin(nu * theta);
    }
    static disk::static_vector<double, 2> grad_phi(RealType r, RealType theta) {
        disk::static_vector<double, 2> g{0.0, 0.0};
        if (r < 1.0e-10) return g;
        RealType A = nu * std::pow(r, nu - 1.0) * std::sin(nu * theta);
        RealType B = nu * std::pow(r, nu - 1.0) * std::cos(nu * theta);
        RealType c = std::cos(theta), s = std::sin(theta);
        g[0] = -A * s - B * c;
        g[1] =  A * c - B * s;
        return g;
    }
    static auto p_exact_at(RealType t) {
        return [t](const point_type& pt) -> RealType {
            auto [r, theta] = polar(pt);
            return -w * std::sin(w * t) * phi(r, theta);
        };
    }
    static auto v_exact_at(RealType t) {
        return [t](const point_type& pt) -> disk::static_vector<double, 2> {
            auto [r, theta] = polar(pt);
            return std::cos(w * t) * grad_phi(r, theta);
        };
    }

    auto dirichlet0() { return null_s_fun; }
    auto p0()         { return null_s_fun; }
    auto q0()         { return v_exact_at(0.0); }

    void prepare(assembler_type& assembler, mesh_type& msh) {
        auto phi_field_fun = [](const point_type& pt) -> RealType {
            auto [r, theta] = polar(pt);
            return phi(r, theta);
        };
        assembler.get_a_bc_conditions().updateDirichletFunction(phi_field_fun, 0);
        assembler.assemble_rhs(msh, null_fun, null_s_fun, true);
        m_rhs_bc_unit = assembler.RHS;
        assembler.get_a_bc_conditions().updateDirichletFunction(null_s_fun, 0);
        assembler.assemble_rhs(msh, null_fun, phi_field_fun, true);
        m_rhs_src_unit = assembler.RHS;
    }

    vector_type rhs(assembler_type&, mesh_type&, RealType t_abs) {
        RealType g = -w * std::sin(w * t_abs);
        RealType f = -w * w * std::cos(w * t_abs);
        return g * m_rhs_bc_unit + f * m_rhs_src_unit;
    }

    auto p_exact(RealType tf) { return p_exact_at(tf); }
    auto q_exact(RealType tf) { return v_exact_at(tf); }
};

// ---------------------------------------------------------------------------
// Problem: mesh + acoustic HHO discretization + explicit operator
// ---------------------------------------------------------------------------

// The constructor reads the mesh and its size range; setup() assembles
// everything (materials, boundary conditions, initial projection, global
// operators, erk_coupling_hho_scheme). Non-copyable: the boundary data
// stored in the assembler refers to `sol`.
template<typename Solution>
class acoustic_problem {
public:
    Solution  sol;
    mesh_type msh;
    RealType  h_max, h_min;
    int       p_global;
    disk::hho_degree_info hho_di;
    std::map<size_t, acoustic_material_data<RealType>> a_material;
    std::map<size_t, elastic_material_data<RealType>>  e_material;
    std::set<size_t> acoustic_bc_face_indexes, interface_face_indexes;
    vector_type x_dof;
    size_t n_scal_dof, n_vec_dof, cell_dof, a_fbs;

    acoustic_problem(const std::string& mesh_fname, Solution s = Solution()) : sol(std::move(s)) {
        load_mesh(msh, mesh_fname);
        h_range(msh, h_max, h_min);
        p_global = static_cast<int>(std::round(h_max / h_min));
        if (p_global < 1) p_global = 1;
    }
    acoustic_problem(const acoustic_problem&) = delete;
    acoustic_problem& operator=(const acoustic_problem&) = delete;

    void setup(const disk::hho_degree_info& di, bool scaled_stabilization_Q) {
        hho_di = di;
        typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
        typedef disk::BoundaryConditions<mesh_type, true>  a_boundary_type;

        for (auto & cell : msh) {
            auto cell_ind = msh.lookup(cell);
            a_material.insert(std::make_pair(cell_ind, acoustic_material_data<RealType>(1.0, 1.0)));
        }

        size_t bc_elastic_id  = 0;
        size_t bc_acoustic_id = 1;
        for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
            auto face  = *face_it;
            auto fc_id = msh.lookup(face);
            disk::boundary_descriptor bi{bc_acoustic_id, true};
            msh.backend_storage()->boundary_info.at(fc_id) = bi;
            acoustic_bc_face_indexes.insert(fc_id);
        }

        auto s_u_fun0    = sol.dirichlet0();
        auto s_v_fun0    = sol.p0();
        auto s_flux_fun0 = sol.q0();

        e_boundary_type e_bnd(msh);
        a_boundary_type a_bnd(msh);
        e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id,  null_fun);
        a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun0);

        std::map<size_t, std::pair<size_t, size_t>> interface_cell_pair_indexes;
        m_assembler = std::make_unique<assembler_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
        assembler_type& assembler = *m_assembler;
        assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
        assembler.set_coupling_stabilization();
        if (scaled_stabilization_Q) assembler.set_scaled_stabilization();

        assembler.assemble_mass(msh);
        assembler.assemble_coupling_terms(msh);

        assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun0, s_flux_fun0);
        assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun0);

        assembler.assemble(msh, null_fun, s_v_fun0, true);
        assembler.LHS += assembler.COUPLING;

        std::set<size_t> elastic_internal_faces, acoustic_internal_faces;
        m_erk = std::make_unique<erk_type>(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING,
            assembler.get_e_n_cells_dof(), assembler.get_a_n_cells_dof(), assembler.get_e_face_dof(), assembler.get_a_face_dof());
        erk_type& erk_an = *m_erk;
        erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(),
            assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
        erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(),
            assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(),
            assembler.get_e_compress(), assembler.get_a_compress(),
            elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
        erk_an.refresh_faces_unknowns(x_dof);

        n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), mesh_type::dimension);
        n_vec_dof  = disk::scalar_basis_size(hho_di.reconstruction_degree(), mesh_type::dimension) - 1;
        cell_dof   = n_scal_dof + n_vec_dof;
        a_fbs      = disk::scalar_basis_size(hho_di.face_degree(), mesh_type::dimension - 1);

        sol.prepare(assembler, msh);
    }

    assembler_type& assembler() { return *m_assembler; }
    erk_type&       erk()       { return *m_erk; }

    vector_type rhs(RealType t_abs) { return sol.rhs(*m_assembler, msh, t_abs); }

    RealType pressure_error(const vector_type& x, RealType tf) {
        return pressure_l2_error(msh, hho_di, x, sol.p_exact(tf));
    }
    RealType velocity_error(const vector_type& x, RealType tf) {
        return velocity_l2_error(msh, hho_di, x, sol.q_exact(tf));
    }

private:
    std::unique_ptr<assembler_type> m_assembler;
    std::unique_ptr<erk_type>       m_erk;
};

// ---------------------------------------------------------------------------
// Multi-level bookkeeping
// ---------------------------------------------------------------------------

// Natural split: thresholds h_max/4, h_max/16, ... while above 2 h_min.
inline std::vector<RealType> natural_thresholds(RealType h_max, RealType h_min) {
    std::vector<RealType> thresholds;
    RealType t_th = h_max / 4.0;
    while (t_th > 2.0 * h_min) { thresholds.push_back(t_th); t_th /= 4.0; }
    return thresholds;
}

// Forced number of levels: thresholds h_max/2^k, k = 1..L-1.
inline std::vector<RealType> forced_thresholds(RealType h_max, int L) {
    std::vector<RealType> thresholds;
    for (int k = 1; k < L; ++k)
        thresholds.push_back(h_max / std::pow(2.0, k));
    return thresholds;
}

// MLTS_FORCE_L when set, natural split otherwise.
inline std::vector<RealType> env_thresholds(RealType h_max, RealType h_min) {
    if (const char* env_l = std::getenv("MLTS_FORCE_L"))
        return forced_thresholds(h_max, std::atoi(env_l));
    return natural_thresholds(h_max, h_min);
}

// Band of a cell: first i with h >= thresholds[i] (h > thresholds[i] when
// strict), L-1 otherwise. A face takes the finest band among its cells
// (assemble_P's union rule generalized to L bands); own_c/own_f are the
// cell/face dofs owned by each band (Dirichlet faces excluded).
struct bands {
    int L;
    std::vector<int> cell_band, face_band;
    std::vector<std::vector<int>> own_c, own_f;
};

template<typename Problem>
bands make_bands(Problem& pb, const std::vector<RealType>& thresholds, bool strict = false) {
    mesh_type& msh = pb.msh;
    bands bd;
    const int L = (int)thresholds.size() + 1;
    bd.L = L;
    bd.cell_band.resize(msh.cells_size());
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            RealType h = diameter(msh, cell);
            int band = L - 1;
            if (strict) {
                for (int i = 0; i < (int)thresholds.size(); ++i) if (h > thresholds[i]) { band = i; break; }
            } else {
                for (int i = 0; i < (int)thresholds.size(); ++i) if (h >= thresholds[i]) { band = i; break; }
            }
            bd.cell_band[ci] = band;
            ci++;
        }
    }
    bd.face_band.assign(msh.faces_size(), 0);
    {
        size_t ci = 0;
        for (auto & cell : msh) {
            for (auto & face : faces(msh, cell)) {
                auto fid = msh.lookup(face);
                if (bd.cell_band[ci] > bd.face_band[fid]) bd.face_band[fid] = bd.cell_band[ci];
            }
            ci++;
        }
    }

    std::vector<size_t> a_compress = pb.assembler().get_a_compress();
    bd.own_c.resize(L); bd.own_f.resize(L);
    for (size_t i = 0; i < msh.cells_size(); ++i) {
        int b = bd.cell_band[i];
        for (size_t d = 0; d < pb.cell_dof; ++d) bd.own_c[b].push_back((int)(i*pb.cell_dof + d));
    }
    for (size_t j = 0; j < msh.faces_size(); ++j) {
        if (pb.acoustic_bc_face_indexes.count(j)) continue;
        int b = bd.face_band[j];
        size_t compressed = a_compress.at(j);
        for (size_t d = 0; d < pb.a_fbs; ++d) bd.own_f[b].push_back((int)(compressed*pb.a_fbs + d));
    }
    return bd;
}

inline void print_bands(const bands& bd) {
    for (int b = 0; b < bd.L; ++b) {
        std::cout << bold << cyan << "      band" << b << ": own_c=" << bd.own_c[b].size()
                  << " own_f=" << bd.own_f[b].size() << reset << std::endl;
    }
}

// Diagonal projector (n_dof x n_dof): 1 at the given cell/face dofs.
inline SparseMatrix<RealType> band_projector(const std::vector<int>& own_c, const std::vector<int>& own_f, size_t n_c_dof, size_t n_dof) {
    std::vector<Eigen::Triplet<RealType>> trip;
    for (int g : own_c) trip.emplace_back(g, g, 1.0);
    for (int g : own_f) trip.emplace_back((int)n_c_dof+g, (int)n_c_dof+g, 1.0);
    SparseMatrix<RealType> P(n_dof, n_dof);
    P.setFromTriplets(trip.begin(), trip.end());
    return P;
}

// Cell dofs followed by the face dofs shifted by n_c_dof, sorted.
inline std::vector<int> merged_dofs(const std::vector<int>& c, const std::vector<int>& f, size_t n_c_dof) {
    std::vector<int> dofs = c;
    for (int g : f) dofs.push_back((int)n_c_dof + g);
    std::sort(dofs.begin(), dofs.end());
    return dofs;
}

// Pell_c/Pell_f: own + every finer band (cumulative, nested P_l), each
// non-terminal level (every level when include_terminal) dilated by
// `overlap_rings` cell/face rings through Kfc -- the Diaz-Grote /
// Almquist-Mehlin (Section 3.2, Table 2: 3 rings for RK4) boundary margin
// RK-LTS stability needs beyond the strict level boundary, reaching a few
// extra dofs on the COARSER side.
inline void pell_sets(erk_type& erk_an, const bands& bd, int overlap_rings,
                      std::vector<std::vector<int>>& Pell_c, std::vector<std::vector<int>>& Pell_f,
                      bool include_terminal = false) {
    const int L = bd.L;
    Pell_c.assign(L, {}); Pell_f.assign(L, {});
    Pell_c[L-1] = bd.own_c[L-1]; Pell_f[L-1] = bd.own_f[L-1];
    for (int b = L-2; b >= 0; --b) {
        Pell_c[b] = Pell_c[b+1]; Pell_c[b].insert(Pell_c[b].end(), bd.own_c[b].begin(), bd.own_c[b].end());
        Pell_f[b] = Pell_f[b+1]; Pell_f[b].insert(Pell_f[b].end(), bd.own_f[b].begin(), bd.own_f[b].end());
    }

    if (overlap_rings > 0) {
        auto& Kfc_g = erk_an.Kfc();
        const int n_dilated = include_terminal ? L : L - 1;
        for (int b = 0; b < n_dilated; ++b) {
            std::set<int> set_c(Pell_c[b].begin(), Pell_c[b].end());
            std::set<int> set_f(Pell_f[b].begin(), Pell_f[b].end());
            for (int ring = 0; ring < overlap_rings; ++ring) {
                std::set<int> new_c, new_f;
                for (int k = 0; k < Kfc_g.outerSize(); ++k) {
                    for (SparseMatrix<RealType>::InnerIterator it(Kfc_g, k); it; ++it) {
                        int face = (int)it.row(), cell = (int)it.col();
                        bool cell_in = set_c.count(cell) > 0;
                        bool face_in = set_f.count(face) > 0;
                        if (cell_in && !face_in) new_f.insert(face);
                        if (face_in && !cell_in) new_c.insert(cell);
                    }
                }
                if (new_c.empty() && new_f.empty()) break;
                set_c.insert(new_c.begin(), new_c.end());
                set_f.insert(new_f.begin(), new_f.end());
            }
            Pell_c[b].assign(set_c.begin(), set_c.end());
            Pell_f[b].assign(set_f.begin(), set_f.end());
        }
    }
}

// Per-level time steps dt_level[i] = cfl * scale[i], scale = (h_max,
// thresholds..., h_min), and sub-step counts p_level[i] =
// round(dt_level[i]/dt_level[i+1]).
inline void level_steps(RealType cfl_factor, RealType h_max, RealType h_min, const std::vector<RealType>& thresholds,
                        std::vector<RealType>& dt_level, std::vector<int>& p_level) {
    const int L = (int)thresholds.size() + 1;
    std::vector<RealType> scale(L);
    dt_level.assign(L, 0.0);
    scale[0] = h_max;
    for (int i = 1; i < L - 1; ++i) scale[i] = thresholds[i-1];
    scale[L-1] = h_min;
    for (int i = 0; i < L; ++i) dt_level[i] = cfl_factor * scale[i];
    p_level.assign(L, 1);
    for (int i = 0; i < L - 1; ++i)
        p_level[i] = std::max(1, (int)std::round(dt_level[i] / dt_level[i+1]));
}

// Number of macro steps over [0, tf] for a target step dt: dt is replaced
// by tf/nt.
inline size_t macro_steps(RealType tf, RealType& dt) {
    size_t nt = std::max<size_t>(1, static_cast<size_t>(std::round(tf / dt)));
    dt = tf / nt;
    return nt;
}

// Taylor coefficients (w0..w3) of a cubic polynomial shifted by a.
inline std::vector<vector_type> taylor_shift(const std::vector<vector_type>& wsrc, RealType a) {
    std::vector<vector_type> out(4);
    RealType a2 = a*a, a3 = a2*a;
    out[0] = wsrc[0] + a*wsrc[1] + (a2/2.0)*wsrc[2] + (a3/6.0)*wsrc[3];
    out[1] = wsrc[1] + a*wsrc[2] + (a2/2.0)*wsrc[3];
    out[2] = wsrc[2] + a*wsrc[3];
    out[3] = wsrc[3];
    return out;
}

inline std::vector<vector_type> zero_weights(size_t n_dof) {
    std::vector<vector_type> w(4);
    for (int j = 0; j < 4; ++j) w[j] = vector_type::Zero(n_dof);
    return w;
}

struct no_hook {
    template<typename... Args> void operator()(Args&&...) const {}
};

// Diaz-Grote-style multilevel recursion over L levels. Level i < L-1 runs
// p_level[i] sub-steps; at each one it recomputes its OWN contribution
// coarse(i, xn, own_w, t0, dtau, shifted) fresh from the CURRENT state
// (shifted = inherited Taylor weights moved to the sub-step start) and
// passes to level i+1:
//   accumulate == true : shifted + own_w (sums, exactly like Diaz-Grote's
//                        v := w + Sum w_k);
//   accumulate == false: own_w alone (the _v2/_v3 kernels fold `shifted`
//                        in through their ancestor_w argument; adding it
//                        again would double-count its order-0 term).
// before/after(i, xn, own_w, dtau) run around the recursion into level
// i+1. Only the terminal level gets genuine RK4 dynamics: fine(xn, w, t,
// len).
template<typename Coarse, typename Fine, typename Before = no_hook, typename After = no_hook>
void multilevel_step(int L, const std::vector<int>& p_level, size_t n_dof, bool accumulate,
                     Coarse&& coarse, Fine&& fine, vector_type& x, RealType tn, RealType dt_macro,
                     Before&& before = Before(), After&& after = After()) {
    std::function<void(int, RealType, RealType, vector_type&, const std::vector<vector_type>&)> recurse;
    recurse = [&](int i, RealType t_start, RealType len, vector_type& xn, const std::vector<vector_type>& w_accum) {
        if (i == L - 1) {
            fine(xn, w_accum, t_start, len);
            return;
        }
        int p_i = p_level[i];
        RealType dtau_i = len / p_i;
        for (int m = 0; m < p_i; ++m) {
            RealType tm = m * dtau_i;
            auto shifted = taylor_shift(w_accum, tm);
            std::vector<vector_type> own_w_i = zero_weights(n_dof);
            coarse(i, xn, own_w_i, t_start + tm, dtau_i, shifted);
            before(i, xn, own_w_i, dtau_i);
            if (accumulate) {
                std::vector<vector_type> combined(4);
                for (int j = 0; j < 4; ++j) combined[j] = shifted[j] + own_w_i[j];
                recurse(i + 1, t_start + tm, dtau_i, xn, combined);
            } else {
                recurse(i + 1, t_start + tm, dtau_i, xn, own_w_i);
            }
            after(i, xn, own_w_i, dtau_i);
        }
    };
    recurse(0, tn, dt_macro, x, zero_weights(n_dof));
}

// Closed-form advance of a level's own dofs over one sub-step by its own
// Taylor weights (skipping dofs in skip_c/skip_f).
inline void advance_own(vector_type& xn, const std::vector<vector_type>& own_w, RealType dtau,
                        const std::vector<int>& own_c, const std::vector<int>& own_f, size_t n_c_dof,
                        const std::set<int>& skip_c, const std::set<int>& skip_f) {
    RealType dt2 = dtau*dtau, dt3 = dt2*dtau, dt4 = dt3*dtau;
    for (int g : own_c) {
        if (skip_c.count(g)) continue;
        xn(g) += own_w[0](g)*dtau + own_w[1](g)*dt2/2.0 + own_w[2](g)*dt3/6.0 + own_w[3](g)*dt4/24.0;
    }
    for (int g : own_f) {
        if (skip_f.count(g)) continue;
        int gg = (int)n_c_dof + g;
        xn(gg) += own_w[0](gg)*dtau + own_w[1](gg)*dt2/2.0 + own_w[2](gg)*dt3/6.0 + own_w[3](gg)*dt4/24.0;
    }
}

// One macro step of the ORIGINAL, unrestricted two-level LTS-RK4 scheme
// (erk_weight_LTS_coarse on Pcoarse, then p sub-steps erk_weight_LTS_fine on
// Pfine) -- the trusted reference of the multi-level validation tests. When
// p == 1 nothing is done, unless fine_when_single_rate (one plain fine step).
template<typename Problem>
void two_level_reference_step(Problem& pb, const SparseMatrix<RealType>& Pcoarse, const SparseMatrix<RealType>& Pfine,
                              int p, vector_type& x, RealType tn, RealType dt_macro, bool fine_when_single_rate = false) {
    erk_type& erk_an = pb.erk();
    const size_t n_dof = x.rows();
    auto x_dof_n = x;
    std::vector<vector_type> wk(4);
    for (int j = 0; j < 4; ++j) { wk[j].resize(n_dof); wk[j].setZero(); }
    erk_an.ZeroFc();
    if (p != 1) {
        vector_type Fn = pb.rhs(tn), Fn12 = pb.rhs(tn + 0.5*dt_macro), Fn1 = pb.rhs(tn + dt_macro);
        erk_an.ZeroFc();
        erk_an.erk_weight_LTS_coarse(x_dof_n, Pcoarse, wk, Fn, Fn12, Fn1, dt_macro);
    }
    if (p != 1 || fine_when_single_rate) {
        RealType dtau = dt_macro / p;
        for (int m = 0; m < p; m++) {
            RealType tm = m*dtau, tmh = (m+0.5)*dtau, tm1 = (m+1.0)*dtau;
            vector_type Fm = pb.rhs(tn+tm), Fmh = pb.rhs(tn+tmh), Fm1 = pb.rhs(tn+tm1);
            erk_an.erk_weight_LTS_fine(x_dof_n, Pfine, wk, Fm, Fmh, Fm1, tm, dtau);
        }
    }
    x = x_dof_n;
}

// Projector on the own dofs of bands [first, last).
inline SparseMatrix<RealType> bands_projector(const bands& bd, int first, int last, size_t n_c_dof, size_t n_dof) {
    std::vector<int> c, f;
    for (int b = first; b < last; ++b) {
        c.insert(c.end(), bd.own_c[b].begin(), bd.own_c[b].end());
        f.insert(f.end(), bd.own_f[b].begin(), bd.own_f[b].end());
    }
    return band_projector(c, f, n_c_dof, n_dof);
}

// max|diff| over each band's own dofs (entries below 1e-12 ignored).
inline void print_band_max_diff(const vector_type& diff, const bands& bd, size_t n_c_dof) {
    std::vector<RealType> max_by_band(bd.L, 0.0);
    std::vector<std::set<int>> band_own_c(bd.L), band_own_f(bd.L);
    for (int b = 0; b < bd.L; ++b) {
        band_own_c[b] = std::set<int>(bd.own_c[b].begin(), bd.own_c[b].end());
        band_own_f[b] = std::set<int>(bd.own_f[b].begin(), bd.own_f[b].end());
    }
    for (int i = 0; i < diff.rows(); ++i) {
        if (diff(i) < 1e-12) continue;
        bool is_face = (size_t)i >= n_c_dof;
        int local = is_face ? (i - (int)n_c_dof) : i;
        for (int b = 0; b < bd.L; ++b) {
            bool in = is_face ? band_own_f[b].count(local) : band_own_c[b].count(local);
            if (in) { if (diff(i) > max_by_band[b]) max_by_band[b] = diff(i); break; }
        }
    }
    for (int b = 0; b < bd.L; ++b)
        std::cout << bold << cyan << "      max|diff| in band" << b << " own dofs: " << std::setprecision(6) << max_by_band[b] << reset << std::endl;
}

// ---------------------------------------------------------------------------
// Multi-level convergence sweep
// ---------------------------------------------------------------------------

// Level split and time steps of one mesh.
struct mlts_levels {
    int L;
    std::vector<RealType> thresholds, dt_level;
    std::vector<int> p_level;
    RealType dt_macro;
    size_t nt;
};

struct mlts_sweep_config {
    std::string title;           // sweep banner
    std::string title_note;      // optional yellow line right after it
    std::string mesh_family;     // directory under mesh_generation/
    std::string mesh_prefix;     // mesh file = <prefix>_k<k>_N<label>.txt
    std::vector<std::string> levels;
    std::string level_column = "N";
    std::string level_tag;       // appended to the "MESH LEVEL N=" banner
    std::string conv_stem;       // lshape/results/<stem>_k_<k><suffix>.txt
    std::string conv_suffix;
    // Level split from the mesh size range (threshold list, L = size+1).
    std::function<std::vector<RealType>(RealType h_max, RealType h_min, int p_global)> thresholds;
    bool strict_bands   = false; // '>' instead of '>=' in the band test
    RealType cfl_factor = 0.05;  // default, MLTS_CFL_FACTOR overrides
    RealType tf         = 0.1;
    bool velocity_error = false; // L2_error_velocity column
    bool cpu_time       = false; // cpu_time_s column (getrusage)
};

inline std::vector<std::string> level_labels(const std::vector<int>& level_N) {
    std::vector<std::string> labels;
    for (int N : level_N) labels.push_back(std::to_string(N));
    return labels;
}

// For every mesh level: problem setup, bands, then
// step = make_stepper(pb, bands) (the integrator's own preparation and
// band report) advances x over one macro step: step(x, tn, levels).
template<typename Solution, typename MakeStepper>
void mlts_sweep(int argc, char **argv, const mlts_sweep_config& cfg, const Solution& sol0, MakeStepper&& make_stepper)
{
    std::cout << std::endl << bold << red << cfg.title << std::endl << std::endl;
    if (!cfg.title_note.empty())
        std::cout << bold << yellow << cfg.title_note << reset << std::endl;

    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    disk::hho_degree_info hho_di = degree_info(sim_data);
    const size_t mesh_k = sim_data.m_k_degree;

    const std::string out_dir = "lshape/results";
    std::filesystem::create_directories(out_dir);
    std::ostringstream conv_fname;
    conv_fname << out_dir << "/" << cfg.conv_stem << "_k_" << mesh_k << cfg.conv_suffix << ".txt";
    std::ofstream conv_log(conv_fname.str());
    conv_log << "# " << cfg.level_column << "  h_max  h_min  p_global  L_levels  L2_error_pressure  wall_time_s"
             << (cfg.velocity_error ? "  L2_error_velocity" : "") << (cfg.cpu_time ? "  cpu_time_s" : "") << "\n";

    for (const std::string& N : cfg.levels) {

        std::cout << bold << red << "\n   ================ MESH LEVEL N=" << N << cfg.level_tag << " ================" << reset << std::endl;

        std::ostringstream mesh_fname;
        mesh_fname << mesh_root << cfg.mesh_family << "/meshes/" << cfg.mesh_prefix << "_k" << mesh_k << "_N" << N << ".txt";
        acoustic_problem<Solution> pb(mesh_fname.str(), sol0);
        const RealType h_max = pb.h_max, h_min = pb.h_min;
        const int p_global = pb.p_global;

        mlts_levels lv;
        lv.thresholds = cfg.thresholds(h_max, h_min, p_global);
        lv.L = (int)lv.thresholds.size() + 1;
        std::cout << bold << cyan << "      n_cells=" << pb.msh.cells_size()
                  << "  h_max=" << h_max << "  h_min=" << h_min
                  << "  p_global=" << p_global << "  L=" << lv.L << " levels" << reset << std::endl;

        pb.setup(hho_di, sim_data.m_scaled_stabilization_Q);
        bands bd = make_bands(pb, lv.thresholds, cfg.strict_bands);
        auto step = make_stepper(pb, bd);

        level_steps(env_real("MLTS_CFL_FACTOR", cfg.cfl_factor), h_max, h_min, lv.thresholds, lv.dt_level, lv.p_level);
        RealType ti = 0.0, tf = cfg.tf;
        lv.dt_macro = lv.dt_level[0];
        lv.nt = macro_steps(tf, lv.dt_macro);
        lv.dt_level[0] = lv.dt_macro;
        std::cout << bold << cyan << "      nt=" << lv.nt << reset << std::endl;

        timecounter tc_run;
        tc_run.tic();
        RealType cpu_start = cpu_seconds();
        for (size_t it = 1; it <= lv.nt; it++) {
            RealType tn = ti + lv.dt_macro * (it - 1);
            step(pb.x_dof, tn, lv);
            if (it % std::max<size_t>(1, lv.nt/10) == 0)
                std::cout << bold << yellow << "         step " << it << "/" << lv.nt << reset << std::endl;
        }
        tc_run.toc();
        RealType cpu_time = cpu_seconds() - cpu_start;

        RealType l2_error   = pb.pressure_error(pb.x_dof, tf);
        RealType v_l2_error = cfg.velocity_error ? pb.velocity_error(pb.x_dof, tf) : 0.0;

        std::cout << bold << red << "      wall time: " << tc_run << " s   ";
        if (cfg.cpu_time) std::cout << "cpu time: " << cpu_time << " s   ";
        std::cout << "L2_error(pressure) = " << std::setprecision(10) << l2_error;
        if (cfg.velocity_error) std::cout << "   L2_error(velocity) = " << v_l2_error;
        std::cout << reset << std::endl;

        // Optional columns are APPENDED so older parsers keep working.
        conv_log << N << " " << std::setprecision(15) << h_max << " " << h_min << " " << p_global << " " << lv.L
                  << " " << l2_error << " " << tc_run.elapsed();
        if (cfg.velocity_error) conv_log << " " << v_l2_error;
        if (cfg.cpu_time) conv_log << " " << cpu_time;
        conv_log << "\n";
        conv_log.flush();
    }

    conv_log.close();
    std::cout << bold << red << "\n   --> " << conv_fname.str() << reset << std::endl << std::endl;
}

} // namespace lts
