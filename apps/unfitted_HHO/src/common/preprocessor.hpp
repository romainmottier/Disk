#pragma once
#ifndef preprocessor_hpp
#define preprocessor_hpp

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>  // for getopt

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/output/cuthho_output.hpp"
#include "meshes.hpp"
#include "refinement.hpp"

#define reset   "\033[0m"
#define red     "\033[31m"
#define bold    "\033[1m"
#define yellow  "\033[33m"
#define blue    "\033[34m"
#define cyan    "\033[36m"

/* Meshes (-M) */
enum class mesh_kind {
    quadrangles = 0,  // Cartesian grid of the square (paper)
    triangles = 1,    // right triangles
    polygons = 2,     // hexagons (FVCA5)
    hexahedra = 3,    // Cartesian grid of the cube
    tetrahedra = 4,   // Kuhn subdivision of the cubes
    polyhedra = 5     // hexagonal prisms
};

inline const char*
mesh_kind_name(size_t m) {
    static const char* names[] = {"quadrangles", "triangles", "polygons (hexagons)", "hexahedra",
                                  "tetrahedra", "polyhedra (hexagonal prisms)"};
    return (m < 6) ? names[m] : "unknown";
}

struct simulation_options {
    size_t prototype = 0;
    size_t mesh = 0;           // -M: mesh_kind
    size_t l_divs = 0;         // -l: space refinement levels
    size_t lv_set = 0;         // -i: 0 circle/sphere, 1 flower (55b), 2 square (2D), 3 flower of ProtoN (2D),
                               //     4 ellipse/ellipsoid, 5 rounded square/cube, 6 annulus/torus
    size_t degree = 0;         // -k: face degree (cell degree k+1)
    size_t int_refsteps = 8;   // -r: the polyline interface has 2^r segments in each cut cell
#ifdef HAVE_INTEL_MKL
    size_t iterative_solver = 2;  // -s: 0 SimplicialLDLT, 1 CG, 2 Pardiso LLT
#else
    size_t iterative_solver = 0;  // -s: 0 SimplicialLDLT, 1 CG, 2 Pardiso LLT (needs MKL)
#endif
    size_t test_case = 0;      // -t
    size_t debug = 0;          // -f: silo files and pairing arrows
    size_t contrast = 0;       // -m: kappa_2 = 10^m kappa_1
    double threshold = 0.3;    // -a: flagging parameter of the ill-cut cells (vartheta)
    double eta_N = -20.0;      // -e: s^N parameter: x >= 0 fixed, x < 0 |x| (k+1)^2 (-1: ProtoN)
    size_t basis = 2;          // -b: cell bases: 0 monomials of the cells, 1 of the subcells (paper), 2 orthonormal
    size_t pairing = 1;        // -P: neighbor of an ill-cut cell: 0 ProtoN (paper), 1 largest subcell
    size_t condensation = 1;   // -c: static condensation: 0 none, 1 local (patches), 2 algebraic
    int    curve = -1;         // -g: interface: -1 exact (implicit), 0 polyline of -r, q >= 2 polynomial of degree q, -2 degree k+2
    size_t cache_uncut = 1;    // -u: compute the local operators of the uncut cells once per shape
    double radius = 1.0/3.0;   // -R: radius of the circle/flower/sphere
    size_t square_p = 2;       // -d: square interface shifted by 0.5*10^-p
    size_t conditioning = 0;   // -q: compute the condition number of the stiffness matrix
    size_t sparsity = 0;       // -z: dump the sparsity pattern of the stiffness matrix
    size_t k_min = 0;          // -K: loop over the degrees k_min..k
    size_t l_min = 0;          // -L: loop over the levels l_min..l
    size_t refinement = 3;     // -x: the cells that the interface does not cut once are refined, at most x times
                               //     (polytopal meshes; 0: never)
    bool   refined_mesh = false; // set by run_mesh: some cells of the meshes are not cut once, -x > 0
};

inline void
simulation_data(int argc, char** argv, simulation_options& opts) {
    int ch;
    optind = 1;
    while ((ch = getopt(argc, argv, "p:M:l:i:k:r:s:t:f:m:a:e:b:P:c:g:u:R:d:q:z:K:L:x:")) != -1) {
        switch (ch) {
            case 'p': opts.prototype = std::atoi(optarg); break;
            case 'M': opts.mesh = std::atoi(optarg); break;
            case 'l': opts.l_divs = std::atoi(optarg); break;
            case 'i': opts.lv_set = std::atoi(optarg); break;
            case 'k': opts.degree = std::atoi(optarg); break;
            case 'r': opts.int_refsteps = std::atoi(optarg); break;
            case 's': opts.iterative_solver = std::atoi(optarg); break;
            case 't': opts.test_case = std::atoi(optarg); break;
            case 'f': opts.debug = std::atoi(optarg); break;
            case 'm': opts.contrast = std::atoi(optarg); break;
            case 'a': opts.threshold = std::atof(optarg); break;
            case 'e': opts.eta_N = std::atof(optarg); break;
            case 'b': opts.basis = std::atoi(optarg); break;
            case 'P': opts.pairing = std::atoi(optarg); break;
            case 'c': opts.condensation = std::atoi(optarg); break;
            case 'g': opts.curve = std::atoi(optarg); break;
            case 'u': opts.cache_uncut = std::atoi(optarg); break;
            case 'R': opts.radius = std::atof(optarg); break;
            case 'd': opts.square_p = std::atoi(optarg); break;
            case 'q': opts.conditioning = std::atoi(optarg); break;
            case 'z': opts.sparsity = std::atoi(optarg); break;
            case 'K': opts.k_min = std::atoi(optarg); break;
            case 'L': opts.l_min = std::atoi(optarg); break;
            case 'x': opts.refinement = std::atoi(optarg); break;
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    std::cout << "    " << bold << red << "SIMULATION DATA:" << reset << blue << bold << std::endl;
    std::cout << "        Prototype                  -p : " << opts.prototype << std::endl;
    std::cout << "        Mesh                       -M : " << opts.mesh << " (" << mesh_kind_name(opts.mesh) << ")" << std::endl;
    std::cout << "        Space refinement level     -l : " << opts.l_divs << std::endl;
    std::cout << "        Polynomial degree          -k : " << opts.degree << "     (Face unknowns)" << std::endl;
    std::cout << "        Interface refinement level -r : " << opts.int_refsteps << std::endl;
    std::cout << "        Test case                  -t : " << opts.test_case << std::endl;
    std::cout << "        Contrast kappa_2/kappa_1   -m : 10^" << opts.contrast << std::endl;
    std::cout << "        Ill-cut threshold          -a : " << opts.threshold << std::endl;
    std::cout << "        eta_N                      -e : ";
    if (opts.eta_N < 0) std::cout << -opts.eta_N << " (k+1)^2" << std::endl;
    else                std::cout << opts.eta_N << std::endl;
    std::cout << "        Cell bases                 -b : " << opts.basis
              << (opts.basis == 0 ? " (monomials of the cells)" : opts.basis == 1 ? " (monomials of the subcells)" : " (orthonormal)")
              << std::endl;
    std::cout << "        Pairing of ill-cut cells   -P : " << opts.pairing
              << (opts.pairing == 0 ? " (ProtoN)" : " (largest subcell)") << std::endl;
    std::cout << "        Static condensation        -c : " << opts.condensation << std::endl;
    std::cout << "        Interface                  -g : ";
    if (opts.curve == -1)      std::cout << "exact (root of the level set)" << std::endl;
    else if (opts.curve == -2) std::cout << "polynomial of degree k+2" << std::endl;
    else if (opts.curve == 0)  std::cout << "0 (polyline of 2^r segments)" << std::endl;
    else                       std::cout << "polynomial of degree " << opts.curve << std::endl;
    std::cout << "        Local refinement           -x : " << opts.refinement
              << (opts.refinement > 0 ? " (cells not cut once, if any)" : " (none)") << std::endl;
    std::cout << "        Debug & Silo files         -f : " << opts.debug << reset << std::endl;
}

/* Representation of the interface for the degree k (option -g) */
inline disk::cut_geometry_options
geometry_options(const simulation_options& opts, size_t k)
{
    disk::cut_geometry_options g;
    g.refine_levels = opts.int_refsteps;
    if (opts.curve == -1)
        g.kind = disk::interface_kind::implicit;
    else if (opts.curve == 0)
        g.kind = disk::interface_kind::polyline;
    else {
        g.kind = disk::interface_kind::polynomial;
        g.degree = (opts.curve < 0) ? k + 2 : size_t(opts.curve);
    }
    return g;
}

template<typename Mesh>
class preprocessor {

    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;

public:
    using cut_mesh_type = disk::cut_mesh<Mesh>;

    static std::shared_ptr<level_set<T, DIM>>
    make_level_set_function(const simulation_options& opts, bool verbose = true) {

        std::ostream null_stream(nullptr);
        std::ostream& out = verbose ? std::cout : null_stream;

        /* interfaces of any dimension, rotated so as not to be aligned with
         * the meshes */
        disk::point<T,DIM> c;
        for (size_t d = 0; d < DIM; d++)
            c[d] = 0.5;
        auto Q = (DIM == 2) ? rotation_matrix<T,DIM>(0.4) : rotation_matrix<T,DIM>(0.5, 0.3, 0.7);
        switch (opts.lv_set) {
            case 4: {
                Eigen::Matrix<T,DIM,1> axes;
                if constexpr (DIM == 2) axes << 0.36, 0.2;
                else                    axes << 0.36, 0.26, 0.18;
                out << "        Level set function         -i : " << (DIM == 2 ? "Ellipse" : "Ellipsoid")
                    << " with semi-axes " << axes.transpose() << std::endl;
                return std::make_shared<superellipse_level_set<T,DIM>>(c, axes, Q, 2);
            }
            case 5: {
                Eigen::Matrix<T,DIM,1> axes = Eigen::Matrix<T,DIM,1>::Constant(0.3);
                out << "        Level set function         -i : Rounded " << (DIM == 2 ? "square" : "cube")
                    << " (sum y_d^4 = 0.3^4)" << std::endl;
                return std::make_shared<superellipse_level_set<T,DIM>>(c, axes, Q, 4);
            }
            case 6: {
                T R = (DIM == 2) ? 0.3 : 0.28, r = (DIM == 2) ? 0.12 : 0.14;
                out << "        Level set function         -i : " << (DIM == 2 ? "Annulus" : "Torus")
                    << " with R = " << R << ", r = " << r << std::endl;
                return std::make_shared<torus_level_set<T,DIM>>(c, R, r, Q);
            }
            default:
                break;
        }

        if constexpr (DIM == 3) {
            switch (opts.lv_set) {
                case 0:
                    out << "        Level set function         -i : Sphere with R = " << opts.radius << std::endl;
                    return std::make_shared<circle_level_set<T, 3>>(opts.radius, c);
                case 1:
                    /* 3D analog of (55b): 4 petals around the axis z */
                    out << "        Level set function         -i : Flower (sphere with 4 petals) with R = "
                        << opts.radius << std::endl;
                    return std::make_shared<star_level_set<T, 3>>(c, opts.radius, 0.03, 4);
                default:
                    throw std::runtime_error("Level set " + std::to_string(opts.lv_set) + " not available in 3D");
            }
        }
        else {
            switch (opts.lv_set) {
                case 0:
                    out << "        Level set function         -i : Circle with R = " << opts.radius << std::endl;
                    return std::make_shared<circle_level_set<T>>(opts.radius, 0.5, 0.5);
                case 1:
                    /* eq. (55b) of the paper: (x-a)^2 + (y-b)^2 - R^2 + c cos(n theta) */
                    out << "        Level set function         -i : Flower with R = " << opts.radius << std::endl;
                    return std::make_shared<flower_level_set<T>>(opts.radius, 0.5, 0.5, 8, -0.03);
                case 3:
                    /* flower of the ProtoN code (- c cos(n theta)), rotated by pi/8 */
                    out << "        Level set function         -i : Flower (ProtoN) with R = " << opts.radius << std::endl;
                    return std::make_shared<flower_level_set<T>>(opts.radius, 0.5, 0.5, 8, 0.03);
                case 2: {
                    T delta = 0.5*std::pow(10.0, -T(opts.square_p));
                    out << "        Level set function         -i : Square with delta = " << delta << std::endl;
                    return std::make_shared<square_level_set<T>>(0.70 + delta, 0.30 - delta, 0.25, 0.75);
                }
                default:
                    throw std::runtime_error("Unsupported level set type: " + std::to_string(opts.lv_set));
            }
        }
    }

    /* Background mesh of the level l: (10*2^l)^2 squares or right triangles,
     * hexagons of the level l+2 of FVCA5, (5*2^l)^3 cubes or Kuhn
     * tetrahedra, prisms on the hexagons of the level l+1 in 5*2^l layers */
    static void
    background_mesh(Mesh& msh, size_t l) {
        size_t N2 = 10 * (size_t(1) << l);
        size_t N3 = 5 * (size_t(1) << l);
        if constexpr (std::is_same_v<Mesh, disk::cartesian_mesh<T,2>>)
            unfitted_meshes::make_cartesian_grid(msh, N2, N2);
        else if constexpr (std::is_same_v<Mesh, disk::simplicial_mesh<T,2>>)
            unfitted_meshes::make_triangle_grid(msh, N2);
        else if constexpr (std::is_same_v<Mesh, disk::generic_mesh<T,2>>)
            unfitted_meshes::make_hexagon_mesh(msh, l + 2);
        else if constexpr (std::is_same_v<Mesh, disk::cartesian_mesh<T,3>>)
            unfitted_meshes::make_hexahedral_grid(msh, N3);
        else if constexpr (std::is_same_v<Mesh, disk::simplicial_mesh<T,3>>)
            unfitted_meshes::make_tetrahedral_grid(msh, N3);
        else {
            disk::generic_mesh<T,2> m2;
            unfitted_meshes::make_hexagon_mesh(m2, l + 1);
            unfitted_meshes::make_prism_mesh(msh, m2, N3);
        }
    }

    /* Background mesh of the level l of the kind -M, refined where the
     * interface does not cut the cells once (-x): a polytopal mesh */
    static void
    refined_background_mesh(Mesh& msh, size_t l, const simulation_options& opts, const level_set<T, DIM>& ls,
                            const disk::cut_geometry_options& geometry) {
        if constexpr (!std::is_same_v<Mesh, disk::generic_mesh<T,DIM>>)
            throw std::logic_error("The locally refined meshes are polytopal meshes");
        else {
            size_t N2 = 10 * (size_t(1) << l);
            size_t N3 = 5 * (size_t(1) << l);
            unfitted_meshes::refinable_mesh<T, DIM> rm;
            auto kind = mesh_kind(opts.mesh);
            if constexpr (DIM == 2) {
                if (kind == mesh_kind::quadrangles || kind == mesh_kind::triangles)
                    unfitted_meshes::refinable_square_grid(rm, N2, kind == mesh_kind::triangles);
                else {
                    disk::generic_mesh<T,2> m2;
                    unfitted_meshes::make_hexagon_mesh(m2, l + 2);
                    unfitted_meshes::refinable_polygonal_mesh(rm, m2);
                }
            }
            else {
                if (kind == mesh_kind::hexahedra || kind == mesh_kind::tetrahedra)
                    unfitted_meshes::refinable_cube_grid(rm, N3, kind == mesh_kind::tetrahedra);
                else {
                    disk::generic_mesh<T,2> m2;
                    unfitted_meshes::make_hexagon_mesh(m2, l + 1);
                    unfitted_meshes::refinable_prism_mesh(rm, m2, N3);
                }
            }
            size_t n0 = rm.cells.size();
            auto st = unfitted_meshes::refine_until_cut_once(rm, msh, ls, geometry, opts.refinement);
            std::cout << "         Local refinement: " << st.refined << " cells refined (up to " << st.levels
                      << " times), " << n0 << " -> " << rm.cells.size() << " cells" << std::endl;
            if (st.unresolved > 0)
                std::cout << "         Warning: level " << l << ": " << st.unresolved
                          << " cells are still not cut once after " << opts.refinement << " refinements" << std::endl;
        }
    }

    /* Do the meshes of the levels of the run have cells that the interface
     * does not cut once? */
    static bool
    needs_refinement(const simulation_options& opts) {
        auto ls = make_level_set_function(opts, false);
        for (size_t l = std::min(opts.l_min, opts.l_divs); l <= opts.l_divs; l++) {
            cut_mesh_type msh;
            background_mesh(msh, l);
            auto bad = disk::cut_once_violations(msh, *ls, geometry_options(opts, opts.degree));
            if (std::find(bad.begin(), bad.end(), true) != bad.end())
                return true;
        }
        return false;
    }

    /* Mesh of the level l cut by the interface, with the pairing of the
     * polynomial extension. */
    static cut_mesh_type
    msh_generation(std::shared_ptr<level_set<T, DIM>> level_set_function, size_t l, const simulation_options& opts,
                   const disk::cut_geometry_options& geometry, bool pairing = true) {

        cut_mesh_type msh;
        if (opts.refined_mesh)
            refined_background_mesh(msh, l, opts, *level_set_function, geometry);
        else
            background_mesh(msh, l);

        auto stats = disk::make_cut_mesh(msh, level_set_function, geometry, T(opts.threshold));
        if (stats.clamped > 0)
            std::cout << "         Warning: level " << l << ": the interface crosses a face twice in " << stats.clamped
                      << " cells (the mesh does not resolve the interface: a cell must be cut once);"
                      << " the interface is projected onto their faces (-x refines these cells)" << std::endl;
        if (stats.fallbacks > 0)
            std::cout << "         Warning: level " << l << ": the curved interface could not be built in "
                      << stats.fallbacks << (DIM == 2 ? " cells" : " faces") << ", which keep the polyline" << std::endl;
        if (stats.split > 0 || stats.subdivided > 0)
            std::cout << "         Level " << l << ": " << stats.cut_cells << " cut cells, " << stats.split
                      << " split into parts, " << stats.subdivided << " with parts integrated by subdivision" << std::endl;
        if (pairing)
            disk::make_polynomial_extension(msh, opts.debug != 0, opts.pairing != 0);

        return msh;
    }

    static void
    msh_debug(cut_mesh_type& msh, const level_set<T, DIM>& level_set_function, size_t l,
              const disk::cut_agglomeration& agglo = disk::cut_agglomeration()) {

        output_mesh_info(msh, level_set_function);

        auto mv = [](const std::string& from, const std::string& to) {
            std::rename(from.c_str(), to.c_str());
        };
        auto suffix = "_l" + std::to_string(l);
        mv("cuthho_meshinfo.silo", "cuthho_meshinfo" + suffix + ".silo");
        mv("interface.3D", "interface" + suffix + ".3D");
        mv("agglo_four.okc", "agglo_four" + suffix + ".okc");
        mv("agglo_five.okc", "agglo_five" + suffix + ".okc");

        write_cells_info(msh, "cells" + suffix + ".txt", agglo);
    }

    /* Plain text description of the cut mesh (for plotting, Fig. 5):
     * "cell id min(x) min(y) [min(z)] max(x) max(y) [max(z)] location cut_type stabilizing_cell agglomeration_group" and
     * "interface cell x0 y0 [z0] x1 y1 [z1] ..." */
    static void
    write_cells_info(const cut_mesh_type& msh, const std::string& filename,
                     const disk::cut_agglomeration& agglo = disk::cut_agglomeration()) {

        std::ofstream ofs(filename);
        for (auto& cl : msh) {
            auto cl_id = offset(msh, cl);
            const auto& ccd = msh.cut_cell_data[cl_id];
            auto pts = points(msh, cl);
            auto lo = pts[0], hi = pts[0];
            for (auto& p : pts)
                for (size_t d = 0; d < DIM; d++) {
                    lo[d] = std::min(lo[d], p[d]);
                    hi[d] = std::max(hi[d], p[d]);
                }
            int stab = ccd.stabilizing_cell.has_value() ? int(ccd.stabilizing_cell.value()) : -1;
            int group = agglo.group_of_cell.empty() ? -1 : agglo.group_of_cell[cl_id];
            ofs << "cell " << cl_id;
            for (size_t d = 0; d < DIM; d++)
                ofs << " " << lo[d];
            for (size_t d = 0; d < DIM; d++)
                ofs << " " << hi[d];
            ofs << " " << int(ccd.loc) << " " << int(ccd.cut) << " " << stab << " " << group << std::endl;
        }
        for (auto& cl : msh) {
            if (!disk::is_cut(msh, cl))
                continue;
            std::vector<typename Mesh::point_type> itf;
            if constexpr (DIM == 2)
                itf = disk::cut_info(msh, cl).poly.polyline;
            else
                disk::cutgeo::polyhedron_interface_samples(disk::cut_info(msh, cl).poly, level_set_function_of(msh), itf);
            ofs << "interface " << offset(msh, cl);
            for (auto& p : itf)
                for (size_t d = 0; d < DIM; d++)
                    ofs << " " << p[d];
            ofs << std::endl;
        }
    }

private:
    static const level_set<T, DIM>&
    level_set_function_of(const cut_mesh_type& msh) { return *msh.level_set_function; }
};

#endif // preprocessor_hpp
