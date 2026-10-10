/*
 * Global assembler of the unfitted HHO method stabilized by polynomial
 * extension (ported from ProtoN's interface_assembler).
 *
 * Global unknowns: cell unknowns first (two copies for the cut cells,
 * negative side first), then face unknowns (two copies for the cut faces),
 * the faces lying on the (Dirichlet) boundary being eliminated.
 */

#pragma once

#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Sparse>

#include "diskpp/mesh/cut_mesh.hpp"
#include "diskpp/bases/cut_bases.hpp"
#include "assembly_index.hpp"
#include "operators.hpp"
#include "agglo_operators.hpp"

template<disk::CutMesh Mesh>
class unfitted_elliptic_interface_assembler {

    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    using face_type = typename Mesh::face_type;
    using cell_type = typename Mesh::cell_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    static constexpr size_t not_assembled = std::numeric_limits<size_t>::max();

    disk::hho_degree_info               m_hdi;
    std::function<T(const point_type&)> m_dirichlet_fun;
    disk::cut_basis_options<T, Mesh::dimension> m_basis;
    std::vector<size_t>                 m_cell_table;
    std::vector<size_t>                 m_face_table;
    size_t                              m_num_cells;
    size_t                              m_num_other_faces;
    size_t                              m_cbs, m_fbs;
    std::vector<Eigen::Triplet<T>>      m_triplets;

    /* Block sparsity pattern (Fig. 8): one block per cell or face unknown,
     * the couplings involving the unknowns of a dependent ill-cut cell being
     * flagged as resulting from the polynomial extension */
    bool                                m_record_pattern = false;
    std::set<std::pair<size_t, size_t>> m_local_blocks, m_extension_blocks;

public:
    Eigen::SparseMatrix<T>  LHS;
    vector_type             RHS;

    unfitted_elliptic_interface_assembler(const Mesh& msh, const disk::hho_degree_info& hdi,
                                          std::function<T(const point_type&)> dirichlet_fun,
                                          const disk::cut_basis_options<T, Mesh::dimension>& basis = {})
        : m_hdi(hdi), m_dirichlet_fun(dirichlet_fun), m_basis(basis)
    {
        m_cbs = disk::cut_scalar_basis<Mesh>::size(m_hdi.cell_degree());
        m_fbs = disk::cut_face_basis<Mesh>::size(m_hdi.face_degree());

        m_num_cells = 0;
        for (auto& cl : msh) {
            m_cell_table.push_back(m_num_cells);
            m_num_cells += disk::is_cut(msh, cl) ? 2 : 1;
        }

        m_num_other_faces = 0;
        for (auto& fc : faces(msh)) {
            if (msh.is_boundary(fc)) {
                m_face_table.push_back(not_assembled);
                continue;
            }
            m_face_table.push_back(m_num_other_faces);
            m_num_other_faces += disk::is_cut(msh, fc) ? 2 : 1;
        }

        auto system_size = m_cbs*m_num_cells + m_fbs*m_num_other_faces;
        LHS = Eigen::SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    size_t num_cell_dofs() const { return m_cbs*m_num_cells; }

    size_t
    face_SOL_offset(const Mesh& msh, const face_type& fc) const {
        return m_num_cells*m_cbs + m_face_table.at(offset(msh, fc))*m_fbs;
    }

    size_t
    cell_SOL_offset(const Mesh& msh, const cell_type& cl, disk::location where) const {
        auto ofs = m_cell_table.at(offset(msh, cl))*m_cbs;
        if (disk::is_cut(msh, cl) && where == disk::location::IN_POSITIVE_SIDE)
            ofs += m_cbs;
        return ofs;
    }

private:
    /* Global indices of a block [cells | faces_1 (| faces_2)] of a cell; the
     * Dirichlet data of the non-assembled face unknowns is stored in dirichlet. */
    void
    add_cell_block(const Mesh& msh, const cell_type& cl, bool two_sides, disk::location single_side,
                   std::vector<assembly_index>& asm_map, std::vector<T>& dirichlet)
    {
        auto cell_LHS_offset = m_cell_table.at(offset(msh, cl))*m_cbs;
        auto ncbs = two_sides ? 2*m_cbs : m_cbs;
        for (size_t i = 0; i < ncbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
            dirichlet.push_back(0.0);
        }

        auto fcs = faces(msh, cl);
        std::vector<disk::location> sides;
        if (two_sides)
            sides = {disk::location::IN_NEGATIVE_SIDE, disk::location::IN_POSITIVE_SIDE};
        else
            sides = {single_side};

        for (size_t s = 0; s < sides.size(); s++) {
            for (auto& fc : fcs) {
                bool dirichlet_face = msh.is_boundary(fc);
                vector_type dir_data = vector_type::Zero(m_fbs);
                size_t face_LHS_offset = 0;
                if (dirichlet_face)
                    dir_data = dirichlet_projection(msh, fc, sides[s]);
                else {
                    /* the second copy of a cut face is the positive side */
                    auto d = (s == 1 && disk::is_cut(msh, fc)) ? m_fbs : 0;
                    face_LHS_offset = face_SOL_offset(msh, fc) + d;
                }
                for (size_t i = 0; i < m_fbs; i++) {
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet_face) );
                    dirichlet.push_back(dir_data(i));
                }
            }
        }
    }

    vector_type
    dirichlet_projection(const Mesh& msh, const face_type& fc, disk::location where) const
    {
        vector_type ret = vector_type::Zero(m_fbs);
        auto loc = disk::locate(msh, fc);
        if (loc != where && loc != disk::location::ON_INTERFACE)
            return ret;

        auto facdeg = m_hdi.face_degree();
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, where, m_basis);
        matrix_type mass = matrix_type::Zero(m_fbs, m_fbs);
        vector_type rhs  = vector_type::Zero(m_fbs);
        auto qps = integrate(msh, fc, 2*facdeg + 2, where);
        for (auto& qp : qps) {
            auto phi = fb.eval_basis(qp.point());
            mass += qp.weight() * phi * phi.transpose();
            rhs  += qp.weight() * phi * m_dirichlet_fun(qp.point());
        }
        return mass.ldlt().solve(rhs);
    }

public:
    std::vector<assembly_index>
    init_asm_map(const Mesh& msh, const disk::cut_subcell& P, std::vector<T>& dirichlet)
    {
        std::vector<assembly_index> asm_map;
        dirichlet.clear();

        auto cl = msh[P.cell];
        add_cell_block(msh, cl, disk::is_cut(msh, cl), P.side, asm_map, dirichlet);
        for (auto dp : P.dependent_cells)
            add_cell_block(msh, msh[dp], true, P.side, asm_map, dirichlet);

        return asm_map;
    }

    void
    assemble(const Mesh& msh, const disk::cut_subcell& P, const matrix_type& lhs, const vector_type& rhs)
    {
        std::vector<T> dirichlet;
        auto asm_map = init_asm_map(msh, P, dirichlet);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());
        assert(asm_map.size() == rhs.rows());

        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble())
                    m_triplets.push_back( Eigen::Triplet<T>(size_t(asm_map[i]), size_t(asm_map[j]), lhs(i,j)) );
                else
                    RHS(size_t(asm_map[i])) -= lhs(i,j)*dirichlet[j];
            }
            RHS(size_t(asm_map[i])) += rhs(i);
        }

        if (m_record_pattern) {
            subcell_layout<Mesh> lay(msh, P, m_hdi);
            for (size_t i = 0; i < lhs.rows(); i++) {
                if (!asm_map[i].assemble())
                    continue;
                for (size_t j = 0; j < lhs.cols(); j++) {
                    if (!asm_map[j].assemble())
                        continue;
                    auto bij = std::make_pair(block_of(asm_map[i]), block_of(asm_map[j]));
                    if (i < lay.current_dofs && j < lay.current_dofs)
                        m_local_blocks.insert(bij);
                    else
                        m_extension_blocks.insert(bij);
                }
            }
        }
    }

    void record_pattern(bool rec) { m_record_pattern = rec; }

    /* Block index of a global unknown: cell blocks first, then face blocks */
    size_t
    block_of(size_t dof) const {
        auto ncd = num_cell_dofs();
        if (dof < ncd)
            return dof / m_cbs;
        return m_num_cells + (dof - ncd) / m_fbs;
    }

    /* Write the block sparsity pattern: "i j flag", flag = 1 if the coupling
     * results from the polynomial extension */
    void
    write_block_pattern(const std::string& filename) const {
        std::ofstream ofs(filename);
        ofs << "# cell_blocks " << m_num_cells << " face_blocks " << m_num_other_faces << std::endl;
        for (auto& b : m_local_blocks)
            if (m_extension_blocks.find(b) == m_extension_blocks.end())
                ofs << b.first << " " << b.second << " 0" << std::endl;
        for (auto& b : m_extension_blocks)
            ofs << b.first << " " << b.second << " 1" << std::endl;
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }

    /* Cell unknowns of the side "where" of the cell cl */
    vector_type
    take_cell_dofs(const Mesh& msh, const cell_type& cl, const vector_type& solution, disk::location where) const
    {
        return solution.segment(cell_SOL_offset(msh, cl, where), m_cbs);
    }
};

/*
 * Global assembler of the method stabilized by cell agglomeration: one block
 * of cell unknowns per macro cell (two if it is cut), no unknown on the faces
 * interior to a macro cell.
 */
template<disk::CutMesh Mesh>
class unfitted_agglomeration_assembler {

    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;
    using face_type = typename Mesh::face_type;
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    static constexpr size_t not_assembled = std::numeric_limits<size_t>::max();

    disk::hho_degree_info               m_hdi;
    std::function<T(const point_type&)> m_dirichlet_fun;
    disk::cut_basis_options<T, Mesh::dimension> m_basis;
    std::vector<size_t>                 m_cell_table;   // per macro cell
    std::vector<size_t>                 m_face_table;   // per mesh face
    size_t                              m_num_cells;
    size_t                              m_num_other_faces;
    size_t                              m_cbs, m_fbs;
    std::vector<Eigen::Triplet<T>>      m_triplets;
    bool                                m_record_pattern = false;
    std::set<std::pair<size_t, size_t>> m_blocks;

public:
    Eigen::SparseMatrix<T>  LHS;
    vector_type             RHS;

    unfitted_agglomeration_assembler(const Mesh& msh, const std::vector<macro_cell<Mesh>>& macros,
                                     const disk::hho_degree_info& hdi,
                                     std::function<T(const point_type&)> dirichlet_fun,
                                     const disk::cut_basis_options<T, Mesh::dimension>& basis = {})
        : m_hdi(hdi), m_dirichlet_fun(dirichlet_fun), m_basis(basis)
    {
        m_cbs = disk::cut_scalar_basis<Mesh>::size(m_hdi.cell_degree());
        m_fbs = disk::cut_face_basis<Mesh>::size(m_hdi.face_degree());

        m_num_cells = 0;
        for (auto& Q : macros) {
            m_cell_table.push_back(m_num_cells);
            m_num_cells += Q.cut ? 2 : 1;
        }

        /* faces carrying unknowns: the boundary faces of the macro cells */
        std::vector<bool> has_unknowns(msh.faces_size(), false);
        for (auto& Q : macros)
            for (auto& fm : Q.faces)
                has_unknowns[offset(msh, fm.first)] = true;

        m_num_other_faces = 0;
        for (auto& fc : faces(msh)) {
            auto fid = offset(msh, fc);
            if (msh.is_boundary(fc) || !has_unknowns[fid]) {
                m_face_table.push_back(not_assembled);
                continue;
            }
            m_face_table.push_back(m_num_other_faces);
            m_num_other_faces += disk::is_cut(msh, fc) ? 2 : 1;
        }

        auto system_size = m_cbs*m_num_cells + m_fbs*m_num_other_faces;
        LHS = Eigen::SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);
    }

    size_t num_cell_dofs() const { return m_cbs*m_num_cells; }

    size_t
    face_SOL_offset(const Mesh& msh, const face_type& fc) const {
        return m_num_cells*m_cbs + m_face_table.at(offset(msh, fc))*m_fbs;
    }

    size_t
    cell_SOL_offset(size_t q, const macro_cell<Mesh>& Q, disk::location where) const {
        auto ofs = m_cell_table.at(q)*m_cbs;
        if (Q.cut && where == disk::location::IN_POSITIVE_SIDE)
            ofs += m_cbs;
        return ofs;
    }

private:
    vector_type
    dirichlet_projection(const Mesh& msh, const face_type& fc, disk::location where) const
    {
        vector_type ret = vector_type::Zero(m_fbs);
        auto loc = disk::locate(msh, fc);
        if (loc != where && loc != disk::location::ON_INTERFACE)
            return ret;

        auto facdeg = m_hdi.face_degree();
        disk::cut_face_basis<Mesh> fb(msh, fc, facdeg, where, m_basis);
        matrix_type mass = matrix_type::Zero(m_fbs, m_fbs);
        vector_type rhs  = vector_type::Zero(m_fbs);
        for (auto& qp : integrate(msh, fc, 2*facdeg + 2, where)) {
            auto phi = fb.eval_basis(qp.point());
            mass += qp.weight() * phi * phi.transpose();
            rhs  += qp.weight() * phi * m_dirichlet_fun(qp.point());
        }
        return mass.ldlt().solve(rhs);
    }

public:
    std::vector<assembly_index>
    init_asm_map(const Mesh& msh, size_t q, const macro_cell<Mesh>& Q, disk::location side, std::vector<T>& dirichlet)
    {
        std::vector<assembly_index> asm_map;
        dirichlet.clear();

        auto cell_LHS_offset = m_cell_table.at(q)*m_cbs;
        auto ncbs = Q.cut ? 2*m_cbs : m_cbs;
        for (size_t i = 0; i < ncbs; i++) {
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
            dirichlet.push_back(0.0);
        }

        std::vector<disk::location> sides;
        if (Q.cut)
            sides = {disk::location::IN_NEGATIVE_SIDE, disk::location::IN_POSITIVE_SIDE};
        else
            sides = {side};

        for (size_t s = 0; s < sides.size(); s++) {
            for (auto& fm : Q.faces) {
                const auto& fc = fm.first;
                bool dirichlet_face = msh.is_boundary(fc);
                vector_type dir_data = vector_type::Zero(m_fbs);
                size_t face_LHS_offset = 0;
                if (dirichlet_face)
                    dir_data = dirichlet_projection(msh, fc, sides[s]);
                else {
                    auto d = (s == 1 && disk::is_cut(msh, fc)) ? m_fbs : 0;
                    face_LHS_offset = face_SOL_offset(msh, fc) + d;
                }
                for (size_t i = 0; i < m_fbs; i++) {
                    asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet_face) );
                    dirichlet.push_back(dir_data(i));
                }
            }
        }

        return asm_map;
    }

    void
    assemble(const Mesh& msh, size_t q, const macro_cell<Mesh>& Q, disk::location side,
             const matrix_type& lhs, const vector_type& rhs)
    {
        std::vector<T> dirichlet;
        auto asm_map = init_asm_map(msh, q, Q, side, dirichlet);
        assert(asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols());

        for (size_t i = 0; i < lhs.rows(); i++) {
            if (!asm_map[i].assemble())
                continue;
            for (size_t j = 0; j < lhs.cols(); j++) {
                if (asm_map[j].assemble()) {
                    m_triplets.push_back( Eigen::Triplet<T>(size_t(asm_map[i]), size_t(asm_map[j]), lhs(i,j)) );
                    if (m_record_pattern)
                        m_blocks.insert( std::make_pair(block_of(asm_map[i]), block_of(asm_map[j])) );
                }
                else
                    RHS(size_t(asm_map[i])) -= lhs(i,j)*dirichlet[j];
            }
            RHS(size_t(asm_map[i])) += rhs(i);
        }
    }

    void
    finalize(void)
    {
        LHS.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }

    void record_pattern(bool rec) { m_record_pattern = rec; }

    size_t
    block_of(size_t dof) const {
        auto ncd = num_cell_dofs();
        if (dof < ncd)
            return dof / m_cbs;
        return m_num_cells + (dof - ncd) / m_fbs;
    }

    void
    write_block_pattern(const std::string& filename) const {
        std::ofstream ofs(filename);
        ofs << "# cell_blocks " << m_num_cells << " face_blocks " << m_num_other_faces << std::endl;
        for (auto& b : m_blocks)
            ofs << b.first << " " << b.second << " 0" << std::endl;
    }

    vector_type
    take_cell_dofs(size_t q, const macro_cell<Mesh>& Q, const vector_type& solution, disk::location where) const
    {
        return solution.segment(cell_SOL_offset(q, Q, where), m_cbs);
    }
};
