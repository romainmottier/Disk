/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

// NewtonRaphson iteration

#pragma once

#include "diskpp/mechanics/NewtonSolver/GenericIteration.hpp"

namespace disk {

namespace mechanics {

/**
 * @brief Newton-Raphson iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template < typename MeshType >
class NewtonIteration : public GenericIteration< MeshType > {
  private:
    typedef typename GenericIteration< MeshType >::mesh_type mesh_type;
    typedef typename GenericIteration< MeshType >::cell_type cell_type;
    typedef typename GenericIteration< MeshType >::scalar_type scalar_type;

    typedef typename GenericIteration< MeshType >::matrix_type matrix_type;
    typedef typename GenericIteration< MeshType >::vector_type vector_type;

    typedef typename GenericIteration< MeshType >::param_type param_type;
    typedef typename GenericIteration< MeshType >::bnd_type bnd_type;
    typedef typename GenericIteration< MeshType >::behavior_type behavior_type;

    typedef typename GenericIteration< MeshType >::elem_type elem_type;
    typedef typename GenericIteration< MeshType >::func_type func_type;

  public:
    NewtonIteration( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                     const MeshDegreeInfo< mesh_type > &degree_infos,
                     const std::shared_ptr< solvers::LinearSolver< scalar_type > > lin_solv,
                     const TimeStep< scalar_type > &current_step )
        : GenericIteration< MeshType >( msh, bnd, rp, degree_infos, lin_solv, current_step ) {}

    AssemblyInfo assemble( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                           const MeshDegreeInfo< mesh_type > &degree_infos,
                           const std::unique_ptr< func_type > &lf,
                           const std::vector< matrix_type > &gradient_precomputed,
                           const std::vector< matrix_type > &stab_precomputed,
                           behavior_type &behavior, StabCoeffManager< scalar_type > &stab_manager,
                           MultiTimeField< scalar_type > &fields ) override {
        elem_type elem;
        AssemblyInfo ai;

        // set RHS to zero
        this->m_assembler.initialize();
        this->m_F_int = 0.0;

        const bool small_def = ( behavior.getDeformation() == SMALL_DEF );

        // Like if it is an implicit scheme
        auto current_time = this->m_time_step.end_time();
        auto depl = fields.getCurrentField( FieldName::DEPL );
        auto depl_faces = fields.getCurrentField( FieldName::DEPL_FACES );

        std::vector< vector_type > resi_cells;

        std::vector< vector_type > acce_cells;
        if ( this->m_dyna.enable() ) {
            acce_cells = fields.getCurrentField( FieldName::ACCE_CELLS );
            resi_cells.reserve( msh.cells_size() );
        }

        auto rlf = this->_getLoad( lf, current_time );

        timecounter tc, ttot;

        ttot.tic();

        for ( auto &cl : msh ) {
            const auto cell_i = msh.lookup( cl );

            const auto huT = depl.at( cell_i );

            const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
            const auto num_cell_dofs = vector_cell_dofs( msh, cell_infos );

            const auto num_tot_dofs = huT.size();
            const auto num_faces_dofs = num_tot_dofs - num_cell_dofs;

            // Gradient Reconstruction
            // std::cout << "Grad" << std::endl;
            tc.tic();
            matrix_type GT = _gradrec( msh, cl, rp, degree_infos, small_def, gradient_precomputed );
            tc.toc();
            ai.m_time_gradrec += tc.elapsed();

            // Mechanical Computation

            tc.tic();
            // std::cout << "Elem" << std::endl;
            elem.compute( msh, cl, bnd, rp, degree_infos, rlf, GT, huT, this->m_time_step, behavior,
                          stab_manager, small_def );

            matrix_type lhs = elem.K_int;
            vector_type rhs = elem.RTF;

            if ( !this->m_dyna.isExplicit() ) {
                this->m_F_int += elem.F_int.squaredNorm();
            } else {
                this->m_F_int += elem.F_int.tail( num_faces_dofs ).squaredNorm();
            }

            tc.toc();
            ai.m_time_elem += tc.elapsed();

            // Stabilisation Contribution
            // std::cout << "Stab" << std::endl;
            tc.tic();
            if ( rp.m_stab ) {
                const auto beta_s = stab_manager.getValue( msh, cl );

                matrix_type stab = beta_s * _stab( msh, cl, rp, degree_infos, stab_precomputed );

                // std::cout << beta_s << std::endl;

                lhs += stab;
                rhs -= stab * huT;
            }
            tc.toc();
            ai.m_time_stab += tc.elapsed();

            // Dynamic contribution
            if ( this->m_dyna.enable() ) {
                this->m_dyna.compute( msh, cl, degree_infos, huT, acce_cells.at( cell_i ),
                                      this->m_time_step );
                if ( !this->m_dyna.isExplicit() ) {
                    lhs.topLeftCorner( num_cell_dofs, num_cell_dofs ) += this->m_dyna.K_iner;
                    rhs.head( num_cell_dofs ) += this->m_dyna.R_iner;
                }
                ai.m_time_dyna += this->m_dyna.time_dyna;
            }

            // std::cout << "R: " << rhs.norm() << std::endl;
            // std::cout << rhs.transpose() << std::endl;

            // Static Condensation
            // std::cout << "StatCond" << std::endl;

            if ( this->m_dyna.isExplicit() ) {
                tc.tic();

                const auto num_cell_dofs = acce_cells.at( cell_i ).size();
                const auto num_tot_dofs = lhs.rows();
                const auto num_faces_dofs = num_tot_dofs - num_cell_dofs;

                this->m_AL[cell_i] = matrix_type::Zero( num_cell_dofs, num_faces_dofs );
                this->m_bL[cell_i] = vector_type::Zero( num_cell_dofs );

                resi_cells.push_back( rhs.head( num_cell_dofs ) );

                tc.toc();
                ai.m_time_statcond += tc.elapsed();

                tc.tic();
                this->m_assembler.assemble_nonlinear(
                    msh, cl, bnd, lhs.bottomRightCorner( num_faces_dofs, num_faces_dofs ),
                    rhs.tail( num_faces_dofs ), depl_faces );
                tc.toc();
                ai.m_time_assembler += tc.elapsed();
            } else {
                tc.tic();

                const auto scnp = make_vector_static_condensation_withMatrix( msh, cl, degree_infos,
                                                                              lhs, rhs, true );

                this->m_AL[cell_i] = std::get< 1 >( scnp );
                this->m_bL[cell_i] = std::get< 2 >( scnp );

                tc.toc();
                ai.m_time_statcond += tc.elapsed();

                const auto &lc = std::get< 0 >( scnp );
                tc.tic();
                this->m_assembler.assemble_nonlinear( msh, cl, bnd, lc.first, lc.second,
                                                      depl_faces );
                tc.toc();
                ai.m_time_assembler += tc.elapsed();
            }
        }

        if ( this->m_dyna.isExplicit() ) {
            fields.setCurrentField( FieldName::RESI_CELLS, resi_cells );
        }

        this->m_F_int = sqrt( this->m_F_int );
        // std::cout << "F_int: " << m_F_int << std::endl;

        ai.m_time_law += elem.time_law;
        ai.m_time_contact += elem.time_contact;
        ai.m_time_load += elem.time_load;
        ai.m_time_rigi += elem.time_rigi;
        ai.m_time_fint += elem.time_fint;

        tc.tic();
        this->m_assembler.impose_neumann_boundary_conditions( msh, bnd );
        this->m_assembler.finalize();
        tc.toc();
        ai.m_time_assembler += tc.elapsed();

        ttot.toc();
        ai.m_time_assembly = ttot.elapsed();
        ai.m_linear_system_size = this->m_assembler.LHS.rows();
        return ai;
    }

    scalar_type postprocess( const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                             const MeshDegreeInfo< mesh_type > &degree_infos,
                             const std::unique_ptr< func_type > &lf,
                             const std::vector< matrix_type > &gradient_precomputed,
                             const std::vector< matrix_type > &stab_precomputed,
                             behavior_type &behavior, StabCoeffManager< scalar_type > &stab_manager,
                             MultiTimeField< scalar_type > &fields ) override {
        timecounter tc;
        tc.tic();

        auto depl_faces = fields.getCurrentField( FieldName::DEPL_FACES );
        auto depl_cells = fields.getCurrentField( FieldName::DEPL_CELLS );
        auto depl = fields.getCurrentField( FieldName::DEPL );

        const auto [dudT, idx] = this->m_assembler.expand_solution_nonlinear(
            msh, bnd, this->m_system_displ, depl_faces );

        // Update cell
        for ( auto &cl : msh ) {
            const auto cell_i = msh.lookup( cl );

            const auto cell_infos = degree_infos.cellDegreeInfo( msh, cl );
            const auto faces_infos = cell_infos.facesDegreeInfo();
            const auto num_faces_dofs = vector_faces_dofs( msh, faces_infos );

            vector_type xdT = vector_type( num_faces_dofs );

            const auto fcs_id = faces_id( msh, cl );
            size_t face_offset = 0;
            for ( size_t face_i = 0; face_i < fcs_id.size(); face_i++ ) {
                const size_t face_id = fcs_id[face_i];
                const auto n_face_dofs = idx( face_id + 1 ) - idx( face_id );

                xdT.segment( face_offset, n_face_dofs ) =
                    dudT.segment( idx( face_id ), n_face_dofs );

                face_offset += n_face_dofs;
            }

            // static decondensation
            const vector_type xT = this->m_bL[cell_i] - this->m_AL[cell_i] * xdT;

            // Update element U^{i+1} = U^i + delta U^i
            depl.at( cell_i ).head( xT.size() ) += xT;
            depl.at( cell_i ).tail( xdT.size() ) += xdT;
            depl_cells.at( cell_i ) += xT;

            // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
            // std::cout << "sol_F" << std::endl;
            // std::cout << xdT.transpose() << std::endl;
            // std::cout << "ft" << std::endl;
            // std::cout << m_bL[cell_i].transpose() << std::endl;
            // std::cout << "sol_T" << std::endl;
            // std::cout << xT.transpose() << std::endl;
            // std::cout << depl.at(cell_i).transpose() << std::endl;
        }
        fields.setCurrentField( FieldName::DEPL, depl );
        fields.setCurrentField( FieldName::DEPL_CELLS, depl_cells );

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i

        auto depl_faces_new = fields.getCurrentField( FieldName::DEPL_FACES );
        int face_i = 0;
        for ( auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++ ) {
            const auto fc = *itor;
            const size_t face_id = msh.lookup( fc );

            depl_faces_new.at( face_i++ ) +=
                dudT.segment( idx( face_id ), idx( face_id + 1 ) - idx( face_id ) );
        }
        fields.setCurrentField( FieldName::DEPL_FACES, depl_faces_new );

        this->m_dyna.postprocess( msh, this->m_time_step, fields );

        tc.toc();
        return tc.elapsed();
    }
};
} // namespace mechanics
} // namespace disk