//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PROJECT_TO_ZERO_DIVERGENCE_WITH_INTERFACE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PROJECT_TO_ZERO_DIVERGENCE_WITH_INTERFACE_HPP

#include <limits>

#include <boost/foreach.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/ref.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/CONSTANT_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/ASSIGN_SIGN_TO_INDEX_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Stencils/ZERO_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Aggregate_Constraints.h"
#include "AGGREGATE_CONSTRAINT_SYSTEM.h"
#include "Build_Interface_Embedding_Subsys.h"
#include "Build_Interface_Regular_Subsys.h"
#include "Build_ZTAZ_Embedding_Subsys.h"
#include "DIVERGENCE_OF_MAC_VECTOR_FIELD.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "Init_ZTAZ_Embedding.h"
#include "INTERFACE_CONSTRAINT_SYSTEM.h"
#include "INTERFACE_INDEX_TRANSFORM.h"
#include "Select_Indys.h"
#include "SYSTEM_SUM.h"

#ifndef PHYSBAM_NO_PETSC
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include "Petsc/Solve_SPD_System_With_ICC_PCG.h"
#include <petsc.h>
#endif // #ifndef PHYSBAM_NO_PETSC

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_BETA_NEGATIVE_OF_INDEX,
    class T_BETA_POSITIVE_OF_INDEX,
    class T_JUMP_P_OF_X_OF_CELL_INDEX,
    class T_JUMP_BETA_GRAD_P_DOT_N_OF_X_AND_N_OF_CELL_INDEX,
    class T_MAC_VECTOR_FIELD
>
void Project_To_Zero_Divergence_With_Interface(
    const unsigned int n_thread,
    const VECTOR<T,D> min_x, const VECTOR<T,D> max_x,
    const MULTI_INDEX_BOUND<D> mac_cell_multi_index_bound,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_BETA_NEGATIVE_OF_INDEX& beta_negative_of_index,
    const T_BETA_POSITIVE_OF_INDEX& beta_positive_of_index,
    const T_JUMP_P_OF_X_OF_CELL_INDEX& jump_p_of_x_of_cell_index,
    const T_JUMP_BETA_GRAD_P_DOT_N_OF_X_AND_N_OF_CELL_INDEX& jump_beta_grad_p_dot_n_of_x_and_n_of_cell_index,
    T_MAC_VECTOR_FIELD mac_vector_field)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef DOMAIN_REGULAR_CROSS_SUBSYS<T,D> REGULAR_SUBSYS_TYPE;
    typedef EMBEDDING_UNSTRUCTURED_SUBSYS<T> EMBEDDING_SUBSYS_TYPE;
    typedef SYSTEM_SUM< boost::mpl::vector2< REGULAR_SUBSYS_TYPE&, EMBEDDING_SUBSYS_TYPE& > > SYSTEM_TYPE;
    typedef INTERFACE_CONSTRAINT_SYSTEM<T,D> CONSTRAINT_SYSTEM_TYPE;

    typedef EMBEDDING_UNSTRUCTURED_SUBSYS<T> ZTAZ_EMBEDDING_SUBSYS_TYPE;
    typedef SYSTEM_SUM< boost::mpl::vector2< REGULAR_SUBSYS_TYPE&, ZTAZ_EMBEDDING_SUBSYS_TYPE& > > ZTAZ_SYSTEM_TYPE;
    typedef AGGREGATE_CONSTRAINT_SYSTEM<T,D> AGGREGATE_CONSTRAINT_SYSTEM_TYPE;

    const MULTI_INDEX_BOUND<D> multi_index_bound = mac_cell_multi_index_bound;
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    const DIVERGENCE_OF_MAC_VECTOR_FIELD< T, D, T_MAC_VECTOR_FIELD >
        divergence_of_mac_vector_field(dx, mac_vector_field);

    REGULAR_SUBSYS_TYPE regular_subsys(multi_index_bound, dx);
    EMBEDDING_SUBSYS_TYPE embedding_subsys;
    SYSTEM_TYPE system(boost::fusion::make_vector(
        boost::ref(regular_subsys),
        boost::ref(embedding_subsys)
    ));
    CONSTRAINT_SYSTEM_TYPE constraint_system;
    ARRAY<T> system_rhs;
    ARRAY<T> constraint_rhs;

    Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT<2>(
        n_thread,
        cell_multi_index_bound,
        Make_Compose_Function(SIGN_FUNCTION(), phi_of_fine_index),
        Make_Assign_Sign_To_Index_Grid_Visitor(
            Make_Array_Wrapper_Function(regular_subsys.sign_of_cell_index)
        ),
        -1 // sign_of_zero
    );

    Build_Interface_Embedding_Subsys(
        n_thread,
        min_x, max_x, multi_index_bound,
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        phi_of_fine_index,
        beta_negative_of_index,
        divergence_of_mac_vector_field,
        beta_positive_of_index,
        divergence_of_mac_vector_field,
        jump_p_of_x_of_cell_index, jump_beta_grad_p_dot_n_of_x_and_n_of_cell_index,
        0.0f, // min_dist_to_vertex
        -1, // sign_of_zero
        embedding_subsys, system_rhs,
        constraint_system, constraint_rhs
    );

    Build_Interface_Regular_Subsys(
        n_thread,
        dx, cell_multi_index_bound,
        beta_negative_of_index,
        divergence_of_mac_vector_field,
        beta_positive_of_index,
        divergence_of_mac_vector_field,
        regular_subsys, As_Array_View(system_rhs)
    );

    const int n_embedding  = embedding_subsys.stencils.Size();
    const int n_virtual    = n_embedding / 2;
    const int n_index      = multi_index_bound.Size() + n_virtual;
    const int n_constraint = constraint_system.stencils.Size();

    // Construct index_transform.
    typedef BOUND_FAST_MEM_FN<
        const int& (HASHTABLE<int,int>::*)( const int& ) const,
        &HASHTABLE<int,int>::Get
    > VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE;
    typedef ARRAY_VIEW<const int> GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE;
    typedef typename Result_Of::MAKE_COMPOSE_FUNCTION<
        SIGN_FUNCTION,
        T_PHI_OF_FINE_INDEX,
        FINE_MULTI_INDEX_FUNCTION<2>,
        MULTI_INDEX_BOUND<D>
    >::type SIGN_OF_GRID_INDEX_TYPE;
    typedef INTERFACE_INDEX_TRANSFORM<
        VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE,
        GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE,
        SIGN_OF_GRID_INDEX_TYPE
    > INDEX_TRANSFORM_TYPE;
    const INDEX_TRANSFORM_TYPE index_transform(
        multi_index_bound.Size(),
        VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE(embedding_subsys.stencil_index_of_index),
        GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE(embedding_subsys.index_of_stencil_index),
        Make_Compose_Function(
            SIGN_FUNCTION(),
            phi_of_fine_index,
            FINE_MULTI_INDEX_FUNCTION<2>(),
            multi_index_bound
        )
    );

    BOUND_FAST_MEM_FN<
        int (INDEX_TRANSFORM_TYPE::*)( int ) const,
        &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
    > grid_index_of_index(index_transform);

    ARRAY<int> indy_index_of_constraint_index(n_constraint, false); // uninit'ed
    AGGREGATE_CONSTRAINT_SYSTEM_TYPE aggregate_constraint_system;
    aggregate_constraint_system.index_of_indy_index.Preallocate(n_embedding / (2 * (1 << D)));
    Select_Indys(
        cell_multi_index_bound,
        Make_Compose_Function(multi_index_bound, grid_index_of_index),
        constraint_system.stencil_index_of_cell_index,
        BOUND_FAST_MEM_FN<
            typename CONSTRAINT_SYSTEM_TYPE::CONST_STENCIL_PROXY_TYPE
                (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
            &CONSTRAINT_SYSTEM_TYPE::Stencil_Proxy
        >(constraint_system),
        PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE(
            index_transform,
            &INDEX_TRANSFORM_TYPE::Index_Is_Virtual
        ),
        0.0f, // min_relative_indy_weight
        As_Array_View(indy_index_of_constraint_index),
        aggregate_constraint_system.index_of_indy_index
    );
    const int n_indy = aggregate_constraint_system.index_of_indy_index.Size();

    aggregate_constraint_system.Init_Indy_Index_Of_Index();
    aggregate_constraint_system.value_of_indy_index.Exact_Resize(n_indy); // init'ed to 0
    aggregate_constraint_system.indyless_stencils.Exact_Resize(n_indy);
    Aggregate_Constraints(
        BOUND_FAST_MEM_FN<
            typename CONSTRAINT_SYSTEM_TYPE::CONST_STENCIL_PROXY_TYPE
                (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
            &CONSTRAINT_SYSTEM_TYPE::Stencil_Proxy
        >(constraint_system),
        As_Const_Array_View(indy_index_of_constraint_index),
        As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
        As_Array_View(aggregate_constraint_system.value_of_indy_index),
        BOUND_FAST_MEM_FN<
            typename AGGREGATE_CONSTRAINT_SYSTEM_TYPE::INDYLESS_STENCIL_PROXY_TYPE
                (AGGREGATE_CONSTRAINT_SYSTEM_TYPE::*)( int ),
            &AGGREGATE_CONSTRAINT_SYSTEM_TYPE::Indyless_Stencil_Proxy
        >(aggregate_constraint_system)
    );
    aggregate_constraint_system.Init_Stencils_Containing_Index();
    ARRAY<T> aggregate_constraint_rhs(n_indy); // init'ed to 0
    for(int constraint_index = 1; constraint_index <= n_constraint; ++constraint_index) {
        const int indy_index = indy_index_of_constraint_index(constraint_index);
        aggregate_constraint_rhs(indy_index) += constraint_rhs(constraint_index);
    }

    ZTAZ_EMBEDDING_SUBSYS_TYPE ztaz_embedding_subsys;
    Init_ZTAZ_Embedding(
        embedding_subsys.stencil_index_of_index,
        As_Const_Array_View(embedding_subsys.index_of_stencil_index),
        aggregate_constraint_system.indy_index_of_index,
        As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
        Make_Constant_Function(ZERO_STENCIL_PROXY<int,T>()), // okay, since all indys are virtual
        ztaz_embedding_subsys.stencil_index_of_index,
        ztaz_embedding_subsys.index_of_stencil_index
    );
    const int n_ztaz_embedding = ztaz_embedding_subsys.stencil_index_of_index.Size();
    ztaz_embedding_subsys.stencils.Exact_Resize(n_ztaz_embedding);
    Build_ZTAZ_Embedding_Subsys(
        Make_Constant_Function(ZERO_STENCIL_PROXY<int,T>()), // okay, since all indys are virtual
        BOUND_FAST_MEM_FN<
            typename EMBEDDING_SUBSYS_TYPE::CONST_STENCIL_PROXY_TYPE
                (EMBEDDING_SUBSYS_TYPE::*)( int ) const,
            &EMBEDDING_SUBSYS_TYPE::Stencil_Proxy
        >(embedding_subsys),
        aggregate_constraint_system.indy_index_of_index,
        As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
        aggregate_constraint_system.stencils_containing_index,
        As_Const_Array_View(aggregate_constraint_system.value_of_indy_index),
        BOUND_FAST_MEM_FN<
            typename AGGREGATE_CONSTRAINT_SYSTEM_TYPE::CONST_INDYLESS_STENCIL_PROXY_TYPE
                (AGGREGATE_CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
            &AGGREGATE_CONSTRAINT_SYSTEM_TYPE::Indyless_Stencil_Proxy
        >(aggregate_constraint_system),
        As_Const_Array_View(ztaz_embedding_subsys.index_of_stencil_index),
        BOUND_FAST_MEM_FN<
            typename ZTAZ_EMBEDDING_SUBSYS_TYPE::STENCIL_PROXY_TYPE
                (ZTAZ_EMBEDDING_SUBSYS_TYPE::*)( int ),
            &ZTAZ_EMBEDDING_SUBSYS_TYPE::Stencil_Proxy_Of_Stencil_Index
        >(ztaz_embedding_subsys)
    );

    // f <- Z^T*(f - A*c)
    ARRAY<T>& ztaz_system_rhs = system_rhs;
    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
        const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
        const T c_value = aggregate_constraint_rhs(indy_index)
                        / aggregate_constraint_system.value_of_indy_index(indy_index);
        system.Apply_Transpose(index, ztaz_system_rhs, -c_value);
    }
    aggregate_constraint_system.Apply_Z_Transpose(ztaz_system_rhs);

    const ZTAZ_SYSTEM_TYPE ztaz_system(boost::fusion::make_vector(
        boost::ref(regular_subsys),
        boost::ref(ztaz_embedding_subsys))
    );

    ARRAY<T> p(n_index); // init'ed to 0

    // Solve for p.
#ifdef PHYSBAM_NO_PETSC
#else // #ifdef PHYSBAM_NO_PETSC
    PHYSBAM_PETSC_CALL_AND_CHKERRQ((
        Petsc::Solve_SPD_System_With_ICC_PCG<T,D>(
            n_thread,
            ztaz_system, ztaz_system_rhs,
            true,                                       // has_constant_vectors_in_null_space
            std::numeric_limits< unsigned int >::max(), // max_iterations
            1e-8f,                                      // relative_tolerance
            std::numeric_limits< float >::min(),        // absolute_tolerance
            false,                                      // print_residuals
            true,                                       // precondition
            p
        )
    ));
#endif // #ifdef PHYSBAM_NO_PETSC

    // p <- c + Z*p
    aggregate_constraint_system.Apply_Z(p);
    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
        const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
        const T c_value = aggregate_constraint_rhs(indy_index)
                        / aggregate_constraint_system.value_of_indy_index(indy_index);
        p(index) += c_value;
    }

    // v <- v + grad(p)
    {
        MULTI_INDEX_BOUND<D> clipped_multi_index_bound = multi_index_bound;
        for(int d = 1; d <= D; ++d) {
            --clipped_multi_index_bound.max_multi_index[d];
            // TODO: MT
            BOOST_FOREACH( const MULTI_INDEX_TYPE multi_index1, clipped_multi_index_bound ) {
                MULTI_INDEX_TYPE multi_index2 = multi_index1;
                ++multi_index2[d];
                MULTI_INDEX_TYPE fine_multi_index = 2 * multi_index1 - 1;
                ++fine_multi_index[d];
                int sign = SIGN_FUNCTION()(phi_of_fine_index(fine_multi_index));
                if(sign == 0)
                    sign = -1;
                const int grid_index1 = multi_index_bound.Linear_Index(multi_index1);
                const int grid_index2 = multi_index_bound.Linear_Index(multi_index2);
                const int index1 = index_transform.Index_Of_Signed_Grid_Index(grid_index1, sign);
                const int index2 = index_transform.Index_Of_Signed_Grid_Index(grid_index2, sign);
                const T dp = (p(index2) - p(index1)) / dx[d];
                const T beta1 = sign == -1 ? beta_negative_of_index(multi_index1) : beta_positive_of_index(multi_index1);
                const T beta2 = sign == -1 ? beta_negative_of_index(multi_index2) : beta_positive_of_index(multi_index2);
                const T beta = (beta1 + beta2) / 2;
                mac_vector_field(d, multi_index2) += beta * dp;
            }
            ++clipped_multi_index_bound.max_multi_index[d];
        }
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PROJECT_TO_ZERO_DIVERGENCE_WITH_INTERFACE_HPP
