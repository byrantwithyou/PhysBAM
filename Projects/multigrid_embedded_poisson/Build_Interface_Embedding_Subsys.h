//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_INTERFACE_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_INTERFACE_EMBEDDING_SUBSYS_HPP

#include <iosfwd>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>

#include "Build_Embedding_Subsys.h"
#include "EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "Init_Cell_Local_Embedding_Interface_System.h"
#include "INTERFACE_CONSTRAINT_SYSTEM.h"
#include "INTERFACE_INDEX_TRANSFORM.h"
#include "POST_EMBEDDING_INIT_INTERFACE_VISITOR.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_SIGN_OF_CELL_INDEX,
    class T_BETA_NEGATIVE_OF_INDEX, class T_F_NEGATIVE_OF_INDEX,
    class T_BETA_POSITIVE_OF_INDEX, class T_F_POSITIVE_OF_INDEX,
    class T_JUMP_U_OF_X_OF_CELL_INDEX, class T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N_OF_CELL_INDEX
>
void
Build_Interface_Embedding_Subsys(
    const unsigned int n_thread,
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_BETA_NEGATIVE_OF_INDEX& beta_negative_of_index,
    const T_F_NEGATIVE_OF_INDEX& f_negative_of_index,
    const T_BETA_POSITIVE_OF_INDEX& beta_positive_of_index,
    const T_F_POSITIVE_OF_INDEX& f_positive_of_index,
    const T_JUMP_U_OF_X_OF_CELL_INDEX& jump_u_of_x_of_cell_index,
    const T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N_OF_CELL_INDEX& jump_beta_grad_u_dot_n_of_x_and_n_of_cell_index,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    EMBEDDING_UNSTRUCTURED_SUBSYS<T>& embedding_subsys,
    ARRAY<T>& system_rhs,
    INTERFACE_CONSTRAINT_SYSTEM<T,D>& constraint_system,
    ARRAY<T>& constraint_rhs,
    std::ostream& lout = PhysBAM::nout)
{
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;

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

    // Construct multi_index_of_index.
    typedef BOUND_FAST_MEM_FN<
        int (INDEX_TRANSFORM_TYPE::*)( int ) const,
        &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
    > GRID_INDEX_OF_INDEX_TYPE;
    typedef COMPOSE_FUNCTION<
        MULTI_INDEX_BOUND<D>,
        GRID_INDEX_OF_INDEX_TYPE
    > MULTI_INDEX_OF_INDEX_TYPE;
    const MULTI_INDEX_OF_INDEX_TYPE multi_index_of_index(
        multi_index_bound,
        GRID_INDEX_OF_INDEX_TYPE(index_transform)
    );

    // Construct index_of_{nega|posi}tive_multi_index.
    typedef COMPOSE_FUNCTION<
        typename INDEX_TRANSFORM_TYPE::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION,
        MULTI_INDEX_BOUND<D>
    > INDEX_OF_SIGNED_MULTI_INDEX_TYPE;
    const INDEX_OF_SIGNED_MULTI_INDEX_TYPE index_of_negative_multi_index(
        index_transform.Index_Of_Signed_Grid_Index_Function(-1),
        multi_index_bound
    );
    const INDEX_OF_SIGNED_MULTI_INDEX_TYPE index_of_positive_multi_index(
        index_transform.Index_Of_Signed_Grid_Index_Function(+1),
        multi_index_bound
    );

    // Construct {nega|posi}tive_index_transform_stencil_proxy_function.
    typedef INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION<
        MULTI_INDEX_OF_INDEX_TYPE,
        INDEX_OF_SIGNED_MULTI_INDEX_TYPE
    > SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE;
    const SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE
        negative_index_transform_stencil_proxy_function(
            multi_index_of_index,
            index_of_negative_multi_index
        );
    const SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE
        positive_index_transform_stencil_proxy_function(
            multi_index_of_index,
            index_of_positive_multi_index
        );

    // Construct system_stencil_proxy_of_index.
    const BOUND_FAST_MEM_FN<
        typename EMBEDDING_UNSTRUCTURED_SUBSYS<T>::STENCIL_PROXY_TYPE
            (EMBEDDING_UNSTRUCTURED_SUBSYS<T>::*)( int ),
        &EMBEDDING_UNSTRUCTURED_SUBSYS<T>::Stencil_Proxy
    > system_stencil_proxy_of_index(embedding_subsys);

    // Construct constraint_stencil_index_of_cell_multi_index.
    typedef BOUND_FAST_MEM_FN<
        const int& (HASHTABLE<int,int>::*)( const int& ) const,
        &HASHTABLE<int,int>::Get
    > CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE;
    typedef COMPOSE_FUNCTION<
        CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE,
        MULTI_INDEX_BOUND<D>
    > CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE;
    const CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE
        constraint_stencil_index_of_cell_multi_index(
            CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE(
                constraint_system.stencil_index_of_cell_index
            ),
            cell_multi_index_bound
        );
    // Construct constraint_stencil_proxy_of_cell_multi_index.
    typedef BOUND_FAST_MEM_FN<
        typename INTERFACE_CONSTRAINT_SYSTEM<T,D>::STENCIL_PROXY_TYPE
            (INTERFACE_CONSTRAINT_SYSTEM<T,D>::*)( int ),
        &INTERFACE_CONSTRAINT_SYSTEM<T,D>::Stencil_Proxy
    > CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE;
    const COMPOSE_FUNCTION<
        CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE,
        CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE
    > constraint_stencil_proxy_of_cell_multi_index(
        CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE(constraint_system),
        constraint_stencil_index_of_cell_multi_index
    );

    Build_Embedding_Subsys(
        n_thread,
        multi_index_bound,
        sign_of_cell_index,
        Make_Post_Embedding_Init_Interface_Visitor(
            multi_index_bound.Size(),
            embedding_subsys, system_rhs,
            constraint_system, constraint_rhs
        ),
        Make_Init_Cell_Local_Embedding_Interface_System_Visitor(
            min_x, max_x, multi_index_bound,
            phi_of_fine_index,
            beta_negative_of_index, f_negative_of_index,
            beta_positive_of_index, f_positive_of_index,
            jump_u_of_x_of_cell_index, jump_beta_grad_u_dot_n_of_x_and_n_of_cell_index,
            min_dist_to_vertex, sign_of_zero,
            Make_Compose_Function(
                negative_index_transform_stencil_proxy_function,
                system_stencil_proxy_of_index,
                index_of_negative_multi_index
            ),
            Make_Compose_Function(
                Make_Array_Wrapper_Function(system_rhs),
                index_of_negative_multi_index
            ),
            Make_Compose_Function(
                positive_index_transform_stencil_proxy_function,
                system_stencil_proxy_of_index,
                index_of_positive_multi_index
            ),
            Make_Compose_Function(
                Make_Array_Wrapper_Function(system_rhs),
                index_of_positive_multi_index
            ),
            Make_Compose_Function(
                negative_index_transform_stencil_proxy_function,
                constraint_stencil_proxy_of_cell_multi_index
            ),
            Make_Compose_Function(
                positive_index_transform_stencil_proxy_function,
                constraint_stencil_proxy_of_cell_multi_index
            ),
            Make_Compose_Function(
                Make_Array_Wrapper_Function(constraint_rhs),
                constraint_stencil_index_of_cell_multi_index
            )
        ),
        embedding_subsys.index_of_stencil_index,
        constraint_system.cell_index_of_stencil_index,
        lout
    );

    constraint_system.Init_Stencils_Containing_Index();
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_INTERFACE_EMBEDDING_SUBSYS_HPP
