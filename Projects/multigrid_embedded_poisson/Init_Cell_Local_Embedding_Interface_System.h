//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP

#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/Pre_Compose_Translate.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FLUX_FUNCTION_TIMES_MULTI_POWERS_OVER_BOUNDARY_FOR_POISSON_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Math/Eval_All_Multilinear_Basis_From_Monomial_Evals.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Evaluate_Cell_Local_Integrations.h"
#include "Init_Cell_Local_Embedding_Domain_System.h"
#include "PRE_COMPOSE_TRANSLATE_X.h"

//#define PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_NEGATIVE_BETA_OF_INDEX,
    class T_NEGATIVE_F_OF_INDEX,
    class T_POSITIVE_BETA_OF_INDEX,
    class T_POSITIVE_F_OF_INDEX,
    class T_JUMP_U_OF_X,
    class T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
    class T_SYSTEM_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_RHS_OF_INDEX,
    class T_CONSTRAINT_STENCIL_PROXY
>
void Init_Cell_Local_Embedding_Interface_System(
    const VECTOR<T,D> min_x, const VECTOR<T,D> max_x,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_NEGATIVE_BETA_OF_INDEX& negative_beta_of_index,
    const T_NEGATIVE_F_OF_INDEX& negative_f_of_index,
    const T_POSITIVE_BETA_OF_INDEX& positive_beta_of_index,
    const T_POSITIVE_F_OF_INDEX& positive_f_of_index,
    const T_JUMP_U_OF_X& jump_u_of_x,
    const T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N& jump_beta_grad_u_dot_n_of_x_and_n,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    const VECTOR<int,D> cell_multi_index,
    const T_SYSTEM_STENCIL_PROXY_OF_INDEX& negative_system_stencil_proxy_of_index,
    const T_SYSTEM_STENCIL_PROXY_OF_INDEX& positive_system_stencil_proxy_of_index,
    const T_SYSTEM_RHS_OF_INDEX& negative_system_rhs_of_index,
    const T_SYSTEM_RHS_OF_INDEX& positive_system_rhs_of_index,
    T_CONSTRAINT_STENCIL_PROXY negative_constraint_stencil_proxy,
    T_CONSTRAINT_STENCIL_PROXY positive_constraint_stencil_proxy,
    T& constraint_rhs)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const VECTOR<T,D> x0 = Multi_Index_X(min_x, max_x, 2 * multi_index_bound - 1, 2 * cell_multi_index);
    const VECTOR<T,D> dx = (max_x - min_x) / (multi_index_bound.max_multi_index - 1);
    const VECTOR<T,D> dx_over_2 = dx / 2;
    const T dv = dx.Product();
    static_cast<void>(dx_over_2);

    VECTOR< T, STATIC_POW_C<3,D>::value > monomials_integrated_over_negative_volume; // init'ed to 0
    VECTOR< T, STATIC_POW_C<2,D>::value > monomials_integrated_over_boundary; // init'ed to 0
    VECTOR< T, (1 << D) > multilinear_basis_integrated_over_negative_volume; // init'ed to 0
    VECTOR< T, (1 << D) > multilinear_basis_integrated_over_boundary; // init'ed to 0
    VECTOR< VECTOR< T, (1 << D) >, (1 << D) >
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_negative_volume; // init'ed to 0

    T jump_u_integrated_over_boundary = 0;
#ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
    T average_jump_beta_grad_u_dot_n = 0;
#else // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
    VECTOR< T, STATIC_POW_C<2,D>::value > jump_beta_grad_u_dot_n_times_monomials_integrated_over_boundary; // init'ed to 0
#endif // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR

    Evaluate_Cell_Local_Integrations(
        dx,
        -1,
        phi_of_fine_index,
        min_dist_to_vertex,
        sign_of_zero,
        Make_Visitor_Sequence(
            Make_Integrate_Function_Over_Boundary_Polytope_Visitor<
                GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,D,2>
            >(
                -1,
                Pre_Compose_Translate(jump_u_of_x, x0),
                jump_u_integrated_over_boundary
            ),
#ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
            Make_Integrate_Flux_Function_Over_Boundary_Polytope_Visitor<
                GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,D,2>
            >(
                -1,
                Make_Pre_Compose_Translate_X(jump_beta_grad_u_dot_n_of_x_and_n, x0),
                average_jump_beta_grad_u_dot_n
            )
#else // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
            Make_Integrate_Flux_Function_Times_Multi_Powers_Over_Boundary_For_Poisson_Polytope_Visitor<2>(
                -1,
                Make_Pre_Compose_Translate_X(jump_beta_grad_u_dot_n_of_x_and_n, x0),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(jump_beta_grad_u_dot_n_times_monomials_integrated_over_boundary),
                    STATIC_MULTI_INDEX_CUBE<D,0,1>()
                )
            )
        ),
#endif // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
        cell_multi_index,
        monomials_integrated_over_negative_volume,
        monomials_integrated_over_boundary,
        multilinear_basis_integrated_over_negative_volume,
        multilinear_basis_integrated_over_boundary,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_negative_volume
    );

    const T negative_volume_measure = monomials_integrated_over_negative_volume[1];
    const T positive_volume_measure = dv - negative_volume_measure;

    const VECTOR< T, (1 << D) > multilinear_basis_integrated_over_positive_volume =
        (dv / (1 << D)) * VECTOR< T, (1 << D) >::All_Ones_Vector()
      - multilinear_basis_integrated_over_negative_volume;
    VECTOR< VECTOR< T, (1 << D) >, (1 << D) >
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_positive_volume; // init'ed to 0
    {
        VECTOR< T, (1 << D) > grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume; // init'ed to 0
        for(int k = 1; k <= (1 << D); ++k) {
            const MULTI_INDEX_TYPE mk = STATIC_MULTI_INDEX_CUBE<D,0,1>::Multi_Index(k);
            for(int d = 1; d <= D; ++d)
                grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[k] +=
                    (1 - 2*mk[d]) * (1 + mk[d]) / (dx[d] * dx[d]);
            grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[k] *=
                dv / (STATIC_POW_C<3,D-1>::value * (1 << mk.Sum()));
        }
        for(int i = 1; i <= (1 << D); ++i) {
            for(int j = 1; j <= (1 << D); ++j) {
                const int k = 1 + ((i - 1) ^ (j - 1));
                grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_positive_volume[i][j] =
                    grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[k]
                  - grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_negative_volume[i][j];
            }
        }
    }

    Init_Cell_Local_Embedding_Domain_System(
        multi_index_bound,
        negative_beta_of_index,
        negative_f_of_index,
        cell_multi_index,
        negative_volume_measure,
        multilinear_basis_integrated_over_negative_volume,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_negative_volume,
        negative_system_stencil_proxy_of_index,
        negative_system_rhs_of_index
    );

    Init_Cell_Local_Embedding_Domain_System(
        multi_index_bound,
        positive_beta_of_index,
        positive_f_of_index,
        cell_multi_index,
        positive_volume_measure,
        multilinear_basis_integrated_over_positive_volume,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_positive_volume,
        positive_system_stencil_proxy_of_index,
        positive_system_rhs_of_index
    );

    const MULTI_INDEX_CUBE<D,0,1> cell_local_multi_index_cube(cell_multi_index);
#ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
    const T boundary_measure = monomials_integrated_over_boundary[1];
    assert(boundary_measure >= 0);
    if(boundary_measure == 0) {
        assert(average_jump_beta_grad_u_dot_n == 0);
        return;
    }
    average_jump_beta_grad_u_dot_n /= boundary_measure;
    for(int i = 1; i <= (1 << D); ++i) {
        const MULTI_INDEX_TYPE multi_index_i = cell_local_multi_index_cube.Multi_Index(i);
        negative_system_rhs_of_index(multi_index_i) -=
            average_jump_beta_grad_u_dot_n * multilinear_basis_integrated_over_boundary[i] / 2;
        positive_system_rhs_of_index(multi_index_i) -=
            average_jump_beta_grad_u_dot_n * multilinear_basis_integrated_over_boundary[i] / 2;
    }
#else // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR
    const VECTOR< T, (1 << D) > jump_beta_grad_u_dot_n_times_multilinear_basis_integrated_over_boundary =
        Eval_All_Multilinear_Basis_From_Monomial_Evals(
            dx_over_2,
            Make_Compose_Function(
                Make_Array_Wrapper_Function(jump_beta_grad_u_dot_n_times_monomials_integrated_over_boundary),
                STATIC_MULTI_INDEX_CUBE<D,0,1>()
            )
        );
    for(int i = 1; i <= (1 << D); ++i) {
        const MULTI_INDEX_TYPE multi_index_i = cell_local_multi_index_cube.Multi_Index(i);
        negative_system_rhs_of_index(multi_index_i) -=
            jump_beta_grad_u_dot_n_times_multilinear_basis_integrated_over_boundary[i] / 2;
        positive_system_rhs_of_index(multi_index_i) -=
            jump_beta_grad_u_dot_n_times_multilinear_basis_integrated_over_boundary[i] / 2;
    }
#endif // #ifdef PHYSBAM_EMBEDDED_POISSON_V2_USE_B_BAR

    CUBE_STENCIL<T,D,0,1> negative_constraint_stencil;
    CUBE_STENCIL<T,D,0,1> positive_constraint_stencil;
    for(int i = 1; i <= (1 << D); ++i) {
        negative_constraint_stencil.values[i-1] = -multilinear_basis_integrated_over_boundary[i];
        positive_constraint_stencil.values[i-1] = +multilinear_basis_integrated_over_boundary[i];
    }
    negative_constraint_stencil_proxy +=
        CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, negative_constraint_stencil);
    positive_constraint_stencil_proxy +=
        CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, positive_constraint_stencil);
    constraint_rhs = jump_u_integrated_over_boundary;
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP
