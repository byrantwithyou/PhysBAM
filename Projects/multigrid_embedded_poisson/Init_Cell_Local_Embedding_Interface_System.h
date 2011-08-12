//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/Pre_Compose_Translate.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FLUX_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FLUX_FUNCTION_TIMES_MULTI_POWERS_OVER_BOUNDARY_FOR_POISSON_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
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

//#define PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_BETA_NEGATIVE_OF_INDEX, class T_F_NEGATIVE_OF_INDEX,
    class T_BETA_POSITIVE_OF_INDEX, class T_F_POSITIVE_OF_INDEX,
    class T_JUMP_U_OF_X, class T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
    class T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_NEGATIVE_RHS_OF_INDEX,
    class T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_POSITIVE_RHS_OF_INDEX,
    class T_CONSTRAINT_NEGATIVE_STENCIL_PROXY,
    class T_CONSTRAINT_POSITIVE_STENCIL_PROXY
>
void Init_Cell_Local_Embedding_Interface_System(
    const VECTOR<T,D> dx,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_BETA_NEGATIVE_OF_INDEX& beta_negative_of_index,
    const T_F_NEGATIVE_OF_INDEX& f_negative_of_index,
    const T_BETA_POSITIVE_OF_INDEX& beta_positive_of_index,
    const T_F_POSITIVE_OF_INDEX& f_positive_of_index,
    const T_JUMP_U_OF_X& jump_u_of_x,
    const T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N& jump_beta_grad_u_dot_n_of_x_and_n,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    const VECTOR<int,D> cell_multi_index,
    const VECTOR<T,D> x0,
    const T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX& system_negative_stencil_proxy_of_index,
    const T_SYSTEM_NEGATIVE_RHS_OF_INDEX& system_negative_rhs_of_index,
    const T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX& system_positive_stencil_proxy_of_index,
    const T_SYSTEM_POSITIVE_RHS_OF_INDEX& system_positive_rhs_of_index,
    T_CONSTRAINT_NEGATIVE_STENCIL_PROXY constraint_negative_stencil_proxy,
    T_CONSTRAINT_POSITIVE_STENCIL_PROXY constraint_positive_stencil_proxy,
    T& constraint_rhs)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

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
#ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
    T average_jump_beta_grad_u_dot_n = 0;
#else // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
    VECTOR< T, STATIC_POW_C<2,D>::value > jump_beta_grad_u_dot_n_times_monomials_integrated_over_boundary; // init'ed to 0
#endif // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR

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
#ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
            Make_Integrate_Flux_Function_Over_Boundary_Polytope_Visitor<
                GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,D,2>
            >(
                -1,
                Make_Pre_Compose_Translate_X(jump_beta_grad_u_dot_n_of_x_and_n, x0),
                average_jump_beta_grad_u_dot_n
            )
#else // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
            Make_Integrate_Flux_Function_Times_Multi_Powers_Over_Boundary_For_Poisson_Polytope_Visitor<2>(
                -1,
                Make_Pre_Compose_Translate_X(jump_beta_grad_u_dot_n_of_x_and_n, x0),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(jump_beta_grad_u_dot_n_times_monomials_integrated_over_boundary),
                    STATIC_MULTI_INDEX_CUBE<D,0,1>()
                )
            )
        ),
#endif // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
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
        beta_negative_of_index,
        f_negative_of_index,
        cell_multi_index,
        negative_volume_measure,
        multilinear_basis_integrated_over_negative_volume,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_negative_volume,
        system_negative_stencil_proxy_of_index,
        system_negative_rhs_of_index
    );

    Init_Cell_Local_Embedding_Domain_System(
        beta_positive_of_index,
        f_positive_of_index,
        cell_multi_index,
        positive_volume_measure,
        multilinear_basis_integrated_over_positive_volume,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_positive_volume,
        system_positive_stencil_proxy_of_index,
        system_positive_rhs_of_index
    );

    const MULTI_INDEX_CUBE<D,0,1> cell_local_multi_index_cube(cell_multi_index);
#ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
    const T boundary_measure = monomials_integrated_over_boundary[1];
    assert(boundary_measure >= 0);
    if(boundary_measure == 0) {
        assert(average_jump_beta_grad_u_dot_n == 0);
        return;
    }
    average_jump_beta_grad_u_dot_n /= boundary_measure;
    for(int i = 1; i <= (1 << D); ++i) {
        const MULTI_INDEX_TYPE multi_index_i = cell_local_multi_index_cube.Multi_Index(i);
        system_negative_rhs_of_index(multi_index_i) -=
            average_jump_beta_grad_u_dot_n * multilinear_basis_integrated_over_boundary[i] / 2;
        system_positive_rhs_of_index(multi_index_i) -=
            average_jump_beta_grad_u_dot_n * multilinear_basis_integrated_over_boundary[i] / 2;
    }
#else // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR
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
        system_negative_rhs_of_index(multi_index_i) -=
            jump_beta_grad_u_dot_n_times_multilinear_basis_integrated_over_boundary[i] / 2;
        system_positive_rhs_of_index(multi_index_i) -=
            jump_beta_grad_u_dot_n_times_multilinear_basis_integrated_over_boundary[i] / 2;
    }
#endif // #ifdef PHYSBAM_MULTIGRID_EMBEDDED_POISSON_USE_B_BAR

    CUBE_STENCIL<T,D,0,1> constraint_negative_stencil;
    CUBE_STENCIL<T,D,0,1> constraint_positive_stencil;
    for(int i = 1; i <= (1 << D); ++i) {
        constraint_negative_stencil.values[i-1] = -multilinear_basis_integrated_over_boundary[i];
        constraint_positive_stencil.values[i-1] = +multilinear_basis_integrated_over_boundary[i];
    }
    constraint_negative_stencil_proxy +=
        CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, constraint_negative_stencil);
    constraint_positive_stencil_proxy +=
        CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, constraint_positive_stencil);
    constraint_rhs = jump_u_integrated_over_boundary;
}

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_BETA_NEGATIVE_OF_INDEX, class T_F_NEGATIVE_OF_INDEX,
    class T_BETA_POSITIVE_OF_INDEX, class T_F_POSITIVE_OF_INDEX,
    class T_JUMP_U_OF_X, class T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
    class T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_NEGATIVE_RHS_OF_INDEX,
    class T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_POSITIVE_RHS_OF_INDEX,
    class T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX,
    class T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX,
    class T_CONSTRAINT_RHS_OF_CELL_INDEX
>
struct INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_VISITOR
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, min_x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, max_x ))
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( typename T_PHI_OF_FINE_INDEX const, phi_of_fine_index ))
        (( typename T_BETA_NEGATIVE_OF_INDEX const, beta_negative_of_index ))
        (( typename T_F_NEGATIVE_OF_INDEX const, f_negative_of_index ))
        (( typename T_BETA_POSITIVE_OF_INDEX const, beta_positive_of_index ))
        (( typename T_F_POSITIVE_OF_INDEX const, f_positive_of_index ))
        (( typename T_JUMP_U_OF_X const, jump_u_of_x ))
        (( typename T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N const, jump_beta_grad_u_dot_n_of_x_and_n ))
        (( /******/ float const, min_dist_to_vertex ))
        (( /******/ int const, sign_of_zero ))
        (( typename T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX const, system_negative_stencil_proxy_of_index ))
        (( typename T_SYSTEM_NEGATIVE_RHS_OF_INDEX const, system_negative_rhs_of_index ))
        (( typename T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX const, system_positive_stencil_proxy_of_index ))
        (( typename T_SYSTEM_POSITIVE_RHS_OF_INDEX const, system_positive_rhs_of_index ))
        (( typename T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX const, constraint_negative_stencil_proxy_of_cell_index ))
        (( typename T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX const, constraint_positive_stencil_proxy_of_cell_index ))
        (( typename T_CONSTRAINT_RHS_OF_CELL_INDEX const, constraint_rhs_of_cell_index ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;
        operator()(cell_multi_index_bound.Multi_Index(cell_linear_index));
    }
    void operator()(const MULTI_INDEX_TYPE cell_multi_index) const
    {
        const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;
        const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;
        const MULTI_INDEX_TYPE fine_cell_multi_index = 2 * cell_multi_index;
        const VECTOR<T,D> x0 = Multi_Index_X(min_x, max_x, fine_multi_index_bound, fine_cell_multi_index);
        const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;
        Init_Cell_Local_Embedding_Interface_System(
            dx,
            phi_of_fine_index,
            beta_negative_of_index, f_negative_of_index,
            beta_positive_of_index, f_positive_of_index,
            jump_u_of_x, jump_beta_grad_u_dot_n_of_x_and_n,
            min_dist_to_vertex, sign_of_zero,
            cell_multi_index, x0,
            system_negative_stencil_proxy_of_index, system_negative_rhs_of_index,
            system_positive_stencil_proxy_of_index, system_positive_rhs_of_index,
            constraint_negative_stencil_proxy_of_cell_index(cell_multi_index),
            constraint_positive_stencil_proxy_of_cell_index(cell_multi_index),
            constraint_rhs_of_cell_index(cell_multi_index)
        );
    }
};

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_BETA_NEGATIVE_OF_INDEX, class T_F_NEGATIVE_OF_INDEX,
    class T_BETA_POSITIVE_OF_INDEX, class T_F_POSITIVE_OF_INDEX,
    class T_JUMP_U_OF_X, class T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
    class T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_NEGATIVE_RHS_OF_INDEX,
    class T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_POSITIVE_RHS_OF_INDEX,
    class T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX,
    class T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX,
    class T_CONSTRAINT_RHS_OF_CELL_INDEX
>
inline INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_VISITOR<
    T, D,
    T_PHI_OF_FINE_INDEX,
    T_BETA_NEGATIVE_OF_INDEX, T_F_NEGATIVE_OF_INDEX,
    T_BETA_POSITIVE_OF_INDEX, T_F_POSITIVE_OF_INDEX,
    T_JUMP_U_OF_X, T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
    T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX,
    T_SYSTEM_NEGATIVE_RHS_OF_INDEX,
    T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX,
    T_SYSTEM_POSITIVE_RHS_OF_INDEX,
    T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX,
    T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX,
    T_CONSTRAINT_RHS_OF_CELL_INDEX
>
Make_Init_Cell_Local_Embedding_Interface_System_Visitor(
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_BETA_NEGATIVE_OF_INDEX& beta_negative_of_index,
    const T_F_NEGATIVE_OF_INDEX& f_negative_of_index,
    const T_BETA_POSITIVE_OF_INDEX& beta_positive_of_index,
    const T_F_POSITIVE_OF_INDEX& f_positive_of_index,
    const T_JUMP_U_OF_X& jump_u_of_x,
    const T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N& jump_beta_grad_u_dot_n_of_x_and_n,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    const T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX& system_negative_stencil_proxy_of_index,
    const T_SYSTEM_NEGATIVE_RHS_OF_INDEX& system_negative_rhs_of_index,
    const T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX& system_positive_stencil_proxy_of_index,
    const T_SYSTEM_POSITIVE_RHS_OF_INDEX& system_positive_rhs_of_index,
    const T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX& constraint_negative_stencil_proxy_of_cell_index,
    const T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX& constraint_positive_stencil_proxy_of_cell_index,
    const T_CONSTRAINT_RHS_OF_CELL_INDEX& constraint_rhs_of_cell_index)
{
    return INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_VISITOR<
        T, D,
        T_PHI_OF_FINE_INDEX,
        T_BETA_NEGATIVE_OF_INDEX, T_F_NEGATIVE_OF_INDEX,
        T_BETA_POSITIVE_OF_INDEX, T_F_POSITIVE_OF_INDEX,
        T_JUMP_U_OF_X, T_JUMP_BETA_GRAD_U_DOT_N_OF_X_AND_N,
        T_SYSTEM_NEGATIVE_STENCIL_PROXY_OF_INDEX,
        T_SYSTEM_NEGATIVE_RHS_OF_INDEX,
        T_SYSTEM_POSITIVE_STENCIL_PROXY_OF_INDEX,
        T_SYSTEM_POSITIVE_RHS_OF_INDEX,
        T_CONSTRAINT_NEGATIVE_STENCIL_PROXY_OF_CELL_INDEX,
        T_CONSTRAINT_POSITIVE_STENCIL_PROXY_OF_CELL_INDEX,
        T_CONSTRAINT_RHS_OF_CELL_INDEX
    >(
        min_x, max_x, multi_index_bound,
        phi_of_fine_index,
        beta_negative_of_index, f_negative_of_index,
        beta_positive_of_index, f_positive_of_index,
        jump_u_of_x, jump_beta_grad_u_dot_n_of_x_and_n,
        min_dist_to_vertex, sign_of_zero,
        system_negative_stencil_proxy_of_index, system_negative_rhs_of_index,
        system_positive_stencil_proxy_of_index, system_positive_rhs_of_index,
        constraint_negative_stencil_proxy_of_cell_index,
        constraint_positive_stencil_proxy_of_cell_index,
        constraint_rhs_of_cell_index
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_INTERFACE_SYSTEM_HPP
