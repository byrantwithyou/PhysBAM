//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_DIRICHLET_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_DIRICHLET_SYSTEM_HPP

#include <Jeffrey_Utilities/Functional/Pre_Compose_Translate.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_FUNCTION_OVER_BOUNDARY_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Math/GAUSSIAN_QUADRATURE_SIMPLEX_RULE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Evaluate_Cell_Local_Integrations.h"
#include "Init_Cell_Local_Embedding_Domain_System.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_BETA_OF_INDEX,
    class T_F_OF_INDEX,
    class T_U_OF_X,
    class T_SYSTEM_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_RHS_OF_INDEX,
    class T_CONSTRAINT_STENCIL_PROXY
>
void Init_Cell_Local_Embedding_Dirichlet_System(
    const VECTOR<T,D> min_x, const VECTOR<T,D> max_x,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    const int domain_sign,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_BETA_OF_INDEX& beta_of_index,
    const T_F_OF_INDEX& f_of_index,
    const T_U_OF_X& u_of_x,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    const VECTOR<int,D> cell_multi_index,
    const T_SYSTEM_STENCIL_PROXY_OF_INDEX& system_stencil_proxy_of_index,
    const T_SYSTEM_RHS_OF_INDEX& system_rhs_of_index,
    T_CONSTRAINT_STENCIL_PROXY constraint_stencil_proxy,
    T& constraint_rhs)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const VECTOR<T,D> dx = (max_x - min_x) / (multi_index_bound.max_multi_index - 1);
    const VECTOR<T,D> x0 = Multi_Index_X(min_x, max_x, 2 * multi_index_bound - 1, 2 * cell_multi_index);

    VECTOR< T, STATIC_POW_C<3,D>::value > monomials_integrated_over_volume; // init'ed to 0
    VECTOR< T, STATIC_POW_C<2,D>::value > monomials_integrated_over_boundary; // init'ed to 0
    VECTOR< T, (1 << D) > multilinear_basis_integrated_over_volume; // init'ed to 0
    VECTOR< T, (1 << D) > multilinear_basis_integrated_over_boundary; // init'ed to 0
    VECTOR< VECTOR< T, (1 << D) >, (1 << D) > grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume; // init'ed to 0

    T u_integrated_over_boundary = 0;

    Evaluate_Cell_Local_Integrations(
        dx,
        domain_sign,
        phi_of_fine_index,
        min_dist_to_vertex,
        sign_of_zero,
        Make_Integrate_Function_Over_Boundary_Polytope_Visitor<
            GAUSSIAN_QUADRATURE_SIMPLEX_RULE<T,D,2>
        >(domain_sign, Pre_Compose_Translate(u_of_x, x0), u_integrated_over_boundary),
        cell_multi_index,
        monomials_integrated_over_volume,
        monomials_integrated_over_boundary,
        multilinear_basis_integrated_over_volume,
        multilinear_basis_integrated_over_boundary,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume
    );

    Init_Cell_Local_Embedding_Domain_System(
        multi_index_bound,
        beta_of_index, f_of_index,
        cell_multi_index,
        monomials_integrated_over_volume[1],
        multilinear_basis_integrated_over_volume,
        grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume,
        system_stencil_proxy_of_index, system_rhs_of_index
    );

    CUBE_STENCIL<T,D,0,1> constraint_stencil;
    for(int i = 1; i <= (1 << D); ++i)
        constraint_stencil.values[i-1] = multilinear_basis_integrated_over_boundary[i];
    constraint_stencil_proxy = CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, constraint_stencil);
    constraint_rhs = u_integrated_over_boundary;
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CELL_LOCAL_EMBEDDING_DIRICHLET_SYSTEM_HPP
