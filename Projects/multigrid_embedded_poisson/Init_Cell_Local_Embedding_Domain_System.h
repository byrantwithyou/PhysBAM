//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_DOMAIN_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_DOMAIN_SYSTEM_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_BETA_OF_INDEX,
    class T_F_OF_INDEX,
    class T_SYSTEM_STENCIL_PROXY_OF_INDEX,
    class T_SYSTEM_RHS_OF_INDEX
>
void Init_Cell_Local_Embedding_Domain_System(
    const T_BETA_OF_INDEX beta_of_index,
    const T_F_OF_INDEX f_of_index,
    const VECTOR<int,D> cell_multi_index,
    const T volume_measure,
    const VECTOR< T, (1 << D) >& multilinear_basis_integrated_over_volume,
    const VECTOR< VECTOR< T, (1 << D) >, (1 << D) >& grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume,
    T_SYSTEM_STENCIL_PROXY_OF_INDEX system_stencil_proxy_of_index,
    T_SYSTEM_RHS_OF_INDEX system_rhs_of_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const MULTI_INDEX_CUBE<D,0,1> cell_local_multi_index_cube(cell_multi_index);

    T average_beta = 0;
    T average_f = 0;
    for(int i = 1; i <= (1 << D); ++i) {
        const MULTI_INDEX_TYPE multi_index = cell_local_multi_index_cube.Multi_Index(i);
        average_beta += beta_of_index(multi_index) * multilinear_basis_integrated_over_volume[i];
        average_f += f_of_index(multi_index) * multilinear_basis_integrated_over_volume[i];
    }
    average_beta /= volume_measure;
    average_f /= volume_measure;

    CUBE_STENCIL<T,D,0,1> cell_local_system_stencil = CUBE_STENCIL<T,D,0,1>::Construct_Zero();
    for(int i = 1; i <= (1 << D); ++i) {
        for(int j = 1; j <= (1 << D); ++j)
            cell_local_system_stencil.values[j-1] =
                average_beta * grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[i][j];
        const MULTI_INDEX_TYPE multi_index_i = cell_local_multi_index_cube.Multi_Index(i);
        system_stencil_proxy_of_index(multi_index_i) +=
            CUBE_STENCIL_PROXY< const CUBE_STENCIL<T,D,0,1> >(cell_multi_index, cell_local_system_stencil);
        system_rhs_of_index(multi_index_i) += average_f * multilinear_basis_integrated_over_volume[i];
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_CELL_LOCAL_EMBEDDING_DOMAIN_SYSTEM_HPP
