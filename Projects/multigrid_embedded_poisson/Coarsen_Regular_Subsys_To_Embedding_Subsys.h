//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_REGULAR_SUBSYS_TO_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_REGULAR_SUBSYS_TO_EMBEDDING_SUBSYS_HPP

#include <cassert>

#include <algorithm>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/Beta_Dv_Over_Dx_Dx.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/CUBE_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/GEOMETRIC_PROLONGATION_STENCIL_PROXY_FUNCTION.h>
#include <Jeffrey_Utilities/Stencils/GEOMETRIC_RESTRICTION_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/Multiply_Stencils.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

//#include "LOCAL_CROSS_SUBSYS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_SIGN_OF_FINE_CELL_INDEX,
    class T_BETA_OF_FINE_CELL_INDEX,
    class T_SIGN_OF_COARSE_CELL_INDEX,
    class T_STENCIL_PROXY_OF_COARSE_INDEX
>
void Coarsen_Regular_Subsys_To_Embedding_Subsys(
    const int sign,
    const VECTOR<T,D> fine_dx,
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound,
    T_SIGN_OF_FINE_CELL_INDEX sign_of_fine_cell_index,
    T_BETA_OF_FINE_CELL_INDEX beta_of_fine_cell_index,
    T_SIGN_OF_COARSE_CELL_INDEX sign_of_coarse_cell_index,
    T_STENCIL_PROXY_OF_COARSE_INDEX stencil_proxy_of_coarse_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    assert(sign == -1 || sign == 1);

    const MULTI_INDEX_BOUND<D> fine_cell_multi_index_bound = fine_multi_index_bound - 1;
    const MULTI_INDEX_BOUND<D> coarse_multi_index_bound = (fine_multi_index_bound + 1) / 2;
    const MULTI_INDEX_BOUND<D> coarse_cell_multi_index_bound = coarse_multi_index_bound - 1;

    const int n_fine_cell = fine_cell_multi_index_bound.Size();
    const int n_coarse_cell = coarse_cell_multi_index_bound.Size();
    for(int coarse_cell_linear_index = 1; coarse_cell_linear_index <= n_coarse_cell; ++coarse_cell_linear_index) {
        if(sign_of_coarse_cell_index(coarse_cell_linear_index) != 0)
            continue;
        const MULTI_INDEX_TYPE coarse_cell_multi_index = coarse_cell_multi_index_bound.Multi_Index(coarse_cell_linear_index);
        const MULTI_INDEX_TYPE fine_base_multi_index = 2 * coarse_cell_multi_index;

#if 0
        LOCAL_CROSS_SUBSYS<T,D,-1,+1> fine_local_subsys = LOCAL_CROSS_SUBSYS<T,D,-1,+1>::Construct_Zero();
        bool fine_local_subsys_is_zero = true;
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE fine_cell_multi_index,
            (MULTI_INDEX_CUBE<D,-1,0>(fine_base_multi_index))
        ) {
            const int fine_cell_linear_index = fine_cell_multi_index_bound.Linear_Index(fine_cell_multi_index);
            assert(1 <= fine_cell_linear_index && fine_cell_linear_index <= n_fine_cell);
            if(sign_of_fine_cell_index(fine_cell_linear_index) != sign)
                continue;
            fine_local_subsys_is_zero = false;
            const T beta = beta_of_fine_cell_index(fine_cell_linear_index);
            fine_local_subsys.Add(Beta_Dv_Over_Dx_Dx(beta, fine_dx), fine_cell_multi_index - fine_base_multi_index);
        }
        if(fine_local_subsys_is_zero)
            continue;

        BOOST_FOREACH(
            const MULTI_INDEX_TYPE coarse_multi_index,
            (MULTI_INDEX_CUBE<D,0,+1>(coarse_cell_multi_index))
        ) {
            // TODO: Need to revisit this...CUBE_STENCIL's might not be enough...
            typedef CUBE_STENCIL<T,D,-1,+1> COARSE_LOCAL_STENCIL_TYPE;
            typedef CUBE_STENCIL_PROXY< COARSE_LOCAL_STENCIL_TYPE > COARSE_LOCAL_STENCIL_PROXY_TYPE;
            COARSE_LOCAL_STENCIL_TYPE coarse_local_stencil = COARSE_LOCAL_STENCIL_TYPE::Construct_Zero();
            COARSE_LOCAL_STENCIL_PROXY_TYPE coarse_local_stencil_proxy(coarse_multi_index, coarse_local_stencil);
            Multiply_Stencils(
                GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>(coarse_multi_index),
                fine_local_subsys.Multi_Index_Proxy_Function(fine_multi_index_bound, fine_base_multi_index),
                GEOMETRIC_PROLONGATION_STENCIL_PROXY_FUNCTION<T,D>(fine_multi_index_bound),
                coarse_local_stencil_proxy
            );
            stencil_proxy_of_coarse_index(coarse_multi_index) += coarse_local_stencil_proxy;
        }
#endif // #if 0

    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_REGULAR_SUBSYS_TO_EMBEDDING_SUBSYS_HPP
