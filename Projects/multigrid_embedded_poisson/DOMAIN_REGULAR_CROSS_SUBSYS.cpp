//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_REGULAR_SUBSYS_BASE.ipp"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
void
DOMAIN_REGULAR_CROSS_SUBSYS<T,D>::
Zero_Stencil(const int linear_index)
{
    STENCIL_TYPE& stencil = stencil_of_index(linear_index);
    const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
    for(MULTI_INDEX_CROSS_ITERATOR<D> it(multi_index); it.Valid(); ++it) {
        const MULTI_INDEX_TYPE other_multi_index = *it;
        const MULTI_INDEX_TYPE multi_offset = other_multi_index - multi_index;
        if(stencil(multi_offset) == 0)
            continue;
        assert(multi_index_bound.Contains(other_multi_index));
        const int other_linear_index = multi_index_bound.Linear_Index(other_multi_index);
        if(other_linear_index == linear_index)
            continue;
        assert(multi_offset != MULTI_INDEX_TYPE());
        STENCIL_TYPE& other_stencil = stencil_of_index(other_linear_index);
        assert(other_stencil(-multi_offset) == stencil(multi_offset));
        other_stencil(-multi_offset) = static_cast<T>(0);
    }
    stencil.Zero();
}

template< class T, int D >
void
DOMAIN_REGULAR_CROSS_SUBSYS<T,D>::
Set_Pure_Neumann_Offset_Grid_BC(
    const MULTI_INDEX_TYPE outside_cell_multi_index,
    const T beta)
{
    assert(!Cell_Multi_Index_Bound().Contains(outside_cell_multi_index));
    const VECTOR<T,D> beta_dv_over_dx_dx = Beta_Dv_Over_Dx_Dx(beta, dx);
    BOOST_FOREACH( const MULTI_INDEX_TYPE cell_multi_offset, (STATIC_MULTI_INDEX_CUBE<D,-1,0>()) ) {
        MULTI_INDEX_TYPE multi_index = outside_cell_multi_index - cell_multi_offset;
        if(!multi_index_bound.Contains(multi_index))
            continue;
        const int linear_index = multi_index_bound.Linear_Index(multi_index);
        STENCIL_TYPE& stencil = stencil_of_index(linear_index);
        for(int d = 1; d <= D; ++d) {
            const int s = 1 + 2 * cell_multi_offset[d];
            multi_index[d] += s;
            if(multi_index_bound.Contains(multi_index)) {
                stencil.Center() += beta_dv_over_dx_dx[d];
                stencil(d,s) -= beta_dv_over_dx_dx[d];
            }
            multi_index[d] -= s;
        }
    }
}

#define EXPLICIT_INSTANTIATION( T, D ) \
template class DOMAIN_REGULAR_CROSS_SUBSYS<T,D>;
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
