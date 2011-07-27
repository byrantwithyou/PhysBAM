//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_STENCIL_CELL_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_STENCIL_CELL_VISITOR_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/Beta_Dv_Over_Dx_Dx.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D, class T_BETA_OF_CELL_INDEX, class T_STENCIL_OF_INDEX >
struct INIT_CROSS_STENCIL_CELL_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INIT_CROSS_STENCIL_CELL_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx ))
        (( typename MULTI_INDEX_BOUND<D> const, cell_multi_index_bound ))
        (( typename T_BETA_OF_CELL_INDEX const, beta_of_cell_index ))
        (( typename T_STENCIL_OF_INDEX const, stencil_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        const T beta = beta_of_cell_index(cell_linear_index);
        const VECTOR<T,D> beta_dv_over_dx_dx = Beta_Dv_Over_Dx_Dx(beta, dx);
        BOOST_FOREACH( const MULTI_INDEX_TYPE multi_index, (MULTI_INDEX_CUBE<D,0,1>(cell_multi_index)) ) {
            const int linear_index = multi_index_bound.Linear_Index(multi_index);
            stencil_of_index(linear_index).Add(beta_dv_over_dx_dx, cell_multi_index - multi_index);
        }
    }
};

template< class T, int D, class T_BETA_OF_CELL_INDEX, class T_STENCIL_OF_INDEX >
inline INIT_CROSS_STENCIL_CELL_VISITOR< T, D, T_BETA_OF_CELL_INDEX, T_STENCIL_OF_INDEX >
Make_Init_Cross_Stencil_Cell_Visitor(
    const VECTOR<T,D>& dx,
    const MULTI_INDEX_BOUND<D>& cell_multi_index_bound,
    const T_BETA_OF_CELL_INDEX& beta_of_cell_index,
    const T_STENCIL_OF_INDEX& stencil_of_index)
{
    return INIT_CROSS_STENCIL_CELL_VISITOR<
        T, D, T_BETA_OF_CELL_INDEX, T_STENCIL_OF_INDEX
    >(dx, cell_multi_index_bound, beta_of_cell_index, stencil_of_index);
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_STENCIL_CELL_VISITOR_HPP
