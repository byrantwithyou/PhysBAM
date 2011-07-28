//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/Beta_Dv_Over_Dx_Dx.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< int D, class T_STENCIL_OF_INDEX >
struct INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR,
        (( typename MULTI_INDEX_BOUND<D> const, cell_multi_index_bound ))
        (( typename T_STENCIL_OF_INDEX const, stencil_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        BOOST_FOREACH( const MULTI_INDEX_TYPE multi_index, (MULTI_INDEX_CUBE<D,0,1>(cell_multi_index)) ) {
            const int linear_index = multi_index_bound.Linear_Index(multi_index);
            stencil_of_index(linear_index).Add(cell_multi_index - multi_index);
        }
    }
};

template< int D, class T_STENCIL_OF_INDEX >
inline INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR< D, T_STENCIL_OF_INDEX >
Make_Init_Cross_Constbeta_Stencil_Cell_Visitor(
    const MULTI_INDEX_BOUND<D>& cell_multi_index_bound,
    const T_STENCIL_OF_INDEX& stencil_of_index)
{
    return INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR<
        D, T_STENCIL_OF_INDEX
    >(cell_multi_index_bound, stencil_of_index);
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR_HPP
