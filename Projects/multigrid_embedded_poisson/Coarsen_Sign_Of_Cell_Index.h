//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_COARSEN_SIGN_OF_CELL_INDEX_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_COARSEN_SIGN_OF_CELL_INDEX_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< int D, class T_SIGN >
void Coarsen_Sign_Of_Cell(
    const MULTI_INDEX_BOUND<D> fine_cell_multi_index_bound,
    const ARRAY_VIEW< const T_SIGN > sign_of_fine_cell_index,
    ARRAY_VIEW< T_SIGN > sign_of_coarse_cell_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const MULTI_INDEX_BOUND<D> coarse_cell_multi_index_bound = fine_cell_multi_index_bound / 2;

    assert(sign_of_fine_cell_index.Size() == fine_cell_multi_index_bound.Size());
    assert(sign_of_coarse_cell_index.Size() == coarse_cell_multi_index_bound.Size());

    for(
        int coarse_cell_linear_index = 1;
        coarse_cell_linear_index <= sign_of_coarse_cell_index.Size();
        ++coarse_cell_linear_index
    ) {
        T_SIGN& coarse_sign = (sign_of_coarse_cell_index(coarse_cell_linear_index) = 0);
        const MULTI_INDEX_TYPE coarse_cell_multi_index =
            coarse_cell_multi_index_bound.Multi_Index(coarse_cell_linear_index);
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE fine_cell_multi_index,
            (MULTI_INDEX_CUBE<D,0,1>(2 * coarse_cell_multi_index - 1))
        ) {
            const int fine_cell_linear_index = fine_cell_multi_index_bound.Linear_Index(fine_cell_multi_index);
            const int fine_sign = sign_of_fine_cell_index(fine_cell_linear_index);
            assert(coarse_sign == 0 || fine_sign == 0 || coarse_sign == fine_sign);
            if((coarse_sign = fine_sign) == 0)
                break;
        }
    }
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_COARSEN_SIGN_OF_CELL_INDEX_HPP
