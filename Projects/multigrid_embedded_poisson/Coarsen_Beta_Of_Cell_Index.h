//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_BETA_OF_CELL_INDEX_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_BETA_OF_CELL_INDEX_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
void Coarsen_Beta_Of_Cell_Index(
    const MULTI_INDEX_BOUND<D> fine_cell_multi_index_bound,
    const ARRAY_VIEW<const T> beta_of_fine_cell_index,
    ARRAY_VIEW<T> beta_of_coarse_cell_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const MULTI_INDEX_BOUND<D> coarse_cell_multi_index_bound = fine_cell_multi_index_bound / 2;

    assert(beta_of_fine_cell_index.Size() == fine_cell_multi_index_bound.Size());
    assert(beta_of_coarse_cell_index.Size() == coarse_cell_multi_index_bound.Size());

    for(
        int coarse_cell_linear_index = 1;
        coarse_cell_linear_index <= beta_of_coarse_cell_index.Size();
        ++coarse_cell_linear_index
    ) {
        T& coarse_beta = (beta_of_coarse_cell_index(coarse_cell_linear_index) = 0);
        const MULTI_INDEX_TYPE coarse_cell_multi_index =
            coarse_cell_multi_index_bound.Multi_Index(coarse_cell_linear_index);
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE fine_cell_multi_index,
            (MULTI_INDEX_CUBE<D,0,1>(2 * coarse_cell_multi_index - 1))
        ) {
            const int fine_cell_linear_index = fine_cell_multi_index_bound.Linear_Index(fine_cell_multi_index);
            coarse_beta += beta_of_fine_cell_index(fine_cell_linear_index);
        }
        coarse_beta /= 1 << D;
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_COARSEN_BETA_OF_CELL_INDEX_HPP
