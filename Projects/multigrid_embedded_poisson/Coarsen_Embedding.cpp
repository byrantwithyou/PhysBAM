//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

#include "Coarsen_Embedding.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< int D >
void Coarsen_Embedding(
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound,
    const HASHTABLE<int,int>& fine_embedding_index_of_linear_index,
    HASHTABLE<int,int>& coarse_embedding_index_of_linear_index,
    ARRAY<int>& coarse_linear_index_of_embedding_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef MULTI_INDEX_BOX_ITERATOR<D> MULTI_INDEX_BOX_ITERATOR_TYPE;

    assert(0 == coarse_embedding_index_of_linear_index.Size());
    assert(0 == coarse_linear_index_of_embedding_index.Size());

    const MULTI_INDEX_BOUND<D> coarse_multi_index_bound = (fine_multi_index_bound + 1) / 2;
    const int n_coarse = coarse_multi_index_bound.Size();

    int n_coarse_embedding = 0;
    BOOST_FOREACH( const MULTI_INDEX_TYPE coarse_multi_index, coarse_multi_index_bound ) {
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE fine_multi_index,
            Multi_Index_Box_Intersect(
                MULTI_INDEX_CUBE<D,-1,+1>(2 * coarse_multi_index - 1),
                fine_multi_index_bound
            )
        ) {
            const int fine_linear_index = fine_multi_index_bound.Linear_Index(fine_multi_index);
            if(fine_embedding_index_of_linear_index.Contains(fine_linear_index)) {
                ++n_coarse_embedding;
                break;
            }
        }
    }

    coarse_embedding_index_of_linear_index.Initialize_New_Table(n_coarse_embedding);
    coarse_linear_index_of_embedding_index.Preallocate(n_coarse_embedding);

    for(int coarse_linear_index = 1; coarse_linear_index <= n_coarse; ++coarse_linear_index) {
        const MULTI_INDEX_TYPE coarse_multi_index = coarse_multi_index_bound.Multi_Index(coarse_linear_index);
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE fine_multi_index,
            Multi_Index_Box_Intersect(
                MULTI_INDEX_CUBE<D,-1,+1>(2 * coarse_multi_index - 1),
                fine_multi_index_bound
            )
        ) {
            const int fine_linear_index = fine_multi_index_bound.Linear_Index(fine_multi_index);
            if(fine_embedding_index_of_linear_index.Contains(fine_linear_index)) {
                assert(
                    coarse_linear_index_of_embedding_index.Size() == 0
                 || coarse_linear_index_of_embedding_index.Last() < coarse_linear_index
                );
                coarse_linear_index_of_embedding_index.Append(coarse_linear_index);
                coarse_embedding_index_of_linear_index.Insert(
                    coarse_linear_index,
                    coarse_linear_index_of_embedding_index.Size()
                );
                break;
            }
        }
    }
}

#define EXPLICIT_INSTANTIATION( D ) \
template void Coarsen_Embedding<D>( \
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound, \
    const HASHTABLE<int,int>& fine_embedding_index_of_linear_index, \
    HASHTABLE<int,int>& coarse_embedding_index_of_linear_index, \
    ARRAY<int>& coarse_linear_index_of_embedding_index);
EXPLICIT_INSTANTIATION( 2 )
EXPLICIT_INSTANTIATION( 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
