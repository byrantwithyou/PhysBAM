//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_ZTAZ_EMBEDDING_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_ZTAZ_EMBEDDING_HPP

#include <cassert>

#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_REGULAR_STENCIL_PROXY_OF_INDY_INDEX >
void
Init_ZTAZ_Embedding(
    const HASHTABLE<int,int>& embedding_index_of_index,
    const ARRAY_VIEW<const int> index_of_embedding_index,
    const HASHTABLE<int,int>& indy_index_of_index,
    const ARRAY_VIEW<const int> index_of_indy_index,
    T_REGULAR_STENCIL_PROXY_OF_INDY_INDEX regular_stencil_proxy_of_indy_index,
    HASHTABLE<int,int>& ztaz_embedding_index_of_index,
    ARRAY<int>& index_of_ztaz_embedding_index)
{
    typedef typename RESULT_OF< T_REGULAR_STENCIL_PROXY_OF_INDY_INDEX ( int ) >::type REGULAR_STENCIL_PROXY_TYPE;
    BOOST_MPL_ASSERT((boost::is_same< typename REGULAR_STENCIL_PROXY_TYPE::INDEX_TYPE, int >));
    typedef typename REGULAR_STENCIL_PROXY_TYPE::INDEX_VALUE_TYPE INDEX_VALUE_TYPE;

    const int n_embedding = index_of_embedding_index.Size();
    const int n_indy = index_of_indy_index.Size();

    int n_ztaz_embedding = n_embedding - n_indy;
    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
        const int index = index_of_indy_index(indy_index);
        assert(embedding_index_of_index.Contains(index));
        BOOST_FOREACH(
            INDEX_VALUE_TYPE const index_value,
            regular_stencil_proxy_of_indy_index(indy_index)
        ) {
            if(
                index_value.value != 0
             && index_value.index != index
             && !embedding_index_of_index.Contains(index_value.index)
            )
                ++n_ztaz_embedding;
        }
    }

    ztaz_embedding_index_of_index.Initialize_New_Table(n_ztaz_embedding);
    index_of_ztaz_embedding_index.Preallocate(n_ztaz_embedding);
    for(int embedding_index = 1; embedding_index <= n_embedding; ++embedding_index) {
        const int index = index_of_embedding_index(embedding_index);
        if(indy_index_of_index.Contains(index))
            continue;
        ztaz_embedding_index_of_index.Insert(index, 0);
        index_of_ztaz_embedding_index.Append(index);
    }
    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
        const int index = index_of_indy_index(indy_index);
        BOOST_FOREACH(
            INDEX_VALUE_TYPE const index_value,
            regular_stencil_proxy_of_indy_index(indy_index)
        ) {
            if(index_value.value == 0 || index_value.index == index)
                continue;
            if(ztaz_embedding_index_of_index.Set(index_value.index, 0)) {
                assert(!embedding_index_of_index.Contains(index_value.index));
                index_of_ztaz_embedding_index.Append(index_value.index);
            }
        }
    }

    assert(ztaz_embedding_index_of_index.Size() <= n_ztaz_embedding);
    assert(index_of_ztaz_embedding_index.Size() <= n_ztaz_embedding);
    n_ztaz_embedding = index_of_ztaz_embedding_index.Size();
    assert(n_ztaz_embedding == ztaz_embedding_index_of_index.Size());
    assert(n_ztaz_embedding == index_of_ztaz_embedding_index.Size());
    std::sort(
        index_of_ztaz_embedding_index.begin(),
        index_of_ztaz_embedding_index.end()
    );
    for(int ztaz_embedding_index = 1; ztaz_embedding_index <= n_ztaz_embedding; ++ztaz_embedding_index) {
        const int index = index_of_ztaz_embedding_index(ztaz_embedding_index);
        ztaz_embedding_index_of_index.Get(index) = ztaz_embedding_index;
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INIT_ZTAZ_EMBEDDING_HPP
