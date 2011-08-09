//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINTS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINTS_HPP

#include <cassert>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T,
    class T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX,
    class T_INDYLESS_AGGREGATE_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX
>
void Aggregate_Constraints(
    T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX constraint_stencil_proxy_of_constraint_index,
    const ARRAY_VIEW<const int> indy_index_of_constraint_index,
    const ARRAY_VIEW<const int> index_of_indy_index,
    ARRAY_VIEW<T> value_of_indy_index,
    T_INDYLESS_AGGREGATE_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_aggregate_constraint_stencil_proxy_of_indy_index)
{
    typedef typename RESULT_OF<
        T_CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_INDEX ( int )
    >::type CONSTRAINT_STENCIL_PROXY_TYPE;
    BOOST_MPL_ASSERT((boost::is_same< typename CONSTRAINT_STENCIL_PROXY_TYPE::INDEX_TYPE, int >));
    BOOST_MPL_ASSERT((boost::is_same< typename CONSTRAINT_STENCIL_PROXY_TYPE::SCALAR_TYPE, T >));

    const int n_constraint = indy_index_of_constraint_index.Size();
    const int n_indy = index_of_indy_index.Size();
    assert(n_indy == index_of_indy_index.Size());
    assert(n_indy == value_of_indy_index.Size());
    static_cast<void>(n_indy);
    assert(value_of_indy_index.Contains_Only(static_cast<T>(0)));

    for(int constraint_index = 1; constraint_index <= n_constraint; ++constraint_index) {
        const int indy_index = indy_index_of_constraint_index(constraint_index);
        const int index = index_of_indy_index(indy_index);
        BOOST_FOREACH(
            typename PHYSBAM_IDENTITY_TYPE(( INDEX_VALUE<int,T> )) const index_value,
            constraint_stencil_proxy_of_constraint_index(constraint_index)
        ) {
            if(index_value.value == 0)
                continue;
            if(index_value.index == index)
                value_of_indy_index(indy_index) += index_value.value;
            else
                indyless_aggregate_constraint_stencil_proxy_of_indy_index(indy_index) += index_value;
        }
    }
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_AGGREGATE_CONSTRAINTS_HPP
