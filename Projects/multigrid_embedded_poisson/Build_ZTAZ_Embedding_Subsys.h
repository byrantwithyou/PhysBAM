//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_ZTAZ_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_ZTAZ_EMBEDDING_SUBSYS_HPP

#include <cassert>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/SCALED_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Detail_Build_ZTAZ_Embedding_Subsys
{

template<
    class T,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_REGULAR_STENCIL_PROXY,
    class T_EMBEDDING_STENCIL_PROXY,
    class T_ZTAZ_STENCIL_PROXY
>
inline void
Multiply_System_Z(
    const HASHTABLE<int,int>& indy_index_of_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    T_REGULAR_STENCIL_PROXY regular_stencil_proxy,
    T_EMBEDDING_STENCIL_PROXY embedding_stencil_proxy,
    const T c,
    T_ZTAZ_STENCIL_PROXY ztaz_stencil_proxy,
    const bool add_regular_z_identity);

} // namespace Detail_Build_ZTAZ_Embedding_Subsys

template<
    class T, int N,
    class T_REGULAR_STENCIL_PROXY_OF_INDEX,
    class T_EMBEDDING_STENCIL_PROXY_OF_INDEX,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX
>
void
Build_ZTAZ_Embedding_Subsys(
    T_REGULAR_STENCIL_PROXY_OF_INDEX regular_stencil_proxy_of_index,
    T_EMBEDDING_STENCIL_PROXY_OF_INDEX embedding_stencil_proxy_of_index,
    const HASHTABLE<int,int>& indy_index_of_index,
    const ARRAY_VIEW<const int> index_of_indy_index,
    const HASHTABLE< int, BOUNDED_LIST<int,N> >& zt_indy_indices_of_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    const ARRAY_VIEW<const int> index_of_ztaz_stencil_index,
    T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX ztaz_stencil_proxy_of_stencil_index)
{
    typedef typename RESULT_OF< T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX ( int ) >::type ZTAZ_STENCIL_PROXY_TYPE;

    const int n_indy = indy_index_of_index.Size();
    assert(n_indy == indy_index_of_index.Size());
    assert(n_indy == index_of_indy_index.Size());
    assert(n_indy == value_of_indy_index.Size());
    static_cast<void>(n_indy);

    const int n_ztaz_stencil = index_of_ztaz_stencil_index.Size();
    for(int ztaz_stencil_index = 1; ztaz_stencil_index <= n_ztaz_stencil; ++ztaz_stencil_index) {
        const int index = index_of_ztaz_stencil_index(ztaz_stencil_index);
        if(indy_index_of_index.Contains(index))
            continue;

        ZTAZ_STENCIL_PROXY_TYPE ztaz_stencil_proxy = ztaz_stencil_proxy_of_stencil_index(ztaz_stencil_index);

        Detail_Build_ZTAZ_Embedding_Subsys::Multiply_System_Z(
            indy_index_of_index,
            value_of_indy_index,
            indyless_constraint_stencil_proxy_of_indy_index,
            regular_stencil_proxy_of_index(index),
            embedding_stencil_proxy_of_index(index),
            static_cast<T>(1),
            ztaz_stencil_proxy,
            false // add_regular_z_identity
        );

        const BOUNDED_LIST<int,N>* const p_zt_indy_indices = zt_indy_indices_of_index.Get_Pointer(index);
        if(!p_zt_indy_indices)
            continue;
        BOOST_FOREACH( const int zt_indy_index, *p_zt_indy_indices ) {
            const int zt_index = index_of_indy_index(zt_indy_index);
            assert(zt_index != index);
            const T zt_value = -indyless_constraint_stencil_proxy_of_indy_index(zt_indy_index)(index)
                             / value_of_indy_index(zt_indy_index);
            if(zt_value == 0)
                continue;
            Detail_Build_ZTAZ_Embedding_Subsys::Multiply_System_Z(
                indy_index_of_index,
                value_of_indy_index,
                indyless_constraint_stencil_proxy_of_indy_index,
                regular_stencil_proxy_of_index(zt_index),
                embedding_stencil_proxy_of_index(zt_index),
                zt_value,
                ztaz_stencil_proxy,
                true // add_regular_z_identity
            );
        }

        assert(ztaz_stencil_proxy.N_Nonzero() != 0);
    }
}

namespace Detail_Build_ZTAZ_Embedding_Subsys
{

template<
    class T,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_SUBSYS_STENCIL_PROXY,
    class T_ZTAZ_STENCIL_PROXY
>
inline void
Multiply_Subsys_Z(
    const HASHTABLE<int,int>& indy_index_of_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    T_SUBSYS_STENCIL_PROXY subsys_stencil_proxy,
    const T c,
    T_ZTAZ_STENCIL_PROXY ztaz_stencil_proxy,
    const bool add_subsys_z_identity)
{
    typedef INDEX_VALUE<int,T> INDEX_VALUE_TYPE;
    BOOST_FOREACH( const INDEX_VALUE_TYPE index_value, subsys_stencil_proxy ) {
        if(index_value.value == 0)
            continue;
        const int* const p_indy_index = indy_index_of_index.Get_Pointer(index_value.index);
        if(p_indy_index) {
            const int indy_index = *p_indy_index;
            ztaz_stencil_proxy += Make_Scaled_Stencil_Proxy(
                indyless_constraint_stencil_proxy_of_indy_index(indy_index),
                -c * index_value.value / value_of_indy_index(indy_index)
            );
        }
        else if(add_subsys_z_identity) {
            ztaz_stencil_proxy += c * index_value;
        }
    }
}

template<
    class T,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_REGULAR_STENCIL_PROXY,
    class T_EMBEDDING_STENCIL_PROXY,
    class T_ZTAZ_STENCIL_PROXY
>
inline void
Multiply_System_Z(
    const HASHTABLE<int,int>& indy_index_of_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    T_REGULAR_STENCIL_PROXY regular_stencil_proxy,
    T_EMBEDDING_STENCIL_PROXY embedding_stencil_proxy,
    const T c,
    T_ZTAZ_STENCIL_PROXY ztaz_stencil_proxy,
    const bool add_regular_z_identity)
{
    BOOST_MPL_ASSERT((boost::is_same< typename T_REGULAR_STENCIL_PROXY::INDEX_TYPE, int >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_REGULAR_STENCIL_PROXY::SCALAR_TYPE, T >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_EMBEDDING_STENCIL_PROXY::INDEX_TYPE, int >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_EMBEDDING_STENCIL_PROXY::SCALAR_TYPE, T >));
    assert(add_regular_z_identity || c == 1);
    if(regular_stencil_proxy.N_Nonzero() != 0)
        Multiply_Subsys_Z(
            indy_index_of_index,
            value_of_indy_index,
            indyless_constraint_stencil_proxy_of_indy_index,
            regular_stencil_proxy,
            c,
            ztaz_stencil_proxy,
            add_regular_z_identity
        );
    Multiply_Subsys_Z(
        indy_index_of_index,
        value_of_indy_index,
        indyless_constraint_stencil_proxy_of_indy_index,
        embedding_stencil_proxy,
        c,
        ztaz_stencil_proxy,
        true
    );
}

} // namespace Detail_Build_ZTAZ_Embedding_Subsys

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_ZTAZ_EMBEDDING_SUBSYS_HPP
