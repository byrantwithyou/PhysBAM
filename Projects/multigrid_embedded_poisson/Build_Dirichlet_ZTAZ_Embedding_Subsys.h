//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DIRICHLET_ZTAZ_EMBEDDING_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DIRICHLET_ZTAZ_EMBEDDING_SUBSYS_HPP

#include <cassert>

#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/SCALED_STENCIL_PROXY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys
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
    const HASHTABLE<int,int>& indy_index_of_linear_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    T_REGULAR_STENCIL_PROXY regular_stencil_proxy,
    T_EMBEDDING_STENCIL_PROXY embedding_stencil_proxy,
    const T c,
    T_ZTAZ_STENCIL_PROXY ztaz_stencil_proxy,
    const bool add_regular_z_identity);

} // namespace Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys

template<
    class T, int D,
    class T_REGULAR_STENCIL_PROXY_OF_INDEX,
    class T_EMBEDDING_STENCIL_PROXY_OF_INDEX,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX
>
void
Build_Dirichlet_ZTAZ_Embedding_Subsys(
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    T_REGULAR_STENCIL_PROXY_OF_INDEX regular_stencil_proxy_of_index,
    T_EMBEDDING_STENCIL_PROXY_OF_INDEX embedding_stencil_proxy_of_index,
    const HASHTABLE<int,int>& constraint_index_of_cell_linear_index,
    const ARRAY_VIEW<const int> indy_index_of_constraint_index,
    const HASHTABLE<int,int>& indy_index_of_linear_index,
    const ARRAY_VIEW<const int> linear_index_of_indy_index,
    const ARRAY_VIEW<const T> value_of_indy_index,
    T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX indyless_constraint_stencil_proxy_of_indy_index,
    const ARRAY_VIEW<const int> linear_index_of_ztaz_stencil_index,
    T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX ztaz_stencil_proxy_of_stencil_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    const int n_constraint = constraint_index_of_cell_linear_index.Size();
    assert(n_constraint == constraint_index_of_cell_linear_index.Size());
    assert(n_constraint == indy_index_of_constraint_index.Size());
    static_cast<void>(n_constraint);
    const int n_indy = indy_index_of_linear_index.Size();
    assert(n_indy == indy_index_of_linear_index.Size());
    assert(n_indy == linear_index_of_indy_index.Size());
    assert(n_indy == value_of_indy_index.Size());
    static_cast<void>(n_indy);

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = multi_index_bound - 1;
    const int n_ztaz_stencil = linear_index_of_ztaz_stencil_index.Size();
    for(int ztaz_stencil_index = 1; ztaz_stencil_index <= n_ztaz_stencil; ++ztaz_stencil_index) {
        const int linear_index = linear_index_of_ztaz_stencil_index(ztaz_stencil_index);
        if(indy_index_of_linear_index.Contains(linear_index))
            continue;

        typename T_ZTAZ_STENCIL_PROXY_OF_STENCIL_INDEX::result_type ztaz_stencil_proxy =
            ztaz_stencil_proxy_of_stencil_index(ztaz_stencil_index);

        Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys::Multiply_System_Z(
            indy_index_of_linear_index,
            value_of_indy_index,
            indyless_constraint_stencil_proxy_of_indy_index,
            regular_stencil_proxy_of_index(linear_index),
            embedding_stencil_proxy_of_index(linear_index),
            static_cast<T>(1),
            ztaz_stencil_proxy,
            false // add_regular_z_identity
        );

        const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
        BOUNDED_LIST< int, (1 << D) > zt_indy_indices;
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE cell_multi_index,
            Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-1,0>(multi_index), cell_multi_index_bound)
        ) {
            const int cell_linear_index = cell_multi_index_bound.Linear_Index(cell_multi_index);
            const int* const p_constraint_index = constraint_index_of_cell_linear_index.Get_Pointer(cell_linear_index);
            if(!p_constraint_index)
                continue;
            const int constraint_index = *p_constraint_index;
            const int indy_index = indy_index_of_constraint_index(constraint_index);
            zt_indy_indices.Append_Unique(indy_index);
        }

        BOOST_FOREACH( const int indy_index, zt_indy_indices ) {
            const int indy_linear_index = linear_index_of_indy_index(indy_index);
            assert(indy_linear_index != linear_index);
            const T zt_value = -indyless_constraint_stencil_proxy_of_indy_index(indy_index)(linear_index)
                             / value_of_indy_index(indy_index);
            if(zt_value == 0)
                continue;
            Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys::Multiply_System_Z(
                indy_index_of_linear_index,
                value_of_indy_index,
                indyless_constraint_stencil_proxy_of_indy_index,
                regular_stencil_proxy_of_index(indy_linear_index),
                embedding_stencil_proxy_of_index(indy_linear_index),
                zt_value,
                ztaz_stencil_proxy,
                true // add_regular_z_identity
            );
        }

        assert(ztaz_stencil_proxy.N_Nonzero() != 0);
    }
}

namespace Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys
{

template<
    class T,
    class T_INDYLESS_CONSTRAINT_STENCIL_PROXY_OF_INDY_INDEX,
    class T_SUBSYS_STENCIL_PROXY,
    class T_ZTAZ_STENCIL_PROXY
>
inline void
Multiply_Subsys_Z(
    const HASHTABLE<int,int>& indy_index_of_linear_index,
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
        const int* const p_indy_index = indy_index_of_linear_index.Get_Pointer(index_value.index);
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
    const HASHTABLE<int,int>& indy_index_of_linear_index,
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
            indy_index_of_linear_index,
            value_of_indy_index,
            indyless_constraint_stencil_proxy_of_indy_index,
            regular_stencil_proxy,
            c,
            ztaz_stencil_proxy,
            add_regular_z_identity
        );
    Multiply_Subsys_Z(
        indy_index_of_linear_index,
        value_of_indy_index,
        indyless_constraint_stencil_proxy_of_indy_index,
        embedding_stencil_proxy,
        c,
        ztaz_stencil_proxy,
        true
    );
}

} // namespace Detail_Build_Dirichlet_ZTAZ_Embedding_Subsys

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DIRICHLET_ZTAZ_EMBEDDING_SUBSYS_HPP
