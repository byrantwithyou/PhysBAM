//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_MULTIPLY_STENCILS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_MULTIPLY_STENCILS_HPP

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/Stencils/SCALED_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template<
    class T_A_STENCIL_PROXY,
    class T_B_STENCIL_PROXY_OF_INDEX,
    class T_RESULT_STENCIL_PROXY
>
inline void Multiply_Stencils(
    const T_A_STENCIL_PROXY a_stencil_proxy,
    const T_B_STENCIL_PROXY_OF_INDEX b_stencil_proxy_of_index,
    T_RESULT_STENCIL_PROXY result_stencil_proxy)
{
    typedef typename T_B_STENCIL_PROXY_OF_INDEX::result_type T_B_STENCIL_PROXY;
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_A_STENCIL_PROXY >));
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_B_STENCIL_PROXY >));
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_RESULT_STENCIL_PROXY >));

    typedef typename T_RESULT_STENCIL_PROXY::SCALAR_TYPE SCALAR_TYPE;
    BOOST_MPL_ASSERT((boost::is_same< typename T_A_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_B_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));

    BOOST_FOREACH( typename T_A_STENCIL_PROXY::reference a, a_stencil_proxy )
        if(a.value != 0)
            result_stencil_proxy += Make_Scaled_Stencil_Proxy(b_stencil_proxy_of_index(a.index), a.value);
}

template<
    class T_A_STENCIL_PROXY,
    class T_B_STENCIL_PROXY_OF_INDEX,
    class T_C_STENCIL_PROXY_OF_INDEX,
    class T_RESULT_STENCIL_PROXY
>
inline void Multiply_Stencils(
    const T_A_STENCIL_PROXY a_stencil_proxy,
    const T_B_STENCIL_PROXY_OF_INDEX b_stencil_proxy_of_index,
    const T_C_STENCIL_PROXY_OF_INDEX c_stencil_proxy_of_index,
    T_RESULT_STENCIL_PROXY result_stencil_proxy)
{
    typedef typename T_B_STENCIL_PROXY_OF_INDEX::result_type T_B_STENCIL_PROXY;
    typedef typename T_C_STENCIL_PROXY_OF_INDEX::result_type T_C_STENCIL_PROXY;
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_A_STENCIL_PROXY >));
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_B_STENCIL_PROXY >));
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_C_STENCIL_PROXY >));
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_RESULT_STENCIL_PROXY >));

    typedef typename T_RESULT_STENCIL_PROXY::SCALAR_TYPE SCALAR_TYPE;
    BOOST_MPL_ASSERT((boost::is_same< typename T_A_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_B_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_C_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));

    BOOST_FOREACH( typename T_A_STENCIL_PROXY::reference a, a_stencil_proxy )
        if(a.value != 0)
            Multiply_Stencils(
                Make_Scaled_Stencil_Proxy(b_stencil_proxy_of_index(a.index), a.value),
                c_stencil_proxy_of_index,
                result_stencil_proxy
            );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_MULTIPLY_STENCILS_HPP
