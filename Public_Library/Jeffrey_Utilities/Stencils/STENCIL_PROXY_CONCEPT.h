//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_CONCEPT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_CONCEPT_HPP

#include <boost/concept/usage.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_STENCIL_PROXY >
struct CONST_STENCIL_PROXY_CONCEPT
{
    typedef typename T_STENCIL_PROXY::INDEX_TYPE INDEX_TYPE;
    typedef typename T_STENCIL_PROXY::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    typedef typename T_STENCIL_PROXY::iterator iterator;
    typedef typename T_STENCIL_PROXY::reference reference;
    BOOST_MPL_ASSERT((boost::is_convertible< reference, INDEX_VALUE_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same<
        typename boost::remove_const< typename boost::remove_reference<
            typename boost::remove_reference<
                typename T_STENCIL_PROXY::reference
            >::type::INDEX_TYPE
        >::type >::type,
        INDEX_TYPE
    >));
    BOOST_MPL_ASSERT((boost::is_same<
        typename boost::remove_const< typename boost::remove_reference<
            typename boost::remove_reference<
                typename T_STENCIL_PROXY::reference
            >::type::VALUE_TYPE
        >::type >::type,
        SCALAR_TYPE
    >));

    BOOST_CONCEPT_USAGE( CONST_STENCIL_PROXY_CONCEPT )
    {
        int n_nonzero = stencil_proxy.N_Nonzero();
        iterator begin_it = stencil_proxy.begin();
        iterator end_it = stencil_proxy.end();

        static_cast< void >(n_nonzero);
        static_cast< void >(begin_it);
        static_cast< void >(end_it);
    }

private:
    T_STENCIL_PROXY stencil_proxy;
};

template< class T_STENCIL_PROXY >
struct STENCIL_PROXY_CONCEPT
    : CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >
{
    typedef typename T_STENCIL_PROXY::INDEX_TYPE INDEX_TYPE;
    typedef typename T_STENCIL_PROXY::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    BOOST_CONCEPT_USAGE( STENCIL_PROXY_CONCEPT )
    {
        stencil_proxy = other_stencil_proxy;
        stencil_proxy += index_value;
        stencil_proxy += other_stencil_proxy;
    }

private:
    T_STENCIL_PROXY stencil_proxy;
    T_STENCIL_PROXY other_stencil_proxy;
    INDEX_VALUE_TYPE index_value;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_CONCEPT_HPP
