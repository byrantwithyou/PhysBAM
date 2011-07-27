//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< class T_BASE >
struct SKIP_ZERO_VALUE_STENCIL_PROXY
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_BASE >));

    typedef typename T_BASE::INDEX_TYPE INDEX_TYPE;
    typedef typename T_BASE::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    typedef T_BASE BASE_TYPE;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        SKIP_ZERO_VALUE_STENCIL_PROXY, (( typename BASE_TYPE, base ))
    )

    template< class T_BASE2 >
    SKIP_ZERO_VALUE_STENCIL_PROXY(const SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE2 >& other,
        typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type* = 0);

    SKIP_ZERO_VALUE_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    template< class T_STENCIL_PROXY >
    SKIP_ZERO_VALUE_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    SKIP_ZERO_VALUE_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    typedef SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< typename BASE_TYPE::iterator > iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

template< class T_BASE >
inline SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >
Make_Skip_Zero_Value_Stencil_Proxy(const T_BASE& base)
{ return SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >(base); }

struct SKIP_ZERO_VALUE_STENCIL_PROXY_FUNCTION
{
    template<class> struct result;
    template< class T_THIS, class T_BASE >
    struct result< T_THIS ( T_BASE ) >
    {
        typedef SKIP_ZERO_VALUE_STENCIL_PROXY<
            typename REMOVE_QUALIFIERS< T_BASE >::type
        > type;
    };

    template< class T_BASE >
    SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >
    operator()(const T_BASE& base)
    { return SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >(base); }
};

//#####################################################################
//#####################################################################

template< class T_BASE >
template< class T_BASE2 >
inline
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
SKIP_ZERO_VALUE_STENCIL_PROXY(const SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE2 >& other,
    typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type*)
    : base(other.base)
{ }

template< class T_BASE >
inline SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >&
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    base += index_value;
    return *this;
}

template< class T_BASE >
template< class T_STENCIL_PROXY >
inline SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >&
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base = stencil_proxy;
    return *this;
}

template< class T_BASE >
template< class T_STENCIL_PROXY >
inline SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >&
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base += stencil_proxy;
    return *this;
}

template< class T_BASE >
inline int
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
N_Nonzero() const
{ return base.N_Nonzero(); }

template< class T_BASE >
inline typename SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::iterator
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
begin() const
{ return iterator(base.begin(), base.end()); }

template< class T_BASE >
inline typename SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::iterator
SKIP_ZERO_VALUE_STENCIL_PROXY< T_BASE >::
end() const
{ return iterator(base.end(), base.end()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_HPP
