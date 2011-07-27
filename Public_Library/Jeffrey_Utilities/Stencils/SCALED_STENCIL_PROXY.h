//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/SCALED_STENCIL_PROXY_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< class T_BASE >
struct SCALED_STENCIL_PROXY
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_BASE >));

    typedef typename T_BASE::INDEX_TYPE INDEX_TYPE;
    typedef typename T_BASE::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    typedef T_BASE BASE_TYPE;

    BASE_TYPE base;
    const SCALAR_TYPE scaling;

    explicit SCALED_STENCIL_PROXY(const BASE_TYPE& base_);
    SCALED_STENCIL_PROXY(const BASE_TYPE& base_, const SCALAR_TYPE scaling_);
    template< class T_BASE2 >
    SCALED_STENCIL_PROXY(const SCALED_STENCIL_PROXY< T_BASE2 >& other,
        typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type* = 0);

    SCALED_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    template< class T_BASE2 >
    SCALED_STENCIL_PROXY& operator=(const SCALED_STENCIL_PROXY< T_BASE2 >& other);
    template< class T_BASE2 >
    SCALED_STENCIL_PROXY& operator+=(const SCALED_STENCIL_PROXY< T_BASE2 >& other);

    template< class T_STENCIL_PROXY >
    SCALED_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    SCALED_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    typedef SCALED_STENCIL_PROXY_ITERATOR<
        typename BASE_TYPE::iterator
    > iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_BASE >
inline SCALED_STENCIL_PROXY< T_BASE >
operator*(const SCALED_STENCIL_PROXY< T_BASE >& stencil_proxy, typename T_BASE::SCALAR_TYPE const scaling)
{ return SCALED_STENCIL_PROXY< T_BASE >(stencil_proxy.base, stencil_proxy.scaling * scaling); }

template< class T_BASE >
inline SCALED_STENCIL_PROXY< T_BASE >
operator*(typename T_BASE::SCALAR_TYPE const scaling, const SCALED_STENCIL_PROXY< T_BASE >& stencil_proxy)
{ return SCALED_STENCIL_PROXY< T_BASE >(stencil_proxy.base, scaling * stencil_proxy.scaling); }

template< class T_BASE >
inline SCALED_STENCIL_PROXY< T_BASE >
operator/(const SCALED_STENCIL_PROXY< T_BASE >& stencil_proxy, typename T_BASE::SCALAR_TYPE const scaling)
{ return SCALED_STENCIL_PROXY< T_BASE >(stencil_proxy.base, stencil_proxy.scaling / scaling); }

//#####################################################################
//#####################################################################

template< class T_STENCIL_PROXY >
struct MAKE_SCALED_STENCIL_PROXY_RESULT
{ typedef SCALED_STENCIL_PROXY< T_STENCIL_PROXY > type; };

template< class T_BASE >
struct MAKE_SCALED_STENCIL_PROXY_RESULT< SCALED_STENCIL_PROXY< T_BASE > >
{ typedef SCALED_STENCIL_PROXY< T_BASE > type; };

template< class T_STENCIL_PROXY >
struct MAKE_SCALED_STENCIL_PROXY_RESULT< T_STENCIL_PROXY& >
    : MAKE_SCALED_STENCIL_PROXY_RESULT< T_STENCIL_PROXY >
{ };

template< class T_STENCIL_PROXY >
struct MAKE_SCALED_STENCIL_PROXY_RESULT< const T_STENCIL_PROXY >
    : MAKE_SCALED_STENCIL_PROXY_RESULT< T_STENCIL_PROXY >
{ };

template< class T_STENCIL_PROXY >
inline SCALED_STENCIL_PROXY< T_STENCIL_PROXY >
Make_Scaled_Stencil_Proxy(
    const T_STENCIL_PROXY& stencil_proxy,
    typename T_STENCIL_PROXY::SCALAR_TYPE const scaling)
{ return SCALED_STENCIL_PROXY< T_STENCIL_PROXY >(stencil_proxy, scaling); }

template< class T_BASE >
inline SCALED_STENCIL_PROXY< T_BASE >
Make_Scaled_Stencil_Proxy(
    const SCALED_STENCIL_PROXY< T_BASE >& stencil_proxy,
    typename T_BASE::SCALAR_TYPE const scaling)
{ return stencil_proxy * scaling; }

//#####################################################################
//#####################################################################

template< class T_BASE >
inline
SCALED_STENCIL_PROXY< T_BASE >::
SCALED_STENCIL_PROXY(const BASE_TYPE& base_)
    : base(base_),
      scaling(1)
{ }

template< class T_BASE >
inline
SCALED_STENCIL_PROXY< T_BASE >::
SCALED_STENCIL_PROXY(const BASE_TYPE& base_, const SCALAR_TYPE scaling_)
    : base(base_),
      scaling(scaling_)
{ }

template< class T_BASE >
template< class T_BASE2 >
inline
SCALED_STENCIL_PROXY< T_BASE >::
SCALED_STENCIL_PROXY(const SCALED_STENCIL_PROXY< T_BASE2 >& other,
    typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type*)
    : base(other.base),
      scaling(other.scaling)
{ }

template< class T_BASE >
inline SCALED_STENCIL_PROXY< T_BASE >&
SCALED_STENCIL_PROXY< T_BASE >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    base += index_value / scaling;
    return *this;
}

template< class T_BASE >
template< class T_BASE2 >
inline SCALED_STENCIL_PROXY< T_BASE >&
SCALED_STENCIL_PROXY< T_BASE >::
operator=(const SCALED_STENCIL_PROXY< T_BASE2 >& other)
{
    BOOST_MPL_ASSERT((boost::is_same< typename T_BASE2::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_BASE2::SCALAR_TYPE, SCALAR_TYPE >));
    base = SCALED_STENCIL_PROXY< T_BASE2 >(other.base, other.scaling / scaling);
    return *this;
}

template< class T_BASE >
template< class T_BASE2 >
inline SCALED_STENCIL_PROXY< T_BASE >&
SCALED_STENCIL_PROXY< T_BASE >::
operator+=(const SCALED_STENCIL_PROXY< T_BASE2 >& other)
{
    BOOST_MPL_ASSERT((boost::is_same< typename T_BASE2::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_BASE2::SCALAR_TYPE, SCALAR_TYPE >));
    base += SCALED_STENCIL_PROXY< T_BASE2 >(other.base, other.scaling / scaling);
    return *this;
}

template< class T_BASE >
template< class T_STENCIL_PROXY >
inline SCALED_STENCIL_PROXY< T_BASE >&
SCALED_STENCIL_PROXY< T_BASE >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base = SCALED_STENCIL_PROXY< T_STENCIL_PROXY >(stencil_proxy, 1/scaling);
}

template< class T_BASE >
template< class T_STENCIL_PROXY >
inline SCALED_STENCIL_PROXY< T_BASE >&
SCALED_STENCIL_PROXY< T_BASE >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base += SCALED_STENCIL_PROXY< T_STENCIL_PROXY >(stencil_proxy, 1/scaling);
}

template< class T_BASE >
inline int
SCALED_STENCIL_PROXY< T_BASE >::
N_Nonzero() const
{ return base.N_Nonzero(); }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY< T_BASE >::iterator
SCALED_STENCIL_PROXY< T_BASE >::
begin() const
{ return iterator(base.begin(), scaling); }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY< T_BASE >::iterator
SCALED_STENCIL_PROXY< T_BASE >::
end() const
{ return iterator(base.end(), scaling); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_HPP
