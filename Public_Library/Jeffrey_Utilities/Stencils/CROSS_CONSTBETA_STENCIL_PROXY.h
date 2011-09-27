//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/utility/enable_if.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CROSS_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_ITERATOR.h>

namespace PhysBAM
{

template< class T_STENCIL, class T >
struct CROSS_CONSTBETA_STENCIL_PROXY
{
    typedef typename T_STENCIL::MULTI_INDEX_TYPE INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;
    typedef typename boost::remove_const< T_STENCIL >::type STENCIL_TYPE;
    typedef T_STENCIL& STENCIL_REFERENCE;

    static const int DIMENSION = T_STENCIL::DIMENSION;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        CROSS_CONSTBETA_STENCIL_PROXY,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR< T, DIMENSION > )) const, beta_dv_over_dx_dx ))
        (( typename INDEX_TYPE const, base_multi_index ))
        (( typename STENCIL_REFERENCE, stencil ))
    )

    template< class T_STENCIL2 >
    CROSS_CONSTBETA_STENCIL_PROXY(const CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL2, T >& other,
        typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type* = 0);

    CROSS_CONSTBETA_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    CROSS_CONSTBETA_STENCIL_PROXY& operator=(const CROSS_CONSTBETA_STENCIL_PROXY& other);
    template< class T_STENCIL_PROXY >
    CROSS_CONSTBETA_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    CROSS_CONSTBETA_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    T Sum() const;

    typedef STENCIL_PROXY_ITERATOR<
        MULTI_INDEX_CROSS_ITERATOR< DIMENSION >,
        CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR< T, DIMENSION >
    > iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_STENCIL, class T >
template< class T_STENCIL2 >
inline
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
CROSS_CONSTBETA_STENCIL_PROXY(const CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL2, T >& other,
    typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type*)
    : beta_dv_over_dx_dx(other.beta_dv_over_dx_dx),
      base_multi_index(other.base_multi_index),
      stencil(other.stencil)
{ }

template< class T_STENCIL, class T >
inline CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >&
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    stencil.Add(beta_dv_over_dx_dx, index_value.index - base_multi_index, index_value.value);
    return *this;
}

template< class T_STENCIL, class T >
inline CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >&
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
operator=(const CROSS_CONSTBETA_STENCIL_PROXY& other)
{ return operator=< CROSS_CONSTBETA_STENCIL_PROXY >(other); }

template< class T_STENCIL, class T >
template< class T_STENCIL_PROXY >
inline CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >&
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil.Assign(beta_dv_over_dx_dx, base_multi_index, stencil_proxy);
    return *this;
}

template< class T_STENCIL, class T >
template< class T_STENCIL_PROXY >
inline CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >&
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil.Add(beta_dv_over_dx_dx, base_multi_index, stencil_proxy);
    return *this;
}

template< class T_STENCIL, class T >
inline int
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
N_Nonzero() const
{ return stencil.N_Nonzero(); }

template< class T_STENCIL, class T >
inline T
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
Sum() const
{ return stencil.Sum(beta_dv_over_dx_dx); }

template< class T_STENCIL, class T >
inline typename CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::iterator
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
begin() const
{
    return iterator(
        MULTI_INDEX_CROSS_ITERATOR< DIMENSION >(base_multi_index, BEGIN_TAG()),
        CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR< T, DIMENSION >(beta_dv_over_dx_dx, &stencil.values[0], BEGIN_TAG())
    );
}

template< class T_STENCIL, class T >
inline typename CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::iterator
CROSS_CONSTBETA_STENCIL_PROXY< T_STENCIL, T >::
end() const
{
    return iterator(
        MULTI_INDEX_CROSS_ITERATOR< DIMENSION >(base_multi_index, END_TAG()),
        CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR< T, DIMENSION >(beta_dv_over_dx_dx, &stencil.values[T_STENCIL::N_VALUES], END_TAG())
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_PROXY_HPP
