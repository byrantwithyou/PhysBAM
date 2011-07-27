//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/range/iterator.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/utility/enable_if.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/PROPAGATE_CONST.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< class T_STENCIL >
struct UNSTRUCTURED_STENCIL_PROXY
{
    typedef typename T_STENCIL::INDEX_TYPE INDEX_TYPE;
    typedef typename T_STENCIL::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;
    typedef typename boost::remove_const< T_STENCIL >::type STENCIL_TYPE;

    typedef typename PROPAGATE_CONST< T_STENCIL, SCALAR_TYPE >::type & SCALAR_REFERENCE;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        UNSTRUCTURED_STENCIL_PROXY, (( typename T_STENCIL&, stencil ))
    )

    template< class T_STENCIL2 >
    UNSTRUCTURED_STENCIL_PROXY(const UNSTRUCTURED_STENCIL_PROXY< T_STENCIL2 >& other,
        typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type* = 0);

    UNSTRUCTURED_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    UNSTRUCTURED_STENCIL_PROXY& operator=(const UNSTRUCTURED_STENCIL_PROXY& other);
    template< class T_STENCIL_PROXY >
    UNSTRUCTURED_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    UNSTRUCTURED_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    SCALAR_REFERENCE operator()(const INDEX_TYPE& index) const;

    typedef typename boost::range_iterator< T_STENCIL >::type iterator;
    typedef iterator const_iterator;
    typedef typename boost::iterator_reference< iterator >::type reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_STENCIL >
template< class T_STENCIL2 >
inline
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
UNSTRUCTURED_STENCIL_PROXY(const UNSTRUCTURED_STENCIL_PROXY< T_STENCIL2 >& other,
    typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type*)
    : stencil(other.stencil)
{ }

template< class T_STENCIL >
inline UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >&
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    stencil += index_value;
    return *this;
}

template< class T_STENCIL >
inline UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >&
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
operator=(const UNSTRUCTURED_STENCIL_PROXY& other)
{ return operator=< UNSTRUCTURED_STENCIL_PROXY >(other); }

template< class T_STENCIL >
template< class T_STENCIL_PROXY >
inline UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >&
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil = stencil_proxy;
    return *this;
}

template< class T_STENCIL >
template< class T_STENCIL_PROXY >
inline UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >&
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil += stencil_proxy;
    return *this;
}

template< class T_STENCIL >
inline int
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
N_Nonzero() const
{ return stencil.values.Size(); }

template< class T_STENCIL >
inline typename UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::SCALAR_REFERENCE
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
operator()(const INDEX_TYPE& index) const
{ return stencil(index); }

template< class T_STENCIL >
inline typename UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::iterator
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
begin() const
{ return stencil.begin(); }

template< class T_STENCIL >
inline typename UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::iterator
UNSTRUCTURED_STENCIL_PROXY< T_STENCIL >::
end() const
{ return stencil.end(); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_PROXY_HPP
