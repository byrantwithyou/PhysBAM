//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_HPP

#include <cassert>

#include <algorithm>

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/LOOP.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_FWD.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM
{

template< class T_INDEX, class T, int MAX_N_NONZERO /*= -1*/ >
struct UNSTRUCTURED_STENCIL
{
    typedef T_INDEX INDEX_TYPE;
    typedef T SCALAR_TYPE;

    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    typedef typename boost::mpl::if_c<
        (MAX_N_NONZERO < 0),
        ARRAY< INDEX_VALUE_TYPE >,
        BOUNDED_LIST< INDEX_VALUE_TYPE, MAX_N_NONZERO >
    >::type VALUES_TYPE;
    VALUES_TYPE values;

    T Sum() const;
    bool Contains(const INDEX_TYPE& index) const;

    UNSTRUCTURED_STENCIL& operator+=(const INDEX_VALUE_TYPE& index_value);

    template< class T_STENCIL_PROXY >
    UNSTRUCTURED_STENCIL& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    UNSTRUCTURED_STENCIL& operator+=(const T_STENCIL_PROXY stencil_proxy);

    /***/ T& operator()(const INDEX_TYPE& index);
    const T& operator()(const INDEX_TYPE& index) const;
    /***/ T& operator()(const INDEX_TYPE& index, T& default_);
    const T& operator()(const INDEX_TYPE& index, const T& default_) const;

    template< class T_VECTOR >
    T Apply(const T_VECTOR& x) const;
    template< class T_VECTOR >
    void Apply_Transpose(T_VECTOR& y, const T x) const;

    typedef typename VALUES_TYPE::/***/ iterator /***/ iterator;
    typedef typename VALUES_TYPE::const_iterator const_iterator;
    /***/ iterator begin();
    const_iterator begin() const;
    /***/ iterator end();
    const_iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline T
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
Sum() const
{
    T sum = 0;
    for(int i = 1; i <= values.Size(); ++i)
        sum += values(i).value;
    return sum;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline bool
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
Contains(const INDEX_TYPE& index) const
{
    const INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index);
    return p != values.end() && p->index == index;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    if(index_value.value != 0) {
        INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index_value.index);
        if(p != values.end() && p->index == index_value.index)
            p->value += index_value.value;
        else
            values.Insert(index_value, 1 + static_cast< int >(p - values.begin()));
    }
    return *this;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
template< class T_STENCIL_PROXY >
inline UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    values.Remove_All();
    //values.Preallocate(stencil_proxy.N_Nonzero());
    return operator+=(stencil_proxy);
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
template< class T_STENCIL_PROXY >
inline UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator+=(const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_FOREACH( const INDEX_VALUE_TYPE index_value, stencil_proxy )
        operator+=(index_value);
    return *this;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline T&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator()(const INDEX_TYPE& index)
{
    INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index);
    assert(p != values.end() && p->index == index);
    return p->value;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline const T&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator()(const INDEX_TYPE& index) const
{
    const INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index);
    assert(p != values.end() && p->index == index);
    return p->value;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline T&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator()(const INDEX_TYPE& index, T& default_)
{
    INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index);
    return p == values.end() || p->index != index ? default_ : p->value;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline const T&
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
operator()(const INDEX_TYPE& index, const T& default_) const
{
    const INDEX_VALUE_TYPE* const p = std::lower_bound(values.begin(), values.end(), index);
    return p == values.end() || p->index != index ? default_ : p->value;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
template< class T_VECTOR >
inline T
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
Apply(const T_VECTOR& x) const
{
    T y = 0;
    for(int i = 1; i <= values.Size(); ++i)
        y += values(i).Apply(x);
    return y;
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
template< class T_VECTOR >
inline void
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
Apply_Transpose(T_VECTOR& y, const T x) const
{
    for(int i = 1; i <= values.Size(); ++i)
        values(i).Apply_Transpose(y, x);
}

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline typename UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::iterator
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
begin()
{ return values.begin(); }

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline typename UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::const_iterator
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
begin() const
{ return values.begin(); }

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline typename UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::iterator
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
end()
{ return values.end(); }

template< class T_INDEX, class T, int MAX_N_NONZERO >
inline typename UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::const_iterator
UNSTRUCTURED_STENCIL< T_INDEX, T, MAX_N_NONZERO >::
end() const
{ return values.end(); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_UNSTRUCTURED_STENCIL_HPP
