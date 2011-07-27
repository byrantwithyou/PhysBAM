//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_VALUE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_VALUE_HPP

#include <boost/call_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/ADD_REFERENCE_ADD_CONST.h>

namespace PhysBAM
{

template< class T_INDEX, class T_VALUE >
struct INDEX_VALUE
{
    typedef T_INDEX INDEX_TYPE;
    typedef T_VALUE VALUE_TYPE;

    T_INDEX index;
    T_VALUE value;

private:
    typedef typename boost::call_traits< INDEX_TYPE >::param_type INDEX_PARAM_TYPE;
    typedef typename boost::call_traits< VALUE_TYPE >::param_type VALUE_PARAM_TYPE;
public:

    INDEX_VALUE()
    { }

    explicit INDEX_VALUE(INDEX_PARAM_TYPE index_)
        : index(index_)
    { }

    INDEX_VALUE(INDEX_PARAM_TYPE index_, VALUE_PARAM_TYPE value_)
        : index(index_),
          value(value_)
    { }

    template< class T_INDEX2, class T_VALUE2 >
    INDEX_VALUE(const INDEX_VALUE< T_INDEX2, T_VALUE2 >& other,
        typename boost::enable_if< boost::mpl::and_<
            boost::is_convertible<
                typename ADD_REFERENCE_ADD_CONST< T_INDEX2 >::type,
                INDEX_TYPE
            >,
            boost::is_convertible<
                typename ADD_REFERENCE_ADD_CONST< T_VALUE2 >::type,
                VALUE_TYPE
            >
        > >::type* = 0)
        : index(other.index),
          value(other.value)
    { }

    template< class T_VECTOR >
    VALUE_TYPE Apply(const T_VECTOR& x) const
    { return value * x(index); }
    template< class T_VECTOR >
    void Apply_Transpose(T_VECTOR& y, VALUE_PARAM_TYPE x) const
    { y(index) += value * x; }
};

template< class T_INDEX, class T_VALUE >
inline bool operator<(
    const INDEX_VALUE< T_INDEX, T_VALUE >& index_value1,
    const INDEX_VALUE< T_INDEX, T_VALUE >& index_value2)
{ return index_value1.index < index_value2.index; }

template< class T_INDEX, class T_VALUE >
inline bool operator<(
    const INDEX_VALUE< T_INDEX, T_VALUE >& index_value,
    typename boost::call_traits< T_INDEX >::param_type index)
{ return index_value.index < index; }

template< class T_INDEX, class T_VALUE >
inline bool operator<(
    typename boost::call_traits< T_INDEX >::param_type index,
    const INDEX_VALUE< T_INDEX, T_VALUE >& index_value)
{ return index < index_value.index; }

template< class T_INDEX, class T_VALUE >
inline INDEX_VALUE< T_INDEX, T_VALUE >
operator*(const INDEX_VALUE< T_INDEX, T_VALUE >& index_value, const T_VALUE value)
{ return INDEX_VALUE< T_INDEX, T_VALUE >(index_value.index, index_value.value * value); }

template< class T_INDEX, class T_VALUE >
inline INDEX_VALUE< T_INDEX, T_VALUE >
operator*(const T_VALUE value, const INDEX_VALUE< T_INDEX, T_VALUE >& index_value)
{ return INDEX_VALUE< T_INDEX, T_VALUE >(index_value.index, value * index_value.value); }

template< class T_INDEX, class T_VALUE >
inline INDEX_VALUE< T_INDEX, T_VALUE >
operator/(const INDEX_VALUE< T_INDEX, T_VALUE >& index_value, const T_VALUE value)
{ return INDEX_VALUE< T_INDEX, T_VALUE >(index_value.index, index_value.value / value); }

//#####################################################################
//#####################################################################

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_VALUE_HPP
