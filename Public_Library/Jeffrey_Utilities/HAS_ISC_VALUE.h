//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_ISC_VALUE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_ISC_VALUE_HPP

#include <boost/type_traits/integral_constant.hpp>
#include <boost/type_traits/is_class.hpp>

namespace PhysBAM
{

namespace Detail_HAS_ISC_VALUE
{

struct YES_TYPE { char _[1]; };
struct  NO_TYPE { char _[2]; };

template< int > struct SFINAE { typedef YES_TYPE type; };
template< class U > typename SFINAE< U::value >::type Has_Isc(int);
template< class U > NO_TYPE Has_Isc(...);

template< class T, bool = boost::is_class<T>::value >
struct HAS_ISC_VALUE_DISPATCH
{
    static const bool value = sizeof( Has_Isc<T>(0) ) == sizeof( YES_TYPE );
    typedef HAS_ISC_VALUE_DISPATCH type;
};

template< class T >
struct HAS_ISC_VALUE_DISPATCH< T, false >
    : boost::false_type
{ };

} // namespace Detail_HAS_ISC_VALUE

template< class T >
struct HAS_ISC_VALUE
    : Detail_HAS_ISC_VALUE::HAS_ISC_VALUE_DISPATCH<T>
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_ISC_VALUE_HPP
