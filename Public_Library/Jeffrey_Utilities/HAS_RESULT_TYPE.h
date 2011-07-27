//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_RESULT_TYPE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_RESULT_TYPE_HPP

#include <boost/type_traits/is_class.hpp>

namespace PhysBAM
{

namespace Detail_HAS_RESULT_TYPE
{

struct YES_TYPE { char _[1]; };
struct  NO_TYPE { char _[2]; };

template< class > struct SFINAE { typedef YES_TYPE type; };
template< class U > typename SFINAE< typename U::RESULT_TYPE >::type Has_Type(int);
template< class U > NO_TYPE Has_Type(...);

template< class T, bool = boost::is_class<T>::value >
struct HAS_RESULT_TYPE_DISPATCH
{
    static const bool value = sizeof( Has_Type<T>(0) ) == sizeof( YES_TYPE );
    typedef HAS_RESULT_TYPE_DISPATCH type;
};

template< class T >
struct HAS_RESULT_TYPE_DISPATCH< T, false >
    : boost::false_type
{ };

} // namespace Detail_HAS_RESULT_TYPE

template< class T >
struct HAS_RESULT_TYPE
    : Detail_HAS_RESULT_TYPE::HAS_RESULT_TYPE_DISPATCH<T>
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_HAS_RESULT_TYPE_HPP
