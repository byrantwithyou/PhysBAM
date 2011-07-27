//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//##################################################################### 
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_QUALIFIERS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_QUALIFIERS_HPP

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_cv.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/type_traits/add_volatile.hpp>

namespace PhysBAM
{

template< class T, class U >
struct PROPAGATE_QUALIFIERS
{ typedef U type; };

template< class T, class U >
struct PROPAGATE_QUALIFIERS< const T, U >
    : boost::add_const<U>
{ };

template< class T, class U >
struct PROPAGATE_QUALIFIERS< volatile T, U >
    : boost::add_volatile<U>
{ };

template< class T, class U >
struct PROPAGATE_QUALIFIERS< const volatile T, U >
    : boost::add_cv<U>
{ };

template< class T, class U >
struct PROPAGATE_QUALIFIERS< T&, U >
    : boost::add_reference<
          typename PROPAGATE_QUALIFIERS<T,U>::type
      >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_QUALIFIERS_HPP
