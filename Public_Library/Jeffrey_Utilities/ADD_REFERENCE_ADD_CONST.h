//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ADD_REFERENCE_ADD_CONST_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ADD_REFERENCE_ADD_CONST_HPP

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_reference.hpp>

namespace PhysBAM
{

template< class T >
struct ADD_REFERENCE_ADD_CONST
    : boost::add_reference<
          typename boost::add_const<T>::type
      >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ADD_REFERENCE_ADD_CONST_HPP
