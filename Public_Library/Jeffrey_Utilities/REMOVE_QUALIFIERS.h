//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_QUALIFIERS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_QUALIFIERS_HPP

#include <boost/type_traits/remove_cv.hpp>
#include <boost/type_traits/remove_reference.hpp>

namespace PhysBAM
{

template< class T >
struct REMOVE_QUALIFIERS
    : boost::remove_cv<
          typename boost::remove_reference<T>::type
      >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_QUALIFIERS_HPP
