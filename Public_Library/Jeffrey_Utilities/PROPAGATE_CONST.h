//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//##################################################################### 
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_CONST_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_CONST_HPP

#include <boost/type_traits/add_const.hpp>

namespace PhysBAM
{

template< class T, class U >
struct PROPAGATE_CONST
{ typedef U type; };

template< class T, class U >
struct PROPAGATE_CONST< const T, U >
    : boost::add_const<U>
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PROPAGATE_CONST_HPP
