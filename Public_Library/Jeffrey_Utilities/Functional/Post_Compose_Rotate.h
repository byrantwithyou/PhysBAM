//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_ROTATE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_ROTATE_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ROTATE_FUNCTION.h>

namespace PhysBAM
{

template< class T_F, class T >
inline COMPOSE_FUNCTION< ROTATE_FUNCTION<T>, T_F >
Post_Compose_Rotate(const T& rot, const T_F& f)
{ return COMPOSE_FUNCTION< ROTATE_FUNCTION<T>, T_F >(ROTATE_FUNCTION<T>(rot), f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_ROTATE_HPP
