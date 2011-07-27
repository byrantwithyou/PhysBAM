//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SCALE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SCALE_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SCALE_FUNCTION.h>

namespace PhysBAM
{

template< class T_F, class T >
inline COMPOSE_FUNCTION< SCALE_FUNCTION<T>, T_F >
Post_Compose_Scale(const T& c, const T_F& f)
{ return COMPOSE_FUNCTION< SCALE_FUNCTION<T>, T_F >(SCALE_FUNCTION<T>(c), f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SCALE_HPP
