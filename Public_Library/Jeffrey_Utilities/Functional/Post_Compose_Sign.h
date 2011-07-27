//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SIGN_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>

namespace PhysBAM
{

template< class T_F >
inline COMPOSE_FUNCTION< SIGN_FUNCTION, T_F >
Post_Compose_Sign(const T_F& f)
{ return COMPOSE_FUNCTION< SIGN_FUNCTION, T_F >(SIGN_FUNCTION(), f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_SIGN_HPP
