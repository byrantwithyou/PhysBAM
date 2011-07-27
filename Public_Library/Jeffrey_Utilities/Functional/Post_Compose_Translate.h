//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_TRANSLATE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_TRANSLATE_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/TRANSLATE_FUNCTION.h>

namespace PhysBAM
{

template< class T_F, class T >
inline COMPOSE_FUNCTION< TRANSLATE_FUNCTION<T>, T_F >
Post_Compose_Translate(const T& x0, const T_F& f)
{ return COMPOSE_FUNCTION< TRANSLATE_FUNCTION<T>, T_F >(TRANSLATE_FUNCTION<T>(x0), f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_POST_COMPOSE_TRANSLATE_HPP
