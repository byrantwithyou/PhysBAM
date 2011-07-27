//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_TRANSLATE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_TRANSLATE_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/TRANSLATE_FUNCTION.h>

namespace PhysBAM
{

template< class T_F, class T >
inline COMPOSE_FUNCTION< T_F, TRANSLATE_FUNCTION<T> >
Pre_Compose_Translate(const T_F& f, const T& x0)
{ return COMPOSE_FUNCTION< T_F, TRANSLATE_FUNCTION<T> >(f, TRANSLATE_FUNCTION<T>(x0)); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_TRANSLATE_HPP
