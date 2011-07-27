//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_NEGATE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_NEGATE_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/NEGATE_FUNCTION.h>

namespace PhysBAM
{

template< class T_F >
inline COMPOSE_FUNCTION< T_F, NEGATE_FUNCTION >
Pre_Compose_Negate(const T_F& f)
{ return COMPOSE_FUNCTION< T_F, NEGATE_FUNCTION >(f, NEGATE_FUNCTION()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_PRE_COMPOSE_NEGATE_HPP
