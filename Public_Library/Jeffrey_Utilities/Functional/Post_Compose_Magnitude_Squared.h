//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MAGNITUDE_SQUARED_OF_RESULT_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MAGNITUDE_SQUARED_OF_RESULT_FUNCTION_HPP

#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/MAGNITUDE_SQUARED_FUNCTION.h>

namespace PhysBAM
{

template< class T_F >
inline COMPOSE_FUNCTION< MAGNITUDE_SQUARED_FUNCTION, T_F >
Post_Compose_Magnitude_Squared(const T_F& f)
{ return COMPOSE_FUNCTION< MAGNITUDE_SQUARED_FUNCTION, T_F >(MAGNITUDE_SQUARED_FUNCTION(), f); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MAGNITUDE_SQUARED_OF_RESULT_FUNCTION_HPP
