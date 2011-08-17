//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_UNWRAP_KRYLOV_VECTOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_UNWRAP_KRYLOV_VECTOR_HPP

#include <boost/type_traits/add_reference.hpp>

#include <Jeffrey_Utilities/ADD_REFERENCE_ADD_CONST.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>

namespace PhysBAM
{

template< class T_VECTOR >
inline T_VECTOR&
Unwrap_Krylov_Vector(T_VECTOR& x)
{ return x; }

template< class T, class T_VECTOR >
inline typename boost::add_reference< T_VECTOR >::type
Unwrap_Krylov_Vector(KRYLOV_VECTOR_WRAPPER< T, T_VECTOR >& x)
{ return x.v; }

template< class T, class T_VECTOR >
inline typename ADD_REFERENCE_ADD_CONST< T_VECTOR >::type
Unwrap_Krylov_Vector(const KRYLOV_VECTOR_WRAPPER< T, T_VECTOR >& x)
{ return x.v; }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_UNWRAP_KRYLOV_VECTOR_HPP
