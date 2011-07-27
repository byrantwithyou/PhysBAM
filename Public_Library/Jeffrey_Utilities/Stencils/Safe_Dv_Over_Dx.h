//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SAFE_DV_OVER_DX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SAFE_DV_OVER_DX_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
inline VECTOR<T,1>
Safe_Dv_Over_Dx(const VECTOR<T,1>& /*dx*/)
{ return 1; }

template< class T >
inline VECTOR<T,2>
Safe_Dv_Over_Dx(const VECTOR<T,2>& dx)
{ return VECTOR<T,2>(dx[2], dx[1]); }

template< class T >
inline VECTOR<T,3>
Safe_Dv_Over_Dx(const VECTOR<T,3>& dx)
{ return VECTOR<T,3>(dx[2]*dx[3], dx[1]*dx[3], dx[1]*dx[2]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SAFE_DV_OVER_DX_HPP
