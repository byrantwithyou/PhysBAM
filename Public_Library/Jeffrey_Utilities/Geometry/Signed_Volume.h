//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Signed_Volume( [ I_{nxn} 0_{nx1} ] ) = 1
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SIGNED_VOLUME_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SIGNED_VOLUME_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Math/Determinant.h>

namespace PhysBAM
{

template< class T >
inline T
Signed_Volume(const VECTOR<T,1>& x1, const VECTOR<T,1>& x2)
{ return Determinant(x1 - x2); }

template< class T >
inline T
Signed_Volume(const VECTOR<T,2>& x1, const VECTOR<T,2>& x2, const VECTOR<T,2>& x3)
{ return Determinant(x1 - x3, x2 - x3); }

template< class T >
inline T
Signed_Volume(const VECTOR<T,3>& x1, const VECTOR<T,3>& x2, const VECTOR<T,3>& x3, const VECTOR<T,3>& x4)
{ return Determinant(x1 - x4, x2 - x4, x3 - x4); }

template< class T >
inline T
Signed_Volume(const VECTOR< VECTOR<T,1>, 2 >& x)
{ return Signed_Volume(x[1], x[2]); }

template< class T >
inline T
Signed_Volume(const VECTOR< VECTOR<T,2>, 3 >& x)
{ return Signed_Volume(x[1], x[2], x[3]); }

template< class T >
inline T
Signed_Volume(const VECTOR< VECTOR<T,3>, 4 >& x)
{ return Signed_Volume(x[1], x[2], x[3], x[4]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_SIGNED_VOLUME_HPP
