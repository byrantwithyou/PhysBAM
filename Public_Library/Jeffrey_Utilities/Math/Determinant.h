//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_DETERMINANT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_DETERMINANT_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
inline T
Determinant(const VECTOR<T,1>& x1)
{ return x1[1]; }

template< class T >
inline T
Determinant(const VECTOR<T,2>& x1, const VECTOR<T,2>& x2)
{ return x1[1] * x2[2] - x1[2] * x2[1]; }

template< class T >
inline T
Determinant(const VECTOR<T,3>& x1, const VECTOR<T,3>& x2, const VECTOR<T,3>& x3)
{
    return x1[1] * (x2[2] * x3[3] - x2[3] * x3[2])
         - x1[2] * (x2[1] * x3[3] - x2[3] * x3[1])
         + x1[3] * (x2[1] * x3[2] - x2[2] * x3[1]);
}

template< class T >
inline T
Determinant(const VECTOR< VECTOR<T,1>, 1 >& x)
{ return Determinant(x[1]); }

template< class T >
inline T
Determinant(const VECTOR< VECTOR<T,2>, 2 >& x)
{ return Determinant(x[1], x[2]); }

template< class T >
inline T
Determinant(const VECTOR< VECTOR<T,3>, 3 >& x)
{ return Determinant(x[1], x[2], x[3]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_DETERMINANT_HPP
