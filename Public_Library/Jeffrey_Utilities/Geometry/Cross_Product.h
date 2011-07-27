//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Cross_Product( [ I_{nxn} ] ) = e_{n+1}
//              ( [ 0_{1xn} ] )
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CROSS_PRODUCT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CROSS_PRODUCT_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
inline VECTOR<T,1>
Cross_Product()
{ return VECTOR<T,1>(1); }

template< class T >
inline VECTOR<T,2>
Cross_Product(const VECTOR<T,2>& x1)
{ return VECTOR<T,2>(-x1[2], x1[1]); }

template< class T >
inline VECTOR<T,3>
Cross_Product(const VECTOR<T,3>& x1, const VECTOR<T,3>& x2)
{ return VECTOR<T,3>(x1[2]*x2[3] - x1[3]*x2[2], x1[3]*x2[1] - x1[1]*x2[3], x1[1]*x2[2] - x1[2]*x2[1]); }

template< class T >
inline VECTOR<T,1>
Cross_Product(const VECTOR< VECTOR<T,1>, 0 >& x)
{ return Cross_Product<T>(); }

template< class T >
inline VECTOR<T,2>
Cross_Product(const VECTOR< VECTOR<T,2>, 1 >& x)
{ return Cross_Product(x[1]); }

template< class T >
inline VECTOR<T,3>
Cross_Product(const VECTOR< VECTOR<T,3>, 2 >& x)
{ return Cross_Product(x[1], x[2]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CROSS_PRODUCT_HPP
