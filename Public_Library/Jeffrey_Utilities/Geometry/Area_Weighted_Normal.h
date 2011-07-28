//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Area_Weighted_Normal( [ I_{nxn} 0_{nx1} ] ) = e_{n+1}
//                     ( [ 0_{1xn} 0_{1x1} ] )
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_AREA_WEIGHTED_NORMAL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_AREA_WEIGHTED_NORMAL_HPP

#include <Jeffrey_Utilities/Geometry/Cross_Product.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
inline VECTOR<T,1>
Area_Weighted_Normal(const VECTOR<T,1>& /*x1*/)
{ return Cross_Product<T>(); }

template< class T >
inline VECTOR<T,2>
Area_Weighted_Normal(const VECTOR<T,2>& x1, const VECTOR<T,2>& x2)
{ return Cross_Product(x1 - x2); }

template< class T >
inline VECTOR<T,3>
Area_Weighted_Normal(const VECTOR<T,3>& x1, const VECTOR<T,3>& x2, const VECTOR<T,3>& x3)
{ return Cross_Product(x1 - x3, x2 - x3); }

template< class T >
inline VECTOR<T,1>
Area_Weighted_Normal(const VECTOR< VECTOR<T,1>, 1 >& x)
{ return Area_Weighted_Normal(x[1]); }

template< class T >
inline VECTOR<T,2>
Area_Weighted_Normal(const VECTOR< VECTOR<T,2>, 2 >& x)
{ return Area_Weighted_Normal(x[1], x[2]); }

template< class T >
inline VECTOR<T,3>
Area_Weighted_Normal(const VECTOR< VECTOR<T,3>, 3 >& x)
{ return Area_Weighted_Normal(x[1], x[2], x[3]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_AREA_WEIGHTED_NORMAL_HPP
