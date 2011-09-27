//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_NORMAL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_NORMAL_HPP

#include <cmath>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Geometry/Area_Weighted_Normal.h>

namespace PhysBAM
{

template< class T >
inline VECTOR<T,1>
Normal(const VECTOR<T,1>& /*x1*/)
{ return VECTOR<T,1>(1); }

template< class T >
inline VECTOR<T,2>
Normal(const VECTOR<T,2>& x1, const VECTOR<T,2>& x2)
{
    const VECTOR<T,2> awn = Area_Weighted_Normal(x1, x2);
    const T a2 = awn.Magnitude_Squared();
    return a2 == 0 ? VECTOR<T,2>() : awn / std::sqrt(a2);
}

template< class T >
inline VECTOR<T,3>
Normal(const VECTOR<T,3>& x1, const VECTOR<T,3>& x2, const VECTOR<T,3>& x3)
{
    const VECTOR<T,3> awn = Area_Weighted_Normal(x1, x2, x3);
    const T a2 = awn.Magnitude_Squared();
    return a2 == 0 ? VECTOR<T,3>() : awn / std::sqrt(a2);
}

template< class T >
inline VECTOR<T,1>
Normal(const VECTOR< VECTOR<T,1>, 1 >& x)
{ return Normal(x[1]); }

template< class T >
inline VECTOR<T,2>
Normal(const VECTOR< VECTOR<T,2>, 2 >& x)
{ return Normal(x[1], x[2]); }

template< class T >
inline VECTOR<T,3>
Normal(const VECTOR< VECTOR<T,3>, 3 >& x)
{ return Normal(x[1], x[2], x[3]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_NORMAL_HPP
