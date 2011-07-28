//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TREFOIL_KNOT_LEVEL_SET_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TREFOIL_KNOT_LEVEL_SET_HPP

#include <cassert>

#include <Jeffrey_Utilities/Level_Sets/Min_Dist_To_Parametric_Curve.h>
#include <Jeffrey_Utilities/Math/Equal_Relative_Tolerance.h>
#include <Jeffrey_Utilities/Parametric_Curves/TREFOIL_KNOT_CURVE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
struct TREFOIL_KNOT_LEVEL_SET
{
    explicit TREFOIL_KNOT_LEVEL_SET(const T minor_radius)
        : m_minor_radius(minor_radius),
          m_cached_x(Cached_X())
    { }

    typedef T result_type;

    T operator()(const VECTOR<T,3> x) const
    {
        const T min_dist = Min_Dist_To_Parametric_Curve_Via_Cached_X<2>(x, TREFOIL_KNOT_CURVE(), m_cached_x);
#if 0
        const T other_min_dist = Min_Dist_To_Parametric_Curve_Via_Subintervals<6>(x, TREFOIL_KNOT_CURVE());
        assert(Equal_Relative_Tolerance<1024*1024>(min_dist, other_min_dist));
#endif // #if 1
        return min_dist - m_minor_radius;
    }

private:
    const T m_minor_radius;
    static const int N_CACHED_X = 48;
    const VECTOR<T,3> (&m_cached_x)[N_CACHED_X];
    static const VECTOR<T,3> (&Init_Cached_X())[N_CACHED_X]
    {
        static VECTOR<T,3> cached_x[N_CACHED_X];
        for(int i = 0; i != N_CACHED_X; ++i)
            cached_x[i] = TREFOIL_KNOT_CURVE()(static_cast<T>(i) / N_CACHED_X);
        return cached_x;
    }
    static const VECTOR<T,3> (&Cached_X())[N_CACHED_X]
    {
        static const VECTOR<T,3> (&cached_x)[N_CACHED_X] = Init_Cached_X();
        return cached_x;
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TREFOIL_KNOT_LEVEL_SET_HPP
