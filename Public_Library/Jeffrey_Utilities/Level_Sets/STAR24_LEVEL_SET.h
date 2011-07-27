//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_STAR24_LEVEL_SET_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_STAR24_LEVEL_SET_HPP

#include <cassert>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
struct STAR24_LEVEL_SET
{
    explicit STAR24_LEVEL_SET(const T min_radius)
        : m_min_radius(min_radius),
          m_max_radius(2 * min_radius)
    { }

    STAR24_LEVEL_SET(const T min_radius, const T max_radius)
        : m_min_radius(min_radius),
          m_max_radius(max_radius)
    { }

    typedef T result_type;

    T operator()(const VECTOR<T,3> x) const
    {
        const int axis = x.Dominant_Axis();
        if(x[axis] == 0)
            return -m_min_radius;
        const VECTOR<T,2> s = x.Remove_Index(axis) / std::abs(x[axis]);
        T s1 = s[1], s2 = s[2];
        assert(-1 <= s1 && s1 <= +1);
        assert(-1 <= s2 && s2 <= +1);
        static const T pi = boost::math::constants::pi<T>();
        s1 = (s1 + static_cast<T>(std::sin(pi * s1 / 2))) / 2;
        s2 = (s2 + static_cast<T>(std::sin(pi * s2 / 2))) / 2;
        const T t = static_cast<T>((1 - std::cos(2 * pi * s1)) * (1 - std::cos(2 * pi * s2)) / 4);
        assert(0 <= t && t <= 1);
        return x.Magnitude() - (m_min_radius + (m_max_radius - m_min_radius) * t);
    }

private:
    const T m_min_radius;
    const T m_max_radius;
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_STAR24_LEVEL_SET_HPP
