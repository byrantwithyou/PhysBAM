//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EQUAL_RELATIVE_TOLERANCE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EQUAL_RELATIVE_TOLERANCE_HPP

#include <cmath>

#include <algorithm>
#include <limits>

namespace PhysBAM
{

template< unsigned int EPS_FACTOR, class T >
inline bool
Equal_Relative_Tolerance(const T& x, const T& y)
{
    static const T eps = std::numeric_limits<T>::epsilon();
    return std::abs(x - y) <= EPS_FACTOR * eps * std::max(std::abs(x), std::abs(y));
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EQUAL_RELATIVE_TOLERANCE_HPP
