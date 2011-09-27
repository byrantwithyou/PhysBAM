//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Example 3, Case I: a = b = 0.40178, m = 2, n = 6
// Example 3, Case II: a = 0.50012563, b = 0.250012563, m = 0, n = 12
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_3_CURVE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_3_CURVE_HPP

#include <cassert>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T >
struct CHEN_STRAIN_EXAMPLE_3_CURVE
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        CHEN_STRAIN_EXAMPLE_3_CURVE,
        (( typename T const, a ))
        (( typename T const, b ))
        (( /******/ int const, m ))
        (( /******/ int const, n ))
    )
public:
    typedef VECTOR<T,2> result_type;
    VECTOR<T,2> operator()(T t) const
    {
        t *= 2 * boost::math::constants::pi<T>();
        const T r = a + b * std::cos(m*t) * std::sin(n*t);
        return VECTOR<T,2>(r * std::cos(t), r * std::sin(t));
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_3_CURVE_HPP
