//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Example 4: a = 0.6012563, b = 0.2401256, n = 4
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_4_CURVE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_4_CURVE_HPP

#include <cassert>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T >
struct CHEN_STRAIN_EXAMPLE_4_CURVE
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        CHEN_STRAIN_EXAMPLE_4_CURVE,
        (( typename T const, a ))
        (( typename T const, b ))
        (( /******/ int const, n ))
    )
public:
    typedef VECTOR<T,2> result_type;
    VECTOR<T,2> operator()(T t) const
    {
        t *= 2 * boost::math::constants::pi<T>();
        const T theta = t + std::sin(n*t);
        const T r = a - b * std::sin(n*t);
        return VECTOR<T,2>(r * std::cos(theta), r * std::sin(theta));
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_CHEN_STRAIN_EXAMPLE_4_CURVE_HPP
