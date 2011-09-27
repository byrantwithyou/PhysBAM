//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// http://en.wikipedia.org/wiki/Trefoil_knot
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_TREFOIL_KNOT_CURVE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_TREFOIL_KNOT_CURVE_HPP

#include <cassert>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>

namespace PhysBAM
{

struct TREFOIL_KNOT_CURVE
{
    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
    { typedef VECTOR< typename REMOVE_QUALIFIERS<T>::type, 3 > type; };

    template< class T >
    VECTOR<T,3> operator()(T t) const
    {
        t *= 2 * boost::math::constants::pi<T>();
        return VECTOR<T,3>(
            (2 + std::cos(3*t)) * std::cos(2*t),
            (2 + std::cos(3*t)) * std::sin(2*t),
            std::sin(3*t)
        );
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_TREFOIL_KNOT_CURVE_HPP
