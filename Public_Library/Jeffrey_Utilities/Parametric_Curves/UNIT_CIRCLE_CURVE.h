//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_UNIT_CIRCLE_CURVE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_UNIT_CIRCLE_CURVE_HPP

#include <cassert>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>

namespace PhysBAM
{

struct UNIT_CIRCLE_CURVE
{
    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
    { typedef VECTOR< typename REMOVE_QUALIFIERS<T>::type, 2 > type; };

    template< class T >
    VECTOR<T,2> operator()(T t) const
    {
        t *= 2 * boost::math::constants::pi<T>();
        return VECTOR<T,2>(std::cos(t), std::sin(t));
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PARAMETRIC_CURVES_UNIT_CIRCLE_CURVE_HPP
