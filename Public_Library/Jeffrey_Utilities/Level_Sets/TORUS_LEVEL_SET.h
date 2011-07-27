//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TORUS_LEVEL_SET_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TORUS_LEVEL_SET_HPP

#include <cmath>

#include <boost/math/special_functions/pow.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T >
struct TORUS_LEVEL_SET
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        TORUS_LEVEL_SET, (( typename T const, minor_radius ))
    )
public:
    typedef T result_type;
    T operator()(const VECTOR<T,3>& x) const
    {
        const T r = std::sqrt(boost::math::pow<2>(x[1]) + boost::math::pow<2>(x[2]));
        return std::sqrt(boost::math::pow<2>(r - 1) + boost::math::pow<2>(x[3])) - minor_radius;
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_TORUS_LEVEL_SET_HPP
