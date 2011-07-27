//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_UNIT_CUBE_LEVEL_SET_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_UNIT_CUBE_LEVEL_SET_HPP

#include <cmath>

#include <boost/math/special_functions/pow.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

struct UNIT_CUBE_LEVEL_SET
{
    template<class> struct result;
    template< class T_THIS, class T_VECTOR >
    struct result< T_THIS ( T_VECTOR ) >
    { typedef typename boost::remove_reference< T_VECTOR >::type::value_type type; };

    template< class T, int D >
    T operator()(const VECTOR<T,D>& x) const
    {
        const T max_abs = x.Max_Abs();
        if(max_abs <= 1)
            return max_abs - 1;
        T sq_dist = 0;
        for(int d = 1; d <= D; ++d)
            sq_dist += std::abs(x[d]) <= 1 ? 0 : boost::math::pow<2>(std::abs(x[d]) - 1);
        return std::sqrt(sq_dist);
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LEVEL_SETS_UNIT_CUBE_LEVEL_SET_HPP
