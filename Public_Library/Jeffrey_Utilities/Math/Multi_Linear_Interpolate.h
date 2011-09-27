//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_HPP

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< int FINE_FACTOR, class T_F_OF_INDEX, int D >
inline typename REMOVE_QUALIFIERS<
    typename RESULT_OF< T_F_OF_INDEX ( VECTOR<int,D> ) >::type
>::type
Multi_Linear_Interpolate(
    T_F_OF_INDEX f_of_index,
    const VECTOR<int,D> fine_multi_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef typename REMOVE_QUALIFIERS<
        typename RESULT_OF< T_F_OF_INDEX ( MULTI_INDEX_TYPE ) >::type
    >::type RESULT_TYPE;
    const MULTI_INDEX_TYPE base_coarse_multi_index = (fine_multi_index - 1) / FINE_FACTOR + 1;
    const MULTI_INDEX_TYPE fine_multi_offset = fine_multi_index - FINE_FACTOR * base_coarse_multi_index + (FINE_FACTOR - 1);
    RESULT_TYPE result = RESULT_TYPE();
    BOOST_FOREACH(
        const MULTI_INDEX_TYPE coarse_multi_offset,
        (STATIC_MULTI_INDEX_CUBE<D,0,1>())
    ) {
        const MULTI_INDEX_TYPE coarse_multi_index = base_coarse_multi_index + coarse_multi_offset;
        const int unscaled_weight = ((1 - coarse_multi_offset) * (FINE_FACTOR - fine_multi_offset) + coarse_multi_offset * fine_multi_offset).Product();
        result += unscaled_weight * f_of_index(coarse_multi_index);
    }
    return result / STATIC_POW_C< FINE_FACTOR, D >::value;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_HPP
