//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_FINE_VERTEX_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_FINE_VERTEX_SIGN_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX >
inline int
Vertex_Sign_Via_Fine_Vertex_Sign(
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    T_SIGN_OF_FINE_INDEX sign_of_fine_index,
    const VECTOR<int,D>& multi_index,
    const int sign_of_zero = -1)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    const MULTI_INDEX_TYPE fine_base_multi_index = FINE_FACTOR * multi_index - (FINE_FACTOR - 1);
    int sign = 2;
    BOOST_FOREACH(
        const MULTI_INDEX_TYPE fine_multi_index,
        Multi_Index_Box_Intersect(
            MULTI_INDEX_CUBE< D, -FINE_FACTOR, +FINE_FACTOR >(fine_base_multi_index),
            FINE_FACTOR * multi_index_bound - (FINE_FACTOR - 1)
        )
    ) {
        const int fine_sign = sign_of_fine_index(fine_multi_index);
        assert(-1 <= fine_sign && fine_sign <= +1);
        if(fine_sign == sign || fine_sign == 0 && sign == sign_of_zero)
            continue;
        if(sign == 2 || sign == 0 && fine_sign == sign_of_zero) {
            sign = fine_sign;
            continue;
        }
        return 0;
    }
    assert(sign != 2);
    return sign;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_FINE_VERTEX_SIGN_HPP
