//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_CELL_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_CELL_SIGN_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_Box_Intersect.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template< int D, class T_SIGN_OF_CELL_INDEX >
inline int
Vertex_Sign_Via_Cell_Sign(
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound,
    T_SIGN_OF_CELL_INDEX sign_of_cell_index,
    const VECTOR<int,D>& multi_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    int sign = 2;
    BOOST_FOREACH(
        const MULTI_INDEX_TYPE cell_multi_index,
        Multi_Index_Box_Intersect(MULTI_INDEX_CUBE<D,-1,0>(multi_index), cell_multi_index_bound)
    ) {
        const int cell_sign = sign_of_cell_index(cell_multi_index);
        assert(-1 <= cell_sign && cell_sign <= +1);
        if(cell_sign == sign)
            continue;
        if(sign == 2 && cell_sign != 0) {
            sign = cell_sign;
            continue;
        }
        return 0;
    }
    assert(sign != 2);
    return sign;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_SIGN_VIA_CELL_SIGN_HPP
