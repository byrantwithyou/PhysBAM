//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_CELL_SIGN_VIA_FINE_VERTEX_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_CELL_SIGN_VIA_FINE_VERTEX_SIGN_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX >
inline int
Cell_Sign_Via_Fine_Vertex_Sign(
    T_SIGN_OF_FINE_INDEX sign_of_fine_index,
    const VECTOR<int,D>& cell_multi_index,
    const int sign_of_zero = -1)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    const MULTI_INDEX_TYPE fine_base_multi_index = FINE_FACTOR * cell_multi_index - (FINE_FACTOR - 1);
    int cell_sign = 2;
    BOOST_FOREACH(
        const MULTI_INDEX_TYPE fine_multi_index,
        (MULTI_INDEX_CUBE< D, 0, FINE_FACTOR >(fine_base_multi_index))
    ) {
        const int sign = sign_of_fine_index(fine_multi_index);
        assert(-1 <= sign && sign <= +1);
        if(sign == cell_sign || (sign == 0 && cell_sign == sign_of_zero))
            continue;
        if(cell_sign == 2 || (cell_sign == 0 && sign == sign_of_zero)) {
            cell_sign = sign;
            continue;
        }
        return 0;
    }
    assert(cell_sign != 2);
    return cell_sign;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_CELL_SIGN_VIA_FINE_VERTEX_SIGN_HPP
