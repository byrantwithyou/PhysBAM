//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CELL2_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CELL2_HPP

#include <cassert>
#include <cmath>

#include <limits>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/Geometry/Divide_Cube2.h>
#include <Jeffrey_Utilities/Level_Sets/Shift_Level_Set_Away_From_Vertices.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_PHI_OF_FINE_INDEX, class T_POLYTOPE_VISITOR >
void
Divide_Cell2(
    const VECTOR<T,D>& dx,
    const VECTOR<int,D>& cell_multi_index,
    T_PHI_OF_FINE_INDEX phi_of_fine_index,
    const T_POLYTOPE_VISITOR& polytope_visitor,
    const float min_dist_to_vertex = 0.0f,
    const int sign_of_zero = -1)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    assert(0 <= min_dist_to_vertex && min_dist_to_vertex <= 0.5f);
    assert(sign_of_zero == -1 || sign_of_zero == +1);

    VECTOR< T, STATIC_POW_C<3,D>::value > phi_of_cube2_vertex; // init'ed to 0
    const MULTI_INDEX_CUBE<D,-1,+1> fine_local_multi_index_cube(2 * cell_multi_index);
    for(int cube2_vertex = 1; cube2_vertex <= STATIC_POW_C<3,D>::value; ++cube2_vertex) {
        const MULTI_INDEX_TYPE fine_multi_index = fine_local_multi_index_cube.Multi_Index(cube2_vertex);
        phi_of_cube2_vertex[cube2_vertex] = phi_of_fine_index(fine_multi_index);
    }

    if(min_dist_to_vertex != 0)
        Shift_Level_Set_Away_From_Vertices(
            STATIC_MULTI_INDEX_CUBE<D,-1,+1>(),
            Make_Array_View(STATIC_POW_C<3,D>::value, &phi_of_cube2_vertex[1]),
            min_dist_to_vertex,
            sign_of_zero
        );

    Divide_Cube2(
        dx,
        Make_Const_Array_View(STATIC_POW_C<3,D>::value, &phi_of_cube2_vertex[1]),
        polytope_visitor,
        sign_of_zero
    );
}

template< class T, int D, class T_PHI_OF_FINE_INDEX, class T_POLYTOPE_VISITOR >
inline void
Divide_Cell2(
    const VECTOR<T,D>& dx,
    const VECTOR<int,D>& cell_multi_index,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const T_POLYTOPE_VISITOR& polytope_visitor,
    const int sign_of_zero)
{
    Divide_Cell2(
        dx,
        cell_multi_index,
        phi_of_fine_index,
        polytope_visitor,
        0.0f, // min_dist_to_vertex
        sign_of_zero
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CELL2_HPP
