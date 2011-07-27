//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_FIND_ALL_SIGNED_VERTICES_VIA_CELL_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_FIND_ALL_SIGNED_VERTICES_VIA_CELL_SIGN_HPP

#include <cassert>

#include <Jeffrey_Utilities/Algorithm/Find_All_If.h>
#include <Jeffrey_Utilities/Grid/VERTEX_HAS_SIGN_VIA_CELL_SIGN.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int D, class T_SIGN_OF_CELL_INDEX >
inline void
Find_All_Signed_Vertices_Via_Cell_Sign(
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const int sign_to_find,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    ARRAY<int>& result)
{
    Find_All_If(
        1, multi_index_bound.Size(),
        Make_Vertex_Has_Sign_Via_Cell_Sign(
            multi_index_bound, sign_to_find, sign_of_cell_index
        ),
        result
    );
}

template< int D, class T_SIGN_OF_CELL_INDEX >
inline void
Find_All_Signed_Vertices_Via_Cell_Sign_MT(
    const unsigned int n_thread,
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const int sign_to_find,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    ARRAY<int>& result)
{
    assert(n_thread >= 1);
    Find_All_If_MT(
        n_thread,
        1, multi_index_bound.Size(),
        Make_Vertex_Has_Sign_Via_Cell_Sign(
            multi_index_bound, sign_to_find, sign_of_cell_index
        ),
        result
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_FIND_ALL_SIGNED_VERTICES_VIA_CELL_SIGN_HPP
