//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_COUNT_SIGNED_CELLS_VIA_CELL_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_COUNT_SIGNED_CELLS_VIA_CELL_SIGN_HPP

#include <cassert>

#include <Jeffrey_Utilities/Algorithm/Count_If.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>

namespace PhysBAM
{

template< class T_SIGN_OF_CELL_INDEX >
inline int
Count_Signed_Cells_Via_Cell_Sign(
    const int n_cell,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const int sign_to_count)
{
    return Count_If(
        1, n_cell,
        Make_Compose_Function(
            Make_Equal_Function(sign_to_count),
            sign_of_cell_index
        )
    );
}

template< class T_SIGN_OF_CELL_INDEX >
inline int
Count_Signed_Cells_Via_Cell_Sign_MT(
    const unsigned int n_thread,
    const int n_cell,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const int sign_to_count)
{
    assert(n_thread >= 1);
    return Count_If_MT(
        n_thread,
        1, n_cell,
        Make_Compose_Function(
            Make_Equal_Function(sign_to_count),
            sign_of_cell_index
        )
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_COUNT_SIGNED_CELLS_VIA_CELL_SIGN_HPP
