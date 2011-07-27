//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_CELL_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_CELL_SIGN_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Grid/Vertex_Sign_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int D, class T_SIGN_OF_CELL_INDEX >
struct VERTEX_HAS_SIGN_VIA_CELL_SIGN
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VERTEX_HAS_SIGN_VIA_CELL_SIGN,
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( /******/ int const, sign_to_count ))
        (( typename T_SIGN_OF_CELL_INDEX const, sign_of_cell_index ))
    )
public:
    typedef bool result_type;
    bool operator()(const int linear_index) const
    {
        const VECTOR<int,D> multi_index = multi_index_bound.Multi_Index(linear_index);
        const int sign = Vertex_Sign_Via_Cell_Sign(
            multi_index_bound - 1, sign_of_cell_index, multi_index
        );
        return sign == sign_to_count;
    }
};

template< int D, class T_SIGN_OF_CELL_INDEX >
inline VERTEX_HAS_SIGN_VIA_CELL_SIGN< D, T_SIGN_OF_CELL_INDEX >
Make_Vertex_Has_Sign_Via_Cell_Sign(
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const int sign_to_count,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index)
{
    return VERTEX_HAS_SIGN_VIA_CELL_SIGN< D, T_SIGN_OF_CELL_INDEX >(
        multi_index_bound, sign_to_count, sign_of_cell_index
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_CELL_SIGN_HPP
