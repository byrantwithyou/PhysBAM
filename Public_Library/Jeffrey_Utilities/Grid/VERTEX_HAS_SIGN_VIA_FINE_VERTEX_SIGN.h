//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Grid/Vertex_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

namespace PhysBAM
{

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX >
struct VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN,
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( /******/ int const, sign_to_count ))
        (( typename T_SIGN_OF_FINE_INDEX const, sign_of_fine_index ))
        (( /******/ int const, sign_of_zero ))
    )
public:
    typedef bool result_type;
    bool operator()(const int linear_index) const
    { return operator()(multi_index_bound.Multi_Index(linear_index)); }
    bool operator()(const VECTOR<int,D>& multi_index) const
    {
        const int sign = Vertex_Sign_Via_Fine_Vertex_Sign< FINE_FACTOR >(
            multi_index_bound, sign_of_fine_index, multi_index, sign_of_zero
        );
        return sign == sign_to_count;
    }
};

template< int FINE_FACTOR, int D, class T_SIGN_OF_FINE_INDEX >
inline VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN< FINE_FACTOR, D, T_SIGN_OF_FINE_INDEX >
Make_Vertex_Has_Sign_Via_Fine_Vertex_Sign(
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const int sign_to_count,
    const T_SIGN_OF_FINE_INDEX& sign_of_fine_index,
    const int sign_of_zero)
{
    return VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN< FINE_FACTOR, D, T_SIGN_OF_FINE_INDEX >(
        multi_index_bound, sign_to_count, sign_of_fine_index, sign_of_zero
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VERTEX_HAS_SIGN_VIA_FINE_VERTEX_SIGN_HPP
