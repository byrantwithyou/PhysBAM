//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIVERGENCE_OF_MAC_VECTOR_FIELD_H
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIVERGENCE_OF_MAC_VECTOR_FIELD_H

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_MAC_VECTOR_FIELD >
inline T
Divergence_Of_MAC_Vector_Field(
    const VECTOR<T,D>& dx,
    T_MAC_VECTOR_FIELD mac_vector_field,
    VECTOR<int,D> multi_index)
{
    T result = 0;
    for(int d = 1; d <= D; ++d) {
        ++multi_index[d];
        T diff = mac_vector_field(d, multi_index);
        --multi_index[d];
        diff -= mac_vector_field(d, multi_index);
        result += diff / dx[d];
    }
    return result;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIVERGENCE_OF_MAC_VECTOR_FIELD_H
