//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D, class T_MAC_VECTOR_FIELD >
struct DIVERGENCE_OF_MAC_VECTOR_FIELD
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        DIVERGENCE_OF_MAC_VECTOR_FIELD,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx ))
        (( typename T_MAC_VECTOR_FIELD const, mac_vector_field ))
    )
public:
    typedef T result_type;
    T operator()(VECTOR<int,D> multi_index) const
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
};

template< class T, int D, class T_MAC_VECTOR_FIELD >
inline DIVERGENCE_OF_MAC_VECTOR_FIELD< T, D, T_MAC_VECTOR_FIELD >
Make_Divergence_Of_MAC_Vector_Field(
    const VECTOR<T,D>& dx,
    const T_MAC_VECTOR_FIELD& mac_vector_field)
{ return DIVERGENCE_OF_MAC_VECTOR_FIELD< T, D, T_MAC_VECTOR_FIELD >(dx, mac_vector_field); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_HPP
