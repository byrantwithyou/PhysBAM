//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Divergence_Of_MAC_Vector_Field.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_MAC_VECTOR_FIELD >
struct DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx ))
        (( typename T_MAC_VECTOR_FIELD const, mac_vector_field ))
    )
public:
    typedef T result_type;
    T operator()(const VECTOR<int,D>& multi_index) const
    { return Divergence_Of_MAC_Vector_Field(dx, mac_vector_field, multi_index); }
};

template< class T, int D, class T_MAC_VECTOR_FIELD >
inline DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION<
    T, D, T_MAC_VECTOR_FIELD
>
Make_Divergence_Of_MAC_Vector_Field_Function(
    const VECTOR<T,D>& dx,
    const T_MAC_VECTOR_FIELD& mac_vector_field)
{
    return DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION<
        T, D, T_MAC_VECTOR_FIELD
    >(dx, mac_vector_field);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DIVERGENCE_OF_MAC_VECTOR_FIELD_FUNCTION_HPP
