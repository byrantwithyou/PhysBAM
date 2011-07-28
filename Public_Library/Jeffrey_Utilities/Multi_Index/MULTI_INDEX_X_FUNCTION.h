//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_MULTI_INDEX_BOX >
struct MULTI_INDEX_X_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        MULTI_INDEX_X_FUNCTION,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, min_x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, max_x ))
        (( typename T_MULTI_INDEX_BOX const, multi_index_box ))
    )
public:

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef VECTOR<T,D> result_type;
    result_type operator()(const int linear_index) const
    { return operator()(multi_index_box.Multi_Index(linear_index)); }
    result_type operator()(const MULTI_INDEX_TYPE& multi_index) const
    { return Multi_Index_X(min_x, max_x, multi_index_box, multi_index); }
};

template< class T, int D, class T_MULTI_INDEX_BOX >
inline MULTI_INDEX_X_FUNCTION< T, D, T_MULTI_INDEX_BOX >
Make_Multi_Index_X_Function(
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const T_MULTI_INDEX_BOX& multi_index_box)
{ return MULTI_INDEX_X_FUNCTION< T, D, T_MULTI_INDEX_BOX >(min_x, max_x, multi_index_box); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_FUNCTION_HPP
