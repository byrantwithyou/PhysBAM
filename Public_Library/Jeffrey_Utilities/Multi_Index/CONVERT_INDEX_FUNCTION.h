//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_CONVERT_INDEX_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_CONVERT_INDEX_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/Convert_Index.h>

namespace PhysBAM
{

template< class T_INDEX_CONVERTER, class T_DEST_INDEX >
struct CONVERT_INDEX_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        CONVERT_INDEX_FUNCTION,
        (( typename T_INDEX_CONVERTER const, index_converter ))
    )
public:
    typedef T_DEST_INDEX result_type;
    template< class T_SRCE_INDEX >
    result_type operator()(const T_SRCE_INDEX& srce_index) const
    { return Convert_Index< result_type >(index_converter, srce_index); }
};

template< class T_DEST_INDEX, class T_INDEX_CONVERTER >
inline CONVERT_INDEX_FUNCTION< T_INDEX_CONVERTER, T_DEST_INDEX >
Make_Convert_Index_Function(const T_INDEX_CONVERTER& index_converter)
{ return CONVERT_INDEX_FUNCTION< T_INDEX_CONVERTER, T_DEST_INDEX >(index_converter); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_CONVERT_INDEX_FUNCTION_HPP
