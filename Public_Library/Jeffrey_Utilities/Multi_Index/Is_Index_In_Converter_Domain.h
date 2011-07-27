//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_IS_INDEX_IN_CONVERTER_DOMAIN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_IS_INDEX_IN_CONVERTER_DOMAIN_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T_INDEX_CONVERTER >
inline bool
Is_Index_In_Converter_Domain(const T_INDEX_CONVERTER& index_converter, const int linear_index)
{ return 1 <= linear_index && linear_index <= index_converter.Size(); }

template< class T_INDEX_CONVERTER, int D >
inline bool
Is_Index_In_Converter_Domain(const T_INDEX_CONVERTER& index_converter, const VECTOR<int,D>& multi_index)
{ return index_converter.Contains(multi_index); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_IS_INDEX_IN_CONVERTER_DOMAIN_HPP
