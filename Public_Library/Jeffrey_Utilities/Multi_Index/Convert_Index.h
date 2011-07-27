//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_CONVERT_INDEX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_CONVERT_INDEX_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Detail_Convert_Index
{

template< class T_DEST_INDEX >
struct CONVERT_INDEX_DISPATCH;

template<>
struct CONVERT_INDEX_DISPATCH< int >
{
    template< class T_INDEX_CONVERTER >
    static int Apply(const T_INDEX_CONVERTER&, const int srce_linear_index)
    { return srce_linear_index; }
    template< class T_INDEX_CONVERTER, int D >
    static int Apply(const T_INDEX_CONVERTER& index_converter, const VECTOR<int,D>& srce_multi_index)
    { return index_converter.Linear_Index(srce_multi_index); }
};

template< int D >
struct CONVERT_INDEX_DISPATCH< VECTOR<int,D> >
{
    template< class T_INDEX_CONVERTER >
    static VECTOR<int,D> Apply(const T_INDEX_CONVERTER& index_converter, const int srce_linear_index)
    { return index_converter.Multi_Index(srce_linear_index); }
    template< class T_INDEX_CONVERTER >
    static VECTOR<int,D> Apply(const T_INDEX_CONVERTER&, const VECTOR<int,D>& srce_multi_index)
    { return srce_multi_index; }
};

} // namespace Detail_Convert_Index

template< class T_DEST_INDEX, class T_INDEX_CONVERTER, class T_SRCE_INDEX >
inline T_DEST_INDEX
Convert_Index(const T_INDEX_CONVERTER& index_converter, const T_SRCE_INDEX& srce_index)
{ return Detail_Convert_Index::CONVERT_INDEX_DISPATCH< T_DEST_INDEX >::Apply(index_converter, srce_index); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_CONVERT_INDEX_HPP
