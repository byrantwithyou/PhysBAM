//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_HPP

#include <boost/config.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

namespace PhysBAM
{

template< class T, int D, class T_MULTI_INDEX_BOX >
inline VECTOR<T,D>
Multi_Index_X(
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const T_MULTI_INDEX_BOX& multi_index_box,
    const VECTOR<int,D>& multi_index)
{
    const VECTOR<int,D> min_multi_index = multi_index_box.Min_Multi_Index();
    const VECTOR<int,D> max_multi_index = multi_index_box.Max_Multi_Index();
#if defined( BOOST_MSVC ) && _MSC_VER <= 1500
    // workaround for MSVC9 ICE
    const VECTOR<T,D> temp = min_x * (max_multi_index - multi_index)
                           + max_x * (multi_index - min_multi_index);
    const VECTOR<T,D> result = temp / (max_multi_index - min_multi_index);
    return result;
#else // #if defined( BOOST_MSVC ) && _MSC_VER <= 1500
    return (min_x * (max_multi_index - multi_index)
          + max_x * (multi_index - min_multi_index))
         / (max_multi_index - min_multi_index);
#endif // #if defined( BOOST_MSVC ) && _MSC_VER <= 1500
}

template< class T, int D, class T_MULTI_INDEX_BOX >
inline VECTOR<T,D>
Multi_Index_X(
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const T_MULTI_INDEX_BOX& multi_index_box,
    const int linear_index)
{ return Multi_Index_X(min_x, max_x, multi_index_box, multi_index_box.Multi_Index(linear_index)); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_X_HPP
