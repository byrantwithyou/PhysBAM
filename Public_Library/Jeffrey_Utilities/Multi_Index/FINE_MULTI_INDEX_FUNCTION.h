//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FINE_MULTI_INDEX_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FINE_MULTI_INDEX_FUNCTION_HPP

#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int FINE_FACTOR >
struct FINE_MULTI_INDEX_FUNCTION
{
    template<class> struct result;
    template< class T_THIS, class T_MULTI_INDEX >
    struct result< T_THIS ( T_MULTI_INDEX ) >
        : REMOVE_QUALIFIERS< T_MULTI_INDEX >
    { };

    template< int D >
    VECTOR<int,D> operator()(const VECTOR<int,D>& coarse_multi_index) const
    { return FINE_FACTOR * coarse_multi_index - (FINE_FACTOR - 1); }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FINE_MULTI_INDEX_FUNCTION_HPP
