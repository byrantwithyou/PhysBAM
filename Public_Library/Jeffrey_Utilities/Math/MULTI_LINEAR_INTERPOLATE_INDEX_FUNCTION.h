//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Math/Multi_Linear_Interpolate.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< int FINE_FACTOR, class T_F_OF_INDEX >
struct MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION,
        (( typename T_F_OF_INDEX const, f_of_index ))
    )
public:
    template<class> struct result;
    template< class T_THIS, class T_MULTI_INDEX >
    struct result< T_THIS ( T_MULTI_INDEX ) >
        : REMOVE_QUALIFIERS<
              typename RESULT_OF< const T_F_OF_INDEX ( T_MULTI_INDEX ) >::type
          >
    { };

    template< int D >
    typename result< const MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION ( VECTOR<int,D> ) >::type
    operator()(const VECTOR<int,D>& fine_multi_index) const
    { return Multi_Linear_Interpolate< FINE_FACTOR >(f_of_index, fine_multi_index); }
};

template< int FINE_FACTOR, class T_F_OF_INDEX >
inline MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION< FINE_FACTOR, T_F_OF_INDEX >
Make_Multi_Linear_Interpolate_Index_Function(const T_F_OF_INDEX& f_of_index)
{ return MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION< FINE_FACTOR, T_F_OF_INDEX >(f_of_index); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_MULTI_LINEAR_INTERPOLATE_INDEX_FUNCTION_HPP
