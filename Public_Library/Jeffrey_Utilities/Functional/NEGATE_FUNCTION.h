//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NEGATE_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NEGATE_FUNCTION_HPP

#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>

namespace PhysBAM
{

struct NEGATE_FUNCTION
{
    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
        : REMOVE_QUALIFIERS<T>
    { };

    template< class T >
    T operator()(const T& x) const
    { return -x; }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NEGATE_FUNCTION_HPP
