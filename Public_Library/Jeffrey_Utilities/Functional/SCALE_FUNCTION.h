//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_SCALE_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_SCALE_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>

namespace PhysBAM
{

template< class T >
struct SCALE_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SCALE_FUNCTION, (( typename T const, c ))
    )
public:

    template<class> struct result;
    template< class T_THIS, class U >
    struct result< T_THIS ( U ) >
        : REMOVE_QUALIFIERS<U>
    { };

    template< class U >
    U operator()(const U& x) const
    { return c * x; }
};

template< class T >
inline SCALE_FUNCTION<T>
Make_Scale_Function(const T& c)
{ return SCALE_FUNCTION<T>(c); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_SCALE_FUNCTION_HPP
