//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_CONSTANT_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_CONSTANT_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T >
struct CONSTANT_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        CONSTANT_FUNCTION, (( typename T const, x ))
    )
public:
    typedef T result_type;
    template< class T0 >
    T operator()(const T0&) const
    { return x; }
    template< class T0, class T1 >
    T operator()(const T0&, const T1&) const
    { return x; }
};

template< class T >
inline CONSTANT_FUNCTION<T>
Make_Constant_Function(const T& x)
{ return CONSTANT_FUNCTION<T>(x); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_CONSTANT_FUNCTION_HPP
