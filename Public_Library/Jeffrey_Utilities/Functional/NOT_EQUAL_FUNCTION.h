//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NOT_EQUAL_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NOT_EQUAL_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

struct NOT_EQUAL_FUNCTION
{
    typedef bool result_type;
    template< class T, class U >
    bool operator()(const T& x, const U& y) const
    { return x != y; }
};

template< class U >
struct NOT_EQUAL1_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        NOT_EQUAL1_FUNCTION, (( typename U const, y ))
    )
public:
    typedef bool result_type;
    template< class T >
    bool operator()(const T& x) const
    { return x != y; }
};

template< class U >
inline NOT_EQUAL1_FUNCTION<U>
Make_Not_Equal_Function(const U& y)
{ return NOT_EQUAL1_FUNCTION<U>(y); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NOT_EQUAL_FUNCTION_HPP
