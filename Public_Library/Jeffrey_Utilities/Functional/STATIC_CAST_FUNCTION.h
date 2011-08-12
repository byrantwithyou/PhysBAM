//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_STATIC_CAST_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_STATIC_CAST_FUNCTION_HPP

namespace PhysBAM
{

template< class T >
struct STATIC_CAST_FUNCTION
{
    typedef T result_type;

    template< class U >
    T operator()(const U& x) const
    { return static_cast<T>(x); }

    template< class U >
    T operator()(U& x) const
    { return static_cast<T>(x); }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_STATIC_CAST_FUNCTION_HPP
