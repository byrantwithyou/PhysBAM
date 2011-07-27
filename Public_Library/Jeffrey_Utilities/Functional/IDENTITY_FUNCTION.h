//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IDENTITY_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IDENTITY_FUNCTION_HPP

#include <boost/type_traits/remove_cv.hpp>

namespace PhysBAM
{

struct IDENTITY_FUNCTION
{
    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
        : boost::remove_cv<T>
    { };

    template< class T >
    const T& operator()(const T& x) const
    { return x; }

    template< class T >
    T& operator()(T& x) const
    { return x; }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IDENTITY_FUNCTION_HPP
