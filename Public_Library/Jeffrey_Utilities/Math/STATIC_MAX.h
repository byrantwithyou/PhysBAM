//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_MAX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_MAX_HPP

namespace PhysBAM
{

template< int M, int N >
struct STATIC_MAX_C
{
    typedef STATIC_MAX_C type;
    static const int value = M < N ? N : M;
};

template< class M, class N >
struct STATIC_MAX
    : STATIC_MAX_C< M::value, N::value >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_MAX_HPP
