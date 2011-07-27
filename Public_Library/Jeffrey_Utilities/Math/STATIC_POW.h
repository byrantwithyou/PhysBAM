//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_POW_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_POW_HPP

namespace PhysBAM
{

template< int B, int P >
struct STATIC_POW_C
{
    typedef STATIC_POW_C type;
    static const int value = (P & 1 ? B : 1) * STATIC_POW_C< STATIC_POW_C< B, P/2 >::value, 2 >::value;
};
template< int B > struct STATIC_POW_C<B,0> { typedef STATIC_POW_C type; static const int value = 1; };
template< int B > struct STATIC_POW_C<B,1> { typedef STATIC_POW_C type; static const int value = B; };
template< int B > struct STATIC_POW_C<B,2> { typedef STATIC_POW_C type; static const int value = B*B; };
template< int B > struct STATIC_POW_C<B,3> { typedef STATIC_POW_C type; static const int value = B*B*B; };

template< class B, class P >
struct STATIC_POW
    : STATIC_POW_C< B::value, P::value >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_STATIC_POW_HPP
