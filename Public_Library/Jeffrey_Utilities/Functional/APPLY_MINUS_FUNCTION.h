//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_MINUS_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_MINUS_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< class T_F, class T_G >
struct APPLY_MINUS_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_MINUS_FUNCTION,
        (( typename T_F const, f ))
        (( typename T_G const, g ))
    )
public:

    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
        : RESULT_OF< const T_F ( T ) >
    { };
    template< class T_THIS, class T, class U >
    struct result< T_THIS ( T, U ) >
        : RESULT_OF< const T_F ( T, U ) >
    { };

    template< class T >
    typename result< const APPLY_MINUS_FUNCTION ( const T& ) >::type
    operator()(const T& x) const
    { return f(x) - g(x); }
    template< class T, class U >
    typename result< const APPLY_MINUS_FUNCTION ( const T&, const U& ) >::type
    operator()(const T& x, const U& y) const
    { return f(x,y) - g(x,y); }
};

template< class T_F, class T_G >
inline APPLY_MINUS_FUNCTION< T_F, T_G >
Make_Apply_Minus_Function(const T_F& f, const T_G& g)
{ return APPLY_MINUS_FUNCTION< T_F, T_G >(f, g); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_MINUS_FUNCTION_HPP
