//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IF_ELSE_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IF_ELSE_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< class T_PREDICATE, class T_TRUE, class T_FALSE >
struct IF_ELSE_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        IF_ELSE_FUNCTION,
        (( typename T_PREDICATE const, predicate ))
        (( typename T_TRUE const, true_ ))
        (( typename T_FALSE const, false_ ))
    )
public:

    template<class> struct result;
    template< class T_THIS, class T >
    struct result< T_THIS ( T ) >
        : RESULT_OF< const T_TRUE ( T ) >
    { };

    template< class T >
    typename RESULT_OF< const T_TRUE ( const T& ) >::type
    operator()(const T& x) const
    { return predicate(x) ? true_(x) : false_(x); }
};

template< class T_PREDICATE, class T_TRUE, class T_FALSE >
inline IF_ELSE_FUNCTION< T_PREDICATE, T_TRUE, T_FALSE >
Make_If_Else_Function(const T_PREDICATE& predicate, const T_TRUE& true_, const T_FALSE& false_)
{ return IF_ELSE_FUNCTION< T_PREDICATE, T_TRUE, T_FALSE >(predicate, true_, false_); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_IF_ELSE_FUNCTION_HPP
