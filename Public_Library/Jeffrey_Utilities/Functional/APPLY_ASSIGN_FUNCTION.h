//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_ASSIGN_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_ASSIGN_FUNCTION_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T_F, class T_G >
struct APPLY_ASSIGN_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_ASSIGN_FUNCTION,
        (( typename T_F const, f ))
        (( typename T_G const, g ))
    )
public:
    typedef void result_type;
    template< class T >
    void operator()(const T& x) const
    { f(x) = g(x); }
};

template< class T_F, class T_G >
inline APPLY_ASSIGN_FUNCTION< T_F, T_G >
Make_Apply_Assign_Function(const T_F& f, const T_G& g)
{ return APPLY_ASSIGN_FUNCTION< T_F, T_G >(f, g); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_APPLY_ASSIGN_FUNCTION_HPP
