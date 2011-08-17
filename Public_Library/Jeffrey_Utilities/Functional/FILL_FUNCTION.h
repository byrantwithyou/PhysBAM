//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_FILL_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_FILL_FUNCTION_HPP

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T >
struct FILL1_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        FILL1_FUNCTION,
        (( typename T const, x ))
    )
public:
    typedef void result_type;
    template< class T_ARRAY >
    void operator()(T_ARRAY& a) const
    { Fill(a, x); }
};

template< class T >
inline FILL1_FUNCTION<T>
Make_Fill_Function(const T& x)
{ return FILL1_FUNCTION<T>(x); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_FILL_FUNCTION_HPP
