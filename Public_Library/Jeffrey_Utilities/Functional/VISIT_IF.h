//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_VISIT_IF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_VISIT_IF_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T_PREDICATE, class T_VISITOR >
struct VISIT_IF
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VISIT_IF,
        (( typename T_PREDICATE const, predicate ))
        (( typename T_VISITOR const, visitor ))
    )
public:
    typedef void result_type;
    template< class T >
    void operator()(const T& x) const
    {
        if(predicate(x))
            visitor(x);
    }
    template< class T >
    void operator()(T& x) const
    {
        if(predicate(x))
            visitor(x);
    }
};

template< class T_PREDICATE, class T_VISITOR >
inline VISIT_IF< T_PREDICATE, T_VISITOR >
Make_Visit_If(const T_PREDICATE& predicate, const T_VISITOR& visitor)
{ return VISIT_IF< T_PREDICATE, T_VISITOR >(predicate, visitor); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_VISIT_IF_HPP
