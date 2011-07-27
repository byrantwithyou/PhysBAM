//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_IF_SIGN_PREDICATE_GRID_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_IF_SIGN_PREDICATE_GRID_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

template< class T_SIGN_PREDICATE, class T_VISITOR >
struct VISIT_IF_SIGN_PREDICATE_GRID_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        VISIT_IF_SIGN_PREDICATE_GRID_VISITOR,
        (( typename T_SIGN_PREDICATE const, sign_predicate ))
        (( typename T_VISITOR const, visitor ))
    )
public:
    typedef void result_type;
    template< class T >
    void operator()(const T& x, const int sign) const
    {
        if(sign_predicate(sign))
            visitor(x);
    }
};

template< class T_SIGN_PREDICATE, class T_VISITOR >
inline VISIT_IF_SIGN_PREDICATE_GRID_VISITOR< T_SIGN_PREDICATE, T_VISITOR >
Make_Visit_If_Sign_Predicate_Grid_Visitor(const T_SIGN_PREDICATE& sign_predicate, const T_VISITOR& visitor)
{ return VISIT_IF_SIGN_PREDICATE_GRID_VISITOR< T_SIGN_PREDICATE, T_VISITOR >(sign_predicate, visitor); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_VISIT_IF_SIGN_PREDICATE_GRID_VISITOR_HPP
