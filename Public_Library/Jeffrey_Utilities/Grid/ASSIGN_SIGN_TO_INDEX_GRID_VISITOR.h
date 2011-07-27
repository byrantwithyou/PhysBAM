//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_ASSIGN_SIGN_TO_INDEX_GRID_VISITOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_ASSIGN_SIGN_TO_INDEX_GRID_VISITOR_HPP

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< class T_SIGN_OF_INDEX >
struct ASSIGN_SIGN_TO_INDEX_GRID_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        ASSIGN_SIGN_TO_INDEX_GRID_VISITOR,
        (( typename T_SIGN_OF_INDEX const, sign_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const int index, const int sign) const
    {
        typedef typename REMOVE_QUALIFIERS<
            typename RESULT_OF< const T_SIGN_OF_INDEX ( int ) >::type
        >::type SIGN_TYPE;
        sign_of_index(index) = static_cast< SIGN_TYPE >(sign);
    }
};

template< class T_SIGN_OF_INDEX >
inline ASSIGN_SIGN_TO_INDEX_GRID_VISITOR< T_SIGN_OF_INDEX >
Make_Assign_Sign_To_Index_Grid_Visitor(const T_SIGN_OF_INDEX& sign_of_index)
{ return ASSIGN_SIGN_TO_INDEX_GRID_VISITOR< T_SIGN_OF_INDEX >(sign_of_index); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GRID_ASSIGN_SIGN_TO_INDEX_GRID_VISITOR_HPP
