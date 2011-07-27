//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARRAY_WRAPPER_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARRAY_WRAPPER_FUNCTION_HPP

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/HAS_ISC_VALUE.h>

namespace PhysBAM
{

template< class T_ARRAY >
struct ARRAY_WRAPPER_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        ARRAY_WRAPPER_FUNCTION, (( typename T_ARRAY&, array ))
    )
public:
    typedef typename boost::mpl::if_<
        boost::is_const< T_ARRAY >,
        typename T_ARRAY::ELEMENT const,
        typename T_ARRAY::ELEMENT /***/
    >::type & result_type;

    template< class T_INT >
    typename boost::enable_if<
        HAS_ISC_VALUE< T_INT >,
        result_type
    >::type operator()(T_INT) const
    { return array(static_cast< typename T_ARRAY::INDEX >(T_INT::value)); }

    result_type operator()(typename T_ARRAY::INDEX const i) const
    { return array(i); }
};

template< class T_ARRAY >
inline ARRAY_WRAPPER_FUNCTION< T_ARRAY >
Make_Array_Wrapper_Function(T_ARRAY& array)
{ return ARRAY_WRAPPER_FUNCTION< T_ARRAY >(array); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARRAY_WRAPPER_FUNCTION_HPP
