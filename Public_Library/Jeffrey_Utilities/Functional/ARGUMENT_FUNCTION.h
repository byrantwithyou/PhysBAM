//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARGUMENT_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARGUMENT_FUNCTION_HPP

#include <boost/function_types/parameter_types.hpp>
#include <boost/mpl/at.hpp>

namespace PhysBAM
{

template< int N >
struct ARGUMENT_FUNCTION;

template<>
struct ARGUMENT_FUNCTION<1>
{
    template< class T_SIGNATURE >
    struct result
        : boost::mpl::at_c<
              typename boost::function_types::parameter_types< T_SIGNATURE >::type,
              0
          >
    { };

    template< class T1 >
    T1& operator()(T1& x1) const
    { return x1; }

    template< class T1 >
    const T1& operator()(const T1& x1) const
    { return x1; }

    template< class T1, class T2 >
    T1& operator()(T1& x1, const T2& /*x2*/) const
    { return x1; }

    template< class T1, class T2 >
    const T1& operator()(const T1& x1, const T2& /*x2*/) const
    { return x1; }
};

template<>
struct ARGUMENT_FUNCTION<2>
{
    template< class T_SIGNATURE >
    struct result
        : boost::mpl::at_c<
              typename boost::function_types::parameter_types< T_SIGNATURE >::type,
              1
          >
    { };

    template< class T1, class T2 >
    T2& operator()(const T1& /*x1*/, T2& x2) const
    { return x2; }

    template< class T1, class T2 >
    const T2& operator()(const T1& /*x1*/, const T2& x2) const
    { return x2; }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_ARGUMENT_FUNCTION_HPP
