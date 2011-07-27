//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ZERO_MPL_VECTOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ZERO_MPL_VECTOR_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector/vector10.hpp>

namespace PhysBAM
{

template< int N >
struct ZERO_MPL_VECTOR_C
    : boost::mpl::push_back<
          typename ZERO_MPL_VECTOR_C< N-1 >::type,
          boost::mpl::int_<0>
      >
{ };

template<>
struct ZERO_MPL_VECTOR_C<0>
{ typedef boost::mpl::vector0<> type; };

template< class N >
struct ZERO_MPL_VECTOR
    : ZERO_MPL_VECTOR_C< N::value >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ZERO_MPL_VECTOR_HPP
