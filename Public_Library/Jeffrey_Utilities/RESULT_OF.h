//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Unfortunately and predictably, PhysBAM arrays don't confirm to the
// Boost.ResultOf protocol.
//##################################################################### 

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_RESULT_OF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_RESULT_OF_HPP

#include <boost/function_types/result_type.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/utility/result_of.hpp>

#include <Jeffrey_Utilities/HAS_CONST_RESULT_TYPE.h>
#include <Jeffrey_Utilities/HAS_RESULT_TYPE.h>

namespace PhysBAM
{

namespace Detail_RESULT_OF
{

template< class T_SIGNATURE, class T_F, bool >
struct RESULT_OF_DISPATCH_ON_HAS_RESULT_TYPE
    : boost::result_of< T_SIGNATURE >
{ };
template< class T_SIGNATURE, class T_F >
struct RESULT_OF_DISPATCH_ON_HAS_RESULT_TYPE< T_SIGNATURE, T_F, true >
{ typedef typename T_F::RESULT_TYPE type; };

template< class T_SIGNATURE, class T_F, bool >
struct RESULT_OF_DISPATCH_ON_HAS_CONST_RESULT_TYPE
    : boost::result_of< T_SIGNATURE >
{ };
template< class T_SIGNATURE, class T_F >
struct RESULT_OF_DISPATCH_ON_HAS_CONST_RESULT_TYPE< T_SIGNATURE, T_F, true >
{ typedef typename T_F::CONST_RESULT_TYPE type; };

template< class T_SIGNATURE, class T_F >
struct RESULT_OF_DISPATCH
    : RESULT_OF_DISPATCH_ON_HAS_RESULT_TYPE<
          T_SIGNATURE, T_F,
          HAS_RESULT_TYPE< T_F >::value
      >
{ };

template< class T_SIGNATURE, class T_F >
struct RESULT_OF_DISPATCH< T_SIGNATURE, const T_F >
    : RESULT_OF_DISPATCH_ON_HAS_CONST_RESULT_TYPE<
          T_SIGNATURE, T_F,
          HAS_CONST_RESULT_TYPE< T_F >::value
      >
{ };

} // namespace Detail_RESULT_OF

template< class T_SIGNATURE >
struct RESULT_OF
    : Detail_RESULT_OF::RESULT_OF_DISPATCH<
          T_SIGNATURE,
          typename boost::function_types::result_type< T_SIGNATURE >::type
      >
{ };

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_RESULT_OF_HPP
