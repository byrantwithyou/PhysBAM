//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_BOUND_FAST_MEM_FN_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_BOUND_FAST_MEM_FN_HPP

#include <boost/call_traits.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

#define PHYSBAM_BOUND_FAST_MEM_FN( X, MEM_PTR ) \
    ::PhysBAM::Detail_BOUND_FAST_MEM_FN::Deduce( MEM_PTR ).Apply< MEM_PTR >( X )
#define PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE( X, MEM_PTR ) \
    ::PhysBAM::Detail_BOUND_FAST_MEM_FN::Deduce( MEM_PTR ).template Apply< MEM_PTR >( X )

namespace PhysBAM
{

template< class T_MEM_PTR, T_MEM_PTR MEM_PTR >
struct BOUND_FAST_MEM_FN;

template< class T, class T_RESULT, class T1, T_RESULT (T::*MEM_PTR)( T1 ) >
struct BOUND_FAST_MEM_FN< T_RESULT (T::*)( T1 ), MEM_PTR >
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BOUND_FAST_MEM_FN, (( typename T&, x ))
    )
public:
    typedef T_RESULT result_type;
    result_type operator()(typename boost::call_traits< T1 >::param_type x1) const
    { return (x.*MEM_PTR)(x1); }
};

template< class T, class T_RESULT, class T1, T_RESULT (T::*MEM_PTR)( T1 ) const >
struct BOUND_FAST_MEM_FN< T_RESULT (T::*)( T1 ) const, MEM_PTR >
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BOUND_FAST_MEM_FN, (( typename T const &, x ))
    )
public:
    typedef T_RESULT result_type;
    result_type operator()(typename boost::call_traits< T1 >::param_type x1) const
    { return (x.*MEM_PTR)(x1); }
};

template< class T, class T_RESULT, class T1, class T2, T_RESULT (T::*MEM_PTR)( T1, T2 ) >
struct BOUND_FAST_MEM_FN< T_RESULT (T::*)( T1, T2 ), MEM_PTR >
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BOUND_FAST_MEM_FN, (( typename T&, x ))
    )
public:
    typedef T_RESULT result_type;
    result_type operator()(
        typename boost::call_traits< T1 >::param_type x1,
        typename boost::call_traits< T2 >::param_type x2) const
    { return (x.*MEM_PTR)(x1, x2); }
};

template< class T, class T_RESULT, class T1, class T2, T_RESULT (T::*MEM_PTR)( T1, T2 ) const >
struct BOUND_FAST_MEM_FN< T_RESULT (T::*)( T1, T2 ) const, MEM_PTR >
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BOUND_FAST_MEM_FN, (( typename T const &, x ))
    )
public:
    typedef T_RESULT result_type;
    result_type operator()(
        typename boost::call_traits< T1 >::param_type x1,
        typename boost::call_traits< T2 >::param_type x2) const
    { return (x.*MEM_PTR)(x1, x2); }
};

namespace Detail_BOUND_FAST_MEM_FN
{

template< class T_MEM_PTR >
struct MAKER
{
    template< T_MEM_PTR MEM_PTR, class T >
    BOUND_FAST_MEM_FN< T_MEM_PTR, MEM_PTR >
    static Apply(T& x)
    { return BOUND_FAST_MEM_FN< T_MEM_PTR, MEM_PTR >(x); }

    template< T_MEM_PTR MEM_PTR, class T >
    BOUND_FAST_MEM_FN< T_MEM_PTR, MEM_PTR >
    static Apply(const T& x)
    { return BOUND_FAST_MEM_FN< T_MEM_PTR, MEM_PTR >(x); }
};

template< class T_MEM_PTR >
inline MAKER< T_MEM_PTR >
Deduce(T_MEM_PTR)
{ return MAKER< T_MEM_PTR >(); }

} // namespace Detail_BOUND_FAST_MEM_FN

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_BOUND_FAST_MEM_FN_HPP
