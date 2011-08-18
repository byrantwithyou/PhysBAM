//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef BOOST_PP_IS_ITERATING

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_COMPOSE_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_COMPOSE_FUNCTION_HPP

#include <boost/function_types/function_type.hpp>
#include <boost/function_types/parameter_types.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/single_view.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/facilities/intercept.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/PROPAGATE_CONST.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

#ifndef PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY
#define PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY 10
#endif // #ifndef PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY

namespace PhysBAM
{

template< class T_F, class T_G >
struct COMPOSE_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        COMPOSE_FUNCTION,
        (( typename T_F const, f ))
        (( typename T_G const, g ))
    )
public:

    template< class T_SIGNATURE >
    struct result
        : RESULT_OF< const T_F (
              typename RESULT_OF<
                  typename boost::function_types::function_type<
                      boost::mpl::joint_view<
                          boost::mpl::single_view< const T_G >,
                          typename boost::function_types::parameter_types< T_SIGNATURE >::type
                      >
                  >::type
              >::type
          ) >
    { };

    template< class T1 >
    typename result< const COMPOSE_FUNCTION ( T1& ) >::type
    operator()(T1& x1) const
    { return f(g(x1)); }

    template< class T1 >
    typename result< const COMPOSE_FUNCTION ( const T1& ) >::type
    operator()(const T1& x1) const
    { return f(g(x1)); }

    template< class T1, class T2 >
    typename result< const COMPOSE_FUNCTION ( T1&, T2& ) >::type
    operator()(T1& x1, T2& x2) const
    { return f(g(x1, x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE_FUNCTION ( T1&, const T2& ) >::type
    operator()(T1& x1, const T2& x2) const
    { return f(g(x1, x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE_FUNCTION ( const T1&, T2& ) >::type
    operator()(const T1& x1, T2& x2) const
    { return f(g(x1, x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE_FUNCTION ( const T1&, const T2& ) >::type
    operator()(const T1& x1, const T2& x2) const
    { return f(g(x1, x2)); }
};

template< class T_F, class T_G1, class T_G2 >
struct COMPOSE2_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        COMPOSE2_FUNCTION,
        (( typename T_F const, f ))
        (( typename T_G1 const, g1 ))
        (( typename T_G2 const, g2 ))
    )
public:

    template< class T_SIGNATURE >
    struct result
        : RESULT_OF< const T_F (
              typename RESULT_OF<
                  typename boost::function_types::function_type<
                      boost::mpl::joint_view<
                          boost::mpl::single_view< const T_G1 >,
                          typename boost::function_types::parameter_types< T_SIGNATURE >::type
                      >
                  >::type
              >::type,
              typename RESULT_OF<
                  typename boost::function_types::function_type<
                      boost::mpl::joint_view<
                          boost::mpl::single_view< const T_G2 >,
                          typename boost::function_types::parameter_types< T_SIGNATURE >::type
                      >
                  >::type
              >::type
          ) >
    { };

    template< class T >
    typename result< const COMPOSE2_FUNCTION ( T& ) >::type
    operator()(T& x) const
    { return f(g1(x), g2(x)); }

    template< class T >
    typename result< const COMPOSE2_FUNCTION ( const T& ) >::type
    operator()(const T& x) const
    { return f(g1(x), g2(x)); }

    template< class T1, class T2 >
    typename result< const COMPOSE2_FUNCTION ( T1&, T2& ) >::type
    operator()(T1& x1, T2& x2) const
    { return f(g1(x1,x2), g2(x1,x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE2_FUNCTION ( T1&, const T2& ) >::type
    operator()(T1& x1, const T2& x2) const
    { return f(g1(x1,x2), g2(x1,x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE2_FUNCTION ( const T1&, T2& ) >::type
    operator()(const T1& x1, T2& x2) const
    { return f(g1(x1,x2), g2(x1,x2)); }

    template< class T1, class T2 >
    typename result< const COMPOSE2_FUNCTION ( const T1&, const T2& ) >::type
    operator()(const T1& x1, const T2& x2) const
    { return f(g1(x1,x2), g2(x1,x2)); }
};

namespace Result_Of
{

template<
    class T_F0, class T_F1,
    BOOST_PP_ENUM_BINARY_PARAMS(
        BOOST_PP_DEC( BOOST_PP_DEC( PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY ) ),
        class T_G, = void BOOST_PP_INTERCEPT
    )
>
struct MAKE_COMPOSE_FUNCTION;

} // namespace Result_Of

#define REPEAT_DATA( z, n, data ) data

#define BOOST_PP_ITERATION_LIMITS ( 2, PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY )
#define BOOST_PP_FILENAME_1 <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include BOOST_PP_ITERATE()

#undef REPEAT_DATA

template< class T_F, class T_G1, class T_G2 >
inline COMPOSE2_FUNCTION< T_F, T_G1, T_G2 >
Make_Compose2_Function(const T_F& f, const T_G1& g1, const T_G2& g2)
{ return COMPOSE2_FUNCTION< T_F, T_G1, T_G2 >(f, g1, g2); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_COMPOSE_FUNCTION_HPP

#else // #ifndef BOOST_PP_IS_ITERATING

#define N BOOST_PP_ITERATION()

namespace Result_Of
{

template< BOOST_PP_ENUM_PARAMS( N, class T_F ) >
struct MAKE_COMPOSE_FUNCTION
#if N < PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY
< BOOST_PP_ENUM_PARAMS( N, T_F ) >
#endif // #if N < PHYSBAM_MAKE_COMPOSE_FUNCTION_MAX_ARITY
{
    typedef BOOST_PP_ENUM_PARAMS( BOOST_PP_DEC( N ), COMPOSE_FUNCTION< T_F ),
            BOOST_PP_CAT( T_F, BOOST_PP_DEC( N ) )
            BOOST_PP_REPEAT( BOOST_PP_DEC( N ), REPEAT_DATA , > ) type;
};

} // namespace Result_Of

template< BOOST_PP_ENUM_PARAMS( N, class T_F ) >
inline typename Result_Of::MAKE_COMPOSE_FUNCTION< BOOST_PP_ENUM_PARAMS( N, T_F ) >::type
Make_Compose_Function( BOOST_PP_ENUM_BINARY_PARAMS( N, const T_F, & f ) )
{
#if N == 2
    return COMPOSE_FUNCTION< T_F0, T_F1 >(f0, f1);
#else // #if N == 2
    return Make_Compose_Function(
        BOOST_PP_ENUM_PARAMS( BOOST_PP_DEC( BOOST_PP_DEC( N ) ), f ),
        Make_Compose_Function(
            BOOST_PP_CAT( f, BOOST_PP_DEC( BOOST_PP_DEC( N ) ) ),
            BOOST_PP_CAT( f, BOOST_PP_DEC( N ) )
        )
    );
#endif // #if N == 2
}

#endif // #ifndef BOOST_PP_IS_ITERATING
