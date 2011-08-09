//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cmath>

#include <boost/math/special_functions/pow.hpp>
#include <boost/throw_exception.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/EXAMPLE_PARAMS.h"

#include "Main_Examples.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< int D, int ID >
struct U_EXAMPLE
{ static const bool is_specialized = false; };

template< int D, int ID >
struct BETA_EXAMPLE
{ static const bool is_specialized = false; };

template< int D, int U_ID, int BETA_ID >
struct F_EXAMPLE
{
    static const bool is_specialized = false;
    template< class T >
    struct _
    {
        static T Apply(const VECTOR<T,D>& /*x_*/)
        {
            BOOST_THROW_EXCEPTION( U_BETA_ERROR() );
            return 0; // might be required by the compiler
        }
    };
};

#define USING_MATH_DECLARATIONS \
    using std::sqrt; \
    using std::cos; using std::sin; \
    using std::exp; using std::log; \
    using boost::math::pow;

#define DEFINE_XYZ( D ) DEFINE_XYZ_ ## D ;
#define DEFINE_XYZ_1 const T x = x_[1]
#define DEFINE_XYZ_2 DEFINE_XYZ_1, y = x_[2]
#define DEFINE_XYZ_3 DEFINE_XYZ_2, z = x_[3]

#define DISABLE_UNUSED_XYZ_WARNING( D ) DISABLE_UNUSED_XYZ_WARNING_ ## D
#define DISABLE_UNUSED_XYZ_WARNING_1 static_cast<void>(x);
#define DISABLE_UNUSED_XYZ_WARNING_2 DISABLE_UNUSED_XYZ_WARNING_1 static_cast<void>(y);
#define DISABLE_UNUSED_XYZ_WARNING_3 DISABLE_UNUSED_XYZ_WARNING_2 static_cast<void>(z);

#define XYZ_STR( D ) XYZ_STR_ ## D
#define XYZ_STR_1 "x"
#define XYZ_STR_2 "x,y"
#define XYZ_STR_3 "x,y,z"

#define APPLY_PREAMBLE( D ) \
    USING_MATH_DECLARATIONS \
    DEFINE_XYZ( D ) \
    DISABLE_UNUSED_XYZ_WARNING( D )

#define DEFINE_U_EXAMPLE( D, ID, Expr, Grad_Expr ) \
template<> struct U_EXAMPLE<D,ID> \
{ \
    static const bool is_specialized = true; \
    template< class T > \
    struct _ \
    { \
        static T Apply(const VECTOR<T,D>& x_) \
        { APPLY_PREAMBLE( D ) return Expr ; } \
        static VECTOR<T,D> Apply_Grad(const VECTOR<T,D>& x_) \
        { APPLY_PREAMBLE( D ) return VECTOR<T,D> Grad_Expr ; } \
    }; \
    static const char str[]; \
}; \
const char U_EXAMPLE<D,ID>::str[] = "(" XYZ_STR( D ) ") -> " #Expr ;

#define DEFINE_BETA_EXAMPLE( D, ID, Expr ) \
template<> struct BETA_EXAMPLE<D,ID> \
{ \
    static const bool is_specialized = true; \
    template< class T > \
    struct _ \
    { \
        static T Apply(const VECTOR<T,D>& x_) \
        { APPLY_PREAMBLE( D ) return Expr ; } \
    }; \
    static const char str[]; \
}; \
const char BETA_EXAMPLE<D,ID>::str[] = "(" XYZ_STR( D ) ") -> " #Expr ;

#define DEFINE_F_EXAMPLE( D, U_ID, BETA_ID, Expr ) \
template<> struct F_EXAMPLE< D, U_ID, BETA_ID > \
{ \
    static const bool is_specialized = true; \
    template< class T > \
    struct _ \
    { \
        static T Apply(const VECTOR<T,D>& x_) \
        { APPLY_PREAMBLE( D ) return Expr ; } \
    }; \
};

DEFINE_U_EXAMPLE( 2, 0, 1, (0,0) )
DEFINE_U_EXAMPLE( 2, 1, x + 2*y, (1,2) )
DEFINE_U_EXAMPLE( 2, 2, pow<2>(x) + 2*pow<2>(y), (2*x, 4*y) )

DEFINE_BETA_EXAMPLE( 2, 0, 1 )
DEFINE_BETA_EXAMPLE( 2, 1, 1 + 2*pow<2>(x) + pow<2>(y) )

DEFINE_F_EXAMPLE( 2, 0, 0, 0 )
DEFINE_F_EXAMPLE( 2, 1, 0, 0 )
DEFINE_F_EXAMPLE( 2, 2, 0, -6 )
DEFINE_F_EXAMPLE( 2, 0, 1, 0 )
DEFINE_F_EXAMPLE( 2, 1, 1, -4*(x + y) )
DEFINE_F_EXAMPLE( 2, 2, 1, -2*(3 + 10*pow<2>(x) + 7*pow<2>(y)) )

DEFINE_U_EXAMPLE( 3, 0, 1, (0,0,0) )
DEFINE_U_EXAMPLE( 3, 1, x + 2*y + 3*z, (1,2,3) )
DEFINE_U_EXAMPLE( 3, 2, pow<2>(x) + 2*pow<2>(y) + 3*pow<2>(z), (2*x, 4*y, 6*z) )

DEFINE_BETA_EXAMPLE( 3, 0, 1 )
DEFINE_BETA_EXAMPLE( 3, 1, 1 + 3*pow<2>(x) + 2*pow<2>(y) + pow<2>(z) )

DEFINE_F_EXAMPLE( 3, 0, 0, 0 )
DEFINE_F_EXAMPLE( 3, 1, 0, 0 )
DEFINE_F_EXAMPLE( 3, 2, 0, -12 )
DEFINE_F_EXAMPLE( 3, 0, 1, 0 )
DEFINE_F_EXAMPLE( 3, 1, 1, -2*(3*x + 4*y + 3*z) )
DEFINE_F_EXAMPLE( 3, 2, 1, -4*(3 + 12*pow<2>(x) + 10*pow<2>(y) + 6*pow<2>(z)) )

#undef DEFINE_U_EXAMPLE
#undef DEFINE_BETA_EXAMPLE
#undef DEFINE_F_EXAMPLE

#undef APPLY_PREAMBLE

#undef XYZ_STR
#undef XYZ_STR_1
#undef XYZ_STR_2
#undef XYZ_STR_3

#undef DISABLE_UNUSED_XYZ_WARNING
#undef DISABLE_UNUSED_XYZ_WARNING_1
#undef DISABLE_UNUSED_XYZ_WARNING_2
#undef DISABLE_UNUSED_XYZ_WARNING_3

#undef DEFINE_XYZ
#undef DEFINE_XYZ_1
#undef DEFINE_XYZ_2
#undef DEFINE_XYZ_3

#undef USING_MATH_DECLARATIONS

namespace
{

template< int D, int ID = 0, bool = U_EXAMPLE<D,ID>::is_specialized >
struct VISIT_U_EXAMPLES_STR_ITERATION;
template< int D, int ID = 0, bool = BETA_EXAMPLE<D,ID>::is_specialized >
struct VISIT_BETA_EXAMPLES_STR_ITERATION;

template< int D, int ID >
struct VISIT_U_EXAMPLES_STR_ITERATION< D, ID, false >
{ static void Apply(const boost::function< void ( unsigned int, const char* ) >& /*visitor*/) { } };

template< int D, int ID >
struct VISIT_U_EXAMPLES_STR_ITERATION< D, ID, true >
{
    static void Apply(const boost::function< void ( unsigned int, const char* ) >& visitor)
    {
        visitor(ID, U_EXAMPLE<D,ID>::str);
        VISIT_U_EXAMPLES_STR_ITERATION< D, ID + 1 >::Apply(visitor);
    }
};

template< int D, int ID >
struct VISIT_BETA_EXAMPLES_STR_ITERATION< D, ID, false >
{ static void Apply(const boost::function< void ( unsigned int, const char* ) >& /*visitor*/) { } };

template< int D, int ID >
struct VISIT_BETA_EXAMPLES_STR_ITERATION< D, ID, true >
{
    static void Apply(const boost::function< void ( unsigned int, const char* ) >& visitor)
    {
        visitor(ID, BETA_EXAMPLE<D,ID>::str);
        VISIT_BETA_EXAMPLES_STR_ITERATION< D, ID + 1 >::Apply(visitor);
    }
};

} // namespace

template< int D >
void Visit_U_Examples_Str(const boost::function< void ( unsigned int, const char* ) >& visitor)
{ VISIT_U_EXAMPLES_STR_ITERATION<D>::Apply(visitor); }

template< int D >
void Visit_Beta_Examples_Str(const boost::function< void ( unsigned int, const char* ) >& visitor)
{ VISIT_BETA_EXAMPLES_STR_ITERATION<D>::Apply(visitor); }

#define EXPLICIT_INSTANTIATION( D ) \
template void Visit_U_Examples_Str<D>(const boost::function< void ( unsigned int, const char* ) >& visitor); \
template void Visit_Beta_Examples_Str<D>(const boost::function< void ( unsigned int, const char* ) >& visitor);
EXPLICIT_INSTANTIATION( 2 )
EXPLICIT_INSTANTIATION( 3 )
#undef EXPLICIT_INSTANTIATION

namespace
{

template< int D, int U_ID = 0, bool = U_EXAMPLE< D, U_ID >::is_specialized >
struct GET_MAIN_EXAMPLE_U_ITERATION;
template< int D, int U_ID, int BETA_ID = 0, bool = BETA_EXAMPLE< D, BETA_ID >::is_specialized >
struct GET_MAIN_EXAMPLE_BETA_ITERATION;

template< int D, int U_ID >
struct GET_MAIN_EXAMPLE_U_ITERATION< D, U_ID, false >
{
    template< class T >
    static void Apply(
        const unsigned int /*u_id*/, const unsigned int /*beta_id*/,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& /*u*/,
        typename EXAMPLE_PARAMS<T,D>::VECTOR_FUNCTION_TYPE& /*grad_u*/,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& /*beta*/,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& /*f*/)
    { }
};

template< int D, int U_ID >
struct GET_MAIN_EXAMPLE_U_ITERATION< D, U_ID, true >
{
    template< class T >
    static void Apply(
        const unsigned int u_id, const unsigned int beta_id,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& u,
        typename EXAMPLE_PARAMS<T,D>::VECTOR_FUNCTION_TYPE& grad_u,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& beta,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& f)
    {
        if(u_id == U_ID) {
            u      = &U_EXAMPLE< D, U_ID >::template _<T>::Apply;
            grad_u = &U_EXAMPLE< D, U_ID >::template _<T>::Apply_Grad;
            GET_MAIN_EXAMPLE_BETA_ITERATION< D, U_ID >::template Apply<T>(beta_id, beta, f);
        }
        else
            GET_MAIN_EXAMPLE_U_ITERATION< D, U_ID + 1 >::template Apply<T>(u_id, beta_id, u, grad_u, beta, f);
    }
};

template< int D, int U_ID, int BETA_ID >
struct GET_MAIN_EXAMPLE_BETA_ITERATION< D, U_ID, BETA_ID, false >
{
    template< class T >
    static void Apply(
        const unsigned int /*beta_id*/,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& /*beta*/,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& /*f*/)
    { }
};

template< int D, int U_ID, int BETA_ID >
struct GET_MAIN_EXAMPLE_BETA_ITERATION< D, U_ID, BETA_ID, true >
{
    template< class T >
    static void Apply(
        const unsigned int beta_id,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& beta,
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& f)
    {
        if(beta_id == BETA_ID) {
            beta = &BETA_EXAMPLE< D,       BETA_ID >::template _<T>::Apply;
            f    = &   F_EXAMPLE< D, U_ID, BETA_ID >::template _<T>::Apply;
        }
        else
            GET_MAIN_EXAMPLE_BETA_ITERATION< D, U_ID, BETA_ID + 1 >::template Apply<T>(beta_id, beta, f);
    }
};

template< int D, int ID = 0, bool = U_EXAMPLE<D,ID>::is_specialized >
struct MAX_U_ID;
template< int D, int ID = 0, bool = BETA_EXAMPLE<D,ID>::is_specialized >
struct MAX_BETA_ID;

template< int D, int ID >
struct MAX_U_ID< D, ID, false >
{ static const int value = ID; };
template< int D, int ID >
struct MAX_U_ID< D, ID, true >
{ static const int value = MAX_U_ID< D, ID + 1 >::value; };

template< int D, int ID >
struct MAX_BETA_ID< D, ID, false >
{ static const int value = ID; };
template< int D, int ID >
struct MAX_BETA_ID< D, ID, true >
{ static const int value = MAX_BETA_ID< D, ID + 1 >::value; };

} // namespace

template< class T, int D >
void
Get_Main_Example(
    const unsigned int u_id, const unsigned int beta_id,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& u,
    typename EXAMPLE_PARAMS<T,D>::VECTOR_FUNCTION_TYPE& grad_u,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& beta,
    typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& f)
{ GET_MAIN_EXAMPLE_U_ITERATION<D>::template Apply<T>(u_id, beta_id, u, grad_u, beta, f); }

#define EXPLICIT_INSTANTIATION( T, D ) \
template void Get_Main_Example<T,D>( \
    const unsigned int u_id, const unsigned int beta_id, \
    EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& u, \
    EXAMPLE_PARAMS<T,D>::VECTOR_FUNCTION_TYPE& grad_u, \
    EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& beta, \
    EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE& f);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
