//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_MULTILINEAR_FROM_MONOMIAL_EVALS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_MULTILINEAR_FROM_MONOMIAL_EVALS_HPP

#include <boost/mpl/vector/vector10_c.hpp>
#include <boost/preprocessor/tuple/rem.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

#define _eval( D, N ) \
    eval_of_monomial( \
        boost::mpl::vector ## D ## _c< \
            int, BOOST_PP_TUPLE_REM_CTOR( D, N ) \
        >() \
    )

template< class T, class T_EVAL_OF_MONOMIAL >
inline T
Eval_Multilinear_From_Monomial_Evals(
    const VECTOR<T,1> x0,
    T_EVAL_OF_MONOMIAL eval_of_monomial)
{ return _eval(1,(1)) - _eval(1,(0)) * x0[1]; }

template< class T, class T_EVAL_OF_MONOMIAL >
inline T
Eval_Multilinear_From_Monomial_Evals(
    const VECTOR<T,2> x0,
    T_EVAL_OF_MONOMIAL eval_of_monomial)
{ return _eval(2,(1,1)) - _eval(2,(1,0)) * x0[2] - _eval(2,(0,1)) * x0[1] + _eval(2,(0,0)) * x0[1] * x0[2]; }

template< class T, class T_EVAL_OF_MONOMIAL >
inline T
Eval_Multilinear_From_Monomial_Evals(
    const VECTOR<T,3> x0,
    T_EVAL_OF_MONOMIAL eval_of_monomial)
{
    return _eval(3,(1,1,1))
         - _eval(3,(1,1,0)) * x0[3] - _eval(3,(1,0,1)) * x0[2] - _eval(3,(0,1,1)) * x0[1]
         + _eval(3,(1,0,0)) * x0[2] * x0[3] + _eval(3,(0,1,0)) * x0[1] * x0[3] + _eval(3,(0,0,1)) * x0[1] * x0[2]
         - _eval(3,(0,0,0)) * x0[1] * x0[2] * x0[3];
}

#undef _eval

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_MULTILINEAR_FROM_MONOMIAL_EVALS_HPP
