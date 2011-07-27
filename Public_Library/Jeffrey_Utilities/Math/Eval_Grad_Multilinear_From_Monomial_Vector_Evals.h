//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_GRAD_MULTILINEAR_FROM_MONOMIAL_VECTOR_EVALS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_GRAD_MULTILINEAR_FROM_MONOMIAL_VECTOR_EVALS_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <boost/mpl/int.hpp>
#include <boost/mpl/vector/vector10_c.hpp>
#include <boost/preprocessor/tuple/rem.hpp>

namespace PhysBAM
{

#define _eval( D, N, d ) \
    eval_of_monomial_vector( \
        boost::mpl::vector ## D ## _c< \
            int, BOOST_PP_TUPLE_REM_CTOR( D, N ) \
        >(), \
        boost::mpl::int_< d >() \
    )

template< class T, class T_EVAL_OF_MONOMIAL_VECTOR >
inline T
Evaluate_Grad_Multilinear_From_Monomial_Vector_Evaluations(
    const VECTOR<T,1>& /*x0*/,
    T_EVAL_OF_MONOMIAL_VECTOR eval_of_monomial_vector)
{ return _eval(1,(0),1); }

template< class T, class T_EVAL_OF_MONOMIAL_VECTOR >
inline T
Evaluate_Grad_Multilinear_From_Monomial_Vector_Evaluations(
    const VECTOR<T,2> x0,
    T_EVAL_OF_MONOMIAL_VECTOR eval_of_monomial_vector)
{
    return _eval(2,(0,1),1) - _eval(2,(0,0),1) * x0[2]
         + _eval(2,(1,0),2) - _eval(2,(0,0),2) * x0[1];
}

template< class T, class T_EVAL_OF_MONOMIAL_VECTOR >
inline T
Evaluate_Grad_Multilinear_From_Monomial_Vector_Evaluations(
    const VECTOR<T,3> x0,
    T_EVAL_OF_MONOMIAL_VECTOR eval_of_monomial_vector)
{
    return _eval(3,(0,1,1),1) - _eval(3,(0,1,0),1) * x0[3] - _eval(3,(0,0,1),1) * x0[2] + _eval(3,(0,0,0),1) * x0[2] * x0[3]
         + _eval(3,(1,0,1),2) - _eval(3,(1,0,0),2) * x0[3] - _eval(3,(0,0,1),2) * x0[1] + _eval(3,(0,0,0),2) * x0[1] * x0[3]
         + _eval(3,(1,1,0),3) - _eval(3,(1,0,0),3) * x0[2] - _eval(3,(0,1,0),3) * x0[1] + _eval(3,(0,0,0),3) * x0[1] * x0[2];
}

#undef _eval

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_GRAD_MULTILINEAR_FROM_MONOMIAL_VECTOR_EVALS_HPP
