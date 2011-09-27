//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Math/Eval_Multilinear_From_Monomial_Evals.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

namespace PhysBAM
{

template< class T, int D, class T_EVAL_OF_MONOMIAL >
inline VECTOR< T, (1 << D) >
Eval_All_Multilinear_Basis_From_Monomial_Evals(
    const VECTOR<T,D> dx,
    T_EVAL_OF_MONOMIAL eval_of_monomial)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    const T dv = (1 << D) * dx.Product();
    VECTOR< T, (1 << D) > basis_evals; // init'ed to 0
    for(int i = 1; i <= 1 << D; ++i) {
        const MULTI_INDEX_TYPE multi_index = 1 - STATIC_MULTI_INDEX_CUBE<D,0,1>::Multi_Index(i);
        const VECTOR<T,D> x0 = (2 * multi_index - 1) * dx;
        const int sign = 1 - 2 * (multi_index.Sum() & 1);
        basis_evals[i] = sign * Eval_Multilinear_From_Monomial_Evals(x0, eval_of_monomial) / dv;
    }
    return basis_evals;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP
