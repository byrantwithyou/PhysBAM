//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_GRAD_MULTILINEAR_BASIS_DOT_GRAD_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_GRAD_MULTILINEAR_BASIS_DOT_GRAD_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Math/Eval_Grad_Multilinear_Dot_Grad_Multilinear_From_Monomial_Evals.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

namespace PhysBAM
{

template< class T, int D, class T_EVAL_OF_MONOMIAL >
inline VECTOR< VECTOR< T, (1 << D) >, (1 << D) >
Eval_All_Grad_Multilinear_Basis_Dot_Grad_Multilinear_Basis_From_Monomial_Evals(
    const VECTOR<T,D> dx,
    T_EVAL_OF_MONOMIAL eval_of_monomial)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    const T dv = (1 << D) * dx.Product();
    const T dv2 = dv * dv;
    VECTOR< VECTOR< T, (1 << D) >, (1 << D) > grad_basis_dot_grad_basis_evals;
    for(int i = 1; i <= 1 << D; ++i) {
        const MULTI_INDEX_TYPE multi_index_i = 1 - STATIC_MULTI_INDEX_CUBE<D,0,1>::Multi_Index(i);
        const VECTOR<T,D> xi = (2 * multi_index_i - 1) * dx;
        const int sign_i = 1 - 2 * (multi_index_i.Sum() & 1);
        for(int j = i; j <= 1 << D; ++j) {
            const MULTI_INDEX_TYPE multi_index_j = 1 - STATIC_MULTI_INDEX_CUBE<D,0,1>::Multi_Index(j);
            const VECTOR<T,D> xj = (2 * multi_index_j - 1) * dx;
            const int sign_j = 1 - 2 * (multi_index_j.Sum() & 1);
            grad_basis_dot_grad_basis_evals[i][j] =
                (sign_i * sign_j)
              * Eval_Grad_Multilinear_Dot_Grad_Multilinear_From_Monomial_Evals(xi, xj, eval_of_monomial)
              / dv2;
            grad_basis_dot_grad_basis_evals[j][i] = grad_basis_dot_grad_basis_evals[i][j];
        }
    }
    return grad_basis_dot_grad_basis_evals;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MATH_EVAL_ALL_GRAD_MULTILINEAR_BASIS_DOT_GRAD_MULTILINEAR_BASIS_FROM_MONOMIAL_EVALS_HPP
