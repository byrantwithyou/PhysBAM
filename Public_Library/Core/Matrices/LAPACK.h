//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPACK
//#####################################################################
#ifndef __LAPACK__
#define __LAPACK__

#include <complex>
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Matrices/MATRIX_BASE.h>
#include <cblas.h>
#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include <lapacke.h>
namespace PhysBAM{

// Overloads

inline int xgelsy(int matrix_layout, int m, int n, int nrhs, double* a,
    int lda, double* b, int ldb, int* jpvt, double rcond, int* rank)
{
    return LAPACKE_dgelsy( matrix_layout, m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank );
}

inline int xgelsy(int matrix_layout, int m, int n, int nrhs, float* a,
    int lda, float* b, int ldb, int* jpvt, float rcond, int* rank)
{
    return LAPACKE_sgelsy( matrix_layout, m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank );
}

inline int xgelsy(int matrix_layout, int m, int n, int nrhs, std::complex<double>* a,
    int lda, std::complex<double>* b, int ldb, int* jpvt, double rcond, int* rank)
{
    return LAPACKE_zgelsy( matrix_layout, m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank );
}

inline int xgelsy(int matrix_layout, int m, int n, int nrhs, std::complex<float>* a,
    int lda, std::complex<float>* b, int ldb, int* jpvt, float rcond, int* rank)
{
    return LAPACKE_cgelsy( matrix_layout, m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank );
}

// Uses

template<class T,class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX1>& M,
    const MATRIX_BASE<T,T_MATRIX2>& rhs,MATRIX_BASE<T,T_MATRIX3>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0)
{
    int m=M.Rows(),n=M.Columns(),nr=rhs.Columns(),ldb=max(m,n);
    assert(m==rhs.Rows());
    assert(n==sol.Rows());
    assert(nr==sol.Columns());
    ARRAY<T> A(m*n);
    ARRAY<T> V(nr*ldb);
    ARRAY<int> jpvt(n);
    for(int j=0;j<n;j++)
        for(int i=0;i<m;i++)
            A(j*m+i)=M(i,j);
    for(int j=0;j<nr;j++)
        for(int i=0;i<m;i++)
            V(j*ldb+i)=rhs(i,j);

    int rnk = 0;
    int ret = xgelsy( LAPACK_COL_MAJOR, m, n, nr, A.Get_Array_Pointer(), m,
        V.Get_Array_Pointer(), ldb, jpvt.Get_Array_Pointer(), rcond, &rnk );
    PHYSBAM_ASSERT(!ret);
    if(rank) *rank=rnk;

    for(int j=0;j<nr;j++)
        for(int i=0;i<n;i++)
            sol(i,j)=V(j*ldb+i);
}

template<class T,class T_MATRIX,class T_VECTOR1,class T_VECTOR2>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX>& M,
    const ARRAY_BASE<T,T_VECTOR1>& rhs,ARRAY_BASE<T,T_VECTOR2>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0)
{
    int m=M.Rows(),n=M.Columns(),ldb=max(m,n);
    assert(m==rhs.Size());
    assert(n==sol.Size());
    ARRAY<T> A(m*n);
    ARRAY<T> V(rhs);
    ARRAY<int> jpvt(n);
    for(int j=0;j<n;j++)
        for(int i=0;i<m;i++)
            A(j*m+i)=M(i,j);

    int rnk = 0;
    int ret = xgelsy( LAPACK_COL_MAJOR, m, n, 1, A.Get_Array_Pointer(), m,
        V.Get_Array_Pointer(), ldb, jpvt.Get_Array_Pointer(), rcond, &rnk );
    PHYSBAM_ASSERT(!ret);
    if(rank) *rank=rnk;

    for(int i=0;i<V.m;i++) sol(i)=V(i);
}

}
#endif
