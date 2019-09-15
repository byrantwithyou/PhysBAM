//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPACK
//#####################################################################
#ifndef __LAPACK__
#define __LAPACK__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Matrices/Make_View.h>
#include <Core/Matrices/MATRIX_BASE.h>
#include <Core/Matrices/MATRIX_VIEW.h>
#include <complex>
#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include <cblas.h>
#include <lapacke.h>
namespace PhysBAM{

//#####################################################################
// Forward Declarations
//#####################################################################

// Least_Squares_Solve_Raw and Least_Squares_Solve work with non-square an singular M.

// Destroys M
template<class T,class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
void Least_Squares_Solve_Raw(MATRIX_BASE<T,T_MATRIX1>& M,
    const MATRIX_BASE<T,T_MATRIX2>& rhs,
    MATRIX_BASE<T,T_MATRIX3>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0);

// Preserves M
template<class T,class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX1>& M,
    const MATRIX_BASE<T,T_MATRIX2>& rhs,
    MATRIX_BASE<T,T_MATRIX3>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0);

// Destroys M
template<class T,class T_MATRIX,class T_ARRAY1,class T_ARRAY2>
void Least_Squares_Solve_Raw(MATRIX_BASE<T,T_MATRIX>& M,
    const ARRAY_BASE<T,T_ARRAY1>& rhs,
    ARRAY_BASE<T,T_ARRAY2>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0);

// Preserves M
template<class T,class T_MATRIX,class T_ARRAY1,class T_ARRAY2>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX>& M,
    const ARRAY_BASE<T,T_ARRAY1>& rhs,
    ARRAY_BASE<T,T_ARRAY2>& sol,
    T rcond=(std::numeric_limits<T>::epsilon()), int * rank=0);

//#####################################################################
// LAPACK Templatizations Wrappers
//#####################################################################

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

//#####################################################################
// Definitions
//#####################################################################

// Destroys M
template<class T,class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
void Least_Squares_Solve_Raw(MATRIX_BASE<T,T_MATRIX1>& M,
    const MATRIX_BASE<T,T_MATRIX2>& rhs,
    MATRIX_BASE<T,T_MATRIX3>& sol,T rcond, int * rank)
{
    assert(M.Rows()==rhs.Rows());
    assert(M.Columns()==sol.Rows());
    assert(rhs.Columns()==sol.Columns());

    MATRIX_MXN<T> TS;
    MATRIX_VIEW<T> VM(M.Derived()),VS(sol.Derived());

    bool b=false;
    if(sol.Rows()<rhs.Rows())
    {
        TS.Resize(rhs.Rows(),rhs.Columns());
        VS.Set(TS);
        b=true;
    }
    VS.m=rhs.Rows();
    VS=rhs;
    VS.m=sol.Rows();

    ARRAY<int> jpvt(VM.n);
    int rnk = 0;
    int ret = xgelsy( LAPACK_COL_MAJOR, VM.m, VM.n, rhs.Columns(), VM.x, VM.s,
        VS.x, VS.s, jpvt.Get_Array_Pointer(), rcond, &rnk );
    PHYSBAM_ASSERT(!ret);
    if(rank) *rank=rnk;
    if(b) sol.Derived()=VS;
}

// Preserves M
template<class T,class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX1>& M,
    const MATRIX_BASE<T,T_MATRIX2>& rhs,
    MATRIX_BASE<T,T_MATRIX3>& sol,T rcond, int * rank)
{
    MATRIX_MXN<T> MM(M);
    Least_Squares_Solve_Raw(MM,rhs,sol,rcond,rank);
}

// Destroys M
template<class T,class T_MATRIX,class T_ARRAY1,class T_ARRAY2>
void Least_Squares_Solve_Raw(MATRIX_BASE<T,T_MATRIX>& M,
    const ARRAY_BASE<T,T_ARRAY1>& rhs,
    ARRAY_BASE<T,T_ARRAY2>& sol,T rcond, int * rank)
{
    assert(M.Rows()==rhs.Size());
    assert(M.Columns()==sol.Size());

    ARRAY<T> TS;
    MATRIX_VIEW<T> VM(M.Derived());
    ARRAY_VIEW<T> VS;

    bool b=true;
    if(sol.Size()>=rhs.Size()) b=Make_View(VS,sol,TS);
    else
    {
        TS.Resize(rhs.Size());
        VS.Set(TS);
    }
    VS.m=rhs.Size();
    VS=rhs;
    VS.m=sol.Size();

    int ldb=max(VM.m,VM.n);
    ARRAY<int> jpvt(VM.n);
    int rnk = 0;
    int ret = xgelsy( LAPACK_COL_MAJOR, VM.m, VM.n, 1, VM.x, VM.s,
        VS.base_pointer, ldb, jpvt.Get_Array_Pointer(), rcond, &rnk );
    PHYSBAM_ASSERT(!ret);
    if(rank) *rank=rnk;
    if(b) sol.Derived()=VS;
}

// Preserves M
template<class T,class T_MATRIX,class T_ARRAY1,class T_ARRAY2>
void Least_Squares_Solve(const MATRIX_BASE<T,T_MATRIX>& M,
    const ARRAY_BASE<T,T_ARRAY1>& rhs,
    ARRAY_BASE<T,T_ARRAY2>& sol,T rcond, int * rank)
{
    MATRIX_MXN<T> MM(M);
    Least_Squares_Solve_Raw(MM,rhs,sol,rcond,rank);
}

}
#endif
