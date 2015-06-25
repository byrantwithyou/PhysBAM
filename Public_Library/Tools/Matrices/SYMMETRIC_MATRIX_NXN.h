//#####################################################################
// Copyright 2005-2007, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_NXN
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_NXN__
#define __SYMMETRIC_MATRIX_NXN__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <cassert>
#include <iostream>
namespace PhysBAM{
template<class T,class GENERATOR> class RANDOM_NUMBERS;
template<class T> class MATRIX_MXN;
template<class T,int d> class VECTOR;

template<class T>
class SYMMETRIC_MATRIX_NXN
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    int n; // size of the n by n matrix
    int size; // number of elements in the matrix: (n*n+n)/2
    T *x; // pointer to the one dimensional data

    SYMMETRIC_MATRIX_NXN()
        :n(0),size(0),x(0)
    {}

    SYMMETRIC_MATRIX_NXN(const int n_input);
    SYMMETRIC_MATRIX_NXN(const SYMMETRIC_MATRIX_NXN<T>& matrix_input);

    ~SYMMETRIC_MATRIX_NXN();

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert((unsigned)i<(unsigned)n && (unsigned)j<=(unsigned)i);return x[((2*n-j-1)*j>>1)+i];}

    const T& Element_Lower(int i,int j) const
    {assert((unsigned)i<(unsigned)n && (unsigned)j<=(unsigned)i);return x[((2*n-j-1)*j>>1)+i];}

    template<class RW> void Read(std::istream& input)
    {delete[] x;
    Read_Binary<RW>(input,n);
    assert(n>=0);
    size=(n*n+n)/2;
    x=0;
    if(n>0){x=new T[size];Read_Binary_Array<RW>(input,x,size);}}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,n);Write_Binary_Array<RW>(output,x,size);}

//#####################################################################
    static SYMMETRIC_MATRIX_NXN<T> Outer_Product(const ARRAY<T>& u);
    SYMMETRIC_MATRIX_NXN<T> Sqr() const;
    void Givens_Conjugate(const int i,const int j,const T c,const T s);
    void Jacobi_Solve_Eigenproblem(ARRAY<VECTOR<int,2> >& givens_pairs,ARRAY<VECTOR<T,2> >& givens_coefficients,const T tolerance=(T)1e-5,
        const int max_iterations=1000000);
    template<class GENERATOR>
    void Maximum_Eigenvalue_Eigenvector_Pair(T& max_eigenvalue,ARRAY<T>& max_eigenvector,RANDOM_NUMBERS<T,GENERATOR>* random_numbers=0,const T tolerance=(T)1e-5,
        const T randomization_decay_factor=(T)0.9,const int max_iterations=1000000);
    void In_Place_Cholesky_Factorization(MATRIX_MXN<T>& L);
    SYMMETRIC_MATRIX_NXN<T>& operator=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator+=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator-=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator*=(const T a);
    SYMMETRIC_MATRIX_NXN<T> operator+(const SYMMETRIC_MATRIX_NXN<T>& A) const;
    SYMMETRIC_MATRIX_NXN<T> operator-(const SYMMETRIC_MATRIX_NXN<T>& A) const;
    SYMMETRIC_MATRIX_NXN<T> operator*(const T a) const;
    ARRAY<T> operator*(const ARRAY<T>& y) const;
    void Set_Identity_Matrix();
    void Set_Zero_Matrix();
    static SYMMETRIC_MATRIX_NXN<T> Identity_Matrix(const int n);
    T Trace() const;
//#####################################################################
};
template<class T>
inline SYMMETRIC_MATRIX_NXN<T> operator*(const T a,const SYMMETRIC_MATRIX_NXN<T>& A)
{return A*a;}

template<class T> std::ostream& operator<<(std::ostream& output_stream,const SYMMETRIC_MATRIX_NXN<T>& A);
}
#endif
