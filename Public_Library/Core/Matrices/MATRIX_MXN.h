//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_MXN
//#####################################################################
#ifndef __MATRIX_MXN__
#define __MATRIX_MXN__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX_BASE.h>
namespace PhysBAM{

template<class T>
class MATRIX_MXN:public MATRIX_BASE<T,MATRIX_MXN<T> >
{
public:
    typedef MATRIX_BASE<T,MATRIX_MXN<T> > BASE;
    using BASE::operator+;using BASE::operator+=;using BASE::operator-;using BASE::operator-=;using BASE::operator*;using BASE::operator*=;
    using BASE::Times_Transpose;using BASE::Transpose_Times;
    typedef T SCALAR;
    int m,n; // size of the m by n matrix
    ARRAY<T> x; // one dimensional data

    MATRIX_MXN()
        :m(0),n(0)
    {}

    MATRIX_MXN(const int m_input,const int n_input)
        :m(m_input),n(n_input),x(m*n)
    {
    }

    explicit MATRIX_MXN(const int n_input)
        :m(n_input),n(n_input),x(m*n)
    {
    }

    MATRIX_MXN(const INITIAL_SIZE m_input,const INITIAL_SIZE n_input)
        :m(Value(m_input)),n(Value(n_input)),x(m*n)
    {
    }

    explicit MATRIX_MXN(const INITIAL_SIZE n_input)
        :m(Value(n_input)),n(Value(n_input)),x(m*n)
    {
    }

    MATRIX_MXN(const MATRIX_MXN<T>& A)
        :m(A.m),n(A.n),x(A.x)
    {
    }

    template<class T2>
    explicit MATRIX_MXN(const MATRIX_MXN<T2>& A)
        :m(A.m),n(A.n),x(A.x)
    {
    }

    template<class T_MATRIX>
    explicit MATRIX_MXN(const MATRIX_BASE<T,T_MATRIX>& A)
        :m(A.Rows()),n(A.Columns()),x(m*n)
    {for(int i=0;i<m;i++) for(int j=0;j<n;j++) (*this)(i,j)=A(i,j);}

    template<int d>
    explicit MATRIX_MXN(const DIAGONAL_MATRIX<T,d>& A)
        :m(d),n(d),x(m*n)
    {for(int i=0;i<d;i++) (*this)(i,i)=A(i,i);}

    template<int d>
    explicit MATRIX_MXN(const SYMMETRIC_MATRIX<T,d>& A)
        :m(A.Rows()),n(A.Columns()),x(m*n)
    {for(int i=0;i<m;i++) for(int j=0;j<=i;j++) (*this)(i,j)=(*this)(j,i)=A.Element_Lower(i,j);}

    template<int d>
    explicit MATRIX_MXN(const UPPER_TRIANGULAR_MATRIX<T,d>& A)
        :m(A.Rows()),n(A.Columns()),x(m*n)
    {for(int i=0;i<m;i++) for(int j=i;j<n;j++) (*this)(i,j)=A(i,j);}

    ~MATRIX_MXN()
    {}

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    void Resize(const int m_new,const int n_new)
    {if(m_new==m && n_new==n) return;
    if(n_new==n){x.Resize(m_new*n);m=m_new;return;}
    ARRAY<T> x_new(m_new*n_new);
    int m1=min(m,m_new),n1=min(n,n_new);for(int i=0;i<m1;i++) for(int j=0;j<n1;j++) x_new(i*n_new+j)=(*this)(i,j);
    x.Exchange(x_new);m=m_new;n=n_new;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);return x(i*n+j);}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);return x(i*n+j);}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<(unsigned)m && (unsigned)j<(unsigned)n;}

    bool operator==(const MATRIX_MXN& A) const
    {if(m!=A.m || n!=A.n) return false;
    for(int i=0;i<m*n;i++) if(x(i)!=A.x(i)) return false;
    return true;}

    bool operator!=(const MATRIX_MXN& A) const
    {return !(*this==A);}

    MATRIX_MXN<T>& operator=(const MATRIX_MXN<T>& A)
    {if(&A!=this){x=A.x;m=A.m;n=A.n;}return *this;}

    template<class T_MATRIX>
    MATRIX_MXN<T>& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {if((void*)&A==(void*)this) return *this;
    x.Resize(A.Rows()*A.Columns());
    m=A.Rows();n=A.Columns();
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) (*this)(i,j)=A(i,j);
    return *this;}

    template<int d>
    MATRIX_MXN<T>& operator=(const DIAGONAL_MATRIX<T,d>& A)
    {x.Resize(A.Rows()*A.Columns());m=n=d;
    for(int k=0;k<m*n;k++) x(k)=0;
    for(int i=0;i<d;i++) (*this)(i,i)=A(i,i);
    return *this;}

    template<int d>
    MATRIX_MXN<T>& operator=(const SYMMETRIC_MATRIX<T,d>& A)
    {x.Resize(A.Rows()*A.Columns());m=n=d;
    for(int i=0;i<m;i++) for(int j=0;j<=i;j++) (*this)(i,j)=(*this)(j,i)=A.Element_Lower(i,j);
    return *this;}

    template<int d>
    MATRIX_MXN<T>& operator=(const UPPER_TRIANGULAR_MATRIX<T,d>& A)
    {x.Resize(d*d);m=n=d;
    for(int j=0;j<d;j++) for(int i=0;i<=j;i++) (*this)(i,j)=A(i,j);
    return *this;}

    void Transpose()
    {*this=Transposed();}

    MATRIX_MXN<T> Transposed() const
    {MATRIX_MXN<T> matrix(n,m);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) matrix(j,i)=(*this)(i,j);
    return matrix;}

    MATRIX_MXN<T> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,2>& v) const
    {assert(n==2);return (*this)*MATRIX<T,1,2>::Cross_Product_Matrix(v).Transposed();}

    MATRIX_MXN<T> Cross_Product_Matrix_Times(const VECTOR<T,2>& v) const
    {assert(m==2);return MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this);}

    MATRIX_MXN<T> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,3>& v) const
    {assert(n==3);MATRIX_MXN<T> matrix(m,3);
    for(int i=0;i<m;i++) matrix.Set_Row(i,v.Cross(VECTOR<T,3>((*this)(i,0),(*this)(i,1),(*this)(i,2))));
    return matrix;}

    MATRIX_MXN<T> Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const
    {assert(m==3);MATRIX_MXN<T> matrix(3,n);
    for(int i=0;i<n;i++) matrix.Set_Column(i,v.Cross(VECTOR<T,3>((*this)(0,i),(*this)(1,i),(*this)(2,i))));
    return matrix;}

    MATRIX_MXN<T> Normal_Equations_Matrix() const
    {MATRIX_MXN<T> result(n);
    for(int j=0;j<n;j++)
        for(int i=j;i<n;i++){
            T a=0;
            for(int k=0;k<m;k++) a+=(*this)(k,i)*(*this)(k,j);
            result(i,j)=result(j,i)=a;}
    return result;}

    ARRAY<T> Normal_Equations_Solve(const ARRAY<T>& b) const
    {MATRIX_MXN<T> A_transpose_A(Normal_Equations_Matrix());ARRAY<T> A_transpose_b(Transpose_Times(b));return A_transpose_A.Cholesky_Solve(A_transpose_b);}

    void Append_Row()
    {x.Resize(++m*n);}
//#####################################################################
};
template<class T>
inline MATRIX_MXN<T> operator*(const T a,const MATRIX_MXN<T>& A)
{return A*a;}
}
#endif
