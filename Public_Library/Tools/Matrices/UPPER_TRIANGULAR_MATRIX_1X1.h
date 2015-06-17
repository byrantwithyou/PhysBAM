//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_1X1
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_1X1__
#define __UPPER_TRIANGULAR_MATRIX_1X1__

#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_0X0.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <cmath>
#include <ostream>
namespace PhysBAM{

using ::std::sqrt;

template<class T> struct is_scalar_BLOCK<UPPER_TRIANGULAR_MATRIX<T,1> >:public is_scalar_BLOCK<T>{};
template<class T> struct is_scalar_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,1> >:public is_scalar_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,1>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,1>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=1};

    T x00;

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(1))
        :x00(T())
    {
        STATIC_ASSERT(sizeof(UPPER_TRIANGULAR_MATRIX)==sizeof(T));assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(1));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,1>& matrix_input)
        :x00(matrix_input.x00)
    {}

    UPPER_TRIANGULAR_MATRIX(const T x11_input)
        :x00(x11_input)
    {}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 1;}

    T& operator()(const int i,const int j)
    {assert(i==0 && j==0);return x00;}

    const T& operator()(const int i,const int j) const
    {assert(i==0 && j==0);return x00;}

    bool Valid_Index(const int i,const int j) const
    {return i==0 && j==0;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return x00==A.x00;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return !(*this==A);}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return UPPER_TRIANGULAR_MATRIX(-x00);}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00+=A.x00;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {x00+=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00-=A.x00;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {x00-=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00);}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,1>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00);}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00+a);}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00);}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,1>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00);}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00-a);}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this=*this*A;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {x00*=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {assert(a!=0);x00/=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(a*x00);}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {assert(a!=0);return UPPER_TRIANGULAR_MATRIX(x00/a);}

    VECTOR<T,1> operator*(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x00*v.x);}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x00);}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,1>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x.x);}

    template<class T_MATRIX>
    T_MATRIX operator*(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==1);return x00*A.Derived();}

    template<class T_MATRIX>
    typename TRANSPOSE<T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Columns()==1);typename TRANSPOSE<T_MATRIX>::TYPE M((INITIAL_SIZE)1,(INITIAL_SIZE)A.Rows());
    for(int j=0;j<A.Rows();j++) M(0,j)+=x00*A(j,0);return M;}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==1);return x00*A.Derived();}

    MATRIX<T,1> Transpose_Times(const UPPER_TRIANGULAR_MATRIX& A) const
    {return MATRIX<T,1>(x00*A.x00);}

    MATRIX<T,1> Transpose_Times(const SYMMETRIC_MATRIX<T,1>& A) const
    {return MATRIX<T,1>(x00*A.x00);}

    UPPER_TRIANGULAR_MATRIX Transpose_Times(const DIAGONAL_MATRIX<T,1>& A) const
    {return A.x.x**this;}

    MATRIX<T,1> Times_Transpose(const UPPER_TRIANGULAR_MATRIX& A) const
    {return MATRIX<T,1>(x00*A.x00);}

    MATRIX<T,1> Times_Transpose(const SYMMETRIC_MATRIX<T,1>& A) const
    {return *this*A;}

    UPPER_TRIANGULAR_MATRIX Times_Transpose(const DIAGONAL_MATRIX<T,1>& A) const
    {return *this*A;}

    SYMMETRIC_MATRIX<T,1> Outer_Product_Matrix() const
    {return SYMMETRIC_MATRIX<T,1>(x00*x00);}

    T Determinant() const
    {return x00;}

    T Trace() const
    {return x00;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {assert(x00!=0);return UPPER_TRIANGULAR_MATRIX(1/x00);}

    VECTOR<T,1> Inverse_Times(const VECTOR<T,1>& b) const
    {return Inverse()*b;}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return UPPER_TRIANGULAR_MATRIX(1);}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX(1);}

    T Max_Abs() const
    {return abs(x00);}

    T Frobenius_Norm() const
    {return abs(x00);}

    T Frobenius_Norm_Squared() const
    {return sqr(x00);}

    MATRIX<T,1> Transposed() const
    {return MATRIX<T,1>(x00);}

//#####################################################################
};
// global functions

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,1> operator*(const DIAGONAL_MATRIX<T,1>& A,const UPPER_TRIANGULAR_MATRIX<T,1>& B) // 3 mults
{return UPPER_TRIANGULAR_MATRIX<T,1>(A.x.x*B.x00);}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,1>& A)
{return output_stream<<"["<<A.x00<<"]";}
//#####################################################################
}
#endif

