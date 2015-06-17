//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_0X0
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_0X0__
#define __UPPER_TRIANGULAR_MATRIX_0X0__

#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <cmath>
#include <ostream>
namespace PhysBAM{

using ::std::sqrt;

template<class T> struct is_scalar_BLOCK<UPPER_TRIANGULAR_MATRIX<T,0> >:public is_scalar_BLOCK<T>{};
template<class T> struct is_scalar_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,0> >:public is_scalar_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,0>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,0>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=0,n=0};

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(0),INITIAL_SIZE nn=INITIAL_SIZE(0))
    {
        assert(mm==INITIAL_SIZE(0) && nn==INITIAL_SIZE(0));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,0>& matrix_input)
    {}

    int Rows() const
    {return 0;}

    int Columns() const
    {return 0;}

    T& operator()(const int i,const int j)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(const int i,const int j) const
    {PHYSBAM_FATAL_ERROR();}

    bool Valid_Index(const int i,const int j) const
    {return false;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return true;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return false;}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,0>& A) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,0>& A) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {return *this;}

    VECTOR<T,0> operator*(const VECTOR<T,0>& v) const
    {return VECTOR<T,0>();}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const
    {return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,0>& A) const
    {return *this;}

    template<class T_MATRIX>
    T_MATRIX operator*(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==0);return A;}

    template<class T_MATRIX>
    typename TRANSPOSE<T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Columns()==0);return A.Derived().Transposed();}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==0);return A;}

    MATRIX<T,0> Transpose_Times(const UPPER_TRIANGULAR_MATRIX& A) const
    {return MATRIX<T,0>();}

    UPPER_TRIANGULAR_MATRIX Transpose_Times(const DIAGONAL_MATRIX<T,0>& A) const
    {return *this;}

    MATRIX<T,0> Transpose_Times(const SYMMETRIC_MATRIX<T,0>& A) const
    {return MATRIX<T,0>();}

    SYMMETRIC_MATRIX<T,0> Outer_Product_Matrix() const
    {return SYMMETRIC_MATRIX<T,0>();}

    T Determinant() const
    {return 1;}

    T Trace() const
    {return 0;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {return *this;}

    VECTOR<T,0> Inverse_Times(const VECTOR<T,0>& b) const
    {return Inverse()*b;}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return *this;}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX();}

    T Max_Abs() const
    {return 0;}

    T Frobenius_Norm() const
    {return 0;}

    T Frobenius_Norm_Squared() const
    {return 0;}

    MATRIX<T,0> Transposed() const
    {return MATRIX<T,0>();}

//#####################################################################
};
// global functions
template<class T,int d>
inline UPPER_TRIANGULAR_MATRIX<T,d> operator*(const T a,const UPPER_TRIANGULAR_MATRIX<T,d>& A)
{return A*a;}

template<class T,int d>
inline UPPER_TRIANGULAR_MATRIX<T,d> operator+(const T a,const UPPER_TRIANGULAR_MATRIX<T,d>& A)
{return A+a;}

template<class T,int d>
inline UPPER_TRIANGULAR_MATRIX<T,d> operator-(const T a,const UPPER_TRIANGULAR_MATRIX<T,d>& A)
{return -A+a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,0> operator*(const DIAGONAL_MATRIX<T,0>& A,const UPPER_TRIANGULAR_MATRIX<T,0>& B)
{return UPPER_TRIANGULAR_MATRIX<T,0>();}

template<class T,int d>
inline UPPER_TRIANGULAR_MATRIX<T,d> operator+(const DIAGONAL_MATRIX<T,d>& A,const UPPER_TRIANGULAR_MATRIX<T,d>& B)
{return B+A;}

template<class T,int d>
inline UPPER_TRIANGULAR_MATRIX<T,d> operator-(const DIAGONAL_MATRIX<T,d>& A,const UPPER_TRIANGULAR_MATRIX<T,d>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,0>& A)
{return output_stream<<"[]";}
//#####################################################################
}
#endif

