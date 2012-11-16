//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX_0X0
//#####################################################################
#ifndef __DIAGONAL_MATRIX_0X0__
#define __DIAGONAL_MATRIX_0X0__

#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

using ::std::log;

template<class T> struct IS_SCALAR_BLOCK<DIAGONAL_MATRIX<T,0> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<DIAGONAL_MATRIX<T,0> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<DIAGONAL_MATRIX<T,0>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class DIAGONAL_MATRIX<T,0>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=0,n=0};

    DIAGONAL_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(0),INITIAL_SIZE nn=INITIAL_SIZE(0))
    {
        assert(mm==INITIAL_SIZE(0) && nn==INITIAL_SIZE(0));
    }

    template<class T2> explicit
    DIAGONAL_MATRIX(const DIAGONAL_MATRIX<T2,0>& matrix_input)
    {}

    explicit DIAGONAL_MATRIX(const VECTOR<T,0>& v)
    {}

    explicit DIAGONAL_MATRIX(const SYMMETRIC_MATRIX<T,0>& matrix_input)
    {}

    explicit DIAGONAL_MATRIX(const MATRIX<T,0>& matrix_input)
    {}

    int Rows() const
    {return 0;}

    int Columns() const
    {return 0;}

    T& operator()(const int i)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(const int i) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator()(const int i,const int j)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(const int i,const int j) const
    {PHYSBAM_FATAL_ERROR();}

    bool Valid_Index(const int i,const int j) const
    {return false;}

    bool operator==(const DIAGONAL_MATRIX& A) const
    {return true;}

    bool operator!=(const DIAGONAL_MATRIX& A) const
    {return false;}

    DIAGONAL_MATRIX operator-() const
    {return *this;}

    DIAGONAL_MATRIX& operator+=(const DIAGONAL_MATRIX& A)
    {return *this;}

    DIAGONAL_MATRIX& operator+=(const T& a)
    {return *this;}

    DIAGONAL_MATRIX& operator-=(const DIAGONAL_MATRIX& A)
    {return *this;}

    DIAGONAL_MATRIX& operator-=(const T& a)
    {return *this;}

    DIAGONAL_MATRIX& operator*=(const T a)
    {return *this;}

    DIAGONAL_MATRIX& operator/=(const T a)
    {return *this;}

    DIAGONAL_MATRIX operator+(const DIAGONAL_MATRIX& A) const
    {return *this;}

    MATRIX<T,0> operator+(const MATRIX<T,0>& A) const
    {return *this;}

    MATRIX<T,0> operator-(const MATRIX<T,0>& A) const
    {return *this;}

    DIAGONAL_MATRIX operator+(const T a) const
    {return *this;}

    DIAGONAL_MATRIX operator-(const DIAGONAL_MATRIX& A) const
    {return *this;}

    DIAGONAL_MATRIX operator-(const T a) const
    {return *this;}

    DIAGONAL_MATRIX operator*(const T a) const
    {return *this;}

    DIAGONAL_MATRIX operator/(const T a) const
    {return *this;}

    VECTOR<T,0> operator*(const VECTOR<T,0>& v) const
    {return VECTOR<T,0>();}

    DIAGONAL_MATRIX operator*(const DIAGONAL_MATRIX& A) const
    {return *this;}

    DIAGONAL_MATRIX& operator*=(const DIAGONAL_MATRIX& A)
    {return *this;}

    DIAGONAL_MATRIX operator/(const DIAGONAL_MATRIX& A) const
    {return *this;}

    T Determinant() const
    {return 1;}

    DIAGONAL_MATRIX Inverse() const
    {return *this;}

    VECTOR<T,0> Solve_Linear_System(const VECTOR<T,0>& v) const
    {return VECTOR<T,0>();}

    VECTOR<T,0> Robust_Solve_Linear_System(const VECTOR<T,0>& v) const
    {return VECTOR<T,0>();}

    DIAGONAL_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    DIAGONAL_MATRIX Cofactor_Matrix() const
    {return *this;}

    T Trace() const
    {return 0;}

    T Dilational() const
    {return 0;}

    T Min() const
    {return FLT_MAX;}

    T Max() const
    {return -FLT_MAX;}

    T Max_Abs() const
    {return 0;}

    static T Inner_Product(const DIAGONAL_MATRIX& A,const DIAGONAL_MATRIX& B)
    {return DIAGONAL_MATRIX();}

    T Inner_Product(const VECTOR<T,0>& a,const VECTOR<T,0>& b) const // inner product with respect to this matrix
    {return 0;}

    T Inverse_Inner_Product(const VECTOR<T,0>& a,const VECTOR<T,0>& b) const // inner product with respect to the inverse of this matrix
    {return 0;}

    T Frobenius_Norm_Squared() const
    {return 0;}

    T Frobenius_Norm() const
    {return 0;}

    bool Positive_Definite() const
    {return true;}

    bool Positive_Semidefinite() const
    {return true;}

    DIAGONAL_MATRIX Positive_Definite_Part() const
    {return *this;}

    DIAGONAL_MATRIX Sqrt() const
    {return *this;}

    DIAGONAL_MATRIX Clamp_Min(const T a) const
    {return *this;}

    DIAGONAL_MATRIX Clamp_Max(const T a) const
    {return *this;}

    DIAGONAL_MATRIX Abs() const
    {return *this;}

    DIAGONAL_MATRIX Sign() const
    {return *this;}

    static DIAGONAL_MATRIX Identity_Matrix()
    {return DIAGONAL_MATRIX();}

    VECTOR<T,0> To_Vector() const
    {return VECTOR<T,0>();}

    template<class T_MATRIX>
    typename PRODUCT<DIAGONAL_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& B) const
    {return *this*B.Transposed();}

    DIAGONAL_MATRIX Times_Transpose(const DIAGONAL_MATRIX& M) const
    {return *this*M;}

    static T Inner_Product_Conjugate(const DIAGONAL_MATRIX& A,const MATRIX<T,0>& Q,const DIAGONAL_MATRIX B)
    {return A.x11*Q.x11*B.x11*Q.x11;}

    SYMMETRIC_MATRIX<T,1> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,1>& v) const
    {return SYMMETRIC_MATRIX<T,1>();}

    template<class RW> void Read(std::istream& input)
    {}

    template<class RW> void Write(std::ostream& output) const
    {}
};
// global functions
template<class T>
inline DIAGONAL_MATRIX<T,0> operator*(const T a,const DIAGONAL_MATRIX<T,0>& A)
{return A*a;}

template<class T>
inline DIAGONAL_MATRIX<T,0> operator+(const T a,const DIAGONAL_MATRIX<T,0>& A)
{return A+a;}

template<class T>
inline DIAGONAL_MATRIX<T,0> operator-(const T a,const DIAGONAL_MATRIX<T,0>& A)
{return -A+a;}

template<class T>
inline MATRIX<T,0> operator+(const MATRIX<T,0>& A,const DIAGONAL_MATRIX<T,0>& B)
{return B+A;}

template<class T>
inline MATRIX<T,0> operator-(const MATRIX<T,0>& A,const DIAGONAL_MATRIX<T,0>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const DIAGONAL_MATRIX<T,0>& A)
{return output_stream<<"[]";}

template<class T>
inline DIAGONAL_MATRIX<T,0> log(const DIAGONAL_MATRIX<T,0>& A)
{return A;}

template<class T>
inline DIAGONAL_MATRIX<T,0> exp(const DIAGONAL_MATRIX<T,0>& A)
{return A;}
//#####################################################################
}
#endif
