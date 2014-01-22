//#####################################################################
// Copyright 2003-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_3X3
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_3X3__
#define __UPPER_TRIANGULAR_MATRIX_3X3__

#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <ostream>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK<UPPER_TRIANGULAR_MATRIX<T,3> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,3> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,3>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,3>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x00,x01,x11,x02,x12,x22;

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
        :x00(T()),x01(T()),x11(T()),x02(T()),x12(T()),x22(T())
    {
        STATIC_ASSERT(sizeof(UPPER_TRIANGULAR_MATRIX)==6*sizeof(T));assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,3>& matrix_input)
        :x00(matrix_input.x00),x01(matrix_input.x01),x11(matrix_input.x11),x02(matrix_input.x02),x12(matrix_input.x12),x22(matrix_input.x22)
    {}

    UPPER_TRIANGULAR_MATRIX(const T x11_input,const T x12_input,const T x22_input,const T x13_input,const T x23_input,const T x33_input)
        :x00(x11_input),x01(x12_input),x11(x22_input),x02(x13_input),x12(x23_input),x22(x33_input)
    {}

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<=(unsigned)j && (unsigned)j<3);return ((T*)this)[((j*(j+1))>>1)+i];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<=(unsigned)j && (unsigned)j<3);return ((const T*)this)[((j*(j+1))>>1)+i];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<=j && (unsigned)j<3;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return x00==A.x00 && x01==A.x01 && x11==A.x11 && x02==A.x02 && x12==A.x12 && x22==A.x22;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return !(*this==A);}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return UPPER_TRIANGULAR_MATRIX(-x00,-x01,-x11,-x02,-x12,-x22);}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00+=A.x00;x01+=A.x01;x11+=A.x11;x02+=A.x02;x12+=A.x12;x22+=A.x22;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {x00+=a;x11+=a;x22+=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00-=A.x00;x01-=A.x01;x11-=A.x11;x02-=A.x02;x12-=A.x12;x22-=A.x22;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {x00-=a;x11-=a;x22-=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00,x01+A.x01,x11+A.x11,x02+A.x02,x12+A.x12,x22+A.x22);}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00,x01,x11+A.x11,x02,x12,x22+A.x22);}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00+a,x01,x11+a,x02,x12,x22+a);}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00,x01-A.x01,x11-A.x11,x02-A.x02,x12-A.x12,x22-A.x22);}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00,x01,x11-A.x11,x02,x12,x22-A.x22);}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00-a,x01,x11-a,x02,x12,x22-a);}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this=*this*A;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {x00*=a;x01*=a;x11*=a;x02*=a;x12*=a;x22*=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x00*=s;x01*=s;x11*=s;x02*=s;x12*=s;x22*=s;return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(a*x00,a*x01,a*x11,a*x02,a*x12,a*x22);}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return UPPER_TRIANGULAR_MATRIX(s*x00,s*x01,s*x11,s*x02,s*x12,s*x22);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const // 6 mults, 3 adds
    {return VECTOR<T,3>(x00*v.x+x01*v.y+x02*v.z,x11*v.y+x12*v.z,x22*v.z);}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x00,x00*A.x01+x01*A.x11,x11*A.x11,x00*A.x02+x01*A.x12+x02*A.x22,x11*A.x12+x12*A.x22,x22*A.x22);}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x.x,x01*A.x.y,x11*A.x.y,x02*A.x.z,x12*A.x.z,x22*A.x.z);}

    template<class T_MATRIX>
    T_MATRIX operator*(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==3);T_MATRIX M(INITIAL_SIZE(3),INITIAL_SIZE(A.Columns()));
    for(int j=0;j<A.Columns();j++) for(int k=0;k<3;k++) for(int i=0;i<=k;i++) M(i,j)+=(*this)(i,k)*A(k,j);return M;}

    template<class T_MATRIX>
    typename TRANSPOSE<T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Columns()==3);typename TRANSPOSE<T_MATRIX>::TYPE M(INITIAL_SIZE(3),INITIAL_SIZE(A.Rows()));
    for(int j=0;j<A.Rows();j++) for(int k=0;k<3;k++) for(int i=0;i<=k;i++) M(i,j)+=(*this)(i,k)*A(j,k);return M;}

    UPPER_TRIANGULAR_MATRIX Times_Transpose(const DIAGONAL_MATRIX<T,3>& A) const
    {return *this*A;}

    MATRIX<T,3> Times_Transpose(const SYMMETRIC_MATRIX<T,3>& A) const
    {return *this*A;}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==3);T_MATRIX M(INITIAL_SIZE(3),INITIAL_SIZE(A.Columns()));
    for(int j=0;j<A.Columns();j++) for(int k=0;k<3;k++) for(int i=0;i<=k;i++) M(k,j)+=(*this)(i,k)*A(i,j);return M;}

    MATRIX<T,3> Transpose_Times(const DIAGONAL_MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x00*A.x.x,x01*A.x.x,x02*A.x.x,0,x11*A.x.y,x12*A.x.y,0,0,x22*A.x.z);}

    MATRIX<T,3> Transpose_Times(const SYMMETRIC_MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x00*A.x00,x01*A.x00+x11*A.x10,x02*A.x00+x12*A.x10+x22*A.x20,x00*A.x10,x01*A.x10+x11*A.x11,x02*A.x10+x12*A.x11+x22*A.x21,x00*A.x20,x01*A.x20+x11*A.x21,x02*A.x20+x12*A.x21+x22*A.x22);}

    MATRIX<T,3> Transpose_Times(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x00*A.x00,x01*A.x00,x02*A.x00,x00*A.x01,x01*A.x01+x11*A.x11,x02*A.x01+x12*A.x11,x00*A.x02,x01*A.x02+x11*A.x12,x02*A.x02+x12*A.x12+x22*A.x22);}

    MATRIX<T,3> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x00*A.x00+x01*A.x01+x02*A.x02,x11*A.x01+x12*A.x02,x22*A.x02,x01*A.x11+x02*A.x12,x11*A.x11+x12*A.x12,x22*A.x12,x02*A.x22,x12*A.x22,x22*A.x22);}

    T Determinant() const
    {return x00*x11*x22;}

    T Trace() const
    {return x00+x11+x22;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {T determinant=x00*x11*x22;assert(determinant!=0);T s=1/determinant;
    return s*UPPER_TRIANGULAR_MATRIX(x11*x22,-x01*x22,x00*x22,x01*x12-x11*x02,-x00*x12,x00*x11);}

    VECTOR<T,3> Solve_Linear_System(const VECTOR<T,3>& b) const
    {return Cofactor_Matrix()*(b/Determinant());}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return UPPER_TRIANGULAR_MATRIX(x11*x22,-x01*x22,x00*x22,x01*x12-x11*x02,-x00*x12,x00*x11);}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX(1,0,1,0,0,1);}

    T Max_Abs() const
    {return maxabs(x00,x01,x11,x02,x12,x22);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {return sqr(x00)+sqr(x01)+sqr(x11)+sqr(x02)+sqr(x12)+sqr(x22);}

    T Simplex_Minimum_Altitude() const
    {return MATRIX<T,3>(*this).Simplex_Minimum_Altitude();}

    MATRIX<T,3> Transposed() const
    {return MATRIX<T,3>(*this).Transposed();}

//#####################################################################
};
// global functions

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator*(const DIAGONAL_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{return UPPER_TRIANGULAR_MATRIX<T,3>(A.x.x*B.x00,A.x.x*B.x01,A.x.y*B.x11,A.x.x*B.x02,A.x.y*B.x12,A.x.z*B.x22);}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,3>& A)
{return output_stream<<"["<<A.x00<<" "<<A.x01<<" "<<A.x02<<" ; 0 "<<A.x11<<" "<<A.x12<<" ; 0 0 "<<A.x22<<"]";}
//#####################################################################
}
#endif
