//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_2X2
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_2X2__
#define __UPPER_TRIANGULAR_MATRIX_2X2__

#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_1X1.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <cmath>
#include <ostream>
namespace PhysBAM{

using ::std::sqrt;

template<class T> struct IS_SCALAR_BLOCK<UPPER_TRIANGULAR_MATRIX<T,2> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,2> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,2>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,2>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=2,n=2};

    union
    {
        struct {T x00,x01,x11;};
        T array[3];
    };

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(2),INITIAL_SIZE nn=INITIAL_SIZE(2))
        :x00(T()),x01(T()),x11(T())
    {
        STATIC_ASSERT(sizeof(UPPER_TRIANGULAR_MATRIX)==3*sizeof(T));assert(mm==INITIAL_SIZE(2) && nn==INITIAL_SIZE(2));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,2>& matrix_input)
        :x00(matrix_input.x00),x01(matrix_input.x01),x11(matrix_input.x11)
    {}

    UPPER_TRIANGULAR_MATRIX(const T x11_input,const T x12_input,const T x22_input)
        :x00(x11_input),x01(x12_input),x11(x22_input)
    {}

    ~UPPER_TRIANGULAR_MATRIX()
    {
        x00.~T();
        x01.~T();
        x11.~T();
    }

    int Rows() const
    {return 2;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<=(unsigned)j && (unsigned)j<2);return array[((j*(j+1))>>1)+i];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<=(unsigned)j && (unsigned)j<2);return array[((j*(j+1))>>1)+i];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<=(unsigned)j && (unsigned)j<2;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return x00==A.x00 && x01==A.x01 && x11==A.x11;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return !(*this==A);}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return UPPER_TRIANGULAR_MATRIX(-x00,-x01,-x11);}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00+=A.x00;x01+=A.x01;x11+=A.x11;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {x00+=a;x11+=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {x00-=A.x00;x01-=A.x01;x11-=A.x11;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {x00-=a;x11-=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00,x01+A.x01,x11+A.x11);}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,2>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00+A.x00,x01,x11+A.x11);}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00+a,x01,x11+a);}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00,x01-A.x01,x11-A.x11);}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,2>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x00-A.x00,x01,x11-A.x11);}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x00-a,x01,x11-a);}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this=*this*A;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {x00*=a;x01*=a;x11*=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x00*=s;x01*=s;x11*=s;return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(a*x00,a*x01,a*x11);}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,2> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x00*v.x+x01*v.y,x11*v.y);}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const // 4 mults, 1 add
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x00,x00*A.x01+x01*A.x11,x11*A.x11);}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,2>& A) const // 3 mults
    {return UPPER_TRIANGULAR_MATRIX(x00*A.x.x,x01*A.x.y,x11*A.x.y);}

    template<class T_MATRIX>
    T_MATRIX operator*(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==2);T_MATRIX M(INITIAL_SIZE(2),INITIAL_SIZE(A.Columns()));
    for(int j=0;j<A.Columns();j++) for(int k=0;k<2;k++) for(int i=0;i<=k;i++) M(i,j)+=(*this)(i,k)*A(k,j);
    return M;}

    template<class T_MATRIX>
    auto Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Columns()==2);decltype(A.Derived().Transposed()) M(INITIAL_SIZE(2),INITIAL_SIZE(A.Rows()));
    for(int j=0;j<A.Rows();j++) for(int k=0;k<2;k++) for(int i=0;i<=k;i++) M(i,j)+=(*this)(i,k)*A(j,k);
    return M;}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {assert(A.Rows()==2);T_MATRIX M(INITIAL_SIZE(2),INITIAL_SIZE(A.Columns()));
    for(int j=0;j<A.Columns();j++) for(int k=0;k<2;k++) for(int i=0;i<=k;i++) M(k,j)+=(*this)(i,k)*A(i,j);
    return M;}

    MATRIX<T,2,3> Times_Transpose(const MATRIX<T,3,2>& A) const
    {return MATRIX<T,2,3>(x00*A.x[0]+x01*A.x[3],x11*A.x[3],x00*A.x[1]+x01*A.x[4],x11*A.x[4],x00*A.x[2]+x01*A.x[5],x11*A.x[5]);}

    UPPER_TRIANGULAR_MATRIX Times_Transpose(const DIAGONAL_MATRIX<T,2>& A) const
    {return *this*A;}

    MATRIX<T,2> Times_Transpose(const UPPER_TRIANGULAR_MATRIX& A) const
    {return MATRIX<T,2>(x00*A.x00+x01*A.x01,x11*A.x01,x01*A.x11,x11*A.x11);}

    MATRIX<T,2> Times_Transpose(const SYMMETRIC_MATRIX<T,2>& A) const
    {return *this*A;}

    MATRIX<T,2> Transpose_Times(const UPPER_TRIANGULAR_MATRIX& A) const
    {return MATRIX<T,2>(x00*A.x00,x01*A.x00,x00*A.x01,x01*A.x01+x11*A.x11);}

    MATRIX<T,2> Transpose_Times(const SYMMETRIC_MATRIX<T,2>& A) const
    {return MATRIX<T,2>(x00*A.x00,x01*A.x00+x11*A.x10,x00*A.x10,x01*A.x10+x11*A.x11);}

    MATRIX<T,2> Transpose_Times(const DIAGONAL_MATRIX<T,2>& A) const
    {return MATRIX<T,2>(x00*A.x.x,x01*A.x.x,0,x11*A.x.y);}

    SYMMETRIC_MATRIX<T,2> Outer_Product_Matrix() const // 4 mults, 1 add
    {return SYMMETRIC_MATRIX<T,2>(x00*x00+x01*x01,x01*x11,x11*x11);}

    T Determinant() const
    {return x00*x11;}

    T Trace() const
    {return x00+x11;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {assert(x00!=0 && x11!=0);T one_over_x11=1/x00,one_over_x22=1/x11;
    return UPPER_TRIANGULAR_MATRIX(one_over_x11,-x01*one_over_x11*one_over_x22,one_over_x22);}

    VECTOR<T,2> Inverse_Times(const VECTOR<T,2>& b) const
    {return Inverse()*b;}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return UPPER_TRIANGULAR_MATRIX(x11,-x01,x00);}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX(1,0,1);}

    T Max_Abs() const
    {return maxabs(x00,x01,x11);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {return sqr(x00)+sqr(x01)+sqr(x11);}

    T Simplex_Minimum_Altitude() const
    {return Determinant()/max(abs(x00),sqrt(sqr(x11)+max(sqr(x01),sqr(x01-x00))));}

    MATRIX<T,2> Transposed() const
    {return MATRIX<T,2>(*this).Transposed();}

//#####################################################################
};
// global functions

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator*(const DIAGONAL_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B) // 3 mults
{return UPPER_TRIANGULAR_MATRIX<T,2>(A.x.x*B.x00,A.x.x*B.x01,A.x.y*B.x11);}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,2>& A)
{return output_stream<<"["<<A.x00<<" "<<A.x01<<" ; 0 "<<A.x11<<"]";}
//#####################################################################
}
#endif

