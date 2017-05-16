//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX
//#####################################################################
#ifndef __DIAGONAL_MATRIX__
#define __DIAGONAL_MATRIX__
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/SCALAR_POLICY.h>
#include <Core/Vectors/VECTOR.h>
#include <Core/Vectors/ZERO_VECTOR.h>

namespace PhysBAM{

template<class T,int d_input> class DIAGONAL_MATRIX;

template<class T,int d> struct IS_SCALAR_BLOCK<DIAGONAL_MATRIX<T,d> >:public IS_SCALAR_BLOCK<T>{};
template<class T,int d> struct IS_SCALAR_VECTOR_SPACE<DIAGONAL_MATRIX<T,d> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,int d,class RW> struct IS_BINARY_IO_SAFE<DIAGONAL_MATRIX<T,d>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T,int d_input>
class DIAGONAL_MATRIX
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {d=d_input,m=d_input,n=d_input};

    VECTOR<T,d> x;

    DIAGONAL_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
    {}

    template<class T2> explicit
    DIAGONAL_MATRIX(const DIAGONAL_MATRIX<T2,d>& matrix_input)
        :x(matrix_input.x)
    {}

    explicit DIAGONAL_MATRIX(const VECTOR<T,d>& v)
        :x(v)
    {}

    explicit DIAGONAL_MATRIX(const T& x0)
        :x(x0)
    {STATIC_ASSERT(d==1);}

    DIAGONAL_MATRIX(const T& x0,const T& x1)
        :x(x0,x1)
    {STATIC_ASSERT(d==2);}

    DIAGONAL_MATRIX(const T& x0,const T& x1,const T& x2)
        :x(x0,x1,x2)
    {STATIC_ASSERT(d==3);}

    DIAGONAL_MATRIX(const T& x0,const T& x1,const T& x2,const T& x3)
        :x(x0,x1,x2,x3)
    {STATIC_ASSERT(d==4);}

    DIAGONAL_MATRIX(const T& x0,const T& x1,const T& x2,const T& x3,const T& x4)
        :x(x0,x1,x2,x3,x4)
    {STATIC_ASSERT(d==5);}

    DIAGONAL_MATRIX(const T& x0,const T& x1,const T& x2,const T& x3,const T& x4,const T& x5)
        :x(x0,x1,x2,x3,x4,x5)
    {STATIC_ASSERT(d==6);}

    int Rows() const
    {return d;}

    int Columns() const
    {return d;}

    T& operator()(const int i)
    {return x(i);}

    const T& operator()(const int i) const
    {return x(i);}

    T& operator()(const int i,const int j)
    {assert(i==j);return x(i);}

    const T& operator()(const int i,const int j) const
    {assert(i==j);return x(i);}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<d && i==j;}

    T First() const
    {return x[0];}

    T Last() const
    {return x[d-1];}

    bool operator==(const DIAGONAL_MATRIX& A) const
    {return x==A.x;}

    bool operator!=(const DIAGONAL_MATRIX& A) const
    {return x!=A.x;}

    DIAGONAL_MATRIX operator-() const
    {return DIAGONAL_MATRIX(-x);}

    DIAGONAL_MATRIX& operator+=(const DIAGONAL_MATRIX& A)
    {x+=A.x;return *this;}

    DIAGONAL_MATRIX& operator-=(const DIAGONAL_MATRIX& A)
    {x-=A.x;return *this;}

    DIAGONAL_MATRIX& operator+=(const T& a)
    {x+=a;return *this;}

    DIAGONAL_MATRIX& operator-=(const T& a)
    {x-=a;return *this;}

    DIAGONAL_MATRIX& operator*=(const T a)
    {x*=a;return *this;}

    DIAGONAL_MATRIX& operator/=(const T a)
    {x/=a;return *this;}

    MATRIX<T,d> operator+(const MATRIX<T,d>& A) const
    {MATRIX<T,d> r(A);
    for(int i=0;i<d;i++) r(i,i)+=x(i);
    return r;}

    MATRIX<T,d> operator-(const MATRIX<T,d>& A) const
    {return *this+-A;}

    DIAGONAL_MATRIX operator+(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x+A.x);}

    DIAGONAL_MATRIX operator-(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x-A.x);}

    DIAGONAL_MATRIX operator+(const T a) const
    {return DIAGONAL_MATRIX(x+a);}

    DIAGONAL_MATRIX operator-(const T a) const
    {return DIAGONAL_MATRIX(x-a);}

    DIAGONAL_MATRIX operator*(const T a) const
    {return DIAGONAL_MATRIX(x*a);}

    DIAGONAL_MATRIX operator/(const T a) const
    {return DIAGONAL_MATRIX(x/a);}

    VECTOR<T,d> operator*(const VECTOR<T,d>& v) const
    {return x*v;}

    ZERO_VECTOR<T,d> operator*(const ZERO_VECTOR<T,d>& v) const
    {return v;}

    DIAGONAL_MATRIX operator*(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x*A.x);}

    DIAGONAL_MATRIX& operator*=(const DIAGONAL_MATRIX& A)
    {x*=A.x;return *this;}

    DIAGONAL_MATRIX operator/(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x/A.x);}

    T Determinant() const
    {return x.Product();}

    DIAGONAL_MATRIX Inverse() const
    {return DIAGONAL_MATRIX(PhysBAM::Inverse(x));}

    VECTOR<T,d> Inverse_Times(const VECTOR<T,d>& v) const
    {return v/x;}

    VECTOR<T,d> Robust_Inverse_Times(const VECTOR<T,d>& v) const
    {VECTOR<T,d> r;
    for(int i=0;i<d;i++) r(i)=Robust_Divide(v(i),x(i));
    return r;}

    template<class T_MATRIX>
    auto Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& B) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(B.Columns()==3);decltype(*this*B.Derived().Transposed()) M((INITIAL_SIZE)B.Columns(),(INITIAL_SIZE)B.Rows());
    for(int k=0;k<B.Rows();k++) for(int i=0;i<B.Columns();i++) M(i,k)=x(i)*B(k,i);
    return M;}

    DIAGONAL_MATRIX Times_Transpose(const DIAGONAL_MATRIX& M) const
    {return *this*M;}

    DIAGONAL_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    template<class T_MATRIX>
    auto Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    MATRIX<T,d> Times_Transpose(const MATRIX<T,d>& A) const
    {MATRIX<T,d> r;
    for(int i=0;i<d;i++) for(int j=0;j<d;j++) r(i,j)=x(i)*A(j,i);
    return r;}

    MATRIX<T,d> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,d>& M) const
    {return (M**this).Transposed();}

    static VECTOR<T,1> Cofactor_Matrix_Helper(const VECTOR<T,1>& v,T& prod){prod=v.x;return VECTOR<T,1>(1);}
    static VECTOR<T,2> Cofactor_Matrix_Helper(const VECTOR<T,2>& v,T& prod){prod=v.Product();return VECTOR<T,2>(v.y,v.x);}
    template<int p> static VECTOR<T,p> Cofactor_Matrix_Helper(const VECTOR<T,p>& v,T& prod)
    {
        T p1,p2;
        VECTOR<T,p/2> v1;
        VECTOR<T,p-p/2> v2;
        v.Split(v1,v2);
        v1=Cofactor_Matrix_Helper(v1,p1);
        v2=Cofactor_Matrix_Helper(v2,p2);
        prod=p1*p2;
        return (v1*p2).Append_Elements(v2*p1);
    }

    DIAGONAL_MATRIX Cofactor_Matrix() const
    {T p;return DIAGONAL_MATRIX(Cofactor_Matrix_Helper(x,p));}

    SYMMETRIC_MATRIX<T,3> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {T yy=sqr(v.y),zz=sqr(v.z),bx=x.y*v.x,cx=x.z*v.x;return SYMMETRIC_MATRIX<T,3>(x.y*zz+x.z*yy,-cx*v.y,-bx*v.z,x.x*zz+cx*v.x,-v.y*v.z*x.x,x.x*yy+bx*v.x);}

    SYMMETRIC_MATRIX<T,2> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x.x*SYMMETRIC_MATRIX<T,2>(sqr(v.y),-v.x*v.y,sqr(v.x));}

    SYMMETRIC_MATRIX<T,1> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,1>& v) const
    {return SYMMETRIC_MATRIX<T,1>();}

    T Trace() const
    {return x.Sum();}

    T Dilational() const
    {return Trace()/d;}

    T Min() const
    {return x.Min();}

    T Max() const
    {return x.Max();}

    T Max_Abs() const
    {return x.Max_Abs();}

    static T Inner_Product(const DIAGONAL_MATRIX& A,const DIAGONAL_MATRIX& B)
    {return A.x.Dot(B.x);}

    T Inner_Product(const VECTOR<T,d>& a,const VECTOR<T,d>& b) const // inner product with respect to this matrix
    {return b.Dot(a*x);}

    T Inverse_Inner_Product(const VECTOR<T,d>& a,const VECTOR<T,d>& b) const // inner product with respect to the inverse of this matrix
    {return b.Dot(a/x);}

    T Frobenius_Norm_Squared() const
    {return x.Magnitude_Squared();}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    bool Positive_Definite() const
    {return x.All_Greater(VECTOR<T,d>());}

    bool Positive_Semidefinite() const
    {return x.All_Greater_Equal(VECTOR<T,d>());}

    DIAGONAL_MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    DIAGONAL_MATRIX Sqrt() const
    {return DIAGONAL_MATRIX(sqrt(x));}

    DIAGONAL_MATRIX Clamp_Min(const T a) const
    {return DIAGONAL_MATRIX(clamp_min(x,a));}

    DIAGONAL_MATRIX Clamp_Max(const T a) const
    {return DIAGONAL_MATRIX(clamp_max(x,a));}

    DIAGONAL_MATRIX Abs() const
    {return DIAGONAL_MATRIX(abs(x));}

    DIAGONAL_MATRIX Sign() const
    {DIAGONAL_MATRIX r;
    for(int i=0;i<d;i++) r.x(i)=sign(x(i));
    return r;}

    static DIAGONAL_MATRIX Identity_Matrix()
    {return DIAGONAL_MATRIX()+1;}

    VECTOR<T,d> To_Vector() const
    {return x;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x);}
//#####################################################################
};
// global functions

template<class T,int d> VECTOR<T,d>
Transpose_Times(const DIAGONAL_MATRIX<T,d>& a,const VECTOR<T,d>& b)
{return a.Transpose_Times(b);}

template<class T,int d> ZERO_VECTOR<T,d>
Transpose_Times(const DIAGONAL_MATRIX<T,d>& a,const ZERO_VECTOR<T,d>& b)
{return b;}

template<class T,int d> DIAGONAL_MATRIX<T,d>
Transpose_Times_Self(const DIAGONAL_MATRIX<T,d>& a)
{return a*a;}

template<class T,int d>
inline DIAGONAL_MATRIX<T,d> operator*(const T a,const DIAGONAL_MATRIX<T,d>& A)
{return A*a;}

template<class T,int d>
inline DIAGONAL_MATRIX<T,d> operator+(const T a,const DIAGONAL_MATRIX<T,d>& A)
{return A+a;}

template<class T,int d>
inline DIAGONAL_MATRIX<T,d> operator-(const T a,const DIAGONAL_MATRIX<T,d>& A)
{return -A+a;}

template<class T,int d>
inline MATRIX<T,d> operator+(const MATRIX<T,d>& A,const DIAGONAL_MATRIX<T,d>& B)
{return B+A;}

template<class T,int d>
inline MATRIX<T,d> operator-(const MATRIX<T,d>& A,const DIAGONAL_MATRIX<T,d>& B)
{return -B+A;}

template<class T,int d>
inline std::istream& operator>>(std::istream& input_stream,DIAGONAL_MATRIX<T,d>& A)
{return input_stream>>A.x;}

template<class T,int d>
inline std::ostream& operator<<(std::ostream& output_stream,const DIAGONAL_MATRIX<T,d>& A)
{return output_stream<<A.x;}

template<class T,int d>
inline DIAGONAL_MATRIX<T,d> log(const DIAGONAL_MATRIX<T,d>& A)
{return DIAGONAL_MATRIX<T,d>(log(A.x));}

template<class T,int d>
inline DIAGONAL_MATRIX<T,d> exp(const DIAGONAL_MATRIX<T,d>& A)
{return DIAGONAL_MATRIX<T,d>(exp(A.x));}
//#####################################################################
}
#endif
