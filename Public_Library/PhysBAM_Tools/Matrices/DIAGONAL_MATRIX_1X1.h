//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX_1X1
//#####################################################################
#ifndef __DIAGONAL_MATRIX_1X1__
#define __DIAGONAL_MATRIX_1X1__

#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

using ::std::log;

template<class T> struct IS_SCALAR_BLOCK<DIAGONAL_MATRIX<T,1> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<DIAGONAL_MATRIX<T,1> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<DIAGONAL_MATRIX<T,1>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class DIAGONAL_MATRIX<T,1>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=1};

    T x11;

    DIAGONAL_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(1))
        :x11(0)
    {
        STATIC_ASSERT(sizeof(DIAGONAL_MATRIX)==sizeof(T));assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(1));
    }

    template<class T2> explicit
    DIAGONAL_MATRIX(const DIAGONAL_MATRIX<T2,1>& matrix_input)
        :x11((T)matrix_input.x11)
    {}

    explicit DIAGONAL_MATRIX(const SYMMETRIC_MATRIX<T,1>& matrix_input)
        :x11(matrix_input.x11)
    {}

    explicit DIAGONAL_MATRIX(const MATRIX<T,1>& matrix_input)
        :x11(matrix_input.x11)
    {}

    DIAGONAL_MATRIX(const T y11)
        :x11(y11)
    {}

    explicit DIAGONAL_MATRIX(const VECTOR<T,1>& v)
        :x11(v.x)
    {}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 1;}

    T& operator()(const int i)
    {assert(i==0);return ((T*)this)[i];}

    const T& operator()(const int i) const
    {assert(i==0);return ((T*)this)[i];}

    T& operator()(const int i,const int j)
    {assert(i==0 && i==j);return ((T*)this)[i];}

    const T& operator()(const int i,const int j) const
    {assert(i==0 && i==j);return ((T*)this)[i];}

    bool Valid_Index(const int i,const int j) const
    {return i==0 && i==j;}

    T First() const
    {return x11;}

    T Last() const
    {return x11;}

    bool operator==(const DIAGONAL_MATRIX& A) const
    {return x11==A.x11;}

    bool operator!=(const DIAGONAL_MATRIX& A) const
    {return x11!=A.x11;}

    DIAGONAL_MATRIX operator-() const
    {return DIAGONAL_MATRIX(-x11);}

    DIAGONAL_MATRIX& operator+=(const DIAGONAL_MATRIX& A)
    {x11+=A.x11;return *this;}

    DIAGONAL_MATRIX& operator+=(const T& a)
    {x11+=a;return *this;}

    DIAGONAL_MATRIX& operator-=(const DIAGONAL_MATRIX& A)
    {x11-=A.x11;return *this;}

    DIAGONAL_MATRIX& operator-=(const T& a)
    {x11-=a;return *this;}

    DIAGONAL_MATRIX& operator*=(const T a)
    {x11*=a;return *this;}

    DIAGONAL_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;return *this;}

    DIAGONAL_MATRIX operator+(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11+A.x11);}

    MATRIX<T,1> operator+(const MATRIX<T,1>& A) const
    {return MATRIX<T,1>(x11+A.x[0],A.x[1],A.x[2]);}

    MATRIX<T,1> operator-(const MATRIX<T,1>& A) const
    {return MATRIX<T,1>(x11-A.x[0],-A.x[1],-A.x[2]);}

    DIAGONAL_MATRIX operator+(const T a) const
    {return DIAGONAL_MATRIX(x11+a);}

    DIAGONAL_MATRIX operator-(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11-A.x11);}

    DIAGONAL_MATRIX operator-(const T a) const
    {return DIAGONAL_MATRIX(x11-a);}

    DIAGONAL_MATRIX operator*(const T a) const
    {return DIAGONAL_MATRIX(a*x11);}

    DIAGONAL_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,1> operator*(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x11*v.x);}

    DIAGONAL_MATRIX operator*(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11*A.x11);}

    DIAGONAL_MATRIX& operator*=(const DIAGONAL_MATRIX& A)
    {return *this=*this*A;}

    DIAGONAL_MATRIX operator/(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11/A.x11);}

    T Determinant() const
    {return x11;}

    DIAGONAL_MATRIX Inverse() const
    {assert(x11!=0);return DIAGONAL_MATRIX(1/x11);}

    VECTOR<T,1> Solve_Linear_System(const VECTOR<T,1>& v) const
    {assert(x11!=0);return VECTOR<T,1>(v.x/x11);}

    VECTOR<T,1> Robust_Solve_Linear_System(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(Robust_Divide(v.x,x11));}

    DIAGONAL_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    DIAGONAL_MATRIX Cofactor_Matrix() const
    {return DIAGONAL_MATRIX(1);}

    T Trace() const
    {return x11;}

    T Dilational() const
    {return Trace();}

    T Min() const
    {return x11;}

    T Max() const
    {return x11;}

    T Max_Abs() const
    {return abs(x11);}

    static T Inner_Product(const DIAGONAL_MATRIX& A,const DIAGONAL_MATRIX& B)
    {return A.x11*B.x11;}

    T Inner_Product(const VECTOR<T,1>& a,const VECTOR<T,1>& b) const // inner product with respect to this matrix
    {return a.x*x11*b.x;}

    T Inverse_Inner_Product(const VECTOR<T,1>& a,const VECTOR<T,1>& b) const // inner product with respect to the inverse of this matrix
    {assert(x11!=0);return a.x/x11*b.x;}

    T Frobenius_Norm_Squared() const
    {return sqr(x11);}

    T Frobenius_Norm() const
    {return abs(x11);}

    bool Positive_Definite() const
    {return x11>0;}

    bool Positive_Semidefinite() const
    {return x11>=0;}

    DIAGONAL_MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    DIAGONAL_MATRIX Sqrt() const
    {return DIAGONAL_MATRIX(sqrt(x11));}

    DIAGONAL_MATRIX Clamp_Min(const T a) const
    {return DIAGONAL_MATRIX(clamp_min(x11,a));}

    DIAGONAL_MATRIX Clamp_Max(const T a) const
    {return DIAGONAL_MATRIX(clamp_max(x11,a));}

    DIAGONAL_MATRIX Abs() const
    {return DIAGONAL_MATRIX(abs(x11));}

    DIAGONAL_MATRIX Sign() const
    {return DIAGONAL_MATRIX(sign(x11));}

    static DIAGONAL_MATRIX Identity_Matrix()
    {return DIAGONAL_MATRIX(1);}

    VECTOR<T,1> To_Vector() const
    {return VECTOR<T,1>(x11);}

    template<class T_MATRIX>
    typename PRODUCT<DIAGONAL_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& B) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(B.Columns()==1);typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)B.Columns(),(INITIAL_SIZE)B.Rows());
    for(int k=0;k<B.Rows();k++) for(int i=0;i<B.Columns();i++) M(i,k)=(*this)(i,i)*B(k,i);return M;}

    MATRIX<T,1> Times_Transpose(const DIAGONAL_MATRIX& M) const
    {return MATRIX<T,1>(x11*M.x11);}

    MATRIX<T,1>Times_Transpose(const SYMMETRIC_MATRIX<T,1>& M) const
    {return MATRIX<T,1>(x11*M.x11);}

    MATRIX<T,1> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,1>& M) const
    {return MATRIX<T,1>(x11*M.x11);}

    SYMMETRIC_MATRIX<T,2> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x11*SYMMETRIC_MATRIX<T,2>(sqr(v.y),-v.x*v.y,sqr(v.x));}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x11);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x11);}

//#####################################################################
    static T Inner_Product_Conjugate(const DIAGONAL_MATRIX& A,const MATRIX<T,1>& Q,const DIAGONAL_MATRIX B);
//#####################################################################
};
// global functions
template<class T>
inline DIAGONAL_MATRIX<T,1> operator*(const T a,const DIAGONAL_MATRIX<T,1>& A)
{return A*a;}

template<class T>
inline DIAGONAL_MATRIX<T,1> operator+(const T a,const DIAGONAL_MATRIX<T,1>& A)
{return A+a;}

template<class T>
inline DIAGONAL_MATRIX<T,1> operator-(const T a,const DIAGONAL_MATRIX<T,1>& A)
{return -A+a;}

template<class T>
inline MATRIX<T,1> operator+(const MATRIX<T,1>& A,const DIAGONAL_MATRIX<T,1>& B)
{return B+A;}

template<class T>
inline MATRIX<T,1> operator-(const MATRIX<T,1>& A,const DIAGONAL_MATRIX<T,1>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const DIAGONAL_MATRIX<T,1>& A)
{return output_stream<<"["<<A.x11<<"]";}

template<class T>
inline DIAGONAL_MATRIX<T,1> log(const DIAGONAL_MATRIX<T,1>& A)
{return DIAGONAL_MATRIX<T,1>(log(A.x11));}

template<class T>
inline DIAGONAL_MATRIX<T,1> exp(const DIAGONAL_MATRIX<T,1>& A)
{return DIAGONAL_MATRIX<T,1>(exp(A.x11));}
//#####################################################################
// Function Inner_Product_Conjugate
//#####################################################################
template<class T> inline T DIAGONAL_MATRIX<T,1>::
Inner_Product_Conjugate(const DIAGONAL_MATRIX<T,1>& A,const MATRIX<T,1>& Q,const DIAGONAL_MATRIX<T,1> B)
{
    return A.x11*Q.x11*B.x11*Q.x11;
}
//#####################################################################
}
#endif
