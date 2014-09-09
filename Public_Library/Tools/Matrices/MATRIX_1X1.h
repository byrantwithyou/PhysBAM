//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_1X1__
#define __MATRIX_1X1__

#include <Tools/Math_Tools/Robust_Arithmetic.h>
#include <Tools/Matrices/MATRIX_BASE.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

using std::log;
using std::exp;
using std::sqrt;

template<class T>
class MATRIX<T,1>:public MATRIX_BASE<T,MATRIX<T,1> >
{
    struct UNUSABLE{};
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=1};
    typedef MATRIX_BASE<T,MATRIX<T,1> > BASE;using BASE::operator*;
    typedef int HAS_UNTYPED_READ_WRITE;

    T x00;

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(1))
        :x00(T())
    {
        assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(1));
    }

    MATRIX(const MATRIX& matrix)
        :x00(matrix.x00)
    {}

    template<class T2>
    explicit MATRIX(const MATRIX<T2,1>& matrix)
        :x00(matrix.x00)
    {}

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
        :x00(A(0,0))
    {
        assert(A.Rows()==1 && A.Columns()==1);
    }

    explicit MATRIX(const T x00)
        :x00(x00)
    {}

    explicit MATRIX(const VECTOR<T,1>& v)
        :x00(v.x)
    {}

    MATRIX(const SYMMETRIC_MATRIX<T,1>& matrix_input)
        :x00(matrix_input.x00)
    {}

    MATRIX(const DIAGONAL_MATRIX<T,1>& matrix_input)
        :x00(matrix_input.x00)
    {}

    MATRIX& operator=(const MATRIX& matrix)
    {x00=matrix.x00;return *this;}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 1;}

    T& operator()(const int i,const int j=1)
    {assert(i==0 && j==0);return x00;}

    const T& operator()(const int i,const int j=1) const
    {assert(i==0 && j==0);return x00;}

    bool Valid_Index(const int i,const int j) const
    {return i==0 && j==0;}

    VECTOR<T,1> Column(const int j) const
    {assert(j==0);return VECTOR<T,1>(x00);}

    void Set_Column(const int j,const VECTOR<T,1>& v)
    {assert(j==0);x00=v.x;}

    VECTOR<T,1> Row(const int j) const
    {assert(j==0);return VECTOR<T,1>(x00);}

    void Set_Row(const int j,const VECTOR<T,1>& v)
    {assert(j==0);x00=v.x;}

    bool operator==(const MATRIX& A) const
    {return x00==A.x00;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,1> Column_Sum() const
    {return Column(0);}

    VECTOR<T,1> Column_Magnitudes() const
    {return VECTOR<T,1>(Column(0).Magnitude());}

    MATRIX Inverse() const
    {assert(x00);return MATRIX(1/x00);}

    VECTOR<T,1> Inverse_Times(const VECTOR<T,1>& b) const
    {return VECTOR<T,1>(b.x/x00);}

    VECTOR<T,1> Robust_Inverse_Times(const VECTOR<T,1>& b) const
    {return VECTOR<T,1>(Robust_Divide(b.x,x00));}

    MATRIX Cofactor_Matrix() const
    {return MATRIX(1);}

    MATRIX Normal_Equations_Matrix() const // 1 mult
    {return MATRIX(sqr(x00));}

    MATRIX operator-() const
    {return MATRIX(-x00);}

    MATRIX operator+(const T a) const
    {return MATRIX(x00+a);}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x00+A.x00);}

    MATRIX operator-(const T a) const
    {return MATRIX(x00-a);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x00-A.x00);}

    MATRIX operator*(const MATRIX& A) const
    {return MATRIX(x00*A.x00);}

    MATRIX operator*(const DIAGONAL_MATRIX<T,1>& A) const
    {return MATRIX(x00*A.x);}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x00);}

    MATRIX operator/(const T a) const
    {return MATRIX(x00/a);}

    VECTOR<T,1> operator*(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x00*v.x);}

    MATRIX& operator+=(const MATRIX& A)
    {x00+=A.x00;return *this;}

    MATRIX& operator+=(const DIAGONAL_MATRIX<T,1>& A)
    {x00+=A.x.x;return *this;}

    MATRIX& operator+=(const SYMMETRIC_MATRIX<T,1>& A)
    {x00+=A.x00;return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {x00-=A.x00;return *this;}

    MATRIX& operator-=(const DIAGONAL_MATRIX<T,1>& A)
    {x00-=A.x.x;return *this;}

    MATRIX& operator-=(const SYMMETRIC_MATRIX<T,1>& A)
    {x00-=A.x00;return *this;}

    MATRIX& operator+=(const T& a)
    {x00+=a;return *this;}

    MATRIX& operator-=(const T& a)
    {x00-=a;return *this;}

    MATRIX& operator*=(const T a)
    {x00*=a;return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {x00*=A.x00;return *this;}

    MATRIX& operator/=(const T a)
    {x00/=a;return *this;}

    VECTOR<T,1> Transpose_Times(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x00*v.x);}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {return *this*A.Derived();}

    MATRIX<T,1> Transpose_Times(const DIAGONAL_MATRIX<T,1>& A) const
    {return *this*A;}

    MATRIX<T,1> Transpose_Times(const SYMMETRIC_MATRIX<T,1>& A) const
    {return *this*A;}

    MATRIX<T,1> Transpose_Times(const UPPER_TRIANGULAR_MATRIX<T,1>& A) const
    {return *this*A;}

    template<class T_MATRIX>
    typename TRANSPOSE<T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {return *this*A.Derived().Transposed();}

    MATRIX Times_Transpose(const DIAGONAL_MATRIX<T,1>& A) const
    {return *this*A;}

    MATRIX Times_Transpose(const SYMMETRIC_MATRIX<T,1>& A) const
    {return *this*A;}

    static MATRIX Outer_Product(const VECTOR<T,1>& u)
    {return MATRIX(u.x*u.x);}

    static MATRIX Outer_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v)
    {return MATRIX(u.x*v.x);}

    MATRIX<T,1>  Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,1>& M) const
    {return *this*M;}

    void Transpose()
    {}

    MATRIX Transposed() const
    {return *this;}

    T Trace() const
    {return x00;}

    T Determinant() const
    {return x00;}

    MATRIX<T,1> Fast_Eigenvalues() const
    {return *this;}

    T Max() const
    {return x00;}

    T Min() const
    {return x00;}

    T Frobenius_Norm() const
    {return abs(x00);}

    T Frobenius_Norm_Squared() const
    {return sqr(x00);}

    T Inner_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v) const
    {return x00*u.x*v.x;}

    T Inverse_Inner_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v) const
    {return u.x*v.x/x00;}

    MATRIX<T,1,2> Times_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x00*MATRIX<T,1,2>::Cross_Product_Matrix(v);}

    MATRIX<T,0,1> Cross_Product_Matrix_Times(const VECTOR<T,1>& v)
    {return MATRIX<T,0,1>();}

    bool Positive_Definite() const
    {return x00>0;}

    bool Positive_Semidefinite() const
    {return x00>=0;}

    MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    MATRIX Diagonal_Part() const
    {return *this;}

    MATRIX Sqrt() const
    {return MATRIX(sqrt(x00));}

    MATRIX Clamp_Min(const T a) const
    {return MATRIX(clamp_min(x00,a));}

    MATRIX Clamp_Max(const T a) const
    {return MATRIX(clamp_max(x00,a));}

    MATRIX Abs() const
    {return MATRIX(abs(x00));}

    MATRIX Sign() const
    {return MATRIX(sign(x00));}

    static MATRIX Identity_Matrix()
    {return MATRIX((T)1);}

    MATRIX Symmetric_Part() const
    {return *this;}

    VECTOR<T,1> To_Vector() const
    {return VECTOR<T,1>(x00);}

    static MATRIX Conjugate(const MATRIX& A,const MATRIX& B)
    {return A*B*A;}

    void Solve_Eigenproblem(MATRIX& D,MATRIX& V) const
    {Fast_Solve_Eigenproblem(D,V);}

    void Fast_Solve_Eigenproblem(MATRIX& D,MATRIX& V) const
    {V.x00=1;D=*this;}
    
    void Fast_Singular_Value_Decomposition(MATRIX& U,MATRIX& D,MATRIX& V) const
    {U.x00=V.x00=1;D=*this;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x00);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x00);}
//#####################################################################
};

template<class T>
inline MATRIX<T,1> operator*(const T a,const MATRIX<T,1>& A)
{return A*a;}

template<class T>
inline MATRIX<T,1> operator+(const T a,const MATRIX<T,1>& A)
{return A+a;}

template<class T>
inline MATRIX<T,1> operator-(const T a,const MATRIX<T,1>& A)
{return MATRIX<T,1>(a-A.x00);}

template<class T>
inline MATRIX<T,1> operator+(const SYMMETRIC_MATRIX<T,1>& A,const MATRIX<T,1>& B)
{return MATRIX<T,1>(A.x00+B.x00);}

template<class T>
inline MATRIX<T,1> operator-(const SYMMETRIC_MATRIX<T,1>& A,const MATRIX<T,1>& B)
{return MATRIX<T,1>(A.x00-B.x00);}

template<class T>
inline MATRIX<T,1> clamp(const MATRIX<T,1>& x,const MATRIX<T,1>& xmin,const MATRIX<T,1>& xmax)
{return MATRIX<T,1>(clamp(x.x00,xmin.x00,xmax.x00));}

template<class T>
inline MATRIX<T,1> clamp_min(const MATRIX<T,1>& x,const MATRIX<T,1>& xmin)
{return MATRIX<T,1>(clamp_min(x.x00,xmin.x00));}

template<class T>
inline MATRIX<T,1> clamp_max(const MATRIX<T,1>& x,const MATRIX<T,1>& xmax)
{return MATRIX<T,1>(clamp_max(x.x00,xmax.x00));}

template<class T>
inline MATRIX<T,1> log(const MATRIX<T,1>& A)
{return MATRIX<T,1>(log(A.x00));}

template<class T>
inline MATRIX<T,1> exp(const MATRIX<T,1>& A)
{return MATRIX<T,1>(exp(A.x00));}

template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,1>& A)
{FILE_UTILITIES::Ignore(input,'[');input>>A.x00;FILE_UTILITIES::Ignore(input,']');return input;}
}
#endif
