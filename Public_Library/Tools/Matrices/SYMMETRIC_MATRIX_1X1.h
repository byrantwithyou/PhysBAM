//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_1X1
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_1X1__
#define __SYMMETRIC_MATRIX_1X1__

#include <Tools/Math_Tools/exchange_sort.h>
#include <Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class T> struct is_scalar_BLOCK<SYMMETRIC_MATRIX<T,1> >:public is_scalar_BLOCK<T>{};
template<class T> struct is_scalar_VECTOR_SPACE<SYMMETRIC_MATRIX<T,1> >:public is_scalar_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,1>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class SYMMETRIC_MATRIX<T,1>
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=1};

    T x00;

    SYMMETRIC_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(1))
        :x00(T())
    {
        assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(1));
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,1>& matrix_input)
        :x00((T)matrix_input.x00)
    {}

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,1>& matrix_input)
        :x00(matrix_input.x.x)
    {}

    explicit SYMMETRIC_MATRIX(const MATRIX<T,1>& matrix_input)
        :x00(matrix_input.x00)
    {}

    SYMMETRIC_MATRIX(const T y00)
        :x00(y00)
    {}

    void From_Matrix(const MATRIX<T,1>& matrix_input)
    {x00=matrix_input(0,0);}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 1;}

    VECTOR<T,1> Column(const int axis) const
    {assert(axis==0);return VECTOR<T,1>(x00);}

    VECTOR<T,1> Row(const int axis) const
    {return Column(axis);}

    T& operator()(int i,int j)
    {assert(i==0 && j==0);return x00;}

    const T& operator()(int i,int j) const
    {assert(i==0 && j==0);return x00;}

    bool Valid_Index(const int i,const int j) const
    {return i==0 && j==0;}

    T& Element_Upper(int i,int j)
    {return (*this)(i,j);}

    const T& Element_Upper(int i,int j) const
    {return (*this)(i,j);}

    T& Element_Lower(int i,int j)
    {return (*this)(i,j);}

    const T& Element_Lower(int i,int j) const
    {return (*this)(i,j);}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {return x00==A.x00;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return !(*this==A);}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(min(v1.x00,v2.x00));}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(max(v1.x00,v2.x00));}

    SYMMETRIC_MATRIX operator-() const
    {return SYMMETRIC_MATRIX(-x00);}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {x00+=A.x00;return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {x00+=a;return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {x00-=A.x00;return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {x00-=a;return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {x00*=a;return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {assert(a!=0);x00/=a;return *this;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x00+A.x00);}

    SYMMETRIC_MATRIX operator+(const T a) const
    {return SYMMETRIC_MATRIX(x00+a);}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x00-A.x00);}

    SYMMETRIC_MATRIX operator-(const T a) const
    {return SYMMETRIC_MATRIX(x00-a);}

    SYMMETRIC_MATRIX operator*(const T a) const
    {return SYMMETRIC_MATRIX(a*x00);}

    SYMMETRIC_MATRIX operator/(const T a) const
    {assert(a!=0);return SYMMETRIC_MATRIX(x00/a);}

    VECTOR<T,1> operator*(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x00*v.x);}

    template<class T_MATRIX>
    typename PRODUCT<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)A.Columns(),(INITIAL_SIZE)A.Rows());A.Add_Times_Transpose(*this,A.Derived(),M);return M;}

    T Determinant() const
    {return x00;}

    SYMMETRIC_MATRIX Inverse() const
    {return SYMMETRIC_MATRIX(1/x00);}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    MATRIX<T,1> Times_Transpose(const MATRIX<T,1>& A) const
    {return *this*A.Transposed();}

    MATRIX<T,1> Times_Transpose(const DIAGONAL_MATRIX<T,1>& A) const
    {return *this*A;}

    MATRIX<T,1> Times_Transpose(const SYMMETRIC_MATRIX<T,1>& A) const
    {return *this*A;}

    MATRIX<T,0,1> Cross_Product_Matrix_Times(const VECTOR<T,1>& v) const
    {return MATRIX<T,0,1>();}

    VECTOR<T,1> Inverse_Times(const VECTOR<T,1>& b) const
    {return Inverse()*b;}

    VECTOR<T,1> Robust_Inverse_Times(const VECTOR<T,1>& b) const
    {return Inverse()*b;}

    T Trace() const
    {return x00;}

    static T Inner_Product(const SYMMETRIC_MATRIX& A,const SYMMETRIC_MATRIX& B)
    {return A.x00*B.x00;}

    T Frobenius_Norm_Squared() const
    {return x00*x00;}

    T Frobenius_Norm() const
    {return abs(x00);}

    SYMMETRIC_MATRIX Cofactor_Matrix() const
    {return SYMMETRIC_MATRIX(1);}

    VECTOR<T,1> Largest_Column() const
    {return VECTOR<T,1>(x00);}

    VECTOR<T,1> Largest_Column_Normalized() const // 5 mults, 2 adds, 1 div, 1 sqrt
    {return VECTOR<T,1>(sign_nonzero(x00));}

    T Max_Abs() const
    {return abs(x00);}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,1>& u)
    {return SYMMETRIC_MATRIX(u.x*u.x);}

    static SYMMETRIC_MATRIX Symmetric_Outer_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v)
    {return SYMMETRIC_MATRIX(2*v.x*u.x);}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX(1);}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {return SYMMETRIC_MATRIX(scale);}

    bool Positive_Definite() const
    {return x00>0;}

    bool Positive_Semidefinite(const T tolerance=(T)1e-7) const
    {return x00>=0;}

    DIAGONAL_MATRIX<T,1> Fast_Eigenvalues() const
    {return DIAGONAL_MATRIX<T,1>(x00);}

    SYMMETRIC_MATRIX Positive_Definite_Part() const
    {return SYMMETRIC_MATRIX(max((T)0,x00));}

    DIAGONAL_MATRIX<T,1> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,1>(x00);}

    VECTOR<T,0> Off_Diagonal_Part() const
    {return VECTOR<T,0>();}

    void Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,1>& eigenvalues,MATRIX<T,1>& eigenvectors) const
    {Solve_Eigenproblem(eigenvalues,eigenvectors);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,1>& A,const MATRIX<T,1>& B) // A^t*B and assume symmetric result
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,1>& A,const UPPER_TRIANGULAR_MATRIX<T,1>& B) // A^t*B and assume symmetric result, 4 mults, 1 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x00);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x00);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x00);}

    void Solve_Eigenproblem(DIAGONAL_MATRIX<T,1>& eigenvalues,MATRIX<T,1>& eigenvectors) const
    {eigenvectors=MATRIX<T,1>(1);eigenvectors=*this;}

    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,1>& A,const DIAGONAL_MATRIX<T,1>& B)
    {return SYMMETRIC_MATRIX<T,1>(sqr(A.x00)*B.x.x);}

    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,1>& A,const SYMMETRIC_MATRIX& B)
    {return SYMMETRIC_MATRIX<T,1>(sqr(A.x00)*B.x00);}

    static SYMMETRIC_MATRIX Conjugate(const DIAGONAL_MATRIX<T,1>& A,const SYMMETRIC_MATRIX& B)
    {return SYMMETRIC_MATRIX<T,1>(sqr(A.x.x)*B.x00);}

    static SYMMETRIC_MATRIX Conjugate(const UPPER_TRIANGULAR_MATRIX<T,1>& A,const SYMMETRIC_MATRIX& B)
    {return SYMMETRIC_MATRIX<T,1>(sqr(A.x00)*B.x00);}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,1>& A,const DIAGONAL_MATRIX<T,1>& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,1>& A,const SYMMETRIC_MATRIX& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,1>& A,const SYMMETRIC_MATRIX& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    MATRIX<T,1> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,1>& A) const
    {return (A**this).Transposed();}

    MATRIX<T,1> operator*(const DIAGONAL_MATRIX<T,1>& A) const
    {return MATRIX<T,1>(x00*A.x.x);}

    MATRIX<T,1> operator*(const UPPER_TRIANGULAR_MATRIX<T,1>& A) const
    {return MATRIX<T,1>(x00*A.x00);}

    SYMMETRIC_MATRIX<T,2> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x00*SYMMETRIC_MATRIX<T,2>(sqr(v.y),-v.x*v.y,sqr(v.x));}

    MATRIX<T,1,2> Times_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x00*MATRIX<T,1,2>::Cross_Product_Matrix(v);}

//#####################################################################
    SYMMETRIC_MATRIX operator+(const DIAGONAL_MATRIX<T,1>& A) const;
//#####################################################################
};
// global functions
template<class T>
inline SYMMETRIC_MATRIX<T,1> operator*(const T a,const SYMMETRIC_MATRIX<T,1>& A) // 4 mults
{return A*a;}

template<class T>
inline SYMMETRIC_MATRIX<T,1> operator+(const T a,const SYMMETRIC_MATRIX<T,1>& A) // 2 adds
{return A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,1> operator-(const T a,const SYMMETRIC_MATRIX<T,1>& A) // 2 adds
{return -A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,1> clamp(const SYMMETRIC_MATRIX<T,1>& x,const SYMMETRIC_MATRIX<T,1>& xmin,const SYMMETRIC_MATRIX<T,1>& xmax)
{return SYMMETRIC_MATRIX<T,1>(clamp(x.x00,xmin.x00,xmax.x00));}

template<class T>
inline SYMMETRIC_MATRIX<T,1> clamp_min(const SYMMETRIC_MATRIX<T,1>& x,const SYMMETRIC_MATRIX<T,1>& xmin)
{return SYMMETRIC_MATRIX<T,1>(clamp_min(x.x00,xmin.x00));}

template<class T>
inline SYMMETRIC_MATRIX<T,1> clamp_max(const SYMMETRIC_MATRIX<T,1>& x,const SYMMETRIC_MATRIX<T,1>& xmax)
{return SYMMETRIC_MATRIX<T,1>(clamp_max(x.x00,xmax.x00));}

template<class T>
inline std::ostream& operator<< (std::ostream& output_stream,const SYMMETRIC_MATRIX<T,1>& A)
{output_stream<<"["<<A.x00<<"]";return output_stream;}

template<class T>
inline SYMMETRIC_MATRIX<T,1> log(const SYMMETRIC_MATRIX<T,1>& A)
{return SYMMETRIC_MATRIX<T,1>(log(A.x00));}

template<class T>
inline SYMMETRIC_MATRIX<T,1> exp(const SYMMETRIC_MATRIX<T,1>& A)
{return SYMMETRIC_MATRIX<T,1>(exp(A.x00));}

//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,1> operator*(const DIAGONAL_MATRIX<T,1>& D,const SYMMETRIC_MATRIX<T,1>& A)
{
    return MATRIX<T,1>(D.x.x*A.x00);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,1> operator*(const UPPER_TRIANGULAR_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return MATRIX<T,1>(A.x00*B.x00);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,1> operator*(const SYMMETRIC_MATRIX<T,1>& A,const MATRIX<T,1>& B)
{
    return MATRIX<T,1>(A.x00*B.x00);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,1> operator*(const SYMMETRIC_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return MATRIX<T,1>(A.x00*B.x00);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,1> SYMMETRIC_MATRIX<T,1>::
operator+(const DIAGONAL_MATRIX<T,1>& A) const
{
    return SYMMETRIC_MATRIX<T,1>(x00+A.x.x);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,1>
operator+(const DIAGONAL_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return B+A;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,1>
operator+(const SYMMETRIC_MATRIX<T,1>& A,const UPPER_TRIANGULAR_MATRIX<T,1>& B)
{
    return MATRIX<T,1>(A.x00+B.x00);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,1>
operator+(const UPPER_TRIANGULAR_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,1>
operator-(const SYMMETRIC_MATRIX<T,1>& A,const UPPER_TRIANGULAR_MATRIX<T,1>& B)
{
    return MATRIX<T,1>(A.x00-B.x00);
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,1>
operator-(const UPPER_TRIANGULAR_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return -B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,1>
operator-(const DIAGONAL_MATRIX<T,1>& A,const SYMMETRIC_MATRIX<T,1>& B)
{
    return SYMMETRIC_MATRIX<T,1>(A.x.x-B.x00);
}
//#####################################################################
}
#endif
