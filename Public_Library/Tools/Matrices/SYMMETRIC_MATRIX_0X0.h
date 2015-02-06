//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_0X0
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_0X0__
#define __SYMMETRIC_MATRIX_0X0__

#include <Tools/Math_Tools/exchange_sort.h>
#include <Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK<SYMMETRIC_MATRIX<T,0> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<SYMMETRIC_MATRIX<T,0> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,0>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class SYMMETRIC_MATRIX<T,0>
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {m=0,n=0};

    SYMMETRIC_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(0),INITIAL_SIZE nn=INITIAL_SIZE(0))
    {
        assert(mm==INITIAL_SIZE(0) && nn==INITIAL_SIZE(0));
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,0>& matrix_input)
    {}

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,0>& matrix_input)
    {}

    explicit SYMMETRIC_MATRIX(const MATRIX<T,0>& matrix_input)
    {}

    void From_Matrix(const MATRIX<T,0>& matrix_input)
    {}

    int Rows() const
    {return 0;}

    int Columns() const
    {return 0;}

    VECTOR<T,0> Column(const int axis) const
    {PHYSBAM_FATAL_ERROR();}

    VECTOR<T,0> Row(const int axis) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator()(int i,int j)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(int i,int j) const
    {PHYSBAM_FATAL_ERROR();}

    bool Valid_Index(const int i,const int j) const
    {return false;}

    T& Element_Upper(int i,int j)
    {PHYSBAM_FATAL_ERROR();}

    const T& Element_Upper(int i,int j) const
    {PHYSBAM_FATAL_ERROR();}

    T& Element_Lower(int i,int j)
    {PHYSBAM_FATAL_ERROR();}

    const T& Element_Lower(int i,int j) const
    {PHYSBAM_FATAL_ERROR();}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {return true;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return false;}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(min(v1.x00,v2.x00));}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(max(v1.x00,v2.x00));}

    SYMMETRIC_MATRIX operator-() const
    {return *this;}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {return *this;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {return *this;}

    SYMMETRIC_MATRIX operator+(const T a) const
    {return *this;}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {return *this;}

    SYMMETRIC_MATRIX operator-(const T a) const
    {return *this;}

    SYMMETRIC_MATRIX operator*(const T a) const
    {return *this;}

    SYMMETRIC_MATRIX operator/(const T a) const
    {return *this;}

    VECTOR<T,0> operator*(const VECTOR<T,0>& v) const
    {return VECTOR<T,0>();}

    template<class T_MATRIX>
    typename PRODUCT<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {return *this*A.Transposed();}

    T Determinant() const
    {return 1;}

    SYMMETRIC_MATRIX Inverse() const
    {return *this;}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    MATRIX<T,0> Times_Transpose(const MATRIX<T,0>& A) const
    {return MATRIX<T,0,0>();}

    MATRIX<T,0> Times_Transpose(const DIAGONAL_MATRIX<T,0>& A) const
    {return MATRIX<T,0,0>();}

    MATRIX<T,0> Times_Transpose(const SYMMETRIC_MATRIX<T,0>& A) const
    {return MATRIX<T,0,0>();}

    MATRIX<T,0,0> Cross_Product_Matrix_Times(const VECTOR<T,0>& v) const
    {return MATRIX<T,0,0>();}

    VECTOR<T,0> Inverse_Times(const VECTOR<T,0>& b) const
    {return Inverse()*b;}

    VECTOR<T,0> Robust_Inverse_Times(const VECTOR<T,0>& b) const
    {return Inverse()*b;}

    T Trace() const
    {return 0;}

    static T Inner_Product(const SYMMETRIC_MATRIX& A,const SYMMETRIC_MATRIX& B)
    {return 0;}

    T Frobenius_Norm_Squared() const
    {return 0;}

    T Frobenius_Norm() const
    {return 0;}

    SYMMETRIC_MATRIX Cofactor_Matrix() const 
    {return *this;}

    VECTOR<T,0> Largest_Column() const
    {return VECTOR<T,0>();}

    VECTOR<T,0> Largest_Column_Normalized() const // 5 mults, 2 adds, 1 div, 1 sqrt
    {return VECTOR<T,0>();}

    T Max_Abs() const
    {return 0;}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,0>& u)
    {return SYMMETRIC_MATRIX();}

    static SYMMETRIC_MATRIX Symmetric_Outer_Product(const VECTOR<T,0>& u,const VECTOR<T,0>& v)
    {return SYMMETRIC_MATRIX();}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX();}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {return SYMMETRIC_MATRIX();}

    bool Positive_Definite() const
    {return true;}

    bool Positive_Semidefinite(const T tolerance=(T)1e-7) const
    {return true;}

    DIAGONAL_MATRIX<T,0> Fast_Eigenvalues() const
    {return *this;}

    SYMMETRIC_MATRIX Positive_Definite_Part() const
    {return *this;}

    DIAGONAL_MATRIX<T,0> Diagonal_Part() const
    {return *this;}

    VECTOR<T,0> Off_Diagonal_Part() const
    {return VECTOR<T,0>();}

    void Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,0>& eigenvalues,MATRIX<T,0>& eigenvectors) const
    {Solve_Eigenproblem(eigenvalues,eigenvectors);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,0>& A,const MATRIX<T,0>& B) // A^t*B and assume symmetric result
    {return SYMMETRIC_MATRIX();}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,0>& A,const UPPER_TRIANGULAR_MATRIX<T,0>& B) // A^t*B and assume symmetric result, 4 mults, 1 adds
    {return SYMMETRIC_MATRIX();}

    template<class RW> void Read(std::istream& input)
    {}

    template<class RW> void Write(std::ostream& output) const
    {}

    void Solve_Eigenproblem(DIAGONAL_MATRIX<T,0>& eigenvalues,MATRIX<T,0>& eigenvectors) const
    {}

    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,0>& A,const DIAGONAL_MATRIX<T,0>& B)
    {return SYMMETRIC_MATRIX<T,0>();}

    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,0>& A,const SYMMETRIC_MATRIX& B)
    {return SYMMETRIC_MATRIX<T,0>();}

    static SYMMETRIC_MATRIX Conjugate(const UPPER_TRIANGULAR_MATRIX<T,0>& A,const SYMMETRIC_MATRIX& B)
    {return SYMMETRIC_MATRIX<T,0>();}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,0>& A,const DIAGONAL_MATRIX<T,0>& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,0>& A,const SYMMETRIC_MATRIX& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,0>& A,const SYMMETRIC_MATRIX& B)
    {return Transpose_Times_With_Symmetric_Result(B*A,A);}

    MATRIX<T,0> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,0>& A) const
    {return (A**this).Transposed();}

    MATRIX<T,0> operator*(const DIAGONAL_MATRIX<T,0>& A) const
    {return MATRIX<T,0>();}

    MATRIX<T,0> operator*(const UPPER_TRIANGULAR_MATRIX<T,0>& A) const
    {return MATRIX<T,0>();}

    SYMMETRIC_MATRIX<T,1> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,1>& v) const
    {return SYMMETRIC_MATRIX<T,1>();}

    MATRIX<T,0,1> Times_Cross_Product_Matrix(const VECTOR<T,1>& v)
    {return MATRIX<T,0,1>();}

//#####################################################################
    SYMMETRIC_MATRIX operator+(const DIAGONAL_MATRIX<T,0>& A) const;
//#####################################################################
};
// global functions
template<class T>
inline SYMMETRIC_MATRIX<T,0> operator*(const T a,const SYMMETRIC_MATRIX<T,0>& A) // 4 mults
{return A*a;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> operator+(const T a,const SYMMETRIC_MATRIX<T,0>& A) // 2 adds
{return A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> operator-(const T a,const SYMMETRIC_MATRIX<T,0>& A) // 2 adds
{return -A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> clamp(const SYMMETRIC_MATRIX<T,0>& x,const SYMMETRIC_MATRIX<T,0>& xmin,const SYMMETRIC_MATRIX<T,0>& xmax)
{return x;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> clamp_min(const SYMMETRIC_MATRIX<T,0>& x,const SYMMETRIC_MATRIX<T,0>& xmin)
{return x;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> clamp_max(const SYMMETRIC_MATRIX<T,0>& x,const SYMMETRIC_MATRIX<T,0>& xmax)
{return x;}

template<class T>
inline std::ostream& operator<< (std::ostream& output_stream,const SYMMETRIC_MATRIX<T,0>& A)
{output_stream<<"[]";return output_stream;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> log(const SYMMETRIC_MATRIX<T,0>& A)
{return A;}

template<class T>
inline SYMMETRIC_MATRIX<T,0> exp(const SYMMETRIC_MATRIX<T,0>& A)
{return A;}

//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,0> operator*(const DIAGONAL_MATRIX<T,0>& D,const SYMMETRIC_MATRIX<T,0>& A)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,0> operator*(const UPPER_TRIANGULAR_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,0> operator*(const SYMMETRIC_MATRIX<T,0>& A,const MATRIX<T,0>& B)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,0> operator*(const SYMMETRIC_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,0> SYMMETRIC_MATRIX<T,0>::
operator+(const DIAGONAL_MATRIX<T,0>& A) const
{
    return SYMMETRIC_MATRIX<T,0>();
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,0>
operator+(const DIAGONAL_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return B+A;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,0>
operator+(const SYMMETRIC_MATRIX<T,0>& A,const UPPER_TRIANGULAR_MATRIX<T,0>& B)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,0>
operator+(const UPPER_TRIANGULAR_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,0>
operator-(const SYMMETRIC_MATRIX<T,0>& A,const UPPER_TRIANGULAR_MATRIX<T,0>& B)
{
    return MATRIX<T,0>();
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,0>
operator-(const UPPER_TRIANGULAR_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return -B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,0>
operator-(const DIAGONAL_MATRIX<T,0>& A,const SYMMETRIC_MATRIX<T,0>& B)
{
    return SYMMETRIC_MATRIX<T,0>();
}
//#####################################################################
}
#endif
