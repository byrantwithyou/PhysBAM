//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX
//#####################################################################
#ifndef __MATRIX__
#define __MATRIX__

#include <Tools/Matrices/MATRIX_0X0.h>
#include <Tools/Matrices/MATRIX_1X1.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_2X3.h>
#include <Tools/Matrices/MATRIX_3X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_4X4.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int m_input,int n_input> // n_input=m_input
class MATRIX:public MATRIX_BASE<T,MATRIX<T,m_input,n_input> >
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    enum WORKAROUND1 {m=m_input,n=n_input,size=n_input*m_input};
    STATIC_ASSERT((!((m>=2 && m<=3 && n>=2 && n<=3) || (m==4 && n==4) || (m==0 && n==0)))); // 0x0, 1x1, 2x2, 2x3, 3x2, 3x3, and 4x4 are specialized
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,m_input,n_input> > BASE;
    using BASE::Frobenius_Norm_Squared;using BASE::operator*;using BASE::Transpose_Times;using BASE::Times_Transpose;using BASE::PLU_Solve;

    T x[m*n+(m*n==0)]; // pointer to the one dimensional data

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(m),INITIAL_SIZE nn=INITIAL_SIZE(n))
    {
        STATIC_ASSERT(sizeof(MATRIX)==(size+!size)*sizeof(T));assert(mm==INITIAL_SIZE(m) && nn==INITIAL_SIZE(n));
        for(int i=0;i<size;i++) x[i]=T();
    }

    MATRIX(const MATRIX& A)
        :BASE()
    {
        for(int i=0;i<size;i++) x[i]=A.x[i];
    }

    explicit MATRIX(const VECTOR<T,size>& column1)
    {
        STATIC_ASSERT(m==1 || n==1);
        for(int i=0;i<column1.Size();i++) x[i]=column1(i);
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(m==A.Rows() && n==A.Columns());for(int j=0;j<n;j++) for(int i=0;i<m;i++) (*this)(i,j)=A(i,j);
    }

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<m);assert((unsigned)j<n);return x[j*m+i];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<m);assert((unsigned)j<n);return x[j*m+i];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<m && (unsigned)j<n;}

    VECTOR<T,m>& Column(const int j)
    {assert((unsigned)j<n);return *(VECTOR<T,m>*)(x+m*j);}

    const VECTOR<T,m>& Column(const int j) const
    {assert((unsigned)j<n);return *(const VECTOR<T,n>*)(x+m*j);}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<size;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    MATRIX& operator=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]=A.x[i];return *this;}

    template<class T_MATRIX>
    MATRIX& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {assert(Rows()==A.Rows() && Columns()==A.Columns());for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) (*this)(i,j)=A(i,j);return *this;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<size;i++) x[i]*=a;return *this;}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]+=A.x[i];return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]-=A.x[i];return *this;}

    MATRIX operator+(const MATRIX& A) const
    {assert(n==A.n && m==A.m);MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]+A.x[i];return matrix;}

    MATRIX operator-(const MATRIX& A) const
    {assert(n==A.n && m==A.m);MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]-A.x[i];return matrix;}

    MATRIX operator-() const
    {MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=-x[i];return matrix;}

    MATRIX operator*(const T a) const
    {MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]*a;return matrix;}

    VECTOR<T,m> operator*(const VECTOR<T,n>& y) const
    {VECTOR<T,m> result;for(int j=0;j<n;j++) for(int i=0;i<m;i++) result(i)+=(*this)(i,j)*y(j);return result;}

    template<int p>
    MATRIX<T,m,p> operator*(const MATRIX<T,n,p>& A) const
    {MATRIX<T,m,p> matrix;for(int j=0;j<p;j++) for(int k=0;k<n;k++) for(int i=0;i<m;i++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,m,n> operator*(const SYMMETRIC_MATRIX<T,n>& A) const
    {MATRIX<T,m,n> matrix;for(int j=0;j<n;j++) for(int k=0;k<n;k++) for(int i=0;i<m;i++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,m,n> operator*(const DIAGONAL_MATRIX<T,n>& A) const
    {MATRIX<T,m,n> matrix;for(int j=0;j<n;j++) for(int i=0;i<m;i++) matrix(i,j)=(*this)(i,j)*A(j,j);return matrix;}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(n==A.m);MATRIX_MXN<T> matrix(m,A.n);for(int j=0;j<A.n;j++) for(int i=0;i<m;i++) for(int k=0;k<n;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,n,m> Transposed() const
    {MATRIX<T,n,m> matrix;for(int i=0;i<m;i++) for(int j=0;j<n;j++) matrix(j,i)=(*this)(i,j);return matrix;}

    VECTOR<T,n> Transpose_Times(const VECTOR<T,m>& y) const
    {VECTOR<T,n> result;for(int j=0;j<n;j++) for(int i=0;i<m;i++) result(j)+=(*this)(i,j)*y(i);return result;}

    template<int p>
    MATRIX<T,n,p> Transpose_Times(const MATRIX<T,m,p>& A) const
    {MATRIX<T,n,p> matrix;for(int j=0;j<p;j++) for(int i=0;i<n;i++) for(int k=0;k<m;k++) matrix(i,j)+=(*this)(k,i)*A(k,j);return matrix;}

    template<int p>
    MATRIX<T,m,p> Times_Transpose(const MATRIX<T,p,n>& A) const
    {MATRIX<T,m,p> matrix;for(int j=0;j<p;j++) for(int i=0;i<m;i++) for(int k=0;k<n;k++) matrix(i,j)+=(*this)(i,k)*A(j,k);return matrix;}

    MATRIX Times_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(n==3);MATRIX matrix;for(int i=0;i<m;i++) matrix.Set_Row(i,VECTOR<T,3>::Cross_Product(VECTOR<T,3>((*this)(i,0),(*this)(i,1),(*this)(i,2)),v));return matrix;}

    MATRIX Times_Cross_Product_Matrix_Transpose(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(n==3);MATRIX matrix;for(int i=0;i<m;i++) matrix.Set_Row(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(i,0),(*this)(i,1),(*this)(i,2))));return matrix;}

    MATRIX Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(m==3);MATRIX matrix;for(int i=0;i<n;i++) matrix.Set_Column(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(0,i),(*this)(1,i),(*this)(2,i))));return matrix;}

    MATRIX Cross_Product_Matrix_Transpose_Times(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(m==3);MATRIX matrix;for(int i=0;i<n;i++) matrix.Set_Column(i,VECTOR<T,3>::Cross_Product(VECTOR<T,3>((*this)(0,i),(*this)(1,i),(*this)(2,i)),v));return matrix;}

    MATRIX<T,m,2> Times_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(n==1);return (*this)*MATRIX<T,1,2>::Cross_Product_Matrix(v);}

    MATRIX<T,m,1> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(n==2);return Times_Transpose(MATRIX<T,1,2>::Cross_Product_Matrix(v));}

    MATRIX<T,1,n> Cross_Product_Matrix_Times(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(m==2);return MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this);}

    SYMMETRIC_MATRIX<T,1> Cross_Product_Matrix_Times_With_Symmetric_Result(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(m==2 && n==1);return SYMMETRIC_MATRIX<T,1>(MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this));}

    MATRIX<T,2,n> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(m==1);return MATRIX<T,1,2>::Cross_Product_Matrix(v).Transpose_Times(*this);}

    static MATRIX<T,1,2> Cross_Product_Matrix(const VECTOR<T,2>& v)
    {STATIC_ASSERT(m==1 && n==2);MATRIX<T,1,2> M;M(0,0)=-v.y;M(0,1)=v.x;return M;}

    SYMMETRIC_MATRIX<T,1> Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(const VECTOR<T,2>& v) const
    {STATIC_ASSERT((m==1 && n==2));return SYMMETRIC_MATRIX<T,1>(Times_Cross_Product_Matrix_Transpose(v));}

    MATRIX<T,m,1> Times_Cross_Product_Matrix(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(n==0);return MATRIX<T,m,1>();}

    MATRIX<T,m,0> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(n==1);return MATRIX<T,m,0>();}

    MATRIX<T,0,n> Cross_Product_Matrix_Times(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(m==1);return MATRIX<T,0,n>();}

    MATRIX<T,1,n> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(m==0);return MATRIX<T,1,n>();}

    SYMMETRIC_MATRIX<T,0> Cross_Product_Matrix_Times_With_Symmetric_Result(const VECTOR<T,1>& v) const
    {STATIC_ASSERT((m==1 && n==0));return SYMMETRIC_MATRIX<T,0>();}

    static MATRIX<T,0,1> Cross_Product_Matrix(const VECTOR<T,1>& v)
    {STATIC_ASSERT((m==0 && n==1));return MATRIX<T,0,1>();}

    MATRIX<T,n> Normal_Equations_Matrix() const
    {MATRIX<T,n> result;for(int j=0;j<n;j++) for(int i=0;i<n;i++) for(int k=0;k<m;k++) result(i,j)+=(*this)(k,i)*(*this)(k,j);return result;}

    VECTOR<T,n> Normal_Equations_Solve(const VECTOR<T,m>& b) const
    {MATRIX<T,n> A_transpose_A(Normal_Equations_Matrix());VECTOR<T,n> A_transpose_b(Transpose_Times(b));return A_transpose_A.Cholesky_Solve(A_transpose_b);}

    T Parallelepiped_Measure() const
    {STATIC_ASSERT(n==1);return sqrt(Frobenius_Norm_Squared());}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,x,m*n);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,x,m*n);}
};

template<class T,int m,int n>
inline MATRIX<T,m,n> operator*(const T a,const MATRIX<T,m,n>& A)
{return A*a;}
}
#endif
