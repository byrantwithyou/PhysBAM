//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Duc Nguyen, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_2X2
//#####################################################################
#ifndef __MATRIX_2X2__
#define __MATRIX_2X2__

#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX_BASE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Core/Vectors/ZERO_VECTOR.h>
namespace PhysBAM{

template<class T>
class MATRIX<T,2>:public MATRIX_BASE<T,MATRIX<T,2> >
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {m=2,n=2};
    typedef MATRIX_BASE<T,MATRIX<T,2> > BASE;using BASE::operator*;using BASE::operator+;using BASE::Transpose_Times;using BASE::Times_Transpose;

    T x[4];

    MATRIX()
    {
        x[0]=T();x[1]=T();x[2]=T();x[3]=T();
    }

    explicit MATRIX(INITIAL_SIZE mm,INITIAL_SIZE nn=INITIAL_SIZE(2))
    {
        assert(mm==INITIAL_SIZE(2) && nn==INITIAL_SIZE(2));
        x[0]=T();x[1]=T();x[2]=T();x[3]=T();
    }

    MATRIX(const MATRIX& matrix_input)
    {
        x[0]=matrix_input.x[0];x[1]=matrix_input.x[1];x[2]=matrix_input.x[2];x[3]=matrix_input.x[3];
    }

    template<class T2> explicit
    MATRIX(const MATRIX<T2,2>& matrix_input)
    {
        x[0]=(T)matrix_input.x[0];x[1]=(T)matrix_input.x[1];x[2]=(T)matrix_input.x[2];x[3]=(T)matrix_input.x[3];
    }

    MATRIX(const SYMMETRIC_MATRIX<T,2>& matrix_input)
    {
        x[0]=matrix_input.x00;x[1]=x[2]=matrix_input.x10;x[3]=matrix_input.x11;
    }

    MATRIX(const DIAGONAL_MATRIX<T,2>& matrix_input)
    {
        x[0]=matrix_input.x.x;x[1]=T();x[2]=T();x[3]=matrix_input.x.y;
    }

    MATRIX(const UPPER_TRIANGULAR_MATRIX<T,2>& matrix_input)
    {
        x[0]=matrix_input.x00;x[2]=matrix_input.x01;x[3]=matrix_input.x11;x[1]=T();
    }

    MATRIX(const T x00,const T x10,const T x01,const T x11)
    {
        x[0]=x00;x[1]=x10;x[2]=x01;x[3]=x11;
    }

    MATRIX(const VECTOR<T,2> & column1,const VECTOR<T,2> & column2)
    {
        x[0]=column1.x;x[1]=column1.y;x[2]=column2.x;x[3]=column2.y;
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==2 && A.Columns()==2);
        x[0]=A(0,0);x[1]=A(1,0);x[2]=A(0,1);x[3]=A(1,1);
    }

    MATRIX& operator=(const MATRIX& matrix_input)
    {
        x[0]=matrix_input.x[0];x[1]=matrix_input.x[1];x[2]=matrix_input.x[2];x[3]=matrix_input.x[3];return *this;
    }

    template<class T_MATRIX>
    MATRIX& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(Rows()==A.Rows() && Columns()==A.Columns());
        x[0]=A(0,0);x[1]=A(1,0);x[2]=A(0,1);x[3]=A(1,1);
        return *this;
    }

    int Rows() const
    {return 2;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<2);assert((unsigned)j<2);return x[i+2*j];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<2);assert((unsigned)j<2);return x[i+2*j];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<2 && (unsigned)j<2;}

    VECTOR<T,2> Column(const int j) const
    {assert((unsigned)j<2);return VECTOR<T,2>(x[2*j],x[2*j+1]);}

    void Set_Column(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<2);x[2*j]=v.x;x[2*j+1]=v.y;}

    void Add_Column(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<2);x[2*j]+=v.x;x[2*j+1]+=v.y;}

    VECTOR<T,2> Row(const int j) const
    {assert((unsigned)j<2);return VECTOR<T,2>(x[j],x[j+2]);}

    void Set_Row(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<2);x[j]=v(0);x[j+2]=v(1);}

    void Add_Row(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<2);x[j]+=v(0);x[j+2]+=v(1);}

    bool operator==(const MATRIX& A) const
    {return x[0]==A.x[0] && x[1]==A.x[1] && x[2]==A.x[2] && x[3]==A.x[3];}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,2> Column_Sum() const
    {return VECTOR<T,2>(x[0]+x[2],x[1]+x[3]);}

    VECTOR<T,2> Normalize_Columns()
    {VECTOR<T,2> a=Column(0),b=Column(1),n(a.Normalize(),b.Normalize());*this=MATRIX<T,2>(a,b);return n;}

    VECTOR<T,2> Column_Magnitudes() const
    {return VECTOR<T,2>(Column(0).Magnitude(),Column(1).Magnitude());}

    MATRIX operator-() const
    {return MATRIX(-x[0],-x[1],-x[2],-x[3]);}

    MATRIX& operator+=(const MATRIX& A)
    {x[0]+=A.x[0];x[1]+=A.x[1];x[2]+=A.x[2];x[3]+=A.x[3];return *this;}

    MATRIX& operator+=(const T a)
    {x[0]+=a;x[3]+=a;return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {x[0]-=A.x[0];x[1]-=A.x[1];x[2]-=A.x[2];x[3]-=A.x[3];return *this;}

    MATRIX& operator-=(const T a)
    {x[0]-=a;x[3]-=a;return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {return *this=*this*A;}

    MATRIX& operator*=(const T a)
    {x[0]*=a;x[1]*=a;x[2]*=a;x[3]*=a;return *this;}

    MATRIX& operator/=(const T a)
    {assert(a!=0);return *this*=(1/a);}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3]);}

    MATRIX operator+(const T a) const
    {return MATRIX(x[0]+a,x[1],x[2],x[3]+a);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3]);}

    MATRIX operator-(const T a) const
    {return MATRIX(x[0]-a,x[1],x[2],x[3]-a);}

    MATRIX operator*(const MATRIX& A) const
    {return MATRIX(x[0]*A.x[0]+x[2]*A.x[1],x[1]*A.x[0]+x[3]*A.x[1],x[0]*A.x[2]+x[2]*A.x[3],x[1]*A.x[2]+x[3]*A.x[3]);}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x[0],a*x[1],a*x[2],a*x[3]);}

    MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return MATRIX(s*x[0],s*x[1],s*x[2],s*x[3]);}

    VECTOR<T,2> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x[0]*v.x+x[2]*v.y,x[1]*v.x+x[3]*v.y);}

    T Determinant() const
    {return x[0]*x[3]-x[1]*x[2];}

    T Parallelepiped_Measure() const
    {return Determinant();}

    void Invert()
    {*this=Inverse();}

    MATRIX Inverse() const
    {T one_over_determinant=1/(x[0]*x[3]-x[1]*x[2]);
    return MATRIX(one_over_determinant*x[3],-one_over_determinant*x[1],-one_over_determinant*x[2],one_over_determinant*x[0]);}

    MATRIX Inverse_Transposed() const
    {return Inverse().Transposed();}

    VECTOR<T,2> Inverse_Times(const VECTOR<T,2>& b) const
    {T one_over_determinant=1/(x[0]*x[3]-x[1]*x[2]);
    return one_over_determinant*VECTOR<T,2>(x[3]*b.x-x[2]*b.y,x[0]*b.y-x[1]*b.x);}

    VECTOR<T,2> Robust_Inverse_Times(const VECTOR<T,2>& b) const
    {T determinant=Determinant();
    VECTOR<T,2> unscaled_result=VECTOR<T,2>(x[3]*b.x-x[2]*b.y,x[0]*b.y-x[1]*b.x);
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;}

    void Transpose()
    {exchange(x[1],x[2]);}

    MATRIX Transposed() const
    {return MATRIX(x[0],x[2],x[1],x[3]);}

    T Trace() const
    {return x[0]+x[3];}

    MATRIX Cofactor_Matrix() const // cheap
    {return MATRIX(x[3],-x[2],-x[1],x[0]);}

    SYMMETRIC_MATRIX<T,2> Normal_Equations_Matrix() const // 6 mults, 3 adds
    {return SYMMETRIC_MATRIX<T,2>(x[0]*x[0]+x[1]*x[1],x[0]*x[2]+x[1]*x[3],x[2]*x[2]+x[3]*x[3]);}

    SYMMETRIC_MATRIX<T,2> Outer_Product_Matrix() const // 6 mults, 3 adds
    {return SYMMETRIC_MATRIX<T,2>(x[0]*x[0]+x[2]*x[2],x[0]*x[1]+x[2]*x[3],x[1]*x[1]+x[3]*x[3]);}

    SYMMETRIC_MATRIX<T,2> Symmetric_Part() const
    {return SYMMETRIC_MATRIX<T,2>(x[0],(T).5*(x[1]+x[2]),x[3]);}

    SYMMETRIC_MATRIX<T,2> Twice_Symmetric_Part() const
    {return SYMMETRIC_MATRIX<T,2>(2*x[0],x[1]+x[2],2*x[3]);}

    DIAGONAL_MATRIX<T,2> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,2>(x[0],x[3]);}

    static MATRIX Transpose(const MATRIX& A)
    {return MATRIX(A.x[0],A.x[2],A.x[1],A.x[3]);}

    static MATRIX Identity_Matrix()
    {return MATRIX(1,0,0,1);}

    static MATRIX Cross_Product_Matrix(const VECTOR<T,1>& v)
    {return MATRIX(0,v.x,-v.x,0);}

    T Antisymmetric_Part_Cross_Product_Vector() const
    {return (T).5*(x[1]-x[2]);}

    static MATRIX Rotation_Matrix(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(c,s,-s,c);}

    static MATRIX Derivative_Rotation_Matrix(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(-s,c,-c,-s);}

    static MATRIX Outer_Product(const VECTOR<T,2>& u,const VECTOR<T,2>& v)
    {return MATRIX(u.x*v.x,u.y*v.x,u.x*v.y,u.y*v.y);}

    static T Inner_Product(const MATRIX& A,const MATRIX& B)
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3];}

    T Frobenius_Norm_Squared() const
    {return sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3]);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    T Max_Abs() const
    {return maxabs(x[0],x[1],x[2],x[3]);}

    MATRIX operator*(const DIAGONAL_MATRIX<T,2>& A) const // 4 mults
    {return MATRIX(x[0]*A.x.x,x[1]*A.x.x,x[2]*A.x.y,x[3]*A.x.y);}

    MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 6 mults, 2 adds
    {return MATRIX(x[0]*A.x00,x[1]*A.x00,x[0]*A.x01+x[2]*A.x11,x[1]*A.x01+x[3]*A.x11);}

    MATRIX operator*(const SYMMETRIC_MATRIX<T,2>& A) const // 8 mults, 4 adds
    {return MATRIX(x[0]*A.x00+x[2]*A.x10,x[1]*A.x00+x[3]*A.x10,x[0]*A.x10+x[2]*A.x11,x[1]*A.x10+x[3]*A.x11);}

    MATRIX Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 6 mults, 2 adds
    {return MATRIX(x[0]*A.x00+x[2]*A.x01,x[1]*A.x00+x[3]*A.x01,x[2]*A.x11,x[3]*A.x11);}

    MATRIX Times_Transpose(const MATRIX& A) const // 8 mults, 4 adds
    {return MATRIX(x[0]*A.x[0]+x[2]*A.x[2],x[1]*A.x[0]+x[3]*A.x[2],x[0]*A.x[1]+x[2]*A.x[3],x[1]*A.x[1]+x[3]*A.x[3]);}

    MATRIX<T,2,3> Times_Transpose(const MATRIX<T,3,2>& A) const // 12 mults, 6 adds
    {return (A.Times_Transpose(*this)).Transposed();}

    MATRIX Transpose_Times(const MATRIX& A) const // 8 mults, 4 adds
    {return MATRIX(x[0]*A.x[0]+x[1]*A.x[1],x[2]*A.x[0]+x[3]*A.x[1],x[0]*A.x[2]+x[1]*A.x[3],x[2]*A.x[2]+x[3]*A.x[3]);}

    VECTOR<T,2> Transpose_Times(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x[0]*v.x+x[1]*v.y,x[2]*v.x+x[3]*v.y);}

    ZERO_VECTOR<T,2> Transpose_Times(const ZERO_VECTOR<T,2>& y) const
    {return ZERO_VECTOR<T,2>();}

    UPPER_TRIANGULAR_MATRIX<T,2> R_From_QR_Factorization() const
    {if(x[1]==0) return UPPER_TRIANGULAR_MATRIX<T,2>(x[0],x[2],x[3]);T c,s;
    if(abs(x[1])>abs(x[0])){T t=-x[0]/x[1];s=1/sqrt(1+t*t);c=s*t;}
    else{T t=-x[1]/x[0];c=1/sqrt(1+t*t);s=c*t;}
    return UPPER_TRIANGULAR_MATRIX<T,2>(c*x[0]-s*x[1],c*x[2]-s*x[3],s*x[2]+c*x[3]);}

    void Indefinite_Polar_Decomposition(MATRIX& Q,SYMMETRIC_MATRIX<T,2>& S) const
    {T x03=x[0]+x[3],cosine,sine;if(x03==0){cosine=0;sine=1;}else{T t=(x[1]-x[2])/x03;cosine=1/sqrt(1+t*t);sine=t*cosine;}
    Q=MATRIX(cosine,sine,-sine,cosine);S=SYMMETRIC_MATRIX<T,2>(Q.x[0]*x[0]+Q.x[1]*x[1],Q.x[0]*x[2]+Q.x[1]*x[3],Q.x[2]*x[2]+Q.x[3]*x[3]);}

    void Singular_Value_Decomposition(MATRIX& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX& V) const
    {MATRIX Q;SYMMETRIC_MATRIX<T,2> S;Indefinite_Polar_Decomposition(Q,S);S.Solve_Eigenproblem(singular_values,V);
    if(singular_values.x.y<0 && abs(singular_values.x.y)>=abs(singular_values.x.x)){
        singular_values=DIAGONAL_MATRIX<T,2>(-singular_values.x.y,-singular_values.x.x);
        Q=-Q;V=MATRIX(V.x[2],V.x[3],-V.x[0],-V.x[1]);}
    U=Q*V;}

    static T Determinant_Differential(const MATRIX& A,const MATRIX& dA)
    {return dA.x[0]*A.x[3]+A.x[0]*dA.x[3]-dA.x[1]*A.x[2]-A.x[1]*dA.x[2];}

    static MATRIX Cofactor_Differential(const MATRIX& dA)
    {return dA.Cofactor_Matrix();}

    T Simplex_Minimum_Altitude() const
    {return Determinant()/sqrt(max(Column(0).Magnitude_Squared(),Column(1).Magnitude_Squared(),(Column(0)-Column(1)).Magnitude_Squared()));}

    static MATRIX Componentwise_Min(const MATRIX& v1,const MATRIX& v2)
    {return MATRIX(min(v1.x[0],v2.x[0]),min(v1.x[1],v2.x[1]),min(v1.x[2],v2.x[2]),min(v1.x[3],v2.x[3]));}

    static MATRIX Componentwise_Max(const MATRIX& v1,const MATRIX& v2)
    {return MATRIX(max(v1.x[0],v2.x[0]),max(v1.x[1],v2.x[1]),max(v1.x[2],v2.x[2]),max(v1.x[3],v2.x[3]));}

    VECTOR<T,1> Contract_Permutation_Tensor() const
    {return VECTOR<T,1>(x[2]-x[1]);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,x,m*n);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,x,m*n);}

//#####################################################################
};
// global functions
template<class T>
inline MATRIX<T,2> operator+(const T a,const MATRIX<T,2>& A)
{return A+a;}

template<class T>
inline MATRIX<T,2> operator-(const T a,const MATRIX<T,2>& A)
{return -A+a;}

template<class T>
inline MATRIX<T,2> operator+(const SYMMETRIC_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x00+B.x[0],A.x10+B.x[1],A.x10+B.x[2],A.x11+B.x[3]);}

template<class T>
inline MATRIX<T,2> operator-(const SYMMETRIC_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x00-B.x[0],A.x10-B.x[1],A.x10-B.x[2],A.x11-B.x[3]);}

template<class T>
inline MATRIX<T,2> operator+(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x00+B.x[0],B.x[1],A.x01+B.x[2],A.x11+B.x[3]);}

template<class T>
inline MATRIX<T,2> operator-(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x00-B.x[0],-B.x[1],A.x01-B.x[2],A.x11-B.x[3]);}

template<class T>
inline MATRIX<T,2> operator*(const T a,const MATRIX<T,2>& A)
{return MATRIX<T,2>(a*A.x[0],a*A.x[1],a*A.x[2],a*A.x[3]);}

template<class T>
inline VECTOR<T,2> operator*(const VECTOR<T,2>& v,const MATRIX<T,2>& A)
{return VECTOR<T,2> (v.x*A.x[0]+v.y*A.x[2],v.x*A.x[1]+v.y*A.x[3]);}

template<class T>
inline MATRIX<T,2> operator*(const DIAGONAL_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x.x*B.x[0],A.x.y*B.x[1],A.x.x*B.x[2],A.x.y*B.x[3]);}

template<class T>
inline MATRIX<T,2> operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const MATRIX<T,2>& B)
{return MATRIX<T,2>(A.x00*B.x[0]+A.x01*B.x[1],A.x11*B.x[1],A.x00*B.x[2]+A.x01*B.x[3],A.x11*B.x[3]);}

template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,2>& A)
{Ignore(input,'[');for(int i=0;i<2;i++){for(int j=0;j<2;j++)input>>A.x[i+j*2];Ignore(input,';');}Ignore(input,']');return input;}
}
#endif
