//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX
//#####################################################################
#ifndef __SYMMETRIC_MATRIX__
#define __SYMMETRIC_MATRIX__
#include <Tools/Matrices/SYMMETRIC_MATRIX_0X0.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_1X1.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
namespace PhysBAM{

template<class T,int d>
class SYMMETRIC_MATRIX
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {m=d,n=d,size=d*(d+1)/2};
    STATIC_ASSERT(d>3);

    T x[size];

    SYMMETRIC_MATRIX()
    {
        for(int i=0;i<size;i++) x[i]=T();
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,d>& matrix_input)
    {
        for(int i=0;i<size;i++) x[i]=(T)matrix_input.x[i];
    }

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,d>& matrix_input)
    {
        for(int i=0;i<size;i++) x[i]=T();
        for(int i=0;i<d;i++) Element_Lower(i,i)=matrix_input.x(i);
    }

    int Rows() const
    {return d;}

    int Columns() const
    {return d;}

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<d && (unsigned)j<d;}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert((unsigned)i<(unsigned)d && (unsigned)j<=(unsigned)i);return x[((2*d-j-1)*j>>1)+i];}

    const T& Element_Lower(int i,int j) const
    {assert((unsigned)i<(unsigned)d && (unsigned)j<=(unsigned)i);return x[((2*d-j-1)*j>>1)+i];}

    VECTOR<T,d> Column(const int axis) const
    {assert((unsigned)axis<(unsigned)d);VECTOR<T,d> v;
    for(int i=0;i<axis;i++) v(i)=Element_Lower(axis,i);
    for(int i=axis;i<d;i++) v(i)=Element_Lower(i,axis);
    return v;}

    VECTOR<T,d> Row(const int axis) const
    {return Column(axis);}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {for(int i=0;i<size;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return !(*this==A);}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=min(v1.x[i],v2.x[i]);return s;}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=max(v1.x[i],v2.x[i]);return s;}

    SYMMETRIC_MATRIX operator-() const
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=-x[i];return s;}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {for(int i=0;i<size;i++) x[i]+=A.x[i];return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {for(int i=0;i<d;i++) Element_Lower(i,i)+=a;return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {for(int i=0;i<size;i++) x[i]-=A.x[i];return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {for(int i=0;i<d;i++) Element_Lower(i,i)-=a;return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {for(int i=0;i<size;i++) x[i]*=a;return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {assert(a!=0);return *this*=1/a;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=x[i]+A.x[i];return s;}

    SYMMETRIC_MATRIX operator+(const T a) const
    {SYMMETRIC_MATRIX s(*this);s+=a;return s;}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=x[i]-A.x[i];return s;}

    SYMMETRIC_MATRIX operator-(const T a) const
    {SYMMETRIC_MATRIX s(*this);s-=a;return s;}

    SYMMETRIC_MATRIX operator*(const T a) const
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=x[i]*a;return s;}

    SYMMETRIC_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,d> operator*(const VECTOR<T,d>& v) const
    {VECTOR<T,d> w;for(int i=0;i<d;i++){
    T z=0;
    for(int j=0;j<i;j++) z+=Element_Lower(i,j)*v(j);
    for(int j=i;j<d;j++) z+=Element_Lower(j,i)*v(j);
    w(i)=z;} return w;}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    T Dilational() const
    {return Trace()/d;}

    SYMMETRIC_MATRIX Deviatoric() const
    {return *this-Dilational();}

    T Trace() const
    {T z=0;for(int i=0;i<d;i++) z+=Element_Lower(i,i);return z;}

    T Max_Abs() const
    {T z=x[0];for(int i=0;i<size;i++) z=max(z,abs(x[i]));return z;}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,d>& u)
    {SYMMETRIC_MATRIX s;for(int i=0;i<d;i++) for(int j=0;j<=i;j++) s.Element_Lower(i,j)=u(i)*u(j);return s;}

    static SYMMETRIC_MATRIX Symmetric_Outer_Product(const VECTOR<T,d>& u,const VECTOR<T,d>& v)
    {SYMMETRIC_MATRIX s;for(int i=0;i<d;i++) for(int j=0;j<=i;j++) s.Element_Lower(i,j)=u(i)*v(j)+v(i)*u(j);return s;}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX()+1;}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {SYMMETRIC_MATRIX s;for(int i=0;i<size;i++) s.x[i]=scale;return s;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,x,size);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,x,size);}
};
// global functions
template<class T,int d>
inline SYMMETRIC_MATRIX<T,d> operator*(const T a,const SYMMETRIC_MATRIX<T,d>& A)
{return A*a;}

template<class T,int d>
inline SYMMETRIC_MATRIX<T,d> operator+(const T a,const SYMMETRIC_MATRIX<T,d>& A)
{return A+a;}

template<class T,int d>
inline SYMMETRIC_MATRIX<T,d> operator-(const T a,const SYMMETRIC_MATRIX<T,d>& A)
{return -A+a;}

//#####################################################################
// Function Symmetric_Times_Transpose
//#####################################################################
template<class T,int m,int n> SYMMETRIC_MATRIX<T,m>
Symmetric_Times_Transpose(const MATRIX<T,m,n>& a,const MATRIX<T,m,n>& b)
{return a.Times_Transpose(b).Twice_Symmetric_Part();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Times_Transpose(const MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b)
{return a.Times_Transpose(b).Twice_Symmetric_Part();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const MATRIX<T,d>& b)
{return a.Times_Transpose(b).Twice_Symmetric_Part();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b)
{return a.Times_Transpose(b).Twice_Symmetric_Part();}

//#####################################################################
// Function Symmetric_Transpose_Times
//#####################################################################
template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const MATRIX<T,d>& b)
{return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b)
{return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int m,int n> SYMMETRIC_MATRIX<T,n>
Symmetric_Transpose_Times(const MATRIX<T,m,n>& a,const MATRIX<T,m,n>& b)
{return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Transpose_Times(const MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b)
{return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int d> VECTOR<T,d>
Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const VECTOR<T,d>& b)
{return a.Transpose_Times(b);}

//#####################################################################
// Function Times_Self_Transpose
//#####################################################################
template<class T,int m,int n> SYMMETRIC_MATRIX<T,m>
Times_Self_Transpose(const MATRIX<T,m,n>& a)
{return a.Outer_Product_Matrix();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Times_Self_Transpose(const SYMMETRIC_MATRIX<T,d>& a)
{return a*a;}

//#####################################################################
// Function Transpose_Times_Self
//#####################################################################
template<class T,int d> SYMMETRIC_MATRIX<T,d>
Transpose_Times_Self(const MATRIX<T,d>& a)
{return a.Normal_Equations_Matrix();}

template<class T,int d> SYMMETRIC_MATRIX<T,d>
Transpose_Times_Self(const SYMMETRIC_MATRIX<T,d>& a)
{return a*a;}

//#####################################################################
// Function Symmetric_Outer_Product
//#####################################################################
template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Outer_Product(const VECTOR<T,d>& u,const VECTOR<T,d>& v)
{return SYMMETRIC_MATRIX<T,d>::Symmetric_Outer_Product(u,v);}

//#####################################################################
// Function Outer_Product
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>
Outer_Product(const TV& u)
{return SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(u);}

template<class T,int d> inline std::ostream&
operator<<(std::ostream& o,const SYMMETRIC_MATRIX<T,d>& A)
{o<<"[";for(int i=0;i<d;i++){for(int j=0;j<d;j++){o<<A(i,j);if(j<d-1) o<<" ";}if(i<d-1) o<<"; ";}o<<"]";return o;}

template<class T,int d> struct SUM<SYMMETRIC_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> > {typedef SYMMETRIC_MATRIX<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<SYMMETRIC_MATRIX<T,d>,SYMMETRIC_MATRIX<T,d> > {typedef SYMMETRIC_MATRIX<T,d> TYPE;};
template<class T,int d> struct NEGATION<SYMMETRIC_MATRIX<T,d> > {typedef SYMMETRIC_MATRIX<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<SYMMETRIC_MATRIX<T,d>,T> {typedef SYMMETRIC_MATRIX<T,d> TYPE;};
}
#endif
