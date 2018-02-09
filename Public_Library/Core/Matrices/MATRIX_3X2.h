//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X2
//#####################################################################
#ifndef __MATRIX_3X2__
#define __MATRIX_3X2__

#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX_BASE.h>
#include <Core/Vectors/VECTOR.h>
#include <Core/Vectors/ZERO_VECTOR.h>
namespace PhysBAM{

template<class T>
class MATRIX<T,3,2>:public MATRIX_BASE<T,MATRIX<T,3,2> >
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,3,2> > BASE;
    enum WORKAROUND1 {m=3,n=2};
    using BASE::operator*;using BASE::Times_Transpose;using BASE::Transpose_Times;

    T x[6];

    MATRIX()
    {
        for(int i=0;i<6;i++) x[i]=T();
    }

    MATRIX(INITIAL_SIZE mm,INITIAL_SIZE nn)
    {
        assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(2));
        for(int i=0;i<6;i++) x[i]=T();
    }

    MATRIX(const MATRIX& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=matrix_input.x[i];
    }

    template<class T2> explicit
    MATRIX(const MATRIX<T2,3,2>& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=(T)matrix_input.x[i];
    }

    MATRIX(const T x00,const T x10,const T x20,const T x01,const T x11,const T x21)
    {
        x[0]=x00;x[1]=x10;x[2]=x20;x[3]=x01;x[4]=x11;x[5]=x21;
    }

    MATRIX(const VECTOR<T,3>& column1,const VECTOR<T,3>& column2)
    {
        x[0]=column1.x;x[1]=column1.y;x[2]=column1.z;x[3]=column2.x;x[4]=column2.y;x[5]=column2.z;
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==3 && A.Columns()==2);
        for(int j=0;j<2;j++) for(int i=0;i<3;i++) (*this)(i,j)=A(i,j);
    }

    MATRIX& operator=(const MATRIX& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=matrix_input.x[i];
        return *this;
    }

    int Rows() const
    {return 3;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<3);assert((unsigned)j<2);return x[i+3*j];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<3);assert((unsigned)j<2);return x[i+3*j];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<3 && (unsigned)j<2;}

    VECTOR<T,3> Column(const int j) const
    {assert((unsigned)j<2);return VECTOR<T,3>(x[3*j],x[3*j+1],x[3*j+2]);}

    void Set_Column(const int j,const VECTOR<T,3>& v)
    {assert((unsigned)j<2);x[3*j]=v.x;x[3*j+1]=v.y;x[3*j+2]=v.z;}

    void Add_Column(const int j,const VECTOR<T,3>& v)
    {assert((unsigned)j<2);x[3*j]+=v.x;x[3*j+1]+=v.y;x[3*j+2]+=v.z;}

    VECTOR<T,2> Row(const int j) const
    {assert((unsigned)j<3);return VECTOR<T,2>(x[j],x[j+3]);}

    void Set_Row(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<3);x[j]=v.x;x[j+3]=v.y;}

    void Add_Row(const int j,const VECTOR<T,2>& v)
    {assert((unsigned)j<3);x[j]+=v.x;x[j+3]+=v.y;}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<6;i++) if(x[i]!=A.x[i]) return false;
    return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,3> Column_Sum() const
    {return VECTOR<T,3>(x[0]+x[3],x[1]+x[4],x[2]+x[5]);}

    MATRIX operator-() const
    {return MATRIX(-x[0],-x[1],-x[2],-x[3],-x[4],-x[5]);}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<6;i++) x[i]+=A.x[i];
    return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<6;i++) x[i]-=A.x[i];
    return *this;}

    MATRIX& operator*=(const MATRIX<T,2>& A)
    {return *this=*this*A;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<6;i++) x[i]*=a;
    return *this;}

    MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;
    for(int i=0;i<6;i++) x[i]*=s;
    return *this;}

    MATRIX operator+(const MATRIX& A) const // 6 adds
    {return MATRIX(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3],x[4]+A.x[4],x[5]+A.x[5]);}

    MATRIX operator-(const MATRIX& A) const // 6 adds
    {return MATRIX(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3],x[4]-A.x[4],x[5]-A.x[5]);}

    MATRIX operator*(const MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[1],x[1]*A.x[0]+x[4]*A.x[1],x[2]*A.x[0]+x[5]*A.x[1],x[0]*A.x[2]+x[3]*A.x[3],x[1]*A.x[2]+x[4]*A.x[3],x[2]*A.x[2]+x[5]*A.x[3]);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(2==A.m);MATRIX_MXN<T> matrix(3,A.n);
    for(int j=0;j<A.n;j++) for(int i=0;i<3;i++) for(int k=0;k<2;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);
    return matrix;}

    MATRIX Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 9 mults, 3 adds
    {return MATRIX(x[0]*A.x00+x[3]*A.x01,x[1]*A.x00+x[4]*A.x01,x[2]*A.x00+x[5]*A.x01,x[3]*A.x11,x[4]*A.x11,x[5]*A.x11);}

    MATRIX Times_Transpose(const MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[2],x[1]*A.x[0]+x[4]*A.x[2],x[2]*A.x[0]+x[5]*A.x[2],x[0]*A.x[1]+x[3]*A.x[3],x[1]*A.x[1]+x[4]*A.x[3],x[2]*A.x[1]+x[5]*A.x[3]);}

    MATRIX<T,3> Times_Transpose(const MATRIX& A) const // 18 mults, 9 adds
    {return MATRIX<T,3>(x[0]*A.x[0]+x[3]*A.x[3],x[1]*A.x[0]+x[4]*A.x[3],x[2]*A.x[0]+x[5]*A.x[3],
                        x[0]*A.x[1]+x[3]*A.x[4],x[1]*A.x[1]+x[4]*A.x[4],x[2]*A.x[1]+x[5]*A.x[4],
                        x[0]*A.x[2]+x[3]*A.x[5],x[1]*A.x[2]+x[4]*A.x[5],x[2]*A.x[2]+x[5]*A.x[5]);}

    MATRIX operator*(const T a) const // 6 mults
    {return MATRIX(a*x[0],a*x[1],a*x[2],a*x[3],a*x[4],a*x[5]);}

    MATRIX operator/(const T a) const // 6 mults, 1 div
    {assert(a!=0);T s=1/a;return MATRIX(s*x[0],s*x[1],s*x[2],s*x[3],s*x[4],s*x[5]);}

    VECTOR<T,3> operator*(const VECTOR<T,2>& v) const // 6 mults, 3 adds
    {return VECTOR<T,3>(x[0]*v.x+x[3]*v.y,x[1]*v.x+x[4]*v.y,x[2]*v.x+x[5]*v.y);}

    UPPER_TRIANGULAR_MATRIX<T,2> R_From_QR_Factorization() const // Gram Schmidt
    {T x_dot_x=Column(0).Magnitude_Squared(),x_dot_y=VECTOR<T,3>::Dot_Product(Column(0),Column(1)),y_dot_y=Column(1).Magnitude_Squared();
    T r00=sqrt(x_dot_x),r01=r00?x_dot_y/r00:0,r11=sqrt(max((T)0,y_dot_y-r01*r01));
    return UPPER_TRIANGULAR_MATRIX<T,2>(r00,r01,r11);}

    SYMMETRIC_MATRIX<T,2> Normal_Equations_Matrix() const // 9 mults, 6 adds
    {return SYMMETRIC_MATRIX<T,2>(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],x[3]*x[0]+x[4]*x[1]+x[5]*x[2],x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);}

    SYMMETRIC_MATRIX<T,3> Outer_Product_Matrix() const // 12 mults, 6 adds
    {return SYMMETRIC_MATRIX<T,3>(x[0]*x[0]+x[3]*x[3],x[0]*x[1]+x[3]*x[4],x[0]*x[2]+x[3]*x[5],x[1]*x[1]+x[4]*x[4],x[1]*x[2]+x[4]*x[5],x[2]*x[2]+x[5]*x[5]);}

    MATRIX<T,2> Transpose_Times(const MATRIX& A) const // 12 mults, 8 adds
    {return MATRIX<T,2>(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5]);}

    MATRIX<T,2,3> Transpose_Times(const MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],
                          x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5],
                          x[0]*A.x[6]+x[1]*A.x[7]+x[2]*A.x[8],x[3]*A.x[6]+x[4]*A.x[7]+x[5]*A.x[8]);}

    MATRIX<T,2,3> Transpose_Times(const SYMMETRIC_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x00+x[1]*A.x10+x[2]*A.x20,x[3]*A.x00+x[4]*A.x10+x[5]*A.x20,
                          x[0]*A.x10+x[1]*A.x11+x[2]*A.x21,x[3]*A.x10+x[4]*A.x11+x[5]*A.x21,
                          x[0]*A.x20+x[1]*A.x21+x[2]*A.x22,x[3]*A.x20+x[4]*A.x21+x[5]*A.x22);}

    MATRIX<T,2,3> Transpose_Times(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x00,x[3]*A.x00,x[0]*A.x01+x[1]*A.x11,x[3]*A.x01+x[4]*A.x11,
                          x[0]*A.x02+x[1]*A.x12+x[2]*A.x22,x[3]*A.x02+x[4]*A.x12+x[5]*A.x22);}

    MATRIX<T,2,3> Transpose_Times(const DIAGONAL_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x.x,x[3]*A.x.x,x[1]*A.x.y,x[4]*A.x.y,x[2]*A.x.z,x[5]*A.x.z);}

    VECTOR<T,2> Transpose_Times(const VECTOR<T,3>& v) const // 6 mults, 4 adds
    {return VECTOR<T,2>(x[0]*v.x+x[1]*v.y+x[2]*v.z,x[3]*v.x+x[4]*v.y+x[5]*v.z);}

    ZERO_VECTOR<T,2> Transpose_Times(const ZERO_VECTOR<T,3>& y) const
    {return ZERO_VECTOR<T,2>();}

    T Max_Abs() const
    {return maxabs(x[0],x[1],x[2],x[3],x[4],x[5]);}

    T Frobenius_Norm_Squared() const
    {return sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3])+sqr(x[4])+sqr(x[5]);}

    MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 9 mults, 3 adds
    {return MATRIX(x[0]*A.x00,x[1]*A.x00,x[2]*A.x00,x[0]*A.x01+x[3]*A.x11,x[1]*A.x01+x[4]*A.x11,x[2]*A.x01+x[5]*A.x11);}

    MATRIX operator*(const SYMMETRIC_MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x00+x[3]*A.x10,x[1]*A.x00+x[4]*A.x10,x[2]*A.x00+x[5]*A.x10,x[0]*A.x10+x[3]*A.x11,x[1]*A.x10+x[4]*A.x11,x[2]*A.x10+x[5]*A.x11);}

    MATRIX operator*(const DIAGONAL_MATRIX<T,2>& A) const // 6 mults
    {return MATRIX(x[0]*A.x.x,x[1]*A.x.x,x[2]*A.x.x,x[3]*A.x.y,x[4]*A.x.y,x[5]*A.x.y);}

    static T Inner_Product(const MATRIX& A,const MATRIX& B) // 6 mults, 5 adds
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5];}

    VECTOR<T,3> Weighted_Normal() const
    {return VECTOR<T,3>::Cross_Product(Column(0),Column(1));}

    MATRIX Cofactor_Matrix() const
    {VECTOR<T,3> normal=Weighted_Normal().Normalized();
    return MATRIX(VECTOR<T,3>::Cross_Product(Column(1),normal),VECTOR<T,3>::Cross_Product(normal,Column(0)));}

    T Parallelepiped_Measure() const
    {return Weighted_Normal().Magnitude();}

    MATRIX<T,2,3> Transposed() const
    {return MATRIX<T,2,3>::Transposed(*this);}

    static MATRIX Outer_Product(const VECTOR<T,3>& u,const VECTOR<T,2>& v)
    {MATRIX result;
    for(int i=0;i<3;i++) for(int j=0;j<2;j++) result(i,j)=u(i)*v(j);
    return result;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,x,m*n);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,x,m*n);}

//#####################################################################
    void Singular_Value_Decomposition(MATRIX& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX<T,2>& V) const;
    void Indefinite_Polar_Decomposition(MATRIX<T,3,2>& Q,SYMMETRIC_MATRIX<T,2>& S) const;
//#####################################################################
};
// global functions
template<class T>
inline MATRIX<T,3,2> operator*(const MATRIX<T,3>& B,const MATRIX<T,3,2>& A) // 18 mults, 12 adds
{return MATRIX<T,3,2>(B.x[0]*A.x[0]+B.x[3]*A.x[1]+B.x[6]*A.x[2],B.x[1]*A.x[0]+B.x[4]*A.x[1]+B.x[7]*A.x[2],B.x[2]*A.x[0]+B.x[5]*A.x[1]+B.x[8]*A.x[2],
    B.x[0]*A.x[3]+B.x[3]*A.x[4]+B.x[6]*A.x[5],B.x[1]*A.x[3]+B.x[4]*A.x[4]+B.x[7]*A.x[5],B.x[2]*A.x[3]+B.x[5]*A.x[4]+B.x[8]*A.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const T a,const MATRIX<T,3,2>& A) // 6 mults
{return MATRIX<T,3,2>(a*A.x[0],a*A.x[1],a*A.x[2],a*A.x[3],a*A.x[4],a*A.x[5]);}

template<class T>
inline VECTOR<T,3> operator*(const VECTOR<T,2>& v,const MATRIX<T,3,2>& A) // 6 mults, 3 adds
{return VECTOR<T,3>(v.x*A.x[0]+v.y*A.x[3],v.x*A.x[1]+v.y*A.x[4],v.x*A.x[2]+v.y*A.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const SYMMETRIC_MATRIX<T,3>& A,const MATRIX<T,3,2>& B) // 18 mults, 12 adds
{return MATRIX<T,3,2>(A.x00*B.x[0]+A.x10*B.x[1]+A.x20*B.x[2],A.x10*B.x[0]+A.x11*B.x[1]+A.x21*B.x[2],A.x20*B.x[0]+A.x21*B.x[1]+A.x22*B.x[2],
    A.x00*B.x[3]+A.x10*B.x[4]+A.x20*B.x[5],A.x10*B.x[3]+A.x11*B.x[4]+A.x21*B.x[5],A.x20*B.x[3]+A.x21*B.x[4]+A.x22*B.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const MATRIX<T,3,2>& B)
{return MATRIX<T,3,2>(A.x00*B.x[0]+A.x01*B.x[1]+A.x02*B.x[2],A.x11*B.x[1]+A.x12*B.x[2],A.x22*B.x[2],
    A.x00*B.x[3]+A.x01*B.x[4]+A.x02*B.x[5],A.x11*B.x[4]+A.x12*B.x[5],A.x22*B.x[5]);}

template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,3,2>& A)
{
    Ignore(input,'[');
    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++) input>>A.x[i+j*3];
        Ignore(input,';');}
    Ignore(input,']');
    return input;
}
}
#endif
