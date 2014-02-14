//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Silvia Salinas-Blemker, Craig Schroeder,
//     Andrew Selle, Tamar Shinar, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X3
//#####################################################################
#ifndef __MATRIX_3X3__
#define __MATRIX_3X3__

#include <Tools/Math_Tools/exchange.h>
#include <Tools/Matrices/MATRIX_BASE.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class MATRIX<T,3>:public MATRIX_BASE<T,MATRIX<T,3> >
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,3> > BASE;
    enum WORKAROUND1 {m=3,n=3};
    using BASE::operator*;using BASE::Times_Transpose;using BASE::Transpose_Times;

    T x[9];

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
    {
        assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));for(int i=0;i<9;i++) x[i]=T();
    }

    MATRIX(const MATRIX& matrix_input)
    {
        for(int i=0;i<9;i++) x[i]=matrix_input.x[i];
    }

    template<class T2> explicit
    MATRIX(const MATRIX<T2,3>& matrix_input)
    {
        for(int i=0;i<9;i++) x[i]=(T)matrix_input.x[i];
    }

    MATRIX(const DIAGONAL_MATRIX<T,3>& matrix_input)
    {
        x[0]=matrix_input.x.x;x[4]=matrix_input.x.y;x[8]=matrix_input.x.z;x[1]=x[2]=x[3]=x[5]=x[6]=x[7]=0;
    }

    MATRIX(const SYMMETRIC_MATRIX<T,3>& matrix_input)
    {
        x[0]=matrix_input.x00;x[1]=x[3]=matrix_input.x10;x[2]=x[6]=matrix_input.x20;x[4]=matrix_input.x11;x[5]=x[7]=matrix_input.x21;x[8]=matrix_input.x22;
    }

    MATRIX(const UPPER_TRIANGULAR_MATRIX<T,3>& matrix_input)
    {
        x[0]=matrix_input.x00;x[3]=matrix_input.x01;x[4]=matrix_input.x11;x[6]=matrix_input.x02;x[7]=matrix_input.x12;x[8]=matrix_input.x22;x[1]=x[2]=x[5]=0;
    }

    MATRIX(const T x00,const T x10,const T x20,const T x01,const T x11,const T x21,const T x02,const T x12,const T x22)
    {
        x[0]=x00;x[1]=x10;x[2]=x20;x[3]=x01;x[4]=x11;x[5]=x21;x[6]=x02;x[7]=x12;x[8]=x22;
    }

    MATRIX(const VECTOR<T,3>& column1,const VECTOR<T,3>& column2,const VECTOR<T,3>& column3)
    {
        x[0]=column1.x;x[1]=column1.y;x[2]=column1.z;x[3]=column2.x;x[4]=column2.y;x[5]=column2.z;x[6]=column3.x;x[7]=column3.y;x[8]=column3.z;
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==3 && A.Columns()==3);for(int j=0;j<3;j++) for(int i=0;i<3;i++) (*this)(i,j)=A(i,j);
    }

    MATRIX& operator=(const MATRIX& matrix_input)
    {
        for(int i=0;i<9;i++) x[i]=matrix_input.x[i];return *this;
    }

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<3);assert((unsigned)j<3);return x[i+3*j];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<3);assert((unsigned)j<3);return x[i+3*j];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<3 && (unsigned)j<3;}

    VECTOR<T,3> Column(const int j) const
    {assert((unsigned)j<3);return VECTOR<T,3>(x[3*j],x[3*j+1],x[3*j+2]);}

    void Set_Column(const int j,const VECTOR<T,3>& v)
    {assert((unsigned)j<3);x[3*j]=v.x;x[3*j+1]=v.y;x[3*j+2]=v.z;}

    VECTOR<T,3> Row(const int j) const
    {assert((unsigned)j<3);return VECTOR<T,3>(x[j],x[j+3],x[j+6]);}

    void Set_Row(const int j,const VECTOR<T,3>& v)
    {assert((unsigned)j<3);x[j]=v(0);x[j+3]=v(1);x[j+6]=v(2);}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<9;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,3> Column_Sum() const
    {return VECTOR<T,3>(x[0]+x[3]+x[6],x[1]+x[4]+x[7],x[2]+x[5]+x[8]);}

    MATRIX operator-() const
    {return MATRIX(-x[0],-x[1],-x[2],-x[3],-x[4],-x[5],-x[6],-x[7],-x[8]);}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<9;i++) x[i]+=A.x[i];return *this;}

    MATRIX& operator+=(const T& a)
    {x[0]+=a;x[4]+=a;x[8]+=a;return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<9;i++) x[i]-=A.x[i];return *this;}

    MATRIX& operator-=(const T& a)
    {x[0]-=a;x[4]-=a;x[8]-=a;return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {return *this=*this*A;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<9;i++) x[i]*=a;return *this;}

    MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;for(int i=0;i<9;i++) x[i]*=s;return *this;}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3],x[4]+A.x[4],x[5]+A.x[5],x[6]+A.x[6],x[7]+A.x[7],x[8]+A.x[8]);}

    MATRIX operator+(const T a) const
    {return MATRIX(x[0]+a,x[1],x[2],x[3],x[4]+a,x[5],x[6],x[7],x[8]+a);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3],x[4]-A.x[4],x[5]-A.x[5],x[6]-A.x[6],x[7]-A.x[7],x[8]-A.x[8]);}

    MATRIX operator-(const T a) const
    {return MATRIX(x[0]-a,x[1],x[2],x[3],x[4]-a,x[5],x[6],x[7],x[8]-a);}

    MATRIX operator*(const MATRIX& A) const // 27 mults, 18 adds
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[1]+x[6]*A.x[2],x[1]*A.x[0]+x[4]*A.x[1]+x[7]*A.x[2],x[2]*A.x[0]+x[5]*A.x[1]+x[8]*A.x[2],
        x[0]*A.x[3]+x[3]*A.x[4]+x[6]*A.x[5],x[1]*A.x[3]+x[4]*A.x[4]+x[7]*A.x[5],x[2]*A.x[3]+x[5]*A.x[4]+x[8]*A.x[5],
        x[0]*A.x[6]+x[3]*A.x[7]+x[6]*A.x[8],x[1]*A.x[6]+x[4]*A.x[7]+x[7]*A.x[8],x[2]*A.x[6]+x[5]*A.x[7]+x[8]*A.x[8]);}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x[0],a*x[1],a*x[2],a*x[3],a*x[4],a*x[5],a*x[6],a*x[7],a*x[8]);}

    MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return MATRIX(s*x[0],s*x[1],s*x[2],s*x[3],s*x[4],s*x[5],s*x[6],s*x[7],s*x[8]);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const // 9 mults, 6 adds
    {return VECTOR<T,3>(x[0]*v.x+x[3]*v.y+x[6]*v.z,x[1]*v.x+x[4]*v.y+x[7]*v.z,x[2]*v.x+x[5]*v.y+x[8]*v.z);}

    VECTOR<T,2> Homogeneous_Times(const VECTOR<T,2>& v) const // assumes w=1 is the 3rd coordinate of v
    {T w=x[2]*v.x+x[5]*v.y+x[8];assert(w!=0);
    if(w==1) return VECTOR<T,2>(x[0]*v.x+x[3]*v.y+x[6],x[1]*v.x+x[4]*v.y+x[7]);
    else{T s=1/w;return VECTOR<T,2>(s*(x[0]*v.x+x[3]*v.y+x[6]),s*(x[1]*v.x+x[4]*v.y+x[7]));}} // rescale so w=1

    VECTOR<T,2> Transform_2X2(const VECTOR<T,2>& v) const // multiplies vector by upper 2x2 of matrix only
    {return VECTOR<T,2>(x[0]*v.x+x[3]*v.y,x[1]*v.x+x[4]*v.y);}

    T Determinant() const // 9 mults, 5 adds
    {return x[0]*(x[4]*x[8]-x[7]*x[5])+x[3]*(x[7]*x[2]-x[1]*x[8])+x[6]*(x[1]*x[5]-x[4]*x[2]);}

    T Parallelepiped_Measure() const
    {return Determinant();}

    void Invert()
    {*this=Inverse();}

    MATRIX Inverse_Transposed() const
    {return Inverse().Transposed();}

    MATRIX Rotation_Only() const
    {return MATRIX(x[0],x[1],0,x[3],x[4],0,0,0,1);}

    VECTOR<T,2> Extract_Translation() const
    {return VECTOR<T,2>(x[6],x[7]);}

    void Set_Translation(const VECTOR<T,2>& v)
    {x[6]=v.x;x[7]=v.y;}

    MATRIX<T,2> Extract_Rotation() const
    {return MATRIX<T,2>(x[0],x[1],x[3],x[4]);}

    static MATRIX From_Linear(const MATRIX<T,2>& M) // Create a homogeneous 3x3 matrix corresponding to a 2x2 transform
    {return MATRIX(M.x[0],M.x[1],0,M.x[2],M.x[3],0,0,0,1);}

    void Transpose()
    {exchange(x[1],x[3]);exchange(x[2],x[6]);exchange(x[5],x[7]);}

    MATRIX Transposed() const
    {return MATRIX(x[0],x[3],x[6],x[1],x[4],x[7],x[2],x[5],x[8]);}

    T Trace() const
    {return x[0]+x[4]+x[8];}

    MATRIX Cofactor_Matrix() const // 18 mults, 9 adds
    {return MATRIX(x[4]*x[8]-x[5]*x[7],x[5]*x[6]-x[3]*x[8],x[3]*x[7]-x[4]*x[6],
                   x[2]*x[7]-x[1]*x[8],x[0]*x[8]-x[2]*x[6],x[1]*x[6]-x[0]*x[7],
                   x[1]*x[5]-x[2]*x[4],x[2]*x[3]-x[0]*x[5],x[0]*x[4]-x[1]*x[3]);}

    SYMMETRIC_MATRIX<T,3> Outer_Product_Matrix() const // 18 mults, 12 adds
    {return SYMMETRIC_MATRIX<T,3>(x[0]*x[0]+x[3]*x[3]+x[6]*x[6],x[1]*x[0]+x[4]*x[3]+x[7]*x[6],x[2]*x[0]+x[5]*x[3]+x[8]*x[6],
                                  x[1]*x[1]+x[4]*x[4]+x[7]*x[7],x[2]*x[1]+x[5]*x[4]+x[8]*x[7],x[2]*x[2]+x[5]*x[5]+x[8]*x[8]);}

    SYMMETRIC_MATRIX<T,3> Normal_Equations_Matrix() const // 18 mults, 12 adds
    {return SYMMETRIC_MATRIX<T,3>(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],x[3]*x[0]+x[4]*x[1]+x[5]*x[2],x[6]*x[0]+x[7]*x[1]+x[8]*x[2],
                                  x[3]*x[3]+x[4]*x[4]+x[5]*x[5],x[6]*x[3]+x[7]*x[4]+x[8]*x[5],x[6]*x[6]+x[7]*x[7]+x[8]*x[8]);}

    SYMMETRIC_MATRIX<T,3> Symmetric_Part() const // 3 mults, 3 adds
    {return SYMMETRIC_MATRIX<T,3>(x[0],(T).5*(x[1]+x[3]),(T).5*(x[2]+x[6]),x[4],(T).5*(x[5]+x[7]),x[8]);}

    SYMMETRIC_MATRIX<T,3> Twice_Symmetric_Part() const // 3 mults, 3 adds
    {return SYMMETRIC_MATRIX<T,3>(2*x[0],x[1]+x[3],x[2]+x[6],2*x[4],x[5]+x[7],2*x[8]);}

    DIAGONAL_MATRIX<T,3> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,3>(x[0],x[4],x[8]);}

    void Normalize_Columns()
    {T magnitude=sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2]));assert(magnitude!=0);T s=1/magnitude;x[0]*=s;x[1]*=s;x[2]*=s;
    magnitude=sqrt(sqr(x[3])+sqr(x[4])+sqr(x[5]));assert(magnitude!=0);s=1/magnitude;x[3]*=s;x[4]*=s;x[5]*=s;
    magnitude=sqrt(sqr(x[6])+sqr(x[7])+sqr(x[8]));assert(magnitude!=0);s=1/magnitude;x[6]*=s;x[7]*=s;x[8]*=s;}

    VECTOR<T,3> Column_Magnitudes() const
    {return VECTOR<T,3>(Column(0).Magnitude(),Column(1).Magnitude(),Column(2).Magnitude());}

    static MATRIX Transpose(const MATRIX& A)
    {return MATRIX(A.x[0],A.x[3],A.x[6],A.x[1],A.x[4],A.x[7],A.x[2],A.x[5],A.x[8]);}

    static MATRIX Translation_Matrix(const VECTOR<T,2>& translation) // treating the 3x3 matrix as a homogeneous transformation on 2d vectors
    {return MATRIX(1,0,0,0,1,0,translation.x,translation.y,1);}

    static MATRIX Identity_Matrix()
    {return MATRIX(1,0,0,0,1,0,0,0,1);}

    static MATRIX Rotation_Matrix_X_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(1,0,0,0,c,s,0,-s,c);}

    static MATRIX Rotation_Matrix_Y_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(c,0,-s,0,1,0,s,0,c);}

    static MATRIX Rotation_Matrix_Z_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(c,s,0,-s,c,0,0,0,1);}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& axis,const T radians)
    {return ROTATION<VECTOR<T,3> >(radians,axis).Rotation_Matrix();}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& rotation)
    {return ROTATION<VECTOR<T,3> >::From_Rotation_Vector(rotation).Rotation_Matrix();}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& x_final,const VECTOR<T,3>& y_final,const VECTOR<T,3>& z_final)
    {return MATRIX(x_final,y_final,z_final);}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& initial_vector,const VECTOR<T,3>& final_vector)
    {return ROTATION<VECTOR<T,3> >::From_Rotated_Vector(initial_vector,final_vector).Rotation_Matrix();}

    static MATRIX Scale_Matrix(const VECTOR<T,2>& scale_vector)
    {return MATRIX(scale_vector.x,0,0,0,scale_vector.y,0,0,0,1);}

    static MATRIX Scale_Matrix(const T scale)
    {return MATRIX(scale,0,0,0,scale,0,0,0,1);}

    static MATRIX Outer_Product(const VECTOR<T,3>& u,const VECTOR<T,3>& v)
    {return MATRIX(u.x*v.x,u.y*v.x,u.z*v.x,u.x*v.y,u.y*v.y,u.z*v.y,u.x*v.z,u.y*v.z,u.z*v.z);}

    static T Inner_Product(const MATRIX& A,const MATRIX& B)
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5]+A.x[6]*B.x[6]+A.x[7]*B.x[7]+A.x[8]*B.x[8];}

    VECTOR<T,3> Antisymmetric_Part_Cross_Product_Vector() const
    {return (T).5*VECTOR<T,3>(x[5]-x[7],x[6]-x[2],x[1]-x[3]);}

    static SYMMETRIC_MATRIX<T,3> Right_Multiply_With_Symmetric_Result(const MATRIX& A,const DIAGONAL_MATRIX<T,3>& B)
    {return SYMMETRIC_MATRIX<T,3>(B.x.x*A.x[0],B.x.x*A.x[1],B.x.x*A.x[2],B.x.y*A.x[4],B.x.y*A.x[5],B.x.z*A.x[8]);}

    T Max_Abs() const
    {return maxabs(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);}

    T Frobenius_Norm_Squared() const
    {return sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3])+sqr(x[4])+sqr(x[5])+sqr(x[6])+sqr(x[7])+sqr(x[8]);}

    MATRIX operator*(const DIAGONAL_MATRIX<T,3>& A) const // 9 mults
    {return MATRIX(x[0]*A.x.x,x[1]*A.x.x,x[2]*A.x.x,x[3]*A.x.y,x[4]*A.x.y,x[5]*A.x.y,x[6]*A.x.z,x[7]*A.x.z,x[8]*A.x.z);}

    MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const // 18 mults, 9 adds
    {return MATRIX(x[0]*A.x00,x[1]*A.x00,x[2]*A.x00,x[0]*A.x01+x[3]*A.x11,x[1]*A.x01+x[4]*A.x11,x[2]*A.x01+x[5]*A.x11,
                          x[0]*A.x02+x[3]*A.x12+x[6]*A.x22,x[1]*A.x02+x[4]*A.x12+x[7]*A.x22,x[2]*A.x02+x[5]*A.x12+x[8]*A.x22);}

    MATRIX operator*(const SYMMETRIC_MATRIX<T,3>& A) const // 27 mults, 18 adds
    {return MATRIX(x[0]*A.x00+x[3]*A.x10+x[6]*A.x20,x[1]*A.x00+x[4]*A.x10+x[7]*A.x20,x[2]*A.x00+x[5]*A.x10+x[8]*A.x20,
                          x[0]*A.x10+x[3]*A.x11+x[6]*A.x21,x[1]*A.x10+x[4]*A.x11+x[7]*A.x21,x[2]*A.x10+x[5]*A.x11+x[8]*A.x21,
                          x[0]*A.x20+x[3]*A.x21+x[6]*A.x22,x[1]*A.x20+x[4]*A.x21+x[7]*A.x22,x[2]*A.x20+x[5]*A.x21+x[8]*A.x22);}

    MATRIX Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const
    {return MATRIX(x[0]*A.x00+x[3]*A.x01+x[6]*A.x02,x[1]*A.x00+x[4]*A.x01+x[7]*A.x02,x[2]*A.x00+x[5]*A.x01+x[8]*A.x02,
                          x[3]*A.x11+x[6]*A.x12,x[4]*A.x11+x[7]*A.x12,x[5]*A.x11+x[8]*A.x12,x[6]*A.x22,x[7]*A.x22,x[8]*A.x22);}

    MATRIX Transpose_Times(const MATRIX& A) const
    {return MATRIX(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],x[6]*A.x[0]+x[7]*A.x[1]+x[8]*A.x[2],
                          x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5],x[6]*A.x[3]+x[7]*A.x[4]+x[8]*A.x[5],
                          x[0]*A.x[6]+x[1]*A.x[7]+x[2]*A.x[8],x[3]*A.x[6]+x[4]*A.x[7]+x[5]*A.x[8],x[6]*A.x[6]+x[7]*A.x[7]+x[8]*A.x[8]);}

    VECTOR<T,3> Transpose_Times(const VECTOR<T,3>& v) const
    {return VECTOR<T,3>(x[0]*v.x+x[1]*v.y+x[2]*v.z,x[3]*v.x+x[4]*v.y+x[5]*v.z,x[6]*v.x+x[7]*v.y+x[8]*v.z);}

    MATRIX Times_Transpose(const MATRIX& A) const
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[3]+x[6]*A.x[6],x[1]*A.x[0]+x[4]*A.x[3]+x[7]*A.x[6],x[2]*A.x[0]+x[5]*A.x[3]+x[8]*A.x[6],
                   x[0]*A.x[1]+x[3]*A.x[4]+x[6]*A.x[7],x[1]*A.x[1]+x[4]*A.x[4]+x[7]*A.x[7],x[2]*A.x[1]+x[5]*A.x[4]+x[8]*A.x[7],
                   x[0]*A.x[2]+x[3]*A.x[5]+x[6]*A.x[8],x[1]*A.x[2]+x[4]*A.x[5]+x[7]*A.x[8],x[2]*A.x[2]+x[5]*A.x[5]+x[8]*A.x[8]);}

    MATRIX Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const // (v*) * (*this)
    {return MATRIX(-v.z*x[1]+v.y*x[2],v.z*x[0]-v.x*x[2],-v.y*x[0]+v.x*x[1],
                   -v.z*x[4]+v.y*x[5],v.z*x[3]-v.x*x[5],-v.y*x[3]+v.x*x[4],
                   -v.z*x[7]+v.y*x[8],v.z*x[6]-v.x*x[8],-v.y*x[6]+v.x*x[7]);}

    MATRIX Cross_Product_Matrix_Transpose_Times(const VECTOR<T,3>& v) const // (v*)^T * (*this)
    {return Cross_Product_Matrix_Times(-v);}

    MATRIX Times_Cross_Product_Matrix(const VECTOR<T,3>& v) const // (*this) * (v*)
    {return MATRIX( x[3]*v.z-x[6]*v.y, x[4]*v.z-x[7]*v.y, x[5]*v.z-x[8]*v.y,
                   -x[0]*v.z+x[6]*v.x,-x[1]*v.z+x[7]*v.x,-x[2]*v.z+x[8]*v.x,
                    x[0]*v.y-x[3]*v.x, x[1]*v.y-x[4]*v.x, x[2]*v.y-x[5]*v.x);}

    MATRIX Times_Cross_Product_Matrix_Transpose(const VECTOR<T,3>& v) const // (*this) * (v*)^T
    {return Times_Cross_Product_Matrix(-v);}

    SYMMETRIC_MATRIX<T,3> Cross_Product_Matrix_Times_With_Symmetric_Result(const VECTOR<T,3>& v) const // (v*) * (*this)
    {return SYMMETRIC_MATRIX<T,3>(-v.z*x[1]+v.y*x[2],v.z*x[0]-v.x*x[2],-v.y*x[0]+v.x*x[1],v.z*x[3]-v.x*x[5],-v.y*x[3]+v.x*x[4],-v.y*x[6]+v.x*x[7]);}

    SYMMETRIC_MATRIX<T,3> Times_Cross_Product_Matrix_With_Symmetric_Result(const VECTOR<T,3>& v) const // (*this) * (v*)
    {return SYMMETRIC_MATRIX<T,3>(x[3]*v.z-x[6]*v.y,x[4]*v.z-x[7]*v.y,x[5]*v.z-x[8]*v.y,-x[1]*v.z+x[7]*v.x,-x[2]*v.z+x[8]*v.x,x[2]*v.y-x[5]*v.x);}

    SYMMETRIC_MATRIX<T,3> Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(const VECTOR<T,3>& v) const // (*this) * (v*)^T
    {return Times_Cross_Product_Matrix_With_Symmetric_Result(-v);}

    static MATRIX Left_Procrustes_Rotation(const MATRIX& A,const MATRIX& B)
    {MATRIX U,V;DIAGONAL_MATRIX<T,3> D;A.Times_Transpose(B).Fast_Singular_Value_Decomposition(U,D,V);return U.Times_Transpose(V);}

    static MATRIX<T,3> Cross_Product_Matrix(const VECTOR<T,3>& v)
    {return MATRIX<T,3>(0,v.z,-v.y,-v.z,0,v.x,v.y,-v.x,0);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,x,m*n);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,x,m*n);}

//#####################################################################
    MATRIX Higham_Iterate(const T tolerance=1e-5,const int max_iterations=20,const bool exit_on_max_iterations=false) const;
    void Fast_Singular_Value_Decomposition(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& singular_values,MATRIX<T,3>& V) const;
    void Fast_Indefinite_Polar_Decomposition(MATRIX<T,3>& Q,SYMMETRIC_MATRIX<T,3>& S) const;
    T Simplex_Minimum_Altitude() const;
    static MATRIX Componentwise_Min(const MATRIX& v1,const MATRIX& v2);
    static MATRIX Componentwise_Max(const MATRIX& v1,const MATRIX& v2);
    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const;
    MATRIX Inverse() const;
    VECTOR<T,3> Solve_Linear_System(const VECTOR<T,3>& b) const; // 33 mults, 17 adds, 1 div
    VECTOR<T,3> Robust_Solve_Linear_System(const VECTOR<T,3>& b) const; // 34 mults, 17 adds, 1 div
    MATRIX Q_From_QR_Factorization() const; // Gram Schmidt
    UPPER_TRIANGULAR_MATRIX<T,3> R_From_QR_Factorization() const; // Gram Schmidt
//#####################################################################
};
// global functions
template<class T>
inline MATRIX<T,3> operator+(const T a,const MATRIX<T,3>& A)
{return A+a;}

template<class T>
inline MATRIX<T,3> operator*(const T a,const MATRIX<T,3>& A)
{return A*a;}

template<class T>
inline MATRIX<T,3> operator-(const T a,const MATRIX<T,3>& A)
{return MATRIX<T,3>(a-A.x[0],-A.x[1],-A.x[2],-A.x[3],a-A.x[4],-A.x[5],-A.x[6],-A.x[7],a-A.x[8]);}

template<class T>
inline VECTOR<T,3> operator*(const VECTOR<T,3>& v,const MATRIX<T,3>& A)
{return VECTOR<T,3>(v.x*A.x[0]+v.y*A.x[1]+v.z*A.x[2],v.x*A.x[3]+v.y*A.x[4]+v.z*A.x[5],v.x*A.x[6]+v.y*A.x[7]+v.z*A.x[8]);}

template<class T>
inline MATRIX<T,3> operator*(const DIAGONAL_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x.x*B.x[0],A.x.y*B.x[1],A.x.z*B.x[2],A.x.x*B.x[3],A.x.y*B.x[4],A.x.z*B.x[5],A.x.x*B.x[6],A.x.y*B.x[7],A.x.z*B.x[8]);}

template<class T>
inline MATRIX<T,3> operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00*B.x[0]+A.x01*B.x[1]+A.x02*B.x[2],A.x11*B.x[1]+A.x12*B.x[2],A.x22*B.x[2],A.x00*B.x[3]+A.x01*B.x[4]+A.x02*B.x[5],
                      A.x11*B.x[4]+A.x12*B.x[5],A.x22*B.x[5],A.x00*B.x[6]+A.x01*B.x[7]+A.x02*B.x[8],A.x11*B.x[7]+A.x12*B.x[8],A.x22*B.x[8]);}

template<class T>
inline MATRIX<T,3> operator*(const SYMMETRIC_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00*B.x[0]+A.x10*B.x[1]+A.x20*B.x[2],A.x10*B.x[0]+A.x11*B.x[1]+A.x21*B.x[2],A.x20*B.x[0]+A.x21*B.x[1]+A.x22*B.x[2],
                      A.x00*B.x[3]+A.x10*B.x[4]+A.x20*B.x[5],A.x10*B.x[3]+A.x11*B.x[4]+A.x21*B.x[5],A.x20*B.x[3]+A.x21*B.x[4]+A.x22*B.x[5],
                      A.x00*B.x[6]+A.x10*B.x[7]+A.x20*B.x[8],A.x10*B.x[6]+A.x11*B.x[7]+A.x21*B.x[8],A.x20*B.x[6]+A.x21*B.x[7]+A.x22*B.x[8]);}

template<class T>
inline MATRIX<T,3> operator+(const SYMMETRIC_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00+B.x[0],A.x10+B.x[1],A.x20+B.x[2],A.x10+B.x[3],A.x11+B.x[4],A.x21+B.x[5],A.x20+B.x[6],A.x21+B.x[7],A.x22+B.x[8]);}

template<class T>
inline MATRIX<T,3> operator-(const SYMMETRIC_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00-B.x[0],A.x10-B.x[1],A.x20-B.x[2],A.x10-B.x[3],A.x11-B.x[4],A.x21-B.x[5],A.x20-B.x[6],A.x21-B.x[7],A.x22-B.x[8]);}

template<class T>
inline MATRIX<T,3> operator+(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00+B.x[0],B.x[1],B.x[2],A.x01+B.x[3],A.x11+B.x[4],B.x[5],A.x02+B.x[6],A.x12+B.x[7],A.x22+B.x[8]);}

template<class T>
inline MATRIX<T,3> operator-(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const MATRIX<T,3>& B)
{return MATRIX<T,3>(A.x00-B.x[0],-B.x[1],-B.x[2],A.x01-B.x[3],A.x11-B.x[4],-B.x[5],A.x02-B.x[6],A.x12-B.x[7],A.x22-B.x[8]);}

template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,3>& A)
{FILE_UTILITIES::Ignore(input,'[');
    for(int i=0;i<3;i++){for(int j=0;j<3;j++) input>>A.x[i+j*3];FILE_UTILITIES::Ignore(input,';');}
FILE_UTILITIES::Ignore(input,']');
return input;}
//#####################################################################
}
#endif

