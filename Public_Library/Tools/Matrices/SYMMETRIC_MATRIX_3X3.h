//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_3X3
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_3X3__
#define __SYMMETRIC_MATRIX_3X3__

#include <Tools/Math_Tools/Robust_Arithmetic.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T> struct HAS_CHEAP_COPY<SYMMETRIC_MATRIX<T,3> > {static const bool value=true;};
template<class T> struct is_scalar_BLOCK<SYMMETRIC_MATRIX<T,3> >:public is_scalar_BLOCK<T>{};
template<class T> struct is_scalar_VECTOR_SPACE<SYMMETRIC_MATRIX<T,3> >:public is_scalar_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,3>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class SYMMETRIC_MATRIX<T,3>
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x00,x10,x20,x11,x21,x22;

    SYMMETRIC_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
        :x00(T()),x10(T()),x20(T()),x11(T()),x21(T()),x22(T())
    {
        STATIC_ASSERT(sizeof(SYMMETRIC_MATRIX)==6*sizeof(T));assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,3>& matrix_input)
        :x00((T)matrix_input.x00),x10((T)matrix_input.x10),x20((T)matrix_input.x20),x11((T)matrix_input.x11),x21((T)matrix_input.x21),x22((T)matrix_input.x22)
    {}

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,3>& matrix_input)
        :x00(matrix_input.x.x),x10(T()),x20(T()),x11(matrix_input.x.y),x21(T()),x22(matrix_input.x.z)
    {}

    SYMMETRIC_MATRIX(const T y00,const T y10,const T y20,const T y11,const T y21,const T y22)
        :x00(y00),x10(y10),x20(y20),x11(y11),x21(y21),x22(y22)
    {}

    void From_Matrix(const MATRIX<T,3>& matrix_input)
    {
        x00=matrix_input(0,0);
        x10=matrix_input(1,0);
        x20=matrix_input(2,0);
        x11=matrix_input(1,1);
        x21=matrix_input(2,1);
        x22=matrix_input(2,2);
    }

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<3 && (unsigned)j<3;}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert((unsigned)i<(unsigned)3 && (unsigned)j<=(unsigned)i);return ((T*)this)[((5-j)*j>>1)+i];}

    const T& Element_Lower(int i,int j) const
    {assert((unsigned)i<(unsigned)3 && (unsigned)j<=(unsigned)i);return ((const T*)this)[((5-j)*j>>1)+i];}

    VECTOR<T,3> Column(const int axis) const
    {assert((unsigned)axis<(unsigned)3);return axis==0?VECTOR<T,3>(x00,x10,x20):axis==1?VECTOR<T,3>(x10,x11,x21):VECTOR<T,3>(x20,x21,x22);}

    VECTOR<T,3> Row(const int axis) const
    {return Column(axis);}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {return x00==A.x00 && x10==A.x10 && x20==A.x20 && x11==A.x11 && x21==A.x21 && x22==A.x22;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return !(*this==A);}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(min(v1.x00,v2.x00),min(v1.x10,v2.x10),min(v1.x20,v2.x20),min(v1.x11,v2.x11),min(v1.x21,v2.x21),min(v1.x22,v2.x22));}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(max(v1.x00,v2.x00),max(v1.x10,v2.x10),max(v1.x20,v2.x20),max(v1.x11,v2.x11),max(v1.x21,v2.x21),max(v1.x22,v2.x22));}

    SYMMETRIC_MATRIX operator-() const
    {return SYMMETRIC_MATRIX(-x00,-x10,-x20,-x11,-x21,-x22);}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {x00+=A.x00;x10+=A.x10;x20+=A.x20;x11+=A.x11;x21+=A.x21;x22+=A.x22;return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {x00+=a;x11+=a;x22+=a;return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {x00-=A.x00;x10-=A.x10;x20-=A.x20;x11-=A.x11;x21-=A.x21;x22-=A.x22;return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {x00-=a;x11-=a;x22-=a;return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {x00*=a;x10*=a;x20*=a;x11*=a;x21*=a;x22*=a;return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x00*=s;x10*=s;x20*=s;x11*=s;x21*=s;x22*=s;return *this;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x00+A.x00,x10+A.x10,x20+A.x20,x11+A.x11,x21+A.x21,x22+A.x22);}

    SYMMETRIC_MATRIX operator+(const T a) const
    {return SYMMETRIC_MATRIX(x00+a,x10,x20,x11+a,x21,x22+a);}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x00-A.x00,x10-A.x10,x20-A.x20,x11-A.x11,x21-A.x21,x22-A.x22);}

    SYMMETRIC_MATRIX operator-(const T a) const
    {return SYMMETRIC_MATRIX(x00-a,x10,x20,x11-a,x21,x22-a);}

    SYMMETRIC_MATRIX operator*(const T a) const
    {return SYMMETRIC_MATRIX(a*x00,a*x10,a*x20,a*x11,a*x21,a*x22);}

    MATRIX<T,3> operator*(const SYMMETRIC_MATRIX& A) const // 27 mults, 18 adds
    {return MATRIX<T,3>(x00*A.x00+x10*A.x10+x20*A.x20,x10*A.x00+x11*A.x10+x21*A.x20,x20*A.x00+x21*A.x10+x22*A.x20,
                        x00*A.x10+x10*A.x11+x20*A.x21,x10*A.x10+x11*A.x11+x21*A.x21,x20*A.x10+x21*A.x11+x22*A.x21,
                        x00*A.x20+x10*A.x21+x20*A.x22,x10*A.x20+x11*A.x21+x21*A.x22,x20*A.x20+x21*A.x21+x22*A.x22);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(3==A.m);MATRIX_MXN<T> matrix(3,A.n);for(int j=0;j<A.n;j++) matrix.Set_Column(j,(*this)*VECTOR<T,3>(A(0,j),A(1,j),A(2,j)));return matrix;}

    template<class T_MATRIX>
    typename PRODUCT<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)A.Columns(),(INITIAL_SIZE)A.Rows());A.Add_Times_Transpose(*this,A.Derived(),M);return M;}

    MATRIX<T,3> Times_Transpose(const SYMMETRIC_MATRIX& M) const // 27 mults, 18 adds
    {return *this*M;}

    MATRIX<T,3> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& M) const
    {return MATRIX<T,3>(x00*M.x00+x10*M.x01+x20*M.x02,x10*M.x00+x11*M.x01+x21*M.x02,x20*M.x00+x21*M.x01+x22*M.x02,x10*M.x11+x20*M.x12,x11*M.x11+x21*M.x12,x21*M.x11+x22*M.x12,x20*M.x22,x21*M.x22,x22*M.x22);}

    MATRIX<T,3> Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const // (v*) * (*this)
    {return MATRIX<T,3>(-v.z*x10+v.y*x20,v.z*x00-v.x*x20,-v.y*x00+v.x*x10,-v.z*x11+v.y*x21,v.z*x10-v.x*x21,-v.y*x10+v.x*x11,-v.z*x21+v.y*x22,v.z*x20-v.x*x22,-v.y*x20+v.x*x21);}

    MATRIX<T,3> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,3>& v) const // (v*)^T * (*this)
    {return MATRIX<T,3>(v.z*x10-v.y*x20,-v.z*x00+v.x*x20,v.y*x00-v.x*x10,v.z*x11-v.y*x21,-v.z*x10+v.x*x21,v.y*x10-v.x*x11,v.z*x21-v.y*x22,-v.z*x20+v.x*x22,v.y*x20-v.x*x21);}

    MATRIX<T,3> Times_Cross_Product_Matrix(const VECTOR<T,3>& v) const // (*this) * (v*)
    {return MATRIX<T,3>(x10*v.z-x20*v.y,x11*v.z-x21*v.y,x21*v.z-x22*v.y,-x00*v.z+x20*v.x,-x10*v.z+x21*v.x,-x20*v.z+x22*v.x,x00*v.y-x10*v.x,x10*v.y-x11*v.x,x20*v.y-x21*v.x);}

    SYMMETRIC_MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return SYMMETRIC_MATRIX(s*x00,s*x10,s*x20,s*x11,s*x21,s*x22);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const
    {return VECTOR<T,3>(x00*v.x+x10*v.y+x20*v.z,x10*v.x+x11*v.y+x21*v.z,x20*v.x+x21*v.y+x22*v.z);}

    T Determinant() const
    {return x00*(x11*x22-x21*x21)+x10*(2*x21*x20-x10*x22)-x20*x11*x20;}

    SYMMETRIC_MATRIX Inverse() const
    {T cofactor00=x11*x22-x21*x21,cofactor01=x21*x20-x10*x22,cofactor02=x10*x21-x11*x20;
    return SYMMETRIC_MATRIX(cofactor00,cofactor01,cofactor02,x00*x22-x20*x20,x10*x20-x00*x21,x00*x11-x10*x10)/(x00*cofactor00+x10*cofactor01+x20*cofactor02);}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    T Dilational() const
    {return ((T)1/3)*Trace();}

    SYMMETRIC_MATRIX Deviatoric() const
    {return *this-Dilational();}

    VECTOR<T,3> Inverse_Times(const VECTOR<T,3>& b) const // 18 mults, 8 adds
    {T cofactor00=x11*x22-x21*x21,cofactor01=x21*x20-x10*x22,cofactor02=x10*x21-x11*x20;
    return SYMMETRIC_MATRIX(cofactor00,cofactor01,cofactor02,x00*x22-x20*x20,x10*x20-x00*x21,x00*x11-x10*x10)*b/(x00*cofactor00+x10*cofactor01+x20*cofactor02);}

    VECTOR<T,3> Robust_Inverse_Times(const VECTOR<T,3>& b) const
    {T cofactor00=x11*x22-x21*x21,cofactor01=x21*x20-x10*x22,cofactor02=x10*x21-x11*x20;
    T determinant=x00*cofactor00+x10*cofactor01+x20*cofactor02;
    VECTOR<T,3> unscaled_result=SYMMETRIC_MATRIX(cofactor00,cofactor01,cofactor02,x00*x22-x20*x20,x10*x20-x00*x21,x00*x11-x10*x10)*b;
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;}

    SYMMETRIC_MATRIX Squared() const
    {return SYMMETRIC_MATRIX(x00*x00+x10*x10+x20*x20,x10*x00+x11*x10+x21*x20,x20*x00+x21*x10+x22*x20,x10*x10+x11*x11+x21*x21,x20*x10+x21*x11+x22*x21,x20*x20+x21*x21+x22*x22);}

    T Trace() const
    {return x00+x11+x22;}

    static T Inner_Product(const SYMMETRIC_MATRIX& A,const SYMMETRIC_MATRIX& B)
    {return A.x00*B.x00+A.x11*B.x11+A.x22*B.x22+2*(A.x10*B.x10+A.x20*B.x20+A.x21*B.x21);}

    T Frobenius_Norm_Squared() const
    {return x00*x00+x11*x11+x22*x22+2*(x10*x10+x20*x20+x21*x21);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    SYMMETRIC_MATRIX Cofactor_Matrix() const// 12 mults, 6 adds
    {return SYMMETRIC_MATRIX(x11*x22-x21*x21,x21*x20-x10*x22,x10*x21-x11*x20,x00*x22-x20*x20,x10*x20-x00*x21,x00*x11-x10*x10);}

    VECTOR<T,3> Largest_Column() const
    {T sqr00=sqr(x00),sqr01=sqr(x10),sqr02=sqr(x20),sqr11=sqr(x11),sqr12=sqr(x21),sqr22=sqr(x22);
    T scale1=sqr00+sqr01+sqr02,scale2=sqr01+sqr11+sqr12,scale3=sqr02+sqr12+sqr22;
    return scale1>scale2?(scale1>scale3?VECTOR<T,3>(x00,x10,x20):VECTOR<T,3>(x20,x21,x22)):(scale2>scale3?VECTOR<T,3>(x10,x11,x21):VECTOR<T,3>(x20,x21,x22));}

    VECTOR<T,3> Largest_Column_Normalized() const // 9 mults, 6 adds, 1 div, 1 sqrt
    {T sqr00=sqr(x00),sqr01=sqr(x10),sqr02=sqr(x20),sqr11=sqr(x11),sqr12=sqr(x21),sqr22=sqr(x22);
    T scale1=sqr00+sqr01+sqr02,scale2=sqr01+sqr11+sqr12,scale3=sqr02+sqr12+sqr22;
    if(scale1>scale2){if(scale1>scale3) return VECTOR<T,3>(x00,x10,x20)/sqrt(scale1);}
    else if(scale2>scale3) return VECTOR<T,3>(x10,x11,x21)/sqrt(scale2);
    if(scale3>0) return VECTOR<T,3>(x20,x21,x22)/sqrt(scale3);else return VECTOR<T,3>(1,0,0);}

    T Max_Abs() const
    {return maxabs(x00,x10,x20,x11,x21,x22);}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,3>& u) // 6 mults
    {return SYMMETRIC_MATRIX(u.x*u.x,u.x*u.y,u.x*u.z,u.y*u.y,u.y*u.z,u.z*u.z);}

    static SYMMETRIC_MATRIX Symmetric_Outer_Product(const VECTOR<T,3>& u,const VECTOR<T,3>& v)
    {return SYMMETRIC_MATRIX(2*v.x*u.x,v.x*u.y+v.y*u.x,v.x*u.z+v.z*u.x,2*v.y*u.y,v.y*u.z+v.z*u.y,2*v.z*u.z);}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX(1,0,0,1,0,1);}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {return SYMMETRIC_MATRIX(scale,scale,scale,scale,scale,scale);}

    bool Positive_Definite() const
    {return x00>0 && x00*x11>x10*x10 && Determinant()>0;}

    bool Positive_Semidefinite(const T tolerance=(T)1e-7) const
    {T scale=Max_Abs();return !scale || (*this+tolerance*scale).Positive_Definite();}

    VECTOR<T,3> First_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,3>& eigenvalues,const T tolerance=1e-5) const
    {T scale=maxabs(eigenvalues.x.x,eigenvalues.x.z),scale_inverse=Robust_Inverse(scale),tiny=tolerance*scale;
    if(eigenvalues.x.x-eigenvalues.x.y>tiny) return ((*this-eigenvalues.x.x)*scale_inverse).Cofactor_Matrix().Largest_Column_Normalized();
    return ((*this-eigenvalues.x.z)*scale_inverse).Cofactor_Matrix().Largest_Column().Unit_Orthogonal_Vector();}

    VECTOR<T,3> Last_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,3>& eigenvalues,const T tolerance=1e-5) const
    {T scale=maxabs(eigenvalues.x.x,eigenvalues.x.z),scale_inverse=Robust_Inverse(scale),tiny=tolerance*scale;
    if(eigenvalues.x.y-eigenvalues.x.z>tiny) return ((*this-eigenvalues.x.z)*scale_inverse).Cofactor_Matrix().Largest_Column_Normalized();
    return ((*this-eigenvalues.x.x)*scale_inverse).Cofactor_Matrix().Largest_Column().Unit_Orthogonal_Vector();}

    SYMMETRIC_MATRIX Positive_Definite_Part() const
    {DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> V;Fast_Solve_Eigenproblem(D,V);D=D.Clamp_Min(0);return Conjugate(V,D);}

    DIAGONAL_MATRIX<T,3> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,3>(x00,x11,x22);}

    VECTOR<T,3> Off_Diagonal_Part() const
    {return VECTOR<T,3>(x21,x20,x10);}

    static SYMMETRIC_MATRIX Multiply_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A*B and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[1]+A.x[6]*B.x[2],A.x[1]*B.x[0]+A.x[4]*B.x[1]+A.x[7]*B.x[2],
                             A.x[2]*B.x[0]+A.x[5]*B.x[1]+A.x[8]*B.x[2],A.x[1]*B.x[3]+A.x[4]*B.x[4]+A.x[7]*B.x[5],
                             A.x[2]*B.x[3]+A.x[5]*B.x[4]+A.x[8]*B.x[5],A.x[2]*B.x[6]+A.x[5]*B.x[7]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Times_Transpose_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A*B^t and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[3]+A.x[6]*B.x[6],A.x[1]*B.x[0]+A.x[4]*B.x[3]+A.x[7]*B.x[6],
                             A.x[2]*B.x[0]+A.x[5]*B.x[3]+A.x[8]*B.x[6],A.x[1]*B.x[1]+A.x[4]*B.x[4]+A.x[7]*B.x[7],
                             A.x[2]*B.x[1]+A.x[5]*B.x[4]+A.x[8]*B.x[7],A.x[2]*B.x[2]+A.x[5]*B.x[5]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Times_Transpose_With_Symmetric_Result(const MATRIX<T,3,2>& A,const MATRIX<T,3,2>& B) // A*B^t and assume symmetric result, 12 mults, 6 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[3],A.x[1]*B.x[0]+A.x[4]*B.x[3],A.x[2]*B.x[0]+A.x[5]*B.x[3],
                             A.x[1]*B.x[1]+A.x[4]*B.x[4],A.x[2]*B.x[1]+A.x[5]*B.x[4],A.x[2]*B.x[2]+A.x[5]*B.x[5]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A^t*B and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2],A.x[3]*B.x[0]+A.x[4]*B.x[1]+A.x[5]*B.x[2],A.x[6]*B.x[0]+A.x[7]*B.x[1]+A.x[8]*B.x[2],
                             A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5],A.x[6]*B.x[3]+A.x[7]*B.x[4]+A.x[8]*B.x[5],A.x[6]*B.x[6]+A.x[7]*B.x[7]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B) // A^t*B and assume symmetric result, 10 mults, 4 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x00,A.x[3]*B.x00,A.x[6]*B.x00,A.x[3]*B.x01+A.x[4]*B.x11,A.x[6]*B.x01+A.x[7]*B.x11,A.x[6]*B.x02+A.x[7]*B.x12+A.x[8]*B.x22);}

    SYMMETRIC_MATRIX<T,3> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {return Cross_Product_Matrix_Times(v).Times_Cross_Product_Matrix_With_Symmetric_Result(-v);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x00,x10,x20,x11,x21,x22);}

    template<class RW>  void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x00,x10,x20,x11,x21,x22);}

//#####################################################################
    MATRIX<T,3> operator*(const DIAGONAL_MATRIX<T,3>& A) const;
    MATRIX<T,3> operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const;
    SYMMETRIC_MATRIX operator+(const DIAGONAL_MATRIX<T,3>& A) const;
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const SYMMETRIC_MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3,2>& A,const DIAGONAL_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    DIAGONAL_MATRIX<T,3> Fast_Eigenvalues() const; // lambda_x > lambda_y > lambda_z
    void Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const;
    void Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const;
private:
    static void Jacobi_Transform(const int sweep,const T threshold,T& app,T& apq,T& aqq,T& arp,T& arq,T& v1p,T& v1q,T& v2p,T& v2q,T& v3p,T& v3q);
//#####################################################################
};
// global functions
template<class T>
inline SYMMETRIC_MATRIX<T,3> operator*(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return A*a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> operator+(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> operator-(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return -A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmin,const SYMMETRIC_MATRIX<T,3>& xmax)
{return SYMMETRIC_MATRIX<T,3>(clamp(x.x00,xmin.x00,xmax.x00),clamp(x.x10,xmin.x10,xmax.x10),clamp(x.x20,xmin.x20,xmax.x20),clamp(x.x11,xmin.x11,xmax.x11),clamp(x.x21,xmin.x21,xmax.x21),clamp(x.x22,xmin.x22,xmax.x22));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp_min(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmin)
{return SYMMETRIC_MATRIX<T,3>(clamp_min(x.x00,xmin.x00),clamp_min(x.x10,xmin.x10),clamp_min(x.x20,xmin.x20),clamp_min(x.x11,xmin.x11),clamp_min(x.x21,xmin.x21),clamp_min(x.x22,xmin.x22));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp_max(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmax)
{return SYMMETRIC_MATRIX<T,3>(clamp_max(x.x00,xmax.x00),clamp_max(x.x10,xmax.x10),clamp_max(x.x20,xmax.x20),clamp_max(x.x11,xmax.x11),clamp_max(x.x21,xmax.x21),clamp_max(x.x22,xmax.x22));}

template<class T>
inline std::ostream& operator<< (std::ostream& output_stream,const SYMMETRIC_MATRIX<T,3>& A)
{output_stream<<"["<<A.x00<<" "<<A.x10<<" "<<A.x20<<" ; "<<A.x10<<" "<<A.x11<<" "<<A.x21<<" ; "<<A.x20<<" "<<A.x21<<" "<<A.x22<<"]";return output_stream;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> log(const SYMMETRIC_MATRIX<T,3>& A)
{DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> Q;A.Fast_Solve_Eigenproblem(D,Q);return A.Conjugate(Q,log(D));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> exp(const SYMMETRIC_MATRIX<T,3>& A)
{DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> Q;A.Fast_Solve_Eigenproblem(D,Q);return A.Conjugate(Q,exp(D));}

//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B) // 27 mults, 12 adds
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const SYMMETRIC_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return SYMMETRIC_MATRIX(A.x.x*A.x.x*B.x00,A.x.y*A.x.x*B.x10,A.x.z*A.x.x*B.x20,A.x.y*A.x.y*B.x11,A.x.z*A.x.y*B.x21,A.x.z*A.x.z*B.x22);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3,2>& A,const DIAGONAL_MATRIX<T,2>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,2>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator*(const DIAGONAL_MATRIX<T,3>& A) const // 9 mults
{
    return MATRIX<T,3>(x00*A.x.x,x10*A.x.x,x20*A.x.x,x10*A.x.y,x11*A.x.y,x21*A.x.y,x20*A.x.z,x21*A.x.z,x22*A.x.z);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const // 18 mults, 9 adds
{
    return MATRIX<T,3>(x00*A.x00,x10*A.x00,x20*A.x00,x00*A.x01+x10*A.x11,x10*A.x01+x11*A.x11,x20*A.x01+x21*A.x11,
                       x00*A.x02+x10*A.x12+x20*A.x22,x10*A.x02+x11*A.x12+x21*A.x22,x20*A.x02+x21*A.x12+x22*A.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3>
operator*(const DIAGONAL_MATRIX<T,3>& D,const SYMMETRIC_MATRIX<T,3>& A) // 9 mults, 
{
    return MATRIX<T,3>(D.x.x*A.x00,D.x.y*A.x10,D.x.z*A.x20,D.x.x*A.x10,D.x.y*A.x11,D.x.z*A.x21,D.x.x*A.x20,D.x.y*A.x21,D.x.z*A.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3>
operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 18 mults, 9 adds
{
    return MATRIX<T,3>(A.x00*B.x00+A.x01*B.x10+A.x02*B.x20,A.x11*B.x10+A.x12*B.x20,A.x22*B.x20,A.x00*B.x10+A.x01*B.x11+A.x02*B.x21,
                       A.x11*B.x11+A.x12*B.x21,A.x22*B.x21,A.x00*B.x20+A.x01*B.x21+A.x02*B.x22,A.x11*B.x21+A.x12*B.x22,A.x22*B.x22);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator+(const DIAGONAL_MATRIX<T,3>& A) const // 3 adds
{
    return SYMMETRIC_MATRIX<T,3>(x00+A.x.x,x10,x20,x11+A.x.y,x21,x22+A.x.z);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3>
operator+(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 3 adds
{
    return B+A;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,3>
operator+(const SYMMETRIC_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{
    return MATRIX<T,3>(A.x00+B.x00,A.x10,A.x20,A.x10+B.x01,A.x11+B.x11,A.x21,A.x20+B.x02,A.x21+B.x12,A.x22+B.x22);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,3>
operator+(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,3>
operator-(const SYMMETRIC_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{
    return MATRIX<T,3>(A.x00-B.x00,A.x10,A.x20,A.x10-B.x01,A.x11-B.x11,A.x21,A.x20-B.x02,A.x21-B.x12,A.x22-B.x22);
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,3>
operator-(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return -B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3>
operator-(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 3 adds
{
    return SYMMETRIC_MATRIX<T,3>(A.x.x-B.x00,-B.x10,-B.x20,A.x.y-B.x11,-B.x21,A.x.z-B.x22);
}
//#####################################################################
}
#endif
