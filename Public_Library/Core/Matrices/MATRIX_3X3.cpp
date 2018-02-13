//#####################################################################
// Copyright 2002-2007, Silvia Salinas-Blemker, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Craig Schroeder, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X3
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/min.h>
#include <Core/Math_Tools/minabs.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Function Higham_Iterate
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Higham_Iterate(const T tolerance,const int max_iterations,const bool exit_on_max_iterations) const
{
    MATRIX<T,3> X=*this;int iterations=0;
    for(;;){
        MATRIX<T,3> Y=(T).5*(X+X.Inverse_Transposed());
        if((X-Y).Max_Abs()<tolerance) return Y;
        X=Y;
        if(++iterations>=max_iterations){
            if(exit_on_max_iterations) PHYSBAM_FATAL_ERROR();
            return X;}}
}
namespace{
template<class T>
VECTOR<T,2> Givens(T a,T b)
{
    return VECTOR<T,2>(a,-b).Normalized();
}
// Apply Givens rotation to matrix (on the left)
template<class T>
void Apply_L(MATRIX<T,3>& A,VECTOR<T,2> u,int i,int j)
{
    auto v=A.Row(i),w=A.Row(j);
    A.Set_Row(i,u.x*v+u.y*w);
    A.Set_Row(j,u.x*w-u.y*v);
}
// Apply Givens rotation to matrix (on the right)
template<class T>
void Apply_R(MATRIX<T,3>& A,VECTOR<T,2> u,int i,int j)
{
    auto v=A.Column(i),w=A.Column(j);
    A.Set_Column(i,u.x*v+u.y*w);
    A.Set_Column(j,u.x*w-u.y*v);
}
// Puts a zero in A(r1,c)
template<class T>
void Zero_L(MATRIX<T,3>& U,MATRIX<T,3>& A,int r0,int r1,int c)
{
    auto g=Givens(A(r0,c),A(r1,c));
    Apply_L(A,g,r1,r0);
    Apply_R(U,g,r1,r0);
    A(r1,c)=0;
}
// Puts a zero in A(r,c1)
template<class T>
void Zero_R(MATRIX<T,3>& V,MATRIX<T,3>& A,int c0,int c1,int r)
{
    auto g=Givens(A(r,c0),A(r,c1));
    Apply_R(A,g,c1,c0);
    Apply_R(V,g,c1,c0);
    A(r,c1)=0;
}
template<class T>
void Zero_Chasing(MATRIX<T,3>& U,MATRIX<T,3>& A,MATRIX<T,3>& V)
{
    // Initial sparsity:  [ * * * ; * * * ; 0 * * ]
    Zero_L(U,A,0,1,0); // [ * * * ; 0 * * ; 0 * * ]
    Zero_R(V,A,1,2,0); // [ * * 0 ; 0 * * ; 0 * * ]
    Zero_L(U,A,1,2,1); // [ * * 0 ; 0 * * ; 0 0 * ]
}
template<class T>
void Swap_Col(MATRIX<T,3>& U,int i,int j)
{
    auto u=U.Column(i),v=U.Column(j);
    U.Set_Column(i,v);
    U.Set_Column(j,-u);
}
template<class T>
void Swap_Row(MATRIX<T,3>& U,int i,int j)
{
    auto u=U.Row(i),v=U.Row(j);
    U.Set_Row(i,v);
    U.Set_Row(j,-u);
}
template<class T>
void Swap(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& D,MATRIX<T,3>& V,int i,int j)
{
    Swap_Col(U,i,j);
    Swap_Col(V,i,j);
    std::swap(D.x(i),D.x(j));
}
template<class T>
void Flip_Sign(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& D,int i,int j)
{
    D.x(i)=-D.x(i);
    D.x(j)=-D.x(j);
    U.Set_Column(i,-U.Column(i));
    U.Set_Column(j,-U.Column(j));
}
template<class T>
void Sort_SV(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& D,MATRIX<T,3>& V)
{
    int i=abs(D.x).Arg_Max();
    if(i>0) Swap(U,D,V,0,i);
    if(abs(D.x.y)<abs(D.x.z)) Swap(U,D,V,1,2);
    if(D.x.x<0) Flip_Sign(U,D,0,D.x.y<0?1:2);
    else if(D.x.y<0) Flip_Sign(U,D,1,2);
}
template<class T>
void Top_Block(MATRIX<T,3>& U,MATRIX<T,3>& B,MATRIX<T,3>& V,DIAGONAL_MATRIX<T,3>& D)
{
    MATRIX<T,2> u,v,a(B(0,0),B(1,0),B(0,1),B(1,1));
    DIAGONAL_MATRIX<T,2> d;
    a.Singular_Value_Decomposition(u,d,v);
    Apply_R(U,u.Column(0),0,1);
    Apply_R(V,v.Column(0),0,1);
    D.x.x=d.x.x;
    D.x.y=d.x.y;
    D.x.z=B(2,2);
    Sort_SV(U,D,V);
}
template<class T>
void Bot_Block(MATRIX<T,3>& U,MATRIX<T,3>& B,MATRIX<T,3>& V,DIAGONAL_MATRIX<T,3>& D)
{
    MATRIX<T,2> u,v,a(B(1,1),B(2,1),B(1,2),B(2,2));
    DIAGONAL_MATRIX<T,2> d;
    a.Singular_Value_Decomposition(u,d,v);
    Apply_R(U,u.Column(0),1,2);
    Apply_R(V,v.Column(0),1,2);
    D.x.x=B(0,0);
    D.x.y=d.x.x;
    D.x.z=d.x.y;
    Sort_SV(U,D,V);
}
template<class T>
bool Check_Tolerances(MATRIX<T,3>& B,MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& D,MATRIX<T,3>& V,T tau)
{
    if(abs(B(1,2))<=tau){
        Top_Block(U,B,V,D);
        return true;}
    if(abs(B(0,1))<=tau){
        Bot_Block(U,B,V,D);
        return true;}
    if(abs(B(1,1))<=tau){
        auto g=Givens(B(1,2),B(2,2)).Orthogonal_Vector();
        Apply_L(B,g,2,1);
        Apply_R(U,g,2,1);
        Top_Block(U,B,V,D);
        return true;}
    if(abs(B(2,2))<=tau){
        Zero_R(V,B,1,2,1);
        Zero_R(V,B,0,2,0);
        Top_Block(U,B,V,D);
        return true;}
    if(abs(B(0,0))<=tau){
        Zero_L(U,B,1,0,1);
        Zero_L(U,B,2,0,2);
        Bot_Block(U,B,V,D);
        return true;}
    return false;
}

template<class T> bool
Lazy_Sort(MATRIX<T,3>& U,MATRIX<T,3>& B,MATRIX<T,3>& V,int i,int j)
{
    T flip_ratio=100;
    if(abs(B(i,i))*flip_ratio>=abs(B(j,j))) return false;
    Swap_Row(B,i,j);
    Swap_Col(U,i,j);
    Swap_Col(B,i,j);
    Swap_Col(V,i,j);
    return true;
}

template<class T> void
Bidiagonalize(MATRIX<T,3>& U,MATRIX<T,3>& B,MATRIX<T,3>& V)
{
    Zero_L(U,B,1,2,0); // [ * * * ; * * * ; 0 * * ]
    Zero_Chasing(U,B,V); // [ * * 0 ; 0 * * ; 0 0 * ]
}
}
//#####################################################################
// Function Singular_Value_Decomposition
//#####################################################################
// U and V rotations, smallest singular value possibly negative
// @techreport{gast:2016:qrsvd,
//   title={Implicit-shifted Symmetric QR Singular Value Decomposition of 3x3 Matrices},
//   author={Gast, T. and Fu, C. and Jiang, C. and Teran, J.},
//   year={2016},
//   institution={University of California Los Angeles}
// }
template<class T> int MATRIX<T,3>::
Singular_Value_Decomposition(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& D,MATRIX<T,3>& V) const
{
    MATRIX<T,3> B=*this;
    V=U=MATRIX<T,3>()+1;
    Bidiagonalize(U,B,V);

    T nu=16*std::numeric_limits<T>::epsilon();
    T tau=nu*std::max((T).5*B.Frobenius_Norm(),(T)1);
    for(int it=0;;it++)
    {
        if(Check_Tolerances(B,U,D,V,tau)) return it;
        if(Lazy_Sort(U,B,V,0,1) | Lazy_Sort(U,B,V,0,2) | Lazy_Sort(U,B,V,1,2)){
            Bidiagonalize(U,B,V);
            if(Check_Tolerances(B,U,D,V,tau)) return it;}

        T al0=B(0,0),al1=B(1,1),al2=B(2,2);
        T be0=B(0,1),be1=B(1,2);
        T ga0=al0*be0,ga1=al1*be1;
        T a0=al1*al1+be0*be0;
        T a1=al2*al2+be1*be1;
        T b0=ga1;
        T d=(T).5*(a0-a1);
        T mu=a1-b0*b0/(d+sign_nonzero(d)*sqrt(d*d+b0*b0));

        auto g=Givens(al0*al0-mu,ga0);
        Apply_R(B,g,1,0);
        Apply_R(V,g,1,0);
        Zero_Chasing(U,B,V);
    }
}
//#####################################################################
// Function Indefinite_Polar_Decomposition
//#####################################################################
template<class T> void MATRIX<T,3>::
Indefinite_Polar_Decomposition(MATRIX<T,3>& Q,SYMMETRIC_MATRIX<T,3>& S) const
{
    MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> D;Singular_Value_Decomposition(U,D,V);
    Q=U.Times_Transpose(V);S=SYMMETRIC_MATRIX<T,3>::Conjugate(V,D);
}
//#####################################################################
// Function Simplex_Minimum_Altitude
//#####################################################################
template<class T> T MATRIX<T,3>::
Simplex_Minimum_Altitude() const
{
    typedef VECTOR<T,3> TV;
    TV X0=Column(0),X1=Column(1),X2=Column(2);
    return minabs(
        TV::Dot_Product(X0,TV::Cross_Product(X1-X0,X2-X0).Normalized()),
        TV::Dot_Product(X1-X0,TV::Cross_Product(X2,X1).Normalized()),
        TV::Dot_Product(X2-X1,TV::Cross_Product(X0,X2).Normalized()),
        TV::Dot_Product(X2,TV::Cross_Product(X0,X1).Normalized()));
}
//#####################################################################
// Function Componentwise_Min
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Componentwise_Min(const MATRIX& v1,const MATRIX& v2)
{
    return MATRIX(min(v1.x[0],v2.x[0]),min(v1.x[1],v2.x[1]),min(v1.x[2],v2.x[2]),min(v1.x[3],v2.x[3]),min(v1.x[4],v2.x[4]),min(v1.x[5],v2.x[5]),
        min(v1.x[6],v2.x[6]),min(v1.x[7],v2.x[7]),min(v1.x[8],v2.x[8]));
}
//#####################################################################
// Function Componentwise_Max
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Componentwise_Max(const MATRIX& v1,const MATRIX& v2)
{
    return MATRIX(max(v1.x[0],v2.x[0]),max(v1.x[1],v2.x[1]),max(v1.x[2],v2.x[2]),max(v1.x[3],v2.x[3]),max(v1.x[4],v2.x[4]),max(v1.x[5],v2.x[5]),
        max(v1.x[6],v2.x[6]),max(v1.x[7],v2.x[7]),max(v1.x[8],v2.x[8]));
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> MATRIX_MXN<T> MATRIX<T,3>::
operator*(const MATRIX_MXN<T>& A) const
{
    assert(3==A.m);
    MATRIX_MXN<T> matrix(3,A.n);
    for(int j=0;j<A.n;j++) for(int i=0;i<3;i++) for(int k=0;k<3;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);
    return matrix;
}
//#####################################################################
// Function Inverse
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Inverse() const
{
    T cofactor00=x[4]*x[8]-x[7]*x[5],cofactor01=x[7]*x[2]-x[1]*x[8],cofactor02=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor00+x[3]*cofactor01+x[6]*cofactor02;
    assert(determinant!=0);
    T s=1/determinant;
    return s*MATRIX(cofactor00,cofactor01,cofactor02,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],x[0]*x[4]-x[3]*x[1]);
}
//#####################################################################
// Function Inverse_Times
//#####################################################################
template<class T> VECTOR<T,3> MATRIX<T,3>::
Inverse_Times(const VECTOR<T,3>& b) const // 33 mults, 17 adds, 1 div
{
    T cofactor00=x[4]*x[8]-x[7]*x[5],cofactor01=x[7]*x[2]-x[1]*x[8],cofactor02=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor00+x[3]*cofactor01+x[6]*cofactor02;
    assert(determinant!=0);
    return MATRIX(cofactor00,cofactor01,cofactor02,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],x[0]*x[4]-x[3]*x[1])*b/determinant;
}
//#####################################################################
// Function Robust_Inverse_Times
//#####################################################################
template<class T> VECTOR<T,3> MATRIX<T,3>::
Robust_Inverse_Times(const VECTOR<T,3>& b) const // 34 mults, 17 adds, 1 div
{
    T cofactor00=x[4]*x[8]-x[7]*x[5],cofactor01=x[7]*x[2]-x[1]*x[8],cofactor02=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor00+x[3]*cofactor01+x[6]*cofactor02;
    VECTOR<T,3> unscaled_result=MATRIX(cofactor00,cofactor01,cofactor02,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],
        x[0]*x[4]-x[3]*x[1])*b;
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;
}
//#####################################################################
// Function Q_From_QR_Factorization
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Q_From_QR_Factorization() const // Gram Schmidt
{
    int k;MATRIX Q=*this;
    T one_over_r00=1/sqrt((sqr(Q.x[0])+sqr(Q.x[1])+sqr(Q.x[2])));
    for(k=0;k<=2;k++) Q.x[k]=one_over_r00*Q.x[k];
    T r01=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];
    Q.x[3]-=r01*Q.x[0];
    Q.x[4]-=r01*Q.x[1];
    Q.x[5]-=r01*Q.x[2];
    T r02=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];
    Q.x[6]-=r02*Q.x[0];
    Q.x[7]-=r02*Q.x[1];
    Q.x[8]-=r02*Q.x[2];
    T one_over_r11=1/sqrt((sqr(Q.x[3])+sqr(Q.x[4])+sqr(Q.x[5])));
    for(k=3;k<=5;k++) Q.x[k]=one_over_r11*Q.x[k];
    T r12=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];
    Q.x[6]-=r12*Q.x[3];
    Q.x[7]-=r12*Q.x[4];
    Q.x[8]-=r12*Q.x[5];
    T one_over_r22=1/sqrt((sqr(Q.x[6])+sqr(Q.x[7])+sqr(Q.x[8])));
    for(k=6;k<=8;k++) Q.x[k]=one_over_r22*Q.x[k];
    return Q;
}
//#####################################################################
// Function R_From_QR_Factorization
//#####################################################################
template<class T> UPPER_TRIANGULAR_MATRIX<T,3> MATRIX<T,3>::
R_From_QR_Factorization() const // Gram Schmidt
{
    int k;
    MATRIX Q=*this;
    UPPER_TRIANGULAR_MATRIX<T,3> R;
    R.x00=sqrt((sqr(Q.x[0])+sqr(Q.x[1])+sqr(Q.x[2])));
    T one_over_r00=1/R.x00;
    for(k=0;k<=2;k++) Q.x[k]=one_over_r00*Q.x[k];
    R.x01=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];
    Q.x[3]-=R.x01*Q.x[0];
    Q.x[4]-=R.x01*Q.x[1];
    Q.x[5]-=R.x01*Q.x[2];
    R.x02=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];
    Q.x[6]-=R.x02*Q.x[0];
    Q.x[7]-=R.x02*Q.x[1];
    Q.x[8]-=R.x02*Q.x[2];
    R.x11=sqrt((sqr(Q.x[3])+sqr(Q.x[4])+sqr(Q.x[5])));
    T one_over_r11=1/R.x11;
    for(k=3;k<=5;k++) Q.x[k]=one_over_r11*Q.x[k];
    R.x12=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];
    Q.x[6]-=R.x12*Q.x[3];
    Q.x[7]-=R.x12*Q.x[4];
    Q.x[8]-=R.x12*Q.x[5];
    R.x22=sqrt((sqr(Q.x[6])+sqr(Q.x[7])+sqr(Q.x[8])));
    return R;
}
//#####################################################################
namespace PhysBAM{
template class MATRIX<float,3>;
template class MATRIX<double,3>;
}
