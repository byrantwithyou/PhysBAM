//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_3X3
//#####################################################################
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/exchange_sort.h>
#include <Tools/Math_Tools/maxabs.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X3.h>
#include <Tools/Matrices/MATRIX_3X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Function Fast_Eigenvalues
//#####################################################################
// lambda_x > lambda_y > lambda_z
// reference: Smith, O. "Eigenvalues of a symmetric 3 x 3 matrix". Commun. ACM 4 (4), p. 168, 1961 (thanks, Gene)
template<class T> DIAGONAL_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Fast_Eigenvalues() const // 24 mults, 20 adds, 1 atan2, 1 sincos, 2 sqrts
{
    if(!IS_SAME<T,double>::value) return DIAGONAL_MATRIX<T,3>(SYMMETRIC_MATRIX<double,3>(*this).Fast_Eigenvalues());
    // now T is double
    T m=(T)one_third*(x00+x11+x22);
    T a00=x00-m,a11=x11-m,a22=x22-m,a12_sqr=x10*x10,a13_sqr=x20*x20,a23_sqr=x21*x21;
    T p=(T)one_sixth*(a00*a00+a11*a11+a22*a22+2*(a12_sqr+a13_sqr+a23_sqr));
    T q=(T).5*(a00*(a11*a22-a23_sqr)-a11*a13_sqr-a22*a12_sqr)+x10*x20*x21;
    T sqrt_p=sqrt(p),disc=p*p*p-q*q;
    T phi=(T)one_third*atan2(sqrt(max((T)0,disc)),q),c=cos(phi),s=sin(phi);
    T sqrt_p_cos=sqrt_p*c,root_three_sqrt_p_sin=(T)root_three*sqrt_p*s;
    DIAGONAL_MATRIX<T,3> lambda(m+2*sqrt_p_cos,m-sqrt_p_cos-root_three_sqrt_p_sin,m-sqrt_p_cos+root_three_sqrt_p_sin);
    lambda.x.Sort();
    return lambda;
}
//#####################################################################
// Function Fast_Eigenvectors
//#####################################################################
namespace{
template<class T> MATRIX<T,3>
Fast_Eigenvectors(const SYMMETRIC_MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& lambda) // 71 mults, 44 adds, 3 divs, 3 sqrts
{
    if(!IS_SAME<T,double>::value) PHYSBAM_FATAL_ERROR();
    // T is now always double

    // flip if necessary so that first eigenvalue is the most different
    bool flipped=false;
    DIAGONAL_MATRIX<T,3> lambda_flip(lambda);
    if(lambda.x.x-lambda.x.y<lambda.x.y-lambda.x.z){ // 2a
        exchange(lambda_flip.x.x,lambda_flip.x.z);
        flipped=true;}

    // get first eigenvector
    VECTOR<T,3> v1=(A-lambda_flip.x.x).Cofactor_Matrix().Largest_Column_Normalized(); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s

    // form basis for orthogonal complement to v1, and reduce A to this space
    VECTOR<T,3> v1_orthogonal=v1.Unit_Orthogonal_Vector(); // 6m+2a+1d+1s (tweak: 5m+1a+1d+1s)
    MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(v1,v1_orthogonal)); // 6m+3a (tweak: 4m+1a)
    SYMMETRIC_MATRIX<T,2> A_reduced=SYMMETRIC_MATRIX<T,2>::Conjugate_With_Transpose(other_v,A); // 21m+12a (tweak: 18m+9a)

    // find third eigenvector from A_reduced, and fill in second via cross product
    VECTOR<T,3> v3=other_v*(A_reduced-lambda_flip.x.z).Cofactor_Matrix().Largest_Column_Normalized(); // 6m+3a + 2a + 5m+2a+1d+1s = 11m+7a+1d+1s (tweak: 10m+6a+1d+1s)
    VECTOR<T,3> v2=VECTOR<T,3>::Cross_Product(v3,v1); // 6m+3a

    // finish
    return flipped?MATRIX<T,3>(v3,v2,-v1):MATRIX<T,3>(v1,v2,v3);
}
}
//#####################################################################
// Function Fast_Solve_Eigenproblem
//#####################################################################
template<class T> void SYMMETRIC_MATRIX<T,3>::
Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const // roughly 95 mults, 64 adds, 3 divs, 5 sqrts, 1 atan2, 1 sincos
{
    if(!IS_SAME<T,double>::value){
        DIAGONAL_MATRIX<double,3> eigenvalues_double;MATRIX<double,3> eigenvectors_double;
        SYMMETRIC_MATRIX<double,3>(*this).Fast_Solve_Eigenproblem(eigenvalues_double,eigenvectors_double);
        eigenvalues=DIAGONAL_MATRIX<T,3>(eigenvalues_double);eigenvectors=MATRIX<T,3>(eigenvectors_double);return;}
    // now T is double
    eigenvalues=Fast_Eigenvalues();
    eigenvectors=Fast_Eigenvectors(*this,eigenvalues);
}
//#####################################################################
// Function Solve_Eigenproblem
//####################################################################
template<class T> void SYMMETRIC_MATRIX<T,3>::
Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const
{
    T a00=x00,a01=x10,a02=x20,a11=x11,a12=x21,a22=x22;
    T v00=1,v01=0,v02=0,v10=0,v11=1,v12=0,v20=0,v21=0,v22=1;
    int sweep;for(sweep=0;sweep<50;sweep++){
        T sum=abs(a01)+abs(a02)+abs(a12);
        if(sum==0) break;
        T threshold=sweep<4?(T)(1./45)*sum:0;
        Jacobi_Transform(sweep,threshold,a00,a01,a11,a02,a12,v00,v01,v10,v11,v20,v21);
        Jacobi_Transform(sweep,threshold,a00,a02,a22,a01,a12,v00,v02,v10,v12,v20,v22);
        Jacobi_Transform(sweep,threshold,a11,a12,a22,a01,a02,v01,v02,v11,v12,v21,v22);}
    assert(sweep<=50);
    eigenvalues=DIAGONAL_MATRIX<T,3>(a00,a11,a22);
    eigenvectors=MATRIX<T,3>(v00,v10,v20,v01,v11,v21,v02,v12,v22);
}
//#####################################################################
// Function Jacobi_Transform
//#####################################################################
template<class T> inline void SYMMETRIC_MATRIX<T,3>::
Jacobi_Transform(const int sweep,const T threshold,T& app,T& apq,T& aqq,T& arp,T& arq,T& v1p,T& v1q,T& v2p,T& v2q,T& v3p,T& v3q)
{
    T epsilon=std::numeric_limits<T>::epsilon();
    T g=100*abs(apq);
    if(sweep > 4 && epsilon*abs(app)>=g && epsilon*abs(aqq)>=g) apq=0;
    else if(abs(apq) > threshold){
        T h=aqq-app,t;
        if(epsilon*abs(h)>=g) t=apq/h;
        else{T theta=(T).5*h/apq;t=1/(abs(theta)+sqrt(1+sqr(theta)));if(theta < 0) t=-t;}
        T cosine=1/sqrt(1+t*t),sine=t*cosine,tau=sine/(1+cosine),da=t*apq;
        app-=da;aqq+=da;apq=0;
        T arp_tmp=arp-sine*(arq+tau*arp);arq=arq+sine*(arp-tau*arq);arp=arp_tmp;
        T v1p_tmp=v1p-sine*(v1q+tau*v1p);v1q=v1q+sine*(v1p-tau*v1q);v1p=v1p_tmp;
        T v2p_tmp=v2p-sine*(v2q+tau*v2p);v2q=v2q+sine*(v2p-tau*v2q);v2p=v2p_tmp;
        T v3p_tmp=v3p-sine*(v3q+tau*v3p);v3q=v3q+sine*(v3p-tau*v3q);v3p=v3p_tmp;}
}
//#####################################################################
namespace PhysBAM{
template class SYMMETRIC_MATRIX<float,3>;
template class SYMMETRIC_MATRIX<double,3>;
}
