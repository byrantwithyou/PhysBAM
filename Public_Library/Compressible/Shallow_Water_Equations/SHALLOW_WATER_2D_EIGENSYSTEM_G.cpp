//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_EIGENSYSTEM_G  
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_EIGENSYSTEM_G.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_2D_EIGENSYSTEM_G<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++){
        T one_over_h=1/U(i)(0),v=U(i)(2)/U(i)(0);
        G(i)(0)=U(i)(2);                                           // h*v
        G(i)(1)=U(i)(1)*U(i)(2)*one_over_h;               // h*u*v
        G(i)(2)=U(i)(2)*v+(T).5*gravity*sqr(U(i)(0));} // h*v^2+.5*g*h^2
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for G(U) at point cell
template<class T> T SHALLOW_WATER_2D_EIGENSYSTEM_G<T>::
Maximum_Magnitude_Eigenvalue(const VECTOR<T,3>& U_cell)
{
    T v=U_cell(2)/U_cell(0);
    T celerity=0;if(U_cell(0) >= 0) celerity=sqrt(gravity*U_cell(0));
    return maxabs(v-celerity,v+celerity);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool SHALLOW_WATER_2D_EIGENSYSTEM_G<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool weakly_hyperbolic=false;

    // eigenvalues on the left - at point i
    T v=U(i)(2)/U(i)(0);
    T celerity=0;if(U(i)(0) >= 0) celerity=sqrt(gravity*U(i)(0));else weakly_hyperbolic=true;
    lambda_left(0)=v-celerity;
    lambda_left(1)=v;
    lambda_left(2)=v+celerity;
        
    // eigenvalues on the right - at point i+1
    v=U(i+1)(2)/U(i+1)(0);
    celerity=0;if(U(i+1)(0) >= 0) celerity=sqrt(gravity*U(i+1)(0));else weakly_hyperbolic=true;
    lambda_right(0)=v-celerity;
    lambda_right(1)=v;
    lambda_right(2)=v+celerity;
        
    // eigenvalues in the center - at flux i
    T h=(T).5*(U(i)(0)+U(i+1)(0)),h_v=(T).5*(U(i)(2)+U(i+1)(2));
    v=h_v/h;
    celerity=0;if(h >= 0) celerity=sqrt(gravity*h);else weakly_hyperbolic=true;
    lambda(0)=v-celerity;
    lambda(1)=v;
    lambda(2)=v+celerity;

    if(weakly_hyperbolic) return false; // loss of hyperbolicity
    else return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void SHALLOW_WATER_2D_EIGENSYSTEM_G<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T h=(T).5*(U(i)(0)+U(i+1)(0)),h_u=(T).5*(U(i)(1)+U(i+1)(1)),h_v=(T).5*(U(i)(2)+U(i+1)(2));
    T one_over_h=1/h,u=h_u*one_over_h,v=h_v*one_over_h;
    T celerity=sqrt(gravity*h),one_over_2_celerity=1/(1*celerity);
                    
    L(0,0)=one_over_2_celerity*(v+celerity);L(0,1)=0;L(0,2)=-one_over_2_celerity;
    L(1,0)=-u;L(1,1)=1;L(1,2)=0;
    L(2,0)=-one_over_2_celerity*(v-celerity);L(2,1)=0;L(2,2)=one_over_2_celerity;
    
    R(0,0)=1;R(0,1)=u;R(0,2)=v-celerity;
    R(1,0)=0;R(1,1)=1;R(1,2)=0;
    R(2,0)=1;R(2,1)=u;R(2,2)=v+celerity;
}  
//#####################################################################
namespace PhysBAM{
template class SHALLOW_WATER_2D_EIGENSYSTEM_G<float>;
template class SHALLOW_WATER_2D_EIGENSYSTEM_G<double>;
}
