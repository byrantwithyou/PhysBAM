//#####################################################################
// Copyright 2003-2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++){
        F(i)(0)=U(i)(0)*U(i)(1); // h*u
        F(i)(1)=(T).5*sqr(U(i)(1))+gravity*eta_ghost(i);} // .5*u^2+g*eta
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int weakly_hyperbolic=0;

    // eigenvalues on the left - at point i
    T u=U(i)(1);
    T celerity=0;if(U(i)(0) >= 0) celerity=sqrt(gravity*U(i)(0));else weakly_hyperbolic=1;
    lambda_left(0)=u-celerity;
    lambda_left(1)=u+celerity;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(1);
    celerity=0;if(U(i+1)(0) >= 0) celerity=sqrt(gravity*U(i+1)(0));else weakly_hyperbolic=1;
    lambda_right(0)=u-celerity;
    lambda_right(1)=u+celerity;
        
    // eigenvalues in the center - at flux i
    T h=(U(i)(0)+U(i+1)(0))/2;
    u=(U(i)(1)+U(i+1)(1))/2;
    celerity=0;if(h >= 0) celerity=sqrt(gravity*h);else weakly_hyperbolic=1;
    lambda(0)=u-celerity;
    lambda(1)=u+celerity;

    if(weakly_hyperbolic) return false; // loss of hyperbolicity
    else return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T h=(U(i)(0)+U(i+1)(0))/2;
    T sqrt_h_over_gravity=sqrt(h/gravity);
    T sqrt_gravity_over_h=0;if(h>min_height) sqrt_gravity_over_h=1/sqrt_h_over_gravity;
                    
    T half_sqrt_h_over_gravity=(T).5*sqrt_h_over_gravity;
    L(0,0)=(T)-.5;L(0,1)=half_sqrt_h_over_gravity;
    L(1,0)=(T).5;L(1,1)=half_sqrt_h_over_gravity;
    
    R(0,0)=-1;R(0,1)=sqrt_gravity_over_h;
    R(1,0)=1;R(1,1)=sqrt_gravity_over_h;
}  
//#####################################################################
template class SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F<double>;
#endif
