//#####################################################################
// Copyright 2002-2007, Ron Fedkiw, Nipun Kwatra, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D_EIGENSYSTEM_G.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void REACTIVE_EULER_2D_EIGENSYSTEM_G<T>::
Flux(const int n,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int j=-3;j<n+3;j++){
        T u=U(j)(1)/U(j)(0),v=U(j)(2)/U(j)(0);
        T Y=U(j)(4)/U(j)(0);
        T p=eos.p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)),Y);
        G(j)(0)=U(j)(2);        // rho*v
        G(j)(1)=U(j)(2)*u;    // rho*v*u
        G(j)(2)=U(j)(2)*v+p;    // rho*v^2+p
        G(j)(3)=(U(j)(3)+p)*v; // (E+p)*u
        G(j)(4)=U(j)(4)*v;} // rho*Y*v
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for G(U) at point cell
template<class T> T REACTIVE_EULER_2D_EIGENSYSTEM_G<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T v=U_cell(2)/U_cell(0);
    T Y=U_cell(4)/U_cell(0);
    T sound_speed=eos.c(U_cell(0),e(U_cell(0),U_cell(1),U_cell(2),U_cell(3)),Y);
    return maxabs(v-sound_speed,v+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for G(U) at flux i and and at points i and i+1
template<class T> bool REACTIVE_EULER_2D_EIGENSYSTEM_G<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point j
    T v=U(j)(2)/U(j)(0);
    T Y=U(j)(4)/U(j)(0);
    T sound_speed=eos.c(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_left(0)=v-sound_speed;
    lambda_left(1)=lambda_left(2)=lambda_left(3)=v;
    lambda_left(4)=v+sound_speed;
        
    // eigenvalues on the right - at point j+1
    v=U(j+1)(2)/U(j+1)(0);
    Y=U(j+1)(4)/U(j+1)(0);
    sound_speed=eos.c(U(j+1)(0),e(U(j+1)(0),U(j+1)(1),U(j+1)(2),U(j+1)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_right(0)=v-sound_speed;
    lambda_right(1)=lambda_right(2)=lambda_right(3)=v;
    lambda_right(4)=v+sound_speed;
        
    // eigenvalues in the center - at flux i
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T E=(U(j)(3)+U(j+1)(3))/2;
    T rho_Y=(U(j)(4)+U(j+1)(4))/2;
    v=rho_v/rho;
    Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    sound_speed=eos.c(rho,internal_energy,Y);
    if(sound_speed == 0) cavitation=1;
    lambda(0)=v-sound_speed;
    lambda(1)=lambda(2)=lambda(3)=v;
    lambda(4)=v+sound_speed;

    if(cavitation) return 0; // loss of hyperbolicity
    else return 1; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for G(U) at flux j between points j and j+1
template<class T> void REACTIVE_EULER_2D_EIGENSYSTEM_G<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux j
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T E=(U(j)(3)+U(j+1)(3))/2;
    T rho_Y=(U(j)(4)+U(j+1)(4))/2;
    T v=rho_v/rho;
    T Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    T sound_speed=eos.c(rho,internal_energy,Y);
        
    T u=rho_u/rho;
    T q2=sqr(u)+sqr(v);
    T h=(E+eos.p(rho,internal_energy,Y))/rho; // capital H=(E+p)/rho
    T b1=eos.p_e(rho,internal_energy,Y)/(rho*sqr(sound_speed));
    T b2=1+b1*(q2-h);
    T z=eos.e_not1-eos.e_not2;
    T b3=b1*Y*z;
        
    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T v_over_2c=v*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b3_over_2=b3/2;
    T b1_over_2_times_u=b1_over_2*u;
    T b1_over_2_times_v=b1_over_2*v;
    T v_times_c=v*sound_speed;
                    
    L(0,0)=b2_over_2+v_over_2c+b3_over_2;
    L(0,1)=-b1_over_2_times_u;
    L(0,2)=-b1_over_2_times_v-one_over_2c;
    L(0,3)=b1_over_2;
    L(0,4)=-b1_over_2*z;
    L(1,0)=1-b2-b3;
    L(1,1)=b1*u;
    L(1,2)=b1*v;
    L(1,3)=-b1;
    L(1,4)=b1*z;
    L(2,0)=-u;
    L(2,1)=1;
    L(2,2)=0;
    L(2,3)=0;
    L(2,4)=0;
    L(3,0)=-Y;
    L(3,1)=0;
    L(3,2)=0;
    L(3,3)=0;
    L(3,4)=1;
    L(4,0)=b2_over_2-v_over_2c+b3_over_2;
    L(4,1)=-b1_over_2_times_u;
    L(4,2)=-b1_over_2_times_v+one_over_2c;
    L(4,3)=b1_over_2;
    L(4,4)=-b1_over_2*z;
    
    R(0,0)=1;
    R(0,1)=u;
    R(0,2)=v-sound_speed;
    R(0,3)=h-v_times_c;
    R(0,4)=Y;
    R(1,0)=1;
    R(1,1)=u;
    R(1,2)=v;
    R(1,3)=h-1/b1;
    R(1,4)=Y;
    R(2,0)=0;
    R(2,1)=1;
    R(2,2)=0;
    R(2,3)=u;
    R(2,4)=0;
    R(3,0)=0;
    R(3,1)=0;
    R(3,2)=0;
    R(3,3)=z;
    R(3,4)=1;
    R(4,0)=1;
    R(4,1)=u;
    R(4,2)=v+sound_speed;
    R(4,3)=h+v_times_c;
    R(4,4)=Y;
}  
//#####################################################################
#if 0 // broken
template class REACTIVE_EULER_2D_EIGENSYSTEM_G<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REACTIVE_EULER_2D_EIGENSYSTEM_G<double>;
#endif
#endif
