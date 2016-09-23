//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/sqr.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void REACTIVE_EULER_2D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++){
        T u=U(i)(1)/U(i)(0),v=U(i)(2)/U(i)(0);
        T Y=U(i)(4)/U(i)(0);
        T p=eos.p(U(i)(0),e(U(i)(0),U(i)(1),U(i)(2),U(i)(3)),Y);
        F(i)(0)=U(i)(1);        // rho*u
        F(i)(1)=U(i)(1)*u+p;    // rho*u^2+p
        F(i)(2)=U(i)(1)*v;    // rho*u*v
        F(i)(3)=(U(i)(3)+p)*u; // (E+p)*u
        F(i)(4)=U(i)(4)*u;} // rho*Y*u
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T REACTIVE_EULER_2D_EIGENSYSTEM_F<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(1)/U_cell(0);
    T Y=U_cell(4)/U_cell(0);
    T sound_speed=eos.c(U_cell(0),e(U_cell(0),U_cell(1),U_cell(2),U_cell(3)),Y);
    return maxabs(u-sound_speed,u+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool REACTIVE_EULER_2D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point i
    T u=U(i)(1)/U(i)(0);
    T Y=U(i)(4)/U(i)(0);
    T sound_speed=eos.c(U(i)(0),e(U(i)(0),U(i)(1),U(i)(2),U(i)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_left(0)=u-sound_speed;
    lambda_left(1)=lambda_left(2)=lambda_left(3)=u;
    lambda_left(4)=u+sound_speed;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(1)/U(i+1)(0);
    Y=U(i+1)(4)/U(i+1)(0);
    sound_speed=eos.c(U(i+1)(0),e(U(i+1)(0),U(i+1)(1),U(i+1)(2),U(i+1)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_right(0)=u-sound_speed;
    lambda_right(1)=lambda_right(2)=lambda_right(3)=u;
    lambda_right(4)=u+sound_speed;
        
    // eigenvalues in the center - at flux i
    T rho=(U(i)(0)+U(i+1)(0))/2;
    T rho_u=(U(i)(1)+U(i+1)(1))/2;
    T rho_v=(U(i)(2)+U(i+1)(2))/2;
    T E=(U(i)(3)+U(i+1)(3))/2;
    T rho_Y=(U(i)(4)+U(i+1)(4))/2;
    u=rho_u/rho;
    Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    sound_speed=eos.c(rho,internal_energy,Y);
    if(sound_speed == 0) cavitation=1;
    lambda(0)=u-sound_speed;
    lambda(1)=lambda(2)=lambda(3)=u;
    lambda(4)=u+sound_speed;

    if(cavitation) return 0; // loss of hyperbolicity
    else return 1; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void REACTIVE_EULER_2D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T rho=(U(i)(0)+U(i+1)(0))/2;
    T rho_u=(U(i)(1)+U(i+1)(1))/2;
    T rho_v=(U(i)(2)+U(i+1)(2))/2;
    T E=(U(i)(3)+U(i+1)(3))/2;
    T rho_Y=(U(i)(4)+U(i+1)(4))/2;
    T u=rho_u/rho;
    T Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    T sound_speed=eos.c(rho,internal_energy,Y);
        
    T v=rho_v/rho;
    T q2=sqr(u)+sqr(v);
    T h=(E+eos.p(rho,internal_energy,Y))/rho; // capital H=(E+p)/rho
    T b1=eos.p_e(rho,internal_energy,Y)/(rho*sqr(sound_speed));
    T b2=1+b1*(q2-h);
    T z=eos.e_not1-eos.e_not2;
    T b3=b1*Y*z;
        
    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T u_over_2c=u*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b3_over_2=b3/2;
    T b1_over_2_times_u=b1_over_2*u;
    T b1_over_2_times_v=b1_over_2*v;
    T u_times_c=u*sound_speed;
                    
    L(0,1)=b2_over_2+u_over_2c+b3_over_2;
    L(0,2)=-b1_over_2_times_u-one_over_2c;
    L(0,3)=-b1_over_2_times_v;
    L(0,4)=b1_over_2;
    L(0,5)=-b1_over_2*z;
    L(1,1)=1-b2-b3;
    L(1,2)=b1*u;
    L(1,3)=b1*v;
    L(1,4)=-b1;
    L(1,5)=b1*z;
    L(2,1)=v;
    L(2,2)=0;
    L(2,3)=-1;
    L(2,4)=0;
    L(2,5)=0;
    L(3,1)=-Y;
    L(3,2)=0;
    L(3,3)=0;
    L(3,4)=0;
    L(3,5)=1;
    L(4,1)=b2_over_2-u_over_2c+b3_over_2;
    L(4,2)=-b1_over_2_times_u+one_over_2c;
    L(4,3)=-b1_over_2_times_v;
    L(4,4)=b1_over_2;
    L(4,5)=-b1_over_2*z;

    R(0,1)=1;
    R(0,2)=u-sound_speed;
    R(0,3)=v;
    R(0,4)=h-u_times_c;
    R(0,5)=Y;
    R(1,1)=1;
    R(1,2)=u;
    R(1,3)=v;
    R(1,4)=h-1/b1;
    R(1,5)=Y;
    R(2,1)=0;
    R(2,2)=0;
    R(2,3)=-1;
    R(2,4)=-v;
    R(2,5)=0;
    R(3,1)=0;
    R(3,2)=0;
    R(3,3)=0;
    R(3,4)=z;
    R(3,5)=1;
    R(4,1)=1;
    R(4,2)=u+sound_speed;
    R(4,3)=v;
    R(4,4)=h+u_times_c;
    R(4,5)=Y;
}  
//#####################################################################
#if 0 // broken
template class REACTIVE_EULER_2D_EIGENSYSTEM_F<float>;
template class REACTIVE_EULER_2D_EIGENSYSTEM_F<double>;
#endif
