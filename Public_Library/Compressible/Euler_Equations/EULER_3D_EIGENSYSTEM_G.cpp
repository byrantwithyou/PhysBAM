//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_G.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// G(U) for j in (-2,n+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_G<T>::
Flux(const int n,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    if(only_pressure_flux){
        for(int j=-3;j<n+3;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=0;
            G(j)(1)=0;
            G(j)(2)=p;
            G(j)(3)=0;
            G(j)(4)=0;}
        return;}

    if(U_clamped){
        for(int j=-3;j<n+3;j++){
            T u=U(j)(1)/U(j)(0),v=U(j)(2)/U(j)(0),w=U(j)(3)/U(j)(0);
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=(*U_clamped)(j)(2);         // rho_clamped*v
            G(j)(1)=(*U_clamped)(j)(2)*u;       // rho_clamped*u*v
            G(j)(2)=(*U_clamped)(j)(2)*v+p;     // rho_clamped*v^2_p
            G(j)(3)=(*U_clamped)(j)(2)*w;       // rho_clamped*v*w
            G(j)(4)=((*U_clamped)(j)(4)+p)*v;}} // (E_from_rho_clamped+p)*v
    else{
        for(int j=-3;j<n+3;j++){
            T u=U(j)(1)/U(j)(0),v=U(j)(2)/U(j)(0),w=U(j)(3)/U(j)(0);
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=U(j)(2);         // rho*v
            G(j)(1)=U(j)(2)*u;       // rho*u*v
            G(j)(2)=U(j)(2)*v+p;     // rho*v^2_p
            G(j)(3)=U(j)(2)*w;       // rho*v*w
            G(j)(4)=(U(j)(4)+p)*v;}} // (E+p)*v
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_3D_EIGENSYSTEM_G<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(only_pressure_flux){
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=0;
            G(j)(1)=0;
            G(j)(2)=p;
            G(j)(3)=0;
            G(j)(4)=0;}
        return;}

    T average_v;
    if(use_standard_average) average_v=(U(face_index)(2)/U(face_index)(0)+U(face_index+1)(2)/U(face_index+1)(0))*(T).5;
    else average_v=(U(face_index)(2)+U(face_index+1)(2))/(U(face_index)(0)+U(face_index+1)(0));

    if(U_clamped){
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=(*U_clamped)(j)(0)*average_v;
            G(j)(1)=(*U_clamped)(j)(1)*average_v;
            G(j)(2)=(*U_clamped)(j)(2)*average_v+p;
            G(j)(3)=(*U_clamped)(j)(3)*average_v;
            G(j)(4)=((*U_clamped)(j)(4)+p)*average_v;}}
    else{
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
            G(j)(0)=U(j)(0)*average_v;
            G(j)(1)=U(j)(1)*average_v;
            G(j)(2)=U(j)(2)*average_v+p;
            G(j)(3)=U(j)(3)*average_v;
            G(j)(4)=(U(j)(4)+p)*average_v;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for G(U) at point cell
template<class T> T EULER_3D_EIGENSYSTEM_G<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T v=U_cell(2)/U_cell(0);
    T sound_speed=eos->c(U_cell(0),e(U_cell(0),U_cell(1),U_cell(2),U_cell(3),U_cell(4)));
    return maxabs(v-sound_speed,v+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for G(U) at flux j and and at points j and j+1
template<class T> bool EULER_3D_EIGENSYSTEM_G<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point j
    T v=U(j)(2)/U(j)(0);
    T sound_speed=eos->c(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3),U(j)(4)));
    if(sound_speed == 0) cavitation=1;
    lambda_left(0)=v-sound_speed;
    lambda_left(1)=lambda_left(2)=lambda_left(3)=v;
    lambda_left(4)=v+sound_speed;
        
    // eigenvalues on the right - at point j+1
    v=U(j+1)(2)/U(j+1)(0);
    sound_speed=eos->c(U(j+1)(0),e(U(j+1)(0),U(j+1)(1),U(j+1)(2),U(j+1)(3),U(j+1)(4)));
    if(sound_speed == 0) cavitation=1;
    lambda_right(0)=v-sound_speed;
    lambda_right(1)=lambda_right(2)=lambda_right(3)=v;
    lambda_right(4)=v+sound_speed;
        
    // eigenvalues in the center - at flux j
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T rho_w=(U(j)(3)+U(j+1)(3))/2;
    T E=(U(j)(4)+U(j+1)(4))/2;
    v=rho_v/rho;
    T internal_energy=e(rho,rho_u,rho_v,rho_w,E);
    sound_speed=eos->c(rho,internal_energy);
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
template<class T> void EULER_3D_EIGENSYSTEM_G<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux j
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T rho_w=(U(j)(3)+U(j+1)(3))/2;
    T E=(U(j)(4)+U(j+1)(4))/2;
    T internal_energy=e(rho,rho_u,rho_v,rho_w,E);
    T sound_speed=eos->c(rho,internal_energy);
    T p=eos->p(rho,internal_energy);
        
    T u=rho_u/rho;
    T v=rho_v/rho;
    T w=rho_w/rho;
    T q2=sqr(u)+sqr(v)+sqr(w);
    T h=(E+p)/rho;
    T b1=eos->p_e(rho,internal_energy)/(rho*sqr(sound_speed));
    T b2=1+b1*(q2-h);
        
    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T v_over_2c=v*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b1_over_2_times_u=b1_over_2*u;
    T b1_over_2_times_v=b1_over_2*v;
    T b1_over_2_times_w=b1_over_2*w;
    T v_times_c=v*sound_speed;
                    
    L(0,0)=b2_over_2+v_over_2c;
    L(0,1)=-b1_over_2_times_u;
    L(0,2)=-b1_over_2_times_v-one_over_2c;
    L(0,3)=-b1_over_2_times_w;
    L(0,4)=b1_over_2;
    L(1,0)=h-q2;
    L(1,1)=u;
    L(1,2)=v;
    L(1,3)=w;
    L(1,4)=-1;
    L(2,0)=u;
    L(2,1)=-1;
    L(2,2)=0;
    L(2,3)=0;
    L(2,4)=0;
    L(3,0)=w;
    L(3,1)=0;
    L(3,2)=0;
    L(3,3)=-1;
    L(3,4)=0;
    L(4,0)=b2_over_2-v_over_2c;
    L(4,1)=-b1_over_2_times_u;
    L(4,2)=-b1_over_2_times_v+one_over_2c;
    L(4,3)=-b1_over_2_times_w;
    L(4,4)=b1_over_2;
    
    R(0,0)=1;
    R(0,1)=u;
    R(0,2)=v-sound_speed;
    R(0,3)=w;
    R(0,4)=h-v_times_c;
    R(1,0)=b1;
    R(1,1)=b1*u;
    R(1,2)=b1*v;
    R(1,3)=b1*w;
    R(1,4)=b1*h-1;
    R(2,0)=0;
    R(2,1)=-1;
    R(2,2)=0;
    R(2,3)=0;
    R(2,4)=-u;
    R(3,0)=0;
    R(3,1)=0;
    R(3,2)=0;
    R(3,3)=-1;
    R(3,4)=-w;
    R(4,0)=1;
    R(4,1)=u;
    R(4,2)=v+sound_speed;
    R(4,3)=w;
    R(4,4)=h+v_times_c;
}  
//#####################################################################
namespace PhysBAM{
template class EULER_3D_EIGENSYSTEM_G<float>;
template class EULER_3D_EIGENSYSTEM_G<double>;
}
