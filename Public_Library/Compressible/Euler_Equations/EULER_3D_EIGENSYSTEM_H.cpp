//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// H(U) for ij in (-2,mn+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_H<T>::
Flux(const int mn,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    if(only_pressure_flux){
        for(int ij=-3;ij<mn+3;ij++){
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=0;
            H(ij)(1)=0;
            H(ij)(2)=0;
            H(ij)(3)=p;
            H(ij)(4)=0;}
        return;}

    if(U_clamped){
        for(int ij=-3;ij<mn+3;ij++){
            T u=U(ij)(1)/U(ij)(0),v=U(ij)(2)/U(ij)(0),w=U(ij)(3)/U(ij)(0);
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=(*U_clamped)(ij)(3);         // rho_clamped*w
            H(ij)(1)=(*U_clamped)(ij)(3)*u;       // rho_clamped*u*w
            H(ij)(2)=(*U_clamped)(ij)(3)*v;       // rho_clamped*v*w
            H(ij)(3)=(*U_clamped)(ij)(3)*w+p;     // rho_clamped*w^2+p
            H(ij)(4)=((*U_clamped)(ij)(4)+p)*w;}} // (E_from_rho_clamped+p)*w
    else{
        for(int ij=-3;ij<mn+3;ij++){
            T u=U(ij)(1)/U(ij)(0),v=U(ij)(2)/U(ij)(0),w=U(ij)(3)/U(ij)(0);
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=U(ij)(3);         // rho*w
            H(ij)(1)=U(ij)(3)*u;       // rho*u*w
            H(ij)(2)=U(ij)(3)*v;       // rho*v*w
            H(ij)(3)=U(ij)(3)*w+p;     // rho*w^2+p
            H(ij)(4)=(U(ij)(4)+p)*w;}} // (E+p)*w
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_3D_EIGENSYSTEM_H<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(only_pressure_flux){
        for(int ij=range.x;ij<range.y;ij++){
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=0;
            H(ij)(1)=0;
            H(ij)(2)=0;
            H(ij)(3)=p;
            H(ij)(4)=0;}
        return;}

    T average_w;
    if(use_standard_average) average_w=(U(face_index)(3)/U(face_index)(0)+U(face_index+1)(3)/U(face_index+1)(0))*(T).5;
    else average_w=(U(face_index)(3)+U(face_index+1)(3))/(U(face_index)(0)+U(face_index+1)(0));

    if(U_clamped){
        for(int ij=range.x;ij<range.y;ij++){
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=(*U_clamped)(ij)(0)*average_w;
            H(ij)(1)=(*U_clamped)(ij)(1)*average_w;
            H(ij)(2)=(*U_clamped)(ij)(2)*average_w;
            H(ij)(3)=(*U_clamped)(ij)(3)*average_w+p;
            H(ij)(4)=((*U_clamped)(ij)(4)+p)*average_w;}}
    else{
        for(int ij=range.x;ij<range.y;ij++){
            T p=eos->p(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
            H(ij)(0)=U(ij)(0)*average_w;
            H(ij)(1)=U(ij)(1)*average_w;
            H(ij)(2)=U(ij)(2)*average_w;
            H(ij)(3)=U(ij)(3)*average_w+p;
            H(ij)(4)=(U(ij)(4)+p)*average_w;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for H(U) at point cell
template<class T> T EULER_3D_EIGENSYSTEM_H<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T w=U_cell(3)/U_cell(0);
    T sound_speed=eos->c(U_cell(0),e(U_cell(0),U_cell(1),U_cell(2),U_cell(3),U_cell(4)));
    return maxabs(w-sound_speed,w+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for H(U) at flux ij and and at points ij and ij+1
template<class T> bool EULER_3D_EIGENSYSTEM_H<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int ij,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point ij
    T w=U(ij)(3)/U(ij)(0);
    T sound_speed=eos->c(U(ij)(0),e(U(ij)(0),U(ij)(1),U(ij)(2),U(ij)(3),U(ij)(4)));
    if(sound_speed == 0) cavitation=1;
    lambda_left(0)=w-sound_speed;
    lambda_left(1)=lambda_left(2)=lambda_left(3)=w;
    lambda_left(4)=w+sound_speed;
        
    // eigenvalues on the right - at point ij+1
    w=U(ij+1)(3)/U(ij+1)(0);
    sound_speed=eos->c(U(ij+1)(0),e(U(ij+1)(0),U(ij+1)(1),U(ij+1)(2),U(ij+1)(3),U(ij+1)(4)));
    if(sound_speed == 0) cavitation=1;
    lambda_right(0)=w-sound_speed;
    lambda_right(1)=lambda_right(2)=lambda_right(3)=w;
    lambda_right(4)=w+sound_speed;
        
    // eigenvalues in the center - at flux ij
    T rho=(U(ij)(0)+U(ij+1)(0))/2;
    T rho_u=(U(ij)(1)+U(ij+1)(1))/2;
    T rho_v=(U(ij)(2)+U(ij+1)(2))/2;
    T rho_w=(U(ij)(3)+U(ij+1)(3))/2;
    T E=(U(ij)(4)+U(ij+1)(4))/2;
    w=rho_w/rho;
    T internal_energy=e(rho,rho_u,rho_v,rho_w,E);
    sound_speed=eos->c(rho,internal_energy);
    if(sound_speed == 0) cavitation=1;
    lambda(0)=w-sound_speed;
    lambda(1)=lambda(2)=lambda(3)=w;
    lambda(4)=w+sound_speed;

    if(cavitation) return 0; // loss of hyperbolicity
    else return 1; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for H(U) at flux ij between points ij and ij+1
template<class T> void EULER_3D_EIGENSYSTEM_H<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int ij,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux ij
    T rho=(U(ij)(0)+U(ij+1)(0))/2;
    T rho_u=(U(ij)(1)+U(ij+1)(1))/2;
    T rho_v=(U(ij)(2)+U(ij+1)(2))/2;
    T rho_w=(U(ij)(3)+U(ij+1)(3))/2;
    T E=(U(ij)(4)+U(ij+1)(4))/2;
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
    T w_over_2c=w*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b1_over_2_times_u=b1_over_2*u;
    T b1_over_2_times_v=b1_over_2*v;
    T b1_over_2_times_w=b1_over_2*w;
    T w_times_c=w*sound_speed;
                    
    L(0,0)=b2_over_2+w_over_2c;
    L(0,1)=-b1_over_2_times_u;
    L(0,2)=-b1_over_2_times_v;
    L(0,3)=-b1_over_2_times_w-one_over_2c;
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
    L(3,0)=v;
    L(3,1)=0;
    L(3,2)=-1;
    L(3,3)=0;
    L(3,4)=0;
    L(4,0)=b2_over_2-w_over_2c;
    L(4,1)=-b1_over_2_times_u;
    L(4,2)=-b1_over_2_times_v;
    L(4,3)=-b1_over_2_times_w+one_over_2c;
    L(4,4)=b1_over_2;
    
    R(0,0)=1;
    R(0,1)=u;
    R(0,2)=v;
    R(0,3)=w-sound_speed;
    R(0,4)=h-w_times_c;
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
    R(3,2)=-1;
    R(3,3)=0;
    R(4,4)=-v;
    R(4,0)=1;
    R(4,1)=u;
    R(4,2)=v;
    R(4,3)=w+sound_speed;
    R(4,4)=h+w_times_c;
}  
//#####################################################################
namespace PhysBAM{
template class EULER_3D_EIGENSYSTEM_H<float>;
template class EULER_3D_EIGENSYSTEM_H<double>;
}
