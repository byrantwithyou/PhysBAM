//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D_EIGENSYSTEM_G  
//##################################################################### 
#ifndef __EULER_2D_EIGENSYSTEM_G__
#define __EULER_2D_EIGENSYSTEM_G__   

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D_EIGENSYSTEM_G:public EULER_EIGENSYSTEM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    typedef EULER_EIGENSYSTEM<GRID<VECTOR<T_input,2> > > BASE;
public:
    using BASE::eos;using BASE::e;

    bool only_pressure_flux;

    EULER_2D_EIGENSYSTEM_G(bool only_pressure_flux_input=false):only_pressure_flux(only_pressure_flux_input)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE; 
//#####################################################################
};   
//#####################################################################
// Function Flux
//#####################################################################
// G(U) for j in (-2,n+3) - 3 ghost cells
template<class T> void EULER_2D_EIGENSYSTEM_G<T>::
Flux(const int n,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    if(only_pressure_flux){
        for(int j=-3;j<n+3;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=0;
            G(j)(1)=0;
            G(j)(2)=p;
            G(j)(3)=0;}
        return;}

    if(U_clamped){
        for(int j=-3;j<n+3;j++){
            T u=U(j)(1)/U(j)(0),v=U(j)(2)/U(j)(0);
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=(*U_clamped)(j)(2);         // rho_clamped*v
            G(j)(1)=(*U_clamped)(j)(2)*u;       // rho_clamped*u*v
            G(j)(2)=(*U_clamped)(j)(2)*v+p;     // rho_clamped*v^2_p
            G(j)(3)=((*U_clamped)(j)(3)+p)*v;}} // (E_from_rho_clamped+p)*v
    else{
        for(int j=-3;j<n+3;j++){
            T u=U(j)(1)/U(j)(0),v=U(j)(2)/U(j)(0);
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=U(j)(2);         // rho*v
            G(j)(1)=U(j)(2)*u;       // rho*u*v
            G(j)(2)=U(j)(2)*v+p;     // rho*v^2_p
            G(j)(3)=(U(j)(3)+p)*v;}} // (E+p)*v
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_2D_EIGENSYSTEM_G<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(only_pressure_flux){
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=0;
            G(j)(1)=0;
            G(j)(2)=p;
            G(j)(3)=0;}
        return;}

    T average_v;
    if(use_standard_average) average_v=(U(face_index)(2)/U(face_index)(0)+U(face_index+1)(2)/U(face_index+1)(0))*(T).5;
    else average_v=(U(face_index)(2)+U(face_index+1)(2))/(U(face_index)(0)+U(face_index+1)(0));

    if(U_clamped){
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=(*U_clamped)(j)(0)*average_v;
            G(j)(1)=(*U_clamped)(j)(1)*average_v;
            G(j)(2)=(*U_clamped)(j)(2)*average_v+p;
            G(j)(3)=((*U_clamped)(j)(3)+p)*average_v;}}
    else{
        for(int j=range.x;j<range.y;j++){
            T p=eos->p(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
            G(j)(0)=U(j)(0)*average_v;
            G(j)(1)=U(j)(1)*average_v;
            G(j)(2)=U(j)(2)*average_v+p;
            G(j)(3)=(U(j)(3)+p)*average_v;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for G(U) at point cell
template<class T> T EULER_2D_EIGENSYSTEM_G<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T v=U_cell(2)/U_cell(0);
    T sound_speed=eos->c(U_cell(0),e(U_cell(0),U_cell(1),U_cell(2),U_cell(3)));
    return maxabs(v-sound_speed,v+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for G(U) at flux j and and at points j and j+1
template<class T> bool EULER_2D_EIGENSYSTEM_G<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool cavitation=false;

    // eigenvalues on the left - at point j
    T v=U(j)(2)/U(j)(0);
    T sound_speed=eos->c(U(j)(0),e(U(j)(0),U(j)(1),U(j)(2),U(j)(3)));
    if(sound_speed == 0) cavitation=true;
    lambda_left(0)=v-sound_speed;
    lambda_left(1)=lambda_left(2)=v;
    lambda_left(3)=v+sound_speed;
        
    // eigenvalues on the right - at point j+1
    v=U(j+1)(2)/U(j+1)(0);
    sound_speed=eos->c(U(j+1)(0),e(U(j+1)(0),U(j+1)(1),U(j+1)(2),U(j+1)(3)));
    if(sound_speed == 0) cavitation=true;
    lambda_right(0)=v-sound_speed;
    lambda_right(1)=lambda_right(2)=v;
    lambda_right(3)=v+sound_speed;
        
    // eigenvalues in the center - at flux j
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T E=(U(j)(3)+U(j+1)(3))/2;
    v=rho_v/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    sound_speed=eos->c(rho,internal_energy);
    if(sound_speed == 0) cavitation=true;
    lambda(0)=v-sound_speed;
    lambda(1)=lambda(2)=v;
    lambda(3)=v+sound_speed;

    return (!cavitation); //cavitation --> loss of hyperbolicity, else well defined eigensystem
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for G(U) at flux j between points j and j+1
template<class T> void EULER_2D_EIGENSYSTEM_G<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux j
    T rho=(U(j)(0)+U(j+1)(0))/2;
    T rho_u=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(2)+U(j+1)(2))/2;
    T E=(U(j)(3)+U(j+1)(3))/2;
    T v=rho_v/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    T sound_speed=eos->c(rho,internal_energy);
    T p=eos->p(rho,internal_energy);
        
    T u=rho_u/rho;
    T q2=sqr(u)+sqr(v);
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
    T v_times_c=v*sound_speed;
                    
    L(0,0)=b2_over_2+v_over_2c;
    L(0,1)=-b1_over_2_times_u;
    L(0,2)=-b1_over_2_times_v-one_over_2c;
    L(0,3)=b1_over_2;
    L(1,0)=h-q2;
    L(1,1)=u;
    L(1,2)=v;
    L(1,3)=-1;
    L(2,0)=-u;
    L(2,1)=1;
    L(2,2)=0;
    L(2,3)=0;
    L(3,0)=b2_over_2-v_over_2c;
    L(3,1)=-b1_over_2_times_u;
    L(3,2)=-b1_over_2_times_v+one_over_2c;
    L(3,3)=b1_over_2;
    
    R(0,0)=1;
    R(0,1)=u;
    R(0,2)=v-sound_speed;
    R(0,3)=h-v_times_c;
    R(1,0)=b1;
    R(1,1)=b1*u;
    R(1,2)=b1*v;
    R(1,3)=b1*h-1;
    R(2,0)=0;
    R(2,1)=1;
    R(2,2)=0;
    R(2,3)=u;
    R(3,0)=1;
    R(3,1)=u;
    R(3,2)=v+sound_speed;
    R(3,3)=h+v_times_c;
}  
//#####################################################################
}
#endif

