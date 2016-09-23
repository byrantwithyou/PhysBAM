//#####################################################################
// Copyright 2007-2009, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/sqr.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Function Flux
//#####################################################################
// H(U) for ij in (-2,mn+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Flux(const int mn,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(U_clamped){
        for(int ij=-3;ij<mn+3;ij++){
            T u=U(ij)(1)/U(ij)(0),v=U(ij)(2)/U(ij)(0),w=U(ij)(3)/U(ij)(0);
            H(ij)(0)=(*U_clamped)(ij)(3);         // rho_clamped*w
            H(ij)(1)=(*U_clamped)(ij)(3)*u;       // rho_clamped*u*w
            H(ij)(2)=(*U_clamped)(ij)(3)*v;       // rho_clamped*v*w
            H(ij)(3)=(*U_clamped)(ij)(3)*w;       // rho_clamped*w^2
            H(ij)(4)=(*U_clamped)(ij)(4)*w;}}     // E_from_rho_clamped*w
    else{
        for(int ij=-3;ij<mn+3;ij++){
            T u=U(ij)(1)/U(ij)(0),v=U(ij)(2)/U(ij)(0),w=U(ij)(3)/U(ij)(0);
            H(ij)(0)=U(ij)(3);         // rho*w
            H(ij)(1)=U(ij)(3)*u;       // rho*u*w
            H(ij)(2)=U(ij)(3)*v;       // rho*v*w
            H(ij)(3)=U(ij)(3)*w;       // rho*w^2
            H(ij)(4)=U(ij)(4)*w;}}     // E*w
}
//#####################################################################
// Function Flux_Divided_By_Velocity
//#####################################################################
// H(U) for ij in (-2,mn+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Flux_Divided_By_Velocity(const int mn,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    H=(U_clamped)?*U_clamped:U;
}
//#####################################################################
// Function Get_Face_Velocity_Component
//#####################################################################
template<class T> T EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U)
{
    T average_w;
    if(use_standard_average) average_w=(U(face_index)(3)/U(face_index)(0)+U(face_index+1)(3)/U(face_index+1)(0))*(T).5;
    else average_w=(U(face_index)(3)+U(face_index+1)(3))/(U(face_index)(0)+U(face_index+1)(0));
    return average_w;
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    T average_w;
    if(use_standard_average) average_w=(U(face_index)(3)/U(face_index)(0)+U(face_index+1)(3)/U(face_index+1)(0))*(T).5;
    else average_w=(U(face_index)(3)+U(face_index+1)(3))/(U(face_index)(0)+U(face_index+1)(0));

    if(U_clamped){
        for(int i=range.x;i<range.y;i++){
            H(i)(0)=(*U_clamped)(i)(0)*average_w;
            H(i)(1)=(*U_clamped)(i)(1)*average_w;
            H(i)(2)=(*U_clamped)(i)(2)*average_w;
            H(i)(3)=(*U_clamped)(i)(3)*average_w;
            H(i)(4)=(*U_clamped)(i)(4)*average_w;}}
    else{
        for(int i=range.x;i<range.y;i++){
            H(i)(0)=U(i)(0)*average_w;
            H(i)(1)=U(i)(1)*average_w;
            H(i)(2)=U(i)(2)*average_w;
            H(i)(3)=U(i)(3)*average_w;
            H(i)(4)=U(i)(4)*average_w;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for H(U) at point cell
template<class T> T EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    return abs(U_cell(3)/U_cell(0));
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for H(U) at flux ij and and at points ij and ij+1
template<class T> bool EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int ij,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point ij
    T w=U(ij)(3)/U(ij)(0);
    lambda_left(0)=lambda_left(1)=lambda_left(2)=lambda_left(3)=lambda_left(4)=w;

    // eigenvalues on the right - at point ij+1
    w=U(ij+1)(3)/U(ij+1)(0);
    lambda_right(0)=lambda_right(1)=lambda_right(2)=lambda_right(3)=lambda_right(4)=w;

    // eigenvalues in the center - at flux ij
    T rho=(U(ij)(0)+U(ij+1)(0))/2;
    T rho_w=(U(ij)(3)+U(ij+1)(3))/2;
    w=rho_w/rho;
    lambda(0)=lambda(1)=lambda(2)=lambda(3)=lambda(4)=w;

    return (!cavitation); //cavitation --> loss of hyperbolicity, else well defined eigensystem
}
template class EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<double>;
template class EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<float>;
}
