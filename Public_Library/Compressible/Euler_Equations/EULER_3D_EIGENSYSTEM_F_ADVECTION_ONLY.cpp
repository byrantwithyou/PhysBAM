//#####################################################################
// Copyright 2007-2009, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/sqr.h>
#include <Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(U_clamped){
        for(int i=-3;i<m+3;i++){
            T one_over_U_1=(T)1/U(i)(0);
            T u=U(i)(1)*one_over_U_1,v=U(i)(2)*one_over_U_1,w=U(i)(3)*one_over_U_1;
            F(i)(0)=(*U_clamped)(i)(1);       // rho_clamped*u
            F(i)(1)=(*U_clamped)(i)(1)*u;     // rho_clamped*u^2
            F(i)(2)=(*U_clamped)(i)(1)*v;     // rho_clamped*u*v
            F(i)(3)=(*U_clamped)(i)(1)*w;     // rho_clamped*u*w
            F(i)(4)=(*U_clamped)(i)(4)*u;}}   // E_from_rho_clamped*u
    else{
        for(int i=-3;i<m+3;i++){
            T one_over_U_1=(T)1/U(i)(0);
            T u=U(i)(1)*one_over_U_1,v=U(i)(2)*one_over_U_1,w=U(i)(3)*one_over_U_1;
            F(i)(0)=U(i)(1);       // rho*u
            F(i)(1)=U(i)(1)*u;     // rho*u^2
            F(i)(2)=U(i)(1)*v;     // rho*u*v
            F(i)(3)=U(i)(1)*w;     // rho*u*w
            F(i)(4)=U(i)(4)*u;}}   // E*u
}
//#####################################################################
// Function Flux_Divided_By_Velocity
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    F=(U_clamped)?*U_clamped:U;
}
//#####################################################################
// Function Get_Face_Velocity_Component
//#####################################################################
template<class T> T EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U)
{
    T average_u;
    if(use_standard_average) average_u=(U(face_index)(1)/U(face_index)(0)+U(face_index+1)(1)/U(face_index+1)(0))*(T).5;
    else average_u=(U(face_index)(1)+U(face_index+1)(1))/(U(face_index)(0)+U(face_index+1)(0));
    return average_u;
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    T average_u;
    if(use_standard_average) average_u=(U(face_index)(1)/U(face_index)(0)+U(face_index+1)(1)/U(face_index+1)(0))*(T).5;
    else average_u=(U(face_index)(1)+U(face_index+1)(1))/(U(face_index)(0)+U(face_index+1)(0));

    if(U_clamped){
        for(int i=range.x;i<range.y;i++){
            F(i)(0)=(*U_clamped)(i)(0)*average_u;
            F(i)(1)=(*U_clamped)(i)(1)*average_u;
            F(i)(2)=(*U_clamped)(i)(2)*average_u;
            F(i)(3)=(*U_clamped)(i)(3)*average_u;
            F(i)(4)=(*U_clamped)(i)(4)*average_u;}}
    else{
        for(int i=range.x;i<range.y;i++){
            F(i)(0)=U(i)(0)*average_u;
            F(i)(1)=U(i)(1)*average_u;
            F(i)(2)=U(i)(2)*average_u;
            F(i)(3)=U(i)(3)*average_u;
            F(i)(4)=U(i)(4)*average_u;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    return abs(U_cell(1)/U_cell(0));
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool cavitation=false;

    // eigenvalues on the left - at point i
    T u=U(i)(1)/U(i)(0);
    lambda_left(0)=lambda_left(1)=lambda_left(2)=lambda_left(3)=lambda_left(4)=u;

    // eigenvalues on the right - at point i+1
    u=U(i+1)(1)/U(i+1)(0);
    lambda_right(0)=lambda_right(1)=lambda_right(2)=lambda_right(3)=lambda_right(4)=u;

    // eigenvalues in the center - at flux i
    T rho=(U(i)(0)+U(i+1)(0))/2;
    T rho_u=(U(i)(1)+U(i+1)(1))/2;
    u=rho_u/rho;
    lambda(0)=lambda(1)=lambda(2)=lambda(3)=lambda(4)=u;

    return (!cavitation); //cavitation --> loss of hyperbolicity, else well defined eigensystem
}
template class EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<double>;
template class EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<float>;
}
