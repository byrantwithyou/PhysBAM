//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{
template<class T,class TV> inline VECTOR<T,TV::m+2>
Stack(T r,TV u,T e)
{
    return u.Prepend(r).Append(e);
}
template<class T,class TV> inline VECTOR<T,TV::m+2>
Stack_Fix(T r,TV u,T e,int a,T p)
{
    u(a)+=p;
    return u.Prepend(r).Append(e);
}
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class TV> void EULER_EIGENSYSTEM<TV>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    if(only_pressure_flux){
        for(int i=-3;i<m+3;i++){
            T p=eos->p(U(i)(0),EULER<TV>::e(U(i)));
            TV_DIMENSION f;
            f(1+a)=p;
            F(i)=f;}
        return;}

    const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_cl=U_clamped?*U_clamped:U;
    for(int i=-3;i<m+3;i++){
        TV u=U(i).template Slice<1,TV::m>()/U(i)(0);
        T p=only_advection?0:eos->p(U(i)(0),EULER<TV>::e(U(i)));
        T u_cl=U_cl(i)(1+a);
        F(i)=Stack_Fix(u_cl,u_cl*u,(U_cl(i)(d-1)+p)*u(a),a,p);}
}
//#####################################################################
// Function Flux_Divided_By_Velocity
//#####################################################################
template<class TV> void EULER_EIGENSYSTEM<TV>::
Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    PHYSBAM_ASSERT(only_advection);
    F=U_clamped?*U_clamped:U;
}
//#####################################################################
// Function Get_Face_Velocity_Component
//#####################################################################
template<class TV> typename TV::SCALAR EULER_EIGENSYSTEM<TV>::
Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U)
{
    PHYSBAM_ASSERT(only_advection);
    T average_u;
    if(use_standard_average) average_u=(U(face_index)(1)/U(face_index)(0)+U(face_index+1)(1)/U(face_index+1)(0))*(T).5;
    else average_u=(U(face_index)(1)+U(face_index+1)(1))/(U(face_index)(0)+U(face_index+1)(0));
    return average_u;
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class TV> void EULER_EIGENSYSTEM<TV>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(only_pressure_flux){
        for(int i=range.x;i<range.y;i++){
            T p=eos->p(U(i)(0),EULER<TV>::e(U(i)));
            TV_DIMENSION f;
            f(1+a)=p;
            F(i)=f;}
        return;}

    T average_u;
    if(use_standard_average) average_u=(U(face_index)(1+a)/U(face_index)(0)+U(face_index+1)(1+a)/U(face_index+1)(0))*(T).5;
    else average_u=(U(face_index)(1+a)+U(face_index+1)(1+a))/(U(face_index)(0)+U(face_index+1)(0));

    const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_cl=U_clamped?*U_clamped:U;
    for(int i=range.x;i<range.y;i++){
        TV_DIMENSION f=U_cl(i)*average_u;
        if(!only_advection){
            T p=eos->p(U(i)(0),EULER<TV>::e(U(i)));
            f(1+a)+=p;
            f(d-1)+=p*average_u;}
        F(i)=f;}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class TV> typename TV::SCALAR EULER_EIGENSYSTEM<TV>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(1+a)/U_cell(0);
    if(only_advection) return abs(u);
    T sound_speed=eos->c(U_cell(0),EULER<TV>::e(U_cell));
    return maxabs(u-sound_speed,u+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class TV> bool EULER_EIGENSYSTEM<TV>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)
{
    bool is_hyperbolic=true;
    auto eig=[this,&is_hyperbolic](TV_DIMENSION UU)
        {
            T u=UU(1+a)/UU(0);
            TV_DIMENSION ev;
            ev.Fill(u);
            if(!only_advection){
                T sound_speed=eos->c(UU(0),EULER<TV>::e(UU));
                if(sound_speed==0) is_hyperbolic=false;
                ev(0)-=sound_speed;
                ev(d-1)+=sound_speed;}
            return ev;
        };

    lambda_left=eig(U(i)); // eigenvalues on the left - at point i
    lambda_right=eig(U(i+1)); // eigenvalues on the right - at point i+1
    lambda=eig((T).5*(U(i)+U(i+1))); // eigenvalues in the center - at flux i

    return is_hyperbolic; // eigensystem is well defined (hyperbolic)
}
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class TV> void EULER_EIGENSYSTEM<TV>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    PHYSBAM_ASSERT(!only_advection);
    // eigensystem in the center - at flux i
    TV_DIMENSION UU=(T).5*(U(i)+U(i+1));
    T rho=UU(0);
    TV rho_u=UU.template Slice<1,TV::m>();
    T E=UU(d-1);
    T internal_energy=EULER<TV>::e(UU);
    T sound_speed=eos->c(rho,internal_energy);
    T p=only_advection?0:eos->p(rho,internal_energy);

    TV u=rho_u/rho;
    T q2=u.Magnitude_Squared();
    T h=(E+p)/rho;
    T b1=eos->p_e(rho,internal_energy)/(rho*sqr(sound_speed));
    T b2=(T).5*(1+b1*(q2-h));

    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T u_over_2c=u(a)*one_over_2c;
    T b1_over_2=b1/2;
    TV u_b1_over_2=b1_over_2*u;
    T u_times_c=u(a)*sound_speed;

    L.Set_Row(0,Stack_Fix(b2+u_over_2c,-u_b1_over_2,b1_over_2,a,-one_over_2c));
    L.Set_Row(d-1,Stack_Fix(b2-u_over_2c,-u_b1_over_2,b1_over_2,a,one_over_2c));
    L.Set_Row(1,Stack(h-q2,u,-(T)1));
    R.Set_Row(0,Stack_Fix((T)1,u,h-u_times_c,a,-sound_speed));
    R.Set_Row(d-1,Stack_Fix((T)1,u,h+u_times_c,a,sound_speed));
    R.Set_Row(1,Stack(b1,b1*u,b1*h-1));
    for(int i=0,j=2;i<TV::m;i++){
        if(i==a) continue;
        TV_DIMENSION l,r;
        l(0)=u(i);
        l(1+i)=-1;
        L.Set_Row(j,l);
        r(1+i)=-1;
        r(d-1)=-u(i);
        R.Set_Row(j,r);
        j++;}
}
//#####################################################################
// Function All_Eigenvalues_Same
//#####################################################################
template<class TV> bool EULER_EIGENSYSTEM<TV>::
All_Eigenvalues_Same()
{
    return only_advection;
}
template class EULER_EIGENSYSTEM<VECTOR<double,1> >;
template class EULER_EIGENSYSTEM<VECTOR<double,2> >;
template class EULER_EIGENSYSTEM<VECTOR<double,3> >;
template class EULER_EIGENSYSTEM<VECTOR<float,1> >;
template class EULER_EIGENSYSTEM<VECTOR<float,2> >;
template class EULER_EIGENSYSTEM<VECTOR<float,3> >;
}
