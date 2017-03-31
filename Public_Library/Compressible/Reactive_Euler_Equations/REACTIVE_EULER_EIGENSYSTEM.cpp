//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER_EIGENSYSTEM.h>
namespace PhysBAM{
template<class T,class TV> inline VECTOR<T,TV::m+3>
Stack(T r,TV u,T e,T y)
{
    return u.Prepend(r).Append(e).Append(y);
}
template<class T,class TV> inline VECTOR<T,TV::m+3>
Stack_Fix(T r,TV u,T e,T y,int a,T p)
{
    u(a)+=p;
    return u.Prepend(r).Append(e).Append(y);
}
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class TV> void REACTIVE_EULER_EIGENSYSTEM<TV>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-3;i<m+3;i++){
        TV_DIMENSION UU=U(i);
        TV u=UU.template Slice<1,TV::m>()/UU(0);
        T Y=UU(d-1),E=UU(d-2);
        T U_a=UU(1+a);
        T p=eos->p(UU(0),EULER<TV>::e(UU.Remove_Index(d-1)),Y);
        F(i)=Stack_Fix(U_a,U_a*u,(E+p)*u(a),Y*u(a),a,p);}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class TV> typename TV::SCALAR REACTIVE_EULER_EIGENSYSTEM<TV>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(1+a)/U_cell(0);
    T Y=U_cell(d-1)/U_cell(0);
    T sound_speed=eos->c(U_cell(0),EULER<TV>::e(U_cell.Remove_Index(d-1)),Y);
    return maxabs(u-sound_speed,u+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class TV> bool REACTIVE_EULER_EIGENSYSTEM<TV>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)
{
    bool is_hyperbolic=true;
    auto eig=[this,&is_hyperbolic](TV_DIMENSION UU)
        {
            T u=UU(1+a)/UU(0);
            T Y=UU(d-1)/UU(0);
            TV_DIMENSION ev;
            ev.Fill(u);
            T sound_speed=eos->c(UU(0),EULER<TV>::e(UU.Remove_Index(d-1)),Y);
            if(sound_speed==0) is_hyperbolic=false;
            ev(0)-=sound_speed;
            ev(d-1)+=sound_speed;
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
template<class TV> void REACTIVE_EULER_EIGENSYSTEM<TV>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    TV_DIMENSION UU=(T).5*(U(i)+U(i+1));
    T rho=UU(0);
    TV rho_u=UU.template Slice<1,TV::m>();
    T E=UU(d-1);
    T rho_Y=UU(d-1);
    T Y=rho_Y/rho;
    T internal_energy=EULER<TV>::e(UU.Remove_Index(d-1));
    T sound_speed=eos->c(rho,internal_energy,Y);
    T p=eos->p(rho,internal_energy,Y);

    TV u=rho_u/rho;
    T q2=u.Magnitude_Squared();
    T h=(E+p)/rho;
    T b1=eos->p_e(rho,internal_energy,Y)/(rho*sqr(sound_speed));
    T b2=(T).5*(1+b1*(q2-h));
    T z=eos->e_not1-eos->e_not2;
    T b3=b1*Y*z;

    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T u_over_2c=u(a)*one_over_2c;
    T b1_over_2=b1/2;
    T b3_over_2=b3/2;
    TV u_b1_over_2=b1_over_2*u;
    T u_times_c=u(a)*sound_speed;

    L.Set_Row(0,Stack_Fix(b2+u_over_2c+b3_over_2,-u_b1_over_2,b1_over_2,-b1_over_2*z,a,-one_over_2c));
    L.Set_Row(d-1,Stack_Fix(b2-u_over_2c+b3_over_2,-u_b1_over_2,b1_over_2,-b1_over_2*z,a,one_over_2c));
    L.Set_Row(1,Stack(1-b2*2-b3,b1*u,-b1,b1*z));
    L.Set_Row(d-2,Stack(-Y,TV(),(T)0,(T)1));
    R.Set_Row(0,Stack_Fix((T)1,u,h-u_times_c,Y,a,-sound_speed));
    R.Set_Row(d-1,Stack_Fix((T)1,u,h+u_times_c,Y,a,sound_speed));
    R.Set_Row(1,Stack((T)1,u,h-1/b1,Y));
    R.Set_Row(d-2,Stack((T)0,TV(),z,(T)1));
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
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<float,1> >;
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<float,2> >;
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<float,3> >;
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<double,1> >;
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<double,2> >;
template class REACTIVE_EULER_EIGENSYSTEM<VECTOR<double,3> >;
}
