//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <cmath>
#include <Deformables/Constitutive_Models/MPM_MATSUOKA_NAKAI.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{
//#####################################################################
// Function Compute_Invariants_Helper
//#####################################################################
template<class T> static void
Compute_Invariants_Helper(const VECTOR<T,2>& sigma,T& I1,T& I2,T& I3)
{
    I1=sigma.Sum();
    I2=sigma.Product();
    I3=T();
}
//#####################################################################
// Function Compute_Invariants_Helper
//#####################################################################
template<class T> static void
Compute_Invariants_Helper(const VECTOR<T,3>& sigma,T& I1,T& I2,T& I3)
{
    I1=sigma.Sum();
    I2=sigma.x*sigma.y+sigma.x*sigma.z+sigma.y*sigma.z;
    I3=sigma.Product();
}
//#####################################################################
// Function Constructor
//#####################################################################
template<class TV> MPM_MATSUOKA_NAKAI<TV>::
MPM_MATSUOKA_NAKAI(typename TV::SCALAR fa,typename TV::SCALAR c):kappa((sqr(sin(fa))-(T)9)/(sqr(sin(fa))-(T)1)),c(c)
{
    if(c!=0) PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Set_Lame_Constants_And_F_Elastic
//#####################################################################
template<class TV> void MPM_MATSUOKA_NAKAI<TV>::
Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe)
{
    D=SYMMETRIC_MATRIX<T,d>::Unit_Matrix(-lambda/(2*d*mu*lambda+4*mu*mu))+SYMMETRIC_MATRIX<T,d>::Identity_Matrix()/(2*mu);

    strain_trial=log(Fe.x);
    tau_trial=D.Inverse_Times(strain_trial);
}
//#####################################################################
// Function Yield_Function
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MATSUOKA_NAKAI<TV>::
Yield_Function() const
{
    T I1,I2,I3;
    Compute_Invariants_Helper(tau_trial,I1,I2,I3);
    return 6*(kappa*I3-I1*I2);
}
//#####################################################################
// Function Yield_Function
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MATSUOKA_NAKAI<TV>::
Yield_Function(const TV& tau) const
{
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    return 6*(kappa*I3-I1*I2);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_MATSUOKA_NAKAI<TV>::
Project_Stress(int max_iterations, T tolerance)
{
    TVP1 x,residual;
    x.Set_Subvector(0,tau_trial);
    residual=Get_Residual(x,strain_trial,D);
    int k=0;
    while(k++<max_iterations && residual.Magnitude_Squared()>tolerance*tolerance){
        x-=Get_Jacobian(x).Inverse_Times(residual);
        residual=Get_Residual(x,strain_trial,D);}
    x.Get_Subvector(0,tau_final);
    return residual.Magnitude_Squared()<=tolerance*tolerance;
}
//#####################################################################
// Function Get_Jacobian
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m+1> MPM_MATSUOKA_NAKAI<TV>::
Get_Jacobian(const TVP1& x) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    TV f_hat_vec=Yield_Function_Derivative(tau);
    MATRIX<T,d+1> ret;
    ret.Set_Submatrix(0,0,D+delta_gamma*Yield_Function_Hessian(tau));
    for(int i=0;i<d;i++){
        ret(i,d)=f_hat_vec(i);
        ret(d,i)=f_hat_vec(i);}

    return ret;
}
//#####################################################################
// Function Yield_Function_Derivative
//#####################################################################
template<class TV> TV MPM_MATSUOKA_NAKAI<TV>::
Yield_Function_Derivative(const TV& tau) const
{
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    TV s=tau-tau.Sum()/d;
    T s_norm_squared=s.Magnitude_Squared();
    TV oo=TV::All_Ones_Vector();
    return 6*kappa*I3*Inverse(tau)+(3*s_norm_squared-8*sqr(I1))*oo+6*I1*tau;
}
//#####################################################################
// Function Yield_Function_Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_MATSUOKA_NAKAI<TV>::
Yield_Function_Hessian(TV tau) const {
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    TV oo=TV::All_Ones_Vector();
    SYMMETRIC_MATRIX<T,d>  op=SYMMETRIC_MATRIX<T,d>::Outer_Product(Inverse(tau));
    DIAGONAL_MATRIX<T,d> ov(sqr(Inverse(tau)));
    SYMMETRIC_MATRIX<T,d> ss=(MATRIX<T,d>::Outer_Product(oo,tau)).Twice_Symmetric_Part();
    SYMMETRIC_MATRIX<T,d> ao=SYMMETRIC_MATRIX<T,d>::Outer_Product(oo);
    return (T)6*(kappa*I3*(op-ov)+ss+I1-3*I1*ao);
}
//#####################################################################
// Function Get_Residual
//#####################################################################
template<class TV> typename MPM_MATSUOKA_NAKAI<TV>::TVP1 MPM_MATSUOKA_NAKAI<TV>::
Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    TV strain_residual=D*tau-strain_trial+delta_gamma*Yield_Function_Derivative(tau);
    return strain_residual.Append(Yield_Function(tau));
}
}
namespace PhysBAM{
template class MPM_MATSUOKA_NAKAI<VECTOR<float,2>>;
template class MPM_MATSUOKA_NAKAI<VECTOR<float,3>>;
template class MPM_MATSUOKA_NAKAI<VECTOR<double,2>>;
template class MPM_MATSUOKA_NAKAI<VECTOR<double,3>>;
}
