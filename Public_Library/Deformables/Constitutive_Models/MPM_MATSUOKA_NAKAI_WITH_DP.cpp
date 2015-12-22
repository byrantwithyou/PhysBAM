#if 0
//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <cmath>
#include <Deformables/Constitutive_Models/MPM_MATSUOKA_NAKAI_WITH_DP.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{
// [1] Mast, Carter M. Modeling landslide-induced flow interactions with structures using the material point method. Diss. University of Washington, 2013.
//#####################################################################
// Function Compute_Invariants_Helper
//#####################################################################
template<class T> static void
Compute_Invariants_Helper(const VECTOR<T,2>& sigma,T& I1,T& I2,T& I3)
{
    I1=sigma.Sum();
    I2=sigma.Product();
    I3=I2;
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
template<class TV> MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
MPM_MATSUOKA_NAKAI_WITH_DP()
    :a0(20),a1(0),a2(0),a3(70),a4(20),a5(50) // Case iii. of [1] pg. 57
{
}
//#####################################################################
// Function Set_Lame_Constants_And_F_Elastic
//#####################################################################
template<class TV> void MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe)
{
    D=SYMMETRIC_MATRIX<T,d>::Unit_Matrix(-lambda/(2*d*mu*lambda+4*mu*mu))+SYMMETRIC_MATRIX<T,d>::Identity_Matrix()/(2*mu);

    strain_trial=log(Fe.x);
    tau_trial=D.Inverse_Times(strain_trial);
}
//#####################################################################
// Function Set_Plastic_Deformation_Lambda
//#####################################################################
template<class TV> void MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Set_Plastic_Deformation_Lambda(T plastic_def)
{
    this->plastic_def=plastic_def;
    phi_F=Get_phi_F(tau_trial,plastic_def);
    phi_cs=Get_phi_cs(plastic_def);
    T sqrsin=sqr(sin(phi_F));
    kappa_F=(sqrsin-9)/(sqrsin-1);
    psi_G=atan(tan(phi_F)-tan(phi_cs));
    T sin_psi_G=sin(psi_G);
    rho_G=sqrt(2.0/3.0)*2*sin_psi_G/(3-sin_psi_G);
}
//#####################################################################
// Function Yield_Function
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Yield_Function(const TV& tau) const
{
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    return 6*(kappa_F*I3-I1*I2);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
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
    if(residual.Magnitude_Squared()<=tolerance*tolerance){
        delta_gamma_final=x(d);
        return true;}
    else{
        delta_gamma_final=0;
        return false;}
}
//#####################################################################
// Function Get_Jacobian
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m+1> MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Get_Jacobian(const TVP1& x) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    MATRIX<T,d+1> ret;
    ret.Set_Submatrix(0,0,D+delta_gamma*Yield_Function_Hessian(tau));
    TV f_hat_vec=Yield_Function_Derivative(tau);
    TV g_hat_vec=Plastic_Flow(tau);
    for(int i=0;i<d;i++){
        ret(i,d)=g_hat_vec(i);
        ret(d,i)=f_hat_vec(i);}

    return ret;
}
//#####################################################################
// Function Yield_Function_Derivative_Helper
//#####################################################################
template<class T> static VECTOR<T,2>
Yield_Function_Derivative_Helper(const VECTOR<T,2>& tau,T kappa_F)
{
    VECTOR<T,2> tr=tau.Reversed();
    return (T)6*(kappa_F*tr-(T)2*tau.Product()-tr*tr);
}
//#####################################################################
// Function Yield_Function_Derivative_Helper
//#####################################################################
template<class T> static VECTOR<T,3>
Yield_Function_Derivative_Helper(const VECTOR<T,3>& tau,T kappa_F)
{
    typedef VECTOR<T,3> TV;
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    VECTOR<T,3> s=tau-tau.Sum()/3;
    T s_norm_squared=s.Magnitude_Squared();
    TV oo=TV::All_Ones_Vector();
    return 6*kappa_F*I3*Inverse(tau)+(3*s_norm_squared-8*sqr(I1))*oo+6*I1*tau;
}
//#####################################################################
// Function Yield_Function_Derivative
//#####################################################################
template<class TV> TV MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Yield_Function_Derivative(const TV& tau) const
{
    return Yield_Function_Derivative_Helper(tau,kappa_F);
}
//#####################################################################
// Function Yield_Function_Hessian_Helper
//#####################################################################
template<class T> static SYMMETRIC_MATRIX<T,2>
Yield_Function_Hessian_Helper(const VECTOR<T,2>& tau,T kappa_F)
{
    return SYMMETRIC_MATRIX<T,2>(-(T)12*tau.y,(T)6*(kappa_F-(T)2*(tau.x+tau.y)),-(T)12*tau.x);
}
//#####################################################################
// Function Yield_Function_Hessian_Helper
//#####################################################################
template<class T> static SYMMETRIC_MATRIX<T,3>
Yield_Function_Hessian_Helper(const VECTOR<T,3>& tau,T kappa_F)
{
    typedef VECTOR<T,3> TV;
    T I1,I2,I3;
    Compute_Invariants_Helper(tau,I1,I2,I3);
    TV oo=TV::All_Ones_Vector();
    SYMMETRIC_MATRIX<T,3>  op=SYMMETRIC_MATRIX<T,3>::Outer_Product(Inverse(tau));
    DIAGONAL_MATRIX<T,3> ov(sqr(Inverse(tau)));
    SYMMETRIC_MATRIX<T,3> ss=(MATRIX<T,3>::Outer_Product(oo,tau)).Twice_Symmetric_Part();
    SYMMETRIC_MATRIX<T,3> ao=SYMMETRIC_MATRIX<T,3>::Outer_Product(oo);
    return (T)6*(kappa_F*I3*(op-ov)+ss+I1-3*I1*ao);
}
//#####################################################################
// Function Yield_Function_Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Yield_Function_Hessian(TV tau) const
{
    return Yield_Function_Hessian_Helper(tau,kappa_F);
}
//#####################################################################
// Function Plastic_Flow
//#####################################################################
template<class TV> TV MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
Plastic_Flow(const TV& tau) const
{
    return (tau-tau.Sum()/d).Normalized()+rho_G; 
}
//#####################################################################
// Function Get_Residual
//#####################################################################
template<class TV> typename MPM_MATSUOKA_NAKAI_WITH_DP<TV>::TVP1 MPM_MATSUOKA_NAKAI_WITH_DP<TV>::
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
template class MPM_MATSUOKA_NAKAI_WITH_DP<VECTOR<float,2>>;
template class MPM_MATSUOKA_NAKAI_WITH_DP<VECTOR<float,3>>;
template class MPM_MATSUOKA_NAKAI_WITH_DP<VECTOR<double,2>>;
template class MPM_MATSUOKA_NAKAI_WITH_DP<VECTOR<double,3>>;
}
#endif
