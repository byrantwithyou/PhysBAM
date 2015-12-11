//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/MPM_DRUCKER_PRAGER_HARDENING.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{

//#####################################################################
// Function Set_Lame_Constants_And_F_Elastic
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER_HARDENING<TV>::
Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe)
{
    this->mu=mu;
    this->lambda=lambda;
    D=SYMMETRIC_MATRIX<T,d>::Unit_Matrix(-lambda/(2*d*mu*lambda+4*mu*mu))+SYMMETRIC_MATRIX<T,d>::Identity_Matrix()/(2*mu);

    strain_trial=log(abs(Fe.x));
    tau_trial=D.Inverse_Times(strain_trial);
}
//#####################################################################
// Function Set_Plastic_Deformation_Lambda
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER_HARDENING<TV>::
Set_Plastic_Deformation_Lambda(T plastic_def)
{
    this->plastic_def=plastic_def;
    T phi_F=M_PI/180*(a0+(a1*plastic_def-a4)*exp(-a3*plastic_def));
    T sin_phi_F=sin(phi_F);
    rho_F=sqrt(2.0/3.0)*2*sin_phi_F/(3-sin_phi_F);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_DRUCKER_PRAGER_HARDENING<TV>::
Project_Stress(int max_iterations, T tolerance)
{
    PHYSBAM_ASSERT(plastic_def>=0 && "Plastic deformation lambda has not been set!");
    if(direct_solution){
        T g=1/(d*lambda+2*mu);
        T b=1/(2*mu);
        T k=strain_trial.Sum();
        TV sh=strain_trial-k/d*TV::All_Ones_Vector();
        T q=sh.Magnitude();
        if(q==0){
           tau_final=TV();
           return true;} 
        else{
            T dg=(g*q+b*(rho_F*k))/g;
            T m=k/(d*g);
            T n=(q-dg)/(b*q);
            tau_final=m*TV::All_Ones_Vector()+n*sh;
            T e=n*q;
            if(e<0) tau_final=TV();
            T e3=(tau_final-tau_final.Sum()/d).Magnitude()+rho_F*tau_final.Sum();
            if(abs(e3)>1e-8) LOG::printf("BAD PARTICLE: yield function value:%g\n",e3);
            return e>=0;}}
    else{
        if(tau_trial.Sum()>0){
            tau_final*=0;
            delta_gamma_final=0;
            PHYSBAM_ASSERT(abs(Yield_Function_Final())<1e-2);
            return true;}
        else{
            TVP1 x,residual;
            x.Set_Subvector(0,tau_trial);
            residual=Get_Residual(x,strain_trial,D);
            int k=0;
            while(k++<max_iterations && residual.Magnitude_Squared()>tolerance*tolerance){
                x-=Get_Jacobian(x).Inverse_Times(residual);
                residual=Get_Residual(x,strain_trial,D);}
            x.Get_Subvector(0,tau_final);
            PHYSBAM_ASSERT(abs(Yield_Function_Final())<1e-2);
            delta_gamma_final=x(d);
            //LOG::printf("residual(%d):%g\n",k,residual.Magnitude());
            //LOG::printf("tau_final=%g\n",tau_final);
            //LOG::printf("F(tau_final)=%g\n",Yield_Function(tau_final));
            //LOG::printf("D*tau_final=%g\n",D*tau_final);
            return residual.Magnitude_Squared()<=tolerance*tolerance;}}
}
//#####################################################################
// Function Get_Jacobian
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m+1> MPM_DRUCKER_PRAGER_HARDENING<TV>::
Get_Jacobian(const TVP1& x) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    TV f_hat_vec=Yield_Function_Derivative(tau);
    TV g_hat_vec=(tau-tau.Sum()/d).Normalized();
    MATRIX<T,d+1> ret;
    ret.Set_Submatrix(0,0,D+delta_gamma*Plastic_Potential_Hessian(tau));
    for(int i=0;i<d;i++){
        ret(i,d)=g_hat_vec(i);
        ret(d,i)=f_hat_vec(i);}

    return ret;
}
//#####################################################################
// Function Yield_Function_Derivative
//#####################################################################
// s:s = (1-(T)1.0/d)*I1*I1-2*I2
template<class TV> TV MPM_DRUCKER_PRAGER_HARDENING<TV>::
Yield_Function_Derivative(const TV& tau) const
{
    return (tau-tau.Sum()/d).Normalized()+rho_F; 
}
//#####################################################################
// Function Plastic_Potential_Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_DRUCKER_PRAGER_HARDENING<TV>::
Plastic_Potential_Hessian(TV tau) const {
    TV s=tau-tau.Sum()/d;
    T s_norm=s.Normalize();
    return ((T)1-SYMMETRIC_MATRIX<T,d>::Unit_Matrix((T)1/d)-SYMMETRIC_MATRIX<T,d>::Outer_Product(s))/s_norm;
}
//#####################################################################
// Function Get_Residual
//#####################################################################
template<class TV> typename MPM_DRUCKER_PRAGER_HARDENING<TV>::TVP1 MPM_DRUCKER_PRAGER_HARDENING<TV>::
Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    TV strain_residual=D*tau-strain_trial+delta_gamma*(tau-tau.Sum()/d).Normalized();
    return strain_residual.Append(Yield_Function(tau));
}
}
namespace PhysBAM{
template class MPM_DRUCKER_PRAGER_HARDENING<VECTOR<float,2>>;
template class MPM_DRUCKER_PRAGER_HARDENING<VECTOR<float,3>>;
template class MPM_DRUCKER_PRAGER_HARDENING<VECTOR<double,2>>;
template class MPM_DRUCKER_PRAGER_HARDENING<VECTOR<double,3>>;
}
