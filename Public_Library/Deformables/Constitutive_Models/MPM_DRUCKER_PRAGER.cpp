//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/MPM_DRUCKER_PRAGER.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{

//#####################################################################
// Function Set_Lame_Constants_And_F_Elastic
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe)
{
    this->mu=mu;
    this->lambda=lambda;
    D=SYMMETRIC_MATRIX<T,d>::Unit_Matrix(-lambda/(2*d*mu*lambda+4*mu*mu))+SYMMETRIC_MATRIX<T,d>::Identity_Matrix()/(2*mu);

    strain_trial=log(Fe.x);
    tau_trial=D.Inverse_Times(strain_trial);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_DRUCKER_PRAGER<TV>::
Project_Stress(int max_iterations, T tolerance)
{
    if(direct_solution){
        T g=1/(d*lambda+2*mu);
        T b=1/(2*mu);
        T k=strain_trial.Sum();
        TV sh=strain_trial-k/d*TV::All_Ones_Vector();
        T q=sh.Magnitude();
        T dg=(g*q+b*(rho*k+g*sigma_Y))/(g+b*d*rho*rho);
        T m=(k-d*rho*dg)/(d*g);
        T n=(q-dg)/(b*q);
        tau_final=m*TV::All_Ones_Vector()+n*sh;

        T e=n*q;
        if(e<0)
            tau_final=TV::Constant_Vector(-0.5*sigma_Y/rho);

        //T e1=e+rho*m*d+sigma_Y;
        //TV e2=D*tau_final-strain_trial+dg*(n*sh/e+rho*TV::All_Ones_Vector());
        T e3=(tau_final-tau_final.Sum()/d).Magnitude()+rho*tau_final.Sum()+sigma_Y;
        if(abs(e3)>1e-8)
            LOG::printf("BAD PARTICLE: yield function value:%g\n",e3);
        //LOG::printf("errors: %P %P\n",e1,e2);
        //LOG::printf("e: %P\n",e);
        //LOG::printf("crit: %P\n",s.Magnitude()+rho*tau.Sum()+sigma_Y);
        return e>=0;
    }else{
        TVP1 x,residual;
        x.Set_Subvector(0,tau_trial);
        residual=Get_Residual(x,strain_trial,D);
        int k=0;
        while(k++<max_iterations && residual.Magnitude_Squared()>tolerance*tolerance){
            x-=Get_Jacobian(x).Inverse_Times(residual);
            residual=Get_Residual(x,strain_trial,D);}
        x.Get_Subvector(0,tau_final);
        //LOG::printf("residual(%d):%g\n",k,residual.Magnitude());
        //LOG::printf("tau_final=%g\n",tau_final);
        //LOG::printf("F(tau_final)=%g\n",Yield_Function(tau_final));
        //LOG::printf("D*tau_final=%g\n",D*tau_final);
        return residual.Magnitude_Squared()<=tolerance*tolerance;
    }
}
//#####################################################################
// Function Get_Jacobian
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m+1> MPM_DRUCKER_PRAGER<TV>::
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
// s:s = (1-(T)1.0/d)*I1*I1-2*I2
template<class TV> TV MPM_DRUCKER_PRAGER<TV>::
Yield_Function_Derivative(const TV& tau) const
{
    return (tau-tau.Sum()/d).Normalized()+rho; 
}
//#####################################################################
// Function Yield_Function_Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_DRUCKER_PRAGER<TV>::
Yield_Function_Hessian(TV tau) const {
    TV s=tau-tau.Sum()/d;
    T s_norm=s.Normalize();
    return ((T)1-SYMMETRIC_MATRIX<T,d>::Unit_Matrix((T)1/d)-SYMMETRIC_MATRIX<T,d>::Outer_Product(s))/s_norm;
}
//#####################################################################
// Function Get_Residual
//#####################################################################
template<class TV> typename MPM_DRUCKER_PRAGER<TV>::TVP1 MPM_DRUCKER_PRAGER<TV>::
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
template class MPM_DRUCKER_PRAGER<VECTOR<float,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<float,3>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,3>>;
}
