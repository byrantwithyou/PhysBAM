#if 0
//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/MPM_SQUARED_DRUCKER_PRAGER.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> void MPM_SQUARED_DRUCKER_PRAGER<TV>::
MPM_SQUARED_DRUCKER_PRAGER(T friction_angle,T cohesion)
    :rho(2*sin(friction_angle)/(sqrt(3)*(3-sin(friction_angle)))),sigma_Y(-2*sqrt(3)*cohesion*cos(friction_angle)/(3-sin(friction_angle)))
{
    PHYSBAM_ASSERT(cohesion>=0);
}
//#####################################################################
// Function Set_Lame_Constants_And_F_Elastic
//#####################################################################
template<class TV> void MPM_SQUARED_DRUCKER_PRAGER<TV>::
Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe)
{
    D=SYMMETRIC_MATRIX<T,d>::Unit_Matrix(-lambda/(2*d*mu*lambda+4*mu*mu))+SYMMETRIC_MATRIX<T,d>::Identity_Matrix()/(2*mu);

    strain_trial=log(Fe.x);
    tau_trial=D.Inverse_Times(strain_trial);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_SQUARED_DRUCKER_PRAGER<TV>::
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
    //LOG::printf("residual(%d):%g\n",k,residual.Magnitude());
    //LOG::printf("tau_final=%g\n",tau_final);
    //LOG::printf("F(tau_final)=%g\n",Yield_Function(tau_final));
    //LOG::printf("D*tau_final=%g\n",D*tau_final);
    return residual.Magnitude_Squared()<=tolerance*tolerance;
}
//#####################################################################
// Function Get_Jacobian
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m+1> MPM_SQUARED_DRUCKER_PRAGER<TV>::
Get_Jacobian(const TVP1& x) const
{
    TV tau;
    x.Get_Subvector(0,tau);
    T delta_gamma=x(d);
    TV f_hat_vec=Yield_Function_Derivative(tau);
    MATRIX<T,d+1> ret;
    ret.Set_Submatrix(0,0,D+delta_gamma*Yield_Function_Hessian());
    for(int i=0;i<d;i++){
        ret(i,d)=f_hat_vec(i);
        ret(d,i)=f_hat_vec(i);}

    return ret;
}
//#####################################################################
// Function Yield_Function_Derivative
//#####################################################################
template<class TV> TV MPM_SQUARED_DRUCKER_PRAGER<TV>::
Yield_Function_Derivative(const TV& tau) const
{
    return (T)2*(tau-tau.Sum()/d)-2*(rho*tau.Sum()+sigma_Y)*rho;
}
//#####################################################################
// Function Yield_Function_Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_SQUARED_DRUCKER_PRAGER<TV>::
Yield_Function_Hessian() const {
    return (T)2-SYMMETRIC_MATRIX<T,d>::Unit_Matrix(2*rho*rho+(T)2/d);
}
//#####################################################################
// Function Get_Residual
//#####################################################################
template<class TV> typename MPM_SQUARED_DRUCKER_PRAGER<TV>::TVP1 MPM_SQUARED_DRUCKER_PRAGER<TV>::
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
template class MPM_SQUARED_DRUCKER_PRAGER<VECTOR<float,2>>;
template class MPM_SQUARED_DRUCKER_PRAGER<VECTOR<float,3>>;
template class MPM_SQUARED_DRUCKER_PRAGER<VECTOR<double,2>>;
template class MPM_SQUARED_DRUCKER_PRAGER<VECTOR<double,3>>;
}
#endif
