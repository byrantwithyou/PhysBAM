//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRUCKER_PRAGER__
#define __MPM_DRUCKER_PRAGER__

#include <cmath>

#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class TV>
class MPM_DRUCKER_PRAGER:public MPM_PLASTICITY_MODEL<TV>
{
public:
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,d+1> TVP1;

    T rho,sigma_Y;
    TV strain_trial,tau_trial,tau_final;
    SYMMETRIC_MATRIX<T,d> D;

    MPM_DRUCKER_PRAGER(T friction_angle,T cohesion):rho(2*sin(friction_angle)/(sqrt(3)*(3-sin(friction_angle)))),sigma_Y(-2*sqrt(3)*cohesion*cos(friction_angle)/(3-sin(friction_angle))){}
    virtual ~MPM_DRUCKER_PRAGER() {}

    virtual void Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe);
    virtual T Yield_Function() const
    {return Yield_Function(tau_trial);}

    T Yield_Function(const TV& tau) const
    {T trace=tau.Sum();return (tau-trace/d).Magnitude()+rho*trace+sigma_Y;}

    virtual bool Project_Stress(int max_iterations,T tolerance);
    virtual TV Get_Updated_Sigma() const
    {return exp(D*tau_final);}

    MATRIX<typename TV::SCALAR,TV::m+1> Get_Jacobian(const TVP1& x) const;
    TV Yield_Function_Derivative(const TV& tau) const;
    TVP1 Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const;
};
}
#endif
