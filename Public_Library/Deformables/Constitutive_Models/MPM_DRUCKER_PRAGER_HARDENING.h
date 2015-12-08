//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRUCKER_PRAGER_HARDENING__
#define __MPM_DRUCKER_PRAGER_HARDENING__

#include <cmath>

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class TV>
class MPM_DRUCKER_PRAGER_HARDENING:public MPM_PLASTICITY_MODEL<TV>
{
public:
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,d+1> TVP1;

    T mu,lambda;
    T rho_F,plastic_def;
    T a0,a1,a3,a4;
    TV strain_trial,tau_trial,tau_final;
    T delta_gamma_final;
    SYMMETRIC_MATRIX<T,d> D;
    bool direct_solution;


    MPM_DRUCKER_PRAGER_HARDENING(T a0, T a1, T a3, T a4, bool direct_solution=true):
        plastic_def(-1),a0(a0),a1(a1),a3(a3),a4(a4),direct_solution(direct_solution){}
    virtual ~MPM_DRUCKER_PRAGER_HARDENING() {}

    void Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe) override;
    T Yield_Function() const override
    {return Yield_Function(tau_trial);}
    T Yield_Function_Final() const override
    {return Yield_Function(tau_final);}

    T Yield_Function(const TV& tau) const
    {T trace=tau.Sum();return (tau-trace/d).Magnitude()+rho_F*trace;}

    bool Project_Stress(int max_iterations,T tolerance) override;
    TV Get_Updated_Sigma() const override
    {return exp(D*tau_final);}

    void Set_Plastic_Deformation_Lambda(T plastic_def) override;
    T Get_Updated_Plastic_Deformation_Lambda() const override
    {return plastic_def+(D*tau_final-strain_trial).Magnitude();}

    MATRIX<typename TV::SCALAR,TV::m+1> Get_Jacobian(const TVP1& x) const;
    TV Yield_Function_Derivative(const TV& tau) const;
    SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Plastic_Potential_Hessian(TV tau) const;
    TVP1 Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const;
};
}
#endif
