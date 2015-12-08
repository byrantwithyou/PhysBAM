//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MATSUOKA_NAKAI_WITH_DP__
#define __MPM_MATSUOKA_NAKAI_WITH_DP__

#include <cmath>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class TV>
class MPM_MATSUOKA_NAKAI_WITH_DP:public MPM_PLASTICITY_MODEL<TV>
{
public:
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,d+1> TVP1;

    T kappa_F,c,rho_G,plastic_def;
    T phi_F, phi_cs, psi_G;
    T a0,a1,a2,a3,a4,a5;
    TV strain_trial,tau_trial,tau_final;
    T delta_gamma_final;
    SYMMETRIC_MATRIX<T,d> D;

    MPM_MATSUOKA_NAKAI_WITH_DP();
    virtual ~MPM_MATSUOKA_NAKAI_WITH_DP(){}

    void Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe) override;
    T Yield_Function() const override
    {return Yield_Function(tau_trial);}
    T Yield_Function(const TV& tau) const;
    T Yield_Function_Final() const override
    {return Yield_Function(tau_final);}

    bool Project_Stress(int max_iterations,T tolerance) override;
    TV Get_Updated_Sigma() const override
    {return exp(D*tau_final);}

    void Set_Plastic_Deformation_Lambda(T plastic_def) override;
    T Get_Updated_Plastic_Deformation_Lambda() const override
    {return plastic_def+delta_gamma_final;}

    MATRIX<typename TV::SCALAR,TV::m+1> Get_Jacobian(const TVP1& x) const;
    TV Yield_Function_Derivative(const TV& tau) const;
    SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Yield_Function_Hessian(TV tau) const;
    TV Plastic_Flow(const TV& tau) const;
    TVP1 Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const;
    T Get_phi_F(const TV& tau, T lambda) const
    {return M_PI/180*(a0+(a1*lambda-a4)*exp(a2*tau.Sum()-a3*lambda));}
    T Get_phi_cs(T lambda) const
    {return M_PI/180*(a0-a4*exp(-a5*lambda));}
};
}
#endif
