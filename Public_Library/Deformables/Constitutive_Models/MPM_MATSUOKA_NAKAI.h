//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MATSUOKA_NAKAI__
#define __MPM_MATSUOKA_NAKAI__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class TV>
class MPM_MATSUOKA_NAKAI:public MPM_PLASTICITY_MODEL<TV>
{
public:
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,d+1> TVP1;

    T kappa,c;
    TV strain_trial,tau_trial,tau_final;
    SYMMETRIC_MATRIX<T,d> D;

    MPM_MATSUOKA_NAKAI(T fa,T c);
    virtual ~MPM_MATSUOKA_NAKAI(){}

    virtual void Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,d>& Fe);
    virtual T Yield_Function() const;
    T Yield_Function(const TV& tau) const;
    virtual T Yield_Function_Final() const
    {return Yield_Function(tau_final);}

    virtual bool Project_Stress(int max_iterations,T tolerance);
    virtual TV Get_Updated_Sigma() const
    {return exp(D*tau_final);}

    MATRIX<typename TV::SCALAR,TV::m+1> Get_Jacobian(const TVP1& x) const;
    SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Yield_Function_Hessian(TV tau) const;
    TV Yield_Function_Derivative(const TV& tau) const;
    TVP1 Get_Residual(const TVP1& x,const TV& strain_trial,const SYMMETRIC_MATRIX<T,d>& D) const;
};
}
#endif
