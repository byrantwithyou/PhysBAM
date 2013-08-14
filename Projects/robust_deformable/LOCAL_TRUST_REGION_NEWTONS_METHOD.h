//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_TRUST_REGION_NEWTONS_METHOD
//#####################################################################
#ifndef __LOCAL_TRUST_REGION_NEWTONS_METHOD__
#define __LOCAL_TRUST_REGION_NEWTONS_METHOD__

#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>

namespace PhysBAM{
template <class T>
struct LOCAL_TRUST_REGION_NEWTONS_METHOD
{
    T eta0;
    T eta1;
    T eta2;
    T sigma1;
    T sigma2;
    T sigma3;
    T use_gradient_descent_failsafe;
    T tr_revise_torelance;
    T tolerance;
    T progress_tolerance;
    int max_iterations;
    T krylov_tolerance;
    bool fail_on_krylov_not_converged;
    int max_krylov_iterations;
    T angle_tolerance;
    bool use_cg;
    bool use_dogleg;
    bool use_trcg;
    LOCAL_TRUST_REGION_NEWTONS_METHOD()
        :eta0((T)1e-4),eta1((T)0.25),eta2((T)0.75),sigma1((T)0.25),sigma2((T)0.5),sigma3((T)4),
        use_gradient_descent_failsafe(true),tr_revise_torelance((T)0.20),tolerance((T)5e-10),
        progress_tolerance((T)5e-10),max_iterations(100),krylov_tolerance((T)1e-10),fail_on_krylov_not_converged(false),
        max_krylov_iterations(100000),angle_tolerance(0),use_cg(false),use_dogleg(false),use_trcg(true)
    {}

    ~LOCAL_TRUST_REGION_NEWTONS_METHOD(){}

    // Minimize F(x)
    bool Trust_Region_Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x);
};
}
#endif
