//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWTONS_METHOD
//#####################################################################
#ifndef __NEWTONS_METHOD__
#define __NEWTONS_METHOD__

#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>

namespace PhysBAM{
template <class T>
struct NEWTONS_METHOD
{
    bool use_golden_section_search;
    bool use_wolfe_search;
    bool use_gradient_descent_failsafe;
    T tolerance;
    T countdown_tolerance;
    int countdown_iterations;
    T progress_tolerance;
    int max_iterations;
    T krylov_tolerance;
    bool fail_on_krylov_not_converged;
    int max_krylov_iterations;
    int max_golden_section_iterations;
    T angle_tolerance;
    bool use_cg;

    NEWTONS_METHOD()
        :use_golden_section_search(false),use_wolfe_search(true),use_gradient_descent_failsafe(true),tolerance((T)5e-10),
        countdown_tolerance(tolerance),countdown_iterations(5),progress_tolerance((T)5e-10),max_iterations(100),
        krylov_tolerance((T)1e-10),fail_on_krylov_not_converged(false),max_krylov_iterations(100000),
        max_golden_section_iterations(10*sizeof(T)),angle_tolerance(0),use_cg(true)
    {}

    ~NEWTONS_METHOD(){}

    // Minimize F(x)
    bool Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x);
};
}
#endif
