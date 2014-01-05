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
    bool use_backtracking;
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
    T max_newton_step_size;
    bool debug;
    bool require_one_iteration;
    bool fixed_tolerance;
    bool finish_before_indefiniteness;
    int iterations_used;

    NEWTONS_METHOD()
        :use_golden_section_search(false),use_wolfe_search(true),use_backtracking(false),use_gradient_descent_failsafe(true),tolerance((T)5e-10),
        countdown_tolerance(tolerance),countdown_iterations(5),progress_tolerance((T)5e-10),max_iterations(100),
        krylov_tolerance((T)1e-10),fail_on_krylov_not_converged(false),max_krylov_iterations(100000),
        max_golden_section_iterations(10*sizeof(T)),angle_tolerance(0),use_cg(true),max_newton_step_size(0),debug(false),
        require_one_iteration(false),fixed_tolerance(false),finish_before_indefiniteness(true),iterations_used(0)
    {}

    ~NEWTONS_METHOD(){}

    // Minimize F(x)
    bool Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x);
    void Make_Downhill_Direction(const KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& dx,const KRYLOV_VECTOR_BASE<T>& grad,T norm_grad);
    T Line_Search(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>& tmp2);
    void Make_Vanilla_Newton()
    {
        use_golden_section_search=false;
        use_wolfe_search=false;
        use_backtracking=false;
        use_gradient_descent_failsafe=false;
        countdown_tolerance=0;
        max_newton_step_size=0;
        fixed_tolerance=true;
        finish_before_indefiniteness=false;
    }
};
}
#endif
