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
    bool use_golden_section_search=false;
    bool use_wolfe_search=true;
    bool use_backtracking=false;
    bool use_gradient_descent_failsafe=true;
    T tolerance=(T)5e-10;
    T countdown_tolerance=tolerance;
    int countdown_iterations=5;
    T progress_tolerance=tolerance;
    int max_iterations=100;
    T krylov_tolerance=(T)1e-10;
    bool fail_on_krylov_not_converged=false;
    int max_krylov_iterations=100000;
    int max_golden_section_iterations=10*sizeof(T);
    T angle_tolerance=(T)1e-2;
    bool use_cg=true;
    bool use_gmres=false;
    T max_newton_step_size=100;
    bool debug=false;
    bool require_one_iteration=false;
    bool fixed_tolerance=false;
    bool finish_before_indefiniteness=true;
    int iterations_used=0;
    bool use_gradient_magnitude_objective=false;

    NEWTONS_METHOD()=default;
    ~NEWTONS_METHOD()=default;

    // Minimize F(x)
    bool Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x,ARRAY<KRYLOV_VECTOR_BASE<T>*>& av);
    void Make_Downhill_Direction(const KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& dx,const KRYLOV_VECTOR_BASE<T>& grad,T norm_grad,KRYLOV_VECTOR_BASE<T>& tmp);
    T Line_Search(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>& tmp2,KRYLOV_VECTOR_BASE<T>& tmp3);
    void Make_Vanilla_Newton();
    void Dump_Parameters() const;
};
}
#endif
