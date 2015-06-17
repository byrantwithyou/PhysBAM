//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/PARAMETRIC_LINE.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <functional>
using namespace PhysBAM;
//#####################################################################
// Function Newtons_Method
//#####################################################################
template <class T> bool NEWTONS_METHOD<T>::
Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x,ARRAY<KRYLOV_VECTOR_BASE<T>*>& av)
{
    KRYLOV_SOLVER<T>::Ensure_Size(av,x,3);
    KRYLOV_VECTOR_BASE<T>& grad=*av.Pop();
    KRYLOV_VECTOR_BASE<T>& dx=*av.Pop();
    KRYLOV_VECTOR_BASE<T>& tm=*av.Pop();
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=finish_before_indefiniteness;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    krylov->relative_tolerance=true;

    bool result=false;

    T last_E=FLT_MAX;
    int local_max_iterations=max_iterations;
    F.Make_Feasible(x);
    for(iterations_used=0;iterations_used<local_max_iterations;iterations_used++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm2_grad=sqrt(sys.Inner_Product(grad,grad));
        T norm_grad=sys.Convergence_Norm(grad);
        if(debug) LOG::printf("GRAD STATS %.16g %.16g %.16g %.16g\n", E, (E-last_E), norm_grad, tolerance);

        if(debug){
            char buff[1000];
            sprintf(buff,"newton %d   %.16g %.16g %.16g %.16g", iterations_used, E, (E-last_E), norm_grad, tolerance);
            PHYSBAM_DEBUG_WRITE_SUBSTEP(buff,1,1);}

        if(norm_grad<tolerance && (iterations_used || !require_one_iteration || !norm_grad)){result=true;break;}
        if(norm_grad<countdown_tolerance){
            result=true;
            local_max_iterations=std::min(local_max_iterations,iterations_used+countdown_iterations);}
        dx*=0;
        tm.Copy(-1,grad);
        T local_krylov_tolerance=std::min((T).5,krylov_tolerance*(T)sqrt(std::max(norm_grad,tolerance)));
        if(fixed_tolerance) local_krylov_tolerance=krylov_tolerance;
        if(!krylov->Solve(sys,dx,tm,av,local_krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
            break;

        if(use_gradient_descent_failsafe) Make_Downhill_Direction(sys,dx,grad,norm2_grad);
        if(max_newton_step_size){
            T norm=sqrt(sys.Inner_Product(dx,dx));
            if(norm>max_newton_step_size) dx*=max_newton_step_size/norm;}

        T a=Line_Search(F,sys,x,dx,grad,tm);
        if(a<=0) break;
        x.Copy(a,dx,x);
        F.Make_Feasible(x);
        last_E=E;}

    av.Append(&tm);
    av.Append(&dx);
    av.Append(&grad);
    return result;
}
//#####################################################################
// Function Make_Downhill_Direction
//#####################################################################
template <class T> void NEWTONS_METHOD<T>::
Make_Downhill_Direction(const KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& dx,const KRYLOV_VECTOR_BASE<T>& grad,T norm_grad)
{
    T norm_dx=sqrt(sys.Inner_Product(dx,dx)),inner=sys.Inner_Product(dx,grad);
    if(inner>angle_tolerance*norm_dx*norm_grad){
        dx*=-1;
        if(debug) LOG::puts("LOOK BACKWARDS");}
    else if(inner>=-angle_tolerance*norm_dx*norm_grad){
        if(debug) LOG::puts("GRADIENT");
        dx.Copy(-norm_dx/norm_grad,grad);}
}
//#####################################################################
// Function Line_Search
//#####################################################################
template <class T> T NEWTONS_METHOD<T>::
Line_Search(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>& tmp2)
{
    T a=1;
    if(use_wolfe_search){
        PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,tmp,&tmp2,&sys);
        if(!LINE_SEARCH<T>::Line_Search_Wolfe_Conditions(pl,0,1,a,(T)1e-4,(T).9)) return 0;}
    else if(use_backtracking){
        PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,tmp,&tmp2,&sys);
        if(!LINE_SEARCH<T>::Line_Search_Backtracking(pl,0,1,a,(T)1e-4)) return 0;}
    else if(use_golden_section_search){
        PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,tmp);
        T phi=(1+sqrt(5))/2;
        a*=phi;
        while(pl(phi*a)<=pl(a)) a*=phi;
        T tau=(T).5*(sqrt((T)5)-1);
        T A=pl(0),B=pl(tau),C=pl(1-tau),D=pl(a),mx=max(A,B,C,D),mn=min(A,B,C,D);
        if(!LINE_SEARCH<T>::Line_Search_Golden_Section(pl,0,a,a,max_golden_section_iterations,(T).25*tolerance)) return 0;
        if(mx-mn<=1e-10*maxabs(mx,mn)){
            if(debug) LOG::printf("OVERRIDE\n");
            a=1;}}

    if(debug) LOG::printf("ALPHA  %g\n", a);
    return a;
}
namespace PhysBAM{
template struct NEWTONS_METHOD<float>;
template struct NEWTONS_METHOD<double>;
}
