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
#include <boost/function.hpp>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template <class T> bool NEWTONS_METHOD<T>::
Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x)
{
    KRYLOV_VECTOR_BASE<T>& grad=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& dx=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tm=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tmp0=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tmp1=*x.Clone_Default();
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    krylov->relative_tolerance=true;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    T phi=(1+sqrt(5))/2;

    bool result=false;

    T last_E=FLT_MAX;
    int local_max_iterations=max_iterations;
    for(int i=0;i<local_max_iterations;i++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad));
        LOG::printf("GRAD STATS %g %g %g\n", E, abs(E-last_E), norm_grad);

        char buff[1000];
        sprintf(buff,"newton %d   %g %g %g", i, E, abs(E-last_E), norm_grad);
        PHYSBAM_DEBUG_WRITE_SUBSTEP(buff,1,1);

        if(norm_grad<tolerance && (i || !norm_grad)){result=true;break;}
        if(norm_grad<countdown_tolerance){
            result=true;
            local_max_iterations=std::min(local_max_iterations,i+countdown_iterations);}
        dx*=0;
        tm.Copy(-1,grad);
        T local_krylov_tolerance=std::min((T).5,krylov_tolerance*(T)sqrt(norm_grad));
        if(!krylov->Solve(sys,dx,tm,av,local_krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
            break;

        sys.Multiply(dx,tm);
        tm.Copy(1,tm,grad);

        static int solve_id=0;
//        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",sys,tmp0,tmp1);
        solve_id++;

        if(use_gradient_descent_failsafe){
            T norm_dx=sqrt(sys.Inner_Product(dx,dx)),inner=sys.Inner_Product(dx,grad);
            if(inner>angle_tolerance*norm_dx*norm_grad){
                dx*=-1;
                LOG::puts("LOOK BACKWARDS");}
            else if(inner>=-angle_tolerance*norm_dx*norm_grad){
                LOG::puts("GRADIENT");
                dx.Copy(-norm_dx/norm_grad,grad);}}

        T n_grad=sqrt(sys.Inner_Product(grad,grad));
        T n_dx=sqrt(sys.Inner_Product(dx,dx));
        T dx_dot_grad=sys.Inner_Product(dx,grad);
        sys.Multiply(dx,tm);
        T dx_H_dx=sys.Inner_Product(dx,tm);
        sys.Multiply(grad,tm);
        T grad_H_grad=sys.Inner_Product(grad,tm);

        T a=1;
        if(use_wolfe_search){
            PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,grad,&tm);
            if(!LINE_SEARCH<T>::Line_Search_Wolfe_Conditions(pl,0,1,a,(T)1e-4,(T).9)){
                result=false;
                break;}}
        else if(use_golden_section_search){
            PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,grad);
            a*=phi;
            while(pl(phi*a)<=pl(a)) a*=phi;
            T tau=(T).5*(sqrt((T)5)-1);
            T A=pl(0),B=pl(tau),C=pl(1-tau),D=pl(a),mx=max(A,B,C,D),mn=min(A,B,C,D);
            if(!LINE_SEARCH<T>::Line_Search_Golden_Section(pl,0,a,a,max_golden_section_iterations,(T).25*tolerance))
                break;
            if(mx-mn<=1e-10*maxabs(mx,mn)){
                LOG::printf("OVERRIDE\n");
                a=1;}}

        LOG::printf("RAW %g %g    %g    %g %g  %g\n", n_grad,n_dx,    dx_dot_grad/n_dx/n_grad,     grad_H_grad/n_grad/n_grad,dx_H_dx/n_dx/n_dx,  a);

        if(a<=0) break;
        x.Copy(a,dx,x);
        F.Make_Feasible(x);
        last_E=E;}

    av.Delete_Pointers_And_Clean_Memory();
    delete &grad;
    delete &dx;
    delete &tm;
    delete &tmp0;
    delete &tmp1;

    return result;
}
namespace PhysBAM{
template struct NEWTONS_METHOD<float>;
template struct NEWTONS_METHOD<double>;
}
