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
// Function Newtons_Method
//#####################################################################
template <class T> bool NEWTONS_METHOD<T>::
Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x)
{
    KRYLOV_VECTOR_BASE<T>& grad=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& dx=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tm=*x.Clone_Default();
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    krylov->relative_tolerance=true;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    bool result=false;

    T last_E=FLT_MAX;
    int local_max_iterations=max_iterations;
    F.Make_Feasible(x);
    for(int i=0;i<local_max_iterations;i++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad));
        if(debug) LOG::printf("GRAD STATS %g %g %g\n", E, (E-last_E), norm_grad);

        if(debug){
            char buff[1000];
            sprintf(buff,"newton %d   %g %g %g", i, E, (E-last_E), norm_grad);
            PHYSBAM_DEBUG_WRITE_SUBSTEP(buff,1,1);}

        if(norm_grad<tolerance && (i || !require_one_iteration || !norm_grad)){result=true;break;}
        if(norm_grad<countdown_tolerance){
            result=true;
            local_max_iterations=std::min(local_max_iterations,i+countdown_iterations);}
        dx*=0;
        tm.Copy(-1,grad);
        T local_krylov_tolerance=std::min((T).5,krylov_tolerance*(T)sqrt(std::max(norm_grad,tolerance)));
        if(!krylov->Solve(sys,dx,tm,av,local_krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
            break;

        if(use_gradient_descent_failsafe) Make_Downhill_Direction(sys,dx,grad,norm_grad);
        if(max_newton_step_size){
            T norm=sqrt(sys.Inner_Product(dx,dx));
            if(norm>max_newton_step_size) dx*=max_newton_step_size/norm;}

        T a=Line_Search(F,x,dx,grad,tm);
        if(a<=0) break;
        x.Copy(a,dx,x);
        F.Make_Feasible(x);
        last_E=E;}

    av.Delete_Pointers_And_Clean_Memory();
    delete &grad;
    delete &dx;
    delete &tm;

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
Line_Search(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>& tmp2)
{
    T a=1;
    if(use_wolfe_search){
        PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,dx,tmp,&tmp2);
        if(!LINE_SEARCH<T>::Line_Search_Wolfe_Conditions(pl,0,1,a,(T)1e-4,(T).9)) return 0;}
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
