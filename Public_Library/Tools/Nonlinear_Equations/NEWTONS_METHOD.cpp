//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/PARAMETRIC_LINE.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
using namespace PhysBAM;
void (*NM_Flush_State)(const char*)=0;
//#####################################################################
// Constructor
//#####################################################################
template <class T> bool NEWTONS_METHOD<T>::
Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x)
{
    KRYLOV_VECTOR_BASE<T>& grad=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& neg_dx=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tm=*x.Clone_Default();
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    T phi=(1+sqrt(5))/2;

    bool result=false;

    char buff[1000];

    T last_E=FLT_MAX;
    for(int i=0;i<max_iterations;i++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad));
//        PHYSBAM_ASSERT(abs(E-last_E)>=progress_tolerance);
        printf("%g %g %g (%g %g)\n", E, abs(E-last_E), norm_grad, progress_tolerance, tolerance);

        sprintf(buff,"newton %d   %g %g %g", i, E, abs(E-last_E), norm_grad);
        if(NM_Flush_State) NM_Flush_State(buff);

        if((abs(E-last_E)<progress_tolerance&&0) || norm_grad<tolerance){result=true;break;}
        neg_dx*=0;
        if(!krylov->Solve(sys,neg_dx,grad,av,krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
            break;

        printf("A   %g\n", sys.Inner_Product(neg_dx,grad));
        if(use_gradient_descent_failsafe){
            T norm_dx=sqrt(sys.Inner_Product(neg_dx,neg_dx)),inner=sys.Inner_Product(neg_dx,grad);
            printf("angle %g\n", inner/(norm_dx*norm_grad));
            if(inner<=-angle_tolerance*norm_dx*norm_grad){
                neg_dx*=-1;
                puts("LOOK BACKWARDS");}
            else if(inner<=angle_tolerance*norm_dx*norm_grad){
                puts("GRADIENT");
                neg_dx.Copy(clamp(norm_dx,(T)1e-3,(T)1e3)/norm_grad,grad);}}

        T n_grad=sqrt(sys.Inner_Product(grad,grad));
        T n_dx=sqrt(sys.Inner_Product(neg_dx,neg_dx));
        T dx_dot_grad=sys.Inner_Product(neg_dx,grad);
        sys.Multiply(neg_dx,tm);
        T dx_H_dx=sys.Inner_Product(neg_dx,tm);
        sys.Multiply(grad,tm);
        T grad_H_grad=sys.Inner_Product(grad,tm);
        printf("RAW %g %g    %g    %g %g\n", n_grad,n_dx,    dx_dot_grad/n_dx/n_grad,     grad_H_grad/n_grad/n_grad,dx_H_dx/n_dx/n_dx);



        T a=-1;
        if(use_golden_section_search){
            PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,neg_dx,grad);
            a*=phi;
            while(pl(phi*a)<=pl(a)) a*=phi;
            T tau=(T).5*(sqrt((T)5)-1);
            T A=pl(0),B=pl(tau),C=pl(1-tau),D=pl(a),mx=max(A,B,C,D),mn=min(A,B,C,D);
            if(!LINE_SEARCH<T>::Line_Search_Golden_Section(pl,a,0,a,max_golden_section_iterations,(T).25*tolerance))
                break;
            if(mx-mn<=1e-10*maxabs(mx,mn)){
                printf("OVERRIDE\n");
                a=-1;}}

        x.Copy(a,neg_dx,x);
        last_E=E;}


    av.Delete_Pointers_And_Clean_Memory();
    delete &grad;
    delete &neg_dx;
    delete &tm;

    return result;
}
namespace PhysBAM{
template struct NEWTONS_METHOD<float>;
template struct NEWTONS_METHOD<double>;
}
