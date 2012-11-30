//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template <class T> bool NEWTONS_METHOD<T>::
Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x)
{
    KRYLOV_VECTOR_BASE<T>& grad=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& neg_dx=*x.Clone_Default();
    MINRES<T> minres;
    KRYLOV_SOLVER<T>& krylov=minres;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    bool result=false;

    T last_E=FLT_MAX;
    for(int i=0;i<max_iterations;i++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad));
        if(abs(E-last_E)<progress_tolerance || norm_grad<tolerance){result=true;break;}

        if(!krylov.Solve(sys,neg_dx,grad,av,krylov_tolerance,0,max_krylov_iterations))
            break;

        if(use_gradient_descent_failsafe){
            T norm_dx=sqrt(sys.Inner_Product(neg_dx,neg_dx)),inner=sys.Inner_Product(neg_dx,grad);
            if(inner<=angle_tolerance*norm_dx*norm_grad){
                neg_dx.Copy(norm_dx/norm_grad,grad);}}


        T a=-1;
        if(use_golden_section_search){
            PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(F,x,neg_dx,grad);
            if(!LINE_SEARCH<T>::Line_Search_Golden_Section(pl,-1,0,a,max_golden_section_iterations,(T).25*tolerance))
                break;}

        x.Copy(a,neg_dx,x);
        last_E=E;}


    av.Delete_Pointers_And_Clean_Memory();
    delete &grad;
    delete &neg_dx;

    return result;
}
template struct NEWTONS_METHOD<float>;
template struct NEWTONS_METHOD<double>;
