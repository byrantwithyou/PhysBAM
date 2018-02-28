//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/MINRES.h>
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
    KRYLOV_SOLVER<T>::Ensure_Size(av,x,3+use_gradient_magnitude_objective);
    KRYLOV_VECTOR_BASE<T>& grad=*av.Pop_Value();
    KRYLOV_VECTOR_BASE<T>& dx=*av.Pop_Value();
    KRYLOV_VECTOR_BASE<T>& tm=*av.Pop_Value();
    KRYLOV_VECTOR_BASE<T>& eff_grad=use_gradient_magnitude_objective?*av.Pop_Value():grad;
    GMRES<T> gmres;
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=finish_before_indefiniteness;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    if(use_gmres) krylov=&gmres;
    krylov->relative_tolerance=true;

    bool result=false;
    if(debug) Dump_Parameters();

    T last_E=FLT_MAX;
    int local_max_iterations=max_iterations;
    F.Make_Feasible(x);
    for(iterations_used=0;iterations_used<local_max_iterations;iterations_used++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad)),norm_eff_grad=norm_grad;
        if(use_gradient_magnitude_objective){
            E=norm_grad*norm_grad/2;
            sys.Multiply(grad,eff_grad,true);
            norm_eff_grad=sqrt(sys.Inner_Product(eff_grad,eff_grad));}
        T conv_norm=sys.Convergence_Norm(eff_grad);
        if(use_gradient_magnitude_objective)
            conv_norm=std::min(conv_norm,sys.Convergence_Norm(grad));

        if(debug){
            LOG::printf("GRAD STATS %.16g %.16g %.16g %.16g\n",E,(E-last_E),conv_norm,tolerance);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("newton %d   %.16g %.16g %.16g %.16g",1,iterations_used,E,(E-last_E),conv_norm,tolerance);}

        if(conv_norm<tolerance && (iterations_used || !require_one_iteration || !conv_norm)){result=true;break;}
        if(conv_norm<countdown_tolerance){
            result=true;
            local_max_iterations=std::min(local_max_iterations,iterations_used+countdown_iterations);}
        dx*=0;
        tm.Copy(-1,grad);
        T local_krylov_tolerance=std::min((T).5,krylov_tolerance*(T)sqrt(std::max(conv_norm,tolerance)));
        if(fixed_tolerance) local_krylov_tolerance=krylov_tolerance;
        if(!krylov->Solve(sys,dx,tm,av,local_krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
            break;

        if(use_gradient_descent_failsafe) Make_Downhill_Direction(sys,dx,eff_grad,norm_eff_grad,tm);
        if(max_newton_step_size){
            T norm=sqrt(sys.Inner_Product(dx,dx));
            if(norm>max_newton_step_size) dx*=max_newton_step_size/norm;}

        T a=Line_Search(F,sys,x,dx,grad,tm,eff_grad);
        if(a<=0) break;

        x.Copy(a,dx,x);
        F.Make_Feasible(x);
        last_E=E;}

    if(use_gradient_magnitude_objective) av.Append(&eff_grad);
    av.Append(&tm);
    av.Append(&dx);
    av.Append(&grad);
    return result;
}
//#####################################################################
// Function Make_Downhill_Direction
//#####################################################################
template <class T> void NEWTONS_METHOD<T>::
Make_Downhill_Direction(const KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& dx,const KRYLOV_VECTOR_BASE<T>& grad,T norm_grad,KRYLOV_VECTOR_BASE<T>& tmp)
{
    T norm_dx=sqrt(sys.Inner_Product(dx,dx)),inner=sys.Inner_Product(dx,grad);
    if(inner>angle_tolerance*norm_dx*norm_grad){
        dx*=-1;
        if(debug) LOG::puts("LOOK BACKWARDS");}
    else if(inner>=-angle_tolerance*norm_dx*norm_grad){
        if(debug) LOG::puts("GRADIENT");
        dx.Copy(-norm_dx/norm_grad,grad);}
}
namespace{
template<class T>
class GRADIENT_SQUARED_OBJECTIVE:public NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>
{
public:
    const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F;
    KRYLOV_SYSTEM_BASE<T>& system;
    KRYLOV_VECTOR_BASE<T>& tmp;

    GRADIENT_SQUARED_OBJECTIVE(
        const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,
        KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& tmp)
        :F(F),system(system),tmp(tmp)
    {}

    ~GRADIENT_SQUARED_OBJECTIVE() override {}

    void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,
        KRYLOV_VECTOR_BASE<T>* g,T* e) const override
    {
        PHYSBAM_ASSERT(!h);
        F.Compute(x,g?&system:0,&tmp,0);
        if(e) *e=system.Inner_Product(tmp,tmp)/2;
        if(g) system.Multiply(tmp,*g,true);
    }
};
}
//#####################################################################
// Function Line_Search
//#####################################################################
template <class T> T NEWTONS_METHOD<T>::
Line_Search(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& dx,KRYLOV_VECTOR_BASE<T>& tmp,KRYLOV_VECTOR_BASE<T>& tmp2,KRYLOV_VECTOR_BASE<T>& tmp3)
{
    T a=1;
    GRADIENT_SQUARED_OBJECTIVE<T> GSO(F,sys,tmp3);
    const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& G=use_gradient_magnitude_objective?GSO:F;

    if(use_wolfe_search || use_backtracking || use_golden_section_search){
        PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)> pl(G,x,dx,tmp,&tmp2,&sys);
        if(use_wolfe_search){
            if(!LINE_SEARCH<T>::Line_Search_Wolfe_Conditions(pl,0,1,a,(T)1e-4,(T).9)) a=0;}
        else if(use_backtracking){
            if(!LINE_SEARCH<T>::Line_Search_Backtracking(pl,0,1,a,(T)1e-4)) a=0;}
        else if(use_golden_section_search){
            T phi=(1+sqrt(5))/2;
            a*=phi;
            while(pl(phi*a)<=pl(a)) a*=phi;
            T tau=(T).5*(sqrt((T)5)-1);
            T A=pl(0),B=pl(tau),C=pl(1-tau),D=pl(a),mx=max(A,B,C,D),mn=min(A,B,C,D);
            if(!LINE_SEARCH<T>::Line_Search_Golden_Section(pl,0,a,a,max_golden_section_iterations,(T).25*tolerance)) return 0;
            if(mx-mn<=1e-10*maxabs(mx,mn)){
                if(debug) LOG::printf("OVERRIDE\n");
                a=1;}}}

    if(debug) LOG::printf("ALPHA  %g\n",a);
    return a;
}
//#####################################################################
// Function Make_Vanilla_Newton
//#####################################################################
template <class T> void NEWTONS_METHOD<T>::
Make_Vanilla_Newton()
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
//#####################################################################
// Function Dump_Parameters
//#####################################################################
template <class T> void NEWTONS_METHOD<T>::
Dump_Parameters() const
{
    auto db=[](const char* s,bool f){if(f) LOG::cout<<" "<<s;};
    auto di=[](const char* s,int i){LOG::printf(" %s=%i",s,i);};
    auto df=[](const char* s,T t){LOG::printf(" %s=%.3g",s,t);};
    LOG::cout<<"flags:";
    db("bt",use_backtracking);
    db("cg",use_cg);
    db("db",debug);
    db("fi",finish_before_indefiniteness);
    db("fk",fail_on_krylov_not_converged);
    db("ft",fixed_tolerance);
    db("gd",use_gradient_descent_failsafe);
    db("gm",use_gmres);
    db("go",use_gradient_magnitude_objective);
    db("gs",use_golden_section_search);
    db("ri",require_one_iteration);
    db("ws",use_wolfe_search);
    di("ci",countdown_iterations);
    di("gi",max_golden_section_iterations);
    di("iu",iterations_used);
    di("ki",max_krylov_iterations);
    di("mi",max_iterations);
    df("at",angle_tolerance);
    df("ct",countdown_tolerance);
    df("kt",krylov_tolerance);
    df("ns",max_newton_step_size);
    df("pt",progress_tolerance);
    df("tl",tolerance);
    LOG::cout<<std::endl;
}
namespace PhysBAM{
template struct NEWTONS_METHOD<float>;
template struct NEWTONS_METHOD<double>;
}
