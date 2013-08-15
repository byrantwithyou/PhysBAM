//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_TR_CONJUGATE_GRADIENT
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Log/LOG.h>
#include <cfloat>
#include <limits>
#include "LOCAL_TR_CONJUGATE_GRADIENT.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LOCAL_KRYLOV_SOLVER<T>::
LOCAL_KRYLOV_SOLVER()
    :print_diagnostics(true),print_residuals(false),relative_tolerance(false),nullspace_tolerance((T)1e-5),iterations_used(0),residual_magnitude_squared(0),
    nullspace_measure(0),restart_iterations(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> LOCAL_KRYLOV_SOLVER<T>::
~LOCAL_KRYLOV_SOLVER()
{
}
//#####################################################################
// Function Ensure_Size
//#####################################################################
template<class T> void LOCAL_KRYLOV_SOLVER<T>::
Ensure_Size(ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,const KRYLOV_VECTOR_BASE<T>& v,int size)
{
    for(int i=0,m=min(size,av.m);i<m;i++)
        av(i)->Resize(v);

    if(size<=av.m) return;
    int old=av.m;
    av.Resize(size);
    for(int i=old;i<av.m;i++)
        av(i)=v.Clone_Default();
}
namespace PhysBAM{
template class LOCAL_KRYLOV_SOLVER<float>;
template class LOCAL_KRYLOV_SOLVER<double>;
}


//#####################################################################
// Destructor
//#####################################################################
template<class T> LOCAL_TR_CONJUGATE_GRADIENT<T>::
~LOCAL_TR_CONJUGATE_GRADIENT()
{}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool LOCAL_TR_CONJUGATE_GRADIENT<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,
    ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,T tolerance,const int min_iterations,const int max_iterations,T trust_region)
{
    Ensure_Size(av,x,3+system.use_preconditioner);
    KRYLOV_VECTOR_BASE<T>& q=*av(0);
    KRYLOV_VECTOR_BASE<T>& s=*av(1);
    KRYLOV_VECTOR_BASE<T>& r=*av(2);
    KRYLOV_VECTOR_BASE<T>& z=*av(2+system.use_preconditioner);
    KRYLOV_VECTOR_BASE<T>& tm=*x.Clone_Default();

    // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);
    T rho_old=(T)FLT_MAX;T convergence_norm=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
            if(print_residuals) LOG::cout<<"restarting cg"<<std::endl;
            r=b;system.Multiply(x,q);r-=q;system.Project(r);}
        // stopping conditions
        system.Project_Nullspace(r);
        convergence_norm=system.Convergence_Norm(r);
        if(relative_tolerance && iterations==0) tolerance*=convergence_norm;
        if(print_residuals) LOG::cout<<convergence_norm<<std::endl;
        if(convergence_norm<=tolerance && (iterations>=min_iterations || convergence_norm<small_number)){
            if(print_diagnostics) LOG::Stat("cg iterations",iterations);
            if(iterations_used) *iterations_used=iterations;
            return true;}
        if(iterations==max_iterations) break;
        // actual iteration
        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(r,z);
        T rho=(T)system.Inner_Product(mr,r);
        if(restart) s=mr;
        else s.Copy(rho/rho_old,s,mr);
        system.Multiply(s,q);
        system.Project(q);
        T s_dot_q=(T)system.Inner_Product(s,q);
        if(s_dot_q<=0){
            if(finish_before_indefiniteness){
                if(iterations==0) x=b;
                if(print_diagnostics) LOG::Stat("cg iterations",iterations);
                if(iterations_used) *iterations_used=iterations;
                return true;}
            LOG::cout<<"CG: matrix appears indefinite or singular, s_dot_q/s_dot_s="
                     <<s_dot_q/(T)system.Inner_Product(s,s)<<std::endl;}
        T alpha=s_dot_q?rho/s_dot_q:(T)FLT_MAX;
        tm=x;
        x.Copy(alpha,s,x);
        //trust_region
        T norm_x=sqrt(system.Inner_Product(x,x));
        if(norm_x>=trust_region){
            //T tm_dot_tm=system.Inner_Product(tm,tm);
             T tm_dot_s=system.Inner_Product(tm,s);
             T s_dot_s=system.Inner_Product(s,s);
             T tau=(-tm_dot_s+trust_region*sqrt(s_dot_s))/s_dot_s;
             if(tau<0){tau=(-tm_dot_s-trust_region*sqrt(s_dot_s))/s_dot_s;}
             printf("trust region: %lf\n",tau);
             x.Copy(tau,s,tm);
             return true;}

        r.Copy(-alpha,q,r);
        rho_old=rho;}

    delete &tm;
    if(print_diagnostics){
        LOG::Stat("cg iterations",iterations);
        LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;}
    if(iterations_used) *iterations_used=iterations;
    return false;
}
//#####################################################################
namespace PhysBAM{
template class LOCAL_TR_CONJUGATE_GRADIENT<float>;
template class LOCAL_TR_CONJUGATE_GRADIENT<double>;
}
