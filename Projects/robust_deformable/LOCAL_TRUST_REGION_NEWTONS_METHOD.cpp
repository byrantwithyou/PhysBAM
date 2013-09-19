//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include "LOCAL_TR_CONJUGATE_GRADIENT.h"
#include "LOCAL_TRUST_REGION_NEWTONS_METHOD.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template <class T> bool LOCAL_TRUST_REGION_NEWTONS_METHOD<T>::
Trust_Region_Newtons_Method(const NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>& F,KRYLOV_SYSTEM_BASE<T>& sys,KRYLOV_VECTOR_BASE<T>& x)
{
    KRYLOV_VECTOR_BASE<T>& grad=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& dx=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& tm=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& pu=*x.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& pk=*x.Clone_Default();
    MINRES<T> minres;
    CONJUGATE_GRADIENT<T> cg;
    LOCAL_TR_CONJUGATE_GRADIENT<T> trcg;
    KRYLOV_SOLVER<T>* krylov=&minres;
    if(use_cg) krylov=&cg;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    bool result=false;

    T last_E=FLT_MAX;
    T trust_region;
    for(int i=0;i<max_iterations;i++){
        T E=0;
        F.Compute(x,&sys,&grad,&E);
        T norm_grad=sqrt(sys.Inner_Product(grad,grad));
        if(i==0){trust_region=norm_grad;}
        printf("%g %g %g (%g %g)\n", E, abs(E-last_E), norm_grad, progress_tolerance, tolerance);
        if((abs(E-last_E)<progress_tolerance&&0) || norm_grad<tolerance){result=true;break;}
        
        if(use_dogleg){
            sys.Multiply(grad,tm);
            T grad_H_grad=sys.Inner_Product(grad,tm);
            pu.Copy((T)-(norm_grad*norm_grad)/grad_H_grad,grad);
            T norm_pu=sqrt(sys.Inner_Product(pu,pu));
            if(norm_pu>=trust_region&&i>0&&trust_region>(T)0.05){pk.Copy(trust_region/norm_pu,pu);printf("unconstrained\n");}
            else{
                dx*=0;
                if(!krylov->Solve(sys,dx,grad,av,krylov_tolerance,0,max_krylov_iterations) && fail_on_krylov_not_converged)
                    break;

                dx*=-1;
                T dx_dx=sys.Inner_Product(dx,dx);
                T n_dx=sqrt(dx_dx);
                if(n_dx==0){pk=pu;}
                if(use_gradient_descent_failsafe){
                    T inner=sys.Inner_Product(dx,grad);
                    printf("angle %g\n", inner/(n_dx*norm_grad));
                    if(-inner<=-angle_tolerance*n_dx*norm_grad){
                        dx*=-1;
                        puts("LOOK BACKWARDS");}}
                if(i==0){trust_region=n_dx;}
                if(n_dx<=trust_region||trust_region<(T)0.05){pk=dx;printf("newton\n");}
                else{
                    dx-=pu;
                    T pu_dot_pbu=sys.Inner_Product(pu,dx);
                    T pbu_pbu=sys.Inner_Product(dx,dx);
                    T tau=(pbu_pbu-pu_dot_pbu+sqrt(pbu_pbu*trust_region*trust_region))/pbu_pbu;
                    if(tau<=(T)1&&tau>=(T)2) tau=(pbu_pbu-pu_dot_pbu-sqrt(pbu_pbu*trust_region*trust_region))/pbu_pbu;
                    if(tau>(T)2){tau=(T)2;}
                    printf("interpolation\n");
                    //printf("tau: %g npu: %g ndx: %g npbu: %g\n",tau,norm_pu,n_dx,sqrt(pbu_pbu));
                    PHYSBAM_ASSERT(tau>=(T)1&&tau<=(T)2);
                    pk.Copy(tau-(T)1,dx,pu);}}
        }
        else{
            pk*=0;
            if(!trcg.Solve(sys,pk,grad,av,krylov_tolerance,0,max_krylov_iterations,trust_region) && fail_on_krylov_not_converged)
                break;
            pk*=-1;
        }

//it's the same algorithm as [Lin and More, 1999].
        T testE;
        T n_pk=sqrt(sys.Inner_Product(pk,pk));
        tm.Copy((T)1,pk,x);
        F.Compute(tm,0,0,&testE);
        sys.Multiply(pk,tm);
        T pk_H_pk=sys.Inner_Product(pk,tm);
        T pk_dot_grad=sys.Inner_Product(pk,grad);

        T actred=E-testE;
        T prered=-(pk_dot_grad+0.5*pk_H_pk);
        T alpha;
        if(actred+pk_dot_grad>=0){alpha=sigma3;}
        else{alpha=max(sigma1,(T)0.5*pk_dot_grad/(actred+pk_dot_grad));}
        T previous_tr=trust_region;//for debug
        if(actred<=eta0*prered){trust_region=min(max(alpha,sigma1)*n_pk,sigma2*trust_region);printf("shrink1: ");} 
        else if(actred<=eta1*prered){trust_region=max(sigma1*trust_region,min(alpha*n_pk,sigma2*trust_region));printf("shrink2: ");}
        else if(actred<=eta2*prered){trust_region=max(sigma1*trust_region,min(alpha*n_pk,sigma3*trust_region));printf("expand or shrink: ");}
        else{trust_region=max(trust_region,min(alpha*n_pk,sigma3*trust_region));printf("expand: ");}
        printf("%g -> %g\n",previous_tr,trust_region);
        printf("rho: %g  n_pk: %g  trustRegion: %g  E-testE: %g\n",actred/prered,n_pk,trust_region,actred);
        //PHYSBAM_ASSERT(trust_region>(T)1e-10);
        if(actred>eta0*prered){x.Copy((T)1,pk,x);E=testE;printf("success\n");}
        else if(trust_region<(T)0.001){x.Copy((T)-0.5,pk,x);F.Compute(x,0,0,&E);printf("backword\n");trust_region=(T)0.01;}//if the trust region is improper small, it will expand again.
        else{printf("trust region faild\n");}

        last_E=E;}

    av.Delete_Pointers_And_Clean_Memory();
    delete &grad;
    delete &dx;
    delete &tm;
    delete &pu;
    delete &pk;

    return result;
}
namespace PhysBAM{
template struct LOCAL_TRUST_REGION_NEWTONS_METHOD<float>;
template struct LOCAL_TRUST_REGION_NEWTONS_METHOD<double>;
}
