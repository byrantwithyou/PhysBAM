//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/Robust_Functions.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_DRUCKER_PRAGER<TV>::
MPM_DRUCKER_PRAGER(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter,T a0,T a1,T a3,T a4)
    :MPM_PLASTICITY_MODEL<TV>(particles,gather_scatter),a0(a0),a1(a1),a3(a3),a4(a4)
{
    particles.Add_Array("plastic_deformation",&plastic_def);
    particles.Add_Array("dp_rho_f",&rho_F);
    particles.Add_Array("dp_cohesion",&sigma_Y);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_DRUCKER_PRAGER<TV>::
~MPM_DRUCKER_PRAGER()
{
}
//#####################################################################
// Function Initialize_Particle
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Initialize_Particles(ARRAY<int>* affected_particles) const
{
    Initialize_Particles(affected_particles,0);
}
//#####################################################################
// Function Initialize_Particle
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Initialize_Particles(ARRAY<int>* affected_particles,T sigma_Y0) const
{
    if(affected_particles)
#pragma omp parallel for
        for(int k=0;k<affected_particles->m;k++){
            int p=(*affected_particles)(k);
            sigma_Y(p)=sigma_Y0;
            Update_Hardening(p,0);}
    else
        if(particles.X.m==0)
            PHYSBAM_WARNING("Adding Drucker-Prager before partices have been added to the system: NO PARTICLES WILL BE INITIALIZED!");
        else
#pragma omp parallel for
            for(int p=0;p<particles.X.m;p++){
                sigma_Y(p)=sigma_Y0;
                Update_Hardening(p,0);}
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_DRUCKER_PRAGER<TV>::
Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,typename TV::SPIN* r_diff,const TV& Fe,
    bool store_hardening,int id) const
{
    typedef typename TV::SPIN TVSPIN;
    T mu=particles.mu(id),lambda=particles.lambda(id);
    T g=1/(d*lambda+2*mu);
    T b=1/(2*mu);
    T beta=-sigma_Y(id)*g/(rho_F(id)*TV::m);
    TV strain_trial=log(abs(Fe))+beta;
    T k=strain_trial.Sum();
    TV sh=strain_trial-k/d;
    T q=sh.Magnitude();
    T r=rho_F(id)*k/g;
    T dg=q+b*r;
    if(dg<=0) return false;
    if(q==0 || k>0){
        strain=sigma_Y(id)?exp(-beta)*TV::All_Ones_Vector():TV::All_Ones_Vector();
        if(dstrain) *dstrain=MATRIX<T,d>();
        if(store_hardening) Update_Hardening(id,(strain_trial-beta).Magnitude());
        return true;} 
    T p=r*b/q;
    TV h_strain=k/d-p*sh;
    strain=exp(h_strain-beta);
    if(store_hardening) Update_Hardening(id,dg);

// What if Fe<0?
    if(!dstrain) return true;
    TV o=TV::All_Ones_Vector();
    DIAGONAL_MATRIX<T,d> DMS(strain);
    TV one_over_Fe=(T)1/Fe;

    DIAGONAL_MATRIX<T,d> D_strain_trial(one_over_Fe);
    TV D_k=one_over_Fe;
    MATRIX<T,d> H1=MATRIX<T,d>::Outer_Product(o,D_k/d);
    TV D_q=D_strain_trial*sh/q;
    TV D_r=rho_F(id)/g*D_k;
    TV D_p=(b*D_r-p*D_q)/q;
    MATRIX<T,d> D_h_strain=(1+p)*H1-MATRIX<T,d>::Outer_Product(sh,D_p)-p*D_strain_trial;
    *dstrain=DMS*D_h_strain;
    if(!r_diff || !r_sum) return true;
    TVSPIN RX_exp,RX_ln;
    for(int i=0;i<TVSPIN::m;i++){
        T x=Fe((i+1)%d),y=Fe((i+2)%d);
        T e1=h_strain((i+1)%d),e2=h_strain((i+2)%d);
        RX_exp(i)=diff_exp_over_diff(e1,e2);
        RX_ln(i)=diff_log_over_diff(x,y);
        (*r_sum)(i)=(strain((i+1)%d)+strain((i+2)%d))/(x+y);}
    *r_diff=-p*RX_exp*RX_ln;
    return true;
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Update_Hardening(int id,T plastic_def_increment) const
{
    plastic_def(id)+=plastic_def_increment;
    T phi_F=pi/180*(a0+(a1*plastic_def(id)-a4)*exp(-a3*plastic_def(id)));
    T sin_phi_F=sin(phi_F);
    rho_F(id)=sqrt(2.0/3.0)*2*sin_phi_F/(3-sin_phi_F);
}
template class MPM_DRUCKER_PRAGER<VECTOR<float,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<float,3>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,3>>;
}
