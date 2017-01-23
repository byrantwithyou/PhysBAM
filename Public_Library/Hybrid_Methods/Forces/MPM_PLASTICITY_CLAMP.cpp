//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_CLAMP.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PLASTICITY_CLAMP<TV>::
MPM_PLASTICITY_CLAMP(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter,
    T theta_c,T theta_s,T max_hardening,T hardening_factor)
    :MPM_PLASTICITY_MODEL<TV>(particles,gather_scatter),theta_c(theta_c),theta_s(theta_s),max_hardening(max_hardening),
    hardening_factor(hardening_factor)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PLASTICITY_CLAMP<TV>::
~MPM_PLASTICITY_CLAMP()
{
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_PLASTICITY_CLAMP<TV>::
Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,
    typename TV::SPIN* r_diff,const TV& Fe,bool store_hardening,int p) const
{
    strain=clamp(Fe,1-theta_c,1+theta_s);
    if(strain==Fe) return false;
    if(store_hardening){
        T det=Fe.Product()/strain.Product()*particles.Fp(p).Determinant();
        T hardening_coeff=exp(min(max_hardening,hardening_factor*(1-det)));
        particles.mu(p)=particles.mu0(p)*hardening_coeff;
        particles.lambda(p)=particles.lambda0(p)*hardening_coeff;}
    // NOTE: updating particles.mu and particles.lambda is not enough; the constitutive models do not use these.

    if(!dstrain) return true;
    VECTOR<bool,TV::m> cb,ct;
    DIAGONAL_MATRIX<T,TV::m> D;
    for(int i=0;i<TV::m;i++){
        bool b=Fe(i)<1-theta_c;
        bool t=Fe(i)>1+theta_s;
        cb(i)=b;
        ct(i)=t;
        D(i)=!b&&!t;}
    *dstrain=D;

    if(!r_diff || !r_sum) return true;
    for(int i=0;i<TV::SPIN::m;i++){
        int j=(i+1)%TV::m,k=(i+2)%TV::m;
        bool b1=cb(j),b2=cb(k);
        bool t1=ct(j),t2=ct(k);
        bool c1=b1||t1,c2=b2||t2;
        T x=Fe(j),y=Fe(k),w=strain(j),z=strain(k);
        if(!c1 && !c2) (*r_diff)(i)=1;
        else if((b1 && b2) || (t1 && t2)) (*r_diff)(i)=0;
        else (*r_diff)(i)=(w-z)/(x-y);
        (*r_sum)(i)=(w+z)/(x+y);}
    return true;
}
template class MPM_PLASTICITY_CLAMP<VECTOR<float,2>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<float,3>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<double,2>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<double,3>>;
}
