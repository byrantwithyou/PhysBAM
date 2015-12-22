//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_CLAMP.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PLASTICITY_CLAMP<TV>::
MPM_PLASTICITY_CLAMP(MPM_PARTICLES<TV>& particles,T theta_c,T theta_s,T max_hardening,T hardening_factor)
    :MPM_PLASTICITY_MODEL<TV>(particles),theta_c(theta_c),theta_s(theta_s),max_hardening(max_hardening),
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
    PHYSBAM_ASSERT(!dstrain);
    strain=clamp(Fe,1-theta_c,1+theta_s);
    if(strain==Fe) return false;
    if(store_hardening){
        T det=Fe.Product()/strain.Product()*particles.Fp(p).Determinant();
        T hardening_coeff=exp(min(max_hardening,hardening_factor*(1-det)));
        particles.mu(p)=particles.mu0(p)*hardening_coeff;
        particles.lambda(p)=particles.lambda0(p)*hardening_coeff;}
    return true;
}
template class MPM_PLASTICITY_CLAMP<VECTOR<float,2>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<float,3>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<double,2>>;
template class MPM_PLASTICITY_CLAMP<VECTOR<double,3>>;
}
