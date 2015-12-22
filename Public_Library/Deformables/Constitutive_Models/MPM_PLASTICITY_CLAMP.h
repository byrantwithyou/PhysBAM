//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_PLASTICITY_CLAMP__
#define __MPM_PLASTICITY_CLAMP__

#include <Tools/Arrays/ATTRIBUTE_ID.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
#include <cmath>
namespace PhysBAM{
const ATTRIBUTE_ID ATTRIBUTE_ID_PLASTIC_DEFORMATION(54);

template<class TV>
class MPM_PLASTICITY_CLAMP:public MPM_PLASTICITY_MODEL<TV>
{
    typedef MPM_PLASTICITY_MODEL<TV> BASE;
public:
    using BASE::particles;

    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;

    T theta_c,theta_s,max_hardening,hardening_factor;

    MPM_PLASTICITY_CLAMP(MPM_PARTICLES<TV>& particles,T theta_c,T theta_s,T max_hardening,T hardening_factor);
    virtual ~MPM_PLASTICITY_CLAMP();

    bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,SYMMETRIC_TENSOR<T,0,TV::m>* ddstrain,
            MATRIX<T,TV::m,TV::SPIN::m>* rdstrain,MATRIX<T,TV::SPIN::m>* rxstrain,
            const TV& Fe,bool store_hardening,int p) const override;
};
}
#endif
