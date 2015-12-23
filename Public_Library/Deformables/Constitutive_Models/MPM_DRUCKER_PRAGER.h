//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRUCKER_PRAGER__
#define __MPM_DRUCKER_PRAGER__

#include <Tools/Arrays/ATTRIBUTE_ID.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
#include <cmath>
namespace PhysBAM{
const ATTRIBUTE_ID ATTRIBUTE_ID_PLASTIC_DEFORMATION(54);
const ATTRIBUTE_ID ATTRIBUTE_ID_DP_RHO_F(55);

template<class TV>
class MPM_DRUCKER_PRAGER:public MPM_PLASTICITY_MODEL<TV>
{
    typedef MPM_PLASTICITY_MODEL<TV> BASE;
public:
    using BASE::particles;using BASE::gather_scatter;

    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;

    mutable ARRAY_VIEW<T> plastic_def;
    mutable ARRAY_VIEW<T> rho_F;
    T a0,a1,a3,a4;

    MPM_DRUCKER_PRAGER(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter,T a0,T a1,T a3,T a4);
    virtual ~MPM_DRUCKER_PRAGER();

    void Initialize_Particles() const override;
    bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,
        typename TV::SPIN* r_diff,const TV& Fe,bool store_hardening,int p) const override;
    void Update_Hardening(int id,T plastic_def_increment) const;
};
}
#endif
