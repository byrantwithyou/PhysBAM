//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRUCKER_PRAGER__
#define __MPM_DRUCKER_PRAGER__

#include <Core/Arrays/ATTRIBUTE_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <cmath>
namespace PhysBAM{
const ATTRIBUTE_ID ATTRIBUTE_ID_DP_RHO_F(55);
const ATTRIBUTE_ID ATTRIBUTE_ID_DP_COHESION(56);

template<class TV>
class MPM_DRUCKER_PRAGER:public MPM_PLASTICITY_MODEL<TV>
{
    typedef MPM_PLASTICITY_MODEL<TV> BASE;
public:
    using BASE::particles;

    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;

    mutable ARRAY_VIEW<T> plastic_def;
    mutable ARRAY_VIEW<T> rho_F;
    mutable ARRAY_VIEW<T> sigma_Y;
    T a0,a1,a3,a4;

    MPM_DRUCKER_PRAGER(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter,T a0,T a1,T a3,T a4);
    virtual ~MPM_DRUCKER_PRAGER();

    void Initialize_Particles(ARRAY<int>* affected_particles) const override;
    void Initialize_Particles(ARRAY<int>* affected_particles,T sigma_Y0) const;
    bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,
        typename TV::SPIN* r_diff,const TV& Fe,bool store_hardening,int p) const override;
    void Update_Hardening(int id,T plastic_def_increment) const;
};
}
#endif
