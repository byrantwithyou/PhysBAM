//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_PLASTICITY_MODEL__
#define __MPM_PLASTICITY_MODEL__

#include <Core/Arrays/ATTRIBUTE_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;
const ATTRIBUTE_ID ATTRIBUTE_ID_PLASTIC_DEFORMATION(54);

template<class TV>
class MPM_PLASTICITY_MODEL
{
public:
    typedef typename TV::SCALAR T;
    MPM_PARTICLES<TV>& particles;
    GATHER_SCATTER<TV>* gather_scatter;
    bool use_implicit;

    MPM_PLASTICITY_MODEL(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter=0);
    MPM_PLASTICITY_MODEL(const MPM_PLASTICITY_MODEL&) = delete;
    void operator=(const MPM_PLASTICITY_MODEL&) = delete;
    virtual ~MPM_PLASTICITY_MODEL();

    virtual void Initialize_Particles(const ARRAY<int>* affected_particles) const;
    virtual bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,
        typename TV::SPIN* r_diff,const TV& Fe,bool store_hardening,int p) const=0;
    virtual void Update_Particles() const;
};
}
#endif
