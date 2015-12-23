//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_PLASTICITY_MODEL__
#define __MPM_PLASTICITY_MODEL__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class GATHER_SCATTER;

template<class TV>
class MPM_PLASTICITY_MODEL:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;
    MPM_PARTICLES<TV>& particles;
    GATHER_SCATTER<TV>* gather_scatter;
    bool use_implicit;

    MPM_PLASTICITY_MODEL(MPM_PARTICLES<TV>& particles,GATHER_SCATTER<TV>* gather_scatter=0);
    virtual ~MPM_PLASTICITY_MODEL();

    virtual void Initialize_Particles() const;
    virtual bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,typename TV::SPIN* r_sum,
        typename TV::SPIN* r_diff,const TV& Fe,bool store_hardening,int p) const=0;
    virtual void Update_Particles() const;
};
}
#endif
