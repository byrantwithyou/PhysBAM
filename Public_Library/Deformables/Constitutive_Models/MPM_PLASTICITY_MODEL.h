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

template<class TV>
class MPM_PLASTICITY_MODEL:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;
    MPM_PARTICLES<TV>& particles;
    bool use_implicit;

    MPM_PLASTICITY_MODEL(MPM_PARTICLES<TV>& particles);
    virtual ~MPM_PLASTICITY_MODEL();

    virtual void Initialize_Particle(int p) const {};
    virtual bool Compute(TV& strain,MATRIX<T,TV::m>* dstrain,SYMMETRIC_TENSOR<T,0,TV::m>* ddstrain,
            MATRIX<T,TV::m,TV::SPIN::m>* rdstrain,MATRIX<T,TV::SPIN::m>* rxstrain,
            const TV& Fe,bool store_hardening,int p) const=0;
    virtual bool Update_Particle(int p) const;
};
}
#endif
