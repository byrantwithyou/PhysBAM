//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTICLE_GRID_WEIGHTS__
#define __PARTICLE_GRID_WEIGHTS__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int d> class SYMMETRIC_MATRIX;

template<class TV>
class PARTICLE_GRID_WEIGHTS
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    struct SCRATCH
    {
        ARRAY<TV_INT> index;
        ARRAY<T> weight;
        ARRAY<TV> gradient;
    };

    bool use_gradient_transfer;
    bool constant_scalar_inertia_tensor;

    PARTICLE_GRID_WEIGHTS(int threads);
    virtual ~PARTICLE_GRID_WEIGHTS();

    virtual void Compute(int p,SCRATCH& scratch,bool want_gradient) const=0;
    virtual void Update(const ARRAY_VIEW<TV>& X)=0;
    virtual T Constant_Scalar_Inverse_Dp() const=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Dp(const TV& X) const=0;
    virtual int Order() const=0;
//#####################################################################
};
}
#endif
