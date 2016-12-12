//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTICLE_GRID_WEIGHTS_SPLINE__
#define __PARTICLE_GRID_WEIGHTS_SPLINE__
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV,int degree>
class PARTICLE_GRID_WEIGHTS_SPLINE:public PARTICLE_GRID_WEIGHTS<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PARTICLE_GRID_WEIGHTS<TV> BASE;
public:
    enum {n=degree+1};
    struct PRECOMPUTE_DATA
    {
        TV_INT base;
        VECTOR<TV,n> w,dw;
    };

    const GRID<TV>& grid;
    ARRAY<PRECOMPUTE_DATA> precompute_data;

    PARTICLE_GRID_WEIGHTS_SPLINE(const GRID<TV>& grid,int threads=1);
    virtual ~PARTICLE_GRID_WEIGHTS_SPLINE();

    void Compute(int p,typename BASE::SCRATCH& scratch,bool want_gradient) const;
    void Update(const ARRAY_VIEW<TV>& X);
    T Constant_Scalar_Inverse_Dp() const;
    SYMMETRIC_MATRIX<T,TV::m> Dp(const TV& X) const;
    int Order() const;
//#####################################################################
};
}
#endif
