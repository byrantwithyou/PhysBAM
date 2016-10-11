//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTICLE_GRID_ITERATOR__
#define __PARTICLE_GRID_ITERATOR__
#include <Core/Vectors/VECTOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_GRID_ITERATOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename PARTICLE_GRID_WEIGHTS<TV>::SCRATCH SCRATCH;
    SCRATCH& scratch;
    int i;

    PARTICLE_GRID_ITERATOR(const PARTICLE_GRID_WEIGHTS<TV>* weights,int p,bool want_gradient,SCRATCH& scratch);
    const TV_INT& Index() const {return scratch.index(i);}
    const T& Weight() const {return scratch.weight(i);}
    const TV& Gradient() const {return scratch.gradient(i);}

    void Next() {i++;}
    bool Valid() const {return i<scratch.index.m;}
};
}
#endif
