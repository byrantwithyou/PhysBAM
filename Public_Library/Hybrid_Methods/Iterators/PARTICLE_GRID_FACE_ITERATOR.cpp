//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_FACE_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_GRID_FACE_ITERATOR<TV>::
PARTICLE_GRID_FACE_ITERATOR(const VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights,
    int p,bool want_gradient,SCRATCH& scratch)
    :scratch(scratch),axis(0),i(0)
{
    for(int i=0;i<TV::m;i++)
        weights(i)->Compute(p,scratch(i),want_gradient);
}
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<float,1> >;
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<float,2> >;
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<float,3> >;
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<double,1> >;
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<double,2> >;
template class PARTICLE_GRID_FACE_ITERATOR<VECTOR<double,3> >;
}
