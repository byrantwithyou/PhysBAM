//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_GRID_ITERATOR<TV>::
PARTICLE_GRID_ITERATOR(const PARTICLE_GRID_WEIGHTS<TV>* weights,int p,bool want_gradient,SCRATCH& scratch)
    :scratch(scratch),i(0)
{
    weights->Compute(p,scratch,want_gradient);
}
template class PARTICLE_GRID_ITERATOR<VECTOR<float,2> >;
template class PARTICLE_GRID_ITERATOR<VECTOR<float,3> >;
template class PARTICLE_GRID_ITERATOR<VECTOR<double,2> >;
template class PARTICLE_GRID_ITERATOR<VECTOR<double,3> >;
}
