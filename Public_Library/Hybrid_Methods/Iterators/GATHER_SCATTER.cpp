//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GATHER_SCATTER<TV>::
GATHER_SCATTER(const ARRAY<int>& simulated_particles,const PARTICLE_GRID_WEIGHTS<TV>* weights)
    :simulated_particles(simulated_particles),weights(weights),threads(weights->thread_scratch.m)
{
}
template class GATHER_SCATTER<VECTOR<float,2> >;
template class GATHER_SCATTER<VECTOR<float,3> >;
template class GATHER_SCATTER<VECTOR<double,2> >;
template class GATHER_SCATTER<VECTOR<double,3> >;
}
