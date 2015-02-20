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
GATHER_SCATTER(const ARRAY<int>& simulated_particles)
    :simulated_particles(simulated_particles),weights(0),threads(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GATHER_SCATTER<TV>::
~GATHER_SCATTER()
{
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void GATHER_SCATTER<TV>::
Set_Weights(const PARTICLE_GRID_WEIGHTS<TV>* weights_input)
{
    weights=weights_input;
    threads=weights->thread_scratch.m;
}
template class GATHER_SCATTER<VECTOR<float,2> >;
template class GATHER_SCATTER<VECTOR<float,3> >;
template class GATHER_SCATTER<VECTOR<double,2> >;
template class GATHER_SCATTER<VECTOR<double,3> >;
}
