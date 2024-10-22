//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_GRID_WEIGHTS<TV>::
PARTICLE_GRID_WEIGHTS(int threads)
    :use_gradient_transfer(false),stencil_width(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_GRID_WEIGHTS<TV>::
~PARTICLE_GRID_WEIGHTS()
{
}
template class PARTICLE_GRID_WEIGHTS<VECTOR<float,1> >;
template class PARTICLE_GRID_WEIGHTS<VECTOR<float,2> >;
template class PARTICLE_GRID_WEIGHTS<VECTOR<float,3> >;
template class PARTICLE_GRID_WEIGHTS<VECTOR<double,1> >;
template class PARTICLE_GRID_WEIGHTS<VECTOR<double,2> >;
template class PARTICLE_GRID_WEIGHTS<VECTOR<double,3> >;
}
