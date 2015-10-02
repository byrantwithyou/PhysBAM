//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_GRID_FORCES<TV>::
PARTICLE_GRID_FORCES(const MPM_FORCE_HELPER<TV>& force_helper)
    :force_helper(force_helper),particles(force_helper.particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_GRID_FORCES<TV>::
~PARTICLE_GRID_FORCES()
{
}
template class PARTICLE_GRID_FORCES<VECTOR<float,1> >;
template class PARTICLE_GRID_FORCES<VECTOR<float,2> >;
template class PARTICLE_GRID_FORCES<VECTOR<float,3> >;
template class PARTICLE_GRID_FORCES<VECTOR<double,1> >;
template class PARTICLE_GRID_FORCES<VECTOR<double,2> >;
template class PARTICLE_GRID_FORCES<VECTOR<double,3> >;
}
