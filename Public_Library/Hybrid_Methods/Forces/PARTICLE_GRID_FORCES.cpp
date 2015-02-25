//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_GRID_FORCES<TV>::
PARTICLE_GRID_FORCES(MPM_PARTICLES<TV>& particles)
    :particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_GRID_FORCES<TV>::
~PARTICLE_GRID_FORCES()
{
}
//#####################################################################
// Function Capture_Stress
//#####################################################################
template<class TV> void PARTICLE_GRID_FORCES<TV>::
Capture_Stress()
{
}
template class PARTICLE_GRID_FORCES<VECTOR<float,2> >;
template class PARTICLE_GRID_FORCES<VECTOR<float,3> >;
template class PARTICLE_GRID_FORCES<VECTOR<double,2> >;
template class PARTICLE_GRID_FORCES<VECTOR<double,3> >;
}
