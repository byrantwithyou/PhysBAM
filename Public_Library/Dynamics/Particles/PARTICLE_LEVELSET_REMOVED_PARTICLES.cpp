//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Vectors/VECTOR.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>::
PARTICLE_LEVELSET_REMOVED_PARTICLES()
    :next(0)
{
    this->Store_Velocity();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>::
~PARTICLE_LEVELSET_REMOVED_PARTICLES()
{}
//#####################################################################
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >;
}
