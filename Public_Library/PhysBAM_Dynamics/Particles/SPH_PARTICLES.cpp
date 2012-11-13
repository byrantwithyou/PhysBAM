//#####################################################################
// Copyright 2004-2008, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/SPH_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SPH_PARTICLES<TV>::
SPH_PARTICLES()
{
    this->Store_Velocity();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SPH_PARTICLES<TV>::
~SPH_PARTICLES()
{}
//#####################################################################
template class SPH_PARTICLES<VECTOR<float,1> >;
template class SPH_PARTICLES<VECTOR<float,2> >;
template class SPH_PARTICLES<VECTOR<float,3> >;
template class SPH_PARTICLES<VECTOR<double,1> >;
template class SPH_PARTICLES<VECTOR<double,2> >;
template class SPH_PARTICLES<VECTOR<double,3> >;
}
