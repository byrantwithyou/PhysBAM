//#####################################################################
// Copyright 2015, Peter Chen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include "SMOKE_PARTICLES.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_PARTICLES<TV>::
SMOKE_PARTICLES()
{
    this->Store_Velocity();
    this->Store_Mass();
    Add_Array(ATTRIBUTE_ID_SMOKE_C,&C);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_PARTICLES<TV>::
~SMOKE_PARTICLES()
{}
//#####################################################################
// Function Initialize_Smoke_Particles
//#####################################################################
static int Initialize_Smoke_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_SMOKE_C,"C");
    return 0;
}
int initialize_deformables_particles=Initialize_Smoke_Particles();
//#####################################################################
template class SMOKE_PARTICLES<VECTOR<float,2> >;
template class SMOKE_PARTICLES<VECTOR<float,3> >;
template class SMOKE_PARTICLES<VECTOR<double,2> >;
template class SMOKE_PARTICLES<VECTOR<double,3> >;
}
