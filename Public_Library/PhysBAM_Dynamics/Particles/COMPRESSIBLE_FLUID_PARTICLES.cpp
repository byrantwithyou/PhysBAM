//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_PARTICLES<TV>::
COMPRESSIBLE_FLUID_PARTICLES()
    :X(0,0),rho(0,0),E(0,0),phi(0,0),grad_phi(0,0),V(0,0)
{
    Add_Array(ATTRIBUTE_ID_X,&X);
    Add_Array(ATTRIBUTE_ID_RHO,&rho);
    Add_Array(ATTRIBUTE_ID_E,&E);
    Add_Array(ATTRIBUTE_ID_PHI,&phi);
    Add_Array(ATTRIBUTE_ID_GRAD_PHI,&grad_phi);
    Add_Array(ATTRIBUTE_ID_V,&V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_PARTICLES<TV>::
~COMPRESSIBLE_FLUID_PARTICLES()
{}
static int Initialize_Compressible_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_VORTICITY,"vorticity");
    Register_Attribute_Name(ATTRIBUTE_ID_E,"E");
    Register_Attribute_Name(ATTRIBUTE_ID_RHO,"rho");
    Register_Attribute_Name(ATTRIBUTE_ID_PHI,"phi");
    Register_Attribute_Name(ATTRIBUTE_ID_GRAD_PHI,"grad_phi");
    Register_Attribute_Name(ATTRIBUTE_ID_AGE,"age");
    Register_Attribute_Name(ATTRIBUTE_ID_MATERIAL_VOLUME,"material_volume");
    Register_Attribute_Name(ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE,"quantized_collision_distance");
    return 0;
}
int initialize_compressible_particles=Initialize_Compressible_Particles();
//#####################################################################
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,1> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,2> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,1> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,2> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,3> >;
#endif
}
