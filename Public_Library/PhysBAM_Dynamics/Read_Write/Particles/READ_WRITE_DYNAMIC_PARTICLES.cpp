//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Particles/READ_WRITE_DEFORMABLES_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Particles/READ_WRITE_RIGIDS_PARTICLES.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
namespace PhysBAM{

void Initialize_Particles()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_VORTICITY,"vorticity");
    Register_Attribute_Name(ATTRIBUTE_ID_E,"E");
    Register_Attribute_Name(ATTRIBUTE_ID_RHO,"rho");
    Register_Attribute_Name(ATTRIBUTE_ID_PHI,"phi");
    Register_Attribute_Name(ATTRIBUTE_ID_GRAD_PHI,"grad_phi");
    Register_Attribute_Name(ATTRIBUTE_ID_AGE,"age");
    Register_Attribute_Name(ATTRIBUTE_ID_MATERIAL_VOLUME,"material_volume");
    Register_Attribute_Name(ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE,"quantized_collision_distance");

    #define READ_WRITE_HELPER(RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<unsigned short>();

    READ_WRITE_HELPER(float);
    #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    READ_WRITE_HELPER(double);
    #endif

    Initialize_Rigids_Particles();
    Initialize_Deformables_Particles();
    Initialize_Geometry_Particle();

}
}
