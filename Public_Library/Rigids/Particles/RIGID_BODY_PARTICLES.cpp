//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
RIGID_BODY_PARTICLES()
    :frame(0,0),twist(0,0),structure_ids(0,0),angular_momentum(0,0),mass(0,0),inertia_tensor(0,0),kinematic(0,0)
{
    Add_Array(ATTRIBUTE_ID_FRAME,&frame);
    Add_Array(ATTRIBUTE_ID_TWIST,&twist);
    Add_Array(ATTRIBUTE_ID_STRUCTURE_IDS,&structure_ids);
    Add_Array(ATTRIBUTE_ID_ANGULAR_MOMENTUM,&angular_momentum);
    Add_Array(ATTRIBUTE_ID_RIGID_MASS,&mass);
    Add_Array(ATTRIBUTE_ID_RIGID_INERTIA_TENSOR,&inertia_tensor);
    Add_Array(ATTRIBUTE_ID_KINEMATIC,&kinematic);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
~RIGID_BODY_PARTICLES()
{
    Delete_All_Particles();
}
static int Initialize_Rigids_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_MASS,"rigid_mass");
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_INERTIA_TENSOR,"rigid_inertia_tensor");
    Register_Attribute_Name(ATTRIBUTE_ID_ANGULAR_MOMENTUM,"angular_momentum");
    Register_Attribute_Name(ATTRIBUTE_ID_KINEMATIC,"kinematic");
    return 0;
}
//#####################################################################
// Resize
//#####################################################################
template<class TV> void RIGID_BODY_PARTICLES<TV>::
Resize(const int new_size)
{
    for(int p=new_size;p<rigid_body.Size();p++)
        if(rigid_body(p)) Remove_Body(p);
    rigid_body.Resize(new_size);
    BASE::Resize(new_size);
}
//#####################################################################
// Function Remove_Body
//#####################################################################
template<class TV> void RIGID_BODY_PARTICLES<TV>::
Remove_Body(const int p)
{
    delete rigid_body(p);
    rigid_body(p)=0;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void RIGID_BODY_PARTICLES<TV>::
Clean_Memory()
{
    for(int p=0;p<this->Size();p++) if(rigid_body(p)) Remove_Body(p);
    rigid_body.Clean_Memory();
    BASE::Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void RIGID_BODY_PARTICLES<TV>::
Delete_All_Particles()
{
    for(int p=0;p<this->Size();p++) if(rigid_body(p)) Remove_Body(p);
    rigid_body.Remove_All();
    BASE::Delete_All_Elements();
}
//#####################################################################
// Function Clone_Helper
//#####################################################################
template<class TV> void RIGID_BODY_PARTICLES<TV>::
Clone_Helper(const RIGID_BODY_PARTICLES& particles)
{
    PHYSBAM_FATAL_ERROR();
}
int initialize_rigids_particles=Initialize_Rigids_Particles();
//#####################################################################
template class RIGID_BODY_PARTICLES<VECTOR<float,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,3> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,3> >;
}
