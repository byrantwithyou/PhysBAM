//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
RIGID_BODY_PARTICLES()
    :frame(0,0),twist(0,0),structure_ids(0,0),angular_momentum(0,0),mass(0,0),inertia_tensor(0,0),kinematic(0,0)
{
    Add_Array("frame",&frame);
    Add_Array("twist",&twist);
    Add_Array("structure_ids",&structure_ids);
    Add_Array("angular_momentum",&angular_momentum);
    Add_Array("rigid_mass",&mass);
    Add_Array("rigid_inertia_tensor",&inertia_tensor);
    Add_Array("kinematic",&kinematic);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
~RIGID_BODY_PARTICLES()
{
    Delete_All_Particles();
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
//#####################################################################
template class RIGID_BODY_PARTICLES<VECTOR<float,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,3> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,3> >;
}
