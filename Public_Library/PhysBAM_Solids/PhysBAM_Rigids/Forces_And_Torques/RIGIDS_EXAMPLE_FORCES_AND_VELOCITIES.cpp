//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
~RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES()
{}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Rigid_Particle_Is_Simulated
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,1> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,2> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,1> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,2> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,3> >;
#endif
