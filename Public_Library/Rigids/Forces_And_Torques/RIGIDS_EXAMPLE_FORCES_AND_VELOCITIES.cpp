//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
~RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES()
{}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_External_Positions(ARRAY_VIEW<FRAME<TV> > rotation,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> bool RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return false;
}
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
namespace PhysBAM{
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,1> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,2> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,3> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,1> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,2> >;
template class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,3> >;
}
