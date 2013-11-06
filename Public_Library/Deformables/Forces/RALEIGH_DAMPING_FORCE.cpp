//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    F+=D_V0;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Lagged_Update_Position_Based_State(const T time)
{
    force.Update_Position_Based_State(time,false);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    force.Add_Implicit_Velocity_Independent_Forces(V,F,coefficient*dt_dv_over_dx/dt*scale,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    force.Enforce_Definiteness(enforce_definiteness_input);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RALEIGH_DAMPING_FORCE<TV>::
Potential_Energy(const T time) const
{
    return -D_V0.Dot(particles.V)*dt/(2*dt_dv_over_dx);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    D_V0.Resize(particles.V.m);
    D_V0.Fill(TV());
    force.Add_Implicit_Velocity_Independent_Forces(particles.V,D_V0,coefficient,time);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    force.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force.Update_Mpi(particle_is_simulated,mpi_solids);
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR RALEIGH_DAMPING_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RALEIGH_DAMPING_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
}
//#####################################################################
namespace PhysBAM{
template class RALEIGH_DAMPING_FORCE<VECTOR<float,2> >;
template class RALEIGH_DAMPING_FORCE<VECTOR<float,3> >;
template class RALEIGH_DAMPING_FORCE<VECTOR<double,2> >;
template class RALEIGH_DAMPING_FORCE<VECTOR<double,3> >;
}
