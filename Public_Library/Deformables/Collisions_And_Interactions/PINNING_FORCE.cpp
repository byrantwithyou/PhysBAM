//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PINNING_FORCE<TV>::
PINNING_FORCE(DEFORMABLE_PARTICLES<TV>& particles,T& dt,T stiffness,T damping_coefficient)
    :LAGGED_FORCE<TV>(particles),stiffness(stiffness),damping_coefficient(damping_coefficient),dt(dt),E(0),H_E(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PINNING_FORCE<TV>::
~PINNING_FORCE()
{
}
//#####################################################################
// Function Lagged_Update_Position_Based_State
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Lagged_Update_Position_Based_State(const T time)
{
    X0=particles.X.Subset(pinned_particles);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
/*
  x = x0 + dt * va;
  xt = xt0 + dt * vta;

  va = (x-x0)/dt;
  vta = (xt-xt0)/dt;

  f = -stiffness*(x-xt)-damping_coefficient*(va-vta);
  f = -stiffness*(x-xt)-damping_coefficient*((x-x0)/dt-(xt-xt0)/dt);
  f = -(stiffness+damping_coefficient/dt)*x + (stiffness*xt-damping_coefficient*(-x0-xt+xt0)/dt);
  f = -tau*x + (stiffness*xt+damping_coefficient*(x0+xt-xt0)/dt);
  f = -tau*(x - sigma);
  E = tau/2*(x - sigma)^2;
 */
    E=0;
    T tau=stiffness+damping_coefficient/dt;
    H_E=tau;
    grad_E.Resize(pinned_particles.m);
    for(int i=0;i<pinned_particles.m;i++){
        TV x=particles.X(pinned_particles(i)),x0=X0(i);
        TV xt=targets(i)(time),xt0=targets(i)(time-dt),sigma;
        if(tau) sigma=(stiffness*xt+damping_coefficient*(x0+xt-xt0)/dt)/tau;
        TV dx=x-sigma;
        E+=tau/2*dx.Magnitude_Squared();
        grad_E(i)=tau*dx;}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    F.Subset(pinned_particles)-=grad_E;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    F.Subset(pinned_particles)-=V.Subset(pinned_particles)*H_E;
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR PINNING_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void PINNING_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR PINNING_FORCE<TV>::
Potential_Energy(const T time) const
{
    return E;
}
//#####################################################################
template class PINNING_FORCE<VECTOR<float,2> >;
template class PINNING_FORCE<VECTOR<float,3> >;
template class PINNING_FORCE<VECTOR<double,2> >;
template class PINNING_FORCE<VECTOR<double,3> >;
}
