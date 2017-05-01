//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Initialize_Phi()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Reflection_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Reflection_Conditions(ARRAY<T,FACE_INDEX<TV::m> >& psi_R,const T time)
{
}
//#####################################################################
// Function Get_Source_Velocities_Masked
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Velocities_Masked(const T time,const ARRAY<bool,FACE_INDEX<TV::m> >& invalid_mask)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Object_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Object_Velocities(LAPLACE_UNIFORM<TV>* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Analytic_Velocities(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Move_Grid_Explicitly
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Move_Grid_Explicitly(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Adjust_Soot_With_Sources
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Adjust_Soot_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Adjust_Density_And_Temperature_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Flame_Speed_Multiplier
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Flame_Speed_Multiplier(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Variable_Surface_Tension
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Surface_Tension(ARRAY<T,TV_INT>& surface_tension,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Variable_Viscosity
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Viscosity(ARRAY<T,TV_INT>& viscosity,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Variable_Vorticity_Confinement
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Vorticity_Confinement(ARRAY<T,TV_INT>& variable_vorticity_confinement,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Divergence
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Divergence(ARRAY<T,TV_INT>& divergence,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_External_Velocity
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_External_Velocity(ARRAY<TV,TV_INT>& V_blend,ARRAY<T,TV_INT>& blend,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Analytic_Velocity
//#####################################################################
template<class TV> TV FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Analytic_Velocity(const TV& location,const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return TV();
}
//#####################################################################
// Function Update_Refinement
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Update_Refinement(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Topology_Changed
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Topology_Changed()
{
}
//#####################################################################
// Function Scalar_Advection_Callback
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Scalar_Advection_Callback(const T dt,const T time)
{
}
//#####################################################################
// Function Modify_Removed_Particles_Before_Advection
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_Before_Advection(const T dt,const T time)
{
}
//#####################################################################
// Function Modify_Removed_Particles_Before_Reincorporation
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_Before_Reincorporation(const T dt,const T time)
{
}
//#####################################################################
// Function Modify_Removed_Particles_After_Reincorporation
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_After_Reincorporation(const T dt,const T time)
{
}
//#####################################################################
// Function Initialize_Fluids_Grids
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Initialize_Fluids_Grids()
{
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Substitute_Coupling_Matrices
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve)
{
}
namespace PhysBAM{
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<double,1> >;
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<double,2> >;
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<double,3> >;
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<float,1> >;
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<float,2> >;
template class FLUIDS_PARAMETERS_CALLBACKS<VECTOR<float,3> >;
}
