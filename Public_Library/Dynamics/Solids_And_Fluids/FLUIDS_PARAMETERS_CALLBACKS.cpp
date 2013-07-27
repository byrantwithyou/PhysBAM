//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
using namespace PhysBAM;
template<class TV> FLUIDS_PARAMETERS_CALLBACKS<TV>::
~FLUIDS_PARAMETERS_CALLBACKS()
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Initialize_Phi()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Reflection_Conditions(T_FACE_ARRAYS_SCALAR& psi_R,const T time)
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Velocities_Masked(const T time,const T_FACE_ARRAYS_BOOL& invalid_mask)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Object_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Object_Velocities(LAPLACE_UNIFORM<TV>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Analytic_Velocities(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Move_Grid_Explicitly(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Adjust_Soot_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Adjust_Density_And_Temperature_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Flame_Speed_Multiplier(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Surface_Tension(ARRAY<T,TV_INT>& surface_tension,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Viscosity(ARRAY<T,TV_INT>& viscosity,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Variable_Vorticity_Confinement(ARRAY<T,TV_INT>& variable_vorticity_confinement,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Divergence(ARRAY<T,TV_INT>& divergence,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_External_Velocity(ARRAY<TV,TV_INT>& V_blend,ARRAY<T,TV_INT>& blend,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> TV FLUIDS_PARAMETERS_CALLBACKS<TV>::
Get_Analytic_Velocity(const TV& location,const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return TV();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Update_Refinement(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Topology_Changed()
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Scalar_Advection_Callback(const T dt,const T time)
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_Before_Advection(const T dt,const T time)
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_Before_Reincorporation(const T dt,const T time)
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Modify_Removed_Particles_After_Reincorporation(const T dt,const T time)
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Initialize_Fluids_Grids()
{
}
template<class TV> void FLUIDS_PARAMETERS_CALLBACKS<TV>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
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
