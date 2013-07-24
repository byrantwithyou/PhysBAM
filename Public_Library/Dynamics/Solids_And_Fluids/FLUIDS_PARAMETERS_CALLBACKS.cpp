//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
using namespace PhysBAM;
template<class T_GRID> FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
~FLUIDS_PARAMETERS_CALLBACKS()
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Initialize_Phi()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Reflection_Conditions(T_FACE_ARRAYS_SCALAR& psi_R,const T time)
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Source_Velocities_Masked(const T time,const T_FACE_ARRAYS_BOOL& invalid_mask)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Object_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Analytic_Velocities(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Move_Grid_Explicitly(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Adjust_Soot_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Adjust_Density_And_Temperature_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Flame_Speed_Multiplier(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Variable_Surface_Tension(T_ARRAYS_SCALAR& surface_tension,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Variable_Viscosity(T_ARRAYS_SCALAR& viscosity,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Variable_Vorticity_Confinement(T_ARRAYS_SCALAR& variable_vorticity_confinement,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Divergence(T_ARRAYS_SCALAR& divergence,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_External_Velocity(ARRAY<TV,TV_INT>& V_blend,T_ARRAYS_SCALAR& blend,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> typename T_GRID::VECTOR_T FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Get_Analytic_Velocity(const TV& location,const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return TV();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Update_Refinement(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Topology_Changed()
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Scalar_Advection_Callback(const T dt,const T time)
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Modify_Removed_Particles_Before_Advection(const T dt,const T time)
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Modify_Removed_Particles_Before_Reincorporation(const T dt,const T time)
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Modify_Removed_Particles_After_Reincorporation(const T dt,const T time)
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Initialize_Fluids_Grids()
{
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class T_GRID> void FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::
Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve)
{
}
namespace PhysBAM{
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<double,1> > >;
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<double,2> > >;
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<double,3> > >;
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<float,1> > >;
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<float,2> > >;
template class FLUIDS_PARAMETERS_CALLBACKS<GRID<VECTOR<float,3> > >;
}
