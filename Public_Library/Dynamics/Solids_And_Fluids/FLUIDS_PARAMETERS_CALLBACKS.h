//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_CALLBACKS
//##################################################################### 
#ifndef __FLUIDS_PARAMETERS_CALLBACKS__
#define __FLUIDS_PARAMETERS_CALLBACKS__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Fluids/PhysBAM_Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class TV> class PARTICLE_LEVELSET_PARTICLES;
template<class TV> class PARTICLE_LEVELSET;
template<class T_GRID> class LAPLACE_UNIFORM;
template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class FLUIDS_PARAMETERS_CALLBACKS:public BOUNDARY_CONDITIONS_CALLBACKS<typename T_GRID::VECTOR_T>
{    
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<bool,FACE_INDEX<TV::m> > T_FACE_ARRAYS_BOOL;
public:

    FLUIDS_PARAMETERS_CALLBACKS()
    {}

    virtual ~FLUIDS_PARAMETERS_CALLBACKS()
    {}

//#####################################################################
    virtual void Initialize_Phi(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Reflection_Conditions(T_FACE_ARRAYS_SCALAR& psi_R,const T time){}
    virtual void Get_Source_Velocities_Masked(const T time,const T_FACE_ARRAYS_BOOL& invalid_mask){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();} // allocate mask and set to true where local reseeding should occur
    virtual void Get_Object_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Analytic_Velocities(const T time) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Set_Dirichlet_Boundary_Conditions(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Move_Grid_Explicitly(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Adjust_Soot_With_Sources(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Adjust_Density_And_Temperature_With_Sources(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Flame_Speed_Multiplier(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Variable_Surface_Tension(T_ARRAYS_SCALAR& surface_tension,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Variable_Viscosity(T_ARRAYS_SCALAR& viscosity,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Variable_Vorticity_Confinement(T_ARRAYS_SCALAR& variable_vorticity_confinement,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Divergence(T_ARRAYS_SCALAR& divergence,const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_External_Velocity(ARRAY<TV,TV_INT>& V_blend,T_ARRAYS_SCALAR& blend,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual TV Get_Analytic_Velocity(const TV& location,const T time) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return TV();}
    virtual void Update_Refinement(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Topology_Changed(){}
    virtual void Scalar_Advection_Callback(const T dt,const T time){}
    virtual void Modify_Removed_Particles_Before_Advection(const T dt,const T time){}
    virtual void Modify_Removed_Particles_Before_Reincorporation(const T dt,const T time){}
    virtual void Modify_Removed_Particles_After_Reincorporation(const T dt,const T time){}
    virtual void Initialize_Fluids_Grids(){}
    virtual void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve){}
//#####################################################################
};
}
#endif
