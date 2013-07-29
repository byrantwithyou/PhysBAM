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
#include <Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class TV> class PARTICLE_LEVELSET_PARTICLES;
template<class TV> class PARTICLE_LEVELSET;
template<class TV> class LAPLACE_UNIFORM;
template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV>
class FLUIDS_PARAMETERS_CALLBACKS:public BOUNDARY_CONDITIONS_CALLBACKS<TV>
{    
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:

    FLUIDS_PARAMETERS_CALLBACKS()
    {}

    virtual ~FLUIDS_PARAMETERS_CALLBACKS();

//#####################################################################
    virtual void Initialize_Phi();
    virtual void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time);
    virtual void Get_Reflection_Conditions(ARRAY<T,FACE_INDEX<TV::m> >& psi_R,const T time);
    virtual void Get_Source_Velocities_Masked(const T time,const ARRAY<bool,FACE_INDEX<TV::m> >& invalid_mask);
    virtual void Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time);
    virtual void Get_Object_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    virtual void Get_Object_Velocities(LAPLACE_UNIFORM<TV>* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    virtual void Get_Analytic_Velocities(const T time) const;
    virtual void Set_Dirichlet_Boundary_Conditions(const T time);
    virtual void Move_Grid_Explicitly(const T time);
    virtual void Adjust_Soot_With_Sources(const T time);
    virtual void Adjust_Density_And_Temperature_With_Sources(const T time);
    virtual void Get_Flame_Speed_Multiplier(const T dt,const T time);
    virtual void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time);
    virtual void Get_Variable_Surface_Tension(ARRAY<T,TV_INT>& surface_tension,const T time);
    virtual void Get_Variable_Viscosity(ARRAY<T,TV_INT>& viscosity,const T time);
    virtual void Get_Variable_Vorticity_Confinement(ARRAY<T,TV_INT>& variable_vorticity_confinement,const T time);
    virtual void Get_Divergence(ARRAY<T,TV_INT>& divergence,const T dt,const T time);
    virtual void Get_External_Velocity(ARRAY<TV,TV_INT>& V_blend,ARRAY<T,TV_INT>& blend,const T time);
    virtual TV Get_Analytic_Velocity(const TV& location,const T time) const;
    virtual void Update_Refinement(const T dt,const T time);
    virtual void Topology_Changed();
    virtual void Scalar_Advection_Callback(const T dt,const T time);
    virtual void Modify_Removed_Particles_Before_Advection(const T dt,const T time);
    virtual void Modify_Removed_Particles_Before_Reincorporation(const T dt,const T time);
    virtual void Modify_Removed_Particles_After_Reincorporation(const T dt,const T time);
    virtual void Initialize_Fluids_Grids();
    virtual void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    virtual void Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve);
//#####################################################################
};
}
#endif
