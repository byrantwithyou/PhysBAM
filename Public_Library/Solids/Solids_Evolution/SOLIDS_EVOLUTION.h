//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EVOLUTION
//#####################################################################
#ifndef __SOLIDS_EVOLUTION__
#define __SOLIDS_EVOLUTION__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Rigids/Rigids_Evolution/RIGIDS_KINEMATIC_EVOLUTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV,bool world_space> class RIGID_BODY_MASS;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class GRID;
template<class TV> class EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
class SOLIDS_EVOLUTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    RIGID_BODY_COLLISIONS<TV>* rigid_body_collisions;
    RIGID_DEFORMABLE_COLLISIONS<TV>* rigid_deformable_collisions;
    T time;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks;
    bool fully_implicit;
protected:
    ARRAY<KRYLOV_VECTOR_BASE<T>*> krylov_vectors;
    GENERALIZED_VELOCITY<TV>& GV_B;
public:
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass_inverse;
private:
    static SOLIDS_EVOLUTION_CALLBACKS<TV> solids_evolution_callbacks_default;
public:
    RIGIDS_KINEMATIC_EVOLUTION<TV> kinematic_evolution;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities;

    SOLIDS_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities);
    SOLIDS_EVOLUTION(const SOLIDS_EVOLUTION&) = delete;
    void operator=(const SOLIDS_EVOLUTION&) = delete;
    virtual ~SOLIDS_EVOLUTION();

    void Set_Solids_Evolution_Callbacks(SOLIDS_EVOLUTION_CALLBACKS<TV>& solids_evolution_callbacks_input)
    {solids_evolution_callbacks=&solids_evolution_callbacks_input;}

//#####################################################################
    virtual bool Use_CFL() const=0;
    virtual void Euler_Step_Position(const T dt,const T time);
    virtual void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids)=0;
    virtual void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids)=0;
    virtual void Initialize_Deformable_Objects(const T frame_rate,const bool restart);
    virtual void Initialize_Rigid_Bodies(const T frame_rate, const bool restart)=0;
    virtual bool Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(const T dt,const T time,int& repulsions_found,int& collisions_found,const bool exit_early=false);
    virtual void Postprocess_Frame(const int frame);
protected:
    void Save_Position(ARRAY<TV>& X,ARRAY<FRAME<TV> >& rigid_frame);
    void Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const FRAME<TV> > rigid_frame);
    //void Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id); // convenience wrapper
    void Clamp_Velocities();
public:
    void Initialize_World_Space_Masses();
    void Restore_Position_After_Hypothetical_Position_Evolution(ARRAY<TV>& X_save,ARRAY<FRAME<TV> >& rigid_frame_save);
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time);
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time);
    void Euler_Step_Position(const T dt,const T time,const int p);
//#####################################################################
};
}
#endif
