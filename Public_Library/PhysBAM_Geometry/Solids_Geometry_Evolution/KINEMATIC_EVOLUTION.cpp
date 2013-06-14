//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KINEMATIC_EVOLUTION<TV>::
KINEMATIC_EVOLUTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigid_body_example_velocities,bool use_kinematic_keyframes_input)
    :rigid_body_collection(rigid_body_collection_input),rigid_body_example_velocities(rigid_body_example_velocities),use_kinematic_keyframes(use_kinematic_keyframes_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KINEMATIC_EVOLUTION<TV>::
~KINEMATIC_EVOLUTION()
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=rigid_body_collection.kinematic_rigid_bodies(i);Set_External_Velocities(twist(p),velocity_time,p);}
    rigid_body_example_velocities.Set_External_Velocities(twist,velocity_time,current_position_time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
    if(rigid_body.is_static){twist=TWIST<TV>();return;}
    if(!use_kinematic_keyframes){Set_Kinematic_Velocities(twist,(T)1e-3,time,id);return;} // Use 1e-3 for backward differencing if the example did not implement velocity calculations.
    RIGID_BODY_STATE<TV> interpolated_state;
    rigid_body.Interpolate_Between_States(kinematic_current_state(rigid_body.particle_index),kinematic_next_state(rigid_body.particle_index),time,interpolated_state);
    twist=interpolated_state.twist;
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id)
{
    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(id);
    if(rigid_body_example_velocities.Set_Kinematic_Velocities(twist,time,id)) return;
    RIGID_BODY_STATE<TV> previous_state,current_state;previous_state.time=time-frame_dt;current_state.time=time;
    int new_id=id;
    rigid_body_example_velocities.Set_Kinematic_Positions(previous_state.frame,previous_state.time,new_id);
    rigid_body_example_velocities.Set_Kinematic_Positions(current_state.frame,current_state.time,new_id);
    rigid_body->Compute_Velocity_Between_States(previous_state,current_state,current_state);
    twist=current_state.twist;
}
//#####################################################################
// Function Get_Current_Kinematic_Keyframes
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Get_Current_Kinematic_Keyframes(const T dt,const T time)
{
    if(!use_kinematic_keyframes) return;
    kinematic_current_state.Remove_All();kinematic_next_state.Remove_All();
    kinematic_current_state.Resize(rigid_body_collection.rigid_body_particles.Size());
    kinematic_next_state.Resize(rigid_body_collection.rigid_body_particles.Size());
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=rigid_body_collection.kinematic_rigid_bodies(i);
        kinematic_current_state(p).time=time;kinematic_next_state(p).time=time+dt;
        rigid_body_example_velocities.Set_Kinematic_Positions(kinematic_current_state(p).frame,time,p);
        rigid_body_example_velocities.Set_Kinematic_Positions(kinematic_next_state(p).frame,time+dt,p);
        Set_Kinematic_Velocities(kinematic_current_state(p).twist,dt,time,p);
        Set_Kinematic_Velocities(kinematic_next_state(p).twist,dt,time+dt,p);}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(FRAME<TV>& frame,const T time,const int id)
{
    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(id);
    int new_id=id;
    int index=rigid_body->particle_index;
    if(!use_kinematic_keyframes){
        rigid_body_example_velocities.Set_Kinematic_Positions(frame,time,new_id);
        return;}
    RIGID_BODY_STATE<TV> interpolated_state;
    rigid_body->Interpolate_Between_States(kinematic_current_state(index),kinematic_next_state(index),time,interpolated_state);
    frame=interpolated_state.frame;
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time)
{
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){
        int p=rigid_body_collection.kinematic_rigid_bodies(i);Set_External_Positions(frame(p),time,p);}
    rigid_body_example_velocities.Set_External_Positions(frame,time);
}
//#####################################################################
// Function Reset_Kinematic_Rigid_Bodies
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Reset_Kinematic_Rigid_Bodies(const T time)
{
    // Move kinematic bodies to their position at given time
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_collection.kinematic_rigid_bodies(i));
        Set_External_Positions(rigid_body.Frame(),time,rigid_body.particle_index);
        Set_External_Velocities(rigid_body.Twist(),time,rigid_body.particle_index);}
}
//#####################################################################
namespace PhysBAM{
template class KINEMATIC_EVOLUTION<VECTOR<float,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,3> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,3> >;
}
