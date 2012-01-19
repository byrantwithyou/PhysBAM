//#####################################################################
// Copyright 2004, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "KINEMATIC_EXAMPLE.h"
using namespace PhysBAM;
template KINEMATIC_EXAMPLE<double>;
template KINEMATIC_EXAMPLE<float>;
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    std::cout << "*** Using parameter " << parameter << std::endl;

    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    last_frame=600;

    if(parameter==1){
        // SPHERE
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(0,5,0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("sphere (coeff 1)");
        Append_Diffuse_Hints(color1);

        // SPHERE
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>((T)2.5,3,0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("sphere (coeff 1)");
        Append_Diffuse_Hints(color1);

        // SPHERE
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(5,4,0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("sphere (coeff 1)");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body("box");
        rigid_body->position=VECTOR_3D<T>(10,10,0);
        rigid_body->orientation=QUATERNION<T>((T).1,VECTOR_3D<T>(1,2,3));
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("box");
        Append_Diffuse_Hints(color1);

        // SPHERE
        position_curve1.Add_Control_Point(0,VECTOR_3D<T>(-5,1,0.1));
        position_curve1.Add_Control_Point(10,VECTOR_3D<T>(10,1,0.1));
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=position_curve1.Value(0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("sphere (coeff 0.5)");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
        Append_Diffuse_Hints(color1);

        position_curve2.Add_Control_Point(0,VECTOR_3D<T>(10,5,0));
        position_curve2.Add_Control_Point(10,VECTOR_3D<T>(10,15,0));
        rigid_body=Initialize_Rigid_Body("box",3);
        rigid_body->position=position_curve2.Value(0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("sphere (coeff 0.5)");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
        Append_Diffuse_Hints(color1);
    }
    else if(parameter==2){
        rigid_body=Initialize_Rigid_Body("box");
        rigid_body->position=VECTOR_3D<T>(10,10,0);
        rigid_body->orientation=QUATERNION<T>((T).1,VECTOR_3D<T>(1,2,3));
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("box");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
        Append_Diffuse_Hints(color1);
    }
    else if(parameter==3){
        T baseboxsize=2;
        T dropboxsize=1.2;
        T plankscale=1;
        T dropheight=70;
        T smallboxsize=0.5;
        T smallboxmass=1;
        T offset = 0.05;

        T stack_epsilon = 0.3;
        T stack_mu = 0.5;

#if 0
        time1=(T)10/24;
        time2=(T)200/24;
        time3=(T)350/24;
        time4=(T)400/24;

        const char *boxfile = "box";

        T boxsize1=1.0;
        T boxsize2=0.5;

        rigid_body=Initialize_Rigid_Body(boxfile,boxsize1);
        rigid_body->position=VECTOR_3D<T>(0,2*baseboxsize+boxsize1,0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1a");
        Append_Diffuse_Hints(color1);

#if 1
        rigid_body=Initialize_Rigid_Body(boxfile,boxsize2);
        rigid_body->position=VECTOR_3D<T>(0,2*baseboxsize+2*boxsize1+boxsize2,0);
        //rigid_body->orientation=QUATERNION<T>(pi/4,VECTOR_3D<T>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1b");
        Append_Diffuse_Hints(color1);
#endif

        rigid_body=Initialize_Rigid_Body(boxfile,baseboxsize);
        rigid_body->Set_Coefficient_Of_Restitution((T).1);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Name("base box");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
        kinematic_body_index=rigid_body_parameters.list.rigid_bodies.m;
        kinematic_body_vel=-5;
        Append_Diffuse_Hints(color1);
#else
        const char *boxfile = "subdivided_box";

        time1=(T)50/24;
        time2=(T)200/24;
        time3=(T)350/24;
        time4=(T)400/24;

        rigid_body=Initialize_Rigid_Body(boxfile,baseboxsize);
        rigid_body->Set_Coefficient_Of_Restitution((T).1);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Name("base box");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
        kinematic_body_index=rigid_body_parameters.list.rigid_bodies.m;
        kinematic_body_vel=5;
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(-0.65,baseboxsize+4,0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1a");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0.65,baseboxsize+5,0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1b");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(offset,baseboxsize+7,0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 2");
        Append_Diffuse_Hints(color1);

#if 0
        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0,baseboxsize+9,offset);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 3");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(-offset,baseboxsize+14,-offset);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 4");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0,baseboxsize+17,0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 5");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0,baseboxsize+19,0);
        rigid_body->orientation=QUATERNION<T>(pi/4,VECTOR_3D<T>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 6");
        Append_Diffuse_Hints(color1);
#endif
#endif

        //for(int i=0;i<rigid_body_parameters.list.rigid_bodies.m;i++) if(!rigid_body_parameters.list.rigid_bodies(i)->is_static) rigid_body_parameters.list.rigid_bodies(i)->is_kinematic=true;
    }

#if 0
    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);
#endif

    Update_Animated_Parameters(0);
    for(int i=0;i<rigid_body_parameters.list.rigid_bodies.m;i++){
        rigid_body=rigid_body_parameters.list.rigid_bodies(i);
        if(rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)){
            Set_External_Position_And_Orientation(rigid_body->position,rigid_body->orientation,0,i);
            Set_External_Velocities(rigid_body->velocity,rigid_body->angular_velocity,0,i);}}

    rigid_body_parameters.list.Set_External_Forces_And_Velocities(*this);
}
//#####################################################################
// Function Compute_State
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Compute_State(RIGID_BODY_STATE_3D<T>& state,const T time)
{
    T baseboxsize=2;
    state.time=time;
    if(time<time1){
        state.position=VECTOR_3D<T>(0,baseboxsize,0);
        state.orientation=QUATERNION<T>();
    }
    else if(time<time2){
        state.position=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time-time1),0);
        state.orientation=QUATERNION<T>();
    }
    else if(time<time3){
        state.position=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);
        state.orientation=QUATERNION<T>();
    }
    else if(time<time4){
        state.position=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);
        state.orientation=QUATERNION<T>(-0.4*(time-time3),VECTOR_3D<T>(0,0,1));
    }
    else{
        state.position=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);
        state.orientation=QUATERNION<T>(-0.4*(time4-time3),VECTOR_3D<T>(0,0,1));
    }
}
//#####################################################################
// Function Compute_Velocities
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Compute_Velocities(const RIGID_BODY<TV>& rigid_body,RIGID_BODY_STATE_3D<T>& state1,RIGID_BODY_STATE_3D<T>& state2)
{
    T one_over_dt=1/(state2.time-state1.time);
    state2.velocity=one_over_dt*(state2.position-state1.position);
    QUATERNION<T> quat=state2.orientation*state1.orientation.Inverse();
    state2.angular_velocity=one_over_dt*quat.Rotation_Vector();
    rigid_body.Update_Angular_Momentum(state2);
}
//#####################################################################
// Function Update_Animated_Parameters
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Update_Animated_Parameters(const T time)
{
    if(parameter==1){
        if(time>2)rigid_body_parameters.list.rigid_bodies(6)->is_kinematic=false;
    }
    else if(parameter==2){
        bool &is_kinematic=rigid_body_parameters.list.rigid_bodies(1)->is_kinematic,new_is_kinematic=(time<3*pi/2);
        if(is_kinematic!=new_is_kinematic){is_kinematic=new_is_kinematic;solid_body_collection.Update_Fragments();} // Chaning is_kinematic requires an Update_Fragments
    }
    else if(parameter==3){
        //rigid_body_parameters.list.rigid_bodies(kinematic_body_index)->is_kinematic=(time>(T)40/24);
        //rigid_body_parameters.list.rigid_bodies(kinematic_body_index)->is_static=!rigid_body_parameters.list.rigid_bodies(kinematic_body_index)->is_kinematic;

        T prev_time=time-(T)1/frame_rate,current_time=time,next_time=time+(T)1/frame_rate;
        RIGID_BODY_STATE_3D<T> prev_state;
        Compute_State(prev_state,prev_time);
        Compute_State(current_state,current_time);
        Compute_State(next_state,next_time);
        Compute_Velocities(*rigid_body_parameters.list.rigid_bodies(kinematic_body_index),prev_state,current_state);
        Compute_Velocities(*rigid_body_parameters.list.rigid_bodies(kinematic_body_index),current_state,next_state);
    }
}    
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Set_External_Velocities(VECTOR_3D<T>& V,VECTOR_3D<T>& omega,const T time,const int id_number)
{
    if(!rigid_body_parameters.list.rigid_bodies(id_number)->is_kinematic) return;
    if(parameter==1){
        if(id_number==5){V=position_curve1.Derivative(time);omega=VECTOR_3D<T>();}
        else if(id_number==6){V=position_curve2.Derivative(time);omega=VECTOR_3D<T>();}
    }
    else if(parameter==2){
        if(id_number==1){
            V=VECTOR_3D<T>(10,-10*sin(time),0);
            omega=VECTOR_3D<T>();
        }
    }
    else if(parameter==3){
        if(id_number==kinematic_body_index){
            RIGID_BODY_STATE_3D<T> interpolated_state;
            rigid_body_parameters.list.rigid_bodies(id_number)->Interpolate_Between_States(current_state,next_state,time,interpolated_state);
            V=interpolated_state.velocity;
            omega=interpolated_state.angular_velocity;
        }
    }
}
//#####################################################################
// Function Set_External_Position_And_Orientation
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Set_External_Position_And_Orientation(VECTOR_3D<T>& X,QUATERNION<T>& orientation,const T time,const int id_number)
{
    if(!rigid_body_parameters.list.rigid_bodies(id_number)->is_kinematic) return;
    if(parameter==1){
        if(id_number==5){X=position_curve1.Value(time);orientation=QUATERNION<T>();}
        else if(id_number==6){X=position_curve2.Value(time);orientation=QUATERNION<T>();}
    }
    else if(parameter==2){
        if(id_number==1){
            X=VECTOR_3D<T>(10*time,10*cos(time)+30,0);
            orientation=QUATERNION<T>();
        }
    }
    else if(parameter==3){
        if(id_number==kinematic_body_index){
            RIGID_BODY_STATE_3D<T> interpolated_state;
            rigid_body_parameters.list.rigid_bodies(id_number)->Interpolate_Between_States(current_state,next_state,time,interpolated_state);
            X=interpolated_state.position;
            orientation=interpolated_state.orientation;
        }
    }
}
