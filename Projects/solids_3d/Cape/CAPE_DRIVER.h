//#####################################################################
// Copyright 2002, 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CAPE_DRIVER
//##################################################################### 
#ifndef __CAPE_DRIVER__
#define __CAPE_DRIVER__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include "CAPE_EXAMPLE.h"
#include "time.h"
#include <Geometry/IMPLICIT_SURFACE_LIST.h>
#include <Geometry/TRIANGULATED_SURFACE_LIST.h>
namespace PhysBAM{

template<class T>
class CAPE_DRIVER:public EXTERNAL_FORCES_AND_VELOCITIES<T,VECTOR_3D<T> >
{
public:
    CAPE_EXAMPLE<T>& example;
    T time;
    T real_time_start,real_time_last;
    TRIANGLE_MESH triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > particles;
    TRIANGULATED_SURFACE<T> triangulated_surface;
    DEFORMABLE_OBJECT<T,VECTOR_3D<T> > deformable_object;
    TRIANGLE_COLLISIONS<T>* triangle_collisions;
    ARRAY<RIGID_BODY<TV>*> rigid_bodies;
    bool collisions_on;
    ARRAY<bool> enforce_collision_velocity; // for external forces and velocities
    ARRAY<VECTOR_3D<T> > collision_normal; // for external forces and velocities
    ARRAY<T> collision_velocity;
    ARRAY<VECTOR_3D<T> > V_save;  // save veolocity in time stepping
    
    TRIANGLE_MESH collision_mesh;
    TRIANGULATED_SURFACE<T> collision_surface;

    CAPE_DRIVER(CAPE_EXAMPLE<T>& example_input)
        :example(example_input),
        triangulated_surface(triangle_mesh,particles),deformable_object(particles),triangle_collisions(0),collisions_on(false),
        collision_surface(collision_mesh,particles)
    {  
        Initialize();
    }

    ~CAPE_DRIVER()
    {delete triangle_collisions;for(int i=0;i<rigid_bodies.m;i++) delete rigid_bodies(i);}

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{  
    // initialize deformable triangulated surface
    example.collision_surface=&collision_surface;
    example.Initialize_Deformable_Object(deformable_object,triangulated_surface);
    
    // initialize collisions
    example.Initialize_Triangulated_Surface_Collisions(collision_surface,triangle_collisions);
    if(!example.allow_intersections && example.check_initial_mesh_for_self_intersection){
        std::cout << "Checking mesh for self intersection... " << std::endl;
        if(collision_surface.Check_For_Self_Intersection(triangle_collisions->small_number)){std::cout << "1) intersections found!" << std::endl;exit(1);}
        else collision_surface.Find_First_Self_Intersection(triangle_collisions->small_number,10);}

    // external forces and velocities
    deformable_object.Set_External_Forces_And_Velocities(*this);
    enforce_collision_velocity.Resize(particles.array_collection->Size());
    collision_normal.Resize(particles.array_collection->Size());
    collision_velocity.Resize(particles.array_collection->Size());
    
    example.cape_collisions->Enable_Constraints(enforce_collision_velocity,collision_normal,collision_velocity);

    // initialize time
    time=example.initial_time+example.restart_step_number/example.frame_rate;

    // restart if necessary
    if(example.restart_step_number){
        deformable_object.CFL();   // needs to be initialized before positions change
        example.Read_Data_Files(triangulated_surface,time,example.restart_step_number);}
    
    // rigid bodies
    example.Initialize_Rigid_Bodies(rigid_bodies);
    for(int k=0;k<rigid_bodies.m;k++) rigid_bodies(k)->implicit_surface->Compute_Cell_Minimum_And_Maximum();
    if(example.restart_step_number){example.Update_Rigid_Body_Positions(rigid_bodies,time);example.Update_Rigid_Body_Velocities(rigid_bodies,time);}

    // output
    if(!example.restart_step_number) example.Write_Data_Files(triangulated_surface,rigid_bodies,time,0);
    else example.Write_Data_Files(triangulated_surface,rigid_bodies,time,1000);
    
    real_time_start=real_time_last=(T)::time(0);
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
void Execute_Main_Program()
{          
    if(collision_surface.Check_For_Self_Intersection(triangle_collisions->small_number)){std::cout << "3) intersections found!" << std::endl;exit(1);}

    T time_per_frame=(T)1/example.frame_rate;int k=example.restart_step_number;
    while(time<example.final_time){
        Advance_One_Frame(time+time_per_frame,++k);
        std::cout << "TIME STEP = " << k << ", TIME = " << " " << time << std::endl << std::endl;
        example.Write_Data_Files(triangulated_surface,rigid_bodies,time,k);
        if(example.perform_self_collision && !example.allow_intersections){
            std::cout << "Checking mesh for self intersection... " << std::endl << std::endl;
            if(collision_surface.Check_For_Self_Intersection(triangle_collisions->small_number)){std::cout << "2) intersections found!" << std::endl;exit(1);}}}
}
//#####################################################################
// Function Advance_One_Frame
//#####################################################################
void Advance_One_Frame(const T end_time,const int frame)
{
    // prepare for force computation
    deformable_object.Update_Position_Based_State(); 
    
    static int total_loops=1,collision_resolved_steps=0;bool done=false;int max_loops=64,substep=0;
    while(!done){substep++;
        if(example.perform_self_collision){
            T time_old=time;
            triangle_collisions->Initialize_Collision_Free_State();bool collisions_resolved=false;
            while(!collisions_resolved){
                std::cout  <<  "                                TOTAL_LOOPS = " << total_loops << std::endl;
                int loops=0;while(loops < total_loops && !done){loops++;
                    T dt=deformable_object.CFL();//if(dt<1e-5)dt=1e-5;
                    if(time+dt >= end_time){dt=end_time-time;done=true;}else if(time+2*dt >= end_time) dt=(T).51*(end_time-time);
                    Advance_One_Time_Step(dt,true);
                    example.Write_Data_Files(triangulated_surface,rigid_bodies,time,frame+1);
                    std::cout << "loops = " << loops << ", dt = " << " " << dt << std::endl;}
                int repulsions=0,collisions=0;bool exit_early=true;if(total_loops == 1) exit_early=false;
                bool self_interaction=triangle_collisions->Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(time-time_old,repulsions,collisions,exit_early);
                if(self_interaction) deformable_object.Update_Position_Based_State();
                if(collisions < 0){collision_resolved_steps=0;total_loops/=2;done=false;time=time_old;} // repeat with less dynamics solver loops
                else{collisions_resolved=true;collision_resolved_steps++;if(collision_resolved_steps%3 == 0 && total_loops < example.max_collision_loops) total_loops*=2;}}}
        else{
            T dt=deformable_object.CFL();if(time+dt >= end_time){dt=end_time-time;done=true;}else if(time+2*dt >= end_time) dt=(T).51*(end_time-time);
            Advance_One_Time_Step(dt,true);
            std::cout << "dt = " << dt << "\n";}
        std::cout << "substep = " << substep << ", time = " << time << std::endl;}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
void Advance_One_Time_Step(const T dt,const bool verbose)
{      
    std::cout<<"maximum speed at beginning of time step: "<<particles.Maximum_Speed(0)<<"\n";
    
    // save velocity for later trapezoidal rule
    V_save.Resize(particles.V.m);   // THIS IS A HACK TO GET THE CAPE EXAMPLE TO WORK
    ARRAY<VECTOR_3D<T> >::copy_up_to(particles.V.array,V_save,particles.array_collection->Size());

    // update V implicitly to time n+1/2
    if(!deformable_object.Backward_Euler_Step_Velocity(dt/2,time,example.cg_tolerance,example.cg_iterations,false)){
        if(verbose) std::cout << "BACKWARD EULER FAILED!" << std::endl;deformable_object.Predictor_Corrector_Integrate_Velocity(time,time+dt/2);} 
   
    // update position
    deformable_object.Euler_Step_Position(dt);
    ARRAY<VECTOR_3D<T> >::Exchange_Arrays(V_save,particles.V.array); // don't care about V_save, but Exchange_Arrays is faster than copy (since it does pointers)

    // rigid body collisions
    example.Update_Rigid_Body_Velocities(rigid_bodies,time+dt);example.Update_Rigid_Body_Positions(rigid_bodies,time+dt);
    example.Read_Fractional_Buddha(time+dt);
    collisions_on=true;Adjust_Nodes_For_Rigid_Body_Collisions(particles.X,particles.V,dt);example.Cape_Set_Positions(time+dt);

    // positions just changed, so update again
    deformable_object.Update_Position_Based_State();

    // finish velocity update
    std::cout<<"maximum speed before Euler step: "<<particles.Maximum_Speed(0)<<"\n";
    deformable_object.Euler_Step_Velocity(dt/2,time);
/*    if(!deformable_object.Backward_Euler_Step_Velocity(dt/2,time,example.cg_tolerance,example.cg_iterations,false,verbose)){
        if(verbose) std::cout << "BACKWARD EULER FAILED!" << std::endl;deformable_object.Predictor_Corrector_Integrate_Velocity(time,time+dt/2);} */
    std::cout<<"maximum speed after Euler step: "<<particles.Maximum_Speed(0)<<"\n";
    time+=dt/2;
    if(!deformable_object.Backward_Euler_Step_Velocity(dt/2,time,example.cg_tolerance,example.cg_iterations,false)){
        if(verbose) std::cout << "BACKWARD EULER FAILED!" << std::endl;deformable_object.Predictor_Corrector_Integrate_Velocity(time,time+dt/2);}  
    time+=dt/2;

    // reset collision velocities to zero
    collisions_on=false;ARRAY<bool>::copy(false,enforce_collision_velocity);
    
    T real_time=(T)::time(0);
    std::cout<<"elapsed time: "<<real_time-real_time_start<<" (+ "<<real_time-real_time_last<<")\n";
    real_time_last=real_time;
}
//#####################################################################
// Function Adjust_Nodes_For_Rigid_Body_Collisions
//#####################################################################
void Adjust_Nodes_For_Rigid_Body_Collisions(ARRAY<VECTOR_3D<T> >& X,ARRAY<VECTOR_3D<T> >& V,const T dt)
{
    // Should be accelerated by using the triangle_hierarchy to eliminate tests on parts of the surface that aren't close to rigid bodies.
    int interactions=0;T depth;
    for(int p=0;p<particles.array_collection->Size();p++) for(int r=0;r<rigid_bodies.m;r++) if(rigid_bodies(r)->Implicit_Surface_Lazy_Inside_And_Value(X(p),depth)){
        depth=max((T)0,-depth)+triangle_collisions->collision_thickness;
        rigid_bodies(r)->Adjust_Point_For_Rigid_Body_Collision(X(p),V(p),depth,dt);
        enforce_collision_velocity(p)=true; // for external forces and velocities
        collision_normal(p)=rigid_bodies(r)->Implicit_Surface_Normal(X(p)); // for external forces and velocities
        collision_velocity(p)=VECTOR_3D<T>::Dot_Product(V(p),collision_normal(p)); // for external forces and velocities
        interactions++;}
    if(interactions) std::cout << "rigid body collisions = " << interactions << std::endl;

    example.cape_collisions->Process_Collisions(dt);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time=1)
{
    if(collisions_on)
        for(int i=0;i<enforce_collision_velocity.m;i++) if(enforce_collision_velocity(i))
            V(i)+=(collision_velocity(i)-VECTOR_3D<T>::Dot_Product(V(i),collision_normal(i)))*collision_normal(i);
    example.Set_External_Velocities(V,time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time=1)
{
    if(collisions_on)
        for(int i=0;i<enforce_collision_velocity.m;i++) if(enforce_collision_velocity(i))
            V(i)-=VECTOR_3D<T>::Dot_Product(V(i),collision_normal(i))*collision_normal(i);
    example.Zero_Out_Enslaved_Velocity_Nodes(V,time);
}
//#####################################################################
};
}
#endif
