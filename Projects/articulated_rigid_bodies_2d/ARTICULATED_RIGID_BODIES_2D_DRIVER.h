//#####################################################################
// Copyright 2004, Rachel Weinstein, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODIES_2D_DRIVER
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODIES_2D_DRIVER__
#define __ARTICULATED_RIGID_BODIES_2D_DRIVER__    

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include "ARTICULATED_RIGID_BODIES_2D_EXAMPLE.h"
#include <Rigid_Bodies/RIGID_BODY_COLLISIONS_2D.h>
#include <Rigid_Bodies/RIGID_BODY_EVOLUTION_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_2D.h>
#include <stdlib.h>
#include <string.h>

namespace PhysBAM{

template<class T,class RW>
class ARTICULATED_RIGID_BODIES_2D_DRIVER : public SOLIDS_FLUIDS_DRIVER_2D<T,RW>
{
public:
    ARTICULATED_RIGID_BODIES_2D_EXAMPLE<T,RW>& example;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies;
    RIGID_BODY_COLLISIONS_2D<T>* collisions;
    RIGID_BODY_EVOLUTION_2D<T>* rigid_body_evolution;
    int tempFrame;

    ARTICULATED_RIGID_BODIES_2D_DRIVER(ARTICULATED_RIGID_BODIES_2D_EXAMPLE<T,RW>& example_input)
        :SOLIDS_FLUIDS_DRIVER_2D<T,RW>(example_input),example(example_input),rigid_bodies(example.rigid_body_list->rigid_bodies),
        collisions(0),rigid_body_evolution(0),tempFrame(1)
    {}

    ~ARTICULATED_RIGID_BODIES_2D_DRIVER()
    {delete collisions;delete rigid_body_evolution;}

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    SOLIDS_FLUIDS_DRIVER_2D<T,RW>::Initialize();

    // initialize rigid bodies
    example.Initialize_Rigid_Bodies();
    example.Initialize_Rigid_Body_Forces();

    if(example.restart) example.Load_Restart_Data(example.restart_frame);
    RIGID_BODY_COLLISIONS_2D<T>::Adjust_Bounding_Boxes(rigid_bodies);
    
    // compute voxel size for spatial partition
    T voxel_size=0;
    if(example.spatial_partition_based_on_scene_size){
        VECTOR_2D<T> scene_box_size=RIGID_BODY_COLLISIONS_2D<T>::Scene_Bounding_Box(rigid_bodies).Size();
        voxel_size=one_third*(scene_box_size.x+scene_box_size.y)/example.spatial_partition_number_of_cells;}
    else if(example.spatial_partition_based_on_object_size){
        if(example.spatial_partition_with_max_size) voxel_size=RIGID_BODY_COLLISIONS_2D<T>::Maximum_Bounding_Box_Size(rigid_bodies);
        else voxel_size=RIGID_BODY_COLLISIONS_2D<T>::Average_Bounding_Box_Size(rigid_bodies);
        voxel_size*=4;} // just to make it bigger
    
    // initialize collisions
    collisions=new RIGID_BODY_COLLISIONS_2D<T>(rigid_bodies,voxel_size);
    RIGID_BODY_COLLISIONS_CALLBACKS<T>* cback=new RIGID_BODY_COLLISIONS_CALLBACKS<T>();
    collisions->Set_Rigid_Body_Collisions_Callbacks(*cback);
    collisions->verbose=example.verbose;
#if 0
    if(example.use_particle_partition){
        collisions->intersections.Use_Particle_Partition(true,example.particle_partition_size,example.particle_partition_size,
                                                                                   example.particle_partition_size);
        if(example.use_particle_partition_center_phi_test) collisions->intersections.Use_Particle_Partition_Center_Phi_Test();}
    if(example.use_triangle_hierarchy){
        collisions->intersections.Use_Triangle_Hierarchy();
        if(example.use_triangle_hierarchy_center_phi_test) collisions->intersections.Use_Triangle_Hierarchy_Center_Phi_Test();
        if(example.use_edge_intersection) collisions->intersections.Use_Edge_Intersection();}
#endif

    // initialize dynamics
    rigid_body_evolution=new RIGID_BODY_EVOLUTION_2D<T>(&example.arb,*collisions); // need to fix to include ARB
    rigid_body_evolution->Set_CFL_Number(example.solids_parameters.cfl);
    rigid_body_evolution->Set_Max_Rotation_Per_Time_Step(example.max_rotation_per_time_step);
    rigid_body_evolution->Set_Max_Linear_Movement_Fraction_Per_Time_Step(example.max_linear_movement_fraction_per_time_step);
    rigid_body_evolution->Set_Minimum_And_Maximum_Time_Step(0,.1/example.frame_rate);
    if(example.artificial_maximum_speed) rigid_body_evolution->Set_Artificial_Maximum_Speed(example.artificial_maximum_speed);

    // precompute normals - faster, but less accurate
    //RIGID_BODY_COLLISIONS_2D<T>::Precompute_Normals(rigid_bodies);

    example.Write_OpenGL_Hints();

    // check for bad initial data
    if (example.extra_verbose) rigid_body_evolution->collisions.Print_Interpenetration_Statistics();
    else if(!collisions->Check_For_Any_Interpenetration()) std::cout << "No initial interpenetration" << std::endl;
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
void Advance_To_Target_Time(const T target_time)
{   T dt;
    collisions->Initialize_Data_Structures(); // Call here to handle possibly changed number of rigid bodies
    bool done=false;for(int substep=1;!done;substep++){
//        example.Update_Animated_Parameters(time);
        T dt_cfl=rigid_body_evolution->CFL(example.verbose_dt);
        dt=dt_cfl;
//        Clamp_Time_Step_With_Target_Time(target_time,dt,done);
        SOLIDS_FLUIDS_EXAMPLE<T>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        std::cout << "substep = " << substep << ", time = " << time << ", dt = " << dt << " (cfl = " << dt_cfl << ")" << std::endl;
        rigid_body_evolution->Advance_One_Time_Step(dt,time);

        //example.Write_Data_Files(tempFrame++);

        //Write_Output_Files(tempFrame++);
        time+=dt;
        //this should be after one time step
        //exit(1);
    }
//    example.Write_Data_Files(frame);
  

}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame)
{
    std::cout<<"frame: "<<frame<<std::endl;
    SOLIDS_FLUIDS_DRIVER_2D<T,RW>::Postprocess_Frame(frame);
    
    if (example.extra_verbose) rigid_body_evolution->collisions.Print_Interpenetration_Statistics();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame)
{
    example.Write_Data_Files(frame);
    example.Write_Articulation_Points(frame);
//    example.Write_Data_Files(tempFrame++);
}
//#####################################################################
};
}
#endif
