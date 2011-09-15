//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODIES_2D_DRIVER
//#####################################################################
#ifndef __RIGID_BODIES_2D_DRIVER__
#define __RIGID_BODIES_2D_DRIVER__    

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include "RIGID_BODIES_2D_EXAMPLE.h"
#include <Rigid_Bodies/RIGID_BODY_COLLISIONS_2D.h>
#include <Rigid_Bodies/RIGID_BODY_EVOLUTION_2D.h>
#include <Solids_And_Fluids/SIMULATION_DRIVER.h>
#include <stdlib.h>
#include <string.h>

namespace PhysBAM{

template<class T>
class RIGID_BODIES_2D_DRIVER : public SIMULATION_DRIVER<T>
{
public:
    RIGID_BODIES_2D_EXAMPLE<T>& example;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies;
    RIGID_BODY_COLLISIONS_2D<T>* collisions;
    RIGID_BODY_EVOLUTION_2D<T>* rigid_body_evolution;

    RIGID_BODIES_2D_DRIVER(RIGID_BODIES_2D_EXAMPLE<T>& example_input)
        :SIMULATION_DRIVER<T>(example_input),example(example_input),rigid_bodies(example.rigid_body_parameters.list.rigid_bodies),
        collisions(0),rigid_body_evolution(0)
    {}

    ~RIGID_BODIES_2D_DRIVER()
    {delete collisions;delete rigid_body_evolution;}

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    SIMULATION_DRIVER<T>::Initialize();

    // initialize rigid bodies
    example.Initialize_Rigid_Bodies();
    example.Initialize_Rigid_Body_Forces();

    if(example.restart) example.Load_Restart_Data(example.restart_frame);
    else Write_Output_Files(current_frame);
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
    rigid_body_evolution=new RIGID_BODY_EVOLUTION_2D<T>(rigid_bodies,*collisions);
    rigid_body_evolution->Set_CFL_Number(example.cfl);
    rigid_body_evolution->Set_Max_Rotation_Per_Time_Step(example.max_rotation_per_time_step);
    rigid_body_evolution->Set_Max_Linear_Movement_Fraction_Per_Time_Step(example.max_linear_movement_fraction_per_time_step);
    rigid_body_evolution->Set_Minimum_And_Maximum_Time_Step(0,.1/example.frame_rate);
    if(example.artificial_maximum_speed) rigid_body_evolution->Set_Artificial_Maximum_Speed(example.artificial_maximum_speed);

    // precompute normals - faster, but less accurate
    //RIGID_BODY_COLLISIONS_2D<T>::Precompute_Normals(rigid_bodies);

    // check for bad initial data
    if (example.extra_verbose) rigid_body_evolution->collisions.Print_Interpenetration_Statistics();
    else if(!collisions->Check_For_Any_Interpenetration()) std::cout << "No initial interpenetration" << std::endl;
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
void Advance_To_Target_Time(const T target_time)
{   
    collisions->Initialize_Data_Structures(); // Call here to handle possibly changed number of rigid bodies
    bool done=false;for(int substep=1;!done;substep++){
        example.Update_Animated_Parameters(time);
        T dt_cfl=rigid_body_evolution->CFL(example.verbose_dt),dt=dt_cfl;
        Clamp_Time_Step_With_Target_Time(target_time,dt,done);
        std::cout << "substep = " << substep << ", time = " << time << ", dt = " << dt << " (cfl = " << dt_cfl << ")" << std::endl;
        rigid_body_evolution->Advance_One_Time_Step(dt,time);
        time+=dt;
    }
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame)
{
    SIMULATION_DRIVER<T>::Postprocess_Frame(frame);
    
    if (example.extra_verbose) rigid_body_evolution->collisions.Print_Interpenetration_Statistics();
    if (example.is_arb) example.Calculate_ARB_Impulse();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame)
{
    example.Write_Data_Files(frame);
}
//#####################################################################
};
}
#endif
