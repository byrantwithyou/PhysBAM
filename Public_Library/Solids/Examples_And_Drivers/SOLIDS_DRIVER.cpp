//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_DRIVER
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Utilities/INTERRUPTS.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_DRIVER<TV>::
SOLIDS_DRIVER(SOLIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input),project_at_frame_boundaries(true),last_dt(0),restart_dt(0),reset_with_restart(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_DRIVER<TV>::
~SOLIDS_DRIVER()
{}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files(example.first_frame);}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Preprocess_Frame(current_frame+1);
        example.solids_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1)
            PHYSBAM_DEBUG_WRITE_SUBSTEP("END Frame %d",example.write_substeps_level,current_frame);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Initialize()
{
    if(example.auto_restart){
        std::string last_frame_file=example.output_directory+"/common/last_frame";
        int last_frame;Read_From_Text_File(last_frame_file,last_frame);
        example.restart=true;example.restart_frame=last_frame;
        LOG::cout<<"Auto Restart from frame "<<last_frame<<" (from file "<<last_frame_file<<")"<<std::endl;}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    example.solid_body_collection.deformable_body_collection.mpi_solids=example.solid_body_collection.deformable_body_collection.mpi_solids;
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();

    example.After_Initialization();

    solids_evolution.time=time;

    if(!solids_parameters.fracture && !example.use_melting) // fracture and melting initialize collisions in Initialize_Bodies
        example.solid_body_collection.deformable_body_collection.Initialize(solids_parameters.triangle_collision_parameters);

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);
}
//#####################################################################
// Function Rigid_Cluster_Fracture
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Rigid_Cluster_Fracture(const T dt_full_advance,const T dt_cfl,const int substep)
{
    Check_For_Interrupts(); // see if keyboard or other interrupts are waiting
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=example.solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    ARRAY<int> active_clusters;

    if(rigid_bindings.callbacks && ((substep-1)%example.solids_parameters.rigid_cluster_fracture_frequency)==0 && rigid_bindings.Size()){
        T dt=min(dt_cfl*example.solids_parameters.rigid_cluster_fracture_frequency,dt_full_advance);
        // TODO update example.fluids_parameters.collision_bodies_affecting_fluid for Deactivate_And_Return_Clusters
        rigid_bindings.Deactivate_And_Return_Clusters(active_clusters);
        example.solid_body_collection.Update_Simulated_Particles();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Before declustered evolution",1);

        rigid_bindings.callbacks->Pre_Advance_Unclustered(dt,time);
        example.solids_evolution->kinematic_evolution.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame,time);
        example.solids_evolution->kinematic_evolution.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist,time,time);
        example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
        solids_evolution.Advance_One_Time_Step_Position(dt,time,true);
        rigid_bindings.callbacks->Post_Advance_Unclustered(dt,time);
        rigid_bindings.callbacks->Compute_New_Clusters_Based_On_Unclustered_Strain();

        PHYSBAM_DEBUG_WRITE_SUBSTEP("After declustered evolution",1);
        NEWMARK_EVOLUTION<TV>& newmark_evolution=dynamic_cast<NEWMARK_EVOLUTION<TV>&>(*example.solids_evolution);
        example.solids_evolution->Restore_Position_After_Hypothetical_Position_Evolution(newmark_evolution.X_save,newmark_evolution.rigid_frame_save);
        rigid_bindings.Reactivate_Bindings(active_clusters);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("After restore",1);

        rigid_bindings.callbacks->Create_New_Clusters();
        example.solid_body_collection.Update_Simulated_Particles();

        Setup_Solids(time,substep); // resetup solids before evolution.
    }
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    T dt_full_advance=target_time-time;

    example.solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        Setup_Solids(time,substep);
        T dt=Compute_Dt(time,target_time,done);
        
        Rigid_Cluster_Fracture(dt_full_advance,dt,substep);

        example.Preprocess_Substep(dt,time);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("solid position update",1);
        Solid_Position_Update(dt,substep);/*S1*/

        PHYSBAM_DEBUG_WRITE_SUBSTEP("solid velocity update",1);
        Solid_Velocity_Update(dt,substep,done);/*S2*/

        example.Postprocess_Substep(dt,time);

        last_dt=restart_dt?restart_dt:dt;time+=last_dt;restart_dt=0;

        PHYSBAM_DEBUG_WRITE_SUBSTEP("END Substep %d",0,substep);}
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Setup_Solids(const T time,const int substep)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
        solids_evolution_callbacks->Self_Collisions_Begin_Callback(time,substep);
        solids_parameters.triangle_collision_parameters.repulsion_pair_update_count=0;
        example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Save_Self_Collision_Free_State();
        if((solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count++)%solids_parameters.triangle_collision_parameters.topological_hierarchy_build_frequency==0){
            LOG::SCOPE scope("hierarchybuild","Building Hierarchy Topology");
            example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Topological_Structure_Of_Hierarchies();}
        solids_parameters.triangle_collision_parameters.self_collision_free_time=time;}

    solids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
    if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects) // TODO - ANDY - why is this needed??? TODO: move this to the right places inside solids evolution 
        example.solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
            solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);
    example.Update_Time_Varying_Material_Properties(time);
}
//#####################################################################
// Function Solid_Position_Update
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Solid_Position_Update(const T dt,const int substep)
{
    Check_For_Interrupts(); // see if keyboard or other interrupts are waiting
    LOG::SCOPE scope("solids position update");

    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.solid_body_collection.deformable_body_collection;

    if((solids_parameters.triangle_collision_parameters.repulsion_pair_update_count++)%solids_parameters.triangle_collision_parameters.repulsion_pair_update_frequency==0){
        example.solid_body_collection.deformable_body_collection.triangle_repulsions.Update_Faces_And_Hierarchies_With_Collision_Free_Positions(&deformable_body_collection.particles.X);
        example.solid_body_collection.deformable_body_collection.triangle_repulsions.Compute_Interaction_Pairs(deformable_body_collection.particles.X);}
    solids_evolution.kinematic_evolution.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame,time);
    solids_evolution.kinematic_evolution.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist,time,time);
    example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
    solids_evolution.Advance_One_Time_Step_Position(dt,time,true);

    if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
        LOG::SCOPE scope("adjust velocity for self repulsion and self collisions");
        int repulsions,collisions_found;
        solids_evolution.Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(dt,time,repulsions,collisions_found,false);
        solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;}
    // Exchange solid positions back to fluid nodes so that they can figure out effective velocity and do collidable advection
    PHYSBAM_DEBUG_WRITE_SUBSTEP("solid position updated",1);
}
//#####################################################################
// Function Solid_Velocity_Update
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Solid_Velocity_Update(const T dt,const int substep,const bool done)
{
    Check_For_Interrupts(); // see if keyboard or other interrupts are waiting
    LOG::SCOPE scope("solids velocity update");

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;
    solids_evolution.Advance_One_Time_Step_Velocity(dt,time,true);
    solids_evolution.time+=dt;
    solids_evolution_callbacks->Postprocess_Solids_Substep(solids_evolution.time,substep);
    solids_evolution_callbacks->Apply_Constraints(dt,solids_evolution.time);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Postprocess_Frame(const int frame)
{
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;

    example.Postprocess_Frame(frame);
    solids_evolution.Postprocess_Frame(frame);

    if(solids_evolution.solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics)
        solids_evolution.rigid_body_collisions->Print_Interpenetration_Statistics();
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR SOLIDS_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    T fluids_dt=FLT_MAX;
    // solids dt
    T solids_dt=FLT_MAX;
    if(solids_evolution.Use_CFL()) solids_dt=min(solids_dt,example.solid_body_collection.CFL(solids_parameters.verbose_dt));
    solids_dt=min(solids_dt,solids_evolution_callbacks->Constraints_CFL());
    if(solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies)
        solids_dt=min(solids_dt,example.solid_body_collection.rigid_body_collection.CFL_Rigid(solids_parameters.rigid_body_evolution_parameters,solids_parameters.verbose_dt));
    solids_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);
    if(example.solid_body_collection.deformable_body_collection.mpi_solids)
        solids_dt=example.solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Min_Global(solids_dt);

    if(example.fixed_dt){fluids_dt=example.fixed_dt;solids_dt=example.fixed_dt;}
    if(example.max_dt){fluids_dt=min(fluids_dt,example.max_dt);solids_dt=min(solids_dt,example.max_dt);}
    T dt=min(fluids_dt,solids_dt);
    LOG::cout<<"dt = solids_dt = "<<dt<<std::endl;
    if(example.abort_when_dt_below && dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(LOG::sprintf("dt too small (%g < %g)",dt,example.abort_when_dt_below));
    done=false;
    SOLIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,solids_parameters.min_dt);
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("writing output files");
    Create_Directory(example.output_directory);
    Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);

    example.Write_Output_Files(frame);

    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Preprocess_Frame(const int frame)
{
    if(example.substeps_delay_frame==frame){example.Set_Write_Substeps_Level(example.substeps_delay_level);output_number=frame-1;}
    example.Preprocess_Frame(frame);
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_DRIVER<VECTOR<float,1> >;
template class SOLIDS_DRIVER<VECTOR<float,2> >;
template class SOLIDS_DRIVER<VECTOR<float,3> >;
template class SOLIDS_DRIVER<VECTOR<double,1> >;
template class SOLIDS_DRIVER<VECTOR<double,2> >;
template class SOLIDS_DRIVER<VECTOR<double,3> >;
}
