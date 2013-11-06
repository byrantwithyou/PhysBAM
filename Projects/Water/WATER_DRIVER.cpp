//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <Tools/Vectors/VECTOR_UTILITIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((WATER_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
WATER_DRIVER(WATER_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_body_collection,example,true),thread_queue(example.thread_queue)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
~WATER_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    
    // initialize collision objects
    kinematic_evolution.Get_Current_Kinematic_Keyframes(1/example.frame_rate,time);
    kinematic_evolution.Set_External_Positions(example.rigid_body_collection.rigid_body_particles.frame,time);
    kinematic_evolution.Set_External_Velocities(example.rigid_body_collection.rigid_body_particles.twist,time,time);

    example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.Particle_Levelset(0).Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    example.face_velocities.Resize(example.mac_grid);

    example.particle_levelset_evolution.time=time;
    example.particle_levelset_evolution.Set_CFL_Number((T).9);

    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    example.incompressible.mpi_grid=example.mpi_grid;
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    example.particle_levelset_evolution.Particle_Levelset(0).mpi_grid=example.mpi_grid;
    if(example.mpi_grid){
        example.boundary=new BOUNDARY_MPI<TV>(example.mpi_grid,example.boundary_scalar);
        example.phi_boundary=new BOUNDARY_MPI<TV>(example.mpi_grid,example.phi_boundary_water);
        example.particle_levelset_evolution.Particle_Levelset(0).last_unique_particle_id=example.mpi_grid->rank*30000000;}
    else{
        example.boundary=&example.boundary_scalar;
        example.phi_boundary=&example.phi_boundary_water;}

    if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<TV> *refine=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<TV>*>(&example.projection)){
        refine->boundary=example.boundary;
        refine->phi_boundary=example.phi_boundary;}
    example.rigid_body_collection.Update_Kinematic_Particles();

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(thread_queue){
        ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T>* threaded_advection_scalar=new ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T>(thread_queue);
        example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(*threaded_advection_scalar);
        example.incompressible.Set_Custom_Advection(*threaded_advection_scalar);
        example.particle_levelset_evolution.Particle_Levelset(0).Set_Thread_Queue(thread_queue);
        example.particle_levelset_evolution.Particle_Levelset(0).levelset.thread_queue=thread_queue;
        if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<TV>* refinement=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<TV>*>(&example.projection))
            refinement->thread_queue=thread_queue;}
    else{
        example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);
        example.incompressible.Set_Custom_Advection(example.advection_scalar);}

    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

    example.particle_levelset_evolution.Particle_Levelset(0).levelset.Set_Custom_Boundary(*example.phi_boundary);
    example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution.Particle_Levelset(0).Use_Removed_Positive_Particles();
    example.particle_levelset_evolution.Particle_Levelset(0).Use_Removed_Negative_Particles();
    example.particle_levelset_evolution.Particle_Levelset(0).Store_Unique_Particle_Id();
    example.particle_levelset_evolution.use_particle_levelset=true;
    example.particle_levelset_evolution.Particle_Levelset(0).levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
    example.particle_levelset_evolution.Particle_Levelset(0).Set_Collision_Distance_Factors(.1,1);

    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.incompressible.projection.elliptic_solver->Set_Relative_Tolerance(1e-8);
    example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations=0;
    example.incompressible.projection.elliptic_solver->pcg.Show_Results();
    example.incompressible.projection.collidable_solver->Use_External_Level_Set(example.particle_levelset_evolution.Particle_Levelset(0).levelset);

    if(example.restart){
        example.Read_Output_Files(example.restart);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.mac_grid.dX.Min(),5);} // compute grid visibility (for advection later)
    else{
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.mac_grid.dX.Min(),5);
        example.Initialize_Phi();
        example.Adjust_Phi_With_Sources(time);
        example.particle_levelset_evolution.Make_Signed_Distance();
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
    if(!example.restart) example.particle_levelset_evolution.Seed_Particles(time);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    
    //add forces
    example.incompressible.gravity=-(T)9.8*TV::Axis_Vector(1-(TV::dimension==1));
    example.incompressible.Set_Body_Force(true);
    example.incompressible.projection.Use_Non_Zero_Divergence(false);
    example.incompressible.projection.elliptic_solver->Solve_Neumann_Regions(true);
    example.incompressible.projection.elliptic_solver->solve_single_cell_neumann_regions=false;
    example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
    example.incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
    example.incompressible.Use_Variable_Vorticity_Confinement(false);
    example.incompressible.Set_Surface_Tension(0);
    example.incompressible.Set_Variable_Surface_Tension(false);
    example.incompressible.Set_Viscosity(0);
    example.incompressible.Set_Variable_Viscosity(false);
    example.incompressible.projection.Set_Density(1e3);

    ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
    example.particle_levelset_evolution.Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,3,0,TV());

    example.Set_Boundary_Conditions(time); // get so CFL is correct
    if(!example.restart) Write_Output_Files(example.first_frame);
}
//#####################################################################
// Run
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Run(RANGE<TV_INT>& domain,const T dt,const T time)
{
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,3,false);
    example.incompressible.boundary->Fill_Ghost_Faces(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);
    LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
    PARTICLE_LEVELSET_UNIFORM<TV>& pls=example.particle_levelset_evolution.Particle_Levelset(0);
    if(pls.use_removed_positive_particles) for(NODE_ITERATOR<TV> iterator(example.mac_grid,domain);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
        for(int p=0;p<particles.Size();p++){
            TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.mac_grid,face_velocities_ghost,X);
            if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=-TV::Axis_Vector(2)*.3; // buoyancy
            particles.V(p)=V;}}
    if(pls.use_removed_negative_particles) for(NODE_ITERATOR<TV> iterator(example.mac_grid,domain);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
        for(int p=0;p<particles.Size();p++) particles.V(p)+=-TV::Axis_Vector(2)*dt*9.8; // ballistic
        for(int p=0;p<particles.Size();p++) particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(example.mac_grid,example.incompressible.force,particles.X(p));} // external forces
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::Time("Calculate Dt");
        example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
        T dt=example.cfl*example.incompressible.CFL(example.face_velocities);dt=min(dt,example.particle_levelset_evolution.CFL(false,false));
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=.5*(target_time-time);}

        // kinematic_update
        LOG::Time("Kinematic Evolution");
        kinematic_evolution.Set_External_Positions(example.rigid_body_collection.rigid_body_particles.frame,time);
        kinematic_evolution.Set_External_Velocities(example.rigid_body_collection.rigid_body_particles.twist,time,time);
        for(int i=0;i<example.rigid_body_collection.kinematic_rigid_bodies.m;i++){
            RIGID_BODY<TV>& rigid_body=example.rigid_body_collection.Rigid_Body(i);            
            rigid_body.Frame().t+=dt*rigid_body.Twist().linear;
            rigid_body.Frame().r=ROTATION<TV>::From_Rotation_Vector(dt*rigid_body.Twist().angular)*rigid_body.Frame().r;rigid_body.Frame().r.Normalize();}
        kinematic_evolution.Set_External_Positions(example.rigid_body_collection.rigid_body_particles.frame,time+dt);

        LOG::Time("Compute Occupied Blocks");
        T maximum_fluid_speed=example.face_velocities.Max_Abs().Max();
        T max_particle_collision_distance=example.particle_levelset_evolution.Particle_Levelset(0).max_collision_distance_factor*example.mac_grid.dX.Max();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.mac_grid.dX.Max(),10);

        LOG::Time("Adjust Phi With Objects");
        ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
        example.incompressible.boundary->Fill_Ghost_Faces(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

        example.Adjust_Phi_With_Objects(time);
        
        //Advect Phi 3.6% (Parallelized)
        LOG::Time("Advect Phi");
        example.phi_boundary_water.Use_Extrapolation_Mode(false);
        example.particle_levelset_evolution.Advance_Levelset(dt);
        example.phi_boundary_water.Use_Extrapolation_Mode(true);

        LOG::Time("Extrapolate Phi Into Objects");
        example.Extrapolate_Phi_Into_Objects(time+dt);
        
        //Advect Particles 12.1% (Parallelized)
        LOG::Time("Step Particles");
        example.particle_levelset_evolution.Particle_Levelset(0).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
        
        //Advect removed particles (Parallelized)
        LOG::Time("Advect Removed Particles");
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
        RANGE<TV_INT> domain(example.mac_grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
        DOMAIN_ITERATOR_THREADED_ALPHA<WATER_DRIVER<TV>,TV>(domain,0).template Run<T,T>(*this,&WATER_DRIVER<TV>::Run,dt,time);

        //Advect Velocities 26% (Parallelized)
        LOG::Time("Advect V");
        example.incompressible.advection->Update_Advection_Equation_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,face_velocities_ghost,*example.incompressible.boundary,dt,time);
        
        //Add Forces 0%
        LOG::Time("Forces");
        example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,dt,time,true,0,example.number_of_ghost_cells);
        kinematic_evolution.Set_External_Velocities(example.rigid_body_collection.rigid_body_particles.twist,time+dt,time+dt);

        //Modify Levelset with Particles 15% (Parallelizedish)
        LOG::Time("Modify Levelset");
        example.particle_levelset_evolution.Particle_Levelset(0).Exchange_Overlap_Particles();
        example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);
        //example.particle_levelset_evolution.Make_Signed_Distance(); //TODO(mlentine) Figure out why this was needed

        //Adjust Phi 0%
        LOG::Time("Adjust Phi");
        example.Adjust_Phi_With_Sources(time+dt);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);

        //Delete Particles 12.5 (Parallelized)
        LOG::Time("Delete Particles");
        example.particle_levelset_evolution.Delete_Particles_Outside_Grid();                                                            //0.1%
        example.particle_levelset_evolution.Particle_Levelset(0).Delete_Particles_In_Local_Maximum_Phi_Cells(1);                           //4.9%
        example.particle_levelset_evolution.Particle_Levelset(0).Delete_Particles_Far_From_Interface(); // uses visibility                 //7.6%
        example.particle_levelset_evolution.Particle_Levelset(0).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt); //2.4%
        
        //Reincorporate Particles 0% (Parallelized)
        LOG::Time("Reincorporate Particles");
        if(example.particle_levelset_evolution.Particle_Levelset(0).use_removed_positive_particles || example.particle_levelset_evolution.Particle_Levelset(0).use_removed_negative_particles)
            example.particle_levelset_evolution.Particle_Levelset(0).Reincorporate_Removed_Particles(1,1,0,true);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);

        //Project 7% (Parallelizedish)
        //LOG::Time("Project");
        {
        LOG::SCOPE *scope=0;
        if(!thread_queue) scope=new LOG::SCOPE("Project");
        LOG::Time("Boundary Conditions");
        example.Set_Boundary_Conditions(time);
        example.incompressible.Set_Dirichlet_Boundary_Conditions(&example.particle_levelset_evolution.phi,0);
        example.projection.p*=dt;
        {
        LOG::SCOPE *scope=0;
        if(!thread_queue) scope=new LOG::SCOPE("Implicit Part");
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
        delete scope;
        }
        LOG::Time("Boundary Condition Face");
        example.projection.p*=(1/dt);
        example.incompressible.boundary->Apply_Boundary_Condition_Face(example.incompressible.grid,example.face_velocities,time+dt);
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
        delete scope;
        }

        //Extrapolate Velocity 7%
        LOG::Time("Extrapolate Velocity");
        ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
        example.particle_levelset_evolution.Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time+dt,8);
        example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,3,0,TV());

        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);        
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        LOG::Time("Reseed");
        if((current_frame-example.first_frame)%1==0){
            example.particle_levelset_evolution.Reseed_Particles(time);
            example.particle_levelset_evolution.Delete_Particles_Outside_Grid();}
        Write_Output_Files(++output_number);
        current_frame++;}
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
namespace PhysBAM{
template class WATER_DRIVER<VECTOR<float,2> >;
template class WATER_DRIVER<VECTOR<float,3> >;
template class WATER_DRIVER<VECTOR<double,2> >;
template class WATER_DRIVER<VECTOR<double,3> >;
}
