//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLS_FSI_DRIVER
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Utilities/INTERRUPTS.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
#include <Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <Dynamics/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <Dynamics/Kang/KANG_POISSON_VISCOSITY.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
//#define DISABLE_LEVELSET_ADVECTION
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PLS_FSI_DRIVER<TV>::
PLS_FSI_DRIVER(PLS_FSI_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FSI_DRIVER<TV>::
~PLS_FSI_DRIVER()
{}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Preprocess_Frame(const int frame)
{
    if(example.substeps_delay_frame==frame){
        example.Set_Write_Substeps_Level(example.substeps_delay_level);
        output_number=frame-1;}
    example.Preprocess_Frame(frame);
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files();}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Preprocess_Frame(current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files();
        else if(example.write_substeps_level!=-1)
            PHYSBAM_DEBUG_WRITE_SUBSTEP("END Frame %d",example.write_substeps_level,current_frame);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Initialize()
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*example.fluids_parameters.collision_bodies_affecting_fluid;

    if(example.auto_restart){
        example.viewer_dir.Read_Last_Frame(0);
        example.restart=true;
        example.restart_frame=example.viewer_dir.frame_stack(0);}
    if(example.restart) current_frame=example.restart_frame;
    else current_frame=0;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    example.fluids_parameters.callbacks=&example;

    *example.fluids_parameters.grid=example.fluids_parameters.grid->Get_MAC_Grid();
    example.fluids_parameters.p_grid=*example.fluids_parameters.grid;
    example.Initialize_Fluids_Grids();
    GRID<TV>& grid=*example.fluids_parameters.grid;

    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=0;
    if(example.use_kang){
        example.kang_poisson_viscosity=new KANG_POISSON_VISCOSITY<TV>(example.fluids_parameters,old_phi);
        example.kang_poisson_viscosity->print_matrix=example.print_matrix;
        example.kang_poisson_viscosity->test_system=example.test_system;}
    else{
        example.fluids_parameters.use_poisson=true;
        coupled_evolution=new SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>(example.solids_parameters,example.solid_body_collection,example,
            example.fluids_parameters,example.solids_fluids_parameters,example.fluid_collection);
        delete example.solids_evolution;
        example.solids_evolution=coupled_evolution;
        example.fluids_parameters.projection=coupled_evolution;
        SOLIDS_EVOLUTION<TV>& solids_evolution=*coupled_evolution;
        solids_evolution.Set_Solids_Evolution_Callbacks(example);}
    
    example.Initialize_Bodies();

    if(example.use_kang)
        example.Set_Boundary_Conditions(example.kang_poisson_viscosity->psi_D,example.kang_poisson_viscosity->psi_N,
            example.kang_poisson_viscosity->psi_D_value,example.kang_poisson_viscosity->psi_N_value);

    example.fluids_parameters.particle_levelset_evolution=new PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>(*example.fluids_parameters.grid,
        collision_bodies_affecting_fluid,example.fluids_parameters.number_of_ghost_cells,false);
    example.fluids_parameters.projection=new PROJECTION_DYNAMICS_UNIFORM<TV>(*example.fluids_parameters.grid,example.fluids_parameters.particle_levelset_evolution->Levelset(0));
    example.fluids_parameters.incompressible=new INCOMPRESSIBLE_UNIFORM<TV>(*example.fluids_parameters.grid,*example.fluids_parameters.projection);
    example.fluids_parameters.phi_boundary=&example.fluids_parameters.phi_boundary_water; // override default
    example.fluids_parameters.phi_boundary_water.Set_Velocity_Pointer(example.fluid_collection.incompressible_fluid_collection.face_velocities);
    example.fluids_parameters.boundary_mac_slip.Set_Phi(example.fluids_parameters.particle_levelset_evolution->phi);

    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=example.fluids_parameters.incompressible;

    Initialize_Fluids_Grids(); // this needs to be here because the arrays have to be resized for multiphase

    example.fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();
    if(coupled_evolution) coupled_evolution->Setup_Boundary_Condition_Collection();
    example.solids_evolution->time=time;

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids();
        example.solids_evolution->time=time=example.Time_At_Frame(example.restart_frame);}

    example.solids_evolution->Initialize_Deformable_Objects(example.frame_rate,example.restart);

    example.solids_evolution->Initialize_Rigid_Bodies(example.frame_rate,example.restart);

    // time
    particle_levelset_evolution->time=time;
    particle_levelset_evolution->Set_CFL_Number(example.fluids_parameters.cfl);

    // sets up the proper wall states
    VECTOR<VECTOR<bool,2>,TV::m> domain_open_boundaries=Complement(example.fluids_parameters.domain_walls);
    if(example.fluids_parameters.phi_boundary) example.fluids_parameters.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(example.fluids_parameters.fluid_boundary) example.fluids_parameters.fluid_boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    example.Initialize_Advection();
    Initialize_Fluids_Grids(); // initialize the valid masks

    // initialize levelset
    particle_levelset_evolution->Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
    particle_levelset_evolution->Set_Levelset_Callbacks(example);

    particle_levelset_evolution->Particle_Levelset(0).levelset.Set_Custom_Boundary(*example.fluids_parameters.phi_boundary);
    particle_levelset_evolution->Bias_Towards_Negative_Particles(example.fluids_parameters.bias_towards_negative_particles);
    if(example.fluids_parameters.use_removed_positive_particles) particle_levelset_evolution->Particle_Levelset(0).Use_Removed_Positive_Particles();
    if(example.fluids_parameters.use_removed_negative_particles) particle_levelset_evolution->Particle_Levelset(0).Use_Removed_Negative_Particles();
    if(example.fluids_parameters.store_particle_ids){
        particle_levelset_evolution->Particle_Levelset(0).Store_Unique_Particle_Id();}
    particle_levelset_evolution->use_particle_levelset=example.fluids_parameters.use_particle_levelset;

    // solid fluid coupling
    particle_levelset_evolution->Particle_Levelset(0).levelset.Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);
    particle_levelset_evolution->Particle_Levelset(0).Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,
        example.fluids_parameters.max_collision_distance_factor);

    // incompressible flow
    incompressible->Set_Custom_Boundary(*example.fluids_parameters.fluid_boundary);

    // set up the initial state
    if(example.restart){
        example.Read_Output_Files_Fluids();
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for advection later)
        particle_levelset_evolution->Set_Seed(2606);
        if(!example.fluids_parameters.write_particles) particle_levelset_evolution->Seed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);}
    else{
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);
        example.Initialize_Phi();
        example.Adjust_Phi_With_Sources(time);
        particle_levelset_evolution->Make_Signed_Distance();
        particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for averaging face velocities to nodes below)
        particle_levelset_evolution->Set_Seed(2606);
        particle_levelset_evolution->Seed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);

        example.Initialize_Velocities();
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);

        int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
        T extrapolation_bandwidth=(T)(extrapolation_ghost_cells-1);
        ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
        particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time,extrapolation_ghost_cells);
//        Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,particle_levelset_evolution->Particle_Levelset(0).levelset,extrapolation_bandwidth);
        if(!example.two_phase) 
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,exchanged_phi_ghost,
                example.fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
    }
}
//#####################################################################
// Function Initialize_Fluids_Grids
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Initialize_Fluids_Grids()
{
    GRID<TV>& grid=*example.fluids_parameters.grid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    *example.fluids_parameters.grid=example.fluids_parameters.grid->Get_MAC_Grid();
    example.fluids_parameters.p_grid=*example.fluids_parameters.grid;
    example.Initialize_Fluids_Grids();
    example.fluid_collection.Initialize_Grids();
    fluids_parameters.particle_levelset_evolution->Initialize_Domain(fluids_parameters.p_grid);
    fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).Set_Band_Width((T)2*fluids_parameters.particle_half_bandwidth);
    fluids_parameters.incompressible->valid_mask.Resize(grid.Domain_Indices(example.fluids_parameters.number_of_ghost_cells),use_init,true);
    fluids_parameters.incompressible->grid=grid.Get_MAC_Grid();
//    p.Resize(p_grid.Domain_Indices(1));p_save_for_projection.Resize(p_grid.Domain_Indices(1));face_velocities_save_for_projection.Resize(p_grid);
    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* slip=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(&solids_evolution))
        slip->Initialize_Grid_Arrays();
}
//#####################################################################
// Function First_Order_Time_Step
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
First_Order_Time_Step(int substep,T dt)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* slip=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(example.solids_evolution);
    GRID<TV>& grid=*fluids_parameters.grid;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;

    if(example.use_kang) old_phi=fluids_parameters.particle_levelset_evolution->Levelset(0).phi;

    example.solid_body_collection.Print_Energy(time,1);
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);  // static occupied blocks
    // swept occupied blocks
    example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,example.fluid_collection.incompressible_fluid_collection.face_velocities);

    if(example.use_pls_evolution_for_structure) Advance_Particles_With_PLS(dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("start step",1);
    if(example.kang_poisson_viscosity && !fluids_parameters.implicit_viscosity){
        face_velocities_scratch=face_velocities;
        PHYSBAM_DEBUG_WRITE_SUBSTEP("explicit viscosity",1);
        example.kang_poisson_viscosity->Apply_Viscosity(face_velocities,dt,fluids_parameters.implicit_viscosity);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("explicit viscosity",1);
        face_velocities_scratch.Exchange(face_velocities,face_velocities_scratch);
        face_velocities_scratch-=face_velocities;
        PHYSBAM_DEBUG_WRITE_SUBSTEP("revert velocity and store difference",1);}

    Advect_Fluid(dt);
    LOG::cout<<"Maximum face velocity (after advect) = ("<<face_velocities.Max_Abs().Magnitude()<<": "<<face_velocities.Max_Abs()<<std::endl;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("advect",1);
    if(example.kang_poisson_viscosity && !fluids_parameters.implicit_viscosity){
        face_velocities+=face_velocities_scratch;
        PHYSBAM_DEBUG_WRITE_SUBSTEP("apply stored explicit viscosity",1);}
    

    example.solids_evolution->example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);
    example.solid_body_collection.Update_Position_Based_State(time+dt,true,true);
//    example.solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<VECTOR<T,2> >*>()->Dump_Curvatures();
    if(slip){
        slip->two_phase=example.two_phase;
        slip->Solve(face_velocities,dt,time,time+dt,false,false);}
    else if(example.kang_poisson_viscosity)
        example.kang_poisson_viscosity->Project_Fluid(face_velocities,dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("pressure solve",1);

    if(slip) slip->Print_Maximum_Velocities(time);

//        if(incompressible) incompressible->projection.p_save_for_projection.Copy(incompressible->projection.p); // save the good pressure for later
    if(!example.use_pls_evolution_for_structure && slip) slip->Euler_Step_Position(dt,time+dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("euler step position",1);

    LOG::Time("extrapolating velocity across interface");
//    Extrapolate_Velocity_Across_Interface(example.face_velocities,particle_levelset_evolution->Particle_Levelset(0).levelset,extrapolation_bandwidth);
    if(!example.two_phase) Extrapolate_Velocity_Across_Interface(time,dt);
        //incompressible->Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,
        //    fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("extrapolate about interface",1);
    incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,face_velocities,time+dt);
    LOG::cout<<"Maximum face velocity = ("<<face_velocities.Max_Abs().Magnitude()<<": "<<face_velocities.Max_Abs()<<std::endl;

    PHYSBAM_DEBUG_WRITE_SUBSTEP("end step",1);
    example.solids_evolution->time+=dt;
    time+=dt;
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Extrapolate_Velocity_Across_Interface(T time,T dt)
{
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    GRID<TV>& grid=*example.fluids_parameters.grid;
    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& slip=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*example.solids_evolution);
    ARRAY<bool,FACE_INDEX<TV::m> >& valid_faces=slip.solved_faces;
    int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
    particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,phi_ghost,0,time+dt,extrapolation_ghost_cells);
    int band_width=extrapolation_ghost_cells-1;
    T delta=band_width*grid.dX.Max();
    for(int axis=0;axis<TV::m;axis++){
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        ARRAY<T,TV_INT> phi_face(face_grid.Domain_Indices(),no_init);
        ARRAYS_ND_BASE<T,TV_INT>& face_velocity=example.fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis);
        ARRAYS_ND_BASE<bool,TV_INT>& fixed_face=valid_faces.Component(axis);
        for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Face_Index();
            T phi1=phi_ghost(iterator.First_Cell_Index()),phi2=phi_ghost(iterator.Second_Cell_Index());
            phi_face(index)=(T).5*(phi1+phi2);
            if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}

        EXTRAPOLATION_UNIFORM<TV,T> extrapolate(face_grid,phi_face,face_velocity,extrapolation_ghost_cells);
        extrapolate.Set_Band_Width((T)band_width);
        extrapolate.Set_Custom_Seed_Done(&fixed_face);
        extrapolate.Extrapolate();}
}
//#####################################################################
// Function Advance_Particles_With_PLS
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advance_Particles_With_PLS(T dt)
{
    ARRAY_VIEW<TV> X=example.solid_body_collection.deformable_body_collection.particles.X;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=example.fluid_collection.incompressible_fluid_collection.face_velocities;
    LINEAR_INTERPOLATION_MAC<TV,T> interp(*example.fluids_parameters.grid);
    for(RUNGEKUTTA<ARRAY_VIEW<TV> > rk(X,example.fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles,dt,0);rk.Valid();rk.Next())
        for(int p=0;p<X.m;p++)
            X(p)+=dt*interp.Clamped_To_Array(face_velocities,X(p));
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=example.solids_evolution->solids_evolution_callbacks;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;

    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        solids_evolution_callbacks->Preprocess_Solids_Substep(time);
        particle_levelset_evolution->Set_Number_Particles_Per_Cell(fluids_parameters.number_particles_per_cell);
        T dt=Compute_Dt(time,target_time,done);
        example.Preprocess_Substep(dt,time);
        example.Update_Fluid_Parameters(dt,time);

        First_Order_Time_Step(substep,dt);

        solids_evolution_callbacks->Postprocess_Solids_Substep(example.solids_evolution->time);
        example.Postprocess_Substep(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("END Substep %d",0,substep);}
}
//#####################################################################
// Function Advect_Fluid
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advect_Fluid(const T dt)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    GRID<TV>& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;

    LOG::SCOPE scalar_scope("SCALAR SCOPE");

    // important to compute ghost velocity values for particles near domain boundaries
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    face_velocities_ghost.Resize(incompressible->grid,fluids_parameters.number_of_ghost_cells,no_init);
    incompressible->boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time+dt,fluids_parameters.number_of_ghost_cells);

    fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Adjust_Phi_With_Objects(time);
    LOG::Time("advecting levelset");
#ifndef DISABLE_LEVELSET_ADVECTION
    particle_levelset_evolution->Advance_Levelset(dt);
#endif // #ifndef DISABLE_LEVELSET_ADVECTION
    fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);
    example.Extrapolate_Phi_Into_Objects(time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after levelset advection",1);
    LOG::Time("advecting particles");
#ifndef DISABLE_LEVELSET_ADVECTION
    particle_levelset_evolution->Advance_Particles(face_velocities_ghost,dt,false);
#endif // #ifndef DISABLE_LEVELSET_ADVECTION
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle advection",1);

    example.Scalar_Advection_Callback(dt,time);

    LOG::Time("updating removed particle velocities");
    example.Modify_Removed_Particles_Before_Advection(dt,time);
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
    PARTICLE_LEVELSET_UNIFORM<TV>& pls=particle_levelset_evolution->Particle_Levelset(0);
    LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
    if(pls.use_removed_positive_particles) for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
        for(int p=0;p<particles.Size();p++){
            TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(grid,face_velocities_ghost,X);
            if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=fluids_parameters.removed_positive_particle_buoyancy_constant*fluids_parameters.gravity.Normalized(); // buoyancy
            particles.V(p)=V;}}
    if(pls.use_removed_negative_particles) for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
        for(int p=0;p<particles.Size();p++) particles.V(p)+=dt*fluids_parameters.gravity; // ballistic
        if(fluids_parameters.use_body_force) for(int p=0;p<particles.Size();p++)
            particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(grid,incompressible->force,particles.X(p));} // external forces

    LOG::Time("updating velocity (explicit part)");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",1);

    int extrapolation_ghost_cells=2*fluids_parameters.number_of_ghost_cells+2;
    T extrapolation_bandwidth=(T)(extrapolation_ghost_cells-1);
    ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
    particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time+dt,extrapolation_ghost_cells);
    if(example.convection_order>1){
        for(RUNGEKUTTA<ARRAY<T,FACE_INDEX<TV::m> > > rk(face_velocities,example.convection_order,dt,time);rk.Valid();rk.Next()){
//            Extrapolate_Velocity_Across_Interface(face_velocities,particle_levelset_evolution->Particle_Levelset(0).levelset,extrapolation_bandwidth);
            if(!example.two_phase)
                incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,exchanged_phi_ghost,
                    fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
            incompressible->Advance_One_Time_Step_Convection(dt,rk.time,face_velocities,face_velocities,fluids_parameters.number_of_ghost_cells);}}
    else incompressible->Advance_One_Time_Step_Convection(dt,time,face_velocities_ghost,face_velocities,fluids_parameters.number_of_ghost_cells);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",1);

    if(fluids_parameters.use_body_force) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.incompressible->force,dt,time);
    incompressible->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution->phi,example.fluids_parameters.number_of_ghost_cells);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after integrate non advection forces",1);

    LOG::Time("updating velocity (explicit part without convection)");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before viscosity",1);
    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(example.solids_evolution))
        coupled_evolution->Apply_Viscosity(face_velocities,dt,time);
    else if(example.kang_poisson_viscosity && fluids_parameters.implicit_viscosity)
        example.kang_poisson_viscosity->Apply_Viscosity(face_velocities,dt,fluids_parameters.implicit_viscosity);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after viscosity",1);

    LOG::Time("effective velocity acceleration structures");
    // revalidate scalars and velocity in body's new position
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before scalar revalidation",1);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);  // static occupied blocks
    collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
    example.Revalidate_Fluid_Scalars(); // uses visibility

    example.Revalidate_Fluid_Velocity(fluid_collection.incompressible_fluid_collection.face_velocities); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after scalar revalidation",1);

    LOG::Time("modifying levelset");
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after filling level set ghost cells and exchanging overlap particles",1);
    particle_levelset_evolution->Modify_Levelset_And_Particles(&face_velocities_ghost);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after modify levelset and particles",1);
    example.Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after revalidate phi",1);

    LOG::Time("adding sources");
    if(example.Adjust_Phi_With_Sources(time+dt)) particle_levelset_evolution->Make_Signed_Distance();
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
    LOG::Time("getting sources");
    ARRAY<bool,TV_INT>* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time+dt);
    if(source_mask){LOG::Time("reseeding sources");particle_levelset_evolution->Reseed_Particles(time+dt,0,source_mask);delete source_mask;}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after adding sources",1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    particle_levelset_evolution->Delete_Particles_Outside_Grid();
    if(fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time+dt);
    LOG::Time("deleting particles in local maxima");
    particle_levelset_evolution->Particle_Levelset(0).Delete_Particles_In_Local_Maximum_Phi_Cells(1);
    LOG::Time("deleting particles far from interface");
    particle_levelset_evolution->Particle_Levelset(0).Delete_Particles_Far_From_Interface(); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete particles far from interface",1);

    LOG::Time("re-incorporating removed particles");
    // TODO: if your particles fall entirely within the grid this shouldn't need ghost cells, but it may need them for MPI
    particle_levelset_evolution->Particle_Levelset(0).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
    if(particle_levelset_evolution->Particle_Levelset(0).use_removed_positive_particles || particle_levelset_evolution->Particle_Levelset(0).use_removed_negative_particles)
        particle_levelset_evolution->Particle_Levelset(0).Reincorporate_Removed_Particles(1,fluids_parameters.removed_particle_mass_scaling,
            fluids_parameters.reincorporate_removed_particle_velocity?&example.fluid_collection.incompressible_fluid_collection.face_velocities:0,true);

    // update ghost phi values
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after scalar update",1);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Postprocess_Frame(const int frame)
{
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;

    example.Postprocess_Phi(time);

    if(particle_levelset_evolution->use_particle_levelset && (frame)%example.fluids_parameters.reseeding_frame_rate==0){
        LOG::Time("Reseeding... ");
        particle_levelset_evolution->Reseed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);
        LOG::Stop_Time();}

    example.Postprocess_Frame(frame);
    example.solids_evolution->Postprocess_Frame(frame);
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR PLS_FSI_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    T fluids_dt=FLT_MAX;
    {T dt_levelset=fluids_parameters.particle_levelset_evolution->CFL(true,false);fluids_dt=min(fluids_dt,dt_levelset);}
    LOG::cout<<"before example cfl clamping:  "<<fluids_dt<<std::endl;

    // solids dt
    T solids_dt=FLT_MAX;
    if(solids_evolution.Use_CFL()) solids_dt=min(solids_dt,example.solid_body_collection.CFL(solids_parameters.verbose_dt));
    solids_dt=min(solids_dt,solids_evolution_callbacks->Constraints_CFL());
    solids_dt=min(solids_dt,example.solid_body_collection.rigid_body_collection.CFL_Rigid(solids_parameters.rigid_body_evolution_parameters,solids_parameters.verbose_dt));

    T dt=min(fluids_dt,solids_dt);
    done=example.Clamp_Time_Step_With_Target_Time(time,target_time,dt);
    LOG::cout<<"fluids_dt = "<<fluids_dt<<", solids_dt = "<<solids_dt<<" dt="<<dt<<std::endl;
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Write_Output_Files()
{
    LOG::SCOPE scope("writing output files");
    example.viewer_dir.Start_Directory(0,example.frame_title);
    example.frame_title="";

    example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Write_Output_Files();
    example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);

    GRID<TV>& grid=*example.fluids_parameters.grid;
    int number_of_positive_particles=0,number_of_negative_particles=0,number_of_removed_positive_particles=0,number_of_removed_negative_particles=0;
    PARTICLE_LEVELSET_UNIFORM<TV>* pls=0;
    pls=&example.fluids_parameters.particle_levelset_evolution->Particle_Levelset(0);
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls->positive_particles(iterator.Cell_Index()))
        number_of_positive_particles+=pls->positive_particles(iterator.Cell_Index())->Size();
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls->negative_particles(iterator.Cell_Index()))
        number_of_negative_particles+=pls->negative_particles(iterator.Cell_Index())->Size();
    LOG::cout<<number_of_positive_particles<<" positive and "<<number_of_negative_particles<<" negative particles "<<std::endl;
    if(pls->use_removed_positive_particles)
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls->removed_positive_particles(iterator.Cell_Index()))
            number_of_removed_positive_particles+=pls->removed_positive_particles(iterator.Cell_Index())->Size();
    if(pls->use_removed_negative_particles)
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls->removed_negative_particles(iterator.Cell_Index()))
            number_of_removed_negative_particles+=pls->removed_negative_particles(iterator.Cell_Index())->Size();
    LOG::cout<<number_of_removed_positive_particles<<" positive and "<<number_of_removed_negative_particles<<" negative removed particles "<<std::endl;

    Write_Time();
    example.viewer_dir.Finish_Directory();
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Delete_Particles_Inside_Objects(const T time)
{
    PARTICLE_LEVELSET_UNIFORM<TV>* particle_levelset=&example.fluids_parameters.particle_levelset_evolution->Particle_Levelset(0);
    Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->positive_particles,PARTICLE_LEVELSET_POSITIVE,time);
    Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->negative_particles,PARTICLE_LEVELSET_NEGATIVE,time);
    if(particle_levelset->use_removed_positive_particles)
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,time);
    if(particle_levelset->use_removed_negative_particles)
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> template<class T_PARTICLES> void PLS_FSI_DRIVER<TV>::
Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    for(NODE_ITERATOR<TV> iterator(*example.fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        BLOCK_UNIFORM<TV> block(*example.fluids_parameters.grid,block_index);
        COLLISION_GEOMETRY_ID body_id;int aggregate_id;
        T_PARTICLES& block_particles=*particles(block_index);
        if(example.fluids_parameters.collision_bodies_affecting_fluid->Occupied_Block(block)){
            for(int k=block_particles.Size()-1;k>=0;k--)
                if(example.fluids_parameters.collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(block_particles.X(k),body_id,aggregate_id))
                    block_particles.Delete_Element(k);}
        example.fluids_parameters.callbacks->Delete_Particles_Inside_Objects(block_particles,particle_type,time);}}
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const LEVELSET<TV>& phi,const T band_width)
{
    GRID<TV>& grid=*example.fluids_parameters.grid;
    EXTRAPOLATION_HIGHER_ORDER<TV,T> eho(grid,phi,20,example.fluids_parameters.number_of_ghost_cells,(int)ceil(band_width));
    eho.Extrapolate_Face([&](const FACE_INDEX<TV::m>& index){return phi.Phi(grid.Face(index))<=0;},face_velocities);
}
//#####################################################################
namespace PhysBAM{
template class PLS_FSI_DRIVER<VECTOR<float,1> >;
template class PLS_FSI_DRIVER<VECTOR<float,2> >;
template class PLS_FSI_DRIVER<VECTOR<float,3> >;
template class PLS_FSI_DRIVER<VECTOR<double,1> >;
template class PLS_FSI_DRIVER<VECTOR<double,2> >;
template class PLS_FSI_DRIVER<VECTOR<double,3> >;
}
