//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Vectors/VECTOR_UTILITIES.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MPI.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <Dynamics/Drivers/PLS_DRIVER.h>
#include <Dynamics/Drivers/PLS_EXAMPLE.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
using namespace PhysBAM;
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_DRIVER<TV>::
PLS_DRIVER(PLS_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_DRIVER<TV>::
~PLS_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;

    if(example.auto_restart){
        example.viewer_dir.Read_Last_Frame(0);
        example.restart=example.viewer_dir.frame_stack(0);}
    else if(example.restart)
        example.viewer_dir.Set(example.restart);
    if(example.restart) example.viewer_dir.Make_Common_Directory(true);
    current_frame=example.restart;
    time=example.Time_At_Frame(current_frame);

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

    VECTOR<VECTOR<bool,2>,TV::m> domain_open_boundaries=Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.particle_levelset_evolution.Levelset_Advection(0).Set_Custom_Advection(example.advection_scalar);
    //example.incompressible.Set_Custom_Advection(example.advection_scalar);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.Particle_Levelset(0).Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    
    example.particle_levelset_evolution.Particle_Levelset(0).number_of_ghost_cells=example.number_of_ghost_cells;
    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);

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

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.Particle_Levelset(0).Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    if(example.restart){
        example.Read_Output_Files();
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
    example.incompressible.gravity=-(T)9.8*TV::Axis_Vector(1-(TV::m==1));
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

    int extrapolation_cells=2*example.number_of_ghost_cells+2;
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(extrapolation_cells));
    example.particle_levelset_evolution.Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,extrapolation_cells);
    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,example.number_of_ghost_cells,0,TV());

    example.Set_Boundary_Conditions(time); // get so CFL is correct
    if(!example.restart) Write_Output_Files();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        T dt=example.incompressible.CFL(example.face_velocities);dt=min(dt,example.particle_levelset_evolution.CFL(false,false));dt=example.cfl*dt;
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=.5*(target_time-time);}
        LOG::cout<<"dt is "<<dt<<std::endl;

        T maximum_fluid_speed=example.face_velocities.Max_Abs().Max();
        T max_particle_collision_distance=example.particle_levelset_evolution.Particle_Levelset(0).max_collision_distance_factor*example.mac_grid.dX.Max();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.mac_grid.dX.Max(),10);

        LOG::Time("Fill ghost");
PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",1);
        ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(example.incompressible.grid,example.number_of_ghost_cells,no_init);
        example.incompressible.boundary->Fill_Ghost_Faces(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

        ARRAY<T,TV_INT> phi_back(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
        example.phi_boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.Particle_Levelset(0).levelset.phi,phi_back,dt,time,example.number_of_ghost_cells);
        LOG::Time("Advect Levelset");
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi",1);
        example.phi_boundary_water.Use_Extrapolation_Mode(false);
        example.particle_levelset_evolution.Advance_Levelset(dt);
        ARRAY<T,TV_INT> phi_ghost(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
        //example.phi_boundary_water.Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,phi_ghost,dt,time,example.number_of_ghost_cells);
        //example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.particle_levelset_evolution.phi,phi_ghost,example.face_velocities,example.phi_boundary_water,dt,time);
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",1);
        LOG::Time("advecting particles");
        example.particle_levelset_evolution.Particle_Levelset(0).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",1);
        
        LOG::Time("updating removed particle velocities");
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
        PARTICLE_LEVELSET_UNIFORM<TV>& pls=example.particle_levelset_evolution.Particle_Levelset(0);
        LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
        if(pls.use_removed_positive_particles) for(NODE_ITERATOR<TV> iterator(example.mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
            for(int p=0;p<particles.Size();p++){
                TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.mac_grid,face_velocities_ghost,X);
                if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=-TV::Axis_Vector(1)*.3; // buoyancy
                particles.V(p)=V;}}
        if(pls.use_removed_negative_particles) for(NODE_ITERATOR<TV> iterator(example.mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
            for(int p=0;p<particles.Size();p++) particles.V(p)+=-TV::Axis_Vector(1)*dt*9.8; // ballistic
            for(int p=0;p<particles.Size();p++) particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(example.mac_grid,example.incompressible.force,particles.X(p));} // external forces
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",1);

        example.particle_levelset_evolution.Make_Signed_Distance();
        example.phi_boundary->Fill_Ghost_Cells(example.mac_grid,pls.levelset.phi,phi_ghost,dt,time,example.number_of_ghost_cells);
        example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
        //example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",1);
        example.Advect_Particles(dt,time);
        example.incompressible.Add_Energy_With_Vorticity(example.face_velocities,example.domain_boundary,dt,time,example.number_of_ghost_cells,&pls.levelset);
        example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,dt,time,false,0,example.number_of_ghost_cells);
        example.boundary->Fill_Ghost_Faces(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",1);

        LOG::Time("modifying levelset");
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        example.particle_levelset_evolution.Particle_Levelset(0).Exchange_Overlap_Particles();
        example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);
        example.particle_levelset_evolution.Make_Signed_Distance();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",1);

        LOG::Time("adding sources");
        example.Adjust_Phi_With_Sources(time+dt);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        LOG::Time("getting sources");

        LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
        example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",1);
        LOG::Time("deleting particles in local maxima");
        example.particle_levelset_evolution.Particle_Levelset(0).Delete_Particles_In_Local_Maximum_Phi_Cells(0);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",1);
        LOG::Time("deleting particles far from interface");
        example.particle_levelset_evolution.Particle_Levelset(0).Delete_Particles_Far_From_Interface(); // uses visibility
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",1);

        LOG::Time("re-incorporating removed particles");
        example.particle_levelset_evolution.Particle_Levelset(0).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",1);
        if(example.particle_levelset_evolution.Particle_Levelset(0).use_removed_positive_particles || example.particle_levelset_evolution.Particle_Levelset(0).use_removed_negative_particles)
            example.particle_levelset_evolution.Particle_Levelset(0).Reincorporate_Removed_Particles(1,1,0,true);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after reincorporate",1);

        // update ghost phi values
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("before boundary",1);
        LOG::Time("getting Neumann and Dirichlet boundary conditions");
        example.Set_Boundary_Conditions(time);
        example.incompressible.Set_Dirichlet_Boundary_Conditions(&example.particle_levelset_evolution.phi,0);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",1);

        example.projection.p*=dt;
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
        example.projection.p*=(1/dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after solve",1);

        example.incompressible.boundary->Apply_Boundary_Condition_Face(example.incompressible.grid,example.face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",1);
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);

        LOG::Time("extrapolating velocity across interface");
        int band_width=example.number_of_ghost_cells+1;
        ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(2*band_width+2));
        example.particle_levelset_evolution.Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time+dt,2*band_width+2);
        example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,band_width,0,TV());
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after extrapolate",1);

        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        example.particle_levelset_evolution.Reseed_Particles(time);
        example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
        Write_Output_Files();
        current_frame++;}
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    example.frame_title=title;
    LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<example.viewer_dir.frame_stack<<", time="<<time<<", frame="<<current_frame<<std::endl;
    Write_Output_Files();
    example.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PLS_DRIVER<TV>::
Write_Output_Files()
{
    example.viewer_dir.Start_Directory(0,example.frame_title);
    example.frame_title="";
    example.Write_Output_Files();
    example.viewer_dir.Finish_Directory();
}
//#####################################################################
namespace PhysBAM{
template class PLS_DRIVER<VECTOR<float,2> >;
template class PLS_DRIVER<VECTOR<float,3> >;
template class PLS_DRIVER<VECTOR<double,2> >;
template class PLS_DRIVER<VECTOR<double,3> >;
}
