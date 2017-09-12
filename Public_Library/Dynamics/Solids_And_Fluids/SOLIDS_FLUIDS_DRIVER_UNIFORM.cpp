//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Levine, Nick Rasmussen, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_UNIFORM
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Utilities/INTERRUPTS.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <Dynamics/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
SOLIDS_FLUIDS_DRIVER_UNIFORM(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>& example)
    :SOLIDS_FLUIDS_DRIVER<TV>(example),example(example),last_dt(0),restart_dt(0),reset_with_restart(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
~SOLIDS_FLUIDS_DRIVER_UNIFORM()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Initialize()
{
    SOLIDS_FLUIDS_DRIVER<TV>::Initialize();
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    example.fluids_parameters.Set_Fluids_Parameters_Callbacks(example);

    bool fluids=Simulate_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());
    if(fluids) example.fluids_parameters.Initialize_Grids(); // this needs to be here because Initialize_Fluid_Evolution and Initialize_Solid_Fluid_Coupling depend on it

    if(Two_Way_Coupled() || (example.fluids_parameters.solid_affects_fluid && example.fluids_parameters.use_slip)) example.Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();

    if(fluids) example.fluids_parameters.Initialize_Fluid_Evolution(example.fluid_collection.incompressible_fluid_collection.face_velocities);
    example.After_Initialization();

    int number_of_regions=example.fluids_parameters.number_of_regions;
    GRID<TV>& grid=*example.fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=example.fluids_parameters.incompressible;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=example.fluids_parameters.particle_levelset_evolution_multiple;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase=example.fluids_parameters.incompressible_multiphase;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*example.fluids_parameters.collision_bodies_affecting_fluid;
    EULER_UNIFORM<TV>* euler=example.fluids_parameters.euler;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;

    // this needs to be before Initialize_Fluids_Grids
    if(example.fluids_parameters.use_sph_for_removed_negative_particles){
        for(int i=0;i<number_of_regions;i++){
            PARTICLE_LEVELSET_UNIFORM<TV>& pls=example.fluids_parameters.particle_levelset_evolution->Particle_Levelset(i);
            pls.template_particles.template Add_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);
            pls.template_removed_particles.template Add_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);}
        example.fluids_parameters.sph_evolution=new SPH_EVOLUTION_UNIFORM<TV>(grid,*incompressible,example.fluids_parameters,particle_levelset_evolution);
        example.fluids_parameters.sph_evolution->Set_SPH_Callbacks(example);
        example.fluids_parameters.fluid_boundary=&example.fluids_parameters.fluid_boundary_water;}

    if(fluids) Initialize_Fluids_Grids(); // this needs to be here because the arrays have to be resized for multiphase

    // initialize MPI
    if(fluids && example.fluids_parameters.compressible) example.Initialize_Advection();
    if(example.fluids_parameters.mpi_grid) example.Initialize_MPI();

    //compressible flow
    if(fluids && example.fluids_parameters.compressible){
        if(example.fluids_parameters.compressible_eos) euler->Set_Custom_Equation_Of_State(*example.fluids_parameters.compressible_eos);
        if(example.fluids_parameters.compressible_conservation_method) euler->Set_Custom_Conservation(*example.fluids_parameters.compressible_conservation_method);
        if(example.fluids_parameters.compressible_boundary) euler->Set_Custom_Boundary(*example.fluids_parameters.compressible_boundary);
        if(example.fluids_parameters.compressible_pressure_boundary) euler->euler_projection.Set_Custom_Pressure_Boundary(*example.fluids_parameters.compressible_pressure_boundary);
        if(example.fluids_parameters.compressible_set_max_time_step) euler->Set_Max_Time_Step(example.fluids_parameters.compressible_max_time_step);
        euler->conservation->Set_Order(example.fluids_parameters.compressible_spatial_order);
        example.fluids_parameters.euler->timesplit=example.fluids_parameters.compressible_timesplit;
        example.fluids_parameters.euler_solid_fluid_coupling_utilities->fluid_affects_solid=example.fluids_parameters.fluid_affects_solid;
        example.fluids_parameters.euler->perform_rungekutta_for_implicit_part=example.fluids_parameters.compressible_perform_rungekutta_for_implicit_part;
        example.fluids_parameters.euler->use_sound_speed_for_cfl=example.fluids_parameters.compressible_use_sound_speed_for_cfl;
        example.fluids_parameters.euler->use_sound_speed_based_dt_multiple_for_cfl=example.fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl;
        example.fluids_parameters.euler->multiplication_factor_for_sound_speed_based_dt=example.fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt;
        example.fluids_parameters.euler->apply_cavitation_correction=example.fluids_parameters.compressible_apply_cavitation_correction;
        euler->euler_projection.elliptic_solver->solve_single_cell_neumann_regions=example.fluids_parameters.solve_single_cell_neumann_regions;
        euler->euler_projection.elliptic_solver->Set_Relative_Tolerance(example.fluids_parameters.compressible_tolerance);
        euler->euler_projection.elliptic_solver->pcg.Set_Maximum_Iterations(example.fluids_parameters.compressible_iterations);
        euler->euler_projection.elliptic_solver->pcg.Show_Results();
        euler->euler_projection.elliptic_solver->pcg.Show_Residuals();}

    if(fluids){
        example.Initialize_Solid_Fluid_Coupling_After_Grid_Initialization();
        example.Initialize_Compressible_Incompressible_Coupling();}
    solids_evolution.time=time;

    if(!solids_parameters.fracture && !example.use_melting) // fracture and melting initialize collisions in Initialize_Bodies
        example.solid_body_collection.deformable_body_collection.Initialize(solids_parameters.triangle_collision_parameters);

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);

    // TODO: this is where a return used to be, so further stuff should be modified to be safe in solids only case
    // NOTHING MORE TO DO FOR SOLIDS SIMS
    if(!fluids) return;

    // fire smoke and water stuff...
    if(example.fluids_parameters.fire) assert(example.fluids_parameters.use_density && example.fluids_parameters.use_temperature); // both assumed to be defined
    // TODO
    if(Two_Way_Coupled()){
//        assert(!example.fluids_parameters.surface_tension && !example.fluids_parameters.variable_surface_tension);
        assert(example.fluids_parameters.use_slip || (!example.fluids_parameters.viscosity && !example.fluids_parameters.variable_viscosity));}

    // time
    if(number_of_regions>=1){
        particle_levelset_evolution->time=time;
        particle_levelset_evolution->Set_CFL_Number(example.fluids_parameters.cfl);}
    if(example.fluids_parameters.compressible) euler->Set_CFL_Number(example.fluids_parameters.cfl);

    // sets up the proper wall states
    example.fluids_parameters.Initialize_Domain_Boundary_Conditions();

    if(!example.fluids_parameters.compressible) example.Initialize_Advection();
    Initialize_Fluids_Grids(); // initialize the valid masks

    // initialize levelset
    if(number_of_regions>=1){
        particle_levelset_evolution->Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
        particle_levelset_evolution->Set_Levelset_Callbacks(example);

        if(number_of_regions>=2) for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Levelset(i).Set_Custom_Boundary(*example.fluids_parameters.phi_boundary_multiphase(i));
        else particle_levelset_evolution->Particle_Levelset(0).levelset.Set_Custom_Boundary(*example.fluids_parameters.phi_boundary);
        particle_levelset_evolution->Bias_Towards_Negative_Particles(example.fluids_parameters.bias_towards_negative_particles);
        for(int i=0;i<number_of_regions;i++){
            if(example.fluids_parameters.use_removed_positive_particles) particle_levelset_evolution->Particle_Levelset(i).Use_Removed_Positive_Particles();
            if(example.fluids_parameters.use_removed_negative_particles) particle_levelset_evolution->Particle_Levelset(i).Use_Removed_Negative_Particles();}
        if(example.fluids_parameters.store_particle_ids){
            for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Particle_Levelset(i).Store_Unique_Particle_Id();}
        particle_levelset_evolution->use_particle_levelset=example.fluids_parameters.use_particle_levelset;}

    // initialize sph
    if(example.fluids_parameters.sph) example.fluids_parameters.sph_evolution->Set_SPH_Callbacks(example);
    // solid fluid coupling
    if(number_of_regions>=2){
        for(int i=0;i<number_of_regions;i++)
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelsets(i)->Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);
        particle_levelset_evolution_multiple->particle_levelset_multiple.Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,
            example.fluids_parameters.max_collision_distance_factor);}
    else if(number_of_regions==1){
        particle_levelset_evolution->Particle_Levelset(0).levelset.Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);
        particle_levelset_evolution->Particle_Levelset(0).Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,
            example.fluids_parameters.max_collision_distance_factor);}

    PHYSBAM_ASSERT(number_of_regions<2 || (!example.fluids_parameters.use_dsd && !example.fluids_parameters.compressible));
    // incompressible flow
    if(incompressible){
        if(number_of_regions>=2){
            incompressible_multiphase->projection.poisson_collidable->Use_Internal_Level_Set(number_of_regions);
            for(int i=0;i<number_of_regions;i++)
                incompressible_multiphase->projection.poisson_collidable->levelset_multiple->levelsets(i)->Set_Custom_Boundary(*example.fluids_parameters.phi_boundary_multiphase(i));}
        incompressible->Set_Custom_Boundary(*example.fluids_parameters.fluid_boundary);
        incompressible->projection.elliptic_solver->Set_Relative_Tolerance(example.fluids_parameters.incompressible_tolerance);
        incompressible->projection.elliptic_solver->pcg.Set_Maximum_Iterations(example.fluids_parameters.incompressible_iterations);
        incompressible->projection.elliptic_solver->pcg.evolution_solver_type=example.fluids_parameters.evolution_solver_type;
        incompressible->projection.elliptic_solver->pcg.cg_restart_iterations=example.fluids_parameters.cg_restart_iterations;
        incompressible->projection.elliptic_solver->pcg.Show_Results();
        if(!example.fluids_parameters.use_reacting_flow && example.fluids_parameters.use_dsd)
            incompressible->projection.Initialize_Dsd(particle_levelset_evolution->Particle_Levelset(0).levelset,example.fluids_parameters.fuel_region);
        if(number_of_regions==1 && example.fluids_parameters.second_order_cut_cell_method)
            incompressible->projection.collidable_solver->Use_External_Level_Set(particle_levelset_evolution->Particle_Levelset(0).levelset);
        if(example.fluids_parameters.use_reacting_flow){
            incompressible_multiphase->projection.poisson_collidable->Set_Jump_Multiphase();
            if(example.fluids_parameters.use_flame_speed_multiplier) incompressible_multiphase->projection.Use_Flame_Speed_Multiplier();
            if(example.fluids_parameters.use_dsd)
                incompressible_multiphase->projection.Initialize_Dsd(particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,example.fluids_parameters.fuel_region);}
        if(number_of_regions==1 && example.fluids_parameters.use_strain){
            incompressible->Use_Strain();
            if(example.fluids_parameters.strain_boundary) incompressible->strain->Set_Custom_Boundary(*example.fluids_parameters.strain_boundary);}
        if(number_of_regions>=2 && example.fluids_parameters.use_multiphase_strain.Count_Matches(0)<number_of_regions){
            incompressible_multiphase->Use_Strain(example.fluids_parameters.use_multiphase_strain);
            if(example.fluids_parameters.strain_boundary){
                for(int i=0;i<number_of_regions;i++) if(incompressible_multiphase->strains(i))
                    incompressible_multiphase->strains(i)->Set_Custom_Boundary(*example.fluids_parameters.strain_boundary);}}}

    // set up Kolmogorov's spectrum
    if(example.fluids_parameters.kolmogorov) example.fluids_parameters.Initialize_Turbulence(time,example.frame_rate);

    // set up the initial state
    if(example.restart){
        example.Read_Output_Files_Fluids(current_frame);
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for advection later)
        if(number_of_regions>=1){
            particle_levelset_evolution->Set_Seed(2606);
            if(!example.fluids_parameters.write_particles) particle_levelset_evolution->Seed_Particles(time);
            particle_levelset_evolution->Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
        if(example.fluids_parameters.fire){
            for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Levelset(i).Compute_Normals(time);}
        if(example.fluids_parameters.use_flame_speed_multiplier) example.Get_Flame_Speed_Multiplier(0,time);
        if(example.fluids_parameters.use_reacting_flow){
            example.Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
            incompressible_multiphase->projection.Update_Phi_And_Move_Velocity_Discontinuity(fluid_collection.incompressible_fluid_collection.face_velocities,particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,time,true);}
        if(example.fluids_parameters.sph || example.fluids_parameters.use_sph_for_removed_negative_particles) example.Initialize_SPH_Particles();
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);}
    else{
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);
        if(example.fluids_parameters.water || example.fluids_parameters.fire || (example.fluids_parameters.compressible && number_of_regions)){
            example.Initialize_Phi();
            example.Adjust_Phi_With_Sources(time);
            if(number_of_regions>=2) particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Project_Levelset();
            if(!example.fluids_parameters.analytic_test) particle_levelset_evolution->Make_Signed_Distance();
            particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);}
        if(example.fluids_parameters.sph || example.fluids_parameters.use_sph_for_removed_negative_particles) example.Initialize_SPH_Particles();

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for averaging face velocities to nodes below)
        if(number_of_regions>=1){
            particle_levelset_evolution->Set_Seed(2606);
            particle_levelset_evolution->Seed_Particles(time);
            particle_levelset_evolution->Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
        if(example.fluids_parameters.use_soot){
            example.fluids_parameters.Initialize_Soot(time);
            example.Adjust_Soot_With_Sources(time);}
        if(example.fluids_parameters.use_density || example.fluids_parameters.use_temperature){
            example.fluids_parameters.Initialize_Density_And_Temperature(time);
            example.Adjust_Density_And_Temperature_With_Sources(time);}

        if(example.fluids_parameters.use_flame_speed_multiplier) example.Get_Flame_Speed_Multiplier(0,time);
        if(example.fluids_parameters.use_reacting_flow) example.Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
        example.Initialize_Velocities();
        if(example.fluids_parameters.compressible){
            example.Initialize_Euler_State();example.fluids_parameters.euler->Invalidate_Ghost_Cells();}
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);
        ARRAY<T,TV_INT> phi_dirichlet_regions;LEVELSET<TV> levelset_dirichlet_regions(grid,phi_dirichlet_regions);
        if(number_of_regions>=2 && example.fluids_parameters.dirichlet_regions.Number_True()>0){
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(example.fluids_parameters.dirichlet_regions,
                levelset_dirichlet_regions,example.fluids_parameters.flood_fill_for_bubbles);
            incompressible_multiphase->levelset_for_dirichlet_regions=&levelset_dirichlet_regions;}
        // TODO: cannot do this at initial time because we do not have two states to use to compute pseudo velocity
        //example.fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions((T)1./example.frame_rate,time); // fictitious dt

        if(example.fluids_parameters.use_reacting_flow)
            incompressible_multiphase->projection.Update_Phi_And_Move_Velocity_Discontinuity(fluid_collection.incompressible_fluid_collection.face_velocities,particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,time,true);
        if(number_of_regions==1){
            int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
            ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
            particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time,extrapolation_ghost_cells);
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,exchanged_phi_ghost,example.fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),
                &collision_bodies_affecting_fluid.face_neighbors_visible);}
        else if(number_of_regions>=2 && example.fluids_parameters.dirichlet_regions.Number_True()>0){
            incompressible_multiphase->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,phi_dirichlet_regions,
                example.fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);}}
    
    if(example.fluids_parameters.compressible){
        T fictitious_dt=(T)1./example.frame_rate;
        if(example.restart) euler->Fill_Ghost_Cells(fictitious_dt,time,3);
        example.fluids_parameters.euler_solid_fluid_coupling_utilities->Fill_Solid_Cells();
        if(number_of_regions){
            int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
            ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
            particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time,extrapolation_ghost_cells);
            euler->Fill_Ghost_Cells(fictitious_dt,time,extrapolation_ghost_cells);
            COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::Extrapolate_Compressible_State_Into_Incompressible_Region(fictitious_dt,time,(T)example.fluids_parameters.number_of_ghost_cells,extrapolation_ghost_cells,*euler->eos,euler->grid,
                exchanged_phi_ghost,fluid_collection.incompressible_fluid_collection.face_velocities,euler->U_ghost,euler->U);}}

    if(number_of_regions==1 && example.fluids_parameters.monitor_mass){
        example.fluids_parameters.mass=particle_levelset_evolution->Levelset_Advection(0).Approximate_Negative_Material();
        LOG::cout<<"Initial Material = "<<example.fluids_parameters.mass<<std::endl;}
    else if(number_of_regions>=2 && example.fluids_parameters.monitor_mass){
        example.fluids_parameters.masses.Resize(number_of_regions);
        for(int i=0;i<number_of_regions;i++){
            T mass_new=particle_levelset_evolution_multiple->Levelset_Advection(i).Approximate_Negative_Material();
            LOG::cout<<"Region "<<i<<": Initial material = "<<mass_new<<std::endl;
            example.fluids_parameters.masses(i)=mass_new;}}
    if(example.fluids_parameters.compressible && example.fluids_parameters.compressible_monitor_conservation_error){
        euler->Compute_Total_Conserved_Quantity(false,(T)0,euler->initial_total_conserved_quantity);

        TV solid_momentum;T solid_kinetic_energy,solid_potential_energy;
        example.solid_body_collection.Compute_Linear_Momentum(solid_momentum);
        example.solid_body_collection.Compute_Energy(time,solid_kinetic_energy,solid_potential_energy);

        for(int i=0;i<TV::m;i++) euler->initial_total_conserved_quantity[i+1]+=solid_momentum[i];
        euler->initial_total_conserved_quantity[TV::m+1]+=(solid_kinetic_energy+solid_potential_energy);}

    if(Simulate_Fluids()){
        if(example.fluids_parameters.compressible){
            example.fluids_parameters.euler_solid_fluid_coupling_utilities->Update_Cut_Out_Grid();
            if(euler->timesplit && euler->thinshell) example.fluids_parameters.euler_solid_fluid_coupling_utilities->Initialize_Collision_Data();}
        if(example.fluids_parameters.analytic_test) example.Get_Analytic_Velocities(time);}
}
//#####################################################################
// Function Rigid_Cluster_Fracture
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
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
        solids_evolution.Advance_One_Time_Step_Position(dt,time,!example.solids_fluids_parameters.mpi_solid_fluid || example.solids_fluids_parameters.mpi_solid_fluid->Solid_Node());
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
// Function Initialize_Fluids_Grids
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Initialize_Fluids_Grids()
{
    GRID<TV>& grid=*example.fluids_parameters.grid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    int number_of_regions=fluids_parameters.number_of_regions;
    fluids_parameters.Initialize_Grids();
    example.fluid_collection.Initialize_Grids();
    if(number_of_regions==1){
        fluids_parameters.particle_levelset_evolution->Initialize_Domain(fluids_parameters.p_grid);
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).Set_Band_Width((T)2*fluids_parameters.particle_half_bandwidth);}
    else if(number_of_regions>=2){
        fluids_parameters.particle_levelset_evolution_multiple->Initialize_Domain(fluids_parameters.p_grid,*fluids_parameters.collision_bodies_affecting_fluid,number_of_regions);
        for(int i=0;i<number_of_regions;i++)
            fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.particle_levelsets(i)->Set_Band_Width(
                (T)2*fluids_parameters.particle_half_bandwidth);}
    if(fluids_parameters.use_soot) fluids_parameters.soot_container.Initialize_Array();
    if(fluids_parameters.use_soot) fluids_parameters.soot_fuel_container.Initialize_Array();
    if(fluids_parameters.use_density) fluids_parameters.density_container.Initialize_Array();
    if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.Initialize_Array();
    if(fluids_parameters.sph || fluids_parameters.use_sph_for_removed_negative_particles) fluids_parameters.sph_evolution->Initialize_Grids(grid);
    if(fluids_parameters.compressible) fluids_parameters.euler->Initialize_Domain(grid);
    if(fluids_parameters.incompressible){
        fluids_parameters.incompressible->Initialize_Grids(grid);
        if(fluids_parameters.use_strain && fluids_parameters.incompressible->strain) fluids_parameters.incompressible->strain->Initialize_Grid(grid);}
    if(fluids_parameters.use_slip) dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(solids_evolution).Initialize_Grid_Arrays();
    if(number_of_regions>=2) for(int i=0;i<fluids_parameters.incompressible_multiphase->strains.m;i++) if(fluids_parameters.incompressible_multiphase->strains(i))
        fluids_parameters.incompressible_multiphase->strains(i)->Initialize_Grid(grid);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Advance_To_Target_Time(const T target_time)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    int number_of_regions=fluids_parameters.number_of_regions;
    const bool fluids=Simulate_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());

    T dt_full_advance=target_time-time;

    example.solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        Setup_Solids(time,substep);
        Setup_Fluids(time);
        T dt=Compute_Dt(time,target_time,done);
        
        if(fluids_parameters.compressible) euler->Save_State(euler->U_save,euler->euler_projection.face_velocities_save,euler->need_to_remove_added_internal_energy_save);
        
        Rigid_Cluster_Fracture(dt_full_advance,dt,substep);

        example.Preprocess_Substep(dt,time);

        if(fluids) example.Update_Fluid_Parameters(dt,time);
        if(fluids_parameters.use_flame_speed_multiplier) example.Get_Flame_Speed_Multiplier(dt,time);
        if(example.use_melting && number_of_regions==1) example.Update_Melting_Substep_Parameters(dt,time);

        if(fluids && Two_Way_Coupled()){
            if(fluids_parameters.compressible) fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(euler->euler_projection.elliptic_solver,euler->euler_projection.face_velocities,dt,time); // Put valid state on coupled faces.
            PHYSBAM_DEBUG_WRITE_SUBSTEP("integrate fluid forces for solid coupling",1);
            Integrate_Fluid_Non_Advection_Forces(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt/2,substep);/*F1*/}

        PHYSBAM_DEBUG_WRITE_SUBSTEP("solid position update",1);
        Solid_Position_Update(dt,substep);/*S1*/

        if(fluids){
            PHYSBAM_DEBUG_WRITE_SUBSTEP("object compatibility",1);
            if(Two_Way_Coupled()){  // Restore time n fluid state
                if(Simulate_Incompressible_Fluids()) incompressible->projection.Restore_After_Projection(example.fluid_collection.incompressible_fluid_collection.face_velocities);
                if(fluids_parameters.compressible) euler->Restore_State(euler->U_save,euler->euler_projection.face_velocities_save,euler->need_to_remove_added_internal_energy_save);}

            if(fluids_parameters.solid_affects_fluid && solids_fluids_parameters.use_leakproof_solve) Advance_Fluid_One_Time_Step_Implicit_Part_For_Object_Compatibility(last_dt,time-last_dt,substep);/*F2*/
            PHYSBAM_DEBUG_WRITE_SUBSTEP("advect fluid",1);
            Advect_Fluid(dt,substep);/*F3*/
            if(fluids_parameters.compressible && !fluids_parameters.use_slip){//slip does one sided interpolation, so dont need it 
                fluids_parameters.euler_solid_fluid_coupling_utilities->Fill_Solid_Cells();}}

        PHYSBAM_DEBUG_WRITE_SUBSTEP("solid velocity update",1);
        Solid_Velocity_Update(dt,substep,done);/*S2*/

        if(fluids){
            PHYSBAM_DEBUG_WRITE_SUBSTEP("project fluid at end of substep",1);
            Advance_Fluid_One_Time_Step_Implicit_Part(done,dt,substep);
            if(fluids_parameters.compressible && (!euler->timesplit || !euler->thinshell)) fluids_parameters.euler_solid_fluid_coupling_utilities->Fill_Solid_Cells();}/*F4*/

        example.Postprocess_Substep(dt,time);

        last_dt=restart_dt?restart_dt:dt;time+=last_dt;restart_dt=0;

        PHYSBAM_DEBUG_WRITE_SUBSTEP("END Substep %d",0,substep);}
}
//#####################################################################
// Function Integrate_Fluid_Non_Advection_Forces
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Integrate_Fluid_Non_Advection_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const int substep)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    int number_of_regions=fluids_parameters.number_of_regions;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=fluids_parameters.particle_levelset_evolution_multiple;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase=fluids_parameters.incompressible_multiphase;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities=fluids_parameters.euler_solid_fluid_coupling_utilities;

    LOG::SCOPE integration_scope("INTEGRATE NON ADVECTION FORCES","integrate non advection forces");

    if(fluids_parameters.use_body_force){
        if(Simulate_Incompressible_Fluids()) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.incompressible->force,dt,time);
        if(fluids_parameters.compressible) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.euler->force,dt,time);}

    LOG::Time("updating velocity (explicit part without convection)"); // TODO: fix vorticity confinement for thin shells
    if(fluids_parameters.analytic_test){time=particle_levelset_evolution->time;example.Get_Analytic_Velocities(particle_levelset_evolution->time+dt);}
    else{
        if(fluids_parameters.compressible){
            euler->Save_State(euler->U_save,euler->euler_projection.face_velocities_save,euler->need_to_remove_added_internal_energy_save);
            euler->Advance_One_Time_Step_Forces(dt,time);
            euler->Fill_Ghost_Cells(dt,time,example.fluids_parameters.number_of_ghost_cells);
            if(!euler->timesplit || !euler->thinshell) euler_solid_fluid_coupling_utilities->Fill_Solid_Cells(); // TODO(kwatra): see if can get rid of this, since one-sided interpolation is used in slip
            euler->Get_Dirichlet_Boundary_Conditions(dt,time);
            euler->euler_projection.Get_Pressure(euler->euler_projection.p_advected);}
        if(Simulate_Incompressible_Fluids()){
            if(fluids_parameters.fluid_affects_solid) incompressible->projection.Set_Up_For_Projection(face_velocities);
            // TODO(kwatra): Check if SPH case is handled properly.
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before viscosity",1);
            if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(fluids_parameters.projection))
                coupled_evolution->Apply_Viscosity(face_velocities,dt,time);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after viscosity",1);
            if(number_of_regions>=2)
                incompressible_multiphase->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution_multiple->phis,
                    &fluids_parameters.pseudo_dirichlet_regions,fluids_parameters.number_of_ghost_cells);
            else if(!fluids_parameters.sph) incompressible->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution->phi,fluids_parameters.number_of_ghost_cells);}
        fluids_parameters.Blend_In_External_Velocity(face_velocities,dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after integrate non advection forces",1);}
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
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
// Function Setup_Fluids
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Setup_Fluids(const T time)
{  
    const T fictitious_dt=(T)1./example.frame_rate;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities=fluids_parameters.euler_solid_fluid_coupling_utilities;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;

    fluids_parameters.collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
    if(fluids_parameters.number_of_regions) particle_levelset_evolution->Set_Number_Particles_Per_Cell(fluids_parameters.number_particles_per_cell);
    if(fluids_parameters.analytic_test) example.Get_Analytic_Velocities(time);
    if(fluids_parameters.compressible && euler->timesplit){
        euler->euler_projection.Compute_One_Over_rho_c_Squared();
        euler->euler_projection.Compute_Density_Weighted_Face_Velocities(fictitious_dt,time,euler->euler_projection.elliptic_solver->psi_N);
        if(euler->thinshell){
            euler_solid_fluid_coupling_utilities->Snapshot_State(euler->U);}}
}
//#####################################################################
// Function Solid_Position_Update
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Solid_Position_Update(const T dt,const int substep)
{
    Check_For_Interrupts(); // see if keyboard or other interrupts are waiting
    LOG::SCOPE scope("solids position update");

    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.solid_body_collection.deformable_body_collection;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;
    GRID<TV>& grid=*fluids_parameters.grid;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    const bool solids=Simulate_Solids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());
    const bool fluids=Simulate_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());

    // MELTING BROKEN!!!! if(example.use_melting){if(number_of_regions){example.Melting_Substep(dt,time);example.Modify_Fluid_For_Melting(dt,time);}else example.Melting_Substep(dt,time);}
    // if(solids_parameters.fracture_evolution) solids_parameters.fracture_evolution->Process_Rigid_Fracture(dt,time,example.solids_evolution,solids_evolution_callbacks);

    if(solids){
        if((solids_parameters.triangle_collision_parameters.repulsion_pair_update_count++)%solids_parameters.triangle_collision_parameters.repulsion_pair_update_frequency==0){
            example.solid_body_collection.deformable_body_collection.triangle_repulsions.Update_Faces_And_Hierarchies_With_Collision_Free_Positions(&deformable_body_collection.particles.X);
            example.solid_body_collection.deformable_body_collection.triangle_repulsions.Compute_Interaction_Pairs(deformable_body_collection.particles.X);}
        solids_evolution.kinematic_evolution.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame,time);
        solids_evolution.kinematic_evolution.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist,time,time);
        example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
        solids_evolution.Advance_One_Time_Step_Position(dt,time,!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());

        if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
            LOG::SCOPE scope("adjust velocity for self repulsion and self collisions");
            int repulsions,collisions_found;
            solids_evolution.Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(dt,time,repulsions,collisions_found,false);
            solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;}}
    else{
        solids_evolution_callbacks->Update_Solids_Parameters(time+dt);
        solids_evolution.kinematic_evolution.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame,time+dt);
        solids_evolution.kinematic_evolution.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist,time+dt,time+dt);
        example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
        dynamic_cast<NEWMARK_EVOLUTION<TV>&>(solids_evolution).Backward_Euler_Step_Velocity_Helper(dt/2,time,time,false);}
    // Exchange solid positions back to fluid nodes so that they can figure out effective velocity and do collidable advection
    if(solids_fluids_parameters.mpi_solid_fluid && Simulate_Fluids()){
        solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(example.solid_body_collection);}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("solid position updated",1);
    if(fluids){
        // MELTING BROKEN !!! if(example.use_melting && !number_of_regions) example.Modify_Fluid_For_Melting(dt,time);
        LOG::Time("rasterize objects");
        if(fluids_parameters.compressible && euler->timesplit && euler->thinshell) collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE,time+dt);
        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt);
        collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        LOG::Time("initializing swept occupied blocks");
        // swept occupied blocks
        if(euler) example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,euler->euler_projection.face_velocities.Max_Abs().Max(),true);
        else if(Simulate_Incompressible_Fluids()) example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,example.fluid_collection.incompressible_fluid_collection.face_velocities.Max_Abs().Max(),true);
        // MELTING BROKEN !!! if(example.use_melting && number_of_regions){ // re-rasterize changed bodies
        // MELTING BROKEN !!!    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
        // MELTING BROKEN !!!    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);}  // static occupied blocks
        
        if(euler && euler->timesplit && euler->thinshell){
            fluids_parameters.euler_solid_fluid_coupling_utilities->Update_Np1_Collision_Data(dt);
            if(fluids_parameters.compressible_rungekutta_order==3) fluids_parameters.euler_solid_fluid_coupling_utilities->Compute_Intermediate_Solid_Position_Data(dt);
            collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);}

        collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);  // static occupied blocks
        PHYSBAM_DEBUG_WRITE_SUBSTEP("body update",1);}
}
//#####################################################################
// Function Project_Fluid
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Project_Fluid(const T dt_projection,const T time_projection,const int substep)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    int number_of_regions=fluids_parameters.number_of_regions;
    GRID<TV>& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=fluids_parameters.particle_levelset_evolution_multiple;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase=fluids_parameters.incompressible_multiphase;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;

    // TODO: time+dt in this function may want to be time when this is used for initial object compatible velocity

    // project advection velocities
    if(fluids_parameters.use_non_zero_divergence){LOG::Time("getting Divergence");
        fluids_parameters.callbacks->Get_Divergence(fluids_parameters.incompressible->projection.divergence,dt_projection,time_projection);}

    if(!fluids_parameters.analytic_test){
        ARRAY<T,TV_INT> phi_for_dirichlet_regions;LEVELSET<TV> levelset_for_dirichlet_regions(grid,phi_for_dirichlet_regions); // for Dirichlet boundaries, surface tension and extrapolation
        LOG::Time("getting Neumann and Dirichlet boundary conditions");
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before get boundary conditions (Project_Fluid)",1);
        // todo: modify Get_Object_Velocities as necessary
        if(!fluids_parameters.compressible||number_of_regions==1)
            fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(incompressible->projection.elliptic_solver,fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection+dt_projection);
        // TODO: for the first timestep, dt_projection will be dt_old and thus zero for one of these calls.  What should we do?
        if(!fluids_parameters.fire && !fluids_parameters.surface_tension && !fluids_parameters.variable_surface_tension && dt_projection && (!fluids_parameters.compressible || number_of_regions==1))
            incompressible->projection.p*=dt_projection; // RESCALE PRESSURE FOR A BETTER INITIAL GUESS!
        if(number_of_regions==1){
            if(fluids_parameters.second_order_cut_cell_method) incompressible->projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
            if(fluids_parameters.surface_tension || fluids_parameters.variable_surface_tension){
                particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
                incompressible->Add_Surface_Tension(particle_levelset_evolution->Particle_Levelset(0).levelset,time_projection+dt_projection);}}
        else if(number_of_regions>=2 && fluids_parameters.dirichlet_regions.Number_True()>0){
            particle_levelset_evolution_multiple->Fill_Levelset_Ghost_Cells(time_projection+dt_projection);
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.dirichlet_regions,levelset_for_dirichlet_regions,
                fluids_parameters.flood_fill_for_bubbles);
            incompressible_multiphase->levelset_for_dirichlet_regions=&levelset_for_dirichlet_regions;

            if(fluids_parameters.second_order_cut_cell_method){
                incompressible_multiphase->projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
                incompressible_multiphase->projection.collidable_solver->Use_External_Level_Set(levelset_for_dirichlet_regions);}
            bool add_surface_tension=false;for(int i=0;i<number_of_regions;i++) for(int j=i+1;j<number_of_regions;j++)
                if(fluids_parameters.surface_tensions(i,j) && (fluids_parameters.dirichlet_regions(i) || fluids_parameters.dirichlet_regions(j))) add_surface_tension=true;
            PHYSBAM_ASSERT(!fluids_parameters.variable_surface_tension); // TODO: handle variable surface tensions here
            if(add_surface_tension) fluids_parameters.incompressible_multiphase->Add_Surface_Tension(levelset_for_dirichlet_regions,time_projection+dt_projection);}

        if(incompressible && incompressible->mpi_grid && fluids_parameters.use_sph_for_removed_negative_particles){
            for(int i=0;i<number_of_regions;i++){
                fluids_parameters.sph_evolution->Move_Particles_Off_Grid_Boundaries(particle_levelset_evolution->Particle_Levelset(i).removed_negative_particles,(T)1e-4*grid.dX.Max());
                Exchange_Ghost_Particles(*particle_levelset_evolution->Particle_Levelset(i).mpi_grid,particle_levelset_evolution->Particle_Levelset(i).template_removed_particles,
                    particle_levelset_evolution->Particle_Levelset(i).removed_negative_particles,example.fluids_parameters.number_of_ghost_cells,particle_levelset_evolution->Particle_Levelset(i));}}

        if(fluids_parameters.use_sph_for_removed_negative_particles && !fluids_parameters.sph_evolution->use_two_way_coupling){
            LOG::Time("updating removed negative particle velocities via sph");
            incompressible->projection.Set_Up_For_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve,true);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before one-way coupled sph solve",1);
            // TODO: check this dt
            fluids_parameters.sph_evolution->Make_Incompressible(particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles,example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after one-way coupled sph solve",1);
            incompressible->projection.Restore_After_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve,true);}

        LOG::Time("solving for the pressure and viscosity");
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before laplace solve (Project_Fluid)",1);
        if(euler){
            if(euler->timesplit && !euler->perform_rungekutta_for_implicit_part){
                LOG::Time("compressible implicit update");
                ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
                euler->Clamp_Internal_Energy(dt_projection,time_projection+dt_projection);
                fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(euler->euler_projection.elliptic_solver,euler->euler_projection.face_velocities,dt_projection,time_projection+dt_projection);
                if(number_of_regions==1){
                    euler->euler_projection.Compute_Density_Weighted_Face_Velocities(dt_projection,time_projection,euler->euler_projection.elliptic_solver->psi_N);
                    incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,incompressible_face_velocities,time_projection+dt_projection);
                    COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::Compute_Compressible_Incompressible_Face_Velocities(euler->grid,
                        incompressible_face_velocities,incompressible->projection.density,particle_levelset_evolution->phi,euler->U,euler->psi,
                        euler->euler_projection.face_velocities);
                    euler->euler_projection.Project(euler->euler_projection.face_velocities,dt_projection,time_projection);}
                else euler->Advance_One_Time_Step_Implicit_Part(dt_projection,time_projection);
                example.Apply_Isobaric_Fix(dt_projection,time_projection);
                euler->Remove_Added_Internal_Energy(dt_projection,time_projection);}
            if(!euler->timesplit&&number_of_regions==1){
                euler->Fill_Ghost_Cells(dt_projection,time_projection,example.fluids_parameters.number_of_ghost_cells);
                fluids_parameters.compressible_incompressible_coupling_utilities->Get_Dirichlet_Boundary_Conditions_For_Incompressible_Region(euler->grid,
                    euler->U_ghost,*(euler->eos),incompressible->projection.density,dt_projection);
                fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(incompressible->projection.elliptic_solver,fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection+dt_projection);
                incompressible->Advance_One_Time_Step_Implicit_Part(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection,
                    fluids_parameters.implicit_viscosity,0,fluids_parameters.use_levelset_viscosity,fluids_parameters.callbacks,fluids_parameters.print_viscosity_matrix);}
            if(euler->timesplit && !euler->perform_rungekutta_for_implicit_part){
                PHYSBAM_DEBUG_WRITE_SUBSTEP("after compressible implicit solve",1);}}
        else if(fluids_parameters.sph){
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before sph solve",1);
            incompressible->projection.Set_Up_For_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after sph solve",1);
            fluids_parameters.sph_evolution->Make_Incompressible(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection);
            incompressible->projection.Restore_After_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve);}
        else if(fluids_parameters.use_sph_for_removed_negative_particles && fluids_parameters.sph_evolution->use_two_way_coupling){
            incompressible->projection.Set_Up_For_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve);
            fluids_parameters.sph_evolution->Copy_Particle_Attributes_From_Array(particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles);
            fluids_parameters.sph_evolution->Set_Up_For_Projection(example.fluid_collection.incompressible_fluid_collection.face_velocities,time);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before two-way coupled sph solve",1);
            // TODO: sph people check this
            incompressible->Advance_One_Time_Step_Implicit_Part(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection,
                fluids_parameters.implicit_viscosity,0,fluids_parameters.use_levelset_viscosity,fluids_parameters.callbacks,fluids_parameters.print_viscosity_matrix);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after two-way coupled sph solve",1);
            fluids_parameters.sph_evolution->Postprocess_Particles(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection);
            fluids_parameters.sph_evolution->Copy_Particle_Attributes_To_Array(particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles);
            incompressible->projection.Restore_After_SPH(example.fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.sph_evolution->use_variable_density_solve);}
        else if(number_of_regions<2 || fluids_parameters.pseudo_dirichlet_regions.Number_True()==0)
            incompressible->Advance_One_Time_Step_Implicit_Part(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection,
                fluids_parameters.implicit_viscosity,0,fluids_parameters.use_levelset_viscosity,fluids_parameters.callbacks,fluids_parameters.print_viscosity_matrix);
        else{
            if(fluids_parameters.second_order_cut_cell_method) PHYSBAM_NOT_IMPLEMENTED(); // TAMAR: can you add the code for this when 2nd order works
            ARRAY<bool,TV_INT> psi_D_old=incompressible->projection.elliptic_solver->psi_D;ARRAY<T,TV_INT> p_old=incompressible->projection.p;
            ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=example.fluid_collection.incompressible_fluid_collection.face_velocities;
            incompressible_multiphase->Set_Dirichlet_Boundary_Conditions(particle_levelset_evolution_multiple->phis,fluids_parameters.pseudo_dirichlet_regions);
            // TODO: check me
            ARRAY<T,FACE_INDEX<TV::m> > air_velocities_save=face_velocities;
            ARRAY<T,TV_INT> phi_for_pseudo_dirichlet_regions;LEVELSET<TV> levelset_for_pseudo_dirichlet_regions(grid,phi_for_pseudo_dirichlet_regions);
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.pseudo_dirichlet_regions,
                levelset_for_pseudo_dirichlet_regions,false);
            incompressible->Extrapolate_Velocity_Across_Interface(face_velocities,phi_for_pseudo_dirichlet_regions);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before pseudo dirichlet solve 1",1);
            incompressible->Advance_One_Time_Step_Implicit_Part(face_velocities,dt_projection,time_projection,fluids_parameters.implicit_viscosity,0,fluids_parameters.use_levelset_viscosity,
                fluids_parameters.callbacks,fluids_parameters.print_viscosity_matrix);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after pseudo dirichlet solve 1",1);
            for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
                if(incompressible_multiphase->projection.elliptic_solver->psi_D(iterator.Cell_Index()))incompressible->projection.p(iterator.Cell_Index())=p_old(iterator.Cell_Index());
            incompressible->projection.elliptic_solver->psi_D=psi_D_old;
            ARRAY<bool,FACE_INDEX<TV::m> > psi_N_old=incompressible_multiphase->projection.elliptic_solver->psi_N;
            for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();int region_1,region_2;T phi_1,phi_2;
                particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
                if((!fluids_parameters.pseudo_dirichlet_regions(region_1) && !incompressible->projection.elliptic_solver->psi_D(cell_1)) ||
                    (!fluids_parameters.pseudo_dirichlet_regions(region_2) && !incompressible->projection.elliptic_solver->psi_D(cell_2)))
                    incompressible_multiphase->projection.elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;
                // TODO: drink me
                else face_velocities(iterator.Axis(),iterator.Face_Index())=air_velocities_save(iterator.Axis(),iterator.Face_Index());}
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before pseudo dirichlet solve 2",1);
            incompressible->Advance_One_Time_Step_Implicit_Part(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt_projection,time_projection,fluids_parameters.implicit_viscosity,0,fluids_parameters.use_levelset_viscosity,fluids_parameters.callbacks,fluids_parameters.print_viscosity_matrix);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after pseudo dirichlet solve 2",1);
            incompressible->projection.elliptic_solver->psi_N=psi_N_old;}
        if(fluids_parameters.use_sph_for_removed_negative_particles) particle_levelset_evolution->Delete_Particles_Outside_Grid();

        PHYSBAM_DEBUG_WRITE_SUBSTEP("after laplace solve (Project_Fluid)",1);
        example.Clamp_Velocities(time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after velocity clamping",1);
        if(!fluids_parameters.fire && !fluids_parameters.surface_tension && !fluids_parameters.variable_surface_tension && dt_projection && (!fluids_parameters.compressible || number_of_regions==1))
            incompressible->projection.p*=(1/dt_projection); // scale pressure back to get a real pressure
        if(fluids_parameters.compressible&&number_of_regions==1) incompressible->projection.p*=incompressible->projection.density;
        if(!fluids_parameters.compressible||number_of_regions==1)
            incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,example.fluid_collection.incompressible_fluid_collection.face_velocities,time_projection+dt_projection);
        if(fluids_parameters.second_order_cut_cell_method) incompressible->projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);

        LOG::Time("extrapolating velocity across interface");
        if(fluids_parameters.use_sph_for_removed_negative_particles && fluids_parameters.sph_evolution->use_two_way_coupling){
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,particle_levelset_evolution->phi,
                fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible,&fluids_parameters.sph_evolution->valid_particle_face_velocities);}
        else if(number_of_regions==1){
            // makes MPI work. TODO: make phi in levelset store 8 by default
            int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
            ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
            particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time_projection+dt_projection,extrapolation_ghost_cells);
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,exchanged_phi_ghost,
                fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
            if(fluids_parameters.use_strain){LOG::Time("extrapolating strain across interface");
                incompressible->strain->Extrapolate_Strain_Across_Interface(particle_levelset_evolution->phi);}

            if(example.fluids_parameters.compressible){
                euler->Fill_Ghost_Cells(dt_projection,time_projection,extrapolation_ghost_cells);
                COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::Extrapolate_Compressible_State_Into_Incompressible_Region(dt_projection,time_projection,(T)example.fluids_parameters.number_of_ghost_cells,extrapolation_ghost_cells,*euler->eos,euler->grid,
                    exchanged_phi_ghost,fluid_collection.incompressible_fluid_collection.face_velocities,euler->U_ghost,euler->U);}}
        else if(number_of_regions>=2){
            if(fluids_parameters.dirichlet_regions.Number_True()>0){
                incompressible_multiphase->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,
                    incompressible_multiphase->levelset_for_dirichlet_regions->phi,fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),
                    &collision_bodies_affecting_fluid.face_neighbors_visible);}
            for(int i=0;i<fluids_parameters.use_multiphase_strain.m;i++) if(fluids_parameters.use_multiphase_strain(i)){
                incompressible_multiphase->strains(i)->Extrapolate_Strain_Across_Interface(particle_levelset_evolution_multiple->phis(i));}}

        if(fluids_parameters.move_grid){fluids_parameters.Move_Grid(fluid_collection.incompressible_fluid_collection.face_velocities,time_projection+dt_projection);Initialize_Fluids_Grids();}}
}
//#####################################################################
// Function Advance_Fluid_One_Time_Step_Implicit_Part_For_Object_Compatibility
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Advance_Fluid_One_Time_Step_Implicit_Part_For_Object_Compatibility(const T dt_projection,const T time_projection,const int substep)
{
    // if(example.fluids_parameters.compressible) PHYSBAM_FATAL_ERROR("This currently doesn't work, as Fill_Solid_Cells is not aware of pseudo-velocities");
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;

    if(fluids_parameters.simulate && fluids_parameters.compressible && (!euler->timesplit || !euler->thinshell))
        fluids_parameters.euler_solid_fluid_coupling_utilities->Fill_Solid_Cells();
    if(Simulate_Incompressible_Fluids()){
        if(fluids_parameters.fluid_affects_solid){
            incompressible->projection.Set_Up_For_Projection(example.fluid_collection.incompressible_fluid_collection.face_velocities);
            incompressible->projection.Exchange_Pressures_For_Projection();}
        if(solids_fluids_parameters.use_leakproof_solve) Project_Fluid(dt_projection,time_projection,substep);
        if(fluids_parameters.fluid_affects_solid) incompressible->projection.Exchange_Pressures_For_Projection();}
}
//#####################################################################
// Function Advect_Fluid
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Advect_Fluid(const T dt,const int substep)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    int number_of_regions=fluids_parameters.number_of_regions;
    GRID<TV>& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=fluids_parameters.particle_levelset_evolution_multiple;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase=fluids_parameters.incompressible_multiphase;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities=fluids_parameters.euler_solid_fluid_coupling_utilities;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;

    LOG::SCOPE scalar_scope("SCALAR SCOPE");

    // TODO: make sure this is in the right place
    ARRAY<bool,FACE_INDEX<TV::m> > region_inside_pseudo_dirichlet_region;
    if(fluids_parameters.pseudo_dirichlet_regions.Number_True()>0){
        region_inside_pseudo_dirichlet_region.Resize(grid);
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int region=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region_Face(iterator.Axis(),iterator.Face_Index());
            region_inside_pseudo_dirichlet_region.Component(iterator.Axis())(iterator.Face_Index())=fluids_parameters.pseudo_dirichlet_regions(region);}}

    // important to compute ghost velocity values for particles near domain boundaries
    ARRAY<T,FACE_INDEX<TV::m> > projected_face_velocities_ghost,face_velocities_ghost;
    ARRAY<T,FACE_INDEX<TV::m> >* advection_face_velocities_ghost=0;
    if(incompressible){
        if(number_of_regions<=1){
            ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
            face_velocities_ghost.Resize(incompressible->grid,example.fluids_parameters.number_of_ghost_cells,false);
            if(Two_Way_Coupled() && solids_fluids_parameters.use_leakproof_solve){
                projected_face_velocities_ghost.Resize(incompressible->grid,example.fluids_parameters.number_of_ghost_cells,false);
                incompressible->boundary->Fill_Ghost_Faces(grid,face_velocities,projected_face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);
                advection_face_velocities_ghost=&projected_face_velocities_ghost;
                incompressible->boundary->Fill_Ghost_Faces(grid,incompressible->projection.face_velocities_save_for_projection,face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);}
            else{
                incompressible->boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);
                advection_face_velocities_ghost=&face_velocities_ghost;}}
        else if(number_of_regions>=2){
            face_velocities_ghost.Resize(incompressible->grid,example.fluids_parameters.number_of_ghost_cells,false);
            if(Two_Way_Coupled() && solids_fluids_parameters.use_leakproof_solve){
                projected_face_velocities_ghost.Resize(incompressible->grid,example.fluids_parameters.number_of_ghost_cells,false);
                particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,
                    particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,projected_face_velocities_ghost,time+dt);
                incompressible->boundary->Fill_Ghost_Faces(grid,incompressible->projection.face_velocities_save_for_projection,face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);
                incompressible->boundary->Fill_Ghost_Faces(grid,projected_face_velocities_ghost,projected_face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);
                advection_face_velocities_ghost=&projected_face_velocities_ghost;}
            else{
                particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,
                    particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,face_velocities_ghost,time+dt);
                incompressible->boundary->Fill_Ghost_Faces(grid,face_velocities_ghost,face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);
                advection_face_velocities_ghost=&face_velocities_ghost;}}}

    if(number_of_regions){
        if(number_of_regions==1) fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
        example.Adjust_Phi_With_Objects(time);
        LOG::Time("advecting levelset");
        if(fluids_parameters.use_reacting_flow) particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Compute_Normals(time);
        particle_levelset_evolution->Advance_Levelset(dt);
        if(number_of_regions==1) fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);
        example.Extrapolate_Phi_Into_Objects(time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after levelset advection",1);
        LOG::Time("advecting particles");
        if(fluids_parameters.analytic_test) particle_levelset_evolution->Advance_Particles(*advection_face_velocities_ghost,dt,fluids_parameters.analytic_test);
        else for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Particle_Levelset(i).Euler_Step_Particles(*advection_face_velocities_ghost,dt,time,true,true,false,fluids_parameters.analytic_test);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle advection",1);}
    if(fluids_parameters.sph){LOG::Time("advecting sph particles");fluids_parameters.sph_evolution->Euler_Step(dt,time);}

    example.Scalar_Advection_Callback(dt,time);

    LOG::Time("updating removed particle velocities");
    example.Modify_Removed_Particles_Before_Advection(dt,time);
    if(number_of_regions==1) fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);
    if(number_of_regions) particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
    if(!fluids_parameters.analytic_test) for(int k=0;k<number_of_regions;k++){
        PARTICLE_LEVELSET_UNIFORM<TV>& pls=particle_levelset_evolution->Particle_Levelset(k);
        LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
        if(pls.use_removed_positive_particles) for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
            for(int p=0;p<particles.Size();p++){
                TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(grid,*advection_face_velocities_ghost,X);
                if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=fluids_parameters.removed_positive_particle_buoyancy_constant*fluids_parameters.gravity.Normalized(); // buoyancy
                particles.V(p)=V;}}
        if(pls.use_removed_negative_particles) for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
            for(int p=0;p<particles.Size();p++) particles.V(p)+=dt*fluids_parameters.gravity; // ballistic
            if(fluids_parameters.use_body_force) for(int p=0;p<particles.Size();p++)
                particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(grid,incompressible->force,particles.X(p));}} // external forces
    if(fluids_parameters.sph) fluids_parameters.sph_evolution->sph_particles.V+=dt*incompressible->gravity; // ballistic
    if(fluids_parameters.use_soot){
        LOG::Time("advecting soot");
        fluids_parameters.Evolve_Soot(dt,time);}
    if(fluids_parameters.use_density || fluids_parameters.use_temperature){
        LOG::Time("advecting density and temperature");fluids_parameters.Evolve_Density_And_Temperature(dt,time);}

    if(fluids_parameters.use_reacting_flow && incompressible_multiphase->projection.dsd){
        LOG::Time("advancing detonation shock dynamics");
        LEVELSET_MULTIPLE<TV>& levelset_multiple=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple;
        levelset_multiple.Compute_Normals();levelset_multiple.Compute_Curvature();
        incompressible_multiphase->projection.dsd->Advance_One_Time_Step(particle_levelset_evolution_multiple->V,dt,time,fluids_parameters.number_of_ghost_cells);}
    else if(incompressible && incompressible->projection.dsd){
        LOG::Time("advancing detonation shock dynamics");
        LEVELSET<TV> &levelset=particle_levelset_evolution->Particle_Levelset(0).levelset;
        levelset.Compute_Normals();levelset.Compute_Curvature();
        incompressible->projection.dsd->Advance_One_Time_Step(particle_levelset_evolution->V,dt,time,fluids_parameters.number_of_ghost_cells);}

    example.Post_Velocity_Advection_Callback(dt,time);

    LOG::Time("updating velocity (explicit part)"); // TODO: fix vorticity confinement for thin shells
    if(fluids_parameters.analytic_test) example.Get_Analytic_Velocities(time+dt);
    else if(!fluids_parameters.compressible||number_of_regions==1){
        ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
        if(number_of_regions>=2){
            if(Two_Way_Coupled() && solids_fluids_parameters.use_leakproof_solve){
                incompressible_multiphase->Advance_One_Time_Step_Convection(dt,time,*advection_face_velocities_ghost,incompressible_multiphase->projection.face_velocities_save_for_projection,&fluids_parameters.pseudo_dirichlet_regions,fluids_parameters.number_of_ghost_cells);
                incompressible->projection.Restore_After_Projection(example.fluid_collection.incompressible_fluid_collection.face_velocities);
                Integrate_Fluid_Non_Advection_Forces(face_velocities,dt,substep);}
            else{
                incompressible_multiphase->Advance_One_Time_Step_Convection(dt,time,*advection_face_velocities_ghost,face_velocities,&fluids_parameters.pseudo_dirichlet_regions,fluids_parameters.number_of_ghost_cells);
                Integrate_Fluid_Non_Advection_Forces(face_velocities,dt,substep);}}
        else if(!fluids_parameters.sph){
            if(Two_Way_Coupled() && solids_fluids_parameters.use_leakproof_solve){
                PHYSBAM_DEBUG_WRITE_SUBSTEP("before forces",1);
                incompressible->Advance_One_Time_Step_Convection(dt,time,*advection_face_velocities_ghost,incompressible->projection.face_velocities_save_for_projection,fluids_parameters.number_of_ghost_cells);
                PHYSBAM_DEBUG_WRITE_SUBSTEP("before restore",1);
                incompressible->projection.Restore_After_Projection(example.fluid_collection.incompressible_fluid_collection.face_velocities);
                Integrate_Fluid_Non_Advection_Forces(face_velocities,dt,substep);
                PHYSBAM_DEBUG_WRITE_SUBSTEP("before convection",1);}
            else{
                PHYSBAM_DEBUG_WRITE_SUBSTEP("before forces",1);
                if(!fluids_parameters.stokes_flow) incompressible->Advance_One_Time_Step_Convection(dt,time,*advection_face_velocities_ghost,face_velocities,fluids_parameters.number_of_ghost_cells);
                Integrate_Fluid_Non_Advection_Forces(face_velocities,dt,substep);}}
        fluids_parameters.Blend_In_External_Velocity(face_velocities,dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after explicit part",1);}

    // compressible update
    if(fluids_parameters.simulate && fluids_parameters.compressible){
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before compressible explicit solve",1);
        euler_solid_fluid_coupling_utilities->Extract_Time_N_Data_For_Explicit_Fluid_Forces();
        LOG::Time("compressible explicit update");
        for(RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR> rk(euler->U,fluids_parameters.compressible_rungekutta_order,dt,time);rk.Valid();){
            if(euler->timesplit && euler->thinshell) euler_solid_fluid_coupling_utilities->Revert_Cells_Near_Interface(rk.substep);
            euler->Fill_Ghost_Cells(dt,time,3);
            if(fluids_parameters.compressible_apply_cavitation_correction){
                fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(euler->euler_cavitation_density.elliptic_solver,euler->euler_projection.face_velocities,dt,rk.time+dt);
                fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(euler->euler_cavitation_internal_energy.elliptic_solver,euler->euler_projection.face_velocities,dt,rk.time+dt);}
            euler->Advance_One_Time_Step_Explicit_Part(dt,rk.time,rk.substep,rk.order);
            for(CELL_ITERATOR<TV> iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
                assert(euler->U(cell_index)(0)>0);assert(EULER<TV>::e(euler->U(cell_index))>0);}
            if(euler->timesplit && euler->perform_rungekutta_for_implicit_part){assert(!euler->thinshell);
                euler->Get_Dirichlet_Boundary_Conditions(dt,rk.time);
                fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(euler->euler_projection.elliptic_solver,euler->euler_projection.face_velocities,dt,rk.time+dt);
                euler->Advance_One_Time_Step_Implicit_Part(dt,rk.time);}
            if(!euler->timesplit || euler->perform_rungekutta_for_implicit_part) example.Apply_Isobaric_Fix(dt,rk.time);
            euler->Remove_Added_Internal_Energy(dt,rk.time);
            rk.Next();
            euler->boundary->Apply_Boundary_Condition(euler->grid,euler->U,time);
            if(rk.substep!=rk.order-1) euler->Clamp_Internal_Energy(dt,rk.time);

            if(euler->timesplit && euler->thinshell){
                PHYSBAM_DEBUG_WRITE_SUBSTEP("before applying FSI update for near-interface cells",1);
                euler_solid_fluid_coupling_utilities->Update_Cells_Near_Interface(dt,rk.order,rk.substep);}
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after compressible explicit rk substep",1);}

        if(euler->timesplit && !euler->perform_rungekutta_for_implicit_part) euler->Get_Dirichlet_Boundary_Conditions(dt,time);
        if(euler->timesplit) euler->euler_projection.Get_Pressure(euler->euler_projection.p_advected);
        if(euler->timesplit && euler->thinshell) euler_solid_fluid_coupling_utilities->Compute_Post_Advected_Variables(); // TODO(jontg): MPI?
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after compressible explicit solve",1);}

    LOG::Time("effective velocity acceleration structures");
    // revalidate scalars and velocity in body's new position
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before scalar revalidation",1);
    collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);  // static occupied blocks
    collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
    example.Revalidate_Fluid_Scalars(); // uses visibility

    if(fluids_parameters.compressible && fluids_parameters.solid_affects_fluid) euler_solid_fluid_coupling_utilities->Update_Cut_Out_Grid();

    if(incompressible) example.Revalidate_Fluid_Velocity(fluid_collection.incompressible_fluid_collection.face_velocities); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after scalar revalidation",1);

    if(number_of_regions){
        if(fluids_parameters.use_reacting_flow){
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Project_Levelset();
            example.Set_Ghost_Density_And_Temperature_Inside_Flame_Core();}

        LOG::Time("modifying levelset");
        particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
        for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Particle_Levelset(i).Exchange_Overlap_Particles();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after filling level set ghost cells and exchanging overlap particles",1);
        if(fluids_parameters.use_sph_for_removed_negative_particles && fluids_parameters.sph_evolution->convert_particles_to_fluid){
            LOG::Time("converting particles to fluid and modifying levelset");
            fluids_parameters.sph_evolution->Modify_Levelset_And_Particles_To_Create_Fluid(time+dt,&face_velocities_ghost);}
        else particle_levelset_evolution->Modify_Levelset_And_Particles(&face_velocities_ghost);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after modify levelset and particles",1);
        example.Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
        if(number_of_regions>=2) particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Project_Levelset();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after revalidate phi",1);

        // fixing the velocity in regions that are no longer pseudo-Dirichlet
        if(fluids_parameters.pseudo_dirichlet_regions.Number_True()>0) for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const TV_INT index=iterator.Face_Index();
            int region1=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region(iterator.First_Cell_Index());
            int region2=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region(iterator.Second_Cell_Index());
            if(region_inside_pseudo_dirichlet_region.Component(axis)(index) && (!fluids_parameters.pseudo_dirichlet_regions(region1) || !fluids_parameters.pseudo_dirichlet_regions(region2)))
                example.fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(index)=particle_levelset_evolution_multiple->V.Component(axis)(index);}

        LOG::Time("adding sources");
        if(example.Adjust_Phi_With_Sources(time+dt)) particle_levelset_evolution->Make_Signed_Distance();
        particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
        LOG::Time("getting sources");
        ARRAY<bool,TV_INT>* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time+dt);
        if(source_mask){LOG::Time("reseeding sources");particle_levelset_evolution->Reseed_Particles(time+dt,0,source_mask);delete source_mask;}
        if(fluids_parameters.sph_evolution){LOG::Time("adding SPH particles for sources");example.Add_SPH_Particles_For_Sources(dt,time+dt);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after adding sources",1);

        LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(fluids_parameters.delete_fluid_inside_objects) fluids_parameters.Delete_Particles_Inside_Objects(time+dt);
        LOG::Time("deleting particles in local maxima");
        if(!fluids_parameters.analytic_test){
            if(number_of_regions==1) particle_levelset_evolution->Particle_Levelset(0).Delete_Particles_In_Local_Maximum_Phi_Cells(1);
            else for(int i=0;i<number_of_regions;i++) if(fluids_parameters.dirichlet_regions(i))
                particle_levelset_evolution->Particle_Levelset(i).Delete_Particles_In_Local_Maximum_Phi_Cells(-1);
            LOG::Time("deleting particles far from interface");
            for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Particle_Levelset(i).Delete_Particles_Far_From_Interface(); // uses visibility
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete particles far from interface",1);}

        LOG::Time("re-incorporating removed particles");
        example.Modify_Removed_Particles_Before_Reincorporation(dt,time+dt);
        for(int i=0;i<number_of_regions;i++){
            // TODO: if your particles fall entirely within the grid this shouldn't need ghost cells, but it may need them for MPI
            particle_levelset_evolution->Particle_Levelset(i).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
            if(particle_levelset_evolution->Particle_Levelset(i).use_removed_positive_particles || particle_levelset_evolution->Particle_Levelset(i).use_removed_negative_particles)
                particle_levelset_evolution->Particle_Levelset(i).Reincorporate_Removed_Particles(1,fluids_parameters.removed_particle_mass_scaling,
                    fluids_parameters.reincorporate_removed_particle_velocity?&example.fluid_collection.incompressible_fluid_collection.face_velocities:0,
                    !fluids_parameters.use_sph_for_removed_negative_particles || !fluids_parameters.sph_evolution->use_two_way_coupling);}
        example.Modify_Removed_Particles_After_Reincorporation(dt,time+dt);}
    else if(fluids_parameters.sph){LOG::Time("adding SPH particles for sources");example.Add_SPH_Particles_For_Sources(dt,time+dt);}

    ARRAY<T,TV_INT> phi_for_dirichlet_regions;LEVELSET<TV> levelset_for_dirichlet_regions(grid,phi_for_dirichlet_regions); // for Dirichlet boundaries, surface tension and extrapolation
    if(number_of_regions>=2){
        particle_levelset_evolution_multiple->Fill_Levelset_Ghost_Cells(time+dt);
        if(fluids_parameters.dirichlet_regions.Number_True()>0)
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.dirichlet_regions,levelset_for_dirichlet_regions,
                fluids_parameters.flood_fill_for_bubbles);
        incompressible_multiphase->levelset_for_dirichlet_regions=&levelset_for_dirichlet_regions;
        if(incompressible_multiphase->mpi_grid && fluids_parameters.flood_fill_for_bubbles) PHYSBAM_NOT_IMPLEMENTED();
        if(!fluids_parameters.use_reacting_flow)
            incompressible_multiphase->projection.poisson_collidable->Update_Internal_Level_Set(particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple);}
    if(fluids_parameters.use_reacting_flow)
        incompressible_multiphase->projection.Update_Phi_And_Move_Velocity_Discontinuity(example.fluid_collection.incompressible_fluid_collection.face_velocities,
            particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple,time+dt);

    // update ghost phi values
    if(number_of_regions>=1) particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after scalar update",1);
    scalar_scope.Pop();
}
//#####################################################################
// Function Solid_Velocity_Update
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Solid_Velocity_Update(const T dt,const int substep,const bool done)
{
    Check_For_Interrupts(); // see if keyboard or other interrupts are waiting
    LOG::SCOPE scope("solids velocity update");

    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    if(Simulate_Solids() || Two_Way_Coupled())
        solids_evolution.Advance_One_Time_Step_Velocity(dt,time,!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());

    // TODO: if this updates positions, it should go in S1.  Otherwise, here.
    solids_evolution.time+=dt;
    // if(solids_parameters.fracture_evolution) solids_parameters.fracture_evolution->Postprocess_Solids_Substep(solids_evolution.time,substep);
    solids_evolution_callbacks->Postprocess_Solids_Substep(solids_evolution.time,substep);
    solids_evolution_callbacks->Apply_Constraints(dt,solids_evolution.time);
}
//#####################################################################
// Function Advance_Fluid_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Advance_Fluid_One_Time_Step_Implicit_Part(const bool done,const T dt,const int substep)
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    GRID<TV>& grid=*fluids_parameters.grid;
    int number_of_regions=fluids_parameters.number_of_regions;
    EULER_UNIFORM<TV>* euler=fluids_parameters.euler;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=fluids_parameters.particle_levelset_evolution_multiple;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase=fluids_parameters.incompressible_multiphase;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible=fluids_parameters.incompressible;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;

    if(fluids_parameters.solid_affects_fluid && !fluids_parameters.fluid_affects_solid && !done && !project_at_frame_boundaries) return;
    assert(!fluids_parameters.fluid_affects_solid || fluids_parameters.solid_affects_fluid);
    if(Two_Way_Coupled()){
        if(fluids_parameters.use_slip)
            dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(solids_evolution).Apply_Pressure(example.fluid_collection.incompressible_fluid_collection.face_velocities,dt,time);
        else
            dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(solids_evolution).Apply_Pressure(dt,time);
        if(incompressible) incompressible->projection.p_save_for_projection.Copy(incompressible->projection.p); // save the good pressure for later

        PHYSBAM_DEBUG_WRITE_SUBSTEP("after apply pressure",1);
        if(fluids_parameters.compressible){
            example.Apply_Isobaric_Fix(dt,time);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after isobaric fix",1);}

        // TODO: this block below is identical with the end of Project_Fluid and could be broken out into a function
        if(fluids_parameters.compressible&&number_of_regions==1) incompressible->projection.p*=fluids_parameters.density;
        if(incompressible){
            incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,example.fluid_collection.incompressible_fluid_collection.face_velocities,time+dt);
        // TODO(kwatra): why do we need Set_Up_Second_Order_Cut_Cell_Method here.
            if(fluids_parameters.second_order_cut_cell_method) incompressible->projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);}

        LOG::Time("extrapolating velocity across interface");
        if(fluids_parameters.use_sph_for_removed_negative_particles && fluids_parameters.sph_evolution->use_two_way_coupling){
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,
                particle_levelset_evolution->phi,fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible,
                &fluids_parameters.sph_evolution->valid_particle_face_velocities);}
        else if(number_of_regions==1){
            // makes MPI work. TODO: make phi in levelset store 8 by default
            int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
            ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
            particle_levelset_evolution->Particle_Levelset(0).levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time+dt,extrapolation_ghost_cells);
            incompressible->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,exchanged_phi_ghost,
                fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
            if(fluids_parameters.use_strain){LOG::Time("extrapolating strain across interface");
                incompressible->strain->Extrapolate_Strain_Across_Interface(particle_levelset_evolution->phi);}
            if(example.fluids_parameters.compressible){
                euler->Fill_Ghost_Cells(dt,time,extrapolation_ghost_cells);
                COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::Extrapolate_Compressible_State_Into_Incompressible_Region(dt,time,(T)example.fluids_parameters.number_of_ghost_cells,extrapolation_ghost_cells,*euler->eos,euler->grid,
                    exchanged_phi_ghost,fluid_collection.incompressible_fluid_collection.face_velocities,euler->U_ghost,euler->U);}
            incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,example.fluid_collection.incompressible_fluid_collection.face_velocities,time+dt);
        }
        else if(number_of_regions>=2){
            if(fluids_parameters.dirichlet_regions.Number_True()>0){
                incompressible_multiphase->Extrapolate_Velocity_Across_Interface(example.fluid_collection.incompressible_fluid_collection.face_velocities,incompressible_multiphase->levelset_for_dirichlet_regions->phi,
                    fluids_parameters.enforce_divergence_free_extrapolation,(T)example.fluids_parameters.number_of_ghost_cells,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
                incompressible_multiphase->boundary->Apply_Boundary_Condition_Face(incompressible->grid,example.fluid_collection.incompressible_fluid_collection.face_velocities,time+dt);}
            for(int i=0;i<fluids_parameters.use_multiphase_strain.m;i++) if(fluids_parameters.use_multiphase_strain(i)){
                incompressible_multiphase->strains(i)->Extrapolate_Strain_Across_Interface(particle_levelset_evolution_multiple->phis(i));}}

        if(fluids_parameters.move_grid){fluids_parameters.Move_Grid(fluid_collection.incompressible_fluid_collection.face_velocities,time+dt);Initialize_Fluids_Grids();}}
    else if(!fluids_parameters.solid_affects_fluid || (done && project_at_frame_boundaries)){ // TODO: use actual velocities rather than effective velocities if one-way coupled && done
        Project_Fluid(dt,time,substep);}

    if(fluids_parameters.compressible){
        if(fluids_parameters.compressible_monitor_conservation_error){
            VECTOR<T,TV::m+2> new_total_conserved_quantity;
            TV solid_momentum;T solid_kinetic_energy,solid_potential_energy;
            euler->Compute_Total_Conserved_Quantity(true,dt,new_total_conserved_quantity);
            example.solid_body_collection.Compute_Linear_Momentum(solid_momentum);
            example.solid_body_collection.Compute_Energy(time,solid_kinetic_energy,solid_potential_energy);

            for(int i=0;i<TV::m;i++) new_total_conserved_quantity[i+1]+=solid_momentum[i];
            new_total_conserved_quantity[TV::m+1]+=(solid_kinetic_energy+solid_potential_energy);
            LOG::cout<<"solid_momentum="<<solid_momentum<<", solid_energy="<<solid_kinetic_energy+solid_potential_energy<<std::endl;

            VECTOR<T,TV::m+2>& initial_total_conserved_quantity=euler->initial_total_conserved_quantity;
            LOG::cout<<"Conserved variable error at time="<<time<<std::endl;
            for(int d=0;d<TV::m+2;d++){
                T total_error=new_total_conserved_quantity(d)-initial_total_conserved_quantity(d);
                T relative_error=(initial_total_conserved_quantity(d)!=0)?total_error/initial_total_conserved_quantity(d):0;
                LOG::cout<<"Conserved variable dimension="<<d<<", initial value="<<initial_total_conserved_quantity(d)<<", new value="<<new_total_conserved_quantity(d)<<
                    ", total error="<<total_error<<", relative error="<<relative_error<<std::endl;}}

        if(fluids_parameters.compressible_log_extremas){
            T max_density=-FLT_MAX,min_density=FLT_MAX;
            TV max_velocity=TV();
            TV_INT min_cell_index,max_cell_index,max_velocity_index;
            for(CELL_ITERATOR<TV> iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
                if(euler->psi(cell_index)){
                    T current_density=euler->U(cell_index)[0];
                    TV current_velocity=EULER<TV>::Get_Velocity(euler->U(cell_index));
                    if(current_density<min_density){min_density=current_density;min_cell_index=cell_index;}
                    if(current_density>max_density){max_density=current_density;max_cell_index=cell_index;}
                    if(current_velocity.Magnitude()>max_velocity.Magnitude()){max_velocity=current_velocity;max_velocity_index=cell_index;}}}
            LOG::cout<<"time="<<time+dt<<", Density extremas: min at "<<min_cell_index<<"="<<min_density<<", max at "<<max_cell_index<<"="<<max_density<<std::endl;
            LOG::cout<<"time="<<time+dt<<", Maximum velocity ="<<max_velocity.Magnitude()<<": "<<max_velocity<<std::endl;}}
    else
        LOG::cout<<"Maximum face velocity = ("<<fluid_collection.incompressible_fluid_collection.face_velocities.Max_Abs().Magnitude()<<": "<<fluid_collection.incompressible_fluid_collection.face_velocities.Max_Abs()<<std::endl;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Postprocess_Frame(const int frame)
{
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple=example.fluids_parameters.particle_levelset_evolution_multiple;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    int number_of_regions=example.fluids_parameters.number_of_regions;

    if(number_of_regions>=1){
        example.Postprocess_Phi(time); // actually changes phi !!!

        if(particle_levelset_evolution->use_particle_levelset && (frame-example.first_frame)%example.fluids_parameters.reseeding_frame_rate==0){
            LOG::Time("Reseeding... ");
            particle_levelset_evolution->Reseed_Particles(time);
            particle_levelset_evolution->Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);
            LOG::Stop_Time();}}

    if(example.fluids_parameters.monitor_mass){
        if(number_of_regions==1){
            T mass_new=particle_levelset_evolution->Levelset_Advection(0).Approximate_Negative_Material();
            LOG::cout<<"Material = "<<mass_new<<" - change = "<<mass_new-example.fluids_parameters.mass<<std::endl;
            example.fluids_parameters.mass=mass_new;}
        else if(number_of_regions>=2){
            for(int i=0;i<example.fluids_parameters.number_of_regions;i++){
                T mass_new=particle_levelset_evolution_multiple->Levelset_Advection(i).Approximate_Negative_Material();
                LOG::cout<<"Region "<<i<<": Material = "<<mass_new<<", change = "<<mass_new-example.fluids_parameters.masses(i)<<std::endl;
                example.fluids_parameters.masses(i)=mass_new;}}}

    SOLIDS_FLUIDS_DRIVER<TV>::Postprocess_Frame(frame);
    solids_evolution.Postprocess_Frame(frame);

    if(solids_evolution.solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics)
        solids_evolution.rigid_body_collisions->Print_Interpenetration_Statistics();
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection=example.fluid_collection;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;
    int number_of_regions=fluids_parameters.number_of_regions;
    bool fluids=Simulate_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());

    T fluids_dt=FLT_MAX;
    if(fluids){
        if(fluids_parameters.use_reacting_flow) for(int i=0;i<fluids_parameters.number_of_regions;i++)
            fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelsets(i)->Compute_Normals(time);
        if((number_of_regions==0 && !fluids_parameters.sph && !fluids_parameters.compressible) || (fluids_parameters.incompressible && fluids_parameters.incompressible->strain))
            fluids_dt=fluids_parameters.cfl*fluids_parameters.incompressible->CFL(fluid_collection.incompressible_fluid_collection.face_velocities);
        else if(number_of_regions>=2 && fluids_parameters.incompressible_multiphase->nonzero_surface_tension)
            fluids_dt=fluids_parameters.cfl*fluids_parameters.incompressible_multiphase->CFL(fluid_collection.incompressible_fluid_collection.face_velocities,fluids_parameters.implicit_viscosity,false);
        else if(fluids_parameters.compressible) fluids_dt=fluids_parameters.euler->cfl_number*fluids_parameters.euler->CFL(time);

        if(number_of_regions>=1){T dt_levelset=fluids_parameters.particle_levelset_evolution->CFL(true,fluids_parameters.analytic_test);fluids_dt=min(fluids_dt,dt_levelset);}
        if(fluids_parameters.sph || fluids_parameters.use_sph_for_removed_negative_particles) fluids_dt=min(fluids_dt,fluids_parameters.cfl*fluids_parameters.sph_evolution->CFL());
        LOG::cout<<"before example cfl clamping:  "<<fluids_dt<<std::endl;
        example.Limit_Dt(fluids_dt,time);
        if(fluids_parameters.mpi_grid && !solids_fluids_parameters.mpi_solid_fluid) fluids_parameters.mpi_grid->Synchronize_Dt(fluids_dt);}

    // solids dt
    T solids_dt=FLT_MAX;
    if(solids_evolution.Use_CFL()) solids_dt=min(solids_dt,example.solid_body_collection.CFL(solids_parameters.verbose_dt));
    solids_dt=min(solids_dt,solids_evolution_callbacks->Constraints_CFL());
    if(solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) solids_dt=min(solids_dt,example.solid_body_collection.rigid_body_collection.CFL_Rigid(solids_parameters.rigid_body_evolution_parameters,solids_parameters.verbose_dt));
    solids_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);
    if(example.solid_body_collection.deformable_body_collection.mpi_solids)
        solids_dt=example.solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Min_Global(solids_dt);
    if(solids_fluids_parameters.mpi_solid_fluid){
        solids_dt=solids_fluids_parameters.mpi_solid_fluid->Reduce_Min(solids_dt);
        fluids_dt=solids_fluids_parameters.mpi_solid_fluid->Reduce_Min(fluids_dt);}

    if(example.fixed_dt){fluids_dt=example.fixed_dt;solids_dt=example.fixed_dt;}
    if(example.max_dt){fluids_dt=min(fluids_dt,example.max_dt);solids_dt=min(solids_dt,example.max_dt);}
    T dt=min(fluids_dt,solids_dt);
    if(Simulate_Fluids())
        LOG::cout<<"fluids_dt = "<<fluids_dt<<", solids_dt = "<<solids_dt<<" dt="<<dt<<std::endl;
    else
        LOG::cout<<"dt = solids_dt = "<<dt<<std::endl;
    if(example.abort_when_dt_below && dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(LOG::sprintf("dt too small (%g < %g)",dt,example.abort_when_dt_below));
    done=false;
    SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,solids_parameters.min_dt);
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("writing output files");
    Create_Directory(example.output_directory);
    Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);

    int number_of_regions=example.fluids_parameters.number_of_regions;

    if(number_of_regions==1) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Write_Output_Files(frame);
    if(number_of_regions==1) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);

    if(number_of_regions>=1){
        GRID<TV>& grid=*example.fluids_parameters.grid;
        for(int i=0;i<number_of_regions;i++){
            int number_of_positive_particles=0,number_of_negative_particles=0,number_of_removed_positive_particles=0,number_of_removed_negative_particles=0;
            PARTICLE_LEVELSET_UNIFORM<TV>* pls=0;
            if(number_of_regions==1) pls=&example.fluids_parameters.particle_levelset_evolution->Particle_Levelset(0);
            else pls=example.fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.particle_levelsets(i);
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
            LOG::cout<<number_of_removed_positive_particles<<" positive and "<<number_of_removed_negative_particles<<" negative removed particles "<<std::endl;}}

    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<float,1> >;
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<float,2> >;
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<float,3> >;
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<double,1> >;
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<double,2> >;
template class SOLIDS_FLUIDS_DRIVER_UNIFORM<VECTOR<double,3> >;
}
