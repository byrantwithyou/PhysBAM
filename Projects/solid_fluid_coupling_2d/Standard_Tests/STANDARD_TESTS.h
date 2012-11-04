//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Falling deformable sphere and smoke
//   2. Falling rigid bodies on smoke and interacting
//   3. Broken
//   4. Deformable closed thin shell on smoke (Broken)
//   5. Cloth falling on smoke (Broken)
//   6. Broken
//   7. Many falling spheres...on smoke (Broken)
//   8. Balloon (Broken)
//   9. Cork rising in a column of water (Broken)
//   10. Broken
//   11. Broken
//   12. Broken
//   13. Broken
//   14. Broken
//   15. Broken
//   16. Broken
//   17. Broken
//   18. Broken
//   19. Broken
//   20. Strange example
//   21. Broken
//   22. Broken
//   23. Broken
//   24. Broken
//   25. Broken
//   26. Broken
//   27. Broken
//   28. Broken
//   29. Broken
//   30. High-mass gravity test
//   31. Advection test
//   32. Rigid body falling under gravity
//   33. Flexible beam test
//   34. Vibrating circle
//   35. Flexible filament
//   36. Simple fluid test (Broken)
//   37. Deformable advection test
//   38. Coupled viscosity test
//   39. Refine circle (Broken)
//   40. Analytic test
//   41. Flow past fixed cylinder
//   42. Sanity test - no advection, no viscosity
//   43. Oscillating disk (Zhao, Freund & Moser)
//   44. Flow past fixed cylinder; shedding
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <fstream>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef ARRAY<T,FACE_INDEX<2> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<bool,FACE_INDEX<2> > T_FACE_ARRAYS_BOOL;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::output_directory;using BASE::last_frame;
    using BASE::frame_rate;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::Add_To_Fluid_Simulation;
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::data_directory;using BASE::Add_Thin_Shell_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    GRID<TV> mattress_grid,mattress_grid_stub;
    int deformable_object_id;
    T velocity_multiplier;
    T mass_multiplier;
    int rigid_body_count;
    ARRAY<int> deformable_object_enslaved_nodes;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
    ARRAY<PAIR<int,TV> > constrained_node_positions;
    int heavy_sphere_index,light_sphere_index;
    T heavy_sphere_drop_time,light_sphere_drop_time,heavy_sphere_initial_height,light_sphere_initial_height;
    T balloon_initial_radius;
    int rigid_body_id;

    T velocity_angle;
    bool flow_particles;
    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_each_matrix;
    bool output_iterators;
    bool use_viscous_forces;
    int circle_refinement;
    T scale_length;
    bool use_solid;
    T fluid_gravity;
    T solid_gravity;
    T solid_width;
    T solid_density;
    int widen_domain;
    T period;
    T max_dt;
    ARRAY<TV> sample_points;
    int beam_elements_width;
    int beam_elements_length;
    int solid_resolution;
    T analytic_solution;
    bool use_viscosity,use_solid_width,use_solid_density;

    GEOMETRY_PARTICLES<TV> debug_particles;
    TRIANGULATED_AREA<T>* test_43_triangulated_area;
    FINITE_VOLUME<TV,2>* finite_volume;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),
        solids_tests(*this,solid_body_collection),deformable_object_id(0),mass_multiplier(1),stiffness_multiplier((T)1),damping_multiplier((T)1),
        bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),rigid_body_id(0),flow_particles(false),run_self_tests(false),print_poisson_matrix(false),
        print_index_map(false),print_matrix(false),print_each_matrix(false),output_iterators(false),circle_refinement(0),scale_length(1),use_solid(false),fluid_gravity((T)9.8),solid_gravity((T)9.8),
        solid_width((T).1111),solid_density(100),widen_domain(0),period(10),max_dt(0),beam_elements_width(2),beam_elements_length(40),solid_resolution(216),
        use_viscosity(false),use_solid_width(false),use_solid_density(false)
    {
        LOG::cout<<std::setprecision(16);
    }

    // Unused callbacks
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    fluids_parameters.incompressible_iterations=3000;
    parse_args->Add("-mass",&mass_multiplier,"value","mass_multiplier");
    parse_args->Add("-cg_iterations",&fluids_parameters.incompressible_iterations,"value","cg iterations");
    parse_args->Add("-slip",&fluids_parameters.use_slip,"use slip");
    parse_args->Add("-viscosity",&fluids_parameters.viscosity,&use_viscosity,"value","fluid viscosity");
    parse_args->Add("-test_system",&run_self_tests,"Run self tests");
    parse_args->Add("-print_poisson_matrix",&print_poisson_matrix,"print_poisson_matrix");
    parse_args->Add("-print_index_map",&print_index_map,"print_index_map");
    parse_args->Add("-print_matrix",&print_matrix,"print_matrix");
    parse_args->Add("-print_each_matrix",&print_each_matrix,"print_each_matrix");
    parse_args->Add("-output_iterators",&output_iterators,"output_iterators");
    parse_args->Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
    parse_args->Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"preconditioner");
    parse_args->Add("-leakproof",&solids_fluids_parameters.use_leakproof_solve,"use leakproof solve");
    parse_args->Add("-use_viscous_forces",&use_viscous_forces,"use_viscous_forces");
    parse_args->Add("-scale_length",&scale_length,"value","scale_length");
    parse_args->Add("-use_solid",&use_solid,"use_solid");
    parse_args->Add("-fluid_gravity",&fluid_gravity,"value","fluid_gravity");
    parse_args->Add("-solid_gravity",&solid_gravity,"value","solid_gravity");
    parse_args->Add("-solid_width",&solid_width,&use_solid_width,"value","solid_width");
    parse_args->Add("-solid_density",&solid_density,&use_solid_density,"value","solid_density");
    parse_args->Add("-widen_domain",&widen_domain,"value","widen_domain");
    parse_args->Add("-period",&period,"value","period");
    parse_args->Add("-max_dt",&max_dt,"value","maximum dt");
    parse_args->Add("-stokes",&fluids_parameters.stokes_flow,"disable advection");
    parse_args->Add("-beam_res_w",&beam_elements_width,"value","beam_elements_width");
    parse_args->Add("-beam_res_h",&beam_elements_length,"value","beam_elements_length");
    parse_args->Add("-solid_resolution",&solid_resolution,"value","solid_resolution");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    last_frame=100;
    frame_rate=24;

    fluids_parameters.cfl=(T).9;
    fluids_parameters.confinement_parameter=(T).04;        
    fluids_parameters.rho_bottom=1;
    fluids_parameters.rho_top=(T).65;
    fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.temperature_container.Set_Cooling_Constant(0);
    solids_parameters.implicit_solve_parameters.cg_iterations=fluids_parameters.incompressible_iterations;
    fluids_parameters.use_coupled_implicit_viscosity=use_viscous_forces;

    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;

    //mattress_grid=GRID<TV>(TV_INT(3,3),RANGE<TV>(TV((T).25,(T).81),TV((T).75,(T)1.05)));
    mattress_grid=GRID<TV>(TV_INT(3,3),RANGE<TV>(TV((T).26,(T).81),TV((T).76,(T)1.04)));
    mattress_grid_stub=GRID<TV>(TV_INT(2,2),RANGE<TV>(TV((T)20.25,(T)-21.05),TV((T)20.75,(T)-21.01)));
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_post_cg_constraints=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=200;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    //solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
    velocity_multiplier=(T)1;

    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_preconditioner_for_slip_system=true;
    solids_fluids_parameters.use_leakproof_solve=false;
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
        
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    // T default_removed_positive_particle_buoyancy_constant=fluids_parameters.removed_positive_particle_buoyancy_constant;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;
    //solid_body_collection.print_residuals=true;

    fluids_parameters.gravity=0;

    switch(test_number){
        case 1:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T)1;
            (*fluids_parameters.grid).Initialize(TV_INT(resolution+1,(int)(1.5*resolution)+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)1.5)));
            fluids_parameters.use_density=fluids_parameters.use_temperature=false;
            break;
        case 2:
            //last_frame=1000;
            (*fluids_parameters.grid).Initialize(TV_INT(resolution+1,(int)(1.5*resolution)+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)1.5)));
            //fluids_parameters.density=(T)300;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T)1000;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
            velocity_multiplier=10;
            fluids_parameters.use_density=fluids_parameters.use_temperature=false;
            break;
        case 3:
        case 4:
        case 5:
        case 6:
            fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,15*resolution+1),RANGE<TV>(TV(0,0),TV(1,1.5)));
            last_frame=1000;
            fluids_parameters.density=(T)1000;
            rigid_body_count=6;
            velocity_multiplier=(T)15;
            fluids_parameters.use_body_force=true;
            break;
        case 7:
            fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,15*resolution+1),RANGE<TV>(TV(0,0),TV(1,1.5)));
            last_frame=1000;
            fluids_parameters.density=(T)1000;
            rigid_body_count=6;
            velocity_multiplier=(T)15;
            fluids_parameters.use_body_force=true;
            break;
        case 8:
            fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,15*resolution+1),RANGE<TV>(TV(0,0),TV(1,1.5)));
            last_frame=1000;
            fluids_parameters.density=(T)10;
            velocity_multiplier=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
            balloon_initial_radius=(T).2;
            light_sphere_drop_time=(T).5;
            heavy_sphere_drop_time=(T)1.5;
            heavy_sphere_initial_height=(T).25;
            light_sphere_initial_height=(T).25;
            fluids_parameters.use_body_force=true;
            break;
        case 9:
            last_frame=5000;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T)1000;
            fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=true;//fluids_parameters.domain_walls[1][1]=true;
            (*fluids_parameters.grid).Initialize(TV_INT(2*resolution+1,5*resolution+1),RANGE<TV>(TV((T)-2,(T)0),TV((T)2,(T)10)));
            fluids_parameters.use_body_force=true;
            break;
        case 10: // flow chamber
        case 11:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,(int)(1.5*resolution)+1),RANGE<TV>(TV(0,0),TV(1,1.5)));
            last_frame=1000;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T)1000;
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            velocity_multiplier=5;
            fluids_parameters.use_body_force=true;
            break;
        case 12:
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            (*fluids_parameters.grid).Initialize(TV_INT(4*resolution+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)4,(T)1)));
            // (*fluids_parameters.grid).Initialize(TV_INT(4*resolution+1,2*resolution+1),RANGE<TV>(TV((T)0,(T)-0.5),TV((T)4,(T)1.5)));
            flow_particles=true;
            fluids_parameters.use_body_force=true;
            
            break;
        case 13:
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            (*fluids_parameters.grid).Initialize(TV_INT(2*resolution+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)));
            velocity_angle=pi/16;
            fluids_parameters.density=(T)1000;
            fluids_parameters.use_body_force=true;
            break;
        case 20:
            last_frame=1000;
            (*fluids_parameters.grid).Initialize(TV_INT((int)(1.5*resolution)+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)));
            //fluids_parameters.density=(T)300;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T)1000;
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;
            velocity_multiplier=10;
            fluids_parameters.use_density=fluids_parameters.use_temperature=false;
            break;
        case 30:
            (*fluids_parameters.grid).Initialize(TV_INT((int)(1.5*resolution)+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)));
            last_frame=1000;
            fluids_parameters.gravity=(T)9.8;
            fluids_parameters.density=(T).1;
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
//                velocity_multiplier=10;
            fluids_parameters.use_density=fluids_parameters.use_temperature=false;
            break;
        case 31:
        case 34:
        case 35:
        case 36:
        case 37:
        case 38:
        case 39:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,(int)(1.5*resolution)+1),RANGE<TV>(TV(0,0),TV(1,1.5)));
            for(int axis=0;axis<TV::dimension;axis++)for(int side=0;side<2;side++) fluids_parameters.domain_walls[axis][side]=false;
            break;
        case 32:
            fluids_parameters.grid->Initialize(TV_INT(resolution,4*resolution),RANGE<TV>(TV(0,0),TV((T).04*scale_length,(T).16*scale_length)),true);
            break;
        case 33:
            fluids_parameters.grid->Initialize(TV_INT(2,1)*resolution,RANGE<TV>(TV(),TV(2,1)),true);
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            break;
        case 40:
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
            {T extend=(T)1/resolution*widen_domain*scale_length;
                (*fluids_parameters.grid).Initialize(TV_INT(resolution+1+2*widen_domain,resolution+1),RANGE<TV>(TV(-extend,(T)0),TV(1+extend,(T)scale_length)));}
            if(!use_viscosity) fluids_parameters.viscosity=100;
            if(!use_solid_density) solid_density=150;
            if(!use_solid_width) solid_width=(T)1/3;
            break;
        case 41:
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            (*fluids_parameters.grid).Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0,0),TV(4,4)));
            break;
        case 42:
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            (*fluids_parameters.grid).Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            break;
        case 43:
            fluids_parameters.grid->Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0,0),TV(1,1)));
            for(int axis=0;axis<TV::dimension;axis++)for(int side=0;side<2;side++) fluids_parameters.domain_walls[axis][side]=false;
            break;
        case 44:
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
            (*fluids_parameters.grid).Initialize(TV_INT((int)(2.5*resolution)+1,resolution+1),RANGE<TV>(TV(-(T)2.5,-(T)3.5),TV(15,(T)3.5)));
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}
        
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);
    if(fluids_parameters.use_slip)
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d_slip",test_number,resolution);
    else
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d",test_number,resolution);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(debug_particles.Size()){
        FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),frame));
        FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles);
        debug_particles.Delete_All_Elements();}

    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution)){
        UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV> iterator_info(*fluids_parameters.collision_bodies_affecting_fluid);
        if(frame==1) evolution->Output_Iterators(this->stream_type,output_directory.c_str(),0);
        evolution->Output_Iterators(this->stream_type,output_directory.c_str(),frame);}
    if(test_number==40){
        T v=fluid_collection.incompressible_fluid_collection.face_velocities(FACE_INDEX<2>(2,fluids_parameters.grid->counts/2));
        if(solid_body_collection.rigid_body_collection.rigid_body_particle.frame.m>=3) v=solid_body_collection.rigid_body_collection.rigid_body_particle.twist(2).linear.y;
        LOG::cout<<"middle-velocity "<<v<<"   error from analytic solution "<<(v/analytic_solution-1)<<std::endl;}
    if(test_number==32){
        LOG::cout<<"solid-velocity "<<solid_body_collection.rigid_body_collection.rigid_body_particle.twist.Last().linear.y<<std::endl;}
    if(test_number==41){
        ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interp;
        for(int i=0;i<sample_points.m;i++){
            TV X=sample_points(i),V;
            for(int d=0;d<V.m;d++)
            V(d)=interp.Clamped_To_Array(fluids_parameters.grid->Get_Face_Grid(d),face_velocities.Component(d),X);
            LOG::cout<<"velocity at "<<X<<" : "<<V<<std::endl;}}
    if(test_number==43){
        ARRAY<T,FACE_INDEX<2> > face_velocities_copy(fluid_collection.incompressible_fluid_collection.face_velocities);
        DEFORMABLE_PARTICLES<TV>& particles=dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(test_43_triangulated_area->particles);
        RANGE<TV_INT> domain_indices=fluids_parameters.grid->Domain_Indices();
        T volume=fluids_parameters.grid->Cell_Size();
        double kinetic_energy=0;
        for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            int simplex=test_43_triangulated_area->Inside(iterator.Location(),1e-4);
            FACE_INDEX<2> face_index=iterator.Full_Index();
            if(simplex!=0){
                TRIANGLE_2D<T> triangle=test_43_triangulated_area->Get_Element(simplex);
                VECTOR<T,3> barycentric_weights=TRIANGLE_2D<T>::Barycentric_Coordinates(iterator.Location(),triangle.X);
                VECTOR<int,3> triangle_indices=test_43_triangulated_area->mesh.elements(simplex);
                T scalar_velocity=0;
                for(int i=0;i<3;i++)
                    scalar_velocity+=barycentric_weights(i)*particles.V(triangle_indices(i))(face_index.axis);
                face_velocities_copy(face_index)=scalar_velocity;
            }
            if(face_index.index(face_index.axis)<=domain_indices.max_corner(face_index.axis)) // avoid double counting wrapped faces
                kinetic_energy+=sqr(face_velocities_copy(face_index));
        }
        kinetic_energy*=.5*volume;
        LOG::cout<<"Total kinetic energy at frame "<<frame<<" is "<<kinetic_energy<<std::endl;

        T potential_energy=0,solid_kinetic_energy=0;
        solid_body_collection.Compute_Energy(frame/frame_rate,solid_kinetic_energy,potential_energy);
        LOG::cout<<"Total potential energy at frame "<<frame<<" is "<<potential_energy<<std::endl;

        ARRAY<T,TV_INT> stream_function(RANGE<TV_INT>(domain_indices.min_corner,domain_indices.max_corner+1));
        stream_function(1,1)=0;
        std::ofstream out(STRING_UTILITIES::string_sprintf("%s/%i/stream_function.dat",output_directory.c_str(),frame).c_str());
        TV dX=fluids_parameters.grid->DX();
        // use nodes so the derivatives are central
        for(int i=0;i<domain_indices.max_corner.x+1;i++){
            for(int j=0;j<domain_indices.max_corner.y+1;j++){
                int count=0;
                T estimate=0;
                if(i>1){ // integrate v_y in the x direction
                    estimate+=stream_function(i-1,j)+dX.x*face_velocities_copy.Component(1)(i-1,j);
                    count++;}
                if(j>1){ // integrate v_x in the y direction
                    estimate+=stream_function(i,j-1)-dX.y*face_velocities_copy.Component(0)(i,j-1);
                    count++;}
                if(count)
                    stream_function(i,j)=estimate/count;
                out<<i<<" "<<j<<" "<<stream_function(i,j)<<std::endl;
            }
            out<<std::endl;
        }
    }
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==40 && solid_body_collection.rigid_body_collection.rigid_body_particle.frame.m>=3) solid_body_collection.rigid_body_collection.rigid_body_particle.frame(2).t=TV(.5,.5)*scale_length;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<constrained_node_positions.m;i++){
        PAIR<int,TV>& node_pair=constrained_node_positions(i);
        V(node_pair.x)=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<deformable_object_enslaved_nodes.m;i++) V(deformable_object_enslaved_nodes(i))=TV();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<constrained_node_positions.m;i++){
        PAIR<int,TV>& node_pair=constrained_node_positions(i);
        X(node_pair.x)=node_pair.y;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
    //if(test_number==43) // use semi-Lagrangian and rasterize structure velocity field onto grid
    //    fluids_parameters.incompressible->Set_Custom_Advection(fluids_parameters.semi_lagrangian);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(fluids_parameters.use_slip){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).use_viscous_forces=use_viscous_forces;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    if(test_number==35)
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(0).Fill(velocity_multiplier);
    
    if(test_number==31 || test_number==37)
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(1).Fill((T).2);

    if(test_number==38)
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(0).Fill((T)1);

    if(fluids_parameters.use_coupled_implicit_viscosity){
        fluid_collection.incompressible_fluid_collection.viscosity.Resize(fluids_parameters.grid->Domain_Indices(1));
        fluid_collection.incompressible_fluid_collection.viscosity.Fill(fluids_parameters.viscosity);
        fluids_parameters.implicit_viscosity=false;}

    if(test_number==42){
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(123);
        ARRAY_VIEW<T> av(fluid_collection.incompressible_fluid_collection.face_velocities.buffer_size,fluid_collection.incompressible_fluid_collection.face_velocities.base_pointer);
        random.Fill_Uniform(av,-1,1);}

    if(test_number==43){
        // fluid velocities
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
            fluid_collection.incompressible_fluid_collection.face_velocities(face_index)=Oscillating_Disk_Domain_Velocity_Sample(iterator.Location())(face_index.axis);
        }

        // structure velocities
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        for(int i=0;i<particles.Size();i++)
            particles.V(i)=Oscillating_Disk_Domain_Velocity_Sample(particles.X(i));
    }
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(test_number>=9) return;
    BASE::Adjust_Density_And_Temperature_With_Sources(RANGE<TV>(TV((T).45,(T)0),TV((T).55,(T).1)),MATRIX<T,3>::Identity_Matrix(),1,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    BASE::Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==10 || test_number==11){
        // set left wall velocities
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=velocity_multiplier;
            psi_N(axis,face_index)=true;}}
    else if(test_number==12){
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=velocity_multiplier;
            psi_N(axis,face_index)=true;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=0;
            psi_N(axis,face_index)=true;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=0;
            psi_N(axis,face_index)=true;}}
    else if(test_number==13){
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(iterator.Location().y>(0.5-sin(velocity_angle)))
                face_velocities(axis,face_index)=velocity_multiplier*cos(velocity_angle);
            else
                face_velocities(axis,face_index)=-velocity_multiplier*cos(velocity_angle);
            if(iterator.Location().y>(0.5-sin(velocity_angle)) || iterator.Location().y<(0.5-sin(velocity_angle)-fluids_parameters.grid->DX()(1)))
                psi_N(axis,face_index)=true;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=-velocity_multiplier*sin(velocity_angle);
            psi_N(axis,face_index)=true;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,1);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(iterator.Location().y<(0.5+sin(velocity_angle)))
                face_velocities(axis,face_index)=-velocity_multiplier*cos(velocity_angle);
            else
                face_velocities(axis,face_index)=velocity_multiplier*cos(velocity_angle);
            if(iterator.Location().y<(0.5+sin(velocity_angle)) || iterator.Location().y>(0.5+sin(velocity_angle)+fluids_parameters.grid->DX()(1)))
                psi_N(axis,face_index)=true;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=velocity_multiplier*sin(velocity_angle);
            psi_N(axis,face_index)=true;}}
    else if(test_number==33){
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            face_velocities(iterator.Full_Index())=velocity_multiplier;psi_N(iterator.Full_Index())=true;}}
    else if(test_number==35 || test_number==36)
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            face_velocities(iterator.Full_Index())=velocity_multiplier;psi_N(iterator.Full_Index())=true;}
    else if(test_number==41){
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            face_velocities(iterator.Full_Index())=1;psi_N(iterator.Full_Index())=true;}}
    else if(test_number==44){
        for(FACE_ITERATOR iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()){
            face_velocities(iterator.Full_Index())=1;psi_N(iterator.Full_Index())=true;}}
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,FACE_INDEX<TV::m> >& outside) PHYSBAM_OVERRIDE
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*fluids_parameters.grid); it.Valid(); it.Next()){
        COLLISION_GEOMETRY_ID body_id;
        if(collision_bodies_affecting_fluid.Inside_Any_Body(it.Location(),body_id)) outside(it.Full_Index())=true;}
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,TV_INT>& outside) PHYSBAM_OVERRIDE
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*fluids_parameters.grid); it.Valid(); it.Next()){
        COLLISION_GEOMETRY_ID body_id;
        if(collision_bodies_affecting_fluid.Inside_Any_Body(it.Location(),body_id)) outside(it.index)=true;}
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE Get_Boundary_Along_Ray(const FACE_INDEX<TV::m>& f1,const FACE_INDEX<TV::m>& f2,T& theta,T& value) PHYSBAM_OVERRIDE
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused;
    COLLISION_GEOMETRY_ID body_id;
    TV X1=fluids_parameters.grid->Axis_X_Face(f1);
    TV X2=fluids_parameters.grid->Axis_X_Face(f2);
    T dx=(X1-X2).Magnitude();
    RAY<TV> ray(SEGMENT_2D<T>(X1,X2));
    bool found=collision_bodies_affecting_fluid.Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
    if(found){
        T th=ray.t_max/dx;
        if(type==BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused || th<theta){
            theta=th;
            value=collision_bodies_affecting_fluid.Object_Velocity(body_id,ray.aggregate_id,ray.Point(ray.t_max))(f1.axis);
            type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::noslip;}}

    VECTOR<RANGE<TV_INT>,TV::dimension> fi=fluids_parameters.grid->Face_Indices(0);
    for(int i=0;i<TV::m;i++){
        if(f2.index(i)<=fi(f2.axis).min_corner(i)-(f2.axis!=i)){
            T th=abs(fluids_parameters.grid->domain.min_corner(i)-X1(i))/dx;
            if(type==BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused || th<theta){
                theta=th;
                value=0;
                type=fluids_parameters.domain_walls[i][0]?BOUNDARY_CONDITIONS_CALLBACKS<TV>::noslip:BOUNDARY_CONDITIONS_CALLBACKS<TV>::free;}}
        else if(f2.index(i)>=fi(f2.axis).max_corner(i)+(f2.axis!=i)){
            T th=abs(fluids_parameters.grid->domain.max_corner(i)-X1(i))/dx;
            if(type==BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused || th<theta){
                theta=th;
                value=0;
                type=fluids_parameters.domain_walls[i][1]?BOUNDARY_CONDITIONS_CALLBACKS<TV>::noslip:BOUNDARY_CONDITIONS_CALLBACKS<TV>::free;}}}

    static VECTOR<T,3> color_map[]={VECTOR<T,3>(1,0,0),VECTOR<T,3>(1,.5,0),VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,.5,0),VECTOR<T,3>(0,1,1),VECTOR<T,3>(1,1,0)};

    if(ARRAY_VIEW<VECTOR<T,3> >* color_attribute=debug_particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR)){
        int p=debug_particles.Add_Element();
        debug_particles.X(p)=X1+theta*(X2-X1);
        (*color_attribute)(p)=color_map[type];}

    return type;
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE Get_Boundary_Along_Ray(const TV_INT& c1,const TV_INT& c2,T& theta,T& value) PHYSBAM_OVERRIDE
{
    PHYSBAM_FATAL_ERROR();
    return BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused;
}
//#####################################################################
// Function Balloon
//#####################################################################
void Balloon()
{
    TV balloon_initial_position=TV((T).5,light_sphere_initial_height);
    int vertices=50;
    deformable_object_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Circle_Deformable_Object(solid_body_collection.deformable_body_collection,vertices,balloon_initial_position,balloon_initial_radius);
    SEGMENTED_CURVE_2D<T>& segmented_curve=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE_2D<T>&>(deformable_object_id);
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,(T)2,true);

    T cut_sphere_radius=(T).1;
    T seperation_distance=balloon_initial_radius+cut_sphere_radius-(T).05;

    SPHERE<TV> analytic_cut_sphere(balloon_initial_position-TV((T)0,seperation_distance),cut_sphere_radius);
    T a=(sqr(seperation_distance)+sqr(balloon_initial_radius)-sqr(cut_sphere_radius))/(2*seperation_distance);
    T height=light_sphere_initial_height-balloon_initial_radius;
    T dist=sqrt(sqr(balloon_initial_radius)-sqr(a));

    ARRAY<int> deletion_list; // List of deleted segments
    ARRAY<bool> is_constrained;
    is_constrained.Resize(segmented_curve.particles.Size());
    is_constrained.Fill(false);

    for(int i=0;i<segmented_curve.mesh.elements.m;i++){
        const SEGMENT_2D<T>& segment=segmented_curve.Get_Element(i);
        int node1,node2;segmented_curve.mesh.elements(i).Get(node1,node2);
        if(analytic_cut_sphere.Lazy_Inside(segment.X.x) || analytic_cut_sphere.Lazy_Inside(segment.X.y)){
            if(analytic_cut_sphere.Lazy_Inside(segment.X.x) && analytic_cut_sphere.Lazy_Inside(segment.X.y))
                deletion_list.Append(i);
            else{
                if(analytic_cut_sphere.Lazy_Inside(segment.X.x)) is_constrained(node1)=true;
                if(analytic_cut_sphere.Lazy_Inside(segment.X.y)) is_constrained(node2)=true;}}}

    segmented_curve.mesh.Delete_Elements(deletion_list);
    ARRAY<int> condensation_mapping;
    segmented_curve.Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);

    // rotate all the points around (.5,0,.5) at a 45 degree angle
    FRAME<TV> frame;
    if(test_number==8) frame=FRAME<TV>(TV(),ROTATION<TV>::From_Angle((T)0));
    else PHYSBAM_NOT_IMPLEMENTED();
    for(int i=0;i<is_constrained.m;i++){
        if(is_constrained(i)){
            TV& position=segmented_curve.particles.X(condensation_mapping(i));
            position=TV(balloon_initial_position.x+(position.x>balloon_initial_position.x?dist:-dist),height);
            position=frame*position;
            constrained_node_positions.Append(PAIR<int,TV>(condensation_mapping(i),position));}
        else if(condensation_mapping(i)){
            TV& position=segmented_curve.particles.X(condensation_mapping(i));
            position=frame*position;}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1:{
//            deformable_circle_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Circle_Deformable_Object(deformable_object,10,TV((T).5,(T)1),(T).25);
//            SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(deformable_object.structures(deformable_circle_id));assert(segmented_curve);
//            segmented_curve->Set_Mass_Of_Particles();
//            /*TRIANGULATED_AREA<T>& triangulated_area=*/solids_tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            /*TRIANGULATED_AREA<T>& triangulated_area=*/
            solids_tests.Create_Mattress(mattress_grid,true,0,2);
            break;}
        case 10:
        case 2:{
            //RIGID_BODY<TV>& rigid_body_square=solids_tests.Add_Rigid_Body("square",(T).1,(T)0);
            T radius=.3;
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",(T)radius,(T)0);
            rigid_body.thin_shell=false;
            rigid_body.Frame()=FRAME<TV>(TV((T).5,(T).5));
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="circle";
            T density=100;
            rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);
            if(test_number==10){
                rigid_body.is_static=true;
                solids_tests.Create_Mattress(mattress_grid_stub,true,0,400);}
            break;}
        case 11:{
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("square",(T).1,(T)0);
            rigid_body.Frame()=FRAME<TV>(TV((T).5,(T).5));
            rigid_body.Frame().r=ROTATION<TV>::From_Angle((T)pi/7);
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="square";
            T density=2000;
            rigid_body.Set_Mass((T)pi*(T).01*(T)density*mass_multiplier);
            rigid_body.is_static=true;
            solids_tests.Create_Mattress(mattress_grid_stub,true,0,400);
            break;}
        case 4:{
             deformable_object_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Circle_Deformable_Object(deformable_body_collection,10,TV((T).5,(T).7),(T).25);
             SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.deformable_geometry.structures(deformable_object_id));
             SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,(T)2000);
             break;}
        case 5:{
            deformable_object_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Deformable_Object(deformable_body_collection,100,TV((T).2,(T).5),TV((T).8,(T).5));
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.deformable_geometry.structures(deformable_object_id));
            THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Set_Mass(segmented_curve,50);
            break;}
        case 6:{
            RIGID_BODY<TV>& rigid_body_cup=solids_tests.Add_Rigid_Body("cup",(T).1,(T)0);
            rigid_body_cup.Frame()=FRAME<TV>(TV((T).5,(T).3));
            rigid_body_cup.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_cup.name="cup";
            rigid_body_cup.Set_Mass(50);
            rigid_body_cup.thin_shell=true;
            break;}
        case 7:{
            for(int i=0;i<rigid_body_count;i++){
                RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",(T).05,(T)0);
                rigid_body.Frame()=FRAME<TV>(TV((T).1+i*(T).1,(T).2+i*(T).1));
                rigid_body.Set_Coefficient_Of_Restitution((T)0);
                rigid_body.name="circle";
                rigid_body.Set_Mass(10);}
            break;}
        case 8:{
            Balloon();
            break;}
        case 9:{
            T solid_density=fluids_parameters.density*mass_multiplier;
            T sphere_mass=sqr((T).4)*(T)pi*solid_density;
            RIGID_BODY<TV>& rigid_body_sphere=solids_tests.Add_Rigid_Body("circle",(T).4,(T)0);
            rigid_body_id=rigid_body_sphere.particle_index;
            rigid_body_sphere.Update_Bounding_Box();
            rigid_body_sphere.Frame().t=TV((T)0,(T)1);
            rigid_body_sphere.Set_Mass(sphere_mass);
            rigid_body_sphere.Set_Coefficient_Of_Restitution((T)1);
            rigid_body_sphere.coefficient_of_friction=(T)1;
            break;}
        case 12:
            {
            // RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("ground",(T).004,(T)0,false,true);
            // rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            // rigid_body.Frame().r=ROTATION<TV>::From_Angle(rotation_angle);
            // rigid_body.Set_Coefficient_Of_Restitution((T)0);
            // rigid_body.name="square");
            // T density=2000;
            // rigid_body.Set_Mass((T)pi*(T).01*(T)density*mass_multiplier);
            // rigid_body.thin_shell=true;
            // rigid_body.is_static=true;
            
            T scale=(T).15;
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",scale,(T)0,true,true);
            rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            // rigid_body.Frame()=FRAME<TV>(TV((T)5.0,(T)0.5));
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="circle_fixed";
            T density=1000;
            rigid_body.Set_Mass(scale*scale*(T)pi*(T).01*(T)density*mass_multiplier);
            // rigid_body.thin_shell=false;
            rigid_body.is_static=true;
            
            RIGID_BODY<TV>& rigid_body_dummy=solids_tests.Add_Rigid_Body("circle",(T).05,(T)0);
            rigid_body_dummy.Frame()=FRAME<TV>(TV((T)10,(T)0));
            rigid_body_dummy.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_dummy.name="circle";
            rigid_body_dummy.Set_Mass(10);
            }
            break;
        case 13:
            {
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("ground",(T).012,(T)0,false,true);
            rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            rigid_body.Frame().r=ROTATION<TV>::From_Angle(velocity_angle);
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="square";
            T density=2000;
            rigid_body.Set_Mass((T)pi*(T).01*(T)density*mass_multiplier);
            rigid_body.thin_shell=true;
            rigid_body.is_static=true;
            
            RIGID_BODY<TV>& rigid_body_dummy=solids_tests.Add_Rigid_Body("circle",(T).05,(T)0);
            rigid_body_dummy.Frame()=FRAME<TV>(TV((T)10,(T)0));
            rigid_body_dummy.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_dummy.name="circle";
            rigid_body_dummy.Set_Mass(10);
            }
            break;
        case 20:{
            //RIGID_BODY<TV>& rigid_body_square=solids_tests.Add_Rigid_Body("square",(T).1,(T)0);
            T radius=.4;
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",radius,(T)0);
            rigid_body.thin_shell=false;
            rigid_body.Frame()=FRAME<TV>(TV((T)1,(T).15));
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="circle";
            T density=100;
            rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);
            rigid_body.Twist().linear.y=5;
            break;}
        case 30:{
            T radius=(T).3;
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",radius,(T)0);
            rigid_body.thin_shell=false;
            rigid_body.Frame()=FRAME<TV>(TV((T).5,(T)1.5));
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="circle";
            T density=1000;
            rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);
            break;}
        case 31: Uniform_Flow_Test(); break;
        case 32: Falling_Rigid_Circle_Test(); break;
        case 33: Flexible_Beam_Test();break;
        case 34: Vibrating_Circle();break;
        case 35: Flexible_Filament_Test();break;
        case 36: Simple_Fluid_Test();break;
        case 37: Deformable_Uniform_Flow_Test();break;
        case 38: Coupled_Viscosity_Test();break;
        case 39: Refine_Circle();break;
        case 40: Analytic_Test();break;
        case 41: Flow_Past_Fixed_Cylinder();break;
        case 42: Sanity_Test_Stokes_No_Viscosity();break;
        case 43: Oscillating_Disk();break;
        case 44: Vortex_Shedding();break;
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
    
    T solid_gravity=(T)9.8;
    switch(test_number){
        case 1:{
//            SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(deformable_body_collection.deformable_geometry.structures(deformable_circle_id));assert(segmented_curve);
//            fluids_parameters.collision_bodies_affecting_fluid->Add_Body(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*segmented_curve));
//            solid_body_collection.Add_Force(Create_Edge_Springs(particles,rigid_body_collection,*segmented_curve));
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            TRIANGULATED_AREA<T>& triangulated_area=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            triangulated_area.Initialize_Hierarchy();
            DEFORMABLES_FORCES<TV>* spring=Create_Edge_Springs(triangulated_area,(T)7/(1+sqrt((T)2)),(T)10);
            solid_body_collection.Add_Force(spring);
            DEFORMABLES_FORCES<TV>* altitude_spring=Create_Altitude_Springs(triangulated_area,(T)7,(T)10);
            solid_body_collection.Add_Force(altitude_spring);
            
            //solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e3,(T).45,(T).01)));
            //DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area.Get_Boundary_Object());
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
            deformable_collisions.object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(deformable_collisions);
            //solids_tests.Add_Ground();
            break;}
        case 2:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            //Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-1),true,true);
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            solids_tests.Add_Ground();
            break;}
        case 10:
        case 11:{
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            TRIANGULATED_AREA<T>& triangulated_area=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            triangulated_area.Initialize_Hierarchy();
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e3,(T).45,(T).01)));
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area.Get_Boundary_Object());
            deformable_collisions.object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(deformable_collisions);
            break;}
        case 4:
        case 5:{
            solids_tests.Add_Ground();
            solid_gravity=(T).75;
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.deformable_geometry.structures(deformable_object_id));
            segmented_curve.Initialize_Hierarchy();
            fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(segmented_curve),0,true);
            solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve));
            break;}
        case 6:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            break;}
        case 7:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            for(int i=0;i<rigid_body_count;i++) Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-i+1));
            break;}
        case 8:{
            SEGMENTED_CURVE_2D<T>& segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE_2D<T>&>(deformable_object_id);
            segmented_curve.Initialize_Hierarchy();
            fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(segmented_curve),0,true);
            T linear_stiffness=stiffness_multiplier*(T).5/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,linear_stiffness,linear_damping)); // were *2 and *10
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Segment_Bending_Springs(segmented_curve,bending_stiffness,bending_damping));
            break;}
        case 9:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_id));
            break;}
        case 12:{
            // solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            // Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            break;}
        case 13:{
            // solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            // Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()-1));
            Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            
            // TRIANGULATED_AREA<T>& triangulated_area=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            // triangulated_area.Initialize_Hierarchy();
            // solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e3,(T).45,(T).01)));
            // DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area.Get_Boundary_Object());
            // deformable_collisions.object.Initialize_Hierarchy();
            // Add_To_Fluid_Simulation(deformable_collisions,true,false);
            break;}
        case 20:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            //Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()-1),true,true);
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            break;}
        case 30:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
            solids_tests.Add_Ground();
            break;}
        case 31: break;
        case 32: break;
        case 33: break;
        case 36: break;
        case 37: break;
        case 38: break;
        default:break;}

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
// Function Uniform_Flow_Test
//#####################################################################
void Uniform_Flow_Test()
{
//    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
//    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    last_frame=1000;
    (*fluids_parameters.grid).Initialize(TV_INT(resolution,(int)(2*resolution)),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    
    T radius=(T).3;
    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",radius,(T)0);
    rigid_body.thin_shell=false;
    rigid_body.Frame()=FRAME<TV>(TV((T).5,(T)1.5));
    rigid_body.Twist().linear.y=(T).2;
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.name="circle";
    T density=100;
    rigid_body.Set_Mass((T)pi*sqr(radius)*(T)density*mass_multiplier);

//    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,(T)9.8));
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
    solids_tests.Add_Ground();
}
//#####################################################################
// Function Deformable_Uniform_Flow_Test
//#####################################################################
void Deformable_Uniform_Flow_Test()
{
//    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
//    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    last_frame=1000;
    (*fluids_parameters.grid).Initialize(TV_INT(resolution,(int)(2*resolution)),RANGE<TV>(TV((T)-.5,(T)-1),TV((T).5,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    
    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Create_Triangulated_Object(data_directory+"/Triangulated_Areas/circle-216.tri2d",RIGID_BODY_STATE<TV>(),true,true,(T).25);

    for(int i=0;i<particles.V.m;i++)
        particles.V(i)=TV(0,.2);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    triangulated_area.Initialize_Hierarchy();
    DEFORMABLES_FORCES<TV>* spring=Create_Edge_Springs(triangulated_area,(T)10/(1+sqrt((T)2)),(T).1);
    solid_body_collection.Add_Force(spring);
    DEFORMABLES_FORCES<TV>* altitude_spring=Create_Altitude_Springs(triangulated_area,(T)10,(T).1,false);
    solid_body_collection.Add_Force(altitude_spring);

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
}
//#####################################################################
// Function Falling_Rigid_Circle_Test
//#####################################################################
void Falling_Rigid_Circle_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    fluids_parameters.use_coupled_implicit_viscosity=true;

    fluids_parameters.gravity=fluid_gravity;
    fluids_parameters.density=(T)1000/sqr(scale_length);
    fluids_parameters.domain_walls[0][0]=true;
    fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    
    T radius=(T).005*scale_length;
    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",radius,(T)0);
    rigid_body.thin_shell=false;
    rigid_body.Frame()=FRAME<TV>(TV((T).02,(T)0.14)*scale_length);
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.name="circle";
    T density=2000/sqr(scale_length);
    rigid_body.Mass()=(T)pi*sqr(radius)*(T)density;
    rigid_body.Inertia_Tensor()(1,1)=rigid_body.Mass()*sqr(radius)/2;

    for(int b=0;b<3;b++){
        Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
        rigid_body_collection.Rigid_Body(b).is_static=true;}

    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity*scale_length));
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));

    PHYSBAM_ASSERT(fluids_parameters.viscosity);
    LOG::cout<<"VISCOSITY "<<fluids_parameters.viscosity<<std::endl;
    // Xiaodong Wang and Wing Kam Liu, Extended immersed boundary method using FEM and RKPM.  Comput. Methods Appl. Mech. Engrg. 193 (2004) 1305–1321
    T aL=.25; // solid diameter divided by channel width
    T analytic_v=(density-fluids_parameters.density)*solid_gravity*scale_length*radius*radius/(4*fluids_parameters.viscosity)*(-log(aL)-(T)0.9157+(T)1.7244*sqr(aL)-1.7302*sqr(sqr(aL)));

    PHYSBAM_ASSERT(sizeof(T)==8);
    LOG::cout<<"Analytic solution for terminal velocity: "<<analytic_v<<std::endl;
}
//#####################################################################
// Function Flexible_Beam_Test
//#####################################################################
void Flexible_Beam_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    //RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(2*resolution),resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.use_body_force=false;
    velocity_multiplier=.25;

//    for(int b=0;b<2;b++){
//        if(b==1) Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
//        rigid_body_collection.Rigid_Body(b).is_static=true;}

    solid_body_collection.Set_CFL_Number(10);
    
    mattress_grid=GRID<TV>(TV_INT(5,20),RANGE<TV>(TV((T).95,(T)0),TV((T)1.05,(T).5)));
    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Create_Mattress(mattress_grid,true,0,200);
    for(int i=0;i<particles.Size();i++)
        if(particles.X(i).y<.0001)
            particles.mass(i)=FLT_MAX;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    //T solid_gravity=(T)9.8;
    //solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
        //deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
    triangulated_area.Initialize_Hierarchy();
    DEFORMABLES_FORCES<TV>* spring=Create_Edge_Springs(triangulated_area,(T)400/(1+sqrt((T)2)),(T)2);
    solid_body_collection.Add_Force(spring);
    DEFORMABLES_FORCES<TV>* altitude_spring=Create_Altitude_Springs(triangulated_area,(T)400,(T)2,false);
    solid_body_collection.Add_Force(altitude_spring);

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
}
//#####################################################################
// Function Vibrating_Circle
//#####################################################################
void Vibrating_Circle()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(2*resolution),resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)1;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    velocity_multiplier=.5;

    solid_body_collection.Set_CFL_Number(10);

    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Create_Triangulated_Object(data_directory+"/Triangulated_Areas/circle-216.tri2d",RIGID_BODY_STATE<TV>(),true,true,(T).25);

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    triangulated_area.Initialize_Hierarchy();
    DEFORMABLES_FORCES<TV>* spring=Create_Edge_Springs(triangulated_area,(T)1,(T).1);
    solid_body_collection.Add_Force(spring);
    DEFORMABLES_FORCES<TV>* altitude_spring=Create_Altitude_Springs(triangulated_area,(T)1,(T).1,false);
    solid_body_collection.Add_Force(altitude_spring);

    MATRIX<T,2> mat((T)1.2,0,0,1/(T)1.2);
    TV shift(1,(T).5);
    for(int i=0;i<particles.X.m;i++) particles.X(i)=mat*particles.X(i)+shift;

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
}
//#####################################################################
// Function Refine_Circle
//#####################################################################
void Refine_Circle()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(2*resolution),resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)1;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    velocity_multiplier=.5;

    solid_body_collection.Set_CFL_Number(10);

    PHYSBAM_ASSERT(circle_refinement);
    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Copy_And_Add_Structure(*TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),(T).25),circle_refinement));
    PHYSBAM_ASSERT(triangulated_area.mesh.Orientations_Consistent());
    particles.mass.Fill((T)100/particles.mass.m);

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e3,(T).45,(T).01)));

    MATRIX<T,2> mat((T)1.2,0,0,1/(T)1.2);
    for(int i=0;i<particles.X.m;i++) particles.X(i)=mat*particles.X(i)+TV(1,(T).5);

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
}
//#####################################################################
// Function Refine_Circle
//#####################################################################
void Analytic_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;

    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.gravity=(T)9.8*scale_length;
    fluids_parameters.density=(T)100/(scale_length*scale_length);
    fluids_parameters.domain_walls[0][0]=true;
    fluids_parameters.domain_walls[0][1]=true;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.use_levelset_viscosity=true;
    velocity_multiplier=.5;
    PHYSBAM_ASSERT(fluids_parameters.use_slip);
    solid_body_collection.Set_CFL_Number(10);

    for(int b=0;b<2;b++){
        Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
        rigid_body_collection.Rigid_Body(b).is_static=true;}

    if(use_solid){
        RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Box(scale_length*TV(solid_width,1.2),5);
        rigid_body.Frame()=FRAME<TV>(scale_length*TV((T).5,.5));
        rigid_body.Set_Coefficient_Of_Restitution((T)0);
        rigid_body.Set_Mass(solid_density*solid_width);

        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity*scale_length));
        Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);}

    T solid_mass=use_solid?solid_body_collection.rigid_body_collection.rigid_body_particle.mass(2):0;
    T rho=fluids_parameters.density;
    TV size=fluids_parameters.grid->domain.Edge_Lengths();
    size.x=(size.x-solid_width)/2;
    analytic_solution=-(solid_mass*solid_gravity+rho*size.x*size.y*fluids_parameters.gravity)*size.x/(2*size.y*fluids_parameters.viscosity);
    LOG::cout<<"analytic_solution "<<analytic_solution<<std::endl;

    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),0));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);
}
//#####################################################################
// Function Flow_Past_Fixed_Cylinder
//#####################################################################
void Flow_Past_Fixed_Cylinder()
{
//    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
//    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
//    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    fluids_parameters.use_coupled_implicit_viscosity=true;

    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.gravity=0;
    fluids_parameters.density=1;
    fluids_parameters.viscosity=(T).1;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.use_levelset_viscosity=true;
//    PHYSBAM_ASSERT(fluids_parameters.use_slip);
    solid_body_collection.Set_CFL_Number(10);

//    for(int b=0;b<2;b++){
//        Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
//        rigid_body_collection.Rigid_Body(b).is_static=true;}

    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Sphere((T).5,1,7);
    rigid_body.Frame()=FRAME<TV>(TV(2,2));
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.is_static=true;
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    sample_points.Append(TV(1.25,2.5));
    sample_points.Append(TV(3,2.25));
    sample_points.Append(TV(2,3));

    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),0));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);
}
//#####################################################################
// Function Flow_Past_Fixed_Cylinder
//#####################################################################
void Vortex_Shedding()
{
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    fluids_parameters.use_coupled_implicit_viscosity=true;

    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.gravity=0;
    fluids_parameters.density=1;
    fluids_parameters.viscosity=(T).01;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.use_levelset_viscosity=true;
    solid_body_collection.Set_CFL_Number(10);

    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Sphere((T).5,1,7);
    rigid_body.Frame()=FRAME<TV>();
    rigid_body.Set_Coefficient_Of_Restitution((T)0);
    rigid_body.is_static=true;
    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body);

    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),0));
    FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);
}
//#####################################################################
// Function Oscillating_Disk_Domain_Velocity_Sample
//#####################################################################
TV Oscillating_Disk_Domain_Velocity_Sample(TV location)
{
    // stream function is psi = psi_0 sin(k_x x) sin (k_y y), psi_0 = 5e-2, k_x=k_y=2pi
    // v_x=psi_0 k_y sin(k_x x) cos(k_y y)
    // v_y=-psi_0 k_x cos(k_x x) sin(k_y y)
    T k_x=two_pi,k_y=two_pi;
    T psi_0=5e-2;
    return TV(psi_0*k_y*sin(k_x*location.x)*cos(k_y*location.y),-psi_0*k_x*cos(k_x*location.x)*sin(k_y*location.y));
}
//#####################################################################
// Function Oscillating_Disk
//#####################################################################
void Oscillating_Disk()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    frame_rate=1000;
    last_frame=(int)frame_rate; // one second;

    (*fluids_parameters.grid).Initialize(TV_INT(resolution,resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)1;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;
    fluids_parameters.domain_walls[1][1]=false;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.cfl=(T).9;

    solid_body_collection.Set_CFL_Number(10);

    PHYSBAM_ASSERT(circle_refinement);
    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Copy_And_Add_Structure(*TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(.5,.5),(T).2),circle_refinement));
    solids_tests.Set_Mass_Of_Particles(triangulated_area,1,false);
    PHYSBAM_ASSERT(triangulated_area.mesh.Orientations_Consistent());

    /*std::string circle_file=data_directory+STRING_UTILITIES::string_sprintf("/Triangulated_Areas/circle-%d.tri2d",solid_resolution);
      TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Create_Triangulated_Object(circle_file,RIGID_BODY_STATE<TV>(),true,true,(T).2);*/

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    //TV shift((T).5,(T).5);
    //for(int i=0;i<particles.X.m;i++){
    //    //particles.X(i)=/*mat*/particles.X(i)+shift;
    //    particles.mass(i)*=.01; // create_triangulated_object sets the density to 100 by default in 2d
    //}
    triangulated_area.Initialize_Hierarchy();
    T poisson_ratio=.475;
    T shear_modulus=1;
    T youngs_modulus=2*(1+poisson_ratio)*shear_modulus;
    finite_volume=Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>(youngs_modulus,poisson_ratio));
    solid_body_collection.Add_Force(finite_volume);

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
    test_43_triangulated_area=&triangulated_area;

    fluids_parameters.viscosity=1e-3;
    fluids_parameters.use_levelset_viscosity=true;
    fluids_parameters.second_order_cut_cell_method=true;

    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    fluids_parameters.fluid_boundary=new BOUNDARY_MAC_GRID_PERIODIC<GRID<TV>,T>;
    fluids_parameters.periodic_boundary.Fill(true);
}
//#####################################################################
// Function Flexible_Filament_Test
//#####################################################################
void Flexible_Filament_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    //RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(3*resolution)+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)3,(T)1)));
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    velocity_multiplier=5;

    solid_body_collection.Set_CFL_Number(10);
    
    GRID<TV> filament_grid(TV_INT(20,5),RANGE<TV>(TV((T).5,(T).475),TV((T)1.25,(T).525)));
    TRIANGULATED_AREA<T>& triangulated_area=solids_tests.Create_Mattress(filament_grid,true,0,200);
    //GRID<VECTOR<T,1> > filament_grid(10,RANGE<VECTOR<T,1> >(VECTOR<T,1>(.5),VECTOR<T,1>(1)));
    //RIGID_GEOMETRY_STATE<TV> state;
    //state.frame.t.x=1;
    //state.frame.t.y=.5;
    //SEGMENTED_CURVE_2D<T>& segmented_curve=solids_tests.Create_Segmented_Curve(filament_grid,state,200);
    for(int i=0;i<particles.Size();i++)
        if(particles.X(i).x<.501)
            particles.mass(i)=FLT_MAX;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    //T solid_gravity=(T)9.8;
    //solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
        //deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
    triangulated_area.Initialize_Hierarchy();
    DEFORMABLES_FORCES<TV>* spring=Create_Edge_Springs(triangulated_area,(T)100/(1+sqrt((T)2)),(T)2);
    spring->compute_half_forces=true;
    solid_body_collection.Add_Force(spring);
    DEFORMABLES_FORCES<TV>* altitude_spring=Create_Altitude_Springs(triangulated_area,(T)100,(T)2,false);
    altitude_spring->compute_half_forces=true;
    solid_body_collection.Add_Force(altitude_spring);

    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area);//.Get_Boundary_Object());
    deformable_collisions.object.Initialize_Hierarchy();
    Add_To_Fluid_Simulation(deformable_collisions);
}
//#####################################################################
// Function Simple_Fluid_Test
//#####################################################################
void Simple_Fluid_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(2*resolution),resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    velocity_multiplier=.25;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Simple_Fluid_Test
//#####################################################################
void Coupled_Viscosity_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    (*fluids_parameters.grid).Initialize(TV_INT((int)(2*resolution),resolution),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)),true);
    fluids_parameters.gravity=(T)0;
    fluids_parameters.density=(T)100;
    fluids_parameters.domain_walls[0][0]=false;
    fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.use_coupled_implicit_viscosity=true;
    velocity_multiplier=.25;
    solids_parameters.use_post_cg_constraints=true;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Sanity_Test_Stokes_No_Viscosity
//#####################################################################
void Sanity_Test_Stokes_No_Viscosity()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->use_collision_face_neighbors=true;
    fluids_parameters.use_coupled_implicit_viscosity=true;

    fluids_parameters.gravity=0;
    fluids_parameters.density=(T)1000/sqr(scale_length);
    fluids_parameters.use_density=fluids_parameters.use_temperature=false;
    fluids_parameters.viscosity=0;
    fluids_parameters.stokes_flow=true;

    for(int b=0;b<3;b++){
        Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(b));
        rigid_body_collection.Rigid_Body(b).is_static=true;}

    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Sphere((T).3,fluids_parameters.density*2,5);
    rigid_body.Frame()=FRAME<TV>(TV((T).5,(T).5));
    rigid_body.name="circle";
    rigid_body.is_static=true;

    Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particle.Size()));
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(max_dt && dt>max_dt) dt=max_dt;
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(max_dt && dt>max_dt) dt=max_dt;
}
//#####################################################################
};
}
#endif
