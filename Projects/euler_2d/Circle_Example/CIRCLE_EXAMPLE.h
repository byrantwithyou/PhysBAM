//#####################################################################
// Copyright 2007-2009 Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CIRCLE_EXAMPLE
//#####################################################################
#ifndef __CIRCLE_EXAMPLE__
#define __CIRCLE_EXAMPLE__

#include <fstream>
#include <iostream>
#include "math.h"

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_ATTENUATION.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>

namespace PhysBAM{

template<class T_input>
class CIRCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
public:
    typedef T_input T;typedef VECTOR<T,2> TV;;typedef GRID<TV> T_GRID;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_FLUID_COLLISION_GEOMETRY_LIST;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::solids_fluids_parameters;
    using BASE::stream_type;using BASE::data_directory;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::solids_evolution;using BASE::Add_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> tests;

    TV_DIMENSION state_inside,state_outside; // // (density,velocity_x,velocity_y,pressure)
    TV shock_center;
    T shock_radius;
    T solid_mass;
    T e_min_for_clamping;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int sphere;
    INTERPOLATION_CURVE<T,TV> motion_curve;
    T solid_gravity;
    bool use_solids_gravity;
    bool transition_to_incompressible;
    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >* eos_smooth_transition;
    bool incompressible;
    bool use_soot;

    bool use_fixed_farfield_boundary;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;

    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool strong_shock;
    bool timesplit;
    bool use_slip;
    bool bucket_walls;
    bool all_walls;
    bool no_walls;
    bool exact;
    bool use_incompressible_gravity;
    T time_start_transition;
    T time_end_transition;
    T one_over_c_incompressible;

    /***************
    example explanation:
    1. Circular higher density region.
    2. Same as 1 + a static circular object.
    3. Same as 1 + a moving circular object.
    4. Same as 1 + a two-way coupled circular object.
    5. Same as 4 + walls on all sides.
    6. Same as 4 + wall on 3 sides.
    7. Static solid circle in a static uniform fluid receives an impulse at t=1
    8. Same as 4, but with a deformable circle
    9. Obrien's shock hitting a wall
    10. Trinity Test
    ***************/

    CIRCLE_EXAMPLE(const STREAM_TYPE stream_type,const bool incompressible_input)
        :BASE(stream_type,0,incompressible_input?fluids_parameters.SMOKE:fluids_parameters.COMPRESSIBLE),
        tests(*this,solid_body_collection),solid_mass((T).0625),rigid_body_collection(solid_body_collection.rigid_body_collection),
        use_solids_gravity(false),transition_to_incompressible(false),
        incompressible(incompressible_input),use_soot(true),use_fixed_farfield_boundary(false),
        run_self_tests(false),print_poisson_matrix(false),print_index_map(false),print_matrix(false),
        print_rhs(false),print_each_matrix(false),output_iterators(false),eno_scheme(2),eno_order(2),
        rk_order(3),cfl_number((T).6),strong_shock(false),timesplit(false),use_slip(false),bucket_walls(false),
        all_walls(false),no_walls(false),exact(false),use_incompressible_gravity(false),time_start_transition((T).5),
        time_end_transition((T).7),one_over_c_incompressible(0)
    {
    }

    virtual ~CIRCLE_EXAMPLE() {}

    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
    parse_args->Add("-eno_order",&eno_order,"eno_order","eno order");
    parse_args->Add("-rk_order",&rk_order,"rk_order","runge kutta order");
    parse_args->Add("-cfl",&cfl_number,"CFL","cfl number");
    parse_args->Add("-strong_shock",&strong_shock,"Use stronger shock with temperature ratio 2900:290, p ratio 1000:1");
    parse_args->Add("-mass",&solid_mass,"solid_mass","the mass of the solid in the simulation");
    parse_args->Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add("-slip",&use_slip,"use slip/spd for coupling");
    parse_args->Add("-bucket_walls",&bucket_walls,"Add walls on bottom, left and right sides");
    parse_args->Add("-all_walls",&all_walls,"Add walls on all sides");
    parse_args->Add("-no_walls",&no_walls,"No walls on all sides");
    parse_args->Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
    parse_args->Add("-transition_to_incompressible",&transition_to_incompressible,"transition to incompressible in a time window");
    parse_args->Add("-time_start_transition",&time_start_transition,"time","time to start transitioning to incompressible flow");
    parse_args->Add("-time_end_transition",&time_end_transition,"time","time to end transitioning to incompressible flow");
    parse_args->Add("-one_over_c_incompressible",&one_over_c_incompressible,"value","one over incompressible sound speed");
    parse_args->Add_Not("-no_soot",&use_soot,"don't advect soot");

    parse_args->Add("-use_fixed_farfield_boundary",&use_fixed_farfield_boundary,"use fixed farfield values for outflow boundaries");

    parse_args->Add("-use_incompressible_gravity",&use_incompressible_gravity,"add gravity on incompressible fluid");
    parse_args->Add("-use_solids_gravity",&use_solids_gravity,"add gravity on solids");

    parse_args->Add("-test_system",&run_self_tests,"run self tests");
    parse_args->Add("-print_poisson_matrix",&print_poisson_matrix,"print poisson matrix");
    parse_args->Add("-print_index_map",&print_index_map,"print index map");
    parse_args->Add("-print_matrix",&print_matrix,"print matrix");
    parse_args->Add("-print_rhs",&print_rhs,"print rhs");
    parse_args->Add("-print_each_matrix",&print_each_matrix,"print each matrix");
    parse_args->Add("-output_iterators",&output_iterators,"output iterators");
    parse_args->Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
    parse_args->Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"use preconditioner");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();

    timesplit=timesplit && !exact;
    //grid
    int cells=resolution;
    T grid_size=(T)1.;
    if(test_number==10) fluids_parameters.grid->Initialize(TV_INT(2,1)*cells+1,RANGE<VECTOR<T,2> >(VECTOR<T,2>((T)-100,(T)0),VECTOR<T,2>((T)100,(T)100))*grid_size);
    else fluids_parameters.grid->Initialize(TV_INT()+cells+1,RANGE<VECTOR<T,2> >(VECTOR<T,2>((T)-1,(T)-1),VECTOR<T,2>((T)1,(T)1))*grid_size);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
    fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
    if (test_number==5||test_number==7){
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
        fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;}
    if (test_number==6){
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;
        fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;}
    if (test_number==9){
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
        fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;}
    if(test_number==10){
        fluids_parameters.domain_walls[1][0]=true;}
    if(bucket_walls){
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
        fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;}
    else if(all_walls){
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;
        fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;}
    else if(no_walls){
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
        fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;}
    //time
    initial_time=(T)0.;last_frame=1000;frame_rate=(T)32.;
    if(strong_shock) frame_rate=2e4;
    if(test_number==10){
        last_frame=1000;frame_rate=1.25e5;}
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    e_min_for_clamping=1e-6;
    if(transition_to_incompressible){
        eos_smooth_transition=new EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >(time_start_transition,time_end_transition,one_over_c_incompressible,false);
        fluids_parameters.compressible_eos=eos_smooth_transition;}
    else fluids_parameters.compressible_eos=new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,e_min_for_clamping);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.use_slip=use_slip;
    fluids_parameters.use_preconditioner_for_slip_system=true;

    if(use_soot){
        fluids_parameters.use_soot=true;
        fluids_parameters.use_fixed_soot_boundary=true;
        fluids_parameters.ambient_soot=(T)0;
        fluids_parameters.soot_boundary=new BOUNDARY_REFLECTION_ATTENUATION<TV,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_soot,(T)1);}

    if(incompressible){
        fluids_parameters.use_vorticity_confinement=false;
        solids_fluids_parameters.use_leakproof_solve=false;
        fluids_parameters.use_body_force=false;
        fluids_parameters.density=(T)1; // not used
        fluids_parameters.use_density=true;
        fluids_parameters.use_temperature=true;
        fluids_parameters.use_fixed_density_boundary=true;
        fluids_parameters.use_fixed_temperature_boundary=true;
        if(use_incompressible_gravity) fluids_parameters.gravity=(T)9.8;
        else fluids_parameters.gravity=(T)0;}

    if(test_number==2||test_number==3) fluids_parameters.solid_affects_fluid=true;
    bool simulate_rigids=(test_number==4||test_number==5||test_number==6||test_number==7||test_number==9);
    if(use_slip) simulate_rigids=true;
    bool simulate_deformable=(test_number==8);
    if(simulate_rigids||simulate_deformable){
        solids_fluids_parameters.use_leakproof_solve=false;
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.use_post_cg_constraints=false;
        if(simulate_deformable){
            solid_body_collection.deformable_body_collection.simulate=true;}
        else{
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            solids_parameters.rigid_body_collision_parameters.use_push_out=true;
            solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;}
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.verbose_dt=true;
        solid_body_collection.print_residuals=false;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
        solids_parameters.implicit_solve_parameters.cg_projection_iterations=0; // TODO: check this
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=40;
        solids_parameters.implicit_solve_parameters.cg_iterations=1000;
    
        solids_parameters.use_rigid_deformable_contact=true;
        solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
        solid_gravity=(T)9.8;}

    motion_curve.Add_Control_Point(0,TV((T).65,(T).11)*grid_size);
    motion_curve.Add_Control_Point((T).04,TV((T).56,(T).02)*grid_size); // .03 was old
    motion_curve.Add_Control_Point((T).6,TV((T)0,(T).02)*grid_size);

    state_inside=TV_DIMENSION((T)1,(T)0,(T)0,(T)1);
    state_outside=TV_DIMENSION((T).125,(T)0,(T)0,(T).1);
    switch(test_number){
        case 7: shock_radius=(T)0; 
                shock_center=TV();
                break;
        case 9: shock_radius=(T).1;
                shock_center=TV(0,(T)-.65);
                strong_shock=true;
                break;
        case 10: shock_radius=(T)6.5;
                 shock_center=TV();
                 break;
        default: shock_radius=(T).4;
                 shock_center=TV();
                 break;}

    if(strong_shock){
        T temperature_inside=(T)2900,temperature_outside=(T)290;
        T p_atm=(T)101325.;
        T p_inside=(T)1000*p_atm,p_outside=p_atm;

        T rho_inside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_inside,temperature_inside);
        T rho_outside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_outside,temperature_outside);

        state_inside=TV_DIMENSION(rho_inside,(T)0,(T)0,p_inside);
        state_outside=TV_DIMENSION(rho_outside,(T)0,(T)0,p_outside);}

    if(test_number==10){
        T p_atm=(T)101325.;
        T temperature_inside=(T)2.62497e8,temperature_outside=(T)290;
        T p_inside=(T)9.41831e10,p_outside=p_atm;

        T rho_inside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_inside,temperature_inside);
        T rho_outside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_outside,temperature_outside);

        state_inside=TV_DIMENSION(rho_inside,(T)0,(T)0,p_inside);
        state_outside=TV_DIMENSION(rho_outside,(T)0,(T)0,p_outside);}

     if(incompressible){
         // Set ambient density and temperature
         EOS_GAMMA<T> eos_temp;
         T rho_outside=state_outside(0);
         T p_outside=state_outside(3);
         T e_outside=eos_temp.e_From_p_And_rho(p_outside,rho_outside);
         T temperature_outside=eos_temp.T(rho_outside,e_outside);
         fluids_parameters.ambient_density=rho_outside;
         fluids_parameters.ambient_temperature=temperature_outside;
         fluids_parameters.density_boundary=new BOUNDARY_REFLECTION_ATTENUATION<TV,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_density,(T).1);
         fluids_parameters.temperature_boundary=new BOUNDARY_REFLECTION_ATTENUATION<TV,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_temperature,(T).1);}

    // Set output directory
    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Circle_Example/Test_%d__Resolution_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
    else output_directory=STRING_UTILITIES::string_sprintf("Circle_Example/Test_%d__Resolution_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
    if(use_slip) output_directory+="_slip";
    if(transition_to_incompressible)  output_directory+="_transition_incompressible";
    if(use_soot) output_directory+="_soot";
    if(use_fixed_farfield_boundary) output_directory+="_fixedFF";
    if(strong_shock) output_directory+="_strong";
    output_directory+=STRING_UTILITIES::string_sprintf("_mass_%f",solid_mass);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    if(incompressible){
        fluids_parameters.Use_Fluid_Coupling_Defaults();
        return;}
    //set custom boundary

    TV far_field_velocity=TV(state_outside(1),state_outside(2));
    if(use_fixed_farfield_boundary){
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(0),state_outside(0),state_outside(0),state_outside(0)),
            T_FACE_VECTOR(state_outside(3),state_outside(3),state_outside(3),state_outside(3)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),true,T_FACE_VECTOR(1,1,1,1),T_FACE_VECTOR_BOOL(true,true,true,true));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(0),state_outside(0),state_outside(0),state_outside(0)),
            T_FACE_VECTOR(state_outside(3),state_outside(3),state_outside(3),state_outside(3)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    if(incompressible) return;

    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& U=fluids_parameters.euler->U;

    fluids_parameters.euler->e_min=e_min_for_clamping;
    if(transition_to_incompressible) fluids_parameters.euler->euler_projection.use_neumann_condition_for_outflow_boundaries=false;

    //initialize grid variables
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho,u_vel,v_vel,p;
        if((grid.X(cell_index)-shock_center).Magnitude()<shock_radius){rho=state_inside(0);u_vel=state_inside(1);v_vel=state_inside(2);p=state_inside(3);}
        else{rho=state_outside(0);u_vel=state_outside(1);v_vel=state_outside(2);p=state_outside(3);}

        U(cell_index)(0)=rho;U(cell_index)(1)=rho*u_vel;U(cell_index)(2)=rho*v_vel;U(cell_index)(3)=rho*(fluids_parameters.euler->eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))*((T).5));}

    // initialize solid_state
    VECTOR<T,T_GRID::dimension+2>& solid_state=fluids_parameters.euler_solid_fluid_coupling_utilities->solid_state;
    T rho=state_outside(0),p=state_outside(3),u_vel=(T)0.,v_vel=(T)0.;
    solid_state(0)=rho;solid_state(1)=rho*u_vel;solid_state(2)=rho*v_vel;solid_state(3)=rho*(fluids_parameters.euler->eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))*((T).5));
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(!incompressible) PHYSBAM_FATAL_ERROR("this shouldn't be called in compressible case");
    if(time>0) return;
    LOG::cout<<"INITIALIZING DENSITY"<<std::endl;
    
    T_GRID& grid=*fluids_parameters.grid;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho;
        if((grid.X(cell_index)-shock_center).Magnitude()<shock_radius) rho=state_inside(0); 
        else rho=state_outside(0);
        fluids_parameters.density_container.density(cell_index)=rho;}
}
//#####################################################################
// Function Adjust_Soot_With_Sources
//#####################################################################
void Adjust_Soot_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(!use_soot) PHYSBAM_FATAL_ERROR("this shouldn't be called in non use_soot case");
    if(time>0) return;
    LOG::cout<<"INITIALIZING SOOT"<<std::endl;
    
    T_GRID& grid=*fluids_parameters.grid;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T soot;
        if((grid.X(cell_index)-shock_center).Magnitude()<shock_radius) soot=1; 
        else soot=0;
        fluids_parameters.soot_container.density(cell_index)=soot;}
}
//#####################################################################
// Function Set_Rigid_Body_Parameters
//#####################################################################
void Set_Rigid_Body_Parameters(int rigid_body_index,TV position,T mass,bool is_kinematic)
{
    rigid_body_collection.Rigid_Body(rigid_body_index).Set_Coefficient_Of_Restitution((T)1);
    rigid_body_collection.Rigid_Body(rigid_body_index).coefficient_of_friction=(T)1;
    rigid_body_collection.rigid_body_particle.frame(rigid_body_index).t=position;
    rigid_body_collection.Rigid_Body(rigid_body_index).Set_Mass(mass);
    rigid_body_collection.Rigid_Body(rigid_body_index).Is_Kinematic()=is_kinematic;
    rigid_body_collection.Rigid_Body(rigid_body_index).simplicial_object->mesh.Initialize_Adjacent_Elements();
}
//#####################################################################
// Function Add_Sphere
//#####################################################################
void Add_Sphere(TV position,bool is_kinematic)
{
    sphere=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies_2D/circle",(T).25,true,true,false);
    Set_Rigid_Body_Parameters(sphere,position,solid_mass,is_kinematic);
}
//#####################################################################
// Function Add_Ground
//#####################################################################
void Add_Ground()
{
    T scale=(T)10;
    RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV::All_Ones_Vector()*scale);
    Set_Rigid_Body_Parameters(ground.particle_index,TV(0,-(T)1-scale*(T).5),(T)1e10,false);
}
//#####################################################################
// Function Add_Wall
//#####################################################################
void Add_Wall()
{
    T x_scale=(T).1,y_scale=(T).6;
    RIGID_BODY<TV>& wall=tests.Add_Analytic_Box(TV(x_scale,y_scale));
    Set_Rigid_Body_Parameters(wall.particle_index,TV(.4+x_scale*(T).5,y_scale*(T).5-1),(T)1e10,false);
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{   
    if(test_number==1) return;

    if(test_number==8){
        TRIANGULATED_AREA<T>& triangulated_area=tests.Create_Triangulated_Object(GRID<TV>(TV_INT()+2,RANGE<TV>(TV()-.1,TV()+.1)),RIGID_BODY_STATE<TV>(FRAME<TV>(motion_curve.Value(0))),(T)1);
        tests.Set_Mass_Of_Particles(triangulated_area,solid_mass/triangulated_area.Total_Size(),false);

        // correct number nodes
        for(int i=0;i<solid_body_collection.deformable_body_collection.deformable_geometry.structures.m;i++) solid_body_collection.deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

        // correct mass
        solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        
        solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_area,(T)1,(T)1.5));
        solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_area,(T)1));

        SEGMENTED_CURVE_2D<T>& segmented_curve=triangulated_area.Get_Boundary_Object();

        solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(segmented_curve);
        deformable_collisions.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(deformable_collisions);}
    else{
        bool is_kinematic=false;
        TV position=motion_curve.Value(0);;
        switch(test_number){
            case 2:is_kinematic=true;
                   Add_Sphere(position,is_kinematic);
                   break;
            case 3:is_kinematic=true;
                   Add_Sphere(position,is_kinematic);
                   break;
            case 4:Add_Sphere(position,is_kinematic);
                   break;
            case 5:Add_Sphere(position,is_kinematic);
                   break;
            case 6:Add_Sphere(position,is_kinematic);
                   break;
            case 7:position=TV(0,0);
                   Add_Sphere(position,is_kinematic);
                   break;
            case 9:Add_Ground();
                   Add_Wall();
                   break;
            case 10:Add_Ground();
                    break;
            default:PHYSBAM_FATAL_ERROR("wrong test number");}

        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);}
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);

    if(use_solids_gravity) solid_body_collection.Add_Force(new GRAVITY<TV>(solid_body_collection.deformable_body_collection.particles,rigid_body_collection,true,true,solid_gravity));

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==sphere) frame.t=motion_curve.Value(time);
    else if (test_number==2 && id==sphere) frame.t=motion_curve.Value(0);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==sphere){twist.linear=motion_curve.Derivative(time);return true;}
    return false;
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==3){
        T_GRID& grid=fluids_parameters.euler->grid;
        TV velocity=rigid_body_collection.rigid_body_particle.twist(sphere).linear;
        T rigid_dt_denominator=abs(velocity.x)/grid.dX.x+abs(velocity.y)/grid.dX.y;
        if(rigid_dt_denominator>(T)1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==7 && time>=(T)1 && time<=(T)1.01){
        wrench(0).linear=(T)10*TV(-1,0);}
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
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_rhs=print_rhs;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
    else if(fluids_parameters.fluid_affects_solid){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(*solids_evolution).print_matrix_rhs_and_solution=print_matrix;}
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(transition_to_incompressible){
        eos_smooth_transition->Set_Current_Time(time);
        //if(time>eos_smooth_transition->t_start_transition) fluids_parameters.euler->euler_projection.Set_Transition_To_Using_Implicit_Pressure(true);
    }
}
/*void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    LAPLACE_UNIFORM<T_GRID>* elliptic_solver=fluids_parameters.euler->euler_projection.elliptic_solver;
    BASE::Set_Dirichlet_Boundary_Conditions(time);
    elliptic_solver->psi_N.Component(0)(TV_INT(1,1))=false;
    elliptic_solver->psi_D(TV_INT(0,1))=true;
    //elliptic_solver->u(TV_INT(0,1))=state_outside(3);

    //for(CELL_ITERATOR iterator(fluids_parameters.euler->grid,1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
     //   TV_INT cell_index=iterator.Cell_Index();
      //  elliptic_solver->u(cell_index)=state_outside(3);}

}*/
//#####################################################################
};
}
#endif
