//#####################################################################
// Copyright 2002-2007 Doug Enright, Ron Fedkiw, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOD_ST
//#####################################################################
#ifndef __SOD_ST__
#define __SOD_ST__

#include <fstream>
#include <iostream>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SOD_ST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,1> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::Add_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    T solid_mass;
    bool write_transparency_output,transition_to_incompressible;
    TV_DIMENSION state_left,state_middle,state_right; // (density,velocity,pressure)
    T middle_state_start_point,right_state_start_point;
    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >* eos_smooth_transition;

    /***************
    example explanation:
    1. Sod shock tube with shock moving to the right.
    2. Sod shock tube with shock moving to the left.
    3. Fluid moving to the left and right at same velocity (similar to test 4 for pistons)
    4. Lax's shock tube problem
    5. Strong shock tube problem
    6. Two symmetric rarefaction waves (same as 3, but different values to match preconditioner paper)
    7. Mach 3 shock test
    8. High mach flow test
    9. Two shocks
    10.Interaction of blast waves (same as bang-bang)
    11. 1 with a 1-D coupled rigid body on the right and no right wall.
    12. 11 with a larger domain.
    13. 11 with walls on both sides
    14. Flat fluid with zero initial velocity, and a 1-D coupled rigid body moving to the right. There should be expansion/comression on the left/right sides of the solid.
    15. 1 with a 1-D deformable beam.
    16. 1d spring-mass system with right particle fixed.
    ***************/

    SOD_ST(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),tests(*this,solid_body_collection),
        rigid_body_collection(solid_body_collection.rigid_body_collection),solid_mass(0),write_transparency_output(false),
        transition_to_incompressible(false),eos_smooth_transition(0)
    {
    }

    virtual ~SOD_ST() {}

    // Unused Callbacks
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Velocities(ARRAY<T,FACE_INDEX<1> >&,ARRAY<bool,FACE_INDEX<1> >&,const T time) PHYSBAM_OVERRIDE {}
    bool Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,ARRAY<T>& rigid_simplex_forces,
        TV& orientation,const T time) PHYSBAM_OVERRIDE {return false;}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-eno_scheme",1,"eno_scheme","eno scheme");
    parse_args->Add_Integer_Argument("-eno_order",2,"eno_order","eno order");
    parse_args->Add_Integer_Argument("-rk_order",3,"rk_order","runge kutta order");
    parse_args->Add_Option_Argument("-noglf","don't use GLF and use LLF for ENO");
    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL","cfl number");
    parse_args->Add_Option_Argument("-cfl_sound_speed","use sound speed based cfl condition");
    parse_args->Add_Double_Argument("-cfl_sound_speed_multiple",(T)0.,"cfl_sound_speed_multiple","multiple of sound speed based cfl. Used if non-zero value set.");
    parse_args->Add_Double_Argument("-solid_mass",(T)1,"solid_mass","the mass of the solid in the simulation");
    parse_args->Add_Option_Argument("-timesplit","split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add_Option_Argument("-implicit_rk","perform runge kutta on the implicit part");
    parse_args->Add_Option_Argument("-all_walls","Add walls on all sides");
    parse_args->Add_Option_Argument("-no_walls","No walls on all sides");
    parse_args->Add_Option_Argument("-exact","output a fully-explicit sim to (output_dir)_exact");
    parse_args->Add_Option_Argument("-write_transparency_output","Akin to Ariente's tests, this allows us to visualize a rigid body as though it were 0-D.");
    parse_args->Add_Option_Argument("-transition_to_incompressible","transition to incompressible in a time window");
    parse_args->Add_Double_Argument("-time_start_transition",(T).5,"time to start transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-time_end_transition",(T).7,"time to end transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-one_over_c_incompressible",(T)0,"incompressible sound speed");
    parse_args->Add_Option_Argument("-apply_cavitation_correction");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    int eno_scheme=parse_args->Get_Integer_Value("-eno_scheme");
    int eno_order=parse_args->Get_Integer_Value("-eno_order");
    int rk_order=parse_args->Get_Integer_Value("-rk_order");
    bool use_glf=!parse_args->Is_Value_Set("-noglf");
    T cfl_number=(T)parse_args->Get_Double_Value("-cfl");
    T use_sound_speed_based_cfl=parse_args->Is_Value_Set("-cfl_sound_speed");
    T multiplication_factor_for_sound_speed_based_dt=(T)parse_args->Get_Double_Value("-cfl_sound_speed_multiple");
    solid_mass=(T)parse_args->Get_Double_Value("-solid_mass");
    bool timesplit=parse_args->Is_Value_Set("-timesplit");
    bool implicit_rk=parse_args->Is_Value_Set("-implicit_rk") && !parse_args->Is_Value_Set("-exact");
    bool all_walls=parse_args->Is_Value_Set("-all_walls");
    bool no_walls=parse_args->Is_Value_Set("-no_walls");
    write_transparency_output=parse_args->Is_Value_Set("-write_transparency_output");
    transition_to_incompressible=parse_args->Is_Value_Set("-transition_to_incompressible");
    T time_start_transition=(T)parse_args->Get_Double_Value("-time_start_transition");
    T time_end_transition=(T)parse_args->Get_Double_Value("-time_end_transition");
    T one_over_c_incompressible=(T)parse_args->Get_Double_Value("-one_over_c_incompressible");

    LOG::cout<<"use_glf="<<use_glf<<std::endl;

    //grid
    if(test_number==12) fluids_parameters.grid->Initialize(resolution,(T)-1,(T)2);
    else fluids_parameters.grid->Initialize(resolution,(T)0,(T)1);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
    if(test_number==1 || test_number==2) fluids_parameters.domain_walls[0][1]=false;
    if(test_number==10){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
    if(test_number==11 || test_number==12) fluids_parameters.domain_walls[0][1]=false;
    if(test_number==13 || test_number==14){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
    if(all_walls){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
    else if(no_walls){fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;}
    //time
    initial_time=(T)0.;last_frame=1500;frame_rate=(T)100.;
    if(test_number==9) last_frame=4000;
    if(test_number==5){frame_rate=(T)5/(T)2.5e-6;last_frame=500;}
    else if(test_number==8){frame_rate=(T)10/(T)1.75e-4;last_frame=1000;}
    else if(test_number==10){frame_rate=(T)10/(T).038;last_frame=500;}
    fluids_parameters.cfl=cfl_number;
    fluids_parameters.compressible_use_sound_speed_for_cfl=use_sound_speed_based_cfl;
    if(multiplication_factor_for_sound_speed_based_dt>0){
        fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl=true;
        fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt=multiplication_factor_for_sound_speed_based_dt;}
    //custom stuff . . .
    if(transition_to_incompressible){
        eos_smooth_transition=new EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >(time_start_transition,time_end_transition,one_over_c_incompressible);
        fluids_parameters.compressible_eos=eos_smooth_transition;}
    else fluids_parameters.compressible_eos=new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(use_glf,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(use_glf,true,false);
    else fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(use_glf,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    //fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_RF<T,T_GRID::dimension+2>;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;
    solid_body_collection.deformable_body_collection.simulate=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;

    if(parse_args->Is_Value_Set("-apply_cavitation_correction")) fluids_parameters.compressible_apply_cavitation_correction=true;

    if(test_number==11 || test_number==12 || test_number==13 || test_number==14 || test_number==15 || test_number==16){
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        if(test_number==15 || test_number==16){
            solid_body_collection.deformable_body_collection.simulate=true;}
        else{
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            solids_parameters.rigid_body_collision_parameters.use_push_out=false;
            solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;}
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.verbose_dt=true;
        solid_body_collection.print_residuals=false;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
        solids_parameters.implicit_solve_parameters.cg_projection_iterations=0; // TODO: check this
        solids_parameters.implicit_solve_parameters.cg_iterations=400;}

    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Sod_ST/Test_%d__Resolution_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x));
    else output_directory=STRING_UTILITIES::string_sprintf("Sod_ST/Test_%d__Resolution_%d_explicit",test_number,(fluids_parameters.grid->counts.x));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";

    if(test_number==11 || test_number==12 || test_number==13 || test_number==14 || test_number==15 || test_number==16){
        output_directory+=STRING_UTILITIES::string_sprintf("_mass_%f",solid_mass);}
    if(transition_to_incompressible) output_directory+="_transition_incompressible";
    if(fluids_parameters.compressible_apply_cavitation_correction) output_directory+="_cavitation";

    middle_state_start_point=0.5;right_state_start_point=0;
    if(test_number==1||test_number==11 || test_number==12 || test_number==13 || test_number==15 || test_number==16){
        state_left=TV_DIMENSION((T)1,(T)0,(T)1);
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);}
    else if(test_number==2){
        state_left=TV_DIMENSION((T).125,(T)0,(T).1);
        state_right=TV_DIMENSION((T)1,(T)0,(T)1);}
    else if(test_number==3){
        state_left=TV_DIMENSION((T)1,(T)-3,(T)1);
        state_right=TV_DIMENSION((T)1,(T)3,(T)1);}
    else if(test_number==4){
        state_left=TV_DIMENSION((T).445,(T).698,(T)3.528);
        state_right=TV_DIMENSION((T).5,(T)0,(T).571);}
    else if(test_number==5){
        state_left=TV_DIMENSION((T)1,(T)0,(T)1e10);
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);}
    else if(test_number==6){
        state_left=TV_DIMENSION((T)1,(T)-2,(T).4);
        state_right=TV_DIMENSION((T)1,(T)2,(T).4);}
    else if(test_number==7){
        state_left=TV_DIMENSION((T)3.857,(T).92,(T)10.333);
        state_right=TV_DIMENSION((T)1,(T)3.55,(T)1);}
    else if(test_number==8){
        state_left=TV_DIMENSION((T)10,(T)2000,(T)500);
        state_right=TV_DIMENSION((T)20,(T)0,(T)500);}
    else if(test_number==9){
        state_left=TV_DIMENSION((T).125,(T)0,(T).1);
        state_middle=TV_DIMENSION((T)1,(T)0,(T)1);
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);
        middle_state_start_point=(T).4;right_state_start_point=(T).6;}
    else if(test_number==10){
        state_left=TV_DIMENSION((T)1,(T)0,(T)1e3);
        state_middle=TV_DIMENSION((T)1,(T)0,(T)1e-2);
        state_right=TV_DIMENSION((T)1,(T)0,(T)1e2);
        middle_state_start_point=(T).1;right_state_start_point=(T).9;}
    else if(test_number==14){
        state_left=TV_DIMENSION((T).125,(T)0,(T).1);
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    if(transition_to_incompressible){
        TV_DIMENSION state_average=(state_left+state_right)*(T).5;
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,T_FACE_VECTOR(state_average(0),state_average(0)),
            T_FACE_VECTOR(state_average(2),state_average(2)),TV_FACE_VECTOR(TV(state_average(1)),TV(state_average(1))),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,T_FACE_VECTOR(state_left(0),state_right(0)),
            T_FACE_VECTOR(state_left(2),state_right(2)),TV_FACE_VECTOR(TV(state_left(1)),TV(state_right(1))),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(test_number<11) return;

    if(test_number==15 || test_number==16){
        SEGMENTED_CURVE<TV>& segmented_curve=tests.Create_Segmented_Curve(GRID<TV>(2,0,.1),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(.7))),(T)1);

        tests.Set_Mass_Of_Particles(segmented_curve,solid_mass/segmented_curve.Total_Size(),false);
        
        // correct number nodes
        for(int i=0;i<solid_body_collection.deformable_body_collection.deformable_geometry.structures.m;i++) solid_body_collection.deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

        // correct mass
        solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

        solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)10,(T)1.5));

        POINT_SIMPLICES_1D<T>& point_simplices_1d=segmented_curve.Get_Boundary_Object();

        for(int i=0;i<point_simplices_1d.mesh.elements.m;i++){
            LOG::cout<<"element "<<i<<"="<<point_simplices_1d.mesh.elements(i)<<", direction="<<point_simplices_1d.mesh.directions(i)<<std::endl;}
        LOG::cout<<"mass="<<solid_body_collection.deformable_body_collection.particles.mass<<std::endl;

        //point_simplices_1d.Update_Point_Simplex_List();

        solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&point_simplices_1d);
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(point_simplices_1d);
        //deformable_collisions.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(deformable_collisions);}
    else{
        T_GRID& grid=*fluids_parameters.grid;
        RANGE<TV> domain=grid.Domain();TV grid_size=domain.Edge_Lengths();
        TV scaling_factor=TV(grid_size.x/(T)6);
        RIGID_BODY<TV>& rect=tests.Add_Analytic_Box(scaling_factor);
        rect.Set_Name("bullet");
        LOG::cout<<"Setting solid mass to "<<solid_mass<<std::endl;
        rect.Set_Mass(solid_mass);
        TV epsilon=TV(.000);
        if(test_number==12) rect.X()=TV((T).7+grid.DX()/(T)2+scaling_factor/(T)2+epsilon);
        else if(test_number==14) {
            epsilon=grid.DX()*(T).1;
            rect.X()=grid.X(TV_INT(0))+scaling_factor/(T)2+epsilon;
            rect.Twist().linear=TV(.2);}
        else rect.X()=grid.X(TV_INT((T).7*TV(grid.counts)))+grid.DX()/(T)2+scaling_factor/(T)2+epsilon;
        rect.Is_Kinematic()=false;

        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS<T>& eos=*fluids_parameters.euler->eos;

    if(transition_to_incompressible) fluids_parameters.euler->euler_projection.use_neumann_condition_for_outflow_boundaries=false;
    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    for(int i=0;i<grid.counts.x;i++){
        T rho=0.,u=0.,p=0.;
        if(grid.Axis_X(i,0) <= middle_state_start_point){rho=state_left(0);u=state_left(1);p=state_left(2);}
        else if(grid.Axis_X(i,0) <= right_state_start_point){rho=state_middle(0);u=state_middle(1);p=state_middle(2);}
        else{rho=state_right(0);u=state_right(1);p=state_right(2);}

        U(i)(0)=rho; U(i)(1)=rho*u; U(i)(2)=rho*(eos.e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}

    if(test_number==11 || test_number==12 || test_number==13 || test_number==14 || test_number==15 || test_number==16){
        VECTOR<T,T_GRID::dimension+2>& solid_state=fluids_parameters.euler_solid_fluid_coupling_utilities->solid_state;
        T rho=(T).125,p=(T).1,u_vel=(T)0.;
        solid_state(0)=rho;solid_state(1)=rho*u_vel;solid_state(2)=rho*(eos.e_From_p_And_rho(p,rho)+(sqr(u_vel))/(T)2.);}
}
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(test_number==16) V(1)=TV();
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(test_number==16) V(1)=TV();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(write_transparency_output){ 
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(output_directory+STRING_UTILITIES::string_sprintf("/U_%d.txt",frame),false,false);
        for(CELL_ITERATOR it(*fluids_parameters.grid);it.Valid();it.Next()){
            if(fluids_parameters.euler->euler_projection.elliptic_solver->psi_D(it.Cell_Index())) continue;
            TV_DIMENSION U_cell=fluids_parameters.euler->U(it.Cell_Index());
            (*output)<<U_cell(0)<<"\t"<<U_cell(1)<<"\t"<<U_cell(2)<<"\t"
                <<EULER<T_GRID>::Get_Velocity(U_cell)<<"\t"<<EULER<T_GRID>::e(U_cell)<<std::endl;}
        delete output;}
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
//#####################################################################
};
}
#endif
