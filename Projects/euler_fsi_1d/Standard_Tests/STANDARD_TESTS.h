//#####################################################################
// Copyright 2009 Nipun Kwatra, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

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
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/HYBRID_SL_ENO_CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_CLAMPED_INTERNAL_ENERGY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,1> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR; typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::solids_evolution;using BASE::Add_To_Fluid_Simulation;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    T solid_mass,solid_position_delta;
    TV_DIMENSION state_left,state_middle,state_right; // (density,velocity,pressure)
    T middle_state_start_point,right_state_start_point;
    T spring_stiffness,spring_overdamping_fraction;
    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >* eos_smooth_transition;
    bool transition_to_incompressible;
    bool simulate_rigids,simulate_deformable;

    T_FACE_ARRAYS_BOOL flux_face;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;

    std::ofstream gnuplot_file_stream;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool use_sound_speed_based_cfl;
    bool multiplication_factor_for_sound_speed_based_dt;
    bool flip_shock;
    bool use_slip;
    bool all_walls;
    bool no_walls;
    T time_start_transition;
    T time_end_transition;
    T one_over_c_incompressible;

    /***************
    example explanation:
    1. Sod shock tube with shock moving to the right.
    2. 1 with a 1-D coupled rigid body on the right, an enlarged domain and no right wall.
    3. 1d spring-mass system in a uniform fluid (Ariente Case A)
    4. 1d spring-mass system in a uniform fluid (Ariente Case B)
    5. 1 with a 1-D coupled rigid body on the right and walls on both sides.
    6. uniform liquid with a 1-D coupled rigid body moving really fast.
    7. Two symmetric rarefaction waves (Fluid moving to the left and right)
    8. Sedov Blast wave (Shu 4.2).
    9. Double rarefaction (Shu 4.3).
    10. 2 but mirrored.
   ***************/

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),tests(*this,solid_body_collection),
        rigid_body_collection(solid_body_collection.rigid_body_collection),solid_mass(1),solid_position_delta(0),eos_smooth_transition(0),
        run_self_tests(false),print_poisson_matrix(false),print_index_map(false),print_matrix(false),print_rhs(false),print_each_matrix(false),
        output_iterators(false),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).5),timesplit(false),
        implicit_rk(false),use_sound_speed_based_cfl(false),multiplication_factor_for_sound_speed_based_dt(false),flip_shock(false),
        use_slip(false),all_walls(false),no_walls(false),time_start_transition((T).5),time_end_transition((T).7),one_over_c_incompressible(0)
    {
    }

    virtual ~STANDARD_TESTS() 
    {
        gnuplot_file_stream.close();
    }

    // Unused Callbacks
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
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
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
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
    parse_args->Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
    parse_args->Add("-eno_order",&eno_order,"eno_order","eno order");
    parse_args->Add("-rk_order",&rk_order,"rk_order","runge kutta order");
    parse_args->Add("-cfl",&cfl_number,"CFL","cfl number");
    parse_args->Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add("-implicit_rk",&implicit_rk,"perform runge kutta on the implicit part");
    parse_args->Add("-cfl_sound_speed",&use_sound_speed_based_cfl,"use sound speed based cfl condition");
    parse_args->Add("-cfl_sound_speed_multiple",&multiplication_factor_for_sound_speed_based_dt,"cfl_sound_speed_multiple","multiple of sound speed based cfl. Used if non-zero value set.");
    parse_args->Add("-solid_mass",&solid_mass,"solid_mass","the mass of the solid in the simulation");
    parse_args->Add("-slip",&use_slip,"use slip/spd for coupling");
    parse_args->Add("-all_walls",&all_walls,"Add walls on all sides");
    parse_args->Add("-no_walls",&no_walls,"No walls on all sides");
    parse_args->Add("-solid_position_delta",&solid_position_delta,"value","Move solid from default positions by specified amount");

    parse_args->Add("-test_system",&run_self_tests,"run self tests");
    parse_args->Add("-print_poisson_matrix",&print_poisson_matrix,"print poisson matrix");
    parse_args->Add("-print_index_map",&print_index_map,"print index map");
    parse_args->Add("-print_matrix",&print_matrix,"print matrix");
    parse_args->Add("-print_rhs",&print_rhs,"print rhs");
    parse_args->Add("-print_each_matrix",&print_each_matrix,"print each matrix");
    parse_args->Add("-output_iterators",&output_iterators,"output iterators");
    parse_args->Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"do not use preconditioner");
    parse_args->Add("-preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"use preconditioner");
    
    parse_args->Add("-flip_shock",&flip_shock,"flip_shock");

    parse_args->Add("-transition_to_incompressible",&transition_to_incompressible,"transition to incompressible in a time window");
    parse_args->Add("-time_start_transition",&time_start_transition,"value","time to start transitioning to incompressible flow");
    parse_args->Add("-time_end_transition",&time_end_transition,"value","time to end transitioning to incompressible flow");
    parse_args->Add("-one_over_c_incompressible",&one_over_c_incompressible,"value","one over incompressible sound speed");

    parse_args->Add("-apply_cavitation_correction",&fluids_parameters.compressible_apply_cavitation_correction,"compressible_apply_cavitation_correction");
    parse_args->Add("-adaptive_time_step",&fluids_parameters.compressible_adaptive_time_step,"compressible_adaptive_time_step");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();

    if(!use_slip && timesplit){LOG::cout<<"*** NON-SLIP VERSION OF SEMI-IMPLICIT FLOW SOLVER NOT AVAILABLE ***"<<std::endl;exit(-1);}
    spring_stiffness=(T)5;

    //grid
    int cells=resolution;
    if(test_number==1 || test_number==2 || test_number==10) fluids_parameters.grid->Initialize(cells+1,(T)-1,(T)3);
    else if(test_number==3||test_number==4) fluids_parameters.grid->Initialize(cells+1,(T)0,(T)20);
    else if(test_number==5) fluids_parameters.grid->Initialize(cells+1,(T)0,(T)3);
    else if(test_number==8) fluids_parameters.grid->Initialize(cells+1,(T)-2,(T)2);
    else if(test_number==9) fluids_parameters.grid->Initialize(cells+1,(T)-1,(T)1);
    else fluids_parameters.grid->Initialize(cells+1,(T)0,(T)1);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid();
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=false;
    if(test_number==5){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
    if(test_number==7 || test_number==8 || test_number==9){fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;}
    if(all_walls){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
    else if(no_walls){fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;}
    //time
    initial_time=(T)0.;last_frame=1500;frame_rate=(T)100.;
    if(test_number==3){frame_rate=(T)2e5;last_frame=3000;}
    if(test_number==4){frame_rate=(T)1e6;last_frame=6000;}
    if(test_number==8){frame_rate=(T)1e5;last_frame=100;}
    if(test_number==9){frame_rate=(T)1e2;last_frame=60;}
    fluids_parameters.cfl=cfl_number;
    fluids_parameters.compressible_use_sound_speed_for_cfl=use_sound_speed_based_cfl;
    if(multiplication_factor_for_sound_speed_based_dt>0){
        fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl=true;
        fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt=multiplication_factor_for_sound_speed_based_dt;}
    //custom stuff . . .
    if(transition_to_incompressible){
        eos_smooth_transition=new EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >(time_start_transition,time_end_transition,one_over_c_incompressible);
        fluids_parameters.compressible_eos=eos_smooth_transition;}
    //else if(fluids_parameters.compressible_apply_cavitation_correction){
    //    T e_min = 1e-14;
    //    T epsilon = e_min*.001;
    //    fluids_parameters.compressible_eos=new EOS_CLAMPED_INTERNAL_ENERGY<T>(*new EOS_GAMMA<T>, e_min, epsilon);
    //}
    else fluids_parameters.compressible_eos=new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else if(eno_scheme==3) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    else fluids_parameters.compressible_conservation_method = new HYBRID_SL_ENO_CONSERVATION<T_GRID,T_GRID::dimension+2>(flux_face,new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false));
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    //fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,T_GRID::dimension+2>;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;
    fluids_parameters.use_slip=use_slip;
    solid_body_collection.deformable_body_collection.simulate=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.use_post_cg_constraints=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;

    fluids_parameters.use_preconditioner_for_slip_system=true;

    simulate_rigids=test_number==2||test_number==5||test_number==6||test_number==10;
    simulate_deformable=test_number==3||test_number==4;

    if(simulate_rigids || simulate_deformable){
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        if(simulate_deformable){
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
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=40;
        solids_parameters.implicit_solve_parameters.cg_iterations=400;}

    middle_state_start_point=(T)0.5;right_state_start_point=0;
    if(test_number==1||test_number==2||test_number==5){
        // middle_state_start_point=(T)-1;right_state_start_point=(T)-1;
        state_left=TV_DIMENSION((T)1,(T)0,(T)1);
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);
        if(test_number==5) middle_state_start_point=(T)1.0;}
    else if(test_number==10){
        // middle_state_start_point=(T)-1;right_state_start_point=(T)-1;
        state_left=TV_DIMENSION((T).125,(T)0,(T).1);
        state_right=TV_DIMENSION((T)1,(T)0,(T)1);}
    else if(test_number==3 || test_number==4){
        middle_state_start_point=(T)1.0;
        state_left=TV_DIMENSION((T)4.,(T)0,(T)1e6);
        state_right=state_left;}
    else if(test_number==6){
        middle_state_start_point=(T)-10;
        right_state_start_point=(T)-10;
        state_right=TV_DIMENSION((T).125,(T)0,(T).1);}
    else if(test_number==7){
        state_left=TV_DIMENSION((T)1,(T)-2,(T).4);
        state_right=TV_DIMENSION((T)1,(T)2,(T).4);}
    else if(test_number==8){
        T_GRID& grid=*fluids_parameters.grid;
        T center_cell_location=grid.X((int)((grid.numbers_of_cells.x+1)*.5)).x;
        T dx=grid.DX().x;
        middle_state_start_point=center_cell_location - dx*.5;
        right_state_start_point=center_cell_location + dx*.5;

        //T internal_energy_usual=1e-12;
        T internal_energy_usual=2e-3;
        //T internal_energy_center_cell=32e5/dx;
        T internal_energy_center_cell=32e3/dx;
        T rho=(T)1;
        T u=(T)0;

        EOS<T>& eos=*fluids_parameters.compressible_eos;
        T pressure_usual=eos.p(rho,internal_energy_usual);
        T pressure_center_cell=eos.p(rho,internal_energy_center_cell);
        state_left=TV_DIMENSION(rho,u,pressure_usual);
        state_right=TV_DIMENSION(rho,u,pressure_usual);
        state_middle=TV_DIMENSION(rho,u,pressure_center_cell);}
    else if(test_number==9){
        middle_state_start_point=(T)0;right_state_start_point=middle_state_start_point;
        state_left=TV_DIMENSION((T)7,(T)-1,(T)0.2);
        state_right=TV_DIMENSION((T)7,(T)1,(T)0.2);
    }

    if(flip_shock){
        exchange(fluids_parameters.domain_walls[0][0],fluids_parameters.domain_walls[0][1]);
        middle_state_start_point=(T)3.-middle_state_start_point;
        exchange(state_left,state_right);}

    if(test_number==3){ 
        solid_mass=(T)6;
        spring_stiffness=1e7;
        spring_overdamping_fraction=0;}
    else if(test_number==4){
        solid_mass=(T).04;
        spring_stiffness=2e7;
        spring_overdamping_fraction=0;}

    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d__Resolution_%d",test_number,(fluids_parameters.grid->counts.x));
    if(timesplit) output_directory+="_semiimplicit";
    else output_directory+="_explicit";
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
    if(use_slip) output_directory+="_slip";
    if(transition_to_incompressible)  output_directory+="_transition_incompressible";
    if(simulate_rigids || simulate_deformable) output_directory+=STRING_UTILITIES::string_sprintf("_mass_%f",solid_mass);
    if(use_sound_speed_based_cfl) output_directory+="_using_sound_speed_cfl";
    if(fluids_parameters.compressible_apply_cavitation_correction) output_directory+="_cavitation";
    if(fluids_parameters.compressible_adaptive_time_step) output_directory+="_adaptive";

    output_directory+=STRING_UTILITIES::string_sprintf("_rk%d_delta%.4f",rk_order,solid_position_delta);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {
    BASE::Parse_Late_Options();
    std::string gnuplot_file=output_directory+"/common/gnuplot_data.dat";
    gnuplot_file_stream.open(gnuplot_file.c_str());
    LOG::cout<<"writing to file "<<gnuplot_file<<std::endl;
}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    T inflow_attenuation=(T).5;
    if(test_number==3||test_number==4) inflow_attenuation=(T)1;
    if(transition_to_incompressible){
        TV_DIMENSION state_average=(state_left+state_right)*(T).5;
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_average(1),state_average(1)),T_FACE_VECTOR(state_average(3),state_average(3)),
            TV_FACE_VECTOR(TV(state_average(2)),TV(state_average(2))),inflow_attenuation,
            VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}

    if(test_number==6){
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_right(1),state_right(1)),T_FACE_VECTOR(state_right(3),state_right(3)),
            TV_FACE_VECTOR(TV(state_right(2)),TV(state_right(2))),inflow_attenuation,
            VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),true,T_FACE_VECTOR((T)0,(T)0),T_FACE_VECTOR_BOOL(true,true));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_left(1),state_right(1)),T_FACE_VECTOR(state_left(3),state_right(3)),
            TV_FACE_VECTOR(TV(state_left(2)),TV(state_right(2))),inflow_attenuation,
            VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),true,T_FACE_VECTOR((T)0,(T)0),T_FACE_VECTOR_BOOL(true,true));}
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(!(simulate_rigids || simulate_deformable)) return;

    if(simulate_deformable){
        SEGMENTED_CURVE<TV>& segmented_curve=tests.Create_Segmented_Curve(GRID<TV>(2,0,1),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(19.5))),(T)1);

        tests.Set_Mass_Of_Particles(segmented_curve,solid_mass/segmented_curve.Total_Size(),false);
        
        // correct number nodes
        for(int i=0;i<solid_body_collection.deformable_body_collection.deformable_geometry.structures.m;i++){
            solid_body_collection.deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();}

        // correct mass
        solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

        solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)spring_stiffness,(T)spring_overdamping_fraction));

        POINT_SIMPLICES_1D<T>& point_simplices_1d=segmented_curve.Get_Boundary_Object();

        LOG::cout<<"point_simplices_1d elements="<<point_simplices_1d.mesh.elements<<", directions=point_simplices_1d.mesh.directions"<<std::endl;
        LOG::cout<<"mass="<<solid_body_collection.deformable_body_collection.particles.mass<<std::endl;
        LOG::cout<<"spring_stiffness="<<spring_stiffness<<", spring_overdamping_fraction="<<spring_overdamping_fraction<<std::endl;

        solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&point_simplices_1d);
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(point_simplices_1d);
        Add_To_Fluid_Simulation(deformable_collisions);}
    else if(simulate_rigids){
        T_GRID& grid=*fluids_parameters.grid;
        RANGE<TV> domain=grid.Domain();TV grid_size=domain.Edge_Lengths();
        TV scaling_factor=grid.DX();
        LOG::cout<<"Setting solid size to "<<scaling_factor<<std::endl;
        RIGID_BODY<TV>& rect=tests.Add_Analytic_Box(scaling_factor);
        rect.name="bullet";
        LOG::cout<<"Setting solid mass to "<<solid_mass<<std::endl;
        rect.Set_Mass(solid_mass);
        if(test_number==5) rect.Frame().t=TV(1.5);
        // else if(test_number==6) rect.Frame().t=TV(.9);
        // else rect.Frame().t=TV(.8);
        rect.Frame().t=grid.X(grid.Cell(TV(1.3),0))+solid_position_delta*grid.DX();
        rect.Is_Kinematic()=false;

        Add_Volumetric_Body_To_Fluid_Simulation(rect,true,true);}

    if(fluids_parameters.use_slip){
        for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
        for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
        for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    fluids_parameters.euler->euler_projection.use_neumann_condition_for_outflow_boundaries=false;

    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS<T>& eos=*fluids_parameters.euler->eos;

    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    // for(int i=0;i<grid.counts.x;i++){
    for(CELL_ITERATOR iter(grid);iter.Valid();iter.Next()){
        TV location=iter.Location();
        T rho=0.,u=0.,p=0.;
        if(location.x <= middle_state_start_point){rho=state_left(1);u=state_left(2);p=state_left(3);}
        else if(location.x <= right_state_start_point){rho=state_middle(1);u=state_middle(2);p=state_middle(3);}
        else{rho=state_right(1);u=state_right(2);p=state_right(3);}

        U(iter.Cell_Index())(1) = rho; U(iter.Cell_Index())(2) = rho*u; U(iter.Cell_Index())(3) = rho*(eos.e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}

    flux_face.Resize(grid.Domain_Indices(3));
    for(FACE_ITERATOR iter(grid,3);iter.Valid();iter.Next())
        flux_face(iter.Full_Index()) = true;
}
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(test_number==3 || test_number==4) V(2)=TV();
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(test_number==3 || test_number==4) V(2)=TV();
}
void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE 
{
    if(test_number==3) dt=min(dt,(T)1e-5);
    else if(test_number==4) dt=min(dt,(T)1e-7);
}
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE 
{
    T position(0),velocity(0);
    if(test_number==2 || test_number==10 || test_number==5){
        RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
        int rigid_body_index=1;
        position=rigid_body_particles.frame(rigid_body_index).t.x;
        velocity=rigid_body_particles.twist(rigid_body_index).linear.x;}
    else if(test_number==3 || test_number==4){
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        int particle_index=1;
        position=particles.X(particle_index).x;
        velocity=particles.V(particle_index).x;}

    LOG::cout<<"writing to file :::::"<<""<<time<<" "<<position<<" "<<velocity<<std::endl;
    gnuplot_file_stream<<""<<time<<" "<<position<<" "<<velocity<<std::endl;
    gnuplot_file_stream.flush();
}
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(fluids_parameters.use_slip){
        SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution);
        evolution.run_self_tests=run_self_tests;
        evolution.output_iterators=output_iterators;
        evolution.print_matrix=print_matrix;
        evolution.print_rhs=print_rhs;
        evolution.print_each_matrix=print_each_matrix;
        evolution.print_poisson_matrix=print_poisson_matrix;
        evolution.print_index_map=print_index_map;}
    else if(fluids_parameters.fluid_affects_solid){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(*solids_evolution).print_matrix_rhs_and_solution=print_matrix;}
}
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==2 || test_number==10 || test_number==5){
        RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
        int rigid_body_index=1;
        T velocity=rigid_body_particles.twist(rigid_body_index).linear.x;
        LOG::cout<<"KINETIC ENERGY = "<<(T).5*solid_mass*velocity*velocity<<std::endl;
    }

    if(transition_to_incompressible){
        eos_smooth_transition->Set_Current_Time(time);
        //if(time>eos_smooth_transition->t_start_transition) fluids_parameters.euler->euler_projection.Set_Transition_To_Using_Implicit_Pressure(true);
    }
    // T_GRID& grid=fluids_parameters.euler->grid;
    // for(FACE_ITERATOR iter(grid,3);iter.Valid();iter.Next())
    //     flux_face(iter.Full_Index()) = true;
}
virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    BASE::Write_Output_Files(frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(fluids_parameters.euler->timesplit){
        ARRAY<bool,TV_INT> irregular_cells(fluids_parameters.euler_solid_fluid_coupling_utilities->near_interface);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",irregular_cells);}
}
//#####################################################################
};
}
#endif
