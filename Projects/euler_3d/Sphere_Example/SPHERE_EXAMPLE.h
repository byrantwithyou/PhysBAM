//#####################################################################
// Copyright 2009 Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_EXAMPLE
//#####################################################################
#ifndef __SPHERE_EXAMPLE__
#define __SPHERE_EXAMPLE__

#include <fstream>
#include <iostream>
#include "math.h"

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class SPHERE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::stream_type;
    using BASE::data_directory;using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::solids_evolution;
    
    TV_DIMENSION state_inside,state_outside; // // (density,velocity_x,velocity_y,velocity_z,pressure)
    T shock_radius;
    T solid_mass;
    TV solid_initial_position;
    int sphere;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >* eos_smooth_transition;
    bool transition_to_incompressible;
    bool simulate_rigids,simulate_deformable;

    bool use_fixed_farfield_boundary;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;


    /***************
    example explanation:
    1. Circular higher density region.
    2. Same as 1 + a two-way coupled circular object.
    ***************/

    SPHERE_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),rigid_body_collection(solid_body_collection.rigid_body_collection)
    {
    }
    
    virtual ~SPHERE_EXAMPLE() {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-eno_scheme",2,"eno_scheme","eno scheme");
    parse_args->Add_Integer_Argument("-eno_order",2,"eno_order","eno order");
    parse_args->Add_Integer_Argument("-rk_order",3,"rk_order","runge kutta order");
    parse_args->Add_Double_Argument("-cfl",(T).6,"CFL","cfl number");
    parse_args->Add_Double_Argument("-mass",(T).0625,"solid_mass","the mass of the solid in the simulation");
    parse_args->Add_Option_Argument("-timesplit","split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add_Option_Argument("-slip","use slip/spd for coupling");
    parse_args->Add_Option_Argument("-exact","output a fully-explicit sim to (output_dir)_exact");

    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-print_poisson_matrix");
    parse_args->Add_Option_Argument("-print_index_map");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-print_rhs");
    parse_args->Add_Option_Argument("-print_each_matrix");
    parse_args->Add_Option_Argument("-output_iterators");
    parse_args->Add_Option_Argument("-no_preconditioner");
    parse_args->Add_Option_Argument("-preconditioner");

    parse_args->Add_Option_Argument("-transition_to_incompressible","transition to incompressible in a time window");
    parse_args->Add_Double_Argument("-time_start_transition",(T).5,"time to start transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-time_end_transition",(T).7,"time to end transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-one_over_c_incompressible",(T)0,"one over incompressible sound speed");
    
    parse_args->Add_Option_Argument("-use_fixed_farfield_boundary","use fixed farfield values for outflow boundaries");
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
    T cfl_number=(T)parse_args->Get_Double_Value("-cfl");
    solid_mass=(T)parse_args->Get_Double_Value("-mass");
    bool timesplit=parse_args->Is_Value_Set("-timesplit") && !parse_args->Is_Value_Set("-exact");
    bool use_slip=parse_args->Is_Value_Set("-slip");
    run_self_tests=parse_args->Is_Value_Set("-test_system");
    print_poisson_matrix=parse_args->Is_Value_Set("-print_poisson_matrix");
    print_index_map=parse_args->Is_Value_Set("-print_index_map");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    print_rhs=parse_args->Is_Value_Set("-print_rhs");
    print_each_matrix=parse_args->Is_Value_Set("-print_each_matrix");
    output_iterators=parse_args->Is_Value_Set("-output_iterators");
    transition_to_incompressible=parse_args->Is_Value_Set("-transition_to_incompressible");
    T time_start_transition=(T)parse_args->Get_Double_Value("-time_start_transition");
    T time_end_transition=(T)parse_args->Get_Double_Value("-time_end_transition");
    T one_over_c_incompressible=(T)parse_args->Get_Double_Value("-one_over_c_incompressible");
    use_fixed_farfield_boundary=parse_args->Is_Value_Set("-use_fixed_farfield_boundary");

    //grid
    int cells=resolution;
    T grid_size=(T)1.;
    fluids_parameters.grid->Initialize(cells,cells,cells,RANGE<TV>(TV((T)-1,(T)-1,(T)-1),TV((T)1,(T)1,(T)1))*grid_size);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[3][1]=false;fluids_parameters.domain_walls[3][2]=false;
    //time
    initial_time=(T)0.;last_frame=1000;frame_rate=(T)100.;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    if(transition_to_incompressible){
        eos_smooth_transition=new EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >(time_start_transition,time_end_transition,one_over_c_incompressible);
        fluids_parameters.compressible_eos=eos_smooth_transition;}
    else fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,1,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.use_slip=use_slip;

    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-no_preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=false;

    if(test_number==2){
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        solids_parameters.use_post_cg_constraints=false;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.rigid_body_collision_parameters.use_push_out=false;
        solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.verbose_dt=true;
        solid_body_collection.print_residuals=false;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
        solids_parameters.implicit_solve_parameters.cg_projection_iterations=0; // TODO: check this
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=40;
        solids_parameters.implicit_solve_parameters.cg_iterations=1000;

        solid_initial_position=TV((T).65,(T).11,(T)0);}

    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Sphere_Example/Test_%d__Resolution_%d_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),
        (fluids_parameters.grid->counts.y),(fluids_parameters.grid->counts.z));
    else output_directory=STRING_UTILITIES::string_sprintf("Sphere_Example/Test_%d__Resolution_%d_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y),
        (fluids_parameters.grid->counts.z));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
    if(use_slip) output_directory+="_slip";
    if(transition_to_incompressible) output_directory+="_transition_incompressible";
    if(use_fixed_farfield_boundary) output_directory+="_fixedFF";
    output_directory+=STRING_UTILITIES::string_sprintf("_mass_%f",solid_mass);

    state_inside=TV_DIMENSION((T)1,(T)0,(T)0,(T)0,(T)1);
    state_outside=TV_DIMENSION((T).125,(T)0,(T)0,(T)0,(T).1);
    shock_radius=(T).4;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    VECTOR<VECTOR<bool,2>,T_GRID::dimension> valid_wall;
    for(int axis=1;axis<=T_GRID::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];

    TV far_field_velocity=TV(state_outside(2),state_outside(3),state_outside(4));
    if(use_fixed_farfield_boundary){
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1)),
            T_FACE_VECTOR(state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,valid_wall,true,T_FACE_VECTOR(1,1,1,1,1,1),T_FACE_VECTOR_BOOL(true,true,true,true,true,true));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1)),
            T_FACE_VECTOR(state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,valid_wall);}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State()
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,5> ,VECTOR<int,3> >& U=fluids_parameters.euler->U;
    EOS<T> *eos = fluids_parameters.euler->eos;

    if(transition_to_incompressible) fluids_parameters.euler->euler_projection.use_neumann_condition_for_outflow_boundaries=false;

    //initialize grid variables
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho,u_vel,v_vel,w_vel,p;
        if(grid.X(cell_index).Magnitude()<shock_radius){rho=state_inside(1);u_vel=state_inside(2);v_vel=state_inside(3);w_vel=state_inside(4);p=state_inside(5);}
        else{rho=state_outside(1);u_vel=state_outside(2);v_vel=state_outside(3);w_vel=state_outside(4);p=state_outside(5);}

        U(cell_index)(1)=rho;U(cell_index)(2)=rho*u_vel;U(cell_index)(3)=rho*v_vel;U(cell_index)(4)=rho*w_vel;
        U(cell_index)(5)=rho*(eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel)+sqr(w_vel))/(T)2.);}
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{   
    if(test_number==1) return;

    sphere=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/sphere",(T).25,true,true,false);
    rigid_body_collection.Rigid_Body(sphere).Set_Coefficient_Of_Restitution((T)1);rigid_body_collection.Rigid_Body(sphere).Set_Coefficient_Of_Friction((T)1);
    rigid_body_collection.rigid_body_particle.X(sphere)=solid_initial_position;rigid_body_collection.Rigid_Body(sphere).Set_Mass(solid_mass);
    rigid_body_collection.Rigid_Body(sphere).Is_Kinematic()=false;
    rigid_body_collection.Rigid_Body(sphere).simplicial_object->mesh.Initialize_Adjacent_Elements();

    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
}
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
