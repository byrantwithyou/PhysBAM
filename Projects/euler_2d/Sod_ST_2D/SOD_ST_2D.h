//#####################################################################
// Copyright 2007 Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOD_ST_2D
//#####################################################################
#ifndef __SOD_ST_2D__
#define __SOD_ST_2D__

#include <fstream>
#include <iostream>
#include "math.h"

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class SOD_ST_2D:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
public:
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    typedef VECTOR<T,2*TV::m> T_FACE_VECTOR;typedef VECTOR<TV,2*TV::m> TV_FACE_VECTOR;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::data_directory;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool exact;

    /***************
    example explanation:
    1. Shock moving in X direction
    2. Shock moving in Y direction
    3. Bullet. 1 with a 2-D coupled rigid body on the right and no right wall.
    ***************/

    SOD_ST_2D(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),tests(stream_type,data_directory,solid_body_collection),
        rigid_body_collection(solid_body_collection.rigid_body_collection),eno_scheme(1),eno_order(2),
        rk_order(3),cfl_number((T).6),timesplit(false),exact(false)
    {
    }
    
    virtual ~SOD_ST_2D() {}

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
    parse_args->Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;

    timesplit=timesplit && !exact;
    //grid
    int cells=resolution;
    T grid_size=(T)10.;
    if(test_number==1 || test_number==3) fluids_parameters.grid->Initialize(TV_INT(cells,cells/5),RANGE<TV>(-TV(grid_size/(T)2,grid_size/(T)10),TV(grid_size/(T)2,grid_size/(T)10)));
    else if(test_number==2) fluids_parameters.grid->Initialize(TV_INT(cells/5,cells),RANGE<TV>(-TV(grid_size/(T)10,grid_size/(T)2),TV(grid_size/(T)10,grid_size/(T)2)));
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
    if(test_number==3) fluids_parameters.domain_walls[0][1]=false;
    //time
    initial_time=(T)0.;last_frame=1000;frame_rate=(T)80.;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    fluids_parameters.compressible_eos=new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_timesplit=timesplit;

    if(test_number==3){
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_collision_parameters.use_push_out=false;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.verbose_dt=true;
        solid_body_collection.print_residuals=false;
        solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
        solids_parameters.implicit_solve_parameters.cg_projection_iterations=0; // TODO: check this
        solids_parameters.implicit_solve_parameters.cg_iterations=400;}

    if(timesplit) output_directory=LOG::sprintf("Sod_ST_2D/Test_%d__Resolution_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
    else output_directory=LOG::sprintf("Sod_ST_2D/Test_%d__Resolution_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    VECTOR<VECTOR<bool,2>,TV::m> valid_wall;
    for(int axis=0;axis<TV::m;axis++) for(int axis_side=0;axis_side<2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];
    fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,T_FACE_VECTOR((T).125,(T).125,(T).125,(T).125),
        T_FACE_VECTOR((T).1,(T).1,(T).1,(T).1),TV_FACE_VECTOR(),(T).5,valid_wall);
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(test_number!=3) return;
    RANGE<TV> domain=fluids_parameters.grid->Domain();TV grid_size=domain.Edge_Lengths();
    TV epsilon=grid_size*(T).1;
    RIGID_BODY<TV>& rect=tests.Add_Analytic_Box(TV(grid_size.x/(T)4,grid_size.y)+epsilon);
    rect.Frame().t=domain.Minimum_Corner()+TV((T).8*grid_size.x,(T).5*grid_size.y);
    rect.Is_Kinematic()=false;

    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);

    VECTOR<T,TV::m+2>& solid_state=fluids_parameters.euler_solid_fluid_coupling_utilities->solid_state;
    EOS_GAMMA<T> *tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);
    T rho=(T).125,p=(T).1,u_vel=(T)0.,v_vel=(T)0.;
    solid_state(0)=rho;solid_state(1)=rho*u_vel;solid_state(2)=rho*v_vel;solid_state(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))/(T)2.);
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    for(CELL_ITERATOR<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho=(T)0.,u_vel=(T)0.,v_vel=(T)0.,p=(T)0.;
        if(test_number==1 || test_number==3){
            if(grid.X(cell_index)(0)<0){rho=(T)1.;p=(T)1.;}
            else{rho=(T).125;p=(T).1;}}
        else if(test_number==2){
            if(grid.X(cell_index)(1)<0){rho=(T)1.;p=(T)1.;}
            else{rho=(T).125;p=(T).1;}}
        U(cell_index)(0)=rho;U(cell_index)(1)=rho*u_vel;U(cell_index)(2)=rho*v_vel;U(cell_index)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))/(T)2.);}
}
//#####################################################################
};
}
#endif
