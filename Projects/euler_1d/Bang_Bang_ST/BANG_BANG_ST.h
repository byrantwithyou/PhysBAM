//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BANG_BANG_ST
//#####################################################################
//
//#####################################################################
// Fedkiw - February 5, 2002
// Enright - September 9, 2003
//#####################################################################
#ifndef __BANG_BANG_ST__
#define __BANG_BANG_ST__

#include <fstream>
#include <iostream>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template <class T_input>
class BANG_BANG_ST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::parse_args;using BASE::resolution;

    BANG_BANG_ST(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE)
    {
    }
        
    virtual ~BANG_BANG_ST() {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-eno_scheme",1,"eno_scheme","eno scheme");
    parse_args->Add_Integer_Argument("-eno_order",2,"eno_order","eno order");
    parse_args->Add_Integer_Argument("-rk_order",3,"rk_order","runge kutta order");
    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL","cfl number");
    parse_args->Add_Option_Argument("-timesplit","split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add_Option_Argument("-implicit_rk","perform runge kutta on the implicit part");
    parse_args->Add_Option_Argument("-exact","output a fully-explicit sim to (output_dir)_exact");
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
    bool timesplit=parse_args->Is_Value_Set("-timesplit") && !parse_args->Is_Value_Set("-exact");
    bool implicit_rk=parse_args->Is_Value_Set("-implicit_rk");

    int cells=20*resolution+1;
    //grid
    fluids_parameters.grid->Initialize(cells,(T)0.,(T)1.);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
    //time
    initial_time=(T)0.;last_frame=10;frame_rate=(T)263;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . .
    fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_rungekutta_order=rk_order;
    //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,T_GRID::dimension+2>;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;

    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Bang_Bang_ST/Test_1__Resolution_%d_semiimplicit",(fluids_parameters.grid->counts.x));
    else output_directory=STRING_UTILITIES::string_sprintf("Bang_Bang_ST/Test_1__Resolution_%d_explicit",(fluids_parameters.grid->counts.x));
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
    fluids_parameters.compressible_boundary= new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,T_FACE_VECTOR((T)1.,(T)1.),T_FACE_VECTOR((T).1,(T).1),
        TV_FACE_VECTOR(),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    for(int i=0;i<grid.counts.x;i++) {
        T rho=(T)0.,u=(T)0.,p=(T)0.;
        if(grid.Axis_X(i,1)<(T).1){rho=(T)1.;p=(T)1000.;} else if(grid.Axis_X(i,1)<(T).9){rho=(T)1.;p=(T).01;} else{rho=(T)1.;p=(T)100.;}
        U(i)(1)=rho;U(i)(2)=rho*u;U(i)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}
}
//#####################################################################
};
}
#endif
