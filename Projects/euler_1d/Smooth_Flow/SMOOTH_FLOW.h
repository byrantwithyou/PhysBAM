//#####################################################################
// Copyright 2007, Jon Gretarsson, Roger Huang, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOOTH_FLOW
//#####################################################################
#ifndef __SMOOTH_FLOW__
#define __SMOOTH_FLOW__

#include <fstream>
#include <iostream>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SMOOTH_FLOW:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::parse_args;using BASE::resolution;

    /***************
    example explanation:
    ***************/

    SMOOTH_FLOW(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE)
    {
    }
    
    virtual ~SMOOTH_FLOW() {}

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
    parse_args->Add_Option_Argument("-cfl_sound_speed","use sound speed based cfl condition");
    parse_args->Add_Double_Argument("-cfl_sound_speed_multiple",(T)0.,"cfl_sound_speed_multiple","multiple of sound speed based cfl. Used if non-zero value set.");
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
    bool use_sound_speed_based_cfl=parse_args->Is_Value_Set("-cfl_sound_speed");
    bool multiplication_factor_for_sound_speed_based_dt=parse_args->Is_Value_Set("-cfl_sound_speed_multiple");

    //grid
    fluids_parameters.grid->Initialize(resolution,(T)0.,(T)1.);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;
    //time
    initial_time=(T)0.;last_frame=50;frame_rate=(T)1000000.;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-3);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,T_GRID::dimension+2>;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;
    fluids_parameters.compressible_use_sound_speed_for_cfl=use_sound_speed_based_cfl;
    if(multiplication_factor_for_sound_speed_based_dt>0){
        fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl=true;
        fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt=multiplication_factor_for_sound_speed_based_dt;}
    fluids_parameters.compressible_monitor_conservation_error=true;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Smooth_Flow/Test_1__Resolution_%d_semiimplicit",(fluids_parameters.grid->counts.x));
    else output_directory=STRING_UTILITIES::string_sprintf("Smooth_Flow/Test_1__Resolution_%d_explicit",(fluids_parameters.grid->counts.x));
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
    fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,VECTOR<T,T_GRID::dimension+2> >();
    fluids_parameters.compressible_pressure_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T>();
    fluids_parameters.euler->euler_projection.elliptic_solver->periodic_boundary(1)=true;
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T>*tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);
    
    T period=1.;
    T gamma=tmp_eos->gamma;
    T one_over_gamma=1/gamma;
    //T rho0=(T)1e-3,p0=(T)1e6,epsilon=(T)1; // actual
    T rho0=(T)1,p0=(T)1e9,epsilon=(T)1e3; // scaling mass (kg->g)
    //T rho0=(T)1e-3,p0=(T)1,epsilon=(T)1e-3;
    //T rho0=(T)1,p0=(T)1,epsilon=(T)1e-3;
    T A=p0/(pow(rho0,gamma));
    for(int i=0;i<grid.counts.x;i++){
        T x=grid.Axis_X(i,0);
        T rho,u,p;
        u=0;
        T p1=(T)60*cos((T)2*(T)pi*x/period)+(T)100*sin((T)4*(T)pi*x/period);
        p=p0+epsilon*p1;
        rho=pow(p/A,one_over_gamma);
        //rho=rho0*pow(p/p0,(T)1./gamma);
        U(i)(1)=rho;U(i)(2)=rho*u;U(i)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);
        LOG::cout<<"i="<<i<<", U="<<U(i)<<std::endl;
    }
}
//#####################################################################
};
}
#endif
