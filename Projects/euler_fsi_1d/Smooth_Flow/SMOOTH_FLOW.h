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

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SMOOTH_FLOW:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,1> >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    typedef VECTOR<T,2*TV::m> T_FACE_VECTOR;typedef VECTOR<TV,2*TV::m> TV_FACE_VECTOR;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::resolution;
    

    /***************
    example explanation:
    ***************/
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool use_sound_speed_based_cfl;
    bool multiplication_factor_for_sound_speed_based_dt;
    bool exact;

    SMOOTH_FLOW(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args,0,fluids_parameters.COMPRESSIBLE),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).5),timesplit(false),
        implicit_rk(false),use_sound_speed_based_cfl(false),multiplication_factor_for_sound_speed_based_dt(false),exact(false)
    {
        parse_args.Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
        parse_args.Add("-eno_order",&eno_order,"eno_order","eno order");
        parse_args.Add("-rk_order",&rk_order,"rk_order","runge kutta order");
        parse_args.Add("-cfl",&cfl_number,"CFL","cfl number");
        parse_args.Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
        parse_args.Add("-implicit_rk",&implicit_rk,"perform runge kutta on the implicit part");
        parse_args.Add("-cfl_sound_speed",&use_sound_speed_based_cfl,"use sound speed based cfl condition");
        parse_args.Add("-cfl_sound_speed_multiple",&multiplication_factor_for_sound_speed_based_dt,"cfl_sound_speed_multiple","multiple of sound speed based cfl. Used if non-zero value set.");
        parse_args.Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
        parse_args.Parse();

        implicit_rk=implicit_rk && !exact;

        //grid
        fluids_parameters.grid->Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box());
        *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
        //time
        initial_time=(T)0.;last_frame=50;frame_rate=(T)1000000.;
        fluids_parameters.cfl=cfl_number;
        //custom stuff . . . 
        fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
        if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,false,false);
        else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,false);
        else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,true);
        fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
        fluids_parameters.compressible_conservation_method->Save_Fluxes();
        fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-3);
        fluids_parameters.compressible_rungekutta_order=rk_order;
        //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,TV::m+2>;
        fluids_parameters.compressible_timesplit=timesplit;
        fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;
        fluids_parameters.compressible_use_sound_speed_for_cfl=use_sound_speed_based_cfl;
        if(multiplication_factor_for_sound_speed_based_dt>0){
            fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl=true;
            fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt=multiplication_factor_for_sound_speed_based_dt;}
        fluids_parameters.compressible_monitor_conservation_error=true;

        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        if(timesplit) output_directory=LOG::sprintf("Smooth_Flow/Test_1__Resolution_%d_semiimplicit",(fluids_parameters.grid->counts.x));
        else output_directory=LOG::sprintf("Smooth_Flow/Test_1__Resolution_%d_explicit",(fluids_parameters.grid->counts.x));
        if(eno_scheme==2) output_directory+="_density_weighted";
        else if(eno_scheme==3) output_directory+="_velocity_weighted";
    }
    
    virtual ~SMOOTH_FLOW() {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    //set custom boundary
    fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<TV,VECTOR<T,TV::m+2> >();
    fluids_parameters.compressible_pressure_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<TV,T>();
    fluids_parameters.euler->euler_projection.elliptic_solver->periodic_boundary(1)=true;
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() override
{
    GRID<TV>& grid=fluids_parameters.euler->grid;
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
        T x=grid.X(VECTOR<int,1>(i)).x;
        T rho,u,p;
        u=0;
        T p1=(T)60*cos((T)2*(T)pi*x/period)+(T)100*sin((T)4*(T)pi*x/period);
        p=p0+epsilon*p1;
        rho=pow(p/A,one_over_gamma);
        //rho=rho0*pow(p/p0,(T)1./gamma);
        U(i)(0)=rho;U(i)(1)=rho*u;U(i)(2)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);
        LOG::cout<<"i="<<i<<", U="<<U(i)<<std::endl;
    }
}
//#####################################################################
};
}
#endif
