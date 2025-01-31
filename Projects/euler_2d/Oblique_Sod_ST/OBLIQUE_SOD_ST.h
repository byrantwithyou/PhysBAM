//#####################################################################
// Copyright 2002-2007 Doug Enright, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OBLIQUE_SOD_ST
//#####################################################################
//
//#####################################################################
// Enright - September 10, 2003
//#####################################################################
#ifndef __OBLIQUE_SOD_ST__
#define __OBLIQUE_SOD_ST__

#include "math.h"
#include <fstream>
#include <iostream>

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
//#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D.h"

namespace PhysBAM{

template<class T_input>
class OBLIQUE_SOD_ST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
public: 
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;using BASE::fluids_parameters;using BASE::stream_type;
    using BASE::resolution;using BASE::user_last_frame;

    T angle;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool exact;

    OBLIQUE_SOD_ST(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args,0,fluids_parameters.COMPRESSIBLE),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).6),timesplit(false),exact(false)
    {
        parse_args.Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
        parse_args.Add("-eno_order",&eno_order,"eno_order","eno order");
        parse_args.Add("-rk_order",&rk_order,"rk_order","runge kutta order");
        parse_args.Add("-cfl",&cfl_number,"CFL","cfl number");
        parse_args.Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
        parse_args.Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
        parse_args.Parse();

        timesplit=timesplit && !exact;

        int scale=6;
        angle=atan((T)1.); //try straight up-down shock . . .
        //grid
        int cells_n=4*resolution,cells_m=scale*cells_n; //need divisibility of (n-1) by 1,2,4 . . .
        T H=(T)1.,L=scale*H;
        fluids_parameters.grid->Initialize(TV_INT(cells_m+1,cells_n+1),RANGE<TV>(TV(),TV(L,H)/sin(angle)));
        *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
        //time
        if(!user_last_frame) last_frame=400;
        if(!this->user_frame_rate) frame_rate=(T)80.;
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
        //conservation_method=new CONSERVATION_ENO_RF<T>;
        //set custom boundary
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D<TV>(angle,true,true);
        fluids_parameters.compressible_timesplit=timesplit;

        if(!this->user_output_directory){
            viewer_dir.output_directory="Oblique_Sod_ST/matlab";        
            if(timesplit) viewer_dir.output_directory=LOG::sprintf("Oblique_Sod_ST/Resolution_%d_%d_semiimplicit",(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
            else viewer_dir.output_directory=LOG::sprintf("Oblique_Sod_ST/Resolution_%d_%d_explicit",(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
            if(eno_scheme==2) viewer_dir.output_directory+="_density_weighted";
            else if(eno_scheme==3) viewer_dir.output_directory+="_velocity_weighted";}
    }
    
    virtual ~OBLIQUE_SOD_ST() {}

//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() override
{
    GRID<TV>& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //amount of shifting to do . .
    T recip_slope=(T)1./tan(angle);
    T initial_x_jump_pos=(T)2.25/sin(angle); //scaled grid location . . .

    //non-oblique for the moment . . .
    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
        T rho=(T)0.,u_vel=(T)0.,v_vel=(T)0.,p=(T)0.;
        T curr_x_jump_pos=initial_x_jump_pos + recip_slope*(j-1)*grid.dX.y;
        if(grid.X(TV_INT(i,0)).x < curr_x_jump_pos) {rho=(T)1.;p=(T)1.;} else {rho=(T).125;p=(T).1;}
        U(i,j)(0)=rho;U(i,j)(1)=rho*u_vel;U(i,j)(2)=rho*v_vel;U(i,j)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))/(T)2.);}
}
//#####################################################################
};
}
#endif
