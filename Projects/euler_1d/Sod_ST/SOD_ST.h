//#####################################################################
// Copyright 2002-2007 Doug Enright, Ron Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOD_ST
//#####################################################################
#ifndef __SOD_ST__
#define __SOD_ST__

#include <fstream>
#include <iostream>

#include <Core/Vectors/VECTOR_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Conservation_Law_Solvers/HYBRID_SL_ENO_CONSERVATION.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SOD_ST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,1> >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    typedef VECTOR<T,2*TV::m> T_FACE_VECTOR;typedef VECTOR<TV,2*TV::m> TV_FACE_VECTOR;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;
    using BASE::user_last_frame;
    
    TV_DIMENSION state_left,state_middle,state_right; // (density,velocity,pressure)
    T middle_state_start_point,right_state_start_point;

    ARRAY<bool,FACE_INDEX<TV::m> > flux_face;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool use_sound_speed_based_cfl;
    bool multiplication_factor_for_sound_speed_based_dt;
    bool exact;

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
    ***************/

    SOD_ST(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
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

        timesplit=timesplit && !exact;
        //grid
        if(test_number==12) fluids_parameters.grid->Initialize(TV_INT(resolution),RANGE<TV>(TV(-1),TV(2)));
        else fluids_parameters.grid->Initialize(TV_INT(resolution),RANGE<TV>(TV(0),TV(1)));
        *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
        if(test_number==1 || test_number==2) fluids_parameters.domain_walls[0][1]=true;
        if(test_number==10){fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;}
        //time
        if(!user_last_frame) last_frame=1500;
        if(!this->user_frame_rate) frame_rate=(T)100.;
        if(test_number==9) if(!user_last_frame) last_frame=4000;
        if(test_number==5){
            frame_rate=(T)5/(T)2.5e-6;
            if(!user_last_frame) last_frame=500;}
        else if(test_number==8){
            frame_rate=(T)10/(T)1.75e-4;
            if(!user_last_frame) last_frame=1000;}
        else if(test_number==10){
            frame_rate=(T)10/(T).038;
            if(!user_last_frame) last_frame=500;}
        fluids_parameters.cfl=cfl_number;
        fluids_parameters.compressible_use_sound_speed_for_cfl=use_sound_speed_based_cfl;
        if(multiplication_factor_for_sound_speed_based_dt>0){
            fluids_parameters.compressible_use_sound_speed_based_dt_multiple_for_cfl=true;
            fluids_parameters.compressible_multiplication_factor_for_sound_speed_based_dt=multiplication_factor_for_sound_speed_based_dt;}
        //custom stuff . . .
        fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
        if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,false,false);
        else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,false);
        else if(eno_scheme==3) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,true);
        else if(eno_scheme==4) fluids_parameters.compressible_conservation_method = new HYBRID_SL_ENO_CONSERVATION<TV,TV::m+2>(flux_face,new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,false,false));
        fluids_parameters.compressible_spatial_order=eno_order;
        fluids_parameters.compressible_conservation_method->Save_Fluxes();
        fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
        fluids_parameters.compressible_rungekutta_order=rk_order;
        //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,TV::m+2>;
        fluids_parameters.compressible_timesplit=timesplit;
        fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;

        if(!this->user_output_directory){
            if(timesplit) viewer_dir.output_directory=LOG::sprintf("Sod_ST/Test_%d_%d_semiimplicit_%d_%1.2f",test_number,(fluids_parameters.grid->counts.x),eno_scheme,cfl_number);
            else viewer_dir.output_directory=LOG::sprintf("Sod_ST/Test_%d__Resolution_%d_explicit",test_number,(fluids_parameters.grid->counts.x));
            if(eno_scheme==2) viewer_dir.output_directory+="_density_weighted";
            else if(eno_scheme==3) viewer_dir.output_directory+="_velocity_weighted";}

        middle_state_start_point=0.5;right_state_start_point=0;
        if(test_number==1){
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
    }

    virtual ~SOD_ST() {}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    //set custom boundary
    fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,T_FACE_VECTOR(state_left(0),state_right(0)),
        T_FACE_VECTOR(state_left(2),state_right(2)),TV_FACE_VECTOR(TV(state_left(1)),TV(state_right(1))),(T).5,Complement(fluids_parameters.domain_walls));
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() override
{
    GRID<TV>& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    for(int i=0;i<grid.counts.x;i++){
        T rho=0.,u=0.,p=0.;
        if(grid.X(TV_INT(i)).x <= middle_state_start_point){rho=state_left(0);u=state_left(1);p=state_left(2);}
        else if(grid.X(TV_INT(i)).x <= right_state_start_point){rho=state_middle(0);u=state_middle(1);p=state_middle(2);}
        else{rho=state_right(0);u=state_right(1);p=state_right(2);}

        U(i)(0) = rho; U(i)(1) = rho*u; U(i)(2) = rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}

    flux_face.Resize(grid.Domain_Indices(3));
    for(FACE_ITERATOR<TV> iter(grid,3);iter.Valid();iter.Next())
        flux_face(iter.Full_Index()) = true;
}
//#####################################################################
};
}
#endif
