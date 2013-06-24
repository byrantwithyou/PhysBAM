//#####################################################################
// Copyright 2007-2008 Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WIND_TUNNEL
//#####################################################################
#ifndef __WIND_TUNNEL__
#define __WIND_TUNNEL__

#include <fstream>
#include <iostream>
#include "math.h"

#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class WIND_TUNNEL:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >,public BOUNDARY<VECTOR<T_input,2>,VECTOR<T_input,4> >
{
public:
    typedef T_input T;typedef VECTOR<T,2> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename COLLISION_BODY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;
    typedef VECTOR<BOUNDARY<TV,TV_DIMENSION>*,2*T_GRID::dimension> T_BOUNDARY_FACE_VECTOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    typedef BOUNDARY<TV,TV_DIMENSION> BASE_BOUNDARY;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::data_directory;using BASE::solid_body_collection;using BASE_BOUNDARY::Boundary;using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    SOLIDS_STANDARD_TESTS<TV> tests;
    T rho_initial,u_vel_initial,v_vel_initial,p_initial,velocity_uniform;
    T wall_thickness;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int square;

    T rho_left,u_left,v_left,p_left;

    // isobaric fix parameters
    bool woodward_fix_entropy,woodward_fix_enthalpy,isobaric_fix_only_6_cells;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool exact;

    /***************
    example explanation:
    1. block on bottom right
    2. block in middle
    3. rigid square in middle
    4. symmetry test. 4 different velocity fields in 4 quadrants.
    5. fluid moving right from wall
    6. fluid moving left from wall
    7. Double Mach Reflection of a Strong Shock
    ***************/

    WIND_TUNNEL(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),tests(stream_type,data_directory,solid_body_collection),rigid_body_collection(solid_body_collection.rigid_body_collection),
        rho_left(0),u_left(0),v_left(0),p_left(0),woodward_fix_entropy(false),woodward_fix_enthalpy(true),isobaric_fix_only_6_cells(true),
        eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).6),timesplit(false),implicit_rk(false),exact(false)
    {
    }
    
    virtual ~WIND_TUNNEL() {}

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
    if(test_number==1) fluids_parameters.grid->Initialize(TV_INT(3,1)*cells,RANGE<TV>(TV(),TV(3,1)),true);
    else if(test_number==7) fluids_parameters.grid->Initialize(TV_INT(4,1)*cells,RANGE<TV>(TV(),TV(4,1)),true);
    else if(test_number==4) fluids_parameters.grid->Initialize(TV_INT(3,3)*cells,RANGE<TV>(TV(),TV(3,3)));
    else fluids_parameters.grid->Initialize(TV_INT(3,1)*cells,RANGE<TV>(TV(),TV(3,1)));
    if(test_number!=1 && test_number!=7) *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
    if(test_number==4){fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;}
    if(test_number==7){fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;}
    //time
    initial_time=(T)0.;last_frame=1000;frame_rate=(T)80.;
    if(test_number==1) last_frame=320;
    if(test_number==7) last_frame=16;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    rho_initial=(T)1.4;u_vel_initial=(T)3.;v_vel_initial=(T)0.;p_initial=(T)1.;
    if(test_number==6){rho_initial=(T)1.4;u_vel_initial=-(T)3.;v_vel_initial=(T)0.;p_initial=(T)1.;}
    if(test_number==4){rho_initial=1;u_vel_initial=(T)0.;v_vel_initial=(T)0.;p_initial=(T)1.;velocity_uniform=3;}
    if(test_number==7){u_vel_initial=(T)0;rho_left=(T)8;u_left=(T)7.145;v_left=(T)-4.125;p_left=(T)116.5;}
    fluids_parameters.compressible_timesplit=timesplit;

    wall_thickness=(T).1;
    if(timesplit)
        output_directory=STRING_UTILITIES::string_sprintf("Wind_Tunnel/Test_%d__Resolution_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
    else
        output_directory=STRING_UTILITIES::string_sprintf("Wind_Tunnel/Test_%d__Resolution_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
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
    TV velocity_initial(u_vel_initial,v_vel_initial);
    VECTOR<VECTOR<bool,2>,T_GRID::dimension> valid_wall;
    for(int axis=0;axis<T_GRID::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];
    valid_wall[2]=VECTOR<bool,2>::Constant_Vector(false);

    if(test_number==7){
        BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>* boundary_euler=
            new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(rho_left,rho_initial,rho_left,rho_initial),T_FACE_VECTOR(p_left,p_initial,p_left,p_initial),
            TV_FACE_VECTOR(TV(u_left,v_left),velocity_initial,TV(u_left,v_left),velocity_initial),(T).5,valid_wall,true,
            T_FACE_VECTOR(1.,0,1.,0),T_FACE_VECTOR_BOOL(true,true,true,true));
        fluids_parameters.compressible_boundary=new BOUNDARY_MULTIPLE_UNIFORM<TV,TV_DIMENSION>(T_BOUNDARY_FACE_VECTOR(boundary_euler,boundary_euler,this,this));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(rho_initial,rho_initial,rho_initial,rho_initial),T_FACE_VECTOR(p_initial,p_initial,p_initial,p_initial),
            TV_FACE_VECTOR(velocity_initial,velocity_initial,velocity_initial,velocity_initial),(T).5,valid_wall,true,
            T_FACE_VECTOR(1.,0,0,0),T_FACE_VECTOR_BOOL(true,true,false,false));}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& U=fluids_parameters.euler->U;
    TV domain_center=fluids_parameters.euler->grid.Domain().Center();

    T rho,u_vel,v_vel,p;
    //initialize grid variables
    for(CELL_ITERATOR<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        rho=rho_initial;u_vel=u_vel_initial;v_vel=v_vel_initial;p=p_initial;
        if(test_number==1 || test_number==2 || test_number==3 || test_number==5 || test_number==6){rho=rho_initial;u_vel=u_vel_initial;v_vel=v_vel_initial;p=p_initial;}
        else if(test_number==4){
            if(location.x>=domain_center.x && location.y>=domain_center.y){ // top right
                u_vel=velocity_uniform;v_vel=0;} // right going velocity
            else if(location.x<domain_center.x && location.y>=domain_center.y){ // top left
                u_vel=0;v_vel=velocity_uniform;} // up going velocity
            else if(location.x<domain_center.x && location.y<domain_center.y){ // bottom left
                u_vel=-velocity_uniform;v_vel=0;} // left going velocity
            else if(location.x>=domain_center.x && location.y<domain_center.y){ // bottom left
                u_vel=0;v_vel=-velocity_uniform;}} // down going velocity
        else if(test_number==7){
            if(location.x>=(T).1667+location.y*(T).577){
                rho=rho_initial;u_vel=u_vel_initial;v_vel=v_vel_initial;p=p_initial;}
            else{
                rho=rho_left;u_vel=u_left;v_vel=v_left;p=p_left;}}
        EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(U,cell_index,rho,TV(u_vel,v_vel),fluids_parameters.euler->eos->e_From_p_And_rho(p,rho));}
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(test_number==1){
        RIGID_BODY<TV>& rect=tests.Add_Analytic_Box(TV((T)2.5,(T).3));
        rect.Frame().t=TV((T)1.85,(T).05); // left most end .6 from the left of the domain (so translating 1.25+.6). Top most end .2 unit from the bottom (total height=.15+.05(translation)=.2)
        rect.is_static=true;}
    else if(test_number==2){
        RIGID_BODY<TV>& rect=tests.Add_Analytic_Box(TV((T)1.8,(T).3));
        rect.Frame().t=TV((T)1.5,(T).35);
        rect.is_static=true;}
    else if(test_number==3){
        square=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies_2D/square",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(square).t=TV((T)1.25,(T).55);
        rigid_body_collection.Rigid_Body(square).Is_Kinematic()=true;}
    else return;

    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);

    fluids_parameters.euler_solid_fluid_coupling_utilities->solid_state=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho_initial,TV(u_vel_initial,v_vel_initial),
        fluids_parameters.euler->eos->e_From_p_And_rho(p_initial,rho_initial));
}
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Set_Dirichlet_Boundary_Conditions(time);
    EULER_UNIFORM<T_GRID >& euler=*((dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID >&>(fluids_parameters)).euler);
    T_FACE_ARRAYS_BOOL& psi_N=euler.euler_projection.elliptic_solver->psi_N;
    T_FACE_ARRAYS_SCALAR& face_velocities=euler.euler_projection.face_velocities;
    RANGE<TV> domain=fluids_parameters.euler->grid.Domain();
    TV domain_center=domain.Center();

    //Set Neumann condition for the .2 unit hight step
    bool inside=false;
    for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){TV location=iterator.Location();
        if(test_number==5) inside=location.x<=domain_center.x;
        else if(test_number==6) inside=location.x>=domain_center.x;
        else if(test_number==4){inside=false;
            inside=inside || (location.x<=domain_center.x+wall_thickness*(T).5 && location.x>=domain_center.x-wall_thickness*(T).5); // middle horizontal
            inside=inside || (location.y<=domain_center.y+wall_thickness*(T).5 && location.y>=domain_center.y-wall_thickness*(T).5); // middle vertical 
            inside=inside || (location.x>=domain_center.x && location.y>=domain.max_corner.y-wall_thickness); // top right
            inside=inside || (location.y>=domain_center.y && location.x<=domain.min_corner.x+wall_thickness); // top left
            inside=inside || (location.x<=domain_center.x && location.y<=domain.min_corner.y+wall_thickness); // bottom left
            inside=inside || (location.y<=domain_center.y && location.x>=domain.max_corner.x-wall_thickness);} // bottom right 
        if(inside){
            euler.psi(iterator.Cell_Index())=false;
            for(int axis=0;axis<T_GRID::dimension;axis++){
                psi_N.Component(axis)(iterator.First_Face_Index(axis))=true;face_velocities.Component(axis)(iterator.First_Face_Index(axis))=0;
                psi_N.Component(axis)(iterator.Second_Face_Index(axis))=true;face_velocities.Component(axis)(iterator.Second_Face_Index(axis))=0;}}}
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_DIMENSION_SCALAR& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time_input,const int number_of_ghost_cells) const
{
    T time=time_input;
    static T last_time_used=0;  // Hack to get around the fact that fluid_parameters.Write_Output_Files doesn't know what time it is...
    if(time) last_time_used=time; else time=last_time_used;
    const EULER_UNIFORM<T_GRID >& euler=*((const_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID >&>(fluids_parameters)).euler);
    assert(test_number==7);
    if(side==3){
        int axis=(side+1)/2;
        int boundary=Boundary(side,region);
        int reflection_times_two;
        if(grid.Is_MAC_Grid())
            reflection_times_two=2*boundary+(side&1?-1:1);
        else reflection_times_two=2*boundary;
        for(CELL_ITERATOR<TV> iterator(grid,region);iterator.Valid();iterator.Next()){
            const TV_INT cell_index=iterator.Cell_Index();const TV location=iterator.Location();
            if(location.x>=(T).1667){
                TV_INT reflected_node=cell_index;reflected_node[axis]=reflection_times_two-cell_index[axis];
                T rho=u_ghost(reflected_node)(1);
                TV velocity=EULER<T_GRID>::Get_Velocity(u_ghost,reflected_node);velocity(axis)*=-1;
                T e=EULER<T_GRID>::e(u_ghost,reflected_node);
                EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(u_ghost,cell_index,rho,velocity,e);}
            else{EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(u_ghost,cell_index,rho_left,TV(u_left,v_left),euler.eos->e_From_p_And_rho(p_left,rho_left));}}}
    else if(side==4){
        T rho=rho_initial,p=p_initial,u_vel=u_vel_initial,v_vel=v_vel_initial;
        for(CELL_ITERATOR<TV> iterator(grid,region);iterator.Valid();iterator.Next()){
            const TV_INT& cell_index=iterator.Cell_Index();const TV location=iterator.Location();
            if(location.y-sqrt((T)3)*(location.x-(T)1./(T)6) + (T)2*(T)10*time <= 0){rho=rho_initial;u_vel=u_vel_initial;v_vel=v_vel_initial;p=p_initial;}
            else{rho=rho_left;u_vel=u_left;v_vel=v_left;p=p_left;}
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(u_ghost,cell_index,rho,TV(u_vel,v_vel),euler.eos->e_From_p_And_rho(p,rho));}}
    else assert(false); // Just to verify that nothing funky is going on with the new boundary class
}
//#####################################################################
// Function Apply_Isobaric_Fix
//#####################################################################
void Woodward_Collela_Fix(bool fix_entropy,bool fix_enthalpy)
{
    EULER_UNIFORM<T_GRID >& euler=*((dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID >&>(fluids_parameters)).euler);
    EOS_GAMMA<T>* gamma_law=dynamic_cast<EOS_GAMMA<T>*>(euler.eos);  // This isobaric fix depends on the EOS being a gamma law
    static TV_INT reference_cell=fluids_parameters.grid->Closest_Node(TV((T).6-fluids_parameters.grid->dX.x,(T).2-fluids_parameters.grid->dX.y));
    const T rho_a=euler.U(reference_cell)(1);
    const T internal_energy_a=euler.e(euler.U(reference_cell));LOG::cout<<internal_energy_a<<std::endl;
    const T P_a=euler.p(euler.eos,euler.U(reference_cell));  assert(P_a>0);
    const T q_a_2=euler.Get_Velocity(euler.U(reference_cell)).Magnitude_Squared();
    for(int i=0;i<4;i++){TV_INT cell_index=reference_cell+TV_INT(i,1);
        const T internal_energy_b=euler.e(euler.U(cell_index));LOG::cout<<internal_energy_b<<std::endl;
        const T gamma=gamma_law->gamma;const T P_b=euler.p(euler.eos,euler.U(cell_index));  assert(P_b>0);
        const TV q_b=euler.Get_Velocity(euler.U(cell_index));
        T rho_b;
        if(fix_entropy) rho_b=rho_a*pow(P_b/P_a,(T)1/gamma);
        else rho_b=euler.U(cell_index)(1);

        T alpha;
        if(fix_enthalpy) alpha=((T).5*q_a_2+(gamma/(gamma-1))*(P_a/pow(rho_a,gamma))*(pow(rho_a,gamma-1)-pow(rho_b,gamma-1)))/((T).5*q_b.Magnitude_Squared());
        else alpha=(T)1;
        LOG::cout<<"alpha="<<((T).5*q_a_2+(gamma/(gamma-1))*(P_a/pow(rho_a,gamma))*(pow(rho_a,gamma-1)-pow(rho_b,gamma-1)))/((T).5*q_b.Magnitude_Squared())<<std::endl;

        euler.Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,rho_b,sqrt(alpha)*q_b,P_b/(rho_b*(gamma-(T)1)));
        assert(euler.p(euler.eos,euler.U(cell_index))>0);}
    for(int i=0;i<2;i++){TV_INT cell_index=reference_cell+TV_INT(i,2);
        const T internal_energy_b=euler.e(euler.U(cell_index));LOG::cout<<internal_energy_b<<std::endl;
        const T gamma=gamma_law->gamma;const T P_b=euler.p(euler.eos,euler.U(cell_index));  assert(P_b>0);
        const TV q_b=euler.Get_Velocity(euler.U(cell_index));
        T rho_b;
        if(fix_entropy) rho_b=rho_a*pow(P_b/P_a,(T)1/gamma);
        else rho_b=euler.U(cell_index)(1);

        T alpha;
        if(fix_enthalpy) alpha=((T).5*q_a_2+(gamma/(gamma-1))*(P_a/pow(rho_a,gamma))*(pow(rho_a,gamma-1)-pow(rho_b,gamma-1)))/((T).5*q_b.Magnitude_Squared());
        else alpha=(T)1;
        LOG::cout<<"alpha="<<((T).5*q_a_2+(gamma/(gamma-1))*(P_a/pow(rho_a,gamma))*(pow(rho_a,gamma-1)-pow(rho_b,gamma-1)))/((T).5*q_b.Magnitude_Squared())<<std::endl;

        euler.Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,rho_b,sqrt(alpha)*q_b,P_b/(rho_b*(gamma-(T)1)));
        assert(euler.p(euler.eos,euler.U(cell_index))>0);}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after applying woodward & colella isobaric fix",0,1);
}
void Fedkiw_Isobaric_Fix(bool fix_only_6_cells)
{
    EULER_UNIFORM<T_GRID >& euler=*((dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID >&>(fluids_parameters)).euler);
    T_FACE_ARRAYS_BOOL& psi_N=euler.euler_projection.elliptic_solver->psi_N;
    EOS_GAMMA<T>* gamma_law=dynamic_cast<EOS_GAMMA<T>*>(euler.eos);  // This isobaric fix depends on the EOS being a gamma law
    const T gamma=gamma_law->gamma;
    // Could do the one-ring and two-ring calculations somewhere before we begin...
    ARRAY<VECTOR<int,1> ,VECTOR<int,2> > ring(euler.grid.Domain_Indices());ring.Fill(VECTOR<int,1>(-1));
    for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){bool is_one_ring=false;
        for(int axis=0;axis<TV::dimension;axis++) is_one_ring=is_one_ring || psi_N(axis,iterator.First_Face_Index(axis)) || psi_N(axis,iterator.Second_Face_Index(axis));
        if(is_one_ring) ring(iterator.Cell_Index())(1)=1;}
    for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){if(ring(iterator.Cell_Index())(1)==1) break;
        bool is_two_ring=false;for(int axis=0;axis<2*TV::dimension;axis++) is_two_ring=is_two_ring || (ring.Valid_Index(iterator.Cell_Neighbor(axis)) && ring(iterator.Cell_Neighbor(axis))(1)==1);
        if(is_two_ring) ring(iterator.Cell_Index())(1)=2;}

    static TV_INT reference_cell=fluids_parameters.grid->Closest_Node(TV((T).6-fluids_parameters.grid->dX.x,(T).2-fluids_parameters.grid->dX.y));
    // Fix the two-ring
    for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){if(ring(iterator.Cell_Index())(1)!=2) continue;
        const TV_INT& cell_index=iterator.Cell_Index();
        if(fix_only_6_cells){
            TV_INT delta=cell_index-reference_cell;
            if(!(delta.y==2 && (delta.x==1 || delta.x==2))) continue;}
        for(int axis=0;axis<2*TV::dimension;axis++){if(!ring.Valid_Index(iterator.Cell_Neighbor(axis)) || ring(iterator.Cell_Neighbor(axis))(1)!=-1) continue;
            const TV_INT& reference_index=iterator.Cell_Neighbor(axis);
            const T new_rho=euler.U(cell_index)(1)*pow(euler.p(euler.eos,euler.U(reference_index))/euler.p(euler.eos,euler.U(cell_index)),(T)1/gamma);
            euler.Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,new_rho,euler.Get_Velocity(euler.U(cell_index)),
                euler.p(euler.eos,euler.U(cell_index))/(new_rho*(gamma-1)));
            break;}}
     
    // Fix the one-ring
    for(CELL_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){if(ring(iterator.Cell_Index())(1)!=1) continue;
        const TV_INT& cell_index=iterator.Cell_Index();
        if(fix_only_6_cells){
            TV_INT delta=cell_index-reference_cell;
            if(!(delta.y==1 && (delta.x>=1 && delta.x<=4))) continue;}
        for(int axis=0;axis<2*TV::dimension;axis++){if(!ring.Valid_Index(iterator.Cell_Neighbor(axis)) || ring(iterator.Cell_Neighbor(axis))(1)!=2) continue;
            const TV_INT& reference_index=iterator.Cell_Neighbor(axis);
            const T new_rho=euler.U(cell_index)(1)*pow(euler.p(euler.eos,euler.U(reference_index))/euler.p(euler.eos,euler.U(cell_index)),(T)1/gamma);
            euler.Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,new_rho,euler.Get_Velocity(euler.U(cell_index)),
                euler.p(euler.eos,euler.U(cell_index))/(new_rho*(gamma-1)));
            break;}}
     PHYSBAM_DEBUG_WRITE_SUBSTEP("after applying isobaric fix",0,1);
}
void Apply_Isobaric_Fix(const T dt,const T time) PHYSBAM_OVERRIDE
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before applying isobaric fix",0,1);
    if(test_number==1){
        Fedkiw_Isobaric_Fix(isobaric_fix_only_6_cells);
        Woodward_Collela_Fix(woodward_fix_entropy,woodward_fix_enthalpy);}
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
void Log_Parameters() const
{
    LOG::SCOPE scope("WIND_TUNNEL parameters");
    BASE::Log_Parameters();
    LOG::cout<<"isobaric_fix_only_6_cells="<<isobaric_fix_only_6_cells<<std::endl;
    LOG::cout<<"woodward_fix_entropy="<<woodward_fix_entropy<<std::endl;
    LOG::cout<<"woodward_fix_enthalpy="<<woodward_fix_enthalpy<<std::endl;
}
//#####################################################################
};
}
#endif
