//#####################################################################
// Copyright 2007 Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PISTON
//#####################################################################
#ifndef __PISTON__
#define __PISTON__

#include <fstream>
#include <iostream>

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Vectors/VECTOR_UTILITIES.h>
#include <Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
#include <Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class PISTON:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >,CONSERVATION_CALLBACKS<T_input>
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,1> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename COLLISION_BODY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;

public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    T u_initial,rho_initial,T_initial,e_initial,p_initial;
    T piston_initial_position,piston_speed,piston_final_position;
    T boundary_set_time;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>* inaccurate_union;
    int piston;
    int piston_width;
    INTERPOLATION_CURVE<T,TV> motion_curve;
    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool exact;

    /***************
    example explanation:
    1. Sod shock hitting an object on the right.
    2. Piston moving right.
    3. Setting rightward velocity at a fixed point.
    4. Air moving to the left on the left of the piston and to the right on the right of the piston. The piston is in the middle.
    ***************/


    PISTON(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.COMPRESSIBLE),rigid_body_collection(solid_body_collection.rigid_body_collection),
        inaccurate_union(0),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).5),timesplit(false),implicit_rk(false),exact(false)
    {
    }

    virtual ~PISTON() {}

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
    implicit_rk=implicit_rk && !exact;

    //grid
    if(test_number==1 || test_number==4)
        fluids_parameters.grid->Initialize(TV_INT(20*resolution+1),RANGE<TV>::Centered_Box()*5);
    else fluids_parameters.grid->Initialize(TV_INT(20*resolution+1),RANGE<TV>::Unit_Box());
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;
    if(test_number==4) fluids_parameters.domain_walls[0][1]=false;
    //time
    initial_time=(T)0.;last_frame=1500;frame_rate=(T)100.;
    fluids_parameters.cfl=cfl_number;;
    //custom stuff . . .
    fluids_parameters.compressible_eos=new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    //fluids_parameters.compressible_conservation_method=new CONSERVATION_ENO_RF<T_GRID,T_GRID::dimension+2>;
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    //fluids_parameters.compressible_conservation_method->Set_Callbacks(this);
    u_initial=0;rho_initial=1;T_initial=300;e_initial=fluids_parameters.compressible_eos->e_From_T_And_rho(T_initial,rho_initial);p_initial=fluids_parameters.compressible_eos->p(rho_initial,e_initial);
    LOG::cout<<"rho_initial="<<rho_initial<<", T_initial="<<T_initial<<", e_initial="<<e_initial<<", p_initial="<<p_initial<<std::endl;
    LOG::cout<<"sound speed="<<fluids_parameters.compressible_eos->c(rho_initial,e_initial)<<std::endl;
    if(test_number==4){rho_initial=1;u_initial=3;p_initial=(T)1.;}
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    motion_curve.Add_Control_Point(0,TV(piston_initial_position));
    motion_curve.Add_Control_Point(last_frame/frame_rate,TV(piston_final_position));

    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Piston/Test_%d__Resolution_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x));
    else output_directory=STRING_UTILITIES::string_sprintf("Piston/Test_%d__Resolution_%d_explicit",test_number,(fluids_parameters.grid->counts.x));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";

    if(test_number==1){
        piston_initial_position=(T)4.8;piston_speed=(T)0.;}
    else if(test_number==2 || test_number==3){
        piston_initial_position=0;piston_speed=(T).1;}
    else if(test_number==4){
        piston_initial_position=(T)0.;piston_speed=(T)0.;piston_width=1;}
    else{LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}
    piston_final_position=piston_initial_position+piston_speed*(last_frame/frame_rate);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    if(test_number==4){
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,T_FACE_VECTOR((T).125,(T).125),T_FACE_VECTOR((T).1,(T).1),
            TV_FACE_VECTOR(),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}
    else{
        TV velocity_initial=TV::All_Ones_Vector()*u_initial;
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,T_FACE_VECTOR(rho_initial,rho_initial),
            T_FACE_VECTOR(p_initial,p_initial),TV_FACE_VECTOR(test_number==4?-velocity_initial:velocity_initial,velocity_initial),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));}
    inaccurate_union=new FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>(fluids_parameters.euler->grid);
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    if(test_number==1)
        for(int i=0;i<grid.counts.x;i++) {
            T rho=(T)0.,u=(T)0.,p=(T)0.;
            if(grid.X(VECTOR<int,1>(i)).x<0) {rho=(T)1.;p=(T)1.;} else {rho=(T).125;p=(T).1;}
            U(i)(0)=rho; U(i)(1)=rho*u; U(i)(2)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}
    else if(test_number==2 || test_number==3)
        for(int i=0;i<grid.counts.x;i++){
            U(i)(0)=rho_initial;U(i)(1)=rho_initial*u_initial;U(i)(2)=rho_initial*(tmp_eos->e_From_T_And_rho(T_initial,rho_initial)+sqr(u_initial)/(T)2.);}
    else if(test_number==4){
        rho_initial=1;u_initial=3;p_initial=(T)1.;
        for(int i=0;i<grid.counts.x;i++){
            TV_INT piston_face_index=grid.Cell(TV(piston_initial_position),0);
            if(i >=piston_face_index[0]) {rho_initial=(T)1.;p_initial=(T)1.;u_initial=(T)3.0;} else {rho_initial=(T)1.;p_initial=1.;u_initial=-(T)3.0;}
            U(i)(0)=rho_initial;U(i)(1)=rho_initial*u_initial;U(i)(2)=rho_initial*(tmp_eos->e_From_p_And_rho(p_initial,rho_initial)+sqr(u_initial)/(T)2.);}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::Set_Dirichlet_Boundary_Conditions(time);
    EULER_UNIFORM<T_GRID>& euler=*((dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID>&>(fluids_parameters)).euler);
    T_FACE_ARRAYS_BOOL& psi_N=euler.euler_projection.elliptic_solver->psi_N;
    T_FACE_ARRAYS_SCALAR& face_velocities=euler.euler_projection.face_velocities;

   T piston_position=piston_initial_position+piston_speed*time;
   if(test_number==1){TV_INT face_index=euler.grid.Cell(TV(piston_position),0);
       psi_N.Component(0)(face_index)=true;face_velocities.Component(0)(face_index)=piston_speed;}
   else if(test_number==2) for(FACE_ITERATOR<TV> iterator(euler.grid);iterator.Valid();iterator.Next()){
       int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();TV location=iterator.Location();
       if(location.x<=piston_position){psi_N.Component(axis)(face_index)=true;face_velocities.Component(axis)(face_index)=piston_speed;}}
   else if(test_number==4){
       TV_INT piston_index=euler.grid.Cell(TV(piston_position),0);
       for(int i=-piston_width;i<=piston_width;i++){
           TV_INT cell_index(piston_index.x+i);
           euler.psi(cell_index)=false;
           for(int axis=0;axis<T_GRID::dimension;axis++){
               psi_N.Component(axis)(euler.grid.First_Face_Index_In_Cell(axis,cell_index))=true;face_velocities.Component(axis)(euler.grid.First_Face_Index_In_Cell(axis,cell_index))=0;
               psi_N.Component(axis)(euler.grid.Second_Face_Index_In_Cell(axis,cell_index))=true;face_velocities.Component(axis)(euler.grid.Second_Face_Index_In_Cell(axis,cell_index))=0;}}}
   if(test_number==3){TV_INT face_index=euler.grid.Cell(TV(piston_initial_position),0);
       psi_N.Component(0)(face_index)=true;face_velocities.Component(0)(face_index)=piston_speed;}

   boundary_set_time=time;
}
void Get_Neumann_Face_Location(const GRID<TV>& grid_1d,const int face_index,T& location) const
{
    location=piston_initial_position+piston_speed*boundary_set_time;
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{   
    if(test_number==1 || test_number==4) return;
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

    POINT_SIMPLICES_1D<T>& point_simplices=*POINT_SIMPLICES_1D<T>::Create();
    POINT_SIMPLEX_MESH& mesh=point_simplices.mesh;
    DEFORMABLE_PARTICLES<TV>& particles=static_cast<DEFORMABLE_PARTICLES<TV>&>(point_simplices.particles);
    particles.Store_Mass();
    mesh.number_nodes=2;mesh.elements.Exact_Resize(2);
    mesh.elements(0).Set(0);mesh.elements(1).Set(1);
    particles.Add_Elements(mesh.number_nodes);
    particles.X(0)=TV(piston_initial_position-(T)10);
    particles.X(1)=TV(piston_initial_position);
    point_simplices.Update_Point_Simplex_List();
    
    rigid_body->Add_Structure(point_simplices);
    rigid_body->Frame().t=TV(piston_initial_position);
    rigid_body->Twist().linear=TV(piston_speed);
    rigid_body->Is_Kinematic()=true;
    piston=rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);

    inaccurate_union->collision_bodies.Add_Bodies(rigid_body_collection);
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);

    VECTOR<T,T_GRID::dimension+2>& solid_state=fluids_parameters.euler_solid_fluid_coupling_utilities->solid_state;
    EOS_GAMMA<T> *tmp_eos=dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);
    T rho=rho_initial,p=p_initial,u_vel=u_initial;
    solid_state(0)=rho;solid_state(1)=rho*u_vel;solid_state(2)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u_vel)/(T)2.);
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if((test_number==2||test_number==3)&&id==piston) frame.t=motion_curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if((test_number==2||test_number==3)&&id==piston){twist.linear=TV(piston_speed);return true;}
    return false;
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==1 || test_number==4) return;
    T_GRID& grid=fluids_parameters.euler->grid;
    TV velocity=TV(piston_speed);
    T piston_dt_denominator=abs(velocity.x)/grid.dX.x;
    if(piston_dt_denominator>1e-8) dt=min(dt,1/piston_dt_denominator);
}
//#####################################################################
};
}
#endif
