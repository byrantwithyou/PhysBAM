//#####################################################################
// Copyright 2008, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_WATER
//#####################################################################
//   1. Two bodies with a bend joint and one way coupling
//   2. Two bodies with a bend joint and two way coupling
//   5. Squirrel falling in a wind tunnel
//   6. Driven Bird that alternates minimizing and maximizing drag
//   7. Octosquid that alternates minimizing and maximizing drag
//   8. Driven Bird that alternates minimizing and maximizing drag with alternative wing joints
//#####################################################################
#ifndef __STANDARD_TESTS_WATER__
#define __STANDARD_TESTS_WATER__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_Evolution/SEARCH_CONTROLLER.h>

#define ANGLE_JOINT_TYPE 1
#define POINT_JOINT_TYPE 2

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef GRID<TV> T_GRID;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef VECTOR<int,3> TV_INT;
    typedef typename TV::SPIN T_SPIN;
public:
    T source_velocity_magnitude;
    SOLIDS_FLUIDS_DRIVER<TV>* driver;
    SOLIDS_STANDARD_TESTS<TV> tests;
    ARRAY<int>* referenced_rigid_particles;
    ARRAY<int>* referenced_particles;
    ARRAY<int>* source_rigid_particles;
    ARRAY<int>* source_particles;
    ARRAY<int>* source_elements;
    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;
    ROTATION<TV> rotation;
    SEARCH_CONTROLLER<T_GRID>* controller;
    ARRAY<ARRAY<int> > bone_hierarchy;
    // optimization
    bool perform_optimization,maximize,use_kinematic_motion,use_deformable,use_embedding;
    //fluids
    SMOKE_STANDARD_TESTS_3D<T_GRID> smoke_tests;
    TV source_velocity;
    RANGE<TV> source;
    MATRIX<T,4> world_to_source;
    RIGID_BODY<TV>* octosquid_body;
    ARRAY<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*> deformable_objects_to_simulate;

    bool use_finite_volume,strain_limit,use_implicit;
    JOINT<TV> *left_wing,*right_wing;
    VECTOR<bool,TV::dimension*2> rigid_body_walls;
    GRAVITY<TV> *solids_source;
    bool use_solids_source;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::frame_rate;using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::resolution;using BASE::parse_args;
    using BASE::Add_To_Fluid_Simulation;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;

    STANDARD_TESTS_WATER(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),source_velocity_magnitude(0),tests(*this,solid_body_collection),
        referenced_rigid_particles(0),collision_manager(0),perform_optimization(false),maximize(false),use_kinematic_motion(false),use_deformable(false),use_embedding(false),
        smoke_tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),
        use_finite_volume(false),strain_limit(false),use_implicit(true),left_wing(0),right_wing(0),use_solids_source(false)
    {
    }

    ~STANDARD_TESTS_WATER()
    {}
    
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {
        if(controller && controller->drag_step) dt=controller->dt_hyp;
        }
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {
        if(controller && controller->hypothetical_step) dt=controller->dt_hyp;
        else dt=1/frame_rate+(T)1e-4;}
    //void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {twist(1).angular=T_SPIN();}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-solve_optimization","Perform Golden Section search along gradient direction");
    parse_args->Add_Option_Argument("-maximize","Maximize rather than minimize the chosen objective function");
    parse_args->Add_Option_Argument("-use_deformable","Wrap all rigid bodies in the body with a mesh");
    parse_args->Add_Option_Argument("-use_embedding","Wrap all rigid bodies in the body with a mesh");
    parse_args->Add_Option_Argument("-use_fem","Use Finite Elements");
    parse_args->Add_Option_Argument("-use_source","");
    parse_args->Add_Double_Argument("-velocity",(T)2);
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    smoke_tests.Initialize(Smoke_Test_Number(test_number),resolution);
    perform_optimization=parse_args->Get_Option_Value("-solve_optimization");
    maximize=parse_args->Get_Option_Value("-maximize");
    use_deformable=parse_args->Get_Option_Value("-use_deformable");
    use_embedding=parse_args->Get_Option_Value("-use_embedding");
    use_finite_volume=parse_args->Get_Option_Value("-use_fem");
    use_solids_source=!parse_args->Get_Option_Value("-use_source");
    source_velocity_magnitude=(T)parse_args->Get_Double_Value("-velocity");
    if(source_velocity_magnitude>10) PHYSBAM_FATAL_ERROR();

    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;

    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=50;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
    solids_parameters.cfl=1000;        
    //solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;

    *fluids_parameters.grid=smoke_tests.grid;
    fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.solve_neumann_regions=false;
    fluids_parameters.fluid_affects_solid=false;
    fluids_parameters.density=(T)100000;
    fluids_parameters.incompressible_iterations=200;
//        fluids_parameters.use_vorticity_confinement=true;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Water/Test_%d_velocity_%f_%s%s%s%s",test_number,source_velocity_magnitude,(maximize?"maximize":"minimize"),(use_solids_source?"_no_source":"_source"),(perform_optimization?"_optimize":""),(use_deformable?"_deformable":""));
    frame_rate=60;

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
            last_frame=360;
            fluids_parameters.grid->Initialize(10*resolution+1,10*resolution+1,15*resolution+1,-5,5,0,10,-(T)7.5,(T)7.5);
            fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=false;
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=true;
            break;
        case 5:
            fluids_parameters.grid->Initialize(10*resolution+1,10*resolution+1,10*resolution+1,-5,5,5,15,-5,5);
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][0]=false;
            fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=true;
            break;
        case 6:
        case 8:
            fluids_parameters.grid->Initialize(10*resolution+1,9*resolution+1,10*resolution+1,-5,5,2,11,-5,5);
            fluids_parameters.domain_walls[1][1]=false;
            fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.incompressible_iterations=400;
            break;
        case 7:
            //fluids_parameters.grid->Initialize(10*resolution+1,15*resolution+1,10*resolution+1,-8,8,0,24,-8,8);
            fluids_parameters.grid->Initialize(10*resolution+1,30*resolution+1,10*resolution+1,-8,8,0,48,-8,8);
            fluids_parameters.domain_walls[1][1]=false;
            //fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=false;
            fluids_parameters.incompressible_iterations=400;
            solids_parameters.implicit_solve_parameters.lanczos_iterations=1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=1000;
            solids_parameters.implicit_solve_parameters.cg_restart_iterations=100;
            fluids_parameters.solve_neumann_regions=false;
            fluids_parameters.density=100000;
            solid_body_collection.print_residuals=true;
            fluids_parameters.fluid_affects_solid=true;
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up fluid
    switch(test_number){
        case 1: 
        case 3:
        case 6:
        case 8: Left_Source();break;
        case 2: 
        case 4: Right_Source();break;
        case 5: Left_Source();break;
        case 7: break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    for(int i=0;i<TV::dimension*2;i+=2){
        if(i==0){rigid_body_walls(i)=fluids_parameters.domain_walls[0][0];rigid_body_walls(i+1)=fluids_parameters.domain_walls[0][1];}
        if(i==2){rigid_body_walls(i)=fluids_parameters.domain_walls[1][0];rigid_body_walls(i+1)=fluids_parameters.domain_walls[1][1];}
        if(i==4){rigid_body_walls(i)=fluids_parameters.domain_walls[2][1];rigid_body_walls(i+1)=fluids_parameters.domain_walls[2][0];}}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
void Set_Driver(SOLIDS_FLUIDS_DRIVER<TV>* driver_input)
{
    driver=driver_input;
}
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{  
    T force_magnitude=2;
    //if(((int)time)%4==2) force_magnitude*=1;
    //else if ((int)time%4==0) force_magnitude*=-1;
    if(((int)time)%3==1) force_magnitude*=1;
    else if ((int)time%3==0) force_magnitude*=-1;
    else force_magnitude=0;
    if(solids_source) solids_source->Set_Gravity(TV(0,force_magnitude,0));
}
bool Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,ARRAY<T>& rigid_simplex_forces,TV& orientation,const T time)
{
    T cycle_time=time*4; //.5 second cycle
    if(!solids_source && test_number==7){
        orientation=octosquid_body->Frame().r.Rotate(TV(0,-1,0));
        T force_magnitude=3;
        T alpha=(cycle_time-7*(((int)cycle_time)/7))/4;
        //if(time<1 || ((int)time)%4==3){
        if(((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==2 || ((int)cycle_time)%7==3) force_magnitude*=-1*sin(alpha*(T)pi);
        //else if ((int)time%4==2){
        //else if ((int)time%4==0) force_magnitude*=-1*sin(alpha*(T)pi);
        else force_magnitude=0;
        if(use_deformable) for(int i=0;i<source_elements->m;i++){deformable_simplices.Append((*source_elements)(i));deformable_simplex_forces.Append(force_magnitude);}
        else{PHYSBAM_FATAL_ERROR();}
        return true;
    }
    return false;
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    twist(1).angular=T_SPIN();    
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=source_velocity[iterator.Axis()];
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(test_number!=7) BASE::Get_Source_Velocities(source,world_to_source,source_velocity);
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(controller && !controller->hypothetical_step) controller->Update_Position_Based_State(fluid_collection.incompressible_fluid_collection.face_velocities,dt,time);
    if(controller && controller->hypothetical_step && controller->drag_step) controller->Preprocess_Drag_Substep(time);
}
//#####################################################################
// Function Smoke_Test_Number
//#####################################################################
static int Smoke_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{}
//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    solids_evolution->rigid_body_collisions->collision_manager=collision_manager;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    RIGID_BODY<TV>* kinematic_body=&solid_body_collection.rigid_body_collection.Rigid_Body(id);
    if(!kinematic_body->Is_Kinematic()) return;
    frame=kinematic_body->Frame();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    return false;
}
//#####################################################################
// Function Set_PD_Targets
//#####################################################################
void Set_PD_Targets(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==6){
        T mod_time=min((T)4,(T)max((T)0,(T)(time+1-5*floor((time+1)/5))));
        T angle=(T)pi*3/8-(T)(1-cos((T)pi*(mod_time/2)))*(T)pi*3/8;

        LOG::cout<<"mod_time = "<<mod_time<<", angle = "<<angle<<std::endl;
        left_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(angle,0,0));
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(-angle,0,0));}
    if(test_number==8){
        T mod_time=(min((T)8,(T)max((T)0,(T)(time+2-10*floor((time+2)/10)))))/(T)2;
        T angle=(T)pi*(T).25-(T)(1-cos((T)pi*mod_time/2))*(T)pi*(T).25;

        LOG::cout<<"mod_time = "<<mod_time<<", angle = "<<angle<<std::endl;
        left_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(angle,0,0));
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(-angle,0,0));}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    if(test_number==6 && !controller->hypothetical_step){
        bool old_min=controller->minimize;
        T mod_time=min((T)4,(T)max((T)0,(T)(time+1-5*floor((time+1)/5))));
        controller->minimize=!(mod_time<4);
        controller->real_dx=(mod_time<=(T)2)?(T)2.5*(T)pi/180:(T)5*(T)pi/180;
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
    if(test_number==8 && !controller->hypothetical_step){
        bool old_min=controller->minimize;
        T mod_time=min((T)8,(T)max((T)0,(T)(time+2-10*floor((time+2)/10))));
        controller->minimize=!(mod_time<4);
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
    if(test_number==7 && controller && !controller->hypothetical_step){
        T cycle_time=4*time;
        controller->drag_direction=octosquid_body->Frame().r.Rotate(TV(0,-1,0));
        bool old_min=controller->minimize;
        controller->minimize=((((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==6)?false:true);
        //controller->minimize=((((int)time)%4==0 || ((int)time)%4==1)?false:true);
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==6) joint.joint_function->Set_k_p(750);
            else joint.joint_function->Set_k_p(1500);}
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
}
//#####################################################################
// Function Swap
//#####################################################################
void Swap(bool& val1,bool& val2){
    bool tmp=val1;val1=val2;val2=tmp;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    controller=new SEARCH_CONTROLLER<T_GRID>(solid_body_collection,driver);
    controller->solve_minimization=perform_optimization;
    controller->minimize=!maximize;
    controller->max_iterations=1;
    controller->incorporate_fluids=true;
    controller->use_projection=true;
    controller->dt_per_search_step=(T).1;
    controller->dt_hyp=(T).00333333333333;
    controller->drag_direction=TV(0,0,-1);
    controller->fluids_parameters=&fluids_parameters;
    //set up rigid_bodies and external forces
    switch(test_number){
        case 1: 
        case 2: 
        case 3:
        case 4:
            controller->drag_direction=TV(0,0,(T)1);
            Cylinder_And_Block();
            break;
        case 5:
            controller->drag_direction=TV(0,(T)1,0);
            Squirrel(controller);break;
        case 6:
            controller->minimize=false;
            controller->threshold=(T)1.5e4;
            controller->min_multiplier=(T)1;
            controller->steady_state_drag=false;
            controller->steady_state_drag_time=(T).5;
            controller->use_projection=false;
            controller->real_dx=(T)5*(T)pi/180;
            controller->dx=(T)10*(T)pi/180;
            controller->drag_direction=TV(0,(T)1,0);
            controller->dt_per_search_step=(T).0333333333;
            Driven_Bird();
            break;
        case 8:
            controller->minimize=false;
            controller->threshold=(T)2.5e3;
            controller->min_multiplier=(T)1;
            controller->steady_state_drag=true;
            controller->steady_state_drag_time=(T).2;
            controller->use_projection=true;
            controller->real_dx=(T)5*(T)pi/180;
            controller->dx=(T)20*(T)pi/180;
            controller->drag_direction=TV(0,(T)1,0);
            controller->dt_per_search_step=(T).0333333333;
            Driven_Bird();
            break;
        case 7:
            controller->drag_direction=TV(0,(T)-1,0);
            controller->real_dx=(T)20*(T)pi/180;
            controller->dx=(T)10*(T)pi/180;
            use_deformable=true;use_embedding=false;
            use_kinematic_motion=false;
            controller->dt_per_search_step=(T).1;
            controller->dt_hyp=(T).0333333;
            controller->use_projection=false;
            controller->steady_state_drag=false;
            Octosquid();
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up joints
    switch(test_number){
        case 1:
        case 2: Joints_From_List(ANGLE_JOINT_TYPE);break;
        case 3:
        case 4: Joints_From_List(POINT_JOINT_TYPE);break;
        case 5:
        case 6:
        case 7:
        case 8: break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(use_deformable){
        bool build_tri=false,build_tet=false;
        TETRAHEDRALIZED_VOLUME<T> *volume_original=TETRAHEDRALIZED_VOLUME<T>::Create();
        TRIANGULATED_SURFACE<T> *surface_original=0,*surface=0;
        if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number)) && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number))){
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);
            if(use_embedding){
                surface_original=TRIANGULATED_SURFACE<T>::Create();
                FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}}
        else if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number))){
            if(use_embedding) build_tri=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);}
        else if(use_embedding && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number))){
            surface_original=TRIANGULATED_SURFACE<T>::Create();
            build_tet=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}

        else{if(use_embedding) build_tri=true;build_tet=true;}

        if(surface_original) for(int i=0;i<surface_original->particles.array_collection->Size();i++) surface_original->particles.X(i)+=TV(0,(T)15,0);
        for(int i=0;i<volume_original->particles.array_collection->Size();i++) volume_original->particles.X(i)+=TV(0,(T)15,0);                
        PHYSBAM_ASSERT(!build_tri); // Make baby jesus cry

        if(build_tri || build_tet){
            RANGE<TV> grid_domain;
            T thickness=(T).2;
            for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                grid_domain.Enlarge_To_Include_Box(solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Box());

            bool read_in_phi_from_file=false;
            GRID<TV> new_grid((T).02,grid_domain.Thickened(2*thickness));
            ARRAY<T,VECTOR<int,3> > new_phi(new_grid.Domain_Indices());new_phi.Fill(1e10);
            LEVELSET_IMPLICIT_OBJECT<TV>* implicit=new LEVELSET_IMPLICIT_OBJECT<TV>(new_grid,new_phi);
            if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number))){
                read_in_phi_from_file=true;
                FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number),*implicit);}
            else{
                for(CELL_ITERATOR iterator(new_grid);iterator.Valid();iterator.Next()){const TV_INT &cell_index=iterator.Cell_Index();
                    //new_phi(cell_index)=-1;
                    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                        new_phi(cell_index)=min(new_phi(cell_index),solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Extended_Phi(iterator.Location())-thickness);}
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number),*implicit);}
            if(build_tet){
                int x_edge=40;
                TV edges=new_grid.Domain().Edge_Lengths();
                TV edge_cells_float=TV((T)x_edge,(T)x_edge*edges.y/edges.x,(T)x_edge*edges.z/edges.x);
                TV_INT edge_cells(edge_cells_float);
                volume_original->Initialize_Cube_Mesh_And_Particles(GRID<TV>(edge_cells,new_grid.Domain()));
                volume_original->Discard_Tetrahedrons_Outside_Implicit_Surface(*implicit);
                volume_original->Discard_Valence_Zero_Particles_And_Renumber();
                volume_original->Update_Number_Nodes();
                if(read_in_phi_from_file) for(int i=0;i<volume_original->particles.array_collection->Size();i++) volume_original->particles.X(i)+=TV(0,(T)15,0);                
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);}
            if(build_tri){
                surface_original=DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(implicit->levelset);  
                surface_original=dynamic_cast<TRIANGULATED_SURFACE<T>*>(surface_original->Append_Particles_And_Create_Copy(particles));
                surface_original->Update_Number_Nodes();
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}
            delete implicit;}

        if(surface_original){
            surface=(TRIANGULATED_SURFACE<T>*)surface_original->Append_Particles_And_Create_Copy(particles);
            surface->Update_Number_Nodes();}
        TETRAHEDRALIZED_VOLUME<T>* volume=(TETRAHEDRALIZED_VOLUME<T>*)volume_original->Append_Particles_And_Create_Copy(particles);
        volume->Update_Number_Nodes();
        volume->Initialize_Hierarchy();
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*volume,1,true);
        particles.Store_Velocity();
        if(surface) deformable_body_collection.deformable_geometry.Add_Structure(surface);
        deformable_body_collection.deformable_geometry.Add_Structure(volume);

        TRIANGULATED_SURFACE<T> *soft_bound_surface=0;
        if(use_embedding){
            assert(surface);
            if(test_number!=6 && test_number!=7) soft_bound_surface=TRIANGULATED_SURFACE<T>::Create(particles);
            ARRAY<TRIPLE<int,int,TV> > bindings; 
            if(false && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number))) FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number),bindings);
            else{
                ARRAY<int> tets;const T tolerance=(T)1e-4;
                ARRAY<int> surface_particles;surface->mesh.elements.Flattened().Get_Unique(surface_particles);
                for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
                    tets.Remove_All();bool got_bind=false;
                    T search_region=tolerance/(T)2;while(tets.m==0){search_region*=(T)2;volume->hierarchy->Intersection_List(particles.X(p),tets,search_region);}
                    for(int tt=0;tt<tets.m;tt++){int t=tets(tt);
                        TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(particles.X(p),particles.X.Subset(volume->mesh.elements(t)));
                        bindings.Append(TRIPLE<int,int,TV>(p,t,bary));got_bind=true;break;}
                    if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(TRIPLE<int,int,TV>(p,0,TV(0,0,0)));}}
                if(false) FILE_UTILITIES::Write_To_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number),bindings);}
            ARRAY<int> lookup;lookup.Resize(particles.array_collection->Size());
            for(int i=0;i<bindings.m;i++){
                if(bindings(i).y==0) continue;
                VECTOR<int,4> nodes=volume->mesh.elements(bindings(i).y);
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,bindings(i).x,nodes,bindings(i).z));
                if(test_number!=6 && test_number!=7){
                    int soft_bound_particle=particles.array_collection->Add_Element_From_Deletion_List();
                    particles.X(soft_bound_particle)=particles.X(bindings(i).x);
                    soft_bindings.Add_Binding(VECTOR<int,2>(soft_bound_particle,bindings(i).x),true);
                    lookup(bindings(i).x)=soft_bound_particle;}
            }

            if(test_number!=6 && test_number!=7){
                soft_bound_surface->mesh.elements.Resize(surface->mesh.elements.m);
                soft_bound_surface->mesh.elements.Flattened()=lookup.Subset(surface->mesh.elements.Flattened());}

            surface->Update_Number_Nodes();
            volume->Update_Number_Nodes();}

        // correct mass
        for(int i=0;i<particles.array_collection->Size();i++){particles.mass(i)=fluids_parameters.density/particles.array_collection->Size()/100;}
        binding_list.Distribute_Mass_To_Parents();
        binding_list.Clear_Hard_Bound_Particles(particles.mass);
        particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
        soft_bindings.Set_Mass_From_Effective_Mass();
    
        T stiffness=(T)1e6;
        T damping=(T)10;
        T restlength_clamp=(T)1e-4;
        T cfl_strain_rate=(T).1;

        if(use_finite_volume){
            FINITE_VOLUME<TV,3>* finite_volume=Create_Finite_Volume(*volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1);
            solid_body_collection.Add_Force(finite_volume);}
        else{
            for(int i=0;i<particles.array_collection->Size();i++){particles.mass(i)=fluids_parameters.density/particles.array_collection->Size()/10;}
            particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
            T linear_stiffness=stiffness,linear_damping=damping;
            LINEAR_SPRINGS<TV>* edge_springs;
            edge_springs=Create_Edge_Springs(*volume,linear_stiffness,linear_damping,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
            edge_springs->Clamp_Restlength(restlength_clamp); 
            solid_body_collection.Add_Force(edge_springs);
            T tet_stiffness=stiffness,tet_damping=damping;
            LINEAR_TET_SPRINGS<T>* tet_springs;
            tet_springs=Create_Tet_Springs(*volume,tet_stiffness,tet_damping,false,(T).1,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
            tet_springs->Clamp_Restlength(restlength_clamp); 
            solid_body_collection.Add_Force(tet_springs);}

        // correct mass
        for(int i=0;i<particles.array_collection->Size();i++){particles.mass(i)=fluids_parameters.density/particles.array_collection->Size()/10;}
        binding_list.Distribute_Mass_To_Parents();
        binding_list.Clear_Hard_Bound_Particles(particles.mass);
        particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
        soft_bindings.Set_Mass_From_Effective_Mass();
 
        //if(use_embedding) solid_body_collection.Add_Force(Create_Edge_Binding_Springs(deformable_body_collection.particles,*soft_bindings.binding_mesh,stiffness,damping));

        if(test_number!=8 && test_number !=6){
            if(use_embedding) deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*soft_bound_surface));
            else deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(volume->Get_Boundary_Object()));}

        // binding the deformable particles to the rigid bodies
        ARRAY<int> particle_array;volume->mesh.elements.Flattened().Get_Unique(particle_array);
        source_elements=new ARRAY<int>();
        for(int p=0;p<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Unbound_Particles_In_Rigid_Body(solid_body_collection.rigid_body_collection.Rigid_Body(p),particle_array);
        
        //for(int i=0;i<soft_bound_surface->mesh.elements.m;i++){int j,k,l;soft_bound_surface->mesh.elements(i).Get(j,k,l);
        /*for(int i=0;i<volume->Get_Boundary_Object().mesh.elements.m;i++){int j,k,l;volume->Get_Boundary_Object().mesh.elements(i).Get(j,k,l);
            if(particles.X(j).y<=14.75&&particles.X(k).y<=14.75&&particles.X(l).y<=14.75&&
               abs(particles.X(j).x)<=1&&abs(particles.X(k).x)<=1&&abs(particles.X(l).x)<=1&&
               abs(particles.X(j).z)<=1&&abs(particles.X(k).z)<=1&&abs(particles.X(l).z)<=1
               ){LOG::cout<<"M_DEBUG Found elem "<<i<<" with positions "<<particles.X(j)<<" and "<<particles.X(k)<<" and "<<particles.X(l)<<" and indices "<<j<<" and "<<k<<" and "<<l<<std::endl;source_elements->Append(i);}}
        for(int i=0;i<volume->Get_Boundary_Object().mesh.elements.m;i++){int j,k,l;volume->Get_Boundary_Object().mesh.elements(i).Get(j,k,l);
            if(particles.X(j).y>=16&&particles.X(k).y>=16&&particles.X(l).y>=16&&
               abs(particles.X(j).x)<=1&&abs(particles.X(k).x)<=1&&abs(particles.X(l).x)<=1&&
               abs(particles.X(j).z)<=1&&abs(particles.X(k).z)<=1&&abs(particles.X(l).z)<=1
               ){LOG::cout<<"M_DEBUG Found elem "<<i<<" with positions "<<particles.X(j)<<" and "<<particles.X(k)<<" and "<<particles.X(l)<<" and indices "<<j<<" and "<<k<<" and "<<l<<std::endl;source_elements->Append(i);}}*/
        for(int i=0;i<volume->Get_Boundary_Object().mesh.elements.m;i++){int j,k,l;volume->Get_Boundary_Object().mesh.elements(i).Get(j,k,l);
            source_elements->Append(i);}

        //deformable_body_collection.deformable_geometry.Add_Structure(&volume->Get_Boundary_Object());

        referenced_particles=new ARRAY<int>();source_particles=new ARRAY<int>();
        //for(int i=0;i<particle_array.m;i++){
        for(int i=0;i<particles.array_collection->Size();i++){
            referenced_particles->Append(i);
            for(int j=0;j<source_rigid_particles->m;j++) if(solid_body_collection.rigid_body_collection.Rigid_Body((*source_rigid_particles)(j)).implicit_object->Inside(particles.X(i))) source_particles->Append(i);
        }

        for(int i=0;i<deformable_objects_to_simulate.m;i++){
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& collision_structure=*deformable_objects_to_simulate(i);
            collision_structure.object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(collision_structure);}
    }

    T friction=(T).8;
    if(test_number!=5) tests.Add_Ground(friction,0,0);
    //solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,referenced_rigid_particles));
    controller->bone_hierarchy=bone_hierarchy;
    if(use_solids_source){
        solids_source=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,source_particles,source_rigid_particles);
        solids_source->Set_Gravity(TV(0,0,0));solid_body_collection.Add_Force(solids_source);}

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
    soft_bindings.Set_Mass_From_Effective_Mass();
    ARRAY<int> walls_added;
    for(int i=0;i<TV::dimension*2;i+=2){
        if(i==0){Swap(fluids_parameters.domain_walls[0][0],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[0][1],rigid_body_walls(i+1));}
        if(i==2){Swap(fluids_parameters.domain_walls[1][0],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[1][1],rigid_body_walls(i+1));}
        if(i==4){Swap(fluids_parameters.domain_walls[2][1],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[2][0],rigid_body_walls(i+1));}}
    fluids_parameters.domain_walls[1][0]=false;
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).5,&walls_added);
    if(test_number!=5) fluids_parameters.domain_walls[1][0]=true;
    for(int i=0;i<TV::dimension*2;i+=2){
        if(i==0){Swap(fluids_parameters.domain_walls[0][0],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[0][1],rigid_body_walls(i+1));}
        if(i==2){Swap(fluids_parameters.domain_walls[1][0],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[1][1],rigid_body_walls(i+1));}
        if(i==4){Swap(fluids_parameters.domain_walls[2][1],rigid_body_walls(i));Swap(fluids_parameters.domain_walls[2][0],rigid_body_walls(i+1));}}
    LOG::cout<<"Walls added with ids: "<<walls_added<<std::endl;
}
//#####################################################################
// Function Cylinder_And_Block
//#####################################################################
void Cylinder_And_Block()
{
    //T friction=(T).8;
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("subdivided_box",(T).95,friction);
    kinematic_body.Is_Kinematic()=1;
    kinematic_body.Frame().t=TV(0,(T)1,0);
    kinematic_body.Update_Bounding_Box();

    RIGID_BODY<TV>& dynamic_body=tests.Add_Rigid_Body("cyllink",1,friction);
    dynamic_body.Frame().t=TV(0,(T)4,0);
    dynamic_body.Update_Bounding_Box();

    static const T volume=(T)pi*4;
    dynamic_body.Set_Mass((T)10*fluids_parameters.density*volume);

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    collision_manager->hash.Insert(PAIR<int,int>(kinematic_body.particle_index,dynamic_body.particle_index));
    collision_manager->hash.Insert(PAIR<int,int>(dynamic_body.particle_index,kinematic_body.particle_index));
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    referenced_rigid_particles->Append(dynamic_body.particle_index);

    Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);
    Add_Volumetric_Body_To_Fluid_Simulation(dynamic_body);
    controller->not_affected_by_fluid.Set(dynamic_body.particle_index);
}
//#####################################################################
// Function Squirrel
//#####################################################################
void Squirrel(SEARCH_CONTROLLER<T_GRID>* controller)
{
    T friction=(T)0;
    ARRAY<ARRAY<int> > bone_hierarchy(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()+5);
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("subdivided_box",(T).45,friction);
    kinematic_body.Is_Kinematic()=1;
    kinematic_body.Frame().t=TV(0,(T)10,0);
    kinematic_body.Update_Bounding_Box();
    Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    for(int i=0;i<4;i++){
        const T radius=(T)1;
        const T angle=(((T)i+(T).5)/(T)2)*(T)pi;
        RIGID_BODY<TV>& dynamic_body=tests.Add_Rigid_Body("cyllink",(T).45,friction);
        dynamic_body.Frame().t=TV(radius*cos(angle),(T)8.5,radius*sin(angle));
        dynamic_body.Update_Bounding_Box();

        static const T volume=(T)pi*(T)4;
        if(i!=1) dynamic_body.Is_Kinematic()=1;
        else{
            dynamic_body.Set_Mass((T)10*fluids_parameters.density*volume);
            referenced_rigid_particles->Append(dynamic_body.particle_index);
            collision_manager->hash.Insert(PAIR<int,int>(kinematic_body.particle_index,dynamic_body.particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(dynamic_body.particle_index,kinematic_body.particle_index));
            Create_Joints(kinematic_body.particle_index,dynamic_body.particle_index,POINT_JOINT_TYPE);}
        bone_hierarchy(kinematic_body.particle_index).Append(dynamic_body.particle_index);
        Add_Volumetric_Body_To_Fluid_Simulation(dynamic_body);
        controller->not_affected_by_fluid.Set(kinematic_body.particle_index);}

    controller->bone_hierarchy=bone_hierarchy;
    controller->root_particle_index=kinematic_body.particle_index;
    controller->Create_All_Clusters(collision_manager);
}
//#####################################################################
// Function Driven_Bird
//#####################################################################
void Driven_Bird()
{
    bone_hierarchy.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()+7);
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    arb.constrain_pd_directions=true;

    RIGID_BODY<TV>& kinematic_body=tests.Add_Analytic_Cylinder((T)2,(T).45,10,10);
    if(true){kinematic_body.Is_Kinematic()=1;}
    else{kinematic_body.Set_Mass((T)40*(T).45*(T).45*(T)2*fluids_parameters.density);}
    kinematic_body.Frame().t=TV(0,(T)7,0);
    kinematic_body.Update_Bounding_Box();
    Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    for(int i=0;i<2;i++){
        ROTATION<TV> joint_rotation=ROTATION<TV>::From_Euler_Angles(0,0,(T)(2*i-3)*(T)pi/2);
        FRAME<TV> object_center_to_wing_center(joint_rotation.Rotate(TV(0,(T)1,(T)0)),joint_rotation);
        RIGID_BODY<TV>& driven_wing=tests.Add_Analytic_Box(TV((T).25,(T)1,(T)2));
        RIGID_BODY<TV>& controlled_wing=tests.Add_Analytic_Box(TV((T).25,(T)2.2,(T)2));
        driven_wing.Frame()=kinematic_body.Frame()*object_center_to_wing_center;driven_wing.Update_Bounding_Box();
        controlled_wing.Frame()=FRAME<TV>(joint_rotation.Rotate(TV(0,(T)1.85,0)),ROTATION<TV>())*driven_wing.Frame();controlled_wing.Update_Bounding_Box();

        static const T volume=(T).25*(T).75*(T)2;
        driven_wing.Set_Mass((T)4*fluids_parameters.density*volume);
        Add_Volumetric_Body_To_Fluid_Simulation(driven_wing);bone_hierarchy(kinematic_body.particle_index).Append(driven_wing.particle_index);
        controller->not_affected_by_fluid.Set(driven_wing.particle_index);
        controlled_wing.Set_Mass((T)5.333*fluids_parameters.density*volume);
        Add_Volumetric_Body_To_Fluid_Simulation(controlled_wing);bone_hierarchy(driven_wing.particle_index).Append(controlled_wing.particle_index);
        controller->not_affected_by_fluid.Set(controlled_wing.particle_index);

        referenced_rigid_particles->Append(driven_wing.particle_index);
        referenced_rigid_particles->Append(controlled_wing.particle_index);
        collision_manager->hash.Insert(PAIR<int,int>(controlled_wing.particle_index,driven_wing.particle_index));
        LOG::cout<<"Adding "<<controlled_wing.particle_index<<" and "<<driven_wing.particle_index<<" to collision_manager"<<std::endl;
        collision_manager->hash.Insert(PAIR<int,int>(driven_wing.particle_index,controlled_wing.particle_index));

        ANGLE_JOINT<TV> *driven_joint=new ANGLE_JOINT<TV>;
        ANGLE_JOINT<TV> *controlled_joint=new ANGLE_JOINT<TV>;
        if(i==0) left_wing=driven_joint; else right_wing=driven_joint;
        { // Create the PD-controlled joint that drives the up-down motion
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*driven_joint,arb);
            driven_joint->impulse_accumulator=arb_impulse_accumulator;
            arb.joint_mesh.Add_Articulation(kinematic_body.particle_index,driven_wing.particle_index,driven_joint);
            driven_joint->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0)));
            driven_joint->Set_Child_To_Joint_Frame(FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0))*object_center_to_wing_center);
            JOINT_FUNCTION<TV>* driven_joint_function=arb.Create_Joint_Function(driven_joint->id_number);
            driven_joint->joint_function->active=true;
            driven_joint_function->Set_k_p(60000);driven_joint_function->Set_Target_Angle(ROTATION<TV>());
        }

        { // Create the controlled joint that min/max's drag
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*controlled_joint,arb);
            controlled_joint->impulse_accumulator=arb_impulse_accumulator;
            arb.joint_mesh.Add_Articulation(driven_wing.particle_index,controlled_wing.particle_index,controlled_joint);
            switch(test_number){
                case 6:
                    controlled_joint->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(0,(T)-.625,0),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0)));
                    controlled_joint->Set_Child_To_Joint_Frame(FRAME<TV>(TV(0,(T)1.225,0),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0)));
                    break;
                case 8:
                    controlled_joint->Set_Parent_To_Joint_Frame(FRAME<TV>(TV((T)-.625,0,0),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,(T)-pi/2)));
                    controlled_joint->Set_Child_To_Joint_Frame(FRAME<TV>(TV((T)1.225,0,0),ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,(T)-pi/2)));
                    break;}
            JOINT_FUNCTION<TV>* controlled_joint_function=arb.Create_Joint_Function(controlled_joint->id_number);
            controlled_joint_function->Set_k_p(7000);controlled_joint_function->Set_Target_Angle(ROTATION<TV>());
            if(use_kinematic_motion){controlled_joint->joint_function->active=true;}
            else{
                controller->objective.Resize(controlled_joint->id_number);
                controller->objective(controlled_joint->id_number)=DRAG;
                controlled_joint->global_post_stabilization=false;controlled_joint->joint_function->active=false;
                for(int i=0;i<T_SPIN::dimension;i++){controlled_joint->control_dof(i)=true;}}
            if(test_number==6){
                if(i==0) controlled_joint->Set_Angle_Constraints(true,(T)-pi/2,(T)pi*3/8);
                else controlled_joint->Set_Angle_Constraints(true,(T)-pi*3/8,(T)pi/2);}
            else{
                if(i==0) controlled_joint->Set_Angle_Constraints(true,(T)-pi/2,(T)0);
                else controlled_joint->Set_Angle_Constraints(true,(T)0,(T)pi/2);}}
    }
}
//#####################################################################
// Function Octosquid
//#####################################################################
void Octosquid()
{
    int num_legs=5;
    int num_segments_per_leg=1;
    bone_hierarchy.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()+num_legs*(num_segments_per_leg*2)-num_legs+1);
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    arb.constrain_pd_directions=true;

    T radius=(T).9,length=(T)2.5,width=(T).25;
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("sphere",radius,friction);
    if(use_kinematic_motion){kinematic_body.Is_Kinematic()=1;}
    else{kinematic_body.Set_Mass((T)40*radius*radius*fluids_parameters.density);}
    kinematic_body.Frame().t=TV(0,(T)15.75,0);
    kinematic_body.Update_Bounding_Box();
    octosquid_body=&kinematic_body;
    if(!use_deformable) Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    source_rigid_particles=new ARRAY<int>;
    source_rigid_particles->Append(kinematic_body.particle_index);

    for(int i=0;i<num_legs;i++){
        T offset=radius-(T).1767767;
        //offset+=.2;
        offset-=(T).45;
        RIGID_BODY<TV>* prev_link=&kinematic_body;
        TV direction=ROTATION<TV>::From_Euler_Angles(0,2*(T)pi/num_legs*(i-1),0).Rotate(TV(0,0,1));
        for(int j=0;j<num_segments_per_leg;j++){
            RIGID_BODY<TV>& tail_link=tests.Add_Analytic_Cylinder(length,width);
            tail_link.Frame()=kinematic_body.Frame();tail_link.Frame().t.y-=radius/2+(T).3;tail_link.Frame().t+=(offset+length)*direction;offset+=length+(T).1767767;
            tail_link.Frame().r=ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction);
            FRAME<TV> J;
            RIGID_BODY<TV>* joint_cover=0;
            if(j==0){
                tail_link.Frame().t-=direction*(T).1;
                FRAME<TV> parent_frame=prev_link->Frame();parent_frame.t.y-=(radius/2+(T).3);
                J=FRAME<TV>((tail_link.Frame().t*(T).4+parent_frame.t*(T).6),ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction));}
            else{
                J=FRAME<TV>((tail_link.Frame().t*(T).5+prev_link->Frame().t*(T).5),ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction));
                if(!use_deformable){
                    joint_cover=&tests.Add_Rigid_Body("sphere",(T).1767767,friction); // sqrt(2)*.25/2
                    joint_cover->Frame().t=J.t;}}

            static const T volume=(T).15*length;
            tail_link.Set_Mass((T).4*fluids_parameters.density*volume);
            if(!use_deformable){
                Add_Volumetric_Body_To_Fluid_Simulation(tail_link);
                controller->not_affected_by_fluid.Set(tail_link.particle_index);
                if(j!=1){Add_Volumetric_Body_To_Fluid_Simulation(*joint_cover);controller->not_affected_by_fluid.Set(joint_cover->particle_index);}}
            bone_hierarchy(prev_link->particle_index).Append(tail_link.particle_index);
            
            referenced_rigid_particles->Append(tail_link.particle_index);
            collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,tail_link.particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,prev_link->particle_index));
            if(j!=1 && !use_deformable){
                collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,joint_cover->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,prev_link->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,joint_cover->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,tail_link.particle_index));}

            ANGLE_JOINT<TV> *joint=new ANGLE_JOINT<TV>;
            { // Create the PD-controlled joint that drives the up-down motion
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
                joint->impulse_accumulator=arb_impulse_accumulator;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,tail_link.particle_index,joint);
                joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                joint->Set_Child_To_Joint_Frame(J.Inverse_Times(tail_link.Frame()));
                controller->objective.Resize(joint->id_number);
                controller->objective(joint->id_number)=DRAG;
                for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
                JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
                joint->global_post_stabilization=false;joint->joint_function->active=false;
                joint_function->Set_k_p(750);joint_function->Set_Target_Angle(ROTATION<TV>());
                T angle=ROTATION<TV>::From_Rotated_Vector(direction,TV(0,-1,0)).Angle();
                //joint->Use_Rotation_Constraint(-abs(angle),abs(angle));
                if(angle<0) joint->Set_Angle_Constraints(true,angle,(T)0);
                else joint->Set_Angle_Constraints(true,(T)0,angle);
    //            if(constrain_joints) driven_joint->Use_Rotation_Constraint((T)-pi/4,(T)pi/4); // TODO(jontg): broken in 2D...
            }
            if(j!=1 && !use_deformable) { // Create the joint that keeps our protective sphere in place
                POINT_JOINT<TV> *cover_joint=new POINT_JOINT<TV>;
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*cover_joint,arb);
                cover_joint->impulse_accumulator=arb_impulse_accumulator;cover_joint->joint_function->active=true;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,joint_cover->particle_index,cover_joint);
                cover_joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                cover_joint->Set_Child_To_Joint_Frame(J.Inverse_Times(joint_cover->Frame()));
                JOINT_FUNCTION<TV>* cover_joint_function=arb.Create_Joint_Function(cover_joint->id_number);
                cover_joint_function->Set_k_p(750);cover_joint_function->Set_Target_Angle(ROTATION<TV>());
            }
            prev_link=&tail_link;
        }
    }
}
//#####################################################################
// Function Initialize_Joint_Between
//#####################################################################
void Initialize_Joint_Between(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,TV up)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    TV x=(cf.t-pf.t).Normalized();
    TV y=up.Projected_Orthogonal_To_Unit_Direction(x).Normalized();
    TV z=TV::Cross_Product(x,y);
    FRAME<TV> J((T).5*(cf.t+pf.t),ROTATION<TV>(MATRIX<T,3>(x,y,z)));
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(pf));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(cf));
}
//#####################################################################
// Function Create_Joints_From_Hierarchy
//#####################################################################
void Create_Joints(int parent_id,int child_id,int joint_type)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    JOINT<TV>* joint=0;
    if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
    if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
    for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
    ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
    joint->impulse_accumulator=arb_impulse_accumulator;

    TV up=(arb.rigid_body_collection.Rigid_Body(parent_id).Frame()*TV(0,0,1)+arb.rigid_body_collection.Rigid_Body(child_id).Frame()*TV(0,0,1))/2;
    Initialize_Joint_Between(joint,arb.rigid_body_collection.Rigid_Body(parent_id),arb.rigid_body_collection.Rigid_Body(child_id),up);
    JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
    joint_function->Set_k_p(10000);joint_function->Set_Target_Angle(ROTATION<TV>());
}
//#####################################################################
// Function Joints_From_List
//#####################################################################
void Joints_From_List(int joint_type)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;

    bone_hierarchy.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    controller->root_particle_index=1;
    controller->objective.Resize(JOINT_ID(Value(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size())-1));
    for(int id(1);id<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        JOINT<TV>* joint=0;
        if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
        if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
        for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
        joint->impulse_accumulator=arb_impulse_accumulator;
        bone_hierarchy(id).Append(id+1);
        arb.joint_mesh.Add_Articulation(id,id+1,joint);
        controller->objective(joint->id_number)=FORCE;
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>());
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.frame(id).t-solid_body_collection.rigid_body_collection.rigid_body_particle.frame(id+1).t,solid_body_collection.rigid_body_collection.rigid_body_particle.frame(id+1).r.Inverse()));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(10000);joint_function->Set_Target_Angle(ROTATION<TV>());}
    
}
//#####################################################################
// Function Left_Source
//#####################################################################
void Left_Source()
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
            source=RANGE<TV>((T)-5,(T)5,(T)0,(T)15,(T)-8,(T)-7);
            world_to_source=MATRIX<T,4>::Identity_Matrix();
            source_velocity=TV((T)0,(T)0,source_velocity_magnitude);
            break;
        case 5:
            source=RANGE<TV>((T)-5,(T)5,(T)4.5,(T)5.5,(T)-5,(T)5);
            world_to_source=MATRIX<T,4>::Identity_Matrix();
            source_velocity=TV((T)0,source_velocity_magnitude,(T)0);
            break;
        case 6:
        case 8:
            source=RANGE<TV>((T)-100,(T)100,(T)-.5,(T)2.5,(T)-100,(T)100);
            world_to_source=MATRIX<T,4>::Identity_Matrix();
            source_velocity=TV((T)0,source_velocity_magnitude,(T)0);
            break;}
}
//#####################################################################
// Function Left_Source
//#####################################################################
void Right_Source()
{Left_Source();}
//#####################################################################
};
}
#endif
