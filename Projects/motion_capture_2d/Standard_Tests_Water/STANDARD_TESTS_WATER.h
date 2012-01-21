//#####################################################################
// Copyright 2008, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_WATER
//#####################################################################
//   1. Two bodies with a bend joint and two-way coupling
//   2. Two bodies with a bend joint and two way coupling, very skinny top object
//   3. Driven Bird
//   4. Three-Link Fish
//   5. Controlled Fish
//   6. Octosquid
//#####################################################################
#ifndef __STANDARD_TESTS_WATER__
#define __STANDARD_TESTS_WATER__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_Evolution/SEARCH_CONTROLLER.h>

#define ANGLE_JOINT_TYPE 1
#define POINT_JOINT_TYPE 2

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef GRID<TV> T_GRID;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef VECTOR<int,2> TV_INT;
    typedef typename TV::SPIN T_SPIN;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
public:
    using BASE::Add_To_Fluid_Simulation;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;
    T source_velocity_magnitude;
    SOLIDS_FLUIDS_DRIVER<TV>* driver;
    SOLIDS_STANDARD_TESTS<TV> tests;
    ARRAY<int>* referenced_rigid_particles;
    ARRAY<int>* referenced_particles;
    ARRAY<int>* source_rigid_particles;
    ARRAY<int>* source_particles;
    ARRAY<ARRAY<int> > bone_hierarchy;
    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;
    ROTATION<TV> rotation;
    SEARCH_CONTROLLER<T_GRID>* controller;
    // optimization
    T initial_angle;
    bool maximize;
    //fluids
    SMOKE_STANDARD_TESTS_2D<T_GRID> smoke_tests;
    TV source_velocity;
    RANGE<TV> source1,source2;
    MATRIX<T,3> world_to_source;
    bool use_gravity,use_kinematic_motion,use_deformable,use_embedding;
    RIGID_BODY<TV>* octosquid_body;
    ARRAY<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*> deformable_objects_to_simulate;

    JOINT<TV> *left_wing,*right_wing;
    GRAVITY<TV> *solids_source;
    bool source_hack,use_solids_source;

    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::frame_rate;using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    STANDARD_TESTS_WATER(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),source_velocity_magnitude(0),tests(*this,solid_body_collection),
        referenced_rigid_particles(0),referenced_particles(0),source_rigid_particles(0),source_particles(0),collision_manager(0),initial_angle(0),maximize(false),
        smoke_tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),
        use_gravity(false),use_kinematic_motion(false),use_deformable(false),use_embedding(false),left_wing(0),right_wing(0),solids_source(0),source_hack(false),use_solids_source(false)
    {
    }

    ~STANDARD_TESTS_WATER()
    {}
    
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    //void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    //void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    //void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    //void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    //void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-initial_angle",(T)0);
    parse_args->Add_Option_Argument("-maximize","Maximize rather than minimize the chosen objective function");
    parse_args->Add_Option_Argument("-gravity","Enable gravity");
    parse_args->Add_Option_Argument("-kinematic","Use the kinematically-driven motion only");
    parse_args->Add_Double_Argument("-velocity",(T)5);
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    smoke_tests.Initialize(Smoke_Test_Number(test_number),resolution);
    initial_angle=(T)parse_args->Get_Double_Value("-initial_angle");
    maximize=parse_args->Get_Option_Value("-maximize");
    use_gravity=parse_args->Get_Option_Value("-gravity");
    use_kinematic_motion=parse_args->Get_Option_Value("-kinematic");
    source_velocity_magnitude=parse_args->Get_Double_Value("-velocity");

    source_hack=false;
    use_solids_source=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=20;
    //solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
    solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    solids_parameters.cfl=100;

    *fluids_parameters.grid=smoke_tests.grid;
    fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.fluid_affects_solid=false;
    fluids_parameters.solve_neumann_regions=false;
    fluids_parameters.density=(T)1000;
    LOG::cout<<"source mag "<<source_velocity_magnitude<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Water/Test_%d_velocity_%d_%s_%s_%s",test_number,(int)source_velocity_magnitude,(initial_angle?STRING_UTILITIES::string_sprintf("_initial_angle_%d",initial_angle).c_str():""),(use_gravity?"_gravity":""),(use_kinematic_motion?"_kinematic":""));
    frame_rate=30;

    switch(test_number){
        case 1:
        case 2:
        case 5:
            last_frame=360;
            fluids_parameters.grid->Initialize(20*resolution+1,10*resolution+1,-11,9,0,10);
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=false;
            fluids_parameters.domain_walls[2][2]=true;
            fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.fluid_affects_solid=false;
            break;
        case 3:
        case 4:
        case 6:
            last_frame=360;
            fluids_parameters.fluid_affects_solid=true;
            //fluids_parameters.grid->Initialize(20*resolution+1,120*resolution+1,-4,4,0,48);
            //fluids_parameters.grid->Initialize(20*resolution+1,30*resolution+1,-4,4,0,12);
            fluids_parameters.grid->Initialize(45*resolution+1,60*resolution+1,-9,9,0,24);
            fluids_parameters.domain_walls[2][2]=false;
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.incompressible_iterations=400;
            solids_parameters.implicit_solve_parameters.lanczos_iterations=1000;
            //solids_parameters.implicit_solve_parameters.cg_iterations=1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=1000;
            solids_parameters.implicit_solve_parameters.cg_restart_iterations=50;
            fluids_parameters.solve_neumann_regions=false;
            fluids_parameters.density=100000;
            solid_body_collection.print_residuals=true;
            //solids_parameters.spectral_analysis=true;
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up fluid
    switch(test_number){
        case 1:
        case 2: 
        case 3:
        case 4: 
        case 5:
        case 6: Left_Source();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{  
    T force_magnitude=10;
    //if(time<1 || ((int)time)%4==3){
    if(((int)time)%4==2) force_magnitude*=1;
    //else if ((int)time%4==2){
    else if ((int)time%4==0) force_magnitude*=-1;
    else force_magnitude=0;
    if(/*false && */solids_source) solids_source->Set_Gravity(TV(0,force_magnitude));
}
//#####################################################################
// Function Set_Driver
//#####################################################################
bool Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,ARRAY<T>& rigid_simplex_forces,TV& orientation,const T time)
{
    if(!source_hack && !solids_source && test_number==6){
        T force_magnitude=10;
        //if(time<1 || ((int)time)%4==3){
        if(((int)time)%4==2) force_magnitude*=1;
        //else if ((int)time%4==2){
        else if ((int)time%4==0) force_magnitude*=-1;
        else force_magnitude=0;
        if(use_deformable){
            //for(int i=0;i<301;i++){deformable_simplices.Append(i);deformable_simplex_forces.Append(force_magnitude);}
            deformable_simplices.Append(66);
            deformable_simplices.Append(68);
            deformable_simplices.Append(72);
            deformable_simplices.Append(74);
            deformable_simplices.Append(76);
            deformable_simplices.Append(80);
            deformable_simplices.Append(83);
            deformable_simplices.Append(86);
            deformable_simplices.Append(82);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);
            deformable_simplex_forces.Append(force_magnitude);}
        else{
            rigid_simplices.Append(PAIR<int,int>(int(1),32));
            rigid_simplices.Append(PAIR<int,int>(int(1),31));
            rigid_simplices.Append(PAIR<int,int>(int(1),30));
            rigid_simplices.Append(PAIR<int,int>(int(1),29));
            rigid_simplex_forces.Append(force_magnitude);
            rigid_simplex_forces.Append(force_magnitude);
            rigid_simplex_forces.Append(force_magnitude);
            rigid_simplex_forces.Append(force_magnitude);}
        return true;
    }
    return false;
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(source_hack && use_deformable && test_number==6 && (((int)velocity_time)%4==2||(int)velocity_time%4==0)){
        //for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){V(i)=TV();}
        for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){if(solid_body_collection.deformable_body_collection.particles.X(i).x<(T).8 && solid_body_collection.deformable_body_collection.particles.X(i).x>-(T).8) V(i)=TV();}
        V(2)=TV();
        V(3)=TV();
        V(4)=TV();
        V(5)=TV();
        V(6)=TV();
        V(7)=TV();
        V(8)=TV();
        V(9)=TV();
        V(10)=TV();
        V(11)=TV();
        V(12)=TV();
    }
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(source_hack && !use_deformable && test_number==6 && (((int)velocity_time)%4==2||(int)velocity_time%4==0)) twist(1)=TWIST<TV>();
}
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    T time=velocity_time;
    if(source_hack && use_deformable && test_number==6){
        T force_magnitude=2;
        if(((int)time)%4==2) force_magnitude*=1;
        else if ((int)time%4==0) force_magnitude*=-1;
        else force_magnitude=0;
        TV force_vector=TV(0,force_magnitude);
        //for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){V(i)=force_vector;}
        for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){if(solid_body_collection.deformable_body_collection.particles.X(i).x<(T).8 && solid_body_collection.deformable_body_collection.particles.X(i).x>-(T).8) V(i)=force_vector;}
        V(2)=force_vector;
        V(3)=force_vector;
        V(4)=force_vector;
        V(5)=force_vector;
        V(6)=force_vector;
        V(7)=force_vector;
        V(8)=force_vector;
        V(9)=force_vector;
        V(10)=force_vector;
        V(11)=force_vector;
        V(12)=force_vector;
    }
}
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    T time=velocity_time;
    if(source_hack && !use_deformable && test_number==6){
        T force_magnitude=2;
        if(((int)time)%4==2) force_magnitude*=1;
        else if ((int)time%4==0) force_magnitude*=-1;
        else force_magnitude=0;
        TV force_vector=TV(0,force_magnitude);
        twist(1).linear=force_vector;
        twist(1).angular*=0;
    }
}
//#####################################################################
// Function Set_Driver
//#####################################################################
void Set_Driver(SOLIDS_FLUIDS_DRIVER<TV>* driver_input)
{
    driver=driver_input;
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
        case 6: return 1;
        default: return 1;}
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
    if(test_number==6) return;
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=source_velocity[iterator.Axis()];
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(test_number<6) BASE::Get_Source_Velocities(source1,world_to_source,source_velocity);
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
    if(!controller || controller->hypothetical_step) return;
    if(test_number==3){
        // Rotate a full period every 4 seconds
        T mod_time=min((T)4,(T)max((T)0,(T)(time-5*floor(time/5))));
        T angle=-(T)(1-cos((T)pi*(mod_time/2)))*(3*(T)pi/8);

        left_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(angle));
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(-angle));}
    else if(test_number==4){
        left_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle((T)-cos(time-(T)pi/2)));
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(cos(time)));}
    else if(test_number==5){
        // Rotate a full period every 4 seconds
        T angle=sin((T)pi*(time/2))*((T)pi/6);
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(angle));}
    else if(!controller && test_number==6){
        T angler=((((int)time)%4==0||((int)time)%4==1)?0:-(T)pi/2);
        T anglel=((((int)time)%4==0||((int)time)%4==1)?0:(T)pi/2);
        right_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(angler));
        left_wing->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(anglel));}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==3 && !controller->hypothetical_step){
        bool old_min=controller->minimize;
        T mod_time=min((T)4,(T)max((T)0,(T)(time-5*floor(time/5))));
        controller->minimize=!(mod_time<(T)2);
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
    if(test_number==6 && controller && !controller->hypothetical_step){
        controller->drag_direction=octosquid_body->Rotation().Rotate(TV(0,-1));
        bool old_min=controller->minimize;
        if(((int)time)%4==1||((int)time)%4==2) LOG::cout<<"M1: Max with time"<<time<<std::endl;
        else LOG::cout<<"M1: Min with time"<<time<<std::endl;
        //controller->minimize=((((int)time)%4==1||((int)time)%4==2)?false:true);
        controller->minimize=((((int)time)%4==0||((int)time)%4==1)?false:true);
        //controller->minimize=((((int)time)%4==0)?false:true);
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    controller=new SEARCH_CONTROLLER<T_GRID>(solid_body_collection,driver);
    controller->solve_minimization=false;
    controller->minimize=!maximize;
    controller->max_iterations=1;
    controller->incorporate_fluids=true;
    controller->use_projection=true;
    controller->fluids_parameters=&fluids_parameters;
    //set up rigid_bodies and external forces
    switch(test_number){
        case 1:
            controller->real_dx=5*(T)pi/180;
            controller->dx=(T)10*(T)pi/180;
            controller->steady_state_drag=true;
            controller->drag_direction=TV((T)1,(T)0);
            controller->dt_per_search_step=(T).1;
            controller->threshold=(T)10;
            Cylinder_And_Block();
            break;
        case 2:
            controller->real_dx=5*(T)pi/180;
            controller->dx=(T)10*(T)pi/180;
            controller->steady_state_drag=true;
            controller->drag_direction=TV((T)1,(T)0);
            controller->dt_per_search_step=(T).1;
            controller->threshold=(T)10;
            Cylinder_And_Block();
            break;
        case 3:
            controller->minimize=false;
            controller->use_projection=false;
            controller->min_multiplier=(T)1e-1;
            controller->real_dx=5*(T)pi/180;
            controller->dx=20*(T)pi/180;
            controller->drag_direction=TV((T)0,(T)1);
            controller->threshold=(T)200;
            controller->dt_per_search_step=(T).0333333;
            controller->dt_hyp=(T).00333333;
            Driven_Bird();
            break;
        case 4:
            controller->drag_direction=TV((T)0,(T)1);
            controller->dt_per_search_step=(T).1;
            Three_Link_Fish();
            break;
        case 5:
            controller->steady_state_drag=false;
            controller->use_projection=false;
            controller->min_multiplier=(T)1;
            controller->real_dx=5*(T)pi/180;
            controller->dx=15*(T)pi/180;
            controller->drag_direction=TV((T)1,(T)0);
            controller->threshold=(T).1;
            controller->dt_per_search_step=(T).0333333;
            controller->dt_hyp=(T).0333333;
            Controlled_Fish();
            break;
        case 6:
            // TODO: what with controller?
            controller->drag_direction=TV(0,(T)1);
            controller->use_projection=false;
            controller->min_multiplier=(T)1;
            controller->real_dx=15*(T)pi/180;
            //controller->dx=20*(T)pi/180;
            controller->dx=15*(T)pi/180;
            controller->dt_per_search_step=(T).1;
            controller->dt_hyp=(T).0333333;
            controller->use_projection=false;
            controller->steady_state_drag=false;
            //controller->steady_state_drag_time=(T).1;
            use_deformable=true;
            Octosquid();
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up joints
    switch(test_number){
        case 1: 
        case 2:
            Joints_From_List(ANGLE_JOINT_TYPE);
            break;
        case 3:
        case 4:
        case 5:
        case 6:
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(use_deformable){
        bool build_tri=false,build_tet=false;
        TRIANGULATED_AREA<T>* volume=TRIANGULATED_AREA<T>::Create(particles);
        SEGMENTED_CURVE_2D<T>* surface=0;
        if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.tri2d",test_number)) && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.curve2d",test_number))){
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.tri2d",test_number),*volume);
            if(use_embedding){
                surface=SEGMENTED_CURVE_2D<T>::Create(particles);
                FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.curve2d",test_number),*surface);}}
        else if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.tri2d",test_number))){
            if(use_embedding) build_tri=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.tri2d",test_number),*volume);}
        else if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.curve2d",test_number))){
            build_tet=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.curve2d",test_number),*surface);}
        else{if(use_embedding) build_tri=true;build_tet=true;} 
        if(build_tri || build_tet){
            RANGE<TV> grid_domain;
            T thickness=(T).1;
            for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++) grid_domain.Enlarge_To_Include_Box(solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Box());
            T_GRID new_grid((T).02,grid_domain.Thickened(2*thickness));
            ARRAY<T,VECTOR<int,2> > new_phi(new_grid.Domain_Indices());new_phi.Fill(1e10);
            for(CELL_ITERATOR iterator(new_grid);iterator.Valid();iterator.Next()){const TV_INT &cell_index=iterator.Cell_Index();
                //new_phi(cell_index)=-1;
                for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                    new_phi(cell_index)=min(new_phi(cell_index),solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Extended_Phi(iterator.Location())-thickness);}
            
            LEVELSET_IMPLICIT_OBJECT<TV> implicit(new_grid,new_phi);
            FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.phi2d",test_number),implicit);
            if(build_tet){
                int x_edge=50;
                TV edges=new_grid.Domain().Edge_Lengths();
                TV edge_cells_float=TV((T)x_edge,(T)x_edge*edges.y/edges.x);
                TV_INT edge_cells(edge_cells_float);
                volume->Initialize_Square_Mesh_And_Particles(T_GRID(edge_cells,new_grid.Domain()));
                volume->Discard_Triangles_Outside_Implicit_Curve(implicit);
                volume->Discard_Valence_Zero_Particles_And_Renumber();
                volume->Update_Number_Nodes();
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.tri2d",test_number),*volume);}
            if(build_tri){
                surface=DUALCONTOUR_2D<T>::Create_Segmented_Curve_From_Levelset(implicit.levelset);  
                surface=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(surface->Append_Particles_And_Create_Copy(particles));
                surface->Update_Number_Nodes();
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_2d_%d.curve2d",test_number),*surface);}}

        if(surface) surface->Update_Number_Nodes();
        volume->Update_Number_Nodes();
        volume->Initialize_Hierarchy();
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*volume,1,true);
        particles.Store_Velocity();
        if(surface) deformable_body_collection.deformable_geometry.Add_Structure(surface);
        deformable_body_collection.deformable_geometry.Add_Structure(volume);

        if(use_embedding){
            assert(surface);
            ARRAY<TRIPLE<int,int,VECTOR<T,3> > > bindings; 
            if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/bindings_2d_%d",test_number))) FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_2d_%d",test_number),bindings);
            else{
                ARRAY<int> tets;const T tolerance=(T)1e-4;
                ARRAY<int> surface_particles;surface->mesh.elements.Flattened().Get_Unique(surface_particles);
                for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
                    tets.Remove_All();volume->hierarchy->Intersection_List(particles.X(p),tets,tolerance);bool got_bind=false;
                    for(int tt=0;tt<tets.m;tt++){int t=tets(tt);
                        VECTOR<T,3> bary=TRIANGLE_2D<T>::Barycentric_Coordinates(particles.X(p),particles.X.Subset(volume->mesh.elements(t)));
                        if(bary.x>-tolerance && bary.y>-tolerance && bary.z>-tolerance && bary.x+bary.y+bary.z<(T)1+tolerance){bindings.Append(TRIPLE<int,int,VECTOR<T,3> >(p,t,bary));got_bind=true;break;}}
                    if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(TRIPLE<int,int,VECTOR<T,3> >(p,0,VECTOR<T,3>(0,0,0)));}}
                FILE_UTILITIES::Write_To_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_2d_%d",test_number),bindings);}
            for(int i=0;i<bindings.m;i++){
                if(bindings(i).y==0) continue;
                VECTOR<int,3> nodes=volume->mesh.elements(bindings(i).y);
                binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,bindings(i).x,nodes,bindings(i).z));
                int soft_bound_particle=particles.array_collection->Add_Element_From_Deletion_List();
                soft_bindings.Add_Binding(VECTOR<int,2>(soft_bound_particle,bindings(i).x),true);}
            surface->Update_Number_Nodes();
            volume->Update_Number_Nodes();}

        // correct mass
        for(int i=0;i<particles.array_collection->Size();i++){particles.mass(i)=fluids_parameters.density/particles.array_collection->Size();}
        binding_list.Distribute_Mass_To_Parents();
        binding_list.Clear_Hard_Bound_Particles(particles.mass);
        particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
        soft_bindings.Set_Mass_From_Effective_Mass();

        T stiffness=(T)1e7;
        T damping=(T)50;
        T cfl_strain_rate=(T).1;
        bool strain_limit=false,use_implicit=false;
        //solid_body_collection.Add_Force(Create_Finite_Volume(*volume,new NEO_HOOKEAN<T,2>(stiffness,(T).45,damping,(T).25),true,(T).1));
        solid_body_collection.Add_Force(Create_Edge_Springs(*volume,stiffness,damping,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit));
        //solid_body_collection.Add_Force(Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)damping));        

        if(use_embedding) deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*surface));
        else deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(volume->Get_Boundary_Object()));

        // binding the deformable particles to the rigid bodies
        ARRAY<int> particle_array;volume->mesh.elements.Flattened().Get_Unique(particle_array);
        for(int p=0;p<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Unbound_Particles_In_Rigid_Body(solid_body_collection.rigid_body_collection.Rigid_Body(p),particle_array);
        for(int i=0;i<volume->Get_Boundary_Object().mesh.elements.m;i++){int j,k;volume->Get_Boundary_Object().mesh.elements(i).Get(j,k);
            if(particles.X(j).y<=6.8&&particles.X(k).y<=6.8) LOG::cout<<"M_DEBUG Found elem "<<i<<" with positions "<<particles.X(j)<<" and "<<particles.X(k)<<" and indices "<<j<<" and "<<k<<std::endl;}
        LOG::cout<<"M_DEBUG volume number "<<volume->Get_Boundary_Object().mesh.elements.m<<std::endl;

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
    tests.Add_Ground(friction,0,0);
    if(use_gravity){
        PHYSBAM_FATAL_ERROR();
        //fluids_parameters.gravity=(T)98;
        GRAVITY<TV>* gravity=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,referenced_particles,referenced_rigid_particles);
        gravity->Set_Gravity(TV(0,(T)-10));solid_body_collection.Add_Force(gravity);}
    if(use_solids_source){
        solids_source=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,source_particles,source_rigid_particles);
        solids_source->Set_Gravity(TV(0,0));solid_body_collection.Add_Force(solids_source);}

    controller->bone_hierarchy=bone_hierarchy;
    //controller=0;

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
    soft_bindings.Set_Mass_From_Effective_Mass();

    fluids_parameters.domain_walls[2][1]=false;
    ARRAY<int> walls_added;
    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).5,&walls_added);
    LOG::cout<<"Walls added with ids: "<<walls_added<<std::endl;
    fluids_parameters.domain_walls[2][1]=true;
}
//#####################################################################
// Function Cylinder_And_Block
//#####################################################################
void Cylinder_And_Block()
{
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("square_refined",(T).95,friction);
    kinematic_body.is_static=true;
    kinematic_body.X()=TV(0,(T)1);
    kinematic_body.Update_Bounding_Box();

    T width=0;
    switch(test_number){
        case 1: width=(T)2; break;
        case 2: width=(T).25; break;
        default: width=(T)2;}
    RIGID_BODY<TV>& dynamic_body=tests.Add_Analytic_Box(TV(width,(T)4));
    dynamic_body.Rotation()=ROTATION<TV>::From_Angle((T)pi*initial_angle/180);
    dynamic_body.X()=TV(0,(T)1)+(T)3*dynamic_body.Rotation().Rotate(TV(0,(T)1));
    dynamic_body.Update_Bounding_Box();

    static const T volume=width*(T)4;
    dynamic_body.Set_Mass((T)40*fluids_parameters.density*volume);

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
// Function Driven_Bird
//#####################################################################
void Driven_Bird()
{
    bone_hierarchy.Resize(7);
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

    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("circle",(T).45,friction);
    if(true){kinematic_body.Is_Kinematic()=1;}
    else{kinematic_body.Set_Mass((T)40*(T).45*(T).45*fluids_parameters.density);}
    kinematic_body.X()=TV(0,(T)7);
    kinematic_body.Update_Bounding_Box();
    Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    for(int i=0;i<2;i++){
        ROTATION<TV> joint_rotation=ROTATION<TV>::From_Angle((T)(2*i-3)*(T)pi/8);
        FRAME<TV> object_center_to_wing_center(joint_rotation.Rotate(TV(0,(T).875)),joint_rotation);

        RIGID_BODY<TV>& driven_wing=tests.Add_Analytic_Box(TV((T).25,(T).75));
        RIGID_BODY<TV>& controlled_wing=tests.Add_Analytic_Box(TV((T).25,(T)1));
        RIGID_BODY<TV>& joint_cover=tests.Add_Rigid_Body("circle",(T).1767767,friction); // sqrt(2)*.25/2
        driven_wing.Set_Frame(kinematic_body.Frame()*object_center_to_wing_center);driven_wing.Update_Bounding_Box();
        controlled_wing.Set_Frame(FRAME<TV>(joint_rotation.Rotate(TV(0,(T)1)),ROTATION<TV>())*driven_wing.Frame());controlled_wing.Update_Bounding_Box();
        joint_cover.Set_Frame(kinematic_body.Frame()*object_center_to_wing_center*FRAME<TV>(TV(0,(T).4375),ROTATION<TV>()));joint_cover.Update_Bounding_Box();

        static const T volume=(T).25*(T).75;
        driven_wing.Set_Mass(4*fluids_parameters.density*volume);
        Add_Volumetric_Body_To_Fluid_Simulation(driven_wing);bone_hierarchy(kinematic_body.particle_index).Append(driven_wing.particle_index);
        controller->not_affected_by_fluid.Set(driven_wing.particle_index);
        controlled_wing.Set_Mass((T)5.333*fluids_parameters.density*volume);
        Add_Volumetric_Body_To_Fluid_Simulation(controlled_wing);bone_hierarchy(driven_wing.particle_index).Append(controlled_wing.particle_index);
        controller->not_affected_by_fluid.Set(controlled_wing.particle_index);
        joint_cover.Set_Mass(4*(T)pi*(T).17678*(T).17678*fluids_parameters.density);
        Add_Volumetric_Body_To_Fluid_Simulation(joint_cover);bone_hierarchy(kinematic_body.particle_index).Append(joint_cover.particle_index);
        controller->not_affected_by_fluid.Set(joint_cover.particle_index);
        
        referenced_rigid_particles->Append(driven_wing.particle_index);
        referenced_rigid_particles->Append(controlled_wing.particle_index);
        referenced_rigid_particles->Append(joint_cover.particle_index);
        collision_manager->hash.Insert(PAIR<int,int>(controlled_wing.particle_index,driven_wing.particle_index));
        LOG::cout<<"Adding "<<controlled_wing.particle_index<<" and "<<driven_wing.particle_index<<" to collision_manager"<<std::endl;
        collision_manager->hash.Insert(PAIR<int,int>(driven_wing.particle_index,controlled_wing.particle_index));
        collision_manager->hash.Insert(PAIR<int,int>(joint_cover.particle_index,driven_wing.particle_index));
        collision_manager->hash.Insert(PAIR<int,int>(driven_wing.particle_index,joint_cover.particle_index));
        collision_manager->hash.Insert(PAIR<int,int>(controlled_wing.particle_index,joint_cover.particle_index));
        collision_manager->hash.Insert(PAIR<int,int>(joint_cover.particle_index,controlled_wing.particle_index));

        POINT_JOINT<TV> *driven_joint=new POINT_JOINT<TV>;
        if(i==1) left_wing=driven_joint; else right_wing=driven_joint;
        POINT_JOINT<TV> *controlled_joint=new POINT_JOINT<TV>;
        { // Create the PD-controlled joint that drives the up-down motion
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*driven_joint,arb);
            driven_joint->impulse_accumulator=arb_impulse_accumulator;driven_joint->joint_function->active=true;
            arb.joint_mesh.Add_Articulation(kinematic_body.particle_index,driven_wing.particle_index,driven_joint);
            driven_joint->Set_Parent_To_Joint_Frame(FRAME<TV>());
            driven_joint->Set_Child_To_Joint_Frame(object_center_to_wing_center);
            JOINT_FUNCTION<TV>* driven_joint_function=arb.Create_Joint_Function(driven_joint->id_number);
            driven_joint_function->Set_k_p(10000);driven_joint_function->Set_Target_Angle(ROTATION<TV>());
        }
        { // Create the joint that keeps our protective sphere in place
            POINT_JOINT<TV> *cover_joint=new POINT_JOINT<TV>;
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*cover_joint,arb);
            cover_joint->impulse_accumulator=arb_impulse_accumulator;cover_joint->joint_function->active=true;
            arb.joint_mesh.Add_Articulation(driven_wing.particle_index,joint_cover.particle_index,cover_joint);
            cover_joint->Set_Child_To_Joint_Frame(FRAME<TV>());
            cover_joint->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(0,(T).4375),ROTATION<TV>()).Inverse());
            JOINT_FUNCTION<TV>* cover_joint_function=arb.Create_Joint_Function(cover_joint->id_number);
            cover_joint_function->Set_k_p(1000);cover_joint_function->Set_Target_Angle(ROTATION<TV>());
        }
        { // Create the controlled joint that min/max's drag
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*controlled_joint,arb);
            controlled_joint->impulse_accumulator=arb_impulse_accumulator;
            arb.joint_mesh.Add_Articulation(driven_wing.particle_index,controlled_wing.particle_index,controlled_joint);
            controlled_joint->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(0,(T).4375),ROTATION<TV>()).Inverse());
            controlled_joint->Set_Child_To_Joint_Frame(FRAME<TV>(TV(0,(T)-.5625),ROTATION<TV>()).Inverse());
            JOINT_FUNCTION<TV>* controlled_joint_function=arb.Create_Joint_Function(controlled_joint->id_number);
            controlled_joint_function->Set_k_p(100000);controlled_joint_function->Set_Target_Angle(ROTATION<TV>());
            if(use_kinematic_motion){controlled_joint->joint_function->active=true;}
            else{
                controller->objective.Resize(controlled_joint->id_number);
                controller->objective(controlled_joint->id_number)=DRAG;
                controlled_joint->global_post_stabilization=false;controlled_joint->joint_function->active=false;
                for(int i=0;i<T_SPIN::dimension;i++){controlled_joint->control_dof(i)=true;}}
            if(i==1) controlled_joint->Use_Rotation_Constraint((T)-pi/2,(T)0);
            else controlled_joint->Use_Rotation_Constraint((T)0,(T)pi/2);
        }
    }
}
//#####################################################################
// Function Controlled_Fish
//#####################################################################
void Controlled_Fish()
{
    bone_hierarchy.Resize(7);
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

    T friction=(T)0;
    T radius=(T).9,length=(T)1.5;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("circle",(T)radius,friction);
    if(true){kinematic_body.Is_Kinematic()=1;}
    else{kinematic_body.Set_Mass((T)40*(T)radius*(T)radius*fluids_parameters.density);}
    kinematic_body.X()=TV((T)0,(T)5);
    kinematic_body.Update_Bounding_Box();
    Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    source_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    source_rigid_particles->Append(kinematic_body.particle_index);
    T offset=radius-(T).1767767;
    RIGID_BODY<TV>* prev_link=&kinematic_body;
    for(int i=0;i<4;i++){
        RIGID_BODY<TV>& tail_link=tests.Add_Analytic_Box(TV(length,(T).15));
        tail_link.Set_Frame(kinematic_body.Frame());tail_link.X()+=TV(offset+length,0);offset+=length+(T).1767767;
        FRAME<TV> J;
        RIGID_BODY<TV>* joint_cover=0;
        if(i==1) J=kinematic_body.Frame();
        else{
            J=FRAME<TV>((tail_link.X()*.5+prev_link->X()*.5),ROTATION<TV>());
            joint_cover=&tests.Add_Rigid_Body("circle",(T).1767767,friction); // sqrt(2)*.25/2
            joint_cover->X()=J.t;}

        static const T volume=(T).15*length;
        tail_link.Set_Mass(4*fluids_parameters.density*volume);
        Add_Volumetric_Body_To_Fluid_Simulation(tail_link);
        controller->not_affected_by_fluid.Set(tail_link.particle_index);
        bone_hierarchy(prev_link->particle_index).Append(tail_link.particle_index);
        
        referenced_rigid_particles->Append(tail_link.particle_index);
        collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,tail_link.particle_index));
        collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,prev_link->particle_index));
        if(i!=1){
            collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,joint_cover->particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,prev_link->particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,joint_cover->particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,tail_link.particle_index));}

        POINT_JOINT<TV> *joint=new POINT_JOINT<TV>;
        { // Create the PD-controlled joint that drives the up-down motion
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
            joint->impulse_accumulator=arb_impulse_accumulator;
            arb.joint_mesh.Add_Articulation(prev_link->particle_index,tail_link.particle_index,joint);
            joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
            joint->Set_Child_To_Joint_Frame(J.Inverse_Times(tail_link.Frame()));
            if(i==1 || use_kinematic_motion) joint->joint_function->active=true;
            else{
                controller->objective.Resize(joint->id_number);
                controller->objective(joint->id_number)=DRAG;
                joint->global_post_stabilization=false;joint->joint_function->active=false;
                for(int i=0;i<T_SPIN::dimension;i++){joint->control_dof(i)=true;}}
            JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
            joint_function->Set_k_p(100000);joint_function->Set_Target_Angle(ROTATION<TV>());
        }
        if(i==1) right_wing=joint;
        else { // Create the joint that keeps our protective sphere in place
            POINT_JOINT<TV> *cover_joint=new POINT_JOINT<TV>;
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*cover_joint,arb);
            cover_joint->impulse_accumulator=arb_impulse_accumulator;cover_joint->joint_function->active=true;
            arb.joint_mesh.Add_Articulation(prev_link->particle_index,joint_cover->particle_index,cover_joint);
            cover_joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
            cover_joint->Set_Child_To_Joint_Frame(J.Inverse_Times(joint_cover->Frame()));
            JOINT_FUNCTION<TV>* cover_joint_function=arb.Create_Joint_Function(cover_joint->id_number);
            cover_joint_function->Set_k_p(100000);cover_joint_function->Set_Target_Angle(ROTATION<TV>());
        }
        prev_link=&tail_link;
    }
}
//#####################################################################
// Function Three_Link_Fish
//#####################################################################
void Three_Link_Fish()
{
    bone_hierarchy.Resize(3);
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

    ROTATION<TV> initial_tail_rotation=ROTATION<TV>::From_Angle((T)-1);

    RIGID_BODY<TV> &central_body=tests.Add_Analytic_Box(TV((T).08,(T).8)),
        &head=tests.Add_Analytic_Box(TV((T).08,(T).8)), &tail=tests.Add_Analytic_Box(TV((T).08,(T).8));
    central_body.X()=TV(0,(T)7);head.X()=TV((T)0,(T)8);
    tail.Set_Frame(FRAME<TV>(TV(0,(T)6.5)+initial_tail_rotation.Rotate(TV(0,(T)-.5)),initial_tail_rotation));
    const T mass=(T)100*(T).08*(T).8*fluids_parameters.density;
    LOG::cout<<"Setting mass = "<<mass<<std::endl;
    central_body.Set_Mass(mass);head.Set_Mass(mass);tail.Set_Mass(mass);
    central_body.Update_Bounding_Box();head.Update_Bounding_Box();tail.Update_Bounding_Box();

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(tail.particle_index);
    referenced_rigid_particles->Append(central_body.particle_index);
    referenced_rigid_particles->Append(head.particle_index);

    controller->root_particle_index=tail.particle_index;
    Add_Volumetric_Body_To_Fluid_Simulation(tail);
    Add_Volumetric_Body_To_Fluid_Simulation(central_body);bone_hierarchy(tail.particle_index).Append(central_body.particle_index);
    Add_Volumetric_Body_To_Fluid_Simulation(head);bone_hierarchy(central_body.particle_index).Append(head.particle_index);

    {  // head joint
        left_wing=new POINT_JOINT<TV>;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*left_wing,arb);
        left_wing->impulse_accumulator=arb_impulse_accumulator;
        arb.joint_mesh.Add_Articulation(central_body.particle_index,head.particle_index,left_wing);
        left_wing->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(0,(T)-.5),ROTATION<TV>()));
        left_wing->Set_Child_To_Joint_Frame(FRAME<TV>(TV(0,(T).5),ROTATION<TV>()));
        JOINT_FUNCTION<TV>* head_joint_function=arb.Create_Joint_Function(left_wing->id_number);
        head_joint_function->Set_k_p(10000);head_joint_function->Set_Target_Angle(ROTATION<TV>());
    }
    {  // tail joint
        right_wing=new POINT_JOINT<TV>;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*right_wing,arb);
        right_wing->impulse_accumulator=arb_impulse_accumulator;
        arb.joint_mesh.Add_Articulation(tail.particle_index,central_body.particle_index,right_wing);
        right_wing->Set_Parent_To_Joint_Frame(FRAME<TV>(TV(0,(T)-.5),ROTATION<TV>()));
        right_wing->Set_Child_To_Joint_Frame(FRAME<TV>(TV(0,(T).5),ROTATION<TV>()));
        JOINT_FUNCTION<TV>* tail_joint_function=arb.Create_Joint_Function(right_wing->id_number);
        tail_joint_function->Set_k_p(10000);tail_joint_function->Set_Target_Angle(initial_tail_rotation.Inverse());
    }
    if(use_kinematic_motion) {left_wing->joint_function->active=true;right_wing->joint_function->active=true;}
    else{
        left_wing->global_post_stabilization=false;left_wing->joint_function->active=false;
        right_wing->global_post_stabilization=false;right_wing->joint_function->active=false;
        for(int i=0;i<T_SPIN::dimension;i++){left_wing->control_dof(i)=true;right_wing->control_dof(i)=true;}}
}
//#####################################################################
// Function Octosquid
//#####################################################################
void Octosquid()
{
    int num_legs=2;
    int num_segments_per_leg=1;
    bone_hierarchy.Resize(num_legs*(num_segments_per_leg*2)-num_legs+1);
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
    arb.constrain_pd_directions=false;

    T radius=(T).9,length=3,width=(T).3;
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("circle",radius,friction);
    if(use_kinematic_motion){kinematic_body.Is_Kinematic()=1;}
    else{kinematic_body.Set_Mass((T)40*radius*radius*fluids_parameters.density);}
    kinematic_body.X()=TV(0,(T)8);
    kinematic_body.Update_Bounding_Box();
    octosquid_body=&kinematic_body;
    if(!use_deformable){Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);controller->not_affected_by_fluid.Set(kinematic_body.particle_index);}
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    source_rigid_particles=new ARRAY<int>;
    source_rigid_particles->Append(kinematic_body.particle_index);

    for(int i=0;i<num_legs;i++){
        T offset=radius-(T).1767767;
        //offset+=(T).2;
        offset-=(T).4;
        RIGID_BODY<TV>* prev_link=&kinematic_body;
        TV direction=ROTATION<TV>::From_Angle(i==1?-(T)pi:0).Rotate(TV(1,0));
        for(int j=0;j<num_segments_per_leg;j++){
            RIGID_BODY<TV>& tail_link=tests.Add_Analytic_Box(TV(length,width));
            tail_link.Set_Frame(kinematic_body.Frame());tail_link.X().y-=radius/2;tail_link.X()+=(offset+length)*direction;offset+=length+.1767767;
            tail_link.Rotation()=ROTATION<TV>::From_Rotated_Vector(TV(1,0),direction); // TODO: hmm
            FRAME<TV> J;
            RIGID_BODY<TV>* joint_cover=0;
            if(j==1){
                tail_link.X()-=direction*.1;
                FRAME<TV> parent_frame=prev_link->Frame();parent_frame.t.y-=radius/2;
                J=FRAME<TV>((tail_link.X()*.6+parent_frame.t*.4),ROTATION<TV>());}
            else{
                J=FRAME<TV>((tail_link.X()*.5+prev_link->X()*.5),ROTATION<TV>());
                if(!use_deformable){
                    joint_cover=&tests.Add_Rigid_Body("circle",(T).1767767,friction); // sqrt(2)*.25/2
                    joint_cover->X()=J.t;}}

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

            POINT_JOINT<TV> *joint=new POINT_JOINT<TV>;
            if(i==1 && j==1) left_wing=joint;
            if(i==2 && j==1) right_wing=joint;
            { // Create the PD-controlled joint that drives the up-down motion
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
                joint->impulse_accumulator=arb_impulse_accumulator;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,tail_link.particle_index,joint);
                joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                joint->Set_Child_To_Joint_Frame(J.Inverse_Times(tail_link.Frame()));
                controller->objective.Resize(joint->id_number);
                controller->objective(joint->id_number)=DRAG;
                joint->global_post_stabilization=false;joint->joint_function->active=false;
                for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
                JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
                //joint_function->Set_k_p(10000);joint_function->Set_Target_Angle(ROTATION<TV>::From_Rotated_Vector(direction,TV(0,-1)));
                joint_function->Set_k_p(10000);joint_function->Set_Target_Angle(ROTATION<TV>());
                T angle=ROTATION<TV>::From_Rotated_Vector(direction,TV(0,-1)).Angle();
                //joint->Use_Rotation_Constraint(-abs(angle),abs(angle));
                if(angle<0) joint->Use_Rotation_Constraint(angle,(T)0);
                else joint->Use_Rotation_Constraint((T)0,angle);
            }
            if(j!=1 && !use_deformable) { // Create the joint that keeps our protective sphere in place
                POINT_JOINT<TV> *cover_joint=new POINT_JOINT<TV>;
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*cover_joint,arb);
                cover_joint->impulse_accumulator=arb_impulse_accumulator;cover_joint->joint_function->active=true;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,joint_cover->particle_index,cover_joint);
                cover_joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                cover_joint->Set_Child_To_Joint_Frame(J.Inverse_Times(joint_cover->Frame()));
                JOINT_FUNCTION<TV>* cover_joint_function=arb.Create_Joint_Function(cover_joint->id_number);
                cover_joint_function->Set_k_p(10000);cover_joint_function->Set_Target_Angle(ROTATION<TV>());
            }
            prev_link=&tail_link;
        }
    }
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
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>());
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.rotation(id+1).Inverse_Rotate(solid_body_collection.rigid_body_collection.rigid_body_particle.X(id)-solid_body_collection.rigid_body_collection.rigid_body_particle.X(id+1)),
                                                  solid_body_collection.rigid_body_collection.rigid_body_particle.rotation(id+1).Inverse()));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint->global_post_stabilization=false;joint->joint_function->active=false;
        joint_function->Set_k_p(1000);joint_function->Set_Target_Angle(ROTATION<TV>());
        controller->objective(joint->id_number)=DRAG;}
}
//#####################################################################
// Function Left_Source
//#####################################################################
void Left_Source()
{
    switch(test_number){
        case 1:
        case 2:
        case 5:
            source1=RANGE<TV>((T)-11.5,(T)-10.5,(T)0,(T)10);
            source2=RANGE<TV>((T)7,(T)8,(T)0,(T)10);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=TV(source_velocity_magnitude,(T)0);
            break;
        case 3:
        case 4:
        case 6:
            source1=RANGE<TV>((T)-4,(T)4,(T)-.5,(T).5);
            source2=RANGE<TV>((T)-4,(T)4,(T)11.5,(T)12.5);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=TV((T)0,source_velocity_magnitude);
            break;
        default: break;}
}
//#####################################################################
// Function Left_Source
//#####################################################################
void Right_Source()
{Left_Source();}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    if(controller) controller->Write(stream_type,output_directory,frame);
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame)
{
    BASE::Read_Output_Files_Solids(frame);
    if(controller) controller->Read(stream_type,output_directory,frame);
}
//#####################################################################
};
}
#endif
