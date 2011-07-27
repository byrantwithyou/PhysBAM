//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
// 1. Deformable sphere on ground
// 2. Deformable torus on ground
// 3. Curtain and ball
// 4. Curtain stretch baseline
// 5. Cloth hanging
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/AXIAL_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    T stiffness_multiplier;
    T damping_multiplier;
    T bending_stiffness_multiplier;
    T bending_damping_multiplier;
    T axial_bending_stiffness_multiplier,axial_bending_damping_multiplier;

    T cloth_cfl;
    int number_side_panels;
    T aspect_ratio;
    T side_length;
    bool fully_explicit;
    
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::frame_rate;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),number_side_panels(40),aspect_ratio((T)1.7),side_length((T)1.0),fully_explicit(false)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Double_Argument("-stiffen_bending",1,"","stiffness multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-dampen_bending",1,"","damping multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-dampen_axial_bending",1,"","axial damping multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-stiffen_axial_bending",1,"","axial stiffness multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-clothcfl",4.,"Cloth CFL");
    parse_args->Add_Integer_Argument("-side_panels",40,"Cloth side panels");
    parse_args->Add_Option_Argument("-fully_explicit","Explicit damping forces");
    parse_args->Add_Option_Argument("-setv","set velocity from positions");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    bending_stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen_bending");
    bending_damping_multiplier=(T)parse_args->Get_Double_Value("-dampen_bending");
    axial_bending_stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen_axial_bending");
    axial_bending_damping_multiplier=(T)parse_args->Get_Double_Value("-dampen_axial_bending");
    cloth_cfl=(T)parse_args->Get_Double_Value("-clothcfl");
    number_side_panels=parse_args->Get_Integer_Value("-side_panels");
    fully_explicit=parse_args->Is_Value_Set("-fully_explicit");
    solids_parameters.set_velocity_from_positions=parse_args->Get_Option_Value("-setv");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Sphere_Test();break;
        case 2: Torus_Test();break;
        case 3: Curtain_And_Ball();break;
        case 4: Curtain_Stretch_Baseline();break;
        case 5: Hanging_Curtain();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    deformable_body_collection.particles.Compute_Auxiliary_Attributes(soft_bindings);
    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // add forces
    switch(test_number){
        case 1:
        case 2:
            Tetrahedralized_Volume_Forces();
            break;
        case 3:{
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            Triangulated_Surface_Forces(true,true,true,linear_stiffness,linear_damping,stiffness_multiplier*8/(1+sqrt((T)2)),damping_multiplier*4);}
            break;
        case 4:
            Triangulated_Surface_Forces(true,false,false,(T).5,(T)0,(T).4,(T).15);
            break;
        case 5:{
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            Triangulated_Surface_Forces(true,true,false,linear_stiffness,linear_damping,stiffness_multiplier*8/(1+sqrt((T)2)),damping_multiplier*4);}
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(test_number!=4 && test_number!=5) solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
    else if(test_number==5) solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,true));

    if(test_number==4) for(int i=1;i<=deformable_body_collection.particles.array_collection->Size();i++) deformable_body_collection.particles.X(i)*=1.2;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
}
//#####################################################################
// Function Sphere_Test
//#####################################################################
void Sphere_Test()
{
    T density=1000;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
    tests.Add_Ground();
}
//#####################################################################
// Function Torus_Test
//#####################################################################
void Torus_Test()
{
    T density=1000;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
    tests.Add_Ground();
}
//#####################################################################
// Function Tetrahedralized_Volume_Forces
//#####################################################################
void Tetrahedralized_Volume_Forces()
{
    for(int i=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        if(tetrahedralized_volume->mesh.elements.m){
            T linear_stiffness=stiffness_multiplier*(T)7.5e3,linear_damping=damping_multiplier*(T).5;
            DEFORMABLES_FORCES<TV>* altitude_springs=Create_Altitude_Springs(*tetrahedralized_volume,(T)linear_stiffness/(1+sqrt((T)2)),(T)linear_damping);
            solid_body_collection.Add_Force(altitude_springs);
            DEFORMABLES_FORCES<TV>* edge_springs=Create_Edge_Springs(*tetrahedralized_volume,(T)linear_stiffness/(1+sqrt((T)2)),(T)linear_damping);
            solid_body_collection.Add_Force(edge_springs);}}
}
//#####################################################################
// Function Curtain_And_Ball
//#####################################################################
void Curtain_And_Ball()
{
    solids_parameters.cfl=(T)50;
    frame_rate=60;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    if(fully_explicit){
        solids_parameters.implicit_solve_parameters.cg_iterations=0;
        solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;}
    else solids_parameters.implicit_solve_parameters.cg_iterations=200;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    tests.Add_Ground();
    RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T).25,(T)0);
    body.X().z=(T).5;
    body.Is_Kinematic()=true;
}
//#####################################################################
// Function Curtain_Stretch_Baseline
//#####################################################################
void Curtain_Stretch_Baseline()
{
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    if(fully_explicit){
        solids_parameters.implicit_solve_parameters.cg_iterations=0;
        solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;}
    else solids_parameters.implicit_solve_parameters.cg_iterations=200;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    
    tests.Create_Cloth_Panel(number_side_panels,side_length,1.0,0);
}
//#####################################################################
// Function Hanging_Curtain
//#####################################################################
void Hanging_Curtain()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    solids_parameters.cfl=(T)50;
    frame_rate=60;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);

    if(test_number==5){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;deformable_body_collection.particles.mass(i+m*(j-1))=FLT_MAX;i=1;j=n;deformable_body_collection.particles.mass(i+m*(j-1))=FLT_MAX;}
}
//#####################################################################
// Function Triangulated_Surface_Forces
//#####################################################################
void Triangulated_Surface_Forces(const bool use_edge,const bool use_bending,const bool use_altitude,const T linear_stiffness,const T linear_damping,const T altitude_stiffness,
    const T altitude_damping)
{
    for(int i=1;TRIANGULATED_SURFACE<T>* triangulated_surface=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        if(use_edge){
            solid_body_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping));} // were *2 and *10
        if(use_bending){
/*            T axial_bending_stiffness=axial_bending_stiffness_multiplier*2/(1+sqrt((T)2)),axial_bending_damping=axial_bending_damping_multiplier*8;
              solid_body_collection.Add_Force(Create_Axial_Bending_Springs(*triangulated_surface,(T).01,axial_bending_stiffness,axial_bending_damping));*/
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));}
        if(use_altitude) solid_body_collection.Add_Force(Create_Altitude_Springs(*triangulated_surface,altitude_stiffness,altitude_damping));}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==3 || test_number==5){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;V(i+m*(j-1))=TV();i=1;j=n;V(i+m*(j-1))=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==3 || test_number==5){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;V(i+m*(j-1))=TV();i=1;j=n;V(i+m*(j-1))=TV();}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==int(2)){
        if(time<2) twist.linear=TV();
        else if(time<(T)3.5) twist.linear=TV(1,(T).5,0);
        else if(time<4) twist.linear=TV(0,(T)-1.5,0);
        else twist.linear=TV(-1.5,0,0);}
    else return false;
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==int(2)){
        if(time<2) frame.t=TV(0,0,(T).5);
        else if(time<(T)3.5) frame.t=TV((time-2),(T).5*(time-2),(T).5);
        else if(time<4) frame.t=TV((T)1.5,(T)((T).75-(T)1.5*(time-(T)3.5)),(T).5);
        else frame.t=TV((T)(1.5-1.5*(time-4)),0,(T).5);}
}
//#####################################################################
// Function Set_Deformable_Particle_Is_Simulated
//#####################################################################
void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
    if(test_number==5){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;particle_is_simulated(i+m*(j-1))=false;i=1;j=n;particle_is_simulated(i+m*(j-1))=false;}
}
//#####################################################################
};
}
#endif
