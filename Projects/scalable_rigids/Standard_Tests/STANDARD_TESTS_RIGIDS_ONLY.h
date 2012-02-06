//#####################################################################
// Copyright 2010, Elliot English, Jon Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_RIGIDS_ONLY
//#####################################################################
// 1. Rings on pegs (concave, simple geometry)
// 2. Bones on pegs (concave, complex geometry) complex geometry will stress functions like particles_in_levelset
// 3. Spheres in a bin (analytics, large scale)
// 4. Drop three balls with different coefficients of restitution
// 5. Stack bricks forming a wall
// 6. Stack bricks forming a wall (modified)
// 7. Cubes, spheres, and rings in a bin
// 8. Bones in a bin
// 9. Push Out
// 10. Stack of bricks angled so that only one edge is in contact between each pair
//#####################################################################
#ifndef __STANDARD_TESTS_RIGIDS_ONLY__
#define __STANDARD_TESTS_RIGIDS_ONLY__
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_RIGIDS_ONLY:public RIGIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    RIGIDS_STANDARD_TESTS<TV> tests;

    int parameter;
    int number_of_rigids;

    //flags for test 7
    bool box_flag;
    bool sphere_flag;
    bool ring_flag;

    typedef RIGIDS_EXAMPLE<TV> BASE;
    using BASE::rigids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::rigid_body_collection;using BASE::rigids_evolution;using BASE::parse_args;using BASE::test_number;
    using BASE::frame_rate;

    STANDARD_TESTS_RIGIDS_ONLY(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,rigid_body_collection), box_flag(true),sphere_flag(true),ring_flag(true)
    {
    }

    ~STANDARD_TESTS_RIGIDS_ONLY()
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
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Integer_Argument("-number_of_rigids",100,"number of rigid bodies in simulation");
    parse_args->Add_Option_Argument("-noanalytic","disable analytic collisions");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");

    parse_args->Add_Option_Argument("-disable_box","disable boxes for test 7");
    parse_args->Add_Option_Argument("-disable_sphere","disable spheres for test 7");
    parse_args->Add_Option_Argument("-disable_ring","disable rings for test 7");

    parse_args->Add_Option_Argument("-use_pgs","use projected-gauss-seidel for solving contact");
    parse_args->Add_Option_Argument("-disable_sp","disable shock propagation");
    parse_args->Add_Option_Argument("-disable_po","disable push out");
    //parse_args->Add_String_Argument("-set_contact_residuals_filename","residuals","set filename used in SOLVE_CONTACT to store residuals");

}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);

    parameter=parse_args->Get_Integer_Value("-parameter");
    number_of_rigids=parse_args->Get_Integer_Value("-number_of_rigids");
    if(parameter) output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);
    rigids_parameters.rigid_body_collision_parameters.use_analytic_collisions=!parse_args->Get_Option_Value("-noanalytic");
    rigid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    rigids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel=parse_args->Is_Value_Set("-use_pgs");
    rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=!parse_args->Is_Value_Set("-disable_sp");
    rigids_parameters.rigid_body_collision_parameters.use_push_out=!parse_args->Is_Value_Set("-disable_po");
    if(rigids_parameters.rigid_body_collision_parameters.use_push_out)
    {
        if(rigids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel)
            rigids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel_push_out=true;
        else
            rigids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
    }
//rigids_parameters.rigid_body_collision_parameters.contact_residuals_filename=parse_args->Get_String_Value("-set_contact_residuals_filename");

    box_flag = !parse_args->Get_Option_Value("-disable_box");
    sphere_flag = !parse_args->Get_Option_Value("-disable_sphere");
    ring_flag = !parse_args->Get_Option_Value("-disable_ring");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    rigids_parameters.cfl=1;
    rigids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;

    //rigids_parameters.rigid_body_collision_parameters.perform_collisions=false;
    rigids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
    rigids_parameters.rigid_body_collision_parameters.contact_proximity=9.8/(frame_rate*frame_rate);

    switch(test_number){
        case 1: Ring_Test();break;
        case 2: Bone_Test();break;
        case 3: Sphere_Test();break;
        case 4: Bounce(0);break;
        case 5: Bricks();break;
        case 6: Brick_Wall(number_of_rigids);break;
        case 7: Wall_Of_Rigids(number_of_rigids,box_flag,sphere_flag,ring_flag);break;
        case 8: Wall_Of_Bones(number_of_rigids);break;
        case 9: Push_Out();break;
        case 10: Angled_Bricks();break;

        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,true));
}
//#####################################################################
// Function Ring_Test
//#####################################################################
void Ring_Test()
{
    last_frame=720;
    T mu=(T)0.6;
    T epsilon=(T)0.3;
    int poles=5;

    LOG::cout<<"Parameter 1: running 200 rings"<<std::endl;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,200,TV(0,30,0));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    Random_Scene_Generator("Rings_Test/ring_revolve",400,11111,random_placement,rigid_body_collection,tests);
    
    for(int i=0;i<rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Mass(10);
        rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,1,3));}

    for(int i=0;i<poles;i++)
        for(int j=0;j<poles;j++){
            // Poles
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,mu);
            rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("pole %d %d",i,j));
            rigid_body.is_static=true;
            rigid_body.X()=TV((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);}

    tests.Add_Ground(mu);
}
//#####################################################################
// Function Bone_Test
//#####################################################################
void Bone_Test()
{
    last_frame=720;
    T mu=(T)0.6;
    T epsilon=(T)0.3;
    int poles=5;

    LOG::cout<<"Parameter 1: running 200 bones"<<std::endl;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,200,TV(0,30,0));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(3);
    Random_Scene_Generator("Bones/Cranium_Resized",400,11111,random_placement,rigid_body_collection,tests);
    
    for(int i=0;i<rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Mass(10);
        rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,1,3));}

    for(int i=0;i<poles;i++)
        for(int j=0;j<poles;j++){
            // Poles
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,mu);
            rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("pole %d %d",i,j));
            rigid_body.is_static=true;
            rigid_body.X()=TV((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);}

    tests.Add_Ground(mu);
}
//#####################################################################
// Function Sphere_Test
//#####################################################################
void Sphere_Test()
{
    last_frame=720;

    TV lower_corner(0,0,0);
    TV upper_corner(10,10,10);
    TV center=(T)0.5*(upper_corner+lower_corner);
    TV size=(upper_corner-lower_corner);
    T coefficient_of_restitution=0;
    T coefficient_of_friction=0;

    STREAM_TYPE stream_type((float)0);

    RIGID_BODY<TV>& left_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    left_wall.X()=VECTOR<T,3>(lower_corner.x,center.y,center.z);
    left_wall.Rotation()=ROTATION<VECTOR<T,3> >(-(T).5*(T)pi,VECTOR<T,3>(0,0,1));
    left_wall.Set_Name("left wall");

    RIGID_BODY<TV>& right_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    right_wall.X()=VECTOR<T,3>(upper_corner.x,center.y,center.z);
    right_wall.Rotation()=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(0,0,1));
    right_wall.Set_Name("right wall");

    RIGID_BODY<TV>& bottom_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.z));
    bottom_wall.X()=VECTOR<T,3>(center.x,lower_corner.y,center.z);
    bottom_wall.Rotation()=ROTATION<TV>(0,TV(0,0,1));
    bottom_wall.Set_Name("bottom wall");

    RIGID_BODY<TV>& front_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    front_wall.X()=VECTOR<T,3>(center.x,center.y,lower_corner.z);
    front_wall.Rotation()=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(1,0,0));
    front_wall.Set_Name("front wall");

    RIGID_BODY<TV>& back_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    back_wall.X()=VECTOR<T,3>(center.x,center.y,upper_corner.z);
    back_wall.Rotation()=ROTATION<VECTOR<T,3> >((T)-.5*(T)pi,VECTOR<T,3>(1,0,0));
    back_wall.Set_Name("back wall");

    //char* object_names[]={"subdivided_box","sphere","New_Bones/Cranium_1","New_Bones/Pelvis_1","New_Bones/Right_Hand"};//,"dumbbell","bowl","spoon"};
    //T object_scales[]={1,1,30,30,30};

    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(min(size(1),size(3))/2.5,30,TV(center(1),1,center(3)));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(0.5);
    Random_Scene_Generator("sphere",500,11111,random_placement,rigid_body_collection,tests);
}
//#####################################################################
// Function Bounce
//#####################################################################
void Bounce(const T angle)
{
    last_frame=240;

    T x_pos[]={-3,0,3},cor[]={(T)1.0,(T).5,0};
    for(int i=0;i<3;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",1,(T).5);
        rigid_body.X()=TV(x_pos[i],5,0);rigid_body.Set_Coefficient_Of_Restitution(cor[i]);
        rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("sphere (cor %g)",cor[i]));}

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",1,(T).5);
    rigid_body.X()=TV(x_pos[2],8,0);rigid_body.Set_Coefficient_Of_Restitution(0);
    rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("sphere"));

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,0,1,1);
    ground.Rotation()=ROTATION<TV>(angle,TV(0,0,1));
}
//#####################################################################
// Function Bricks
//#####################################################################
void Bricks()
{
    last_frame=240;
    //rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
    //rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
    //rigids_parameters.use_trapezoidal_rule_for_velocities=true;

    int width=1;
    int height=parameter;
    T offset=0.2;

    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,1);
            rigid_body.X()=TV(2.0*x+offset*(y%2),2.0*y-1.0,0);
            rigid_body.Set_Coefficient_Of_Restitution(0);
            
        }

    tests.Add_Ground((T).5,0,1,1);
}
//#####################################################################
// Function Brick_Wall
//#####################################################################
void Brick_Wall(int number)
{
    last_frame=500;

    int width=(int)sqrt(number)+1;
    int height=width;
    T offset=1;

    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,1);
            rigid_body.X()=TV(2.0*x+offset*(y%2),2.0*y-1.0,0);
            rigid_body.Set_Coefficient_Of_Restitution(0);
        }

    tests.Add_Ground((T).5,0,1,1);

}
//#####################################################################
// Function Construct_Wall
//#####################################################################
void Construct_Wall(TV lower_corner, TV upper_corner, TV center, TV size, T coefficient_of_restitution, T coefficient_of_friction)
{
    STREAM_TYPE stream_type((float)0);

    RIGID_BODY<TV>& left_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    left_wall.X()=VECTOR<T,3>(lower_corner.x,center.y,center.z);
    left_wall.Rotation()=ROTATION<VECTOR<T,3> >(-(T).5*(T)pi,VECTOR<T,3>(0,0,1));
    left_wall.Set_Name("left wall");

    RIGID_BODY<TV>& right_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    right_wall.X()=VECTOR<T,3>(upper_corner.x,center.y,center.z);
    right_wall.Rotation()=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(0,0,1));
    right_wall.Set_Name("right wall");

    RIGID_BODY<TV>& bottom_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.z));
    bottom_wall.X()=VECTOR<T,3>(center.x,lower_corner.y,center.z);
    bottom_wall.Rotation()=ROTATION<TV>(0,TV(0,0,1));
    bottom_wall.Set_Name("bottom wall");

    RIGID_BODY<TV>& front_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    front_wall.X()=VECTOR<T,3>(center.x,center.y,lower_corner.z);
    front_wall.Rotation()=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(1,0,0));
    front_wall.Set_Name("front wall");

    RIGID_BODY<TV>& back_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    back_wall.X()=VECTOR<T,3>(center.x,center.y,upper_corner.z);
    back_wall.Rotation()=ROTATION<VECTOR<T,3> >((T)-.5*(T)pi,VECTOR<T,3>(1,0,0));
    back_wall.Set_Name("back wall");
}

//#####################################################################
// Function Wall_Of_Rigids
//#####################################################################
void Wall_Of_Rigids(int number, bool useBox, bool useSphere, bool useRing)
{
    last_frame=400;

    TV lower_corner(0,0,0);
    TV upper_corner(20,20,20);
    TV center=(T)0.5*(upper_corner+lower_corner);
    TV size=(upper_corner-lower_corner);
    T coefficient_of_restitution=0;
    T coefficient_of_friction=0;

    Construct_Wall(lower_corner, upper_corner, center, size, coefficient_of_restitution, coefficient_of_friction);

    int height = (int)(number/10) + 10;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(min(size(1),size(3))/2.5,height,TV(center(1),1,center(3)));

    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(0.5);

    int split = ((useBox?1:0) + (useSphere?1:0) + (useRing?1:0));
    if (split==0)split=1;
    int count = number/split + number%split;

    if (useBox) Random_Scene_Generator("subdivided_box",count,11111,random_placement,rigid_body_collection,tests);
    if (useSphere) Random_Scene_Generator("sphere",count,11111,random_placement,rigid_body_collection,tests);
    if (useRing) Random_Scene_Generator("Rings_Test/ring_revolve",count,11111,random_placement,rigid_body_collection,tests);

}

//#####################################################################
// Function Wall_Of_Bones
//#####################################################################
void Wall_Of_Bones(int number)
{
    last_frame=720;

    TV lower_corner(0,0,0);
    TV upper_corner(20,20,20);
    TV center=(T)0.5*(upper_corner+lower_corner);
    TV size=(upper_corner-lower_corner);
    T coefficient_of_restitution=0;
    T coefficient_of_friction=0;

    Construct_Wall(lower_corner, upper_corner, center, size, coefficient_of_restitution, coefficient_of_friction);

    LOG::cout<<"Test 8: running "<< number <<" bones in enclosed space"<<std::endl;

    int height = (int)(number/10)+10;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(min(size(1),size(3))/2.5,height,TV(center(1),1,center(3)));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(30);


    ARRAY<std::string> object_names(3);
    object_names(1)= "New_Bones/Cranium_1";
    object_names(2)= "New_Bones/Pelvis_1";
    object_names(3)= "New_Bones/Right_Hand";

    int count = number/object_names.m + number%object_names.m;

    for (int i=1; i<=object_names.m; i++)
        Random_Scene_Generator(object_names(i), count, 11111, random_placement, rigid_body_collection, tests);

}
//#####################################################################
// Function Push_Out
//#####################################################################
void Push_Out()
{
    last_frame=240;
    //rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
    //rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
    //rigids_parameters.use_trapezoidal_rule_for_velocities=true;

    //int width=1;
    int height=parameter;
    //T offset=0;

    for(int y=0;y<height;y++)
    {
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,1);
        rigid_body.X()=TV(0,0.5+(y-1)*1.5,0);
        rigid_body.Set_Coefficient_Of_Restitution(0);
    }

    tests.Add_Ground((T).5,0,1,1);
}
//#####################################################################
// Function Angled_Bricks
//#####################################################################
void Angled_Bricks()
{
    last_frame=240;
    //rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
    //rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
    //rigids_parameters.use_trapezoidal_rule_for_velocities=true;

    int width=1;
    int height=parameter;
    T angle=0.1*pi;
    T angle_offset=(sin(angle)+cos(angle));

    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,0);
            rigid_body.X()=TV(x*2.0,(2*angle_offset)*y-angle_offset,0);
            rigid_body.Rotation()=ROTATION<TV>((y%2)?angle:-angle,VECTOR<T,3>(0,0,1));
            rigid_body.Set_Coefficient_Of_Restitution(0);
        }

    tests.Add_Ground(0,0,0,1);
}
//#####################################################################
};
}
#endif
