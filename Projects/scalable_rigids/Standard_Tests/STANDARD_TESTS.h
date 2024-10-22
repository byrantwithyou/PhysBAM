//#####################################################################
// Copyright 2010, Elliot English, Jon Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
// 1. Rings on pegs (concave, simple geometry)
// 2. Bones on pegs (concave, complex geometry) complex geometry will stress functions like particles_in_levelset
// 3. Spheres in a bin (analytics, large scale)
// 4. Drop three balls with different coefficients of restitution
// 5. Stack bricks forming a wall
//#####################################################################1
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"


namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    int parameter;

    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::viewer_dir;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;
    using BASE::user_last_frame;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),parameter(0)
    {
        parse_args.Add("-parameter",&parameter,"value","parameter used by multiple tests to change the parameters of the test");
        parse_args.Add_Not("-noanalytic",&solids_parameters.rigid_body_collision_parameters.use_analytic_collisions,"disable analytic collisions");
        parse_args.Add("-print_energy",&solid_body_collection.rigid_body_collection.print_energy,"print energy statistics");
        parse_args.Parse();

        tests.data_directory=data_directory;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
    }

    ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;

    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=true;

    switch(test_number){
        case 1: Ring_Test();break;
        case 2: Bone_Test();break;
        case 3: Sphere_Test();break;
        case 4: Bounce(0);break;
        case 5: Bricks();break;
        default: PHYSBAM_FATAL_ERROR(LOG::sprintf("Unrecognized test number %d",test_number));}

    tests.Add_Gravity();
}
//#####################################################################
// Function Ring_Test
//#####################################################################
void Ring_Test()
{
    if(!user_last_frame) last_frame=720;
    T mu=(T)0.6;
    T epsilon=(T)0.3;
    int poles=5;

    LOG::cout<<"Parameter 1: running 200 rings"<<std::endl;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,200,TV(0,30,0));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    Random_Scene_Generator("Rings_Test/ring_revolve",400,11111,random_placement,solid_body_collection.rigid_body_collection,tests);
    
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.coefficient_of_friction=mu;
        rigid_body.Set_Mass(10);
        rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,1,3));}

    for(int i=0;i<poles;i++)
        for(int j=0;j<poles;j++){
            // Poles
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,mu);
            rigid_body.name=LOG::sprintf("pole %d %d",i,j);
            rigid_body.is_static=true;
            rigid_body.Frame().t=TV((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);}

    tests.Add_Ground(mu);
}
//#####################################################################
// Function Bone_Test
//#####################################################################
void Bone_Test()
{
    if(!user_last_frame) last_frame=720;
    T mu=(T)0.6;
    T epsilon=(T)0.3;
    int poles=5;

    LOG::cout<<"Parameter 1: running 200 bones"<<std::endl;
    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,200,TV(0,30,0));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(3);
    Random_Scene_Generator("Bones/Cranium_Resized",400,11111,random_placement,solid_body_collection.rigid_body_collection,tests);
    
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.coefficient_of_friction=mu;
        rigid_body.Set_Mass(10);
        rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,1,3));}

    for(int i=0;i<poles;i++)
        for(int j=0;j<poles;j++){
            // Poles
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,mu);
            rigid_body.name=LOG::sprintf("pole %d %d",i,j);
            rigid_body.is_static=true;
            rigid_body.Frame().t=TV((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);}

    tests.Add_Ground(mu);
}
//#####################################################################
// Function Sphere_Test
//#####################################################################
void Sphere_Test()
{
    if(!user_last_frame) last_frame=720;

    TV lower_corner(0,0,0);
    TV upper_corner(10,10,10);
    TV center=(T)0.5*(upper_corner+lower_corner);
    TV size=(upper_corner-lower_corner);
    T coefficient_of_restitution=0;
    T coefficient_of_friction=0;

    STREAM_TYPE stream_type((float)0);

    RIGID_BODY<TV>& left_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    left_wall.Frame().t=VECTOR<T,3>(lower_corner.x,center.y,center.z);
    left_wall.Frame().r=ROTATION<VECTOR<T,3> >(-(T).5*(T)pi,VECTOR<T,3>(0,0,1));
    left_wall.name="left wall";

    RIGID_BODY<TV>& right_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.y,size.z));
    right_wall.Frame().t=VECTOR<T,3>(upper_corner.x,center.y,center.z);
    right_wall.Frame().r=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(0,0,1));
    right_wall.name="right wall";

    RIGID_BODY<TV>& bottom_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.z));
    bottom_wall.Frame().t=VECTOR<T,3>(center.x,lower_corner.y,center.z);
    bottom_wall.Frame().r=ROTATION<TV>(0,TV(0,0,1));
    bottom_wall.name="bottom wall";

    RIGID_BODY<TV>& front_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    front_wall.Frame().t=VECTOR<T,3>(center.x,center.y,lower_corner.z);
    front_wall.Frame().r=ROTATION<VECTOR<T,3> >((T).5*(T)pi,VECTOR<T,3>(1,0,0));
    front_wall.name="front wall";

    RIGID_BODY<TV>& back_wall=tests.Add_Ground(coefficient_of_friction,0,coefficient_of_restitution,max(size.x,size.y));
    back_wall.Frame().t=VECTOR<T,3>(center.x,center.y,upper_corner.z);
    back_wall.Frame().r=ROTATION<VECTOR<T,3> >((T)-.5*(T)pi,VECTOR<T,3>(1,0,0));
    back_wall.name="back wall";

    //char* object_names[]={"subdivided_box","sphere","New_Bones/Cranium_1","New_Bones/Pelvis_1","New_Bones/Right_Hand"};//,"dumbbell","bowl","spoon"};
    //T object_scales[]={1,1,30,30,30};

    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(min(size(1),size(3))/2.5,30,TV(center(1),1,center(3)));
    random_placement.Set_Max_Orientation_Angle((T)0.2);
    random_placement.Set_Speed_Range(0,1);
    random_placement.Set_Angular_Speed_Range(0,1);
    random_placement.Set_Fixed_Scale(0.5);
    Random_Scene_Generator("sphere",500,11111,random_placement,solid_body_collection.rigid_body_collection,tests);
}
//#####################################################################
// Function Bounce
//#####################################################################
void Bounce(const T angle)
{
    if(!user_last_frame) last_frame=240;

    T x_pos[]={-3,0,3},cor[]={(T)1.0,(T).5,0};
    for(int i=0;i<3;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",1,(T).5);
        rigid_body.Frame().t=TV(x_pos[i],5,0);rigid_body.Set_Coefficient_Of_Restitution(cor[i]);
        rigid_body.name=LOG::sprintf("sphere (cor %g)",cor[i]);}

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",1,(T).5);
    rigid_body.Frame().t=TV(x_pos[2],8,0);rigid_body.Set_Coefficient_Of_Restitution(0);
    rigid_body.name=LOG::sprintf("sphere");

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,0,1,1);
    ground.Frame().r=ROTATION<TV>(angle,TV(0,0,1));
}
//#####################################################################
// Function Bricks
//#####################################################################
void Bricks()
{
    solids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel=true;

    if(!user_last_frame) last_frame=240;

    int width=1;
    int height=1;
    T offset=0;

    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,1);
            rigid_body.Frame().t=TV(2.0*x+offset*(y%2),2.0*y-1.0,0);
            rigid_body.Set_Coefficient_Of_Restitution(0);
        }

    tests.Add_Ground((T).5,0,1,1);
}
//#####################################################################
};
}
#endif
