//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Kinematic box stack
//   2. pile of boxes and circles
//   3. cluster examples
//   4. Pyramid of boxes
//   5. Stacked boxes
//   6. Partition Test
//   7. Contact Test 1
//   8. Contact Test 2
//   9. Simple Collision Test
//   10. Collision Test
//   11. Pushout Test

// TODOs
// 3. restarts
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
#include <Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    int width, height, num_bodies;

    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::frame_rate;
    using BASE::data_directory;using BASE::last_frame;using BASE::viewer_dir;using BASE::stream_type;
    using BASE::user_last_frame;
    
    SOLIDS_STANDARD_TESTS<TV> tests;

    int kinematic_body_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),width(2),height(4),num_bodies(6),tests(stream_type_input,data_directory,solid_body_collection)
    {
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        //solids_parameters.rigid_body_collision_parameters.use_shock_propagation=false;
        solids_parameters.cfl=1;
        parse_args.Add("-width",&width,"value","number of stacks");
        parse_args.Add("-height",&height,"value","height of each stack");
        parse_args.Add("-num_bodies",&num_bodies,"value","number of total bodies");
        parse_args.Parse();

        tests.data_directory=data_directory;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
    }

    ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) override
{
    if(time>(T)2 && test_number==3){
        static bool deleted=false;
        if(!deleted){
            //solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Set_Binding_Active(4,false);
            //solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Delete_Binding(4,solids_parameters.collision_body_list);
            deleted=true;}}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) override
{
    if(test_number==11){
        RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
        solids_parameters.rigid_body_collision_parameters.perform_collisions=false;
        solids_parameters.rigid_body_collision_parameters.collision_iterations=0;
        solids_parameters.rigid_body_collision_parameters.contact_iterations=0;
        collisions.Set_Shock_Propagation_Iterations(0);}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    switch(test_number){
      case 1: Kinematic();break;
      case 2: Test_Example();break;
      case 3: Cluster();break;
      case 4: Pyramid_Of_Boxes(); break;
      case 5: Stacked_Boxes(); break;
      case 6: Partition_Test(); break;
      case 7: Contact_Test_1(); break;
      case 8: Contact_Test_2(); break;
      case 9: Simple_Collision_Test(); break;
      case 10: Collision_Test(); break;
      case 11: Pushout_Test(); break;
      default: PHYSBAM_FATAL_ERROR(LOG::sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Kinematic
//#####################################################################
void Kinematic()
{
    RIGID_BODY<TV>* rigid_body=0;
    T baseboxsize=2;
    T boxsize1=(T)1.0;
    T boxsize2=(T)0.5;
    T smallboxmass=1;
    T stack_epsilon=(T)0.3;
    T stack_mu=(T)0.5;
    const char* boxfile="square";

    rigid_body=&tests.Add_Rigid_Body(boxfile,boxsize1,stack_mu);
    rigid_body->Frame().t=TV(0,2*baseboxsize+boxsize1);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->name="stack box 1a";

    rigid_body=&tests.Add_Rigid_Body(boxfile,boxsize2,stack_mu);
    rigid_body->Frame().t=TV(0,2*baseboxsize+2*boxsize1+boxsize2);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->name="stack box 1b";

    //boxfile="square_refined";
    rigid_body=&tests.Add_Rigid_Body(boxfile,baseboxsize,stack_mu);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->name="base box";
    solid_body_collection.rigid_body_collection.rigid_body_particles.kinematic(rigid_body->particle_index)=true;
    kinematic_body_id=rigid_body->particle_index;

    TV t0(0,baseboxsize),t1(0,baseboxsize+(T)17.5);
    ROTATION<TV> r0,r1(ROTATION<TV>::From_Angle((T)0.8));
    curve.Add_Control_Point((T).5,FRAME<TV>(t0,r0));
    curve.Add_Control_Point((T)4,FRAME<TV>(t1,r0));
    curve.Add_Control_Point((T)6,FRAME<TV>(t1,r0));
    curve.Add_Control_Point((T)8,FRAME<TV>(t1,r1));

    rigid_body->Frame()=curve.Value(0);
    rigid_body->Twist()=curve.Derivative(0);

    tests.Add_Ground((T).3,0,1);

    if(!user_last_frame) last_frame=(int)(15*frame_rate);

    // add forces
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Test_Example
//#####################################################################
void Test_Example()
{
    RECTANGULAR_RANDOM_PLACEMENT<TV> random_placement(RANGE<TV>(TV(-20,1),TV(20,100)));
    Random_Scene_Generator("circle", 100, 1234, random_placement,solid_body_collection.rigid_body_collection,tests);
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++)
        solid_body_collection.rigid_body_collection.Rigid_Body(i).Set_Coefficient_Of_Restitution((T)0.5);

    RIGID_BODY<TV>* rigid_body0=&tests.Add_Rigid_Body("circle",(T)1,(T).1);
    rigid_body0->Frame().r=ROTATION<TV>::From_Angle((T)pi/5);
    rigid_body0->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body0->name="circle";

    RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("circle",(T)2,(T).1);
    rigid_body1->Frame().t=TV(0,120);
    rigid_body1->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body1->name="circle";

    tests.Add_Ground(1,-20);
    //if(!user_last_frame) last_frame=(int)(10*frame_rate);
    if(!user_last_frame) last_frame=200;

    // add forces
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Pyramid_Of_Boxes
//#####################################################################
void Pyramid_Of_Boxes()
{
    int height=21;

    T current_x;
    T first_x=0;
    for (int i=1; i < height; i++) {
        current_x=first_x;
        for (int j= 1; j <= i ; j++) {
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square_refined",(T)2, (T).1);
            rigid_body->Frame().t=TV((T)current_x, (height+1-i)*5+80);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
            current_x += 5;}
        first_x -= 2.5;}

    for(int i=0;i<2;i++){
        RIGID_BODY<TV>* rigid_body0=&tests.Add_Rigid_Body("square_refined",(T)10,(T).1);
        rigid_body0->Frame().t=TV(first_x-10,i*20);
        rigid_body0->is_static=true;
        rigid_body0->name="left_square";

        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square_refined",(T)10,(T).1);
        rigid_body1->Frame().t=TV(current_x+7.5,i*20);
        rigid_body1->is_static=true;
        rigid_body1->name="right_square";}

    tests.Add_Ground(1, -10);
    if(!user_last_frame) last_frame=400;
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Stacked_Boxes
//#####################################################################
void Stacked_Boxes() {
    int height=50;

    for (int i=1; i <= height; i++) {
        for (int j= 1; j < 20 ; j++) {
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square_refined",(T)2, (T).1);
            rigid_body->Frame().t=TV(5*j - 50, 5*i+10);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);}}

    for (int i=0; i < height/4; i++) {
        RIGID_BODY<TV>* rigid_body0=&tests.Add_Rigid_Body("square_refined",(T)10,(T).1);
        rigid_body0->Frame().t=TV(-58,i*20);
        rigid_body0->is_static=true;
        rigid_body0->name="left_square";
    
        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square_refined",(T)10,(T).1);
        rigid_body1->Frame().t=TV(58,i*20);
        rigid_body1->is_static=true;
        rigid_body1->name="right_square";}

    tests.Add_Ground(1, -10);
    if(!user_last_frame) last_frame=250;
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Partition_Test
//#####################################################################
void Partition_Test() {
    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square", (T)1, (T).1);
    rigid_body->Frame().t=TV(0,0);
    rigid_body->Twist().linear=TV(1,0);

    RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square", (T)1, (T).1);
    rigid_body1->Frame().t=TV(10,3);
    rigid_body1->Twist().linear=TV(-1,0);

    if(!user_last_frame) last_frame=200;
}
//#####################################################################
// Function Cluster
//#####################################################################
void Cluster()
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;

    RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>* cluster_fracture_callbacks=new RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>(solid_body_collection.rigid_body_collection,rigid_bindings);
    rigid_bindings.callbacks=cluster_fracture_callbacks;
    cluster_fracture_callbacks->allowed_strain=(T).1;
    

    const char* boxfile="square_refined";

    // add ground first
    tests.Add_Ground((T).1,(T)-4);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_1->Frame().t=TV(0,2);
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->name="box1";

    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_2->Frame().t=TV(2,3);
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->name="box2";

    RIGID_BODY<TV>* rigid_body_3=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_3->Frame().t=TV(4,3);
    rigid_body_3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_3->name="box3";

    RIGID_BODY<TV>* rigid_body_4=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_4->Frame().t=TV(4,5);
    rigid_body_4->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_4->name="box4";

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    children.Append(rigid_body_1->particle_index);
    children.Append(rigid_body_2->particle_index);
    children.Append(rigid_body_3->particle_index);
    children.Append(rigid_body_4->particle_index);
    int cluster_particle=rigid_bindings.Add_Binding(children);
    RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster";

    // make duplicate object non-clustered
    RIGID_BODY<TV>* rigid_body_cluster_test=new RIGID_BODY<TV>(solid_body_collection.rigid_body_collection,true);
    rigid_body_cluster_test->name="clustertest";
    rigid_body_cluster_test->Frame()=FRAME<TV>(TV((T)15,0))*rigid_body_cluster->Frame();
    SEGMENTED_CURVE_2D<T>* segmented_curve=SEGMENTED_CURVE_2D<T>::Create();
    segmented_curve->mesh.elements=rigid_body_cluster->simplicial_object->mesh.elements;
    segmented_curve->particles.Add_Elements(rigid_body_cluster->simplicial_object->particles.Size());
    segmented_curve->mesh.number_nodes=segmented_curve->particles.Size();
    for(int i=0;i<segmented_curve->particles.Size();i++){
        segmented_curve->particles.X(i)=rigid_body_cluster->simplicial_object->particles.X(i);
        LOG::cout<<"particles.X("<<i<<")"<<" is "<<segmented_curve->particles.X(i)<<std::endl;}
    rigid_body_cluster_test->Add_Structure(*segmented_curve);
    solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body_cluster_test);
    for(int i=0;i<rigid_body_cluster->simplicial_object->mesh.elements.m;i++){
        LOG::cout<<"mesh "<<i<<" -> "<<rigid_body_cluster->simplicial_object->mesh.elements(i)<<std::endl;
    }
    MASS_PROPERTIES<TV> mp(*rigid_body_cluster->simplicial_object,false);mp.Set_Mass(rigid_body_cluster->Mass());
    rigid_body_cluster_test->Mass()=mp.Mass();
    rigid_body_cluster_test->Inertia_Tensor()=DIAGONAL_MATRIX<T,1>(mp.Inertia_Tensor().x00);
    LOG::cout<<"tensor sucks "<<rigid_body_cluster_test->Mass()<<std::endl;
    LOG::cout<<"old mass "<<rigid_body_cluster->Mass()<<std::endl;
    PHYSBAM_ASSERT(rigid_body_cluster_test->Rigid_Mass().Valid());

    // add gravity
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_3->particle_index);
    referenced_rigid_particles->Append(rigid_body_4->particle_index);
    referenced_rigid_particles->Append(rigid_body_cluster_test->particle_index);
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));

    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++){
        LOG::cout<<"Rigid body index "<<i<<solid_body_collection.rigid_body_collection.Rigid_Body(i).name<<std::endl;
        LOG::cout<<"    FRAME  "<<solid_body_collection.rigid_body_collection.Rigid_Body(i).Frame()<<std::endl;
        LOG::cout<<"    MASS  "<<solid_body_collection.rigid_body_collection.rigid_body_particles.mass(i)<<std::endl;
        LOG::cout<<"    INERTIA TENSOR  "<<solid_body_collection.rigid_body_collection.rigid_body_particles.inertia_tensor(i)<<std::endl;
        LOG::cout<<"    TWIST  "<<solid_body_collection.rigid_body_collection.rigid_body_particles.twist(i)<<std::endl;
        LOG::cout<<"    cluster parent "<<((i<=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.binding_index.m && solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.binding_index(i).Size()>0)?solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.binding_index(i)(1).x:0)<<std::endl;
    }                                                 
}
//#####################################################################
// Function Contact_Test_1
//#####################################################################
void Contact_Test_1() {
    for (int i=0;i<width;i++) {
        for (int j=0;j<height;j++) {
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
            rigid_body->Frame().t=TV(i*4,j*2);}}

    tests.Add_Ground(1, -1);
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Contact_Test_2
//#####################################################################
void Contact_Test_2() {
    if (height>width)
        height=width;

    for (int i=0;i<height;i++) {
        for (int j=0;j<width;j++) {
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
            rigid_body->Frame().t=TV((T)j*2+i*0.5,i*2);}
        width--;}

    tests.Add_Ground(1, -1);
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection));
}
//#####################################################################
// Function Simple_Collision_Test
//#####################################################################
void Simple_Collision_Test() {
    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square", (T)1, (T).1);
    rigid_body->Frame().t=TV(0,0);
    rigid_body->Twist().linear=TV(2,0);

    RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square", (T)1, (T).1);
    rigid_body1->Frame().t=TV(10,3);
    rigid_body1->Twist().linear=TV(-1,0);

    RIGID_BODY<TV>* rigid_body2=&tests.Add_Rigid_Body("square", (T)1, (T).1);
    rigid_body2->Frame().t=TV(10,0);
    rigid_body2->Twist().linear=TV(-1,0);

    if(!user_last_frame) last_frame=200;
}
//#####################################################################
// Function Collision_Test
//#####################################################################
void Collision_Test() {
    for (int i=0;i<height;i++) {
        RIGID_BODY<TV>* left_body=&tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
        RIGID_BODY<TV>* right_body=&tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
        left_body->Frame().t=TV(0,i*3);
        right_body->Frame().t=TV((width+1)*3,i*3);
        left_body->Set_Coefficient_Of_Restitution((T)1.0);
        right_body->Set_Coefficient_Of_Restitution((T)1.0);
        left_body->is_static=true;
        right_body->is_static=true;

        for (int j=1;j<=width;j++) {
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
            rigid_body->Frame().t=TV(j*3,i*3);
            rigid_body->Twist().linear=TV(4*(.5-j%2),0);
            rigid_body->Set_Coefficient_Of_Restitution((T)1.0);}}

    if(!user_last_frame) last_frame=300;
}
//#####################################################################
// Function Pushout_Test
//#####################################################################
void Pushout_Test() {
    for (int i=0;i<num_bodies;i++) {
        RIGID_BODY<TV>* rigid_body= &tests.Add_Rigid_Body("square_refined",(T)1,(T).1);
        rigid_body->Frame().t=TV((T)1.75*i,(T)1.75*i);}
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    if(test_number==1 && id==kinematic_body_id) twist=curve.Derivative(time);
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    if(test_number==1 && id==kinematic_body_id) frame=curve.Value(time);
}
//#####################################################################
};
}
#endif
