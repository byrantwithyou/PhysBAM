//#####################################################################
// Copyright 2007-2008, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. drop three balls with different coefficients of restitution
//   2. drop three nonconvex objects with different coefficients of restitution
//   3. five billiard balls in a row impacted with a sixth ball
//   4. drop spinning ball, which transfers spin into velocity
//   5. two balls slide with different spins
//   6. <none>
//   7. big and little block balance on a seesaw (no contact)
//   8. big and little block balance on a seesaw
//   9. stack of blocks on a plank upset by falling block (no contact)
//  10. stack of blocks on a plank upset by falling block
//  11. <none>
//  12. rings fall on a set of pegs
//  13. bones fall on a set of pegs
//  14. suspended cubes
//  15. <none>
//  16. ether test
//  17. wind drag test
//  18. friction and force propogation test
//  19. suspended cube
//  20. <none>
//  21. sphere rolling down incline - analytic test
//  22. spinning top
//  23. spinning cylinder
//  24. drop three balls with different coefficients of restitution on inclined plane
//  25. <none>
//  26. cluster cluster testing
//  27. cluster cluster testing with a kinematic body
//  28. push out with cylinder
//  29. Boxes falling with removed rigid bodies
//  30. Collision with kinematic rigid body
//  31. Sanity Tests
//  32. Drop Cubes
//  33. Collision/Contact pairs sticking test
//  34. Kinematically Deforming Sphere
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Forces_And_Torques/RIGID_ETHER_DRAG.h>
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/KINEMATIC_COLLISION_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_CLUSTER_CONSTITUENT_ID.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#include "../RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    T small_block_mass;
    int parameter;
    INTERPOLATION_CURVE<T,TV> curve_left,curve_right,curve_plank,curve;
    int curve_left_id,curve_right_id,curve_plank_id,torus_id;
    ARRAY<int> rigid_body_particles_with_ether_drag; // test 16

    KINEMATIC_COLLISION_BODY<TV>* deforming_sphere;

    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;
    bool print_matrix;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;
    using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::stream_type;using BASE::parse_args;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),small_block_mass(1),
        parameter(0),collision_manager(0),print_matrix(false)
    {
    }

    ~STANDARD_TESTS()
    {
        if(test_number==21){
            T angle=parameter*(T)pi/40;
            T y=solid_body_collection.rigid_body_collection.rigid_body_particles.frame(0).t.y,yanalytic=-(T)5./14*(T)9.8*sqr((T)last_frame/24)*sqr(sin(angle));
            T ke=solid_body_collection.rigid_body_collection.Rigid_Body(0).Kinetic_Energy(),pe=y*solid_body_collection.rigid_body_collection.rigid_body_particles.mass(0)*(T)9.8;
            T omega=solid_body_collection.rigid_body_collection.rigid_body_particles.twist(0).angular.z,speed=solid_body_collection.rigid_body_collection.rigid_body_particles.twist(0).linear.Magnitude();
            LOG::cout<<"ERROR in y for angle "<<angle<<"  ("<<(angle*180/(T)pi)<<" degrees):  "<<(1-y/yanalytic)<<std::endl;
            LOG::cout<<"SLIPPAGE  "<<omega/speed<<std::endl;LOG::cout<<"PE = "<<pe<<"    KE = "<<ke<<std::endl;}
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
    {
        dynamic_cast<NEWMARK_EVOLUTION<TV>&>(*solids_evolution).print_matrix=print_matrix;
    }
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}

    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

    void Register_Options() PHYSBAM_OVERRIDE
    {
        BASE::Register_Options();
        parse_args->Add("-small_block_mass",&small_block_mass,"value","mass for small blocks in plank test");
        parse_args->Add("-parameter",&parameter,"value","parameter used by multiple tests to change the parameters of the test");
        parse_args->Add_Not("-noanalytic",&solids_parameters.rigid_body_collision_parameters.use_analytic_collisions,"disable analytic collisions");
        parse_args->Add("-print_energy",&solid_body_collection.rigid_body_collection.print_energy,"print energy statistics");
        parse_args->Add("-print_matrix",&print_matrix,"Print Krylov matrix");
    }
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        tests.data_directory=data_directory;
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);

        if(small_block_mass!=1) output_directory+=STRING_UTILITIES::string_sprintf("_m%g",small_block_mass);
        if(parameter) output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);
    }
    
    void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}

//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    collisions.Set_Push_Out_Level_Iterations(1);
    if(test_number==26||test_number==27) collisions.collision_manager=collision_manager;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==14){
        if(id==curve_left_id) frame.t=curve_left.Value(time);
        else if(id==curve_right_id) frame.t=curve_right.Value(time);
        else if(id==curve_plank_id){frame.t=curve_plank.Value(time);frame.r=ROTATION<TV>((T)pi/2,TV(0,1,0));}}
    if(test_number==30) frame.t=curve.Value(time);
    if(test_number==34 && id==deforming_sphere->particle_index)
        Build_Deforming_Sphere(deforming_sphere,frame,time,true,false);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==14){
        if(id==0) twist.linear=curve_left.Derivative(time);
        else if(id==1) twist.linear=curve_right.Derivative(time);
        else if(id==2) twist.linear=curve_plank.Derivative(time);}
    else if(test_number==30) twist.linear=curve.Derivative(time);
    else if(test_number==34)
        Build_Deforming_Sphere(deforming_sphere,deforming_sphere->Frame(),time,false,true);
    else return false;
    return true;
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==19){
        bool &is_static=solid_body_collection.rigid_body_collection.Rigid_Body(4).is_static,new_is_static=time<=(T).5;
        if(is_static!=new_is_static){is_static=new_is_static;solid_body_collection.rigid_body_collection.Update_Simulated_Particles();}} // Changing is_static requires an Update_Simulated_Particles
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==31){
        if(parameter==4 && time>=(T).4999 && time<=(T).5001){
            solid_body_collection.rigid_body_collection.Reactivate_Body(1);
            solid_body_collection.rigid_body_collection.Update_Simulated_Particles();}
        if(parameter==3 && time>=(T).4999 && time<=(T).5001){
            solid_body_collection.rigid_body_collection.rigid_body_particles.Remove_Body(0);
            solid_body_collection.rigid_body_collection.Update_Simulated_Particles();}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;

    switch(test_number){
      case 1: Bounce(0);break;
      case 2: Nonconvex_Bounce();break;
      case 3: Billiard_Balls();break;
      case 4: Bounce_Friction();break;
      case 5: Rolling_Friction();break;
      case 7: solids_parameters.rigid_body_collision_parameters.contact_iterations=0; //fall through
      case 8: Seesaw();break;
      case 9: solids_parameters.rigid_body_collision_parameters.contact_iterations=0; //fall through
      case 10: Plank();break;
      case 12: Ring_Test();break;
      case 13: Bone_Test();break;
      case 14: Suspended_Cubes();break;
      case 16: Ether_Test();break;
      case 17: Wind_Test();break;
      case 18: Force_Propogation();break;
      case 19: Suspended_Cube();break;
      case 21: Sphere_Incline();break;
      case 22: Spinning_Top();break;
      case 23: Spinning_Cylinder();break;
      case 24: Bounce((T)pi/6);break;
      case 26: Cluster_Cluster_Interaction();break;
      case 27: Cluster_Cluster_Interaction_Kinematic();break;
      case 28: Push_Out();break;
      case 29: Removed_Bodies();break;
      case 30: Kinematic_Collision();break;
      case 31: Sphere_Sanity_Tests();break;
      case 32: Drop_Cubes();break;
      case 33: Collision_Contact_Pairs_Test();break;
      case 34: Deforming_Sphere();break;
      case 35: Cluster_Fracture();break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
    
    // add forces
    if(test_number==16){
        rigid_body_particles_with_ether_drag.Append(1);
        solid_body_collection.Add_Force(new RIGID_ETHER_DRAG<TV>(solid_body_collection.rigid_body_collection,&rigid_body_particles_with_ether_drag,(T).5,(T).5));}
    else if(test_number==17){
        for(int i=2;i<=3;i++){
            PHYSBAM_NOT_IMPLEMENTED("RIGID_BODY_BASIC_FORCES is obsolete");
//             RIGID_BODY_BASIC_FORCES<TV>* basic_forces=0; // do rigid body i
//             basic_forces->Use_Wind_Drag();
//             basic_forces->Set_Wind_Density((T).1);
//             basic_forces->Use_Constant_Wind(0,TV(1,0,0));
        }}
    if(test_number!=23 && test_number!=26 && test_number!=27 && test_number!=28) solid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,true));
}
//#####################################################################
// Function Bounce
//#####################################################################
void Bounce(const T angle)
{
    last_frame=240;
    //T x_pos[]={-3,0,3},cor[]={(T)1.0,(T).5,0};
    for(int i=0;i<100;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",1,(T).5);
        rigid_body.Frame().t=TV(0,5*i,0);rigid_body.Set_Coefficient_Of_Restitution(.1);
        rigid_body.name=STRING_UTILITIES::string_sprintf("sphere %i",i);}

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,0,1,1);
    ground.Frame().r=ROTATION<TV>(angle,TV(0,0,1));
}
//#####################################################################
// Function Nonconvex_Bounce
//#####################################################################
void Nonconvex_Bounce()
{
    last_frame=240;
    T x_pos[]={-3,0,3},cor[]={(T)1.0,(T).6,0},z_angle[]={(T).05,(T).02,0};
    for(int i=0;i<3;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("New_Bones/Pelvis_1",20,(T).5);
        rigid_body.Set_Mass(1);
        rigid_body.Frame().t=TV(x_pos[i],5,0);
        rigid_body.Frame().r=ROTATION<TV>(z_angle[i],TV(0,0,1))*ROTATION<TV>((T)-0.5,TV(1,0,0))*ROTATION<TV>((T)pi,TV(0,1,0));
        rigid_body.Set_Coefficient_Of_Restitution(cor[i]);
        rigid_body.name=STRING_UTILITIES::string_sprintf("bone (cor %g)",cor[i]);}

    tests.Add_Ground((T).5,0,1);
}
//#####################################################################
// Function Billiard_Balls
//#####################################################################
void Billiard_Balls()
{
    last_frame=240;
    T mu=(T)0.3,cor=(T)1.0;
    int number=5;
    for (int k=0;k<=number;k++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",(T)1,mu);
        rigid_body.Frame().t=TV((T)(k?2*(k-(number+1)/2):-30),1,0);
        if(!k) rigid_body.Twist().linear=TV(20,0,0);
        rigid_body.name=STRING_UTILITIES::string_sprintf("sphere%d",k);
        rigid_body.Set_Coefficient_Of_Restitution(cor);}

    tests.Add_Ground((T)1,0,0).Set_Coefficient_Of_Rolling_Friction(1);
}
//#####################################################################
// Function Bounce_Friction
//#####################################################################
void Bounce_Friction()
{
    last_frame=360;
    RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",(T)1,(T).8);
    sphere.Frame().t=TV(-4,10,0);
    sphere.Twist().angular=TV(-10,0,0);
    sphere.Update_Angular_Momentum();
    sphere.Set_Coefficient_Of_Restitution((T)0.5);
    sphere.Set_Coefficient_Of_Rolling_Friction((T)0.01);
    sphere.name="spinning sphere";

    tests.Add_Ground((T)1,0,1).Set_Coefficient_Of_Rolling_Friction(1);
}
//#####################################################################
// Function Rolling_Friction
//#####################################################################
void Rolling_Friction()
{
    last_frame=240;
    RIGID_BODY<TV>* rigid_body=0;

    if(parameter==0){
        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).1);
        rigid_body->Frame().t=TV(-2,1,45);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Coefficient_Of_Rolling_Friction((T)0.01);
        rigid_body->name="sphere (mu=0.05)";

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).1);
        rigid_body->Frame().t=TV(2,1,45);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Twist().angular=TV(0,0,-10);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Coefficient_Of_Rolling_Friction((T)0.2);
        rigid_body->name="sphere (mu=0.2)";}
    else{
        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).05);
        rigid_body->Frame().t=TV(-8,1,25);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->name="sphere (mu=0.05)";

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).02);
        rigid_body->Frame().t=TV(8,1,25);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->name="sphere (mu=0.02)";

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).5);
        rigid_body->Frame().t=TV(-2,10,0);
        rigid_body->Twist().angular=TV(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->name="spinning sphere";

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).05);
        rigid_body->Frame().t=TV(2,10,0);
        rigid_body->Twist().angular=TV(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->name="low friction spinning sphere";}

    tests.Add_Ground(1,0,1).Set_Coefficient_Of_Rolling_Friction(1);
}
//#####################################################################
// Function Seesaw
//#####################################################################
void Seesaw()
{
    last_frame=480;
    RIGID_BODY<TV>* rigid_body=0;
    T mu=(T)0.2;
    bool fromrest=true;

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T)0.9,(T).15);
    rigid_body->Set_Mass(20);
    rigid_body->Frame().t=TV((T)3.5,fromrest?(T)2.4:20,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->name="box1";

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T)0.5,mu);
    rigid_body->Set_Mass(5);
    rigid_body->Frame().t=TV(-4,fromrest?(T)2:(T)3,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->name="box2";

    rigid_body=&tests.Add_Rigid_Body("plank",1,1);
    rigid_body->Frame().t=TV(0,fromrest?(T)1.25:(T)1.5,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Mass(10);
    rigid_body->name="plank";

    rigid_body=&tests.Add_Rigid_Body("Rings_Test/cylinder_revolve",1,mu);
    rigid_body->Frame().t=TV(0,(T)0.5,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
    rigid_body->name="cylinder";
    rigid_body->is_static=true;
    rigid_body->Set_Coefficient_Of_Restitution(1);

    tests.Add_Ground(1,0,1);
}
//#####################################################################
// Function Plank
//#####################################################################
void Plank()
{
    last_frame=480;
    RIGID_BODY<TV>* rigid_body=0;
    T baseboxsize=2;
    T dropboxsize=(T)1.2;
    T plankscale=1;
    T dropheight=70;
    T smallboxsize=(T)0.5;
    T smallboxmass=small_block_mass;
    T offset=(T)0.05;
    const char* boxfile=parameter?"box":"subdivided_box";
    const char* plankfile=parameter?"unsubdivided_plank":"plank";
    T stack_epsilon=(T)0.3;
    T stack_mu=(T)0.5;
    const char* name[]={"1a","1b","2","3","4","5","6"};
    T x_pos[]={-(T)0.65,(T)0.65,offset,0,-offset,0,0};
    T y_pos[]={4,5,7,9,14,17,19};
    T z_pos[]={0,0,0,offset,-offset,0,0};

    for(int i=0;i<7;i++){
        rigid_body=&tests.Add_Rigid_Body(boxfile,smallboxsize,stack_mu);
        rigid_body->Frame().t=TV(x_pos[i],baseboxsize+y_pos[i],z_pos[i]);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->name=STRING_UTILITIES::string_sprintf("stack box %s",name[i]);}
    rigid_body->Frame().r=ROTATION<TV>((T)pi/4,TV(0,1,0));

    rigid_body=&tests.Add_Rigid_Body(boxfile,1,1);
    rigid_body->Frame().t=TV(baseboxsize+plankscale*5-dropboxsize,2*baseboxsize+(T)0.5*plankscale+dropboxsize+dropheight,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Mass(100);
    rigid_body->name="drop box";

    rigid_body=&tests.Add_Rigid_Body(plankfile,plankscale,1);
    rigid_body->Frame().t=TV(baseboxsize,2*baseboxsize+(T)0.25,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->Set_Mass(10);
    rigid_body->name="plank";

    rigid_body=&tests.Add_Rigid_Body(boxfile,baseboxsize,1);
    rigid_body->Frame().t=TV(0,baseboxsize,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->coefficient_of_friction=(T)0.3;
    rigid_body->name="base box";
    rigid_body->is_static=true;

    tests.Add_Ground((T).3,0,1);
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

    if(parameter==0){
        LOG::cout<<"Parameter 0: running 500 rings"<<std::endl;
        CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,300,TV(0,30,0));
        random_placement.Set_Max_Orientation_Angle((T)0.2);
        random_placement.Set_Speed_Range(0,1);
        random_placement.Set_Angular_Speed_Range(0,1);
        Random_Scene_Generator("Rings_Test/ring_revolve",500,12345,random_placement,solid_body_collection.rigid_body_collection,tests);}
    else{
        LOG::cout<<"Parameter 1: running 1000 rings"<<std::endl;
        CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,600,TV(0,30,0));
        random_placement.Set_Max_Orientation_Angle((T)0.2);
        random_placement.Set_Speed_Range(0,1);
        random_placement.Set_Angular_Speed_Range(0,1);
        Random_Scene_Generator("Rings_Test/ring_revolve",1000,11111,random_placement,solid_body_collection.rigid_body_collection,tests);}

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
            rigid_body.name=STRING_UTILITIES::string_sprintf("pole %d %d",i,j);
            rigid_body.is_static=true;
            rigid_body.Frame().t=TV((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);}

    tests.Add_Ground(mu);
}
//#####################################################################
// Function Bone_Test
//#####################################################################
void Bone_Test()
{
    last_frame=240;
    RIGID_BODY<TV>* rigid_body=0;

    ARRAY<std::string>filenames;
    for(int i=0;i<300;i++) filenames.Append("New_Bones/Cranium_1");
    for(int i=0;i<150;i++) filenames.Append("New_Bones/Pelvis_1");
    for(int i=0;i<50;i++) filenames.Append("New_Bones/Left_Femur_2");

    T mu=(T)0.4;
    T epsilon=(T)0.3;
    T rolling_friction=(T)0.1;

    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(1,40,TV(0,3,0));
    random_placement.Set_Fixed_Scale(3);
    Random_Scene_Generator(filenames,12345,random_placement,solid_body_collection.rigid_body_collection,tests);
    T mass_scale=10;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Mass(rigid_body.Mass()*mass_scale);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.coefficient_of_friction=mu;
        rigid_body.Set_Coefficient_Of_Rolling_Friction(rolling_friction);

        if(rigid_body.name.find("Cranium"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(2,2,2));
        else if(rigid_body.name.find("Left_Femur"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,4,1));
        else if(rigid_body.name.find("Pelvis"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,3,2));}

    rigid_body=&tests.Add_Rigid_Body("funnel_thicker_revolve",(T)0.4,1);
    rigid_body->Frame().t=TV(0,2,0);
    rigid_body->name="funnel";
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,3,3));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    T yy=(T)0.2;
    T ss=(T)0.25;
    T xx=ss*(T)4.75;

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->Frame().t=TV(xx,yy,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->name="plank 1";
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->Frame().t=TV(0,yy,xx);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->name="plank 2";
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->Frame().t=TV(-xx,yy,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->name="plank 3";
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->Frame().t=TV(0,yy,-xx);
    rigid_body->Frame().r=ROTATION<TV>(3*(T)pi/2,TV(0,1,0))*rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->name="plank 4";
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Ground(mu);
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,1,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(rolling_friction);
}
//#####################################################################
// Function Suspended_Cubes
//#####################################################################
void Suspended_Cubes()
{
    last_frame=240;
    T mu=(T)0.4,height=5,cor=(T).5;
    RIGID_BODY<TV>& ram1=tests.Add_Rigid_Body("subdivided_box",1,mu);
    ram1.Set_Coefficient_Of_Restitution(cor);
    ram1.Is_Kinematic()=true;
    curve_left_id=ram1.particle_index;
    curve_left.Add_Control_Point((T)0,TV(-7,height,0));
    curve_left.Add_Control_Point((T)1,TV(-6+(T)1e-4,height,(T)0));
    curve_left.Add_Control_Point((T)5,TV(-6+(T)1e-4,height,(T)0));
    curve_left.Add_Control_Point((T)6,TV(-7,height,(T)0));
    ram1.Frame().t=curve_left.Value(0);
    ram1.Twist().linear=curve_left.Derivative(0);

    RIGID_BODY<TV>& ram2=tests.Add_Rigid_Body("subdivided_box",1,mu);
    ram2.Set_Coefficient_Of_Restitution(cor);
    ram2.Is_Kinematic()=true;
    curve_right_id=ram2.particle_index;
    curve_right.Add_Control_Point((T)0,TV(7,height,0));
    curve_right.Add_Control_Point((T)1,TV(6-(T)1e-4,height,(T)0));
    curve_right.Add_Control_Point((T)5,TV(6-(T)1e-4,height,(T)0));
    curve_right.Add_Control_Point((T)6,TV(7,height,(T)0));
    ram2.Frame().t=curve_right.Value(0);
    ram2.Twist().linear=curve_right.Derivative(0);

    RIGID_BODY<TV>& plank=tests.Add_Rigid_Body("plank",1,mu);
    plank.Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    plank.Set_Coefficient_Of_Restitution(cor);
    plank.Is_Kinematic()=true;
    curve_plank_id=plank.particle_index;
    curve_plank.Add_Control_Point((T)0,TV(0,height-(T)1.25,0));
    curve_plank.Add_Control_Point((T)2,TV(0,height-(T)1.25,(T)0));
    curve_plank.Add_Control_Point((T)3,TV(0,1,(T)0));
    curve_plank.Add_Control_Point((T)4,TV(0,1,(T)4));
    curve_plank.Add_Control_Point((T)6,TV(0,1,(T)4));
    plank.Frame().t=curve_plank.Value(0);
    plank.Twist().linear=curve_plank.Derivative(0);

    for(int i=0;i<5;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,mu);
        rigid_body.Frame().t=TV((T)(2*i-4),height,(T)0);
        rigid_body.Set_Coefficient_Of_Restitution(cor);
        rigid_body.Set_Mass(1);}

    tests.Add_Ground((T).3,0,cor);
}
//#####################################################################
// Function Ether_Test
//#####################################################################
void Ether_Test()
{
    last_frame=240;
    RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("box",(T)1,(T)1);
    box1.Frame().t=TV(-2,20,0);
    box1.Twist().angular=TV(4,4,4);
    box1.Update_Angular_Momentum();
    box1.name="Normal";
    RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("box",(T)1,(T)1);
    box2.Frame().t=TV(2,20,0);
    box2.Twist().angular=TV(4,4,4);
    box2.Update_Angular_Momentum();
    box2.name="Ether";
    tests.Add_Ground();
}
//#####################################################################
// Function Wind_Test
//#####################################################################
void Wind_Test()
{
    last_frame=240;
    RIGID_BODY<TV>& plank1=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank1.Frame().t=TV(-8,20,0);
    plank1.Twist().angular=TV(4,0,0);
    plank1.Update_Angular_Momentum();
    plank1.name="Normal";
    RIGID_BODY<TV>& plank2=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank2.Frame().t=TV(0,20,0);
    plank2.Frame().r=ROTATION<TV>((T)pi/2,TV(0,0,1));
    plank2.Twist().angular=TV(4,0,0);
    plank2.Update_Angular_Momentum();
    plank2.name="Wind";
    RIGID_BODY<TV>& plank3=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank3.Frame().t=TV(8,20,0);
    plank3.Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    plank3.Twist().angular=TV(4,0,0);
    plank3.Update_Angular_Momentum();
    plank3.name="Wind";
    tests.Add_Ground();
}
//#####################################################################
// Function Force_Propogation
//#####################################################################
void Force_Propogation()
{
    last_frame=240;
    RIGID_BODY<TV>* rigid_body=0;
    T ground_mu=(T)0.5;
    T sliding_block_mu=(T)0.5;
    T top_block_mu=(T)2.0;
    T epsilon=1;
    T initial_speed=20;

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,sliding_block_mu);
    rigid_body->Frame().t=TV(0,1,-4);
    rigid_body->Twist().linear=TV(initial_speed,0,0);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,sliding_block_mu);
    rigid_body->Frame().t=TV(0,1,0);
    rigid_body->Twist().linear=TV(initial_speed,0,0);
    rigid_body->Mass()*=2;

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,sliding_block_mu);
    rigid_body->Frame().t=TV(0,1,4);
    rigid_body->Twist().linear=TV(initial_speed,0,0);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,top_block_mu);
    rigid_body->Frame().t=TV(0,3,4);
    rigid_body->Twist().linear=TV(initial_speed,0,0);

    tests.Add_Ground(ground_mu,0,epsilon);
}
//#####################################################################
// Function Suspended_Cube
//#####################################################################
void Suspended_Cube()
{
    last_frame=240;
    // parameters
    T cylinder_radius=4,cube_scale=1,theta=(T)pi/3,m=tan(theta);
    T x=cube_scale+cylinder_radius,y=m*(cylinder_radius+cube_scale)+cylinder_radius*sqrt(m*m+1);
    T friction=1;

    RIGID_BODY<TV>& inclined_plane1=tests.Add_Ground((T).3,0);
    inclined_plane1.Frame().r=ROTATION<TV>(theta,TV(0,0,1));
    inclined_plane1.name="inclined_plane1";
    RIGID_BODY<TV>& inclined_plane2=tests.Add_Ground((T).3,0);
    inclined_plane2.Frame().r=ROTATION<TV>(-theta,TV(0,0,1));
    inclined_plane2.name="inclined_plane2";
    RIGID_BODY<TV>& cylinder1=tests.Add_Rigid_Body("cylinder",cylinder_radius,friction);
    cylinder1.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
    cylinder1.Frame().t=TV(x,y,0);
    RIGID_BODY<TV>& cylinder2=tests.Add_Rigid_Body("cylinder",cylinder_radius,friction);
    cylinder2.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
    cylinder2.Frame().t=TV(-x,y,0);
    RIGID_BODY<TV>& cube=tests.Add_Rigid_Body("subdivided_box",cube_scale,friction);
    cube.Frame().t=TV(0,y,0);
}
//#####################################################################
// Function Sphere_Incline
//#####################################################################
void Sphere_Incline()
{
    last_frame=120;
    T angle=parameter*(T)pi/40;
    RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",1,(T)1e10);
    sphere.Set_Coefficient_Of_Restitution(0);

    RIGID_BODY<TV>& ground=tests.Add_Ground((T)1e10,-1/cos(angle),(T)0);
    ground.Frame().r=ROTATION<TV>(angle,TV(0,0,1));

    // analytic solution: y = -5/14 g t^2 sin(angle)^2
}
//#####################################################################
// Function Spinning_Top
//#####################################################################
void Spinning_Top()
{
    last_frame=240;
    RIGID_BODY<TV>* rigid_body=0;

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->Frame().t=TV(-6,10,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-40,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->name="spinning, falling top";

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->Frame().t=TV(-2,1.5,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-40,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->name="spinning top";

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->Frame().t=TV(2,1.5,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-10,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->name="spinning slow top";

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->Frame().t=TV(6,1.5,0);
    rigid_body->Frame().r=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,0,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->name="stationary (control) top";

    tests.Add_Ground(1,0,1).Set_Coefficient_Of_Rolling_Friction(1);
}
//#####################################################################
// Function Spinning_Cylinder
//#####################################################################
void Spinning_Cylinder()
{
    last_frame=7200;
    RIGID_BODY<TV>* rigid_body=0;
    rigid_body=&tests.Add_Rigid_Body("Rings_Test/medium_cylinder",(T)1,(T).5);
    rigid_body->Twist().angular=TV(1,1,1);
    rigid_body->Update_Angular_Momentum();
}
//#####################################################################
// Function Cluster_Cluster_Interaction
//#####################################################################
void Cluster_Cluster_Interaction()
{
    last_frame=240;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=true;
    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;

    const char *boxfile="subdivided_box";

    // add ground first
    tests.Add_Ground(1,-4,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_1->Frame().t=TV(0,8,0);
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->name="box1";

    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_2->Frame().t=TV(2,9,0);
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->name="box2";

    RIGID_BODY<TV>* rigid_body_3=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_3->Frame().t=TV(4,9,0);
    rigid_body_3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_3->name="box3";

    RIGID_BODY<TV>* rigid_body_4=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_4->Frame().t=TV(4,11,0);
    rigid_body_4->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_4->name="box4";

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_1;
    children_1.Append(rigid_body_1->particle_index);
    children_1.Append(rigid_body_2->particle_index);
    children_1.Append(rigid_body_3->particle_index);
    children_1.Append(rigid_body_4->particle_index);
    int cluster_particle=rigid_bindings.Add_Binding(children_1);
    RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster1";

    RIGID_BODY<TV>* rigid_body_5=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_5->Frame().t=TV(4.5,4,0);
    rigid_body_5->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_5->name="box5";

    RIGID_BODY<TV>* rigid_body_6=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_6->Frame().t=TV(6.5,5,0);
    rigid_body_6->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_6->name="box6";

    RIGID_BODY<TV>* rigid_body_7=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_7->Frame().t=TV(8.5,5,0);
    rigid_body_7->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_7->name="box7";

    RIGID_BODY<TV>* rigid_body_8=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_8->Frame().t=TV(8.5,7,0);
    rigid_body_8->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_8->name="box8";

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_2;
    children_2.Append(rigid_body_5->particle_index);
    children_2.Append(rigid_body_6->particle_index);
    children_2.Append(rigid_body_7->particle_index);
    children_2.Append(rigid_body_8->particle_index);
    cluster_particle=rigid_bindings.Add_Binding(children_2);
    rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster2";

    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_3->particle_index,rigid_body_5->particle_index));
    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_5->particle_index,rigid_body_3->particle_index));

    // add gravity
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_3->particle_index);
    referenced_rigid_particles->Append(rigid_body_4->particle_index);
    referenced_rigid_particles->Append(rigid_body_5->particle_index);
    referenced_rigid_particles->Append(rigid_body_6->particle_index);
    referenced_rigid_particles->Append(rigid_body_7->particle_index);
    referenced_rigid_particles->Append(rigid_body_8->particle_index);
    solid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));
    //solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Cluster_Cluster_Interaction_Kinematic
//#####################################################################
void Cluster_Cluster_Interaction_Kinematic()
{
    last_frame=240;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=true;
    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;

    const char *boxfile="subdivided_box";
    const char *longboxfile="short_plank";
    const char *spherefile="sphere_66k";
    
    // add ground first
    tests.Add_Ground(1,-4,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Rigid_Body(longboxfile,(T)2,(T)0);
    rigid_body_1->is_static=true;
    rigid_body_1->Frame().t=TV(0,4,0);
    rigid_body_1->Frame().r=ROTATION<TV>((T)half_pi,TV(0,0,1));
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->name="longbox1";
    
    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Rigid_Body(longboxfile,(T)3,(T)0);
    rigid_body_2->is_static=true;
    rigid_body_2->Frame().t=TV(0,7,0);
    rigid_body_2->Frame().r=ROTATION<TV>((T)half_pi,TV(0,1,0));
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->name="longbox2";

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_1;
    children_1.Append(rigid_body_1->particle_index);
    children_1.Append(rigid_body_2->particle_index);
    int cluster_particle=rigid_bindings.Add_Binding(children_1);
    RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster1";

    RIGID_BODY<TV>* rigid_body_5=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_5->Frame().t=TV(2,11,0);
    rigid_body_5->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_5->name="box1";

    RIGID_BODY<TV>* rigid_body_6=&tests.Add_Rigid_Body(spherefile,(T)2,(T)0);
    rigid_body_6->Frame().t=TV(-2,11,0);
    rigid_body_6->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_6->name="sphere1";

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_2;
    children_2.Append(rigid_body_5->particle_index);
    children_2.Append(rigid_body_6->particle_index);
    cluster_particle=rigid_bindings.Add_Binding(children_2);
    rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster2";

    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_2->particle_index,rigid_body_5->particle_index));
    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_5->particle_index,rigid_body_2->particle_index));

    // add gravity
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_5->particle_index);
    referenced_rigid_particles->Append(rigid_body_6->particle_index);
    solid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));
    //solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Push_Out
//#####################################################################
void Push_Out()
{
    last_frame=240;
    solids_parameters.rigid_body_collision_parameters.use_push_out_rotation=false;

    tests.Add_Ground(1,-1,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Analytic_Cylinder(5,1);
    rigid_body_1->Frame().r=ROTATION<TV>::From_Euler_Angles(0,0,0);
    //rigid_body_1->Frame().t=TV(0,(T).65,0);
    rigid_body_1->Frame().t=TV(0,(T).8,0);
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->name="cyl1";

    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Analytic_Cylinder(10,(T).7);
    rigid_body_2->Frame().r=ROTATION<TV>::From_Euler_Angles(-(T)pi/100,7*(T)pi/16,0);
    rigid_body_2->Frame().t=TV((T)1.5,(T).1,-(T).9);
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->name="cyl2";

    RIGID_BODY<TV>* rigid_body_3=&tests.Add_Analytic_Cylinder(10,(T).7);
    rigid_body_3->Frame().r=ROTATION<TV>::From_Euler_Angles((T)pi/16,-7*(T)pi/16,0);
    rigid_body_3->Frame().t=TV(2,(T).3,(T).6);
    rigid_body_3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_3->name="cyl3";

    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_3->particle_index);
    solid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Removed_Bodies
//#####################################################################
void Removed_Bodies()
{
    last_frame=240;
    ARRAY<RIGID_BODY<TV>*> bodies;
    for(int i=0;i<20;i++){
        bodies.Append(&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5));
        bodies(i)->Frame().t=TV(0,i,0);
        bodies(i)->Set_Coefficient_Of_Restitution((T).5);
        bodies(i)->name=STRING_UTILITIES::string_sprintf("box-%i",i);}
    for(int i=1;i<=20;i+=2) solid_body_collection.rigid_body_collection.rigid_body_particles.Remove_Body(bodies(i)->particle_index);
    tests.Add_Ground();
}
//#####################################################################
// Function Kinematic_Collision
//#####################################################################
void Kinematic_Collision()
{
    last_frame=240;
    ARRAY<RIGID_BODY<TV>*> bodies;
    for(int i=0;i<2;i++){
        bodies.Append(&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5));
        bodies(i)->Frame().t=TV(0,2*i-1,0);
        bodies(i)->Set_Coefficient_Of_Restitution((T).5);
        bodies(i)->name=STRING_UTILITIES::string_sprintf("box-%i",i);}

    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("sphere_66k",(T).5,(T).5);
    rigid_body->Set_Coefficient_Of_Restitution((T).5);
    rigid_body->name="sphere";
    rigid_body->Is_Kinematic()=true;
    curve.Add_Control_Point(0,TV(0,3,-15));
    curve.Add_Control_Point(5,TV(0,3,15));
    rigid_body->Frame().t=curve.Value(0);

    tests.Add_Ground();
}
//#####################################################################
// Function Drop_Cubes
//#####################################################################
void Drop_Cubes()
{
    last_frame=240;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;

    tests.Add_Ground(1,-1,1);

    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    for(int i=0;i<1;i++){
        RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5);
        children.Append(rigid_body->particle_index);
        rigid_body->Frame().t=TV(0,2,0);
        rigid_body->name=STRING_UTILITIES::string_sprintf("box-%i",i);}

    int cluster_particle=rigid_bindings.Add_Binding(children);
    RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->name="cluster";
    
    tests.Add_Ground(0,0,0);

    solid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,(ARRAY<int>*)0));
}
//#####################################################################
// Function Build_Deforming_Sphere
//#####################################################################
void Build_Deforming_Sphere(KINEMATIC_COLLISION_BODY<TV>* sphere,FRAME<TV>& frame,T time,bool update_positions,bool update_velocities)
{
    typedef VECTOR<int,3> TV_INT;

    T rate=1;
    T radius=rate*time+1;
    int n=20;
    
    TV_INT cells(n,n,n);
    RANGE<TV> range(TV(-radius*1.1,-radius*1.1,-radius*1.1),TV(radius*1.1,radius*1.1,radius*1.1));
    
    if(update_positions)
    {
        frame=FRAME<TV>(TV(1,-2,0));

        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=(LEVELSET_IMPLICIT_OBJECT<TV>*)sphere->implicit_object->object_space_implicit_object;
        sphere->Initialize_Implicit_Object_Levelset(cells,range);

        for(NODE_ITERATOR<TV> iterator(levelset->levelset.grid);iterator.Valid();iterator.Next())
            levelset->levelset.phi(iterator.Node_Index())=iterator.Location().Magnitude()-radius;
        levelset->levelset.Compute_Cell_Minimum_And_Maximum();
    }
    if(update_velocities)
    {
        sphere->velocity_grid->Initialize(cells,range);
        sphere->velocity_field->Resize(sphere->velocity_grid->Domain_Indices());
        for(NODE_ITERATOR<TV> iterator(*sphere->velocity_grid);iterator.Valid();iterator.Next())
            sphere->velocity_field->operator()(iterator.Node_Index())=rate*iterator.Location().Normalized();
    }
}
//#####################################################################
// Function Deforming_Sphere
//#####################################################################
void Deforming_Sphere()
{
    last_frame=240;
    deforming_sphere=new KINEMATIC_COLLISION_BODY<TV>(solid_body_collection.rigid_body_collection,true,new GRID<TV>,new ARRAY<TV,TV_INT>);
    Build_Deforming_Sphere(deforming_sphere,deforming_sphere->Frame(),0,true,true);
    solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(deforming_sphere);
    
    tests.Add_Ground((T).5,0,1);

    int n=6;
    for(int i=-n;i<=n;i++)
        for(int j=-n;j<=n;j++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5);
            rigid_body.coefficient_of_friction=0.2;
            rigid_body.Frame().t=TV(i*2.1,2,j*2.1);
        }
}
//#####################################################################
// Function Cluster_Fracture
//#####################################################################
void Cluster_Fracture()
{
    int num_bodies=24;
    std::ifstream istream("Standard_Tests/bodies_cluster.txt");
    TRIANGULATED_SURFACE<T>* surface=0;    
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    //RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;    
    //RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>* cluster_fracture_callbacks=new RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>(rigid_body_collection,rigid_bindings);
    //cluster_fracture_callbacks->allowed_compresive_strain=-.001;
    //rigid_bindings.callbacks=cluster_fracture_callbacks;
    for(int i=0;i<num_bodies;i++){
        RIGID_BODY_PARTICLES<TV>& particles=solid_body_collection.rigid_body_collection.rigid_body_particles;    
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(solid_body_collection.rigid_body_collection);
        children.Append(rigid_body.particle_index);
        referenced_rigid_particles->Append(rigid_body.particle_index);
        T radius;TV center;/*char dummy;*/istream>>center.x>>center.y>>center.z>>radius;//>>dummy;
        particles.frame(rigid_body.particle_index)=FRAME<TV>(center);
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=DIAGONAL_MATRIX<T,TV::SPIN::m>();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),radius);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        if(!surface) surface=TESSELLATION::Generate_Triangles(sphere);
        rigid_body.Add_Structure(*surface);  
        //if(i==1) rigid_body.is_static=true;
        solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        solid_body_collection.rigid_body_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        solid_body_collection.rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    istream.close();
    tests.Add_Ground((T).5,0,0);
    
    //int cluster_particle=rigid_bindings.Add_Binding(children);
    //RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    //rigid_body_cluster->name="cluster");
    
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Sphere_Sanity_Tests
//#####################################################################
void Sphere_Sanity_Tests()
{
    last_frame=24;
    int n=3;
    if(parameter==2) n=0;
    for(int i=0;i<n;i++) tests.Add_Rigid_Body("sphere",1,1).Frame().t.x=(T)3*i;
    if(parameter==3) solid_body_collection.rigid_body_collection.rigid_body_particles.Remove_Body(2);
    if(parameter==4) solid_body_collection.rigid_body_collection.Deactivate_Body(1);
    if(parameter==4) solid_body_collection.rigid_body_collection.Rigid_Body(1).is_static=true;
}
//#####################################################################
// Function Collision_Contact_Pairs_Test
//#####################################################################
void Collision_Contact_Pairs_Test()
{
    RIGID_BODY<TV>& rigid_body_1=tests.Add_Rigid_Body("box",1,(T).5,true);
    rigid_body_1.Frame().t=TV(0,(T)1,0);
    RIGID_BODY<TV>& rigid_body_2=tests.Add_Rigid_Body("sphere",(T).5,(T).5,true);
    rigid_body_2.Frame().t=TV(0,(T)3.5,0);
    T initial_vel=-24+9.8/72;
    rigid_body_2.Twist().linear=TV(0,initial_vel,0);
    rigid_body_2.Set_Coefficient_Of_Restitution(100);
    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,0,1);
    (void)ground;
}
//#####################################################################
};
}
#endif
