//#####################################################################
// Copyright 2007-2008, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_RIGIDS_ONLY
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
//  35. Analytic Spheres
//  36. Cluster Fracture
//  37. Break LSV
//#####################################################################
#ifndef __STANDARD_TESTS_RIGIDS_ONLY__
#define __STANDARD_TESTS_RIGIDS_ONLY__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/KINEMATIC_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_CLUSTER_CONSTITUENT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include <fstream>
#include "../RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_RIGIDS_ONLY:public RIGIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T_input,3> TV;
public:
    RIGIDS_STANDARD_TESTS<TV> tests;

    T small_block_mass;
    int parameter;
    INTERPOLATION_CURVE<T,TV> curve_left,curve_right,curve_plank,curve;
    int curve_left_id,curve_right_id,curve_plank_id,torus_id;
    ARRAY<int> rigid_body_particles_with_ether_drag; // test 16

    KINEMATIC_COLLISION_BODY<GRID<TV> >* deforming_sphere;

    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;

    typedef RIGIDS_EXAMPLE<TV> BASE;
    using BASE::rigids_parameters;using BASE::rigid_body_collection;using BASE::rigids_evolution;using BASE::test_number;
    using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::stream_type;using BASE::parse_args;

    STANDARD_TESTS_RIGIDS_ONLY(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,rigid_body_collection),collision_manager(0)
    {
    }

    ~STANDARD_TESTS_RIGIDS_ONLY()
    {if(test_number==21){
        T angle=parameter*(T)pi/40;
        T y=rigid_body_collection.rigid_body_particle.X(1).y,yanalytic=-(T)5./14*(T)9.8*sqr((T)last_frame/24)*sqr(sin(angle));
        T ke=rigid_body_collection.Rigid_Body(1).Kinetic_Energy(),pe=y*rigid_body_collection.rigid_body_particle.mass(1)*(T)9.8;
        T omega=rigid_body_collection.rigid_body_particle.twist(1).angular.z,speed=rigid_body_collection.rigid_body_particle.twist(1).linear.Magnitude();
        LOG::cout<<"ERROR in y for angle "<<angle<<"  ("<<(angle*180/(T)pi)<<" degrees):  "<<(1-y/yanalytic)<<std::endl;
        LOG::cout<<"SLIPPAGE  "<<omega/speed<<std::endl;LOG::cout<<"PE = "<<pe<<"    KE = "<<ke<<std::endl;}}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {
        RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;
        RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>* cluster_fracture=dynamic_cast<RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>*>(rigid_bindings.callbacks);
        if(!cluster_fracture) return;
        for(typename HASHTABLE<int,typename RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::FRACTURE_DATA>::ITERATOR iterator(cluster_fracture->fracture_data);iterator.Valid();iterator.Next()){
            typename RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::FRACTURE_DATA& data=iterator.Data();
            typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_bindings.reverse_bindings.Get(iterator.Key());
            for(int i=data.connections.m;i>=1;i--){
                const VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2>& edge=data.connections(i);
                RIGID_BODY<TV> &child_1=rigid_body_collection.Rigid_Body(cluster.children(edge[1])),&child_2=rigid_body_collection.Rigid_Body(cluster.children(edge[2]));
                VECTOR<int,2> hash_index=VECTOR<int,2>(cluster.children(edge[1]),cluster.children(edge[2])).Sorted();
                T minX=min(child_1.X()(1),child_2.X()(2));
                T& decay=cluster_fracture->decay_rate.Get_Or_Insert(hash_index,0);decay=time-minX;}}
    }

    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

    void Register_Options() PHYSBAM_OVERRIDE
    {
        BASE::Register_Options();
        parse_args->Add_Double_Argument("-small_block_mass",1,"mass for small blocks in plank test");
        parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
        parse_args->Add_Option_Argument("-noanalytic","disable analytic collisions");
        parse_args->Add_Option_Argument("-print_energy","print energy statistics");
        parse_args->Add_Option_Argument("-combined_collisions","use combined collisions and contact");
        parse_args->Add_Option_Argument("-test_combined_system","perform test system for combined collisions");
    }
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);

        small_block_mass=(T)parse_args->Get_Double_Value("-small_block_mass");
        if(small_block_mass!=1) output_directory+=STRING_UTILITIES::string_sprintf("_m%g",small_block_mass);
        parameter=parse_args->Get_Integer_Value("-parameter");
        if(parameter) output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);
        rigids_parameters.rigid_body_collision_parameters.use_analytic_collisions=!parse_args->Get_Option_Value("-noanalytic");
        rigid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
        rigids_parameters.rigid_body_collision_parameters.use_combined_collisions=parse_args->Is_Value_Set("-combined_collisions");
        rigids_parameters.rigid_body_collision_parameters.test_combined_system=parse_args->Is_Value_Set("-test_combined_system");
    }
    
    void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}

//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*rigids_evolution->rigid_body_collisions;
    collisions.Set_Push_Out_Level_Iterations(1);
    if(test_number==26||test_number==27) collisions.collision_manager=collision_manager;
    if(test_number==35||test_number==36||test_number==37){
        int level=1;
        rigids_parameters.rigid_body_collision_parameters.collision_iterations=level;
        collisions.Set_Collision_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.contact_iterations=level;
        collisions.Set_Contact_Level_Iterations(level);
        collisions.Set_Contact_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
        collisions.Set_Shock_Propagation_Iterations(level);
        collisions.Set_Shock_Propagation_Level_Iterations(level);
        collisions.Set_Shock_Propagation_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
        rigids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        collisions.Set_Push_Out_Iterations(level);
        collisions.Set_Push_Out_Level_Iterations(level);
        collisions.Set_Push_Out_Pair_Iterations(level*2);
        collisions.Register_Analytic_Collisions();}
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
        Build_Deforming_Sphere(deforming_sphere,frame.t,frame.r,time,true,false);
    if(test_number==36)
        frame=FRAME<TV>(TV(-100+100*time,125,-100+100*time),ROTATION<TV>());
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==14){
        if(id==1) twist.linear=curve_left.Derivative(time);
        else if(id==int(2)) twist.linear=curve_right.Derivative(time);
        else if(id==int(3)) twist.linear=curve_plank.Derivative(time);}
    else if(test_number==30) twist.linear=curve.Derivative(time);
    else if(test_number==34)
        Build_Deforming_Sphere(deforming_sphere,deforming_sphere->X(),deforming_sphere->Rotation(),time,false,true);
    else return false;
    return true;
}
//#####################################################################
// Function Update_Rigids_Parameters
//#####################################################################
void Update_Rigids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==19){
        bool &is_static=rigid_body_collection.Rigid_Body(5).is_static,new_is_static=time<=(T).5;
        if(is_static!=new_is_static){is_static=new_is_static;rigid_body_collection.Update_Simulated_Particles();}} // Changing is_static requires an Update_Simulated_Particles
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==31){
        if(parameter==4 && time>=(T).4999 && time<=(T).5001){
            rigid_body_collection.rigid_geometry_collection.Reactivate_Geometry(2);
            rigid_body_collection.Update_Simulated_Particles();}
        if(parameter==3 && time>=(T).4999 && time<=(T).5001){
            rigid_body_collection.rigid_body_particle.Remove_Body(1);
            rigid_body_collection.Update_Simulated_Particles();}}
}
//#####################################################################
// Function Set_Rigid_Particle_Is_Simulated
//#####################################################################
void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
    if(test_number==31 && parameter==6) particle_is_simulated(2)=false;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    rigids_parameters.cfl=1;
    rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
    rigids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;

    switch(test_number){
      case 1: Bounce(0);break;
      case 2: Nonconvex_Bounce();break;
      case 3: Billiard_Balls();break;
      case 4: Bounce_Friction();break;
      case 5: Rolling_Friction();break;
      case 7: rigids_parameters.rigid_body_collision_parameters.contact_iterations=0; //fall through
      case 8: Seesaw();break;
      case 9: rigids_parameters.rigid_body_collision_parameters.contact_iterations=0; //fall through
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
      case 35: Analytic_Contact();break;
      case 36: Cluster_Fracture();break;
      case 37: Break_Levelset();break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
    
    // add forces
    if(test_number==16){
        rigid_body_particles_with_ether_drag.Append(2);
        rigid_body_collection.Add_Force(new RIGID_ETHER_DRAG<GRID<TV> >(rigid_body_collection,&rigid_body_particles_with_ether_drag,(T).5,(T).5));}
    else if(test_number==17){
        for(int i=2;i<=3;i++){
            PHYSBAM_NOT_IMPLEMENTED("RIGID_BODY_BASIC_FORCES is obsolete");
//             RIGID_BODY_BASIC_FORCES<TV>* basic_forces=0; // do rigid body i
//             basic_forces->Use_Wind_Drag();
//             basic_forces->Set_Wind_Density((T).1);
//             basic_forces->Use_Constant_Wind(0,TV(1,0,0));
        }}
    if(test_number!=23 && test_number!=26 && test_number!=27 && test_number!=28) rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,true));
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

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,0,1,1);
    ground.Rotation()=ROTATION<TV>(angle,TV(0,0,1));
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
        rigid_body.X()=TV(x_pos[i],5,0);
        rigid_body.Rotation()=ROTATION<TV>(z_angle[i],TV(0,0,1))*ROTATION<TV>((T)-0.5,TV(1,0,0))*ROTATION<TV>((T)pi,TV(0,1,0));
        rigid_body.Set_Coefficient_Of_Restitution(cor[i]);
        rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("bone (cor %g)",cor[i]));}

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
        rigid_body.X()=TV((T)(k?2*(k-(number+1)/2):-30),1,0);
        if(!k) rigid_body.Twist().linear=TV(20,0,0);
        rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("sphere%d",k));
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
    sphere.X()=TV(-4,10,0);
    sphere.Twist().angular=TV(-10,0,0);
    sphere.Update_Angular_Momentum();
    sphere.Set_Coefficient_Of_Restitution((T)0.5);
    sphere.Set_Coefficient_Of_Rolling_Friction((T)0.01);
    sphere.Set_Name("spinning sphere");

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
        rigid_body->X()=TV(-2,1,45);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Coefficient_Of_Rolling_Friction((T)0.01);
        rigid_body->Set_Name("sphere (mu=0.05)");

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).1);
        rigid_body->X()=TV(2,1,45);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Twist().angular=TV(0,0,-10);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Coefficient_Of_Rolling_Friction((T)0.2);
        rigid_body->Set_Name("sphere (mu=0.2)");}
    else{
        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).05);
        rigid_body->X()=TV(-8,1,25);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Set_Name("sphere (mu=0.05)");

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).02);
        rigid_body->X()=TV(8,1,25);
        rigid_body->Twist().linear=TV(0,0,-10);
        rigid_body->Set_Name("sphere (mu=0.02)");

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).5);
        rigid_body->X()=TV(-2,10,0);
        rigid_body->Twist().angular=TV(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("spinning sphere");

        rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).05);
        rigid_body->X()=TV(2,10,0);
        rigid_body->Twist().angular=TV(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Name("low friction spinning sphere");}

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
    rigid_body->X()=TV((T)3.5,fromrest?(T)2.4:20,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Name("box1");

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T)0.5,mu);
    rigid_body->Set_Mass(5);
    rigid_body->X()=TV(-4,fromrest?(T)2:(T)3,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Name("box2");

    rigid_body=&tests.Add_Rigid_Body("plank",1,1);
    rigid_body->X()=TV(0,fromrest?(T)1.25:(T)1.5,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Mass(10);
    rigid_body->Set_Name("plank");

    rigid_body=&tests.Add_Rigid_Body("Rings_Test/cylinder_revolve",1,mu);
    rigid_body->X()=TV(0,(T)0.5,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Name("cylinder");
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
        rigid_body->X()=TV(x_pos[i],baseboxsize+y_pos[i],z_pos[i]);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("stack box %s",name[i]));}
    rigid_body->Rotation()=ROTATION<TV>((T)pi/4,TV(0,1,0));

    rigid_body=&tests.Add_Rigid_Body(boxfile,1,1);
    rigid_body->X()=TV(baseboxsize+plankscale*5-dropboxsize,2*baseboxsize+(T)0.5*plankscale+dropboxsize+dropheight,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.05);
    rigid_body->Set_Mass(100);
    rigid_body->Set_Name("drop box");

    rigid_body=&tests.Add_Rigid_Body(plankfile,plankscale,1);
    rigid_body->X()=TV(baseboxsize,2*baseboxsize+(T)0.25,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->Set_Mass(10);
    rigid_body->Set_Name("plank");

    rigid_body=&tests.Add_Rigid_Body(boxfile,baseboxsize,1);
    rigid_body->X()=TV(0,baseboxsize,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->Set_Coefficient_Of_Friction((T)0.3);
    rigid_body->Set_Name("base box");
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
        Random_Scene_Generator("Rings_Test/ring_revolve",500,12345,random_placement,rigid_body_collection,tests);}
    else{
        LOG::cout<<"Parameter 1: running 1000 rings"<<std::endl;
        CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(20,600,TV(0,30,0));
        random_placement.Set_Max_Orientation_Angle((T)0.2);
        random_placement.Set_Speed_Range(0,1);
        random_placement.Set_Angular_Speed_Range(0,1);
        Random_Scene_Generator("Rings_Test/ring_revolve",1000,11111,random_placement,rigid_body_collection,tests);}

    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Mass(10);
        rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,1,3));}

    for(int i=1;i<=poles;i++)
        for(int j=1;j<=poles;j++){
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
    last_frame=240;
    RIGID_BODY<TV>* rigid_body=0;

    ARRAY<std::string>filenames;
    for(int i=1;i<=300;i++) filenames.Append("New_Bones/Cranium_1");
    for(int i=1;i<=150;i++) filenames.Append("New_Bones/Pelvis_1");
    for(int i=1;i<=50;i++) filenames.Append("New_Bones/Left_Femur_2");

    T mu=(T)0.4;
    T epsilon=(T)0.3;
    T rolling_friction=(T)0.1;

    CYLINDRICAL_RANDOM_PLACEMENT<TV> random_placement(1,40,TV(0,3,0));
    random_placement.Set_Fixed_Scale(3);
    Random_Scene_Generator(filenames,12345,random_placement,rigid_body_collection,tests);
    T mass_scale=10;
    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        rigid_body.Set_Mass(rigid_body.Mass()*mass_scale);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Coefficient_Of_Rolling_Friction(rolling_friction);

        if(rigid_body.name.find("Cranium"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(2,2,2));
        else if(rigid_body.name.find("Left_Femur"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,4,1));
        else if(rigid_body.name.find("Pelvis"))
            rigid_body.simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,3,2));}

    rigid_body=&tests.Add_Rigid_Body("funnel_thicker_revolve",(T)0.4,1);
    rigid_body->X()=TV(0,2,0);
    rigid_body->Set_Name("funnel");
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(3,3,3));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    T yy=(T)0.2;
    T ss=(T)0.25;
    T xx=ss*(T)4.75;

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->X()=TV(xx,yy,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->Set_Name("plank 1");
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->X()=TV(0,yy,xx);
    rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->Set_Name("plank 2");
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->X()=TV(-xx,yy,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->Set_Name("plank 3");
    rigid_body->is_static=true;
    rigid_body->simplicial_object->Set_Desired_Particle_Partition_Size(VECTOR<int,3>(1,10,1));
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);

    rigid_body=&tests.Add_Rigid_Body("plank",ss,1);
    rigid_body->X()=TV(0,yy,-xx);
    rigid_body->Rotation()=ROTATION<TV>(3*(T)pi/2,TV(0,1,0))*rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,0,1));
    rigid_body->Set_Name("plank 4");
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
    ram1.X()=curve_left.Value(0);
    ram1.Twist().linear=curve_left.Derivative(0);

    RIGID_BODY<TV>& ram2=tests.Add_Rigid_Body("subdivided_box",1,mu);
    ram2.Set_Coefficient_Of_Restitution(cor);
    ram2.Is_Kinematic()=true;
    curve_right_id=ram2.particle_index;
    curve_right.Add_Control_Point((T)0,TV(7,height,0));
    curve_right.Add_Control_Point((T)1,TV(6-(T)1e-4,height,(T)0));
    curve_right.Add_Control_Point((T)5,TV(6-(T)1e-4,height,(T)0));
    curve_right.Add_Control_Point((T)6,TV(7,height,(T)0));
    ram2.X()=curve_right.Value(0);
    ram2.Twist().linear=curve_right.Derivative(0);

    RIGID_BODY<TV>& plank=tests.Add_Rigid_Body("plank",1,mu);
    plank.Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0));
    plank.Set_Coefficient_Of_Restitution(cor);
    plank.Is_Kinematic()=true;
    curve_plank_id=plank.particle_index;
    curve_plank.Add_Control_Point((T)0,TV(0,height-(T)1.25,0));
    curve_plank.Add_Control_Point((T)2,TV(0,height-(T)1.25,(T)0));
    curve_plank.Add_Control_Point((T)3,TV(0,1,(T)0));
    curve_plank.Add_Control_Point((T)4,TV(0,1,(T)4));
    curve_plank.Add_Control_Point((T)6,TV(0,1,(T)4));
    plank.X()=curve_plank.Value(0);
    plank.Twist().linear=curve_plank.Derivative(0);

    for(int i=0;i<5;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,mu);
        rigid_body.X()=TV((T)(2*i-4),height,(T)0);
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
    box1.X()=TV(-2,20,0);
    box1.Twist().angular=TV(4,4,4);
    box1.Update_Angular_Momentum();
    box1.Set_Name("Normal");
    RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("box",(T)1,(T)1);
    box2.X()=TV(2,20,0);
    box2.Twist().angular=TV(4,4,4);
    box2.Update_Angular_Momentum();
    box2.Set_Name("Ether");
    tests.Add_Ground();
}
//#####################################################################
// Function Wind_Test
//#####################################################################
void Wind_Test()
{
    last_frame=240;
    RIGID_BODY<TV>& plank1=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank1.X()=TV(-8,20,0);
    plank1.Twist().angular=TV(4,0,0);
    plank1.Update_Angular_Momentum();
    plank1.Set_Name("Normal");
    RIGID_BODY<TV>& plank2=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank2.X()=TV(0,20,0);
    plank2.Rotation()=ROTATION<TV>((T)pi/2,TV(0,0,1));
    plank2.Twist().angular=TV(4,0,0);
    plank2.Update_Angular_Momentum();
    plank2.Set_Name("Wind");
    RIGID_BODY<TV>& plank3=tests.Add_Rigid_Body("plank",(T)1,(T)1);
    plank3.X()=TV(8,20,0);
    plank3.Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0));
    plank3.Twist().angular=TV(4,0,0);
    plank3.Update_Angular_Momentum();
    plank3.Set_Name("Wind");
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
    rigid_body->X()=TV(0,1,-4);
    rigid_body->Twist().linear=TV(initial_speed,0,0);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,sliding_block_mu);
    rigid_body->X()=TV(0,1,0);
    rigid_body->Twist().linear=TV(initial_speed,0,0);
    rigid_body->Mass()*=2;

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,sliding_block_mu);
    rigid_body->X()=TV(0,1,4);
    rigid_body->Twist().linear=TV(initial_speed,0,0);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,top_block_mu);
    rigid_body->X()=TV(0,3,4);
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
    inclined_plane1.Rotation()=ROTATION<TV>(theta,TV(0,0,1));
    inclined_plane1.Set_Name("inclined_plane1");
    RIGID_BODY<TV>& inclined_plane2=tests.Add_Ground((T).3,0);
    inclined_plane2.Rotation()=ROTATION<TV>(-theta,TV(0,0,1));
    inclined_plane2.Set_Name("inclined_plane2");
    RIGID_BODY<TV>& cylinder1=tests.Add_Rigid_Body("cylinder",cylinder_radius,friction);
    cylinder1.Rotation()=ROTATION<TV>((T)pi/2,TV(1,0,0));
    cylinder1.X()=TV(x,y,0);
    RIGID_BODY<TV>& cylinder2=tests.Add_Rigid_Body("cylinder",cylinder_radius,friction);
    cylinder2.Rotation()=ROTATION<TV>((T)pi/2,TV(1,0,0));
    cylinder2.X()=TV(-x,y,0);
    RIGID_BODY<TV>& cube=tests.Add_Rigid_Body("subdivided_box",cube_scale,friction);
    cube.X()=TV(0,y,0);
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
    ground.Rotation()=ROTATION<TV>(angle,TV(0,0,1));

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
    rigid_body->X()=TV(-6,10,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-40,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("spinning, falling top");

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->X()=TV(-2,1.5,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-40,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("spinning top");

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->X()=TV(2,1.5,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,-10,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("spinning slow top");

    rigid_body=&tests.Add_Rigid_Body("cone",(T)1,(T).5);
    rigid_body->X()=TV(6,1.5,0);
    rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,0,1));
    rigid_body->Twist().angular=TV(0,0,0);
    rigid_body->Update_Angular_Momentum();
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("stationary (control) top");

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
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=true;
    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;

    const char *boxfile="subdivided_box";

    // add ground first
    tests.Add_Ground(1,-4,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_1->X()=TV(0,8,0);
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->Set_Name("box1");

    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_2->X()=TV(2,9,0);
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->Set_Name("box2");

    RIGID_BODY<TV>* rigid_body_3=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_3->X()=TV(4,9,0);
    rigid_body_3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_3->Set_Name("box3");

    RIGID_BODY<TV>* rigid_body_4=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_4->X()=TV(4,11,0);
    rigid_body_4->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_4->Set_Name("box4");

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_1;
    children_1.Append(rigid_body_1->particle_index);
    children_1.Append(rigid_body_2->particle_index);
    children_1.Append(rigid_body_3->particle_index);
    children_1.Append(rigid_body_4->particle_index);
    int cluster_particle=rigid_bindings.Add_Binding(children_1);
    RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster1");

    RIGID_BODY<TV>* rigid_body_5=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_5->X()=TV(4.5,4,0);
    rigid_body_5->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_5->Set_Name("box5");

    RIGID_BODY<TV>* rigid_body_6=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_6->X()=TV(6.5,5,0);
    rigid_body_6->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_6->Set_Name("box6");

    RIGID_BODY<TV>* rigid_body_7=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_7->X()=TV(8.5,5,0);
    rigid_body_7->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_7->Set_Name("box7");

    RIGID_BODY<TV>* rigid_body_8=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_8->X()=TV(8.5,7,0);
    rigid_body_8->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_8->Set_Name("box8");

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_2;
    children_2.Append(rigid_body_5->particle_index);
    children_2.Append(rigid_body_6->particle_index);
    children_2.Append(rigid_body_7->particle_index);
    children_2.Append(rigid_body_8->particle_index);
    cluster_particle=rigid_bindings.Add_Binding(children_2);
    rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster2");

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
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
    //solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Cluster_Cluster_Interaction_Kinematic
//#####################################################################
void Cluster_Cluster_Interaction_Kinematic()
{
    last_frame=240;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=true;
    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;

    const char *boxfile="subdivided_box";
    const char *longboxfile="short_plank";
    const char *spherefile="sphere_66k";
    
    // add ground first
    tests.Add_Ground(1,-4,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Rigid_Body(longboxfile,(T)2,(T)0);
    rigid_body_1->is_static=true;
    rigid_body_1->X()=TV(0,4,0);
    rigid_body_1->Rotation()=ROTATION<TV>((T)half_pi,TV(0,0,1));
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->Set_Name("longbox1");
    
    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Rigid_Body(longboxfile,(T)3,(T)0);
    rigid_body_2->is_static=true;
    rigid_body_2->X()=TV(0,7,0);
    rigid_body_2->Rotation()=ROTATION<TV>((T)half_pi,TV(0,1,0));
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->Set_Name("longbox2");

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_1;
    children_1.Append(rigid_body_1->particle_index);
    children_1.Append(rigid_body_2->particle_index);
    int cluster_particle=rigid_bindings.Add_Binding(children_1);
    RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster1");

    RIGID_BODY<TV>* rigid_body_5=&tests.Add_Rigid_Body(boxfile,(T)1,(T)0);
    rigid_body_5->X()=TV(2,11,0);
    rigid_body_5->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_5->Set_Name("box1");

    RIGID_BODY<TV>* rigid_body_6=&tests.Add_Rigid_Body(spherefile,(T)2,(T)0);
    rigid_body_6->X()=TV(-2,11,0);
    rigid_body_6->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_6->Set_Name("sphere1");

    // make clustered object
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children_2;
    children_2.Append(rigid_body_5->particle_index);
    children_2.Append(rigid_body_6->particle_index);
    cluster_particle=rigid_bindings.Add_Binding(children_2);
    rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster2");

    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_2->particle_index,rigid_body_5->particle_index));
    collision_manager->hash.Insert(PAIR<int,int>(rigid_body_5->particle_index,rigid_body_2->particle_index));

    // add gravity
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_5->particle_index);
    referenced_rigid_particles->Append(rigid_body_6->particle_index);
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
    //solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Push_Out
//#####################################################################
void Push_Out()
{
    last_frame=240;
    rigids_parameters.rigid_body_collision_parameters.use_push_out_rotation=false;

    tests.Add_Ground(1,-1,1);

    RIGID_BODY<TV>* rigid_body_1=&tests.Add_Analytic_Cylinder(5,1);
    rigid_body_1->Rotation()=ROTATION<TV>::From_Euler_Angles(0,0,0);
    //rigid_body_1->X()=TV(0,(T).65,0);
    rigid_body_1->X()=TV(0,(T).8,0);
    rigid_body_1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_1->Set_Name("cyl1");

    RIGID_BODY<TV>* rigid_body_2=&tests.Add_Analytic_Cylinder(10,(T).7);
    rigid_body_2->Rotation()=ROTATION<TV>::From_Euler_Angles(-(T)pi/100,7*(T)pi/16,0);
    rigid_body_2->X()=TV((T)1.5,(T).1,-(T).9);
    rigid_body_2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_2->Set_Name("cyl2");

    RIGID_BODY<TV>* rigid_body_3=&tests.Add_Analytic_Cylinder(10,(T).7);
    rigid_body_3->Rotation()=ROTATION<TV>::From_Euler_Angles((T)pi/16,-7*(T)pi/16,0);
    rigid_body_3->X()=TV(2,(T).3,(T).6);
    rigid_body_3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body_3->Set_Name("cyl3");

    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(rigid_body_1->particle_index);
    referenced_rigid_particles->Append(rigid_body_2->particle_index);
    referenced_rigid_particles->Append(rigid_body_3->particle_index);
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Removed_Bodies
//#####################################################################
void Removed_Bodies()
{
    last_frame=240;
    ARRAY<RIGID_BODY<TV>*> bodies;
    for(int i=1;i<=20;i++){
        bodies.Append(&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5));
        bodies(i)->X()=TV(0,i,0);
        bodies(i)->Set_Coefficient_Of_Restitution((T).5);
        bodies(i)->Set_Name(STRING_UTILITIES::string_sprintf("box-%i",i));}
    for(int i=1;i<=20;i+=2) rigid_body_collection.rigid_body_particle.Remove_Body(bodies(i)->particle_index);
    tests.Add_Ground();
}
//#####################################################################
// Function Kinematic_Collision
//#####################################################################
void Kinematic_Collision()
{
    last_frame=240;
    ARRAY<RIGID_BODY<TV>*> bodies;
    for(int i=1;i<=2;i++){
        bodies.Append(&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5));
        bodies(i)->X()=TV(0,2*i-1,0);
        bodies(i)->Set_Coefficient_Of_Restitution((T).5);
        bodies(i)->Set_Name(STRING_UTILITIES::string_sprintf("box-%i",i));}

    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("sphere_66k",(T).5,(T).5);
    rigid_body->Set_Coefficient_Of_Restitution((T).5);
    rigid_body->Set_Name("sphere");
    rigid_body->Is_Kinematic()=true;
    curve.Add_Control_Point(0,TV(0,3,-15));
    curve.Add_Control_Point(5,TV(0,3,15));
    rigid_body->X()=curve.Value(0);

    tests.Add_Ground();
}
//#####################################################################
// Function Drop_Cubes
//#####################################################################
void Drop_Cubes()
{
    last_frame=240;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;

    tests.Add_Ground(1,-1,1);

    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    for(int i=0;i<1;i++){
        RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5);
        children.Append(rigid_body->particle_index);
        rigid_body->X()=TV(0,2,0);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("box-%i",i));}

    int cluster_particle=rigid_bindings.Add_Binding(children);
    RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster");
    
    tests.Add_Ground(0,0,0);

    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,(ARRAY<int>*)0));
}
//#####################################################################
// Function Build_Deforming_Sphere
//#####################################################################
    void Build_Deforming_Sphere(KINEMATIC_COLLISION_BODY<GRID<TV> >* sphere,TV& X,ROTATION<TV>& rotation,T time,bool update_positions,bool update_velocities)
{
    typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;

    T rate=1;
    T radius=rate*time+1;
    int n=20;
    
    TV_INT cells(n,n,n);
    RANGE<TV> range(TV(-radius*1.1,-radius*1.1,-radius*1.1),TV(radius*1.1,radius*1.1,radius*1.1));
    
    if(update_positions)
    {
        X=TV(1,-2,0);
        rotation=ROTATION<TV>();

        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=(LEVELSET_IMPLICIT_OBJECT<TV>*)sphere->implicit_object->object_space_implicit_object;
        sphere->Initialize_Implicit_Object_Levelset(cells,range);

        for(NODE_ITERATOR iterator(levelset->levelset.grid);iterator.Valid();iterator.Next())
            levelset->levelset.phi(iterator.Node_Index())=iterator.Location().Magnitude()-radius;
        levelset->levelset.Compute_Cell_Minimum_And_Maximum();
    }
    if(update_velocities)
    {
        sphere->velocity_grid->Initialize(cells,range);
        sphere->velocity_field->Resize(sphere->velocity_grid->Domain_Indices());
        for(NODE_ITERATOR iterator(*sphere->velocity_grid);iterator.Valid();iterator.Next())
            sphere->velocity_field->operator()(iterator.Node_Index())=rate*iterator.Location().Normalized();
    }
}
//#####################################################################
// Function Deforming_Sphere
//#####################################################################
void Deforming_Sphere()
{
    last_frame=240;
    deforming_sphere=new KINEMATIC_COLLISION_BODY<GRID<TV> >(rigid_body_collection,true,new GRID<TV>,new typename KINEMATIC_COLLISION_BODY<GRID<TV> >::T_ARRAYS_TV);
    Build_Deforming_Sphere(deforming_sphere,deforming_sphere->X(),deforming_sphere->Rotation(),0,true,true);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(deforming_sphere);
    
    tests.Add_Ground((T).5,0,1);

    int n=6;
    for(int i=-n;i<=n;i++)
        for(int j=-n;j<=n;j++)
        {
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T)1,(T).5);
            rigid_body.Set_Coefficient_Of_Friction(0.2);
            rigid_body.X()=TV(i*2.1,2,j*2.1);
        }
}
//#####################################################################
// Function Analytic_Contact
//#####################################################################
void Analytic_Contact()
{
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;    
    int num_bodies=2406;
    std::ifstream istream("Standard_Tests/bodies.txt");
    TRIANGULATED_SURFACE<T>* surface=0;    
    for(int i=1;i<=num_bodies;i++){
        RIGID_BODY_PARTICLES<TV>& particles=rigid_body_collection.rigid_body_particle;    
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        T radius;TV center;char dummy;istream>>center.x>>center.y>>center.z>>radius>>dummy;
        particles.X(rigid_body.particle_index)=center;
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),radius);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        if(!surface) surface=TESSELLATION::Generate_Triangles(sphere);
        rigid_body.Add_Structure(*surface);  
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    istream.close();
    tests.Add_Ground((T).5,0,0);
    RIGID_GRAVITY<TV> *gravity=new RIGID_GRAVITY<TV>(rigid_body_collection,true);
    rigid_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Cluster_Fracture
//#####################################################################
void Cluster_Fracture()
{
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;    
    int num_bodies=1000;
    std::ifstream istream("Standard_Tests/bodies_cluster.txt");
    TRIANGULATED_SURFACE<T>* surface=0;    
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    RIGID_BODY_PARTICLES<TV>& particles=rigid_body_collection.rigid_body_particle;    
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=false;
    RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>* cluster_fracture_callbacks=new RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>(rigid_body_collection,rigid_bindings);
    cluster_fracture_callbacks->allowed_strain=.001;
    rigid_bindings.callbacks=cluster_fracture_callbacks;
    for(int i=1;i<=num_bodies;i++){
        T radius;TV center;/*char dummy;*/istream>>center.x>>center.y>>center.z>>radius;//>>dummy;
        if(center.y<100) continue;
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        children.Append(rigid_body.particle_index);
        referenced_rigid_particles->Append(rigid_body.particle_index);
        particles.X(rigid_body.particle_index)=center;
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),radius);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        surface=TESSELLATION::Generate_Triangles(sphere);
        rigid_body.Add_Structure(*surface);  
        if(rigid_body.particle_index==251) rigid_body.is_static=true;
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    istream.close();
    tests.Add_Ground((T).5,0,0);
    
    int cluster_particle=rigid_bindings.Add_Binding(children);
    RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
    rigid_body_cluster->Set_Name("cluster");
    
    {
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        rigid_body.Is_Kinematic()=true;
        particles.X(rigid_body.particle_index)=TV(-100,125,-100);
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),10);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        surface=TESSELLATION::Generate_Triangles(sphere);
        rigid_body.Add_Structure(*surface);  
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;
    }

    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Cluster_Fracture
//#####################################################################
void Break_Levelset()
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    bool use_clustering=true;
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    RIGID_BODY_PARTICLES<TV>& particles=rigid_body_collection.rigid_body_particle;    
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.collide_constituent_bodies=false;
    if(use_clustering){
        RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>* cluster_fracture_callbacks=new RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>(rigid_body_collection,rigid_bindings);
        cluster_fracture_callbacks->allowed_strain=1;
        rigid_bindings.callbacks=cluster_fracture_callbacks;}
    GRID<TV>& grid=*new GRID<TV>;
    ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >;
    LEVELSET_3D<GRID<TV> > levelset(grid,phi);
    data_directory="../../../Public_Data/Archives/";
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Rigid_Bodies/sphere_66k.phi",data_directory.c_str()),levelset);
    TV_INT counts=TV_INT::All_Ones_Vector()*10;counts(2)=10;
    GRID<TV> body_grid(counts,levelset.grid.Domain());
    TRIANGULATED_SURFACE<T>* surface=0;    
    for(typename GRID<TV>::CELL_ITERATOR iterator(body_grid);iterator.Valid();iterator.Next()){
        TV_INT cell=levelset.grid.Cell(iterator.Location(),0);
        if(!phi.Valid_Index(cell)) continue;
        if(phi(cell)>0) continue;
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        children.Append(rigid_body.particle_index);
        referenced_rigid_particles->Append(rigid_body.particle_index);
        particles.X(rigid_body.particle_index)=iterator.Location();
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        BOX<TV> box(body_grid.dX*-.5,body_grid.dX*.5);
        //SPHERE<TV> box(TV(),body_grid.min_dX*.5);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(box));
        //rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(box));
        if(!surface) surface=TESSELLATION::Generate_Triangles(box);
        rigid_body.Add_Structure(*surface);  
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    tests.Add_Ground((T).5,-2,0);
    
    if(use_clustering){
        int cluster_particle=rigid_bindings.Add_Binding(children);
        RIGID_BODY<TV>* rigid_body_cluster=&rigid_body_collection.Rigid_Body(cluster_particle);
        rigid_body_cluster->Set_Name("cluster");}
    
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Sphere_Sanity_Tests
//#####################################################################
void Sphere_Sanity_Tests()
{
    last_frame=24;
    int n=3;
    if(parameter==2) n=0;
    for(int i=1;i<=n;i++) tests.Add_Rigid_Body("sphere",1,1).X().x=(T)3*i;
    if(parameter==3) rigid_body_collection.rigid_body_particle.Remove_Body(3);
    if(parameter==4) rigid_body_collection.rigid_geometry_collection.Deactivate_Geometry(2);
}
//#####################################################################
// Function Collision_Contact_Pairs_Test
//#####################################################################
void Collision_Contact_Pairs_Test()
{
    RIGID_BODY<TV>& rigid_body_1=tests.Add_Rigid_Body("box",1,(T).5,true);
    rigid_body_1.X()=TV(0,(T)1,0);
    RIGID_BODY<TV>& rigid_body_2=tests.Add_Rigid_Body("sphere",(T).5,(T).5,true);
    rigid_body_2.X()=TV(0,(T)3.5,0);
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
