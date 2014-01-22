//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Point joint between 2 blocks
//   2. Rigid joint between 2 blocks
//   3. Multiple rigid joints
//   4. PD example
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Joints/ANGLE_JOINT.h>
#include <Rigids/Joints/JOINT_FUNCTION.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <Rigids/Joints/RIGID_JOINT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::parse_args;using BASE::test_number;using BASE::Set_External_Positions;using BASE::data_directory; // silence -Woverloaded-virtual

    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    int njoints;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),njoints(6)
    {
        fluids_parameters.simulate=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=1;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.use_rigid_deformable_contact=false;
    }

    virtual ~STANDARD_TESTS()
    {
        delete arb;
    }

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-njoints",&njoints,"num","number of joints to use");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Iterative_Tolerance((T)1e-4);
    arb->Set_Contact_Level_Iterations(5);
    arb->Set_Shock_Propagation_Level_Iterations(5);
    arb->Set_Poststabilization_Iterations(5);
    arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);

    ROTATION<TV> rotation=ROTATION<TV>::From_Euler_Angles(VECTOR<T,1>((T)30*(T)pi/(T)180));
    switch(test_number){
      case 1: Point_Joint();break;
      case 2: Rigid_Joint();break;
      case 3: Multiple_Rigid_Joints();break;
      case 4: PD_Example();break;
      case 5:{ // cluster with rigid block of same size
          Large_Cluster_Square(FRAME<TV>(TV(0,2),rotation));
          RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("subdivided_box",2,(T).5);
          rigid_body->Frame()=FRAME<TV>(TV(10,2),rotation);
          break;}
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    tests.Add_Ground(.5,-2,1);

    // add forces
    tests.Add_Gravity();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Point_Joint
//#####################################################################
void Point_Joint()
{
    last_frame=96;
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>& rigid_body0=tests.Add_Rigid_Body("square_refined",1,(T).5);
    RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("square_refined",1,(T).5);
    rigid_body0.Frame().t=TV(0,2);
    rigid_body1.Frame().t=TV(2,4);
    rigid_body0.name="parent";
    rigid_body1.name="child";

    joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(rigid_body0.particle_index,rigid_body1.particle_index,joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1)));
}
//#####################################################################
// Function Rigid_Joint
//#####################################################################
void Rigid_Joint()
{
    last_frame=96;
    RIGID_BODY<TV>& rigid_body0=tests.Add_Rigid_Body("square_refined",1,(T).5);
    RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("square_refined",1,(T).5);
    rigid_body0.Frame().t=TV(0,2);
    rigid_body1.Frame().t=TV((T)2.5,2);
    rigid_body0.name="parent";
    rigid_body1.name="child";

    RIGID_JOINT<TV>* joint=new RIGID_JOINT<TV>();arb->joint_mesh.Add_Articulation(rigid_body0.particle_index,rigid_body1.particle_index,joint);
    joint->Set_Prismatic_Component_Translation(TV((T).5,0));
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1)));
}
//#####################################################################
// Function Rigid_Joint
//#####################################################################
void Multiple_Rigid_Joints()
{
    last_frame=480;
    JOINT<TV>* joint=0;
    for(int i=0;i<8;i++){
        RIGID_BODY<TV>&rigid_body=tests.Add_Rigid_Body("square_refined",1,(T)0);
        rigid_body.Frame().t=TV((T)2*i,(T)2*i);
        rigid_body.Set_Coefficient_Of_Restitution(1);
        if(i>0){
            joint=new RIGID_JOINT<TV>();arb->joint_mesh.Add_Articulation(rigid_body.particle_index-1,rigid_body.particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1)));}}
}
//#####################################################################
// Function Rigid_Joint
//#####################################################################
void PD_Example()
{
    last_frame=480;
    JOINT<TV>*joint=0;
    RIGID_BODY<TV>*parent_body=NULL,*child_body=NULL;
    T cheight=(T)0;
    T k_p=(T)1000;

    // Create first body
    parent_body=&tests.Add_Rigid_Body("square_refined",(T).2,(T).5);
    parent_body->Frame().t=TV(cheight,0);
    parent_body->name="parent";
    parent_body->Set_Mass(5);
    parent_body->is_static=true;

    // Add children and joints
    T desired_angle=2*(T)pi/(njoints+1);
    for(int i=0;i<njoints;i++){
        cheight+=1.25;
        child_body=&tests.Add_Rigid_Body("square_refined",(T).2,(T).5);
        child_body->Frame().t=TV(cheight,0);
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);

        joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(child_body->particle_index-1,child_body->particle_index,joint);
        JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*parent_body,*child_body);
        joint->Set_Joint_Function(jfunc);
        joint->joint_function->Set_k_p(k_p);
        joint->joint_function->Set_Target_Angle(ROTATION<TV>::From_Angle(desired_angle));
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0)));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0)));

        // Swap indices
        parent_body=child_body;
        child_body=NULL;
        joint=NULL;}
    
    arb->Use_PD_Actuators();
    arb->global_post_stabilization=true;
}
//#####################################################################
// Function Large_Cluster_Square
//#####################################################################
int Large_Cluster_Square(FRAME<TV>shift_frame,T scale=1)
{
#if 0
    last_frame=480;
    // CLUSTER 1
    ARRAY<RIGID_BODY<TV>*>& bodies=*new ARRAY<RIGID_BODY<TV>*>(4);
    int count=0;
    for(int i=0;i<4;i++){
//        bodies(i)->name=STRING_UTILITIES::string_sprintf("child::%d",bodies(i)));}
        bodies(i)=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        bodies(i)->name=STRING_UTILITIES::string_sprintf("child::%d",bodies(i)));
        solids_parameters.collision_body_list.Add_Body(bodies(i));}
    for(int i=-1;i<=1;i+=2) for(int j=-1;j<=1;j+=2) bodies(++count)->Frame()=shift_frame*FRAME<TV>(TV((T)j,(T)i));
    tests.Add_Gravity();
    int cluster_id=rigid_body_particles.Add_Cluster_Body(&bodies);
    rigid_body_particles.Rigid_Body(cluster_id).name="combo_square");
    return cluster_id;
#endif
    return int();
}
//#####################################################################
};}

#endif
