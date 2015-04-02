//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Kinematic box
//   2. Kinematic box stack
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,1> >
{
    typedef T_input T;typedef VECTOR<T,1> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::parse_args;using BASE::test_number;using BASE::Set_External_Positions;using BASE::data_directory; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;

    int kinematic_body_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(stream_type,data_directory,solid_body_collection)
    {
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=1;
        output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    }

    virtual ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    switch(test_number){
      case 1: Kinematic();break;
      case 2: Kinematic();break;
      default: PHYSBAM_FATAL_ERROR(LOG::sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Kinematic
//#####################################################################
void Kinematic()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    RIGID_BODY<TV>* rigid_body=0;
    T baseboxsize=2;
    T stack_mu=(T)0.5;

    if(test_number!=1){
        T stack_epsilon=(T)0.3;
        T smallboxmass=1;
        T boxsize1=(T)1.0;
        T boxsize2=(T)0.5;
        rigid_body=&tests.Add_Analytic_Box(VECTOR<T,1>(boxsize1));
        rigid_body->coefficient_of_friction=stack_mu;
        rigid_body->Frame().t=TV(2*baseboxsize+boxsize1);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->name="stack box 1a";

        rigid_body=&tests.Add_Analytic_Box(VECTOR<T,1>(boxsize2));
        rigid_body->coefficient_of_friction=stack_mu;
        rigid_body->Frame().t=TV(2*baseboxsize+2*boxsize1+boxsize2);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->name="stack box 1b";}

    //boxfile="square_refined";
    rigid_body=&tests.Add_Analytic_Box(VECTOR<T,1>(baseboxsize));
    rigid_body->coefficient_of_friction=stack_mu;
    rigid_body->Set_Coefficient_Of_Restitution((T)0.1);
    rigid_body->name="base box";
    rigid_body->Is_Kinematic()=true;
    kinematic_body_id=rigid_body->particle_index;

    TV t0(baseboxsize),t1(baseboxsize+(T)17.5);
    ROTATION<TV> r0;
    curve.Add_Control_Point((T).5,FRAME<TV>(t0,r0));
    curve.Add_Control_Point((T)4,FRAME<TV>(t1,r0));
    curve.Add_Control_Point((T)6,FRAME<TV>(t1,r0));
    curve.Add_Control_Point((T)8,FRAME<TV>(t1,r0));

    rigid_body->Frame()=curve.Value(0);
    rigid_body->Twist()=curve.Derivative(0);

    last_frame=(int)(15*frame_rate);

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // add forces
    tests.Add_Gravity();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==1 && id==kinematic_body_id) twist=curve.Derivative(time);
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==1 && id==kinematic_body_id) frame=curve.Value(time);
}
//#####################################################################
};
}
#endif
