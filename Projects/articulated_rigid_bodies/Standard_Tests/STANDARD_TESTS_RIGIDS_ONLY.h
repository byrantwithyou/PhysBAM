//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_RIGIDS_ONLY
//   1. Point joint between 2 blocks
//#####################################################################
#ifndef __STANDARD_TESTS_RIGIDS_ONLY__
#define __STANDARD_TESTS_RIGIDS_ONLY__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_EXAMPLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_RIGIDS_ONLY:public RIGIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef RIGIDS_EXAMPLE<TV> BASE;
    typedef typename TV::SPIN T_SPIN;

    using BASE::output_directory;using BASE::rigids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::last_frame;
    using BASE::stream_type;using BASE::frame_rate;using BASE::rigid_body_collection;using BASE::test_number;using BASE::parse_args;

    RIGIDS_STANDARD_TESTS<TV> tests;

    STANDARD_TESTS_RIGIDS_ONLY(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,rigid_body_collection)
    {
        rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        rigids_parameters.cfl=(T).9;
        rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
    }

    ~STANDARD_TESTS_RIGIDS_ONLY()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {dt=1;}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    void Update_Rigids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

    void Register_Options() PHYSBAM_OVERRIDE
    {
        BASE::Register_Options();
        parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    }
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
        rigid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    }
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>* rigid_body1=0,*rigid_body2=0;

    ARTICULATED_RIGID_BODY<TV>& arb=rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    if(test_number<=4){
        rigid_body1=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        rigid_body2=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        rigid_body1->Frame().t=TV(0,2,0);
        rigid_body1->Set_Name("parent");
        rigid_body2->Set_Name("child");}

    switch(test_number){
        case 1: // point joint
            last_frame=96;
            rigid_body2->Frame().t=TV(0,4,2);
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,-1)));
            break;
        case 2: // rigid with prismatic translation
            last_frame=96;
            rigid_body2->Frame().t=TV((T)2.5,2,0);
            joint=new RIGID_JOINT<TV>();((RIGID_JOINT<TV>*)joint)->Set_Prismatic_Component_Translation(TV((T).5,0,0));
            arb.joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
            break;
        case 3: // hinge
            last_frame=96;
            rigid_body2->Frame().t=TV(2,4,0);
            joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1),ROTATION<TV>((T)pi/2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            break;
        case 4: // twist
            last_frame=96;
            rigid_body2->Frame().t=TV((T)2.1,2,0);
            joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.1,0,0)));
            rigid_body2->Angular_Momentum()=TV(15,0,0);
            rigid_body2->Update_Angular_Velocity();
            break;
        case 5: // multiple rigid contraints
            last_frame=240;
            for(int i=0;i<8;i++){
                RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
                rigid_body.Frame().t=TV((T)2*i,(T)2*i,(T)2*i);
                rigid_body.Set_Coefficient_Of_Restitution((T).9);
                if(i>0){
                    joint=new RIGID_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body.particle_index-1,rigid_body.particle_index,joint);
                    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
                    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,-1)));}}
            break;
        case 6: // closed loop
            last_frame=240;
            tests.Make_Lathe_Chain(FRAME<TV>(TV(0,10,0)));
            break;
        default:PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    tests.Add_Ground((T).5,-2,1);

    // add forces
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,true));
}

};
}
#endif
