//#####################################################################
// Copyright 2004-2008, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BRIDGE_EXAMPLE
//##################################################################### 
#ifndef __BRIDGE_EXAMPLE__
#define __BRIDGE_EXAMPLE__

#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Joints/JOINT.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

class PARSE_ARGS;

template<class T_input>
class BRIDGE_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::viewer_dir;using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::data_directory;
    using BASE::stream_type;using BASE::restart;using BASE::last_frame;using BASE::frame_rate;
    using BASE::test_number;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Position_Nodes; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
    
    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY<TV> *box1,*box2;
    int start_rolling,num_rolling_frames,start_drop; // allow the bridge to settle until start_drop
    int num_rungs;
    int selection;
    bool use_rigid_deformable_evolution_old;

    BRIDGE_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~BRIDGE_EXAMPLE();

    //#####################################################################
    void Initialize_Bodies() override;
    void Update_Solids_Parameters(const T time) override;
    void Preprocess_Frame(const int frame) override;
private:
    RIGID_BODY<TV>& Add_Rigid_Body(const std::string& filename,const T scale,const T cof,const T cor,const FRAME<TV>& frame,const std::string& name,const T mass_scale=0);
    JOINT<TV>& Add_Joint(const int parent_id,const int child_id,JOINT<TV>* joint,const FRAME<TV>& joint_to_child_frame);
    void Make_Bridge();
    void Make_Lathe_Chains();
    void Make_Blocks();
    void Point_Constraint_With_2_Blocks(TV shift);
//#####################################################################
};
}
#endif
