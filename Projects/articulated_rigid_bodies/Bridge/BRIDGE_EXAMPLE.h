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
    using BASE::output_directory;using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::write_last_frame;using BASE::data_directory;
    using BASE::stream_type;using BASE::restart;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::restart_frame;using BASE::frame_rate;
    using BASE::test_number;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Position_Nodes; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY<TV> *box1,*box2;
    int start_rolling,num_rolling_frames,start_drop; // allow the bridge to settle until start_drop
    int num_rungs;
    int selection;
    bool use_rigid_deformable_evolution_old;

    BRIDGE_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~BRIDGE_EXAMPLE();

    // Unused callbacks
    void Postprocess_Frame(const int frame) override {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) override {}

    void After_Initialization() override {BASE::After_Initialization();}
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
