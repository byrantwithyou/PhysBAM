//#####################################################################
// Copyright 2004-2008, Craig Schroeder, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAGNETS_EXAMPLE
//##################################################################### 
#ifndef __MAGNETS_EXAMPLE__
#define __MAGNETS_EXAMPLE__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
namespace PhysBAM{

template<class T_input>
class MAGNETS_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_EXAMPLE<TV> BASE;
public:
    using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::stream_type;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities; // silence -Woverloaded-virtual

    RIGID_BODY<TV>* shelf00,*shelf01,*shelf10,*shelf11;
    int current_frame,start_move,end_move;
    T increment;
    int selection;
    RANDOM_NUMBERS<T> rg;
    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    int parent_id;

    MAGNETS_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~MAGNETS_EXAMPLE();

    // Unused callbacks
    void Preprocess_Frame(const int frame) override {}
    void Postprocess_Frame(const int frame) override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Update_Solids_Parameters(const T time) override {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) override {}

    void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
    void Initialize_Bodies() override;
    void Drop_Letter(std::string letter,int parent_id,TV start,ROTATION<TV> orient,bool stop_in_middle);
    void Random_Letters(int levels,T up_shift);
    std::string Get_Random_Letter();
//#####################################################################
};
}
#include "MAGNETS_EXAMPLE.cpp"
#endif
