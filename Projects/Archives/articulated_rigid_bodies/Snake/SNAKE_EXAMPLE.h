//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SNAKE_EXAMPLE
//##################################################################### 
#ifndef __SNAKE_EXAMPLE__
#define __SNAKE_EXAMPLE__

#include <Articulated_Rigid_Bodies/BEND_FUNCTION_JOINT.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW>
class SNAKE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
    typedef VECTOR<T,3> TV;
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;

    ARTICULATED_RIGID_BODY<TV>* arb;
    BEND_FUNCTION_JOINT<TV>* bend_function_joint;
    int current_frame;
    bool flip_frame;

    SNAKE_EXAMPLE();
    ~SNAKE_EXAMPLE();
    void Initialize_Bodies();
    void Make_Snake_Chain(TV start_point,TV direction,TV shift,QUATERNION<T> orient,int& num_joints,int& num_bodies,std::string name);
    void Postprocess_Frame(const int frame);
    void Preprocess_Frame(const int frame);
};
}
#endif
