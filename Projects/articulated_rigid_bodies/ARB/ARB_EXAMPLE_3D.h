//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARB_EXAMPLE_3D
//##################################################################### 
#ifndef __ARB_EXAMPLE_3D__
#define __ARB_EXAMPLE_3D__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T>
class ARB_EXAMPLE_3D:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;typedef VECTOR<T,3> TV;
    using BASE::fluids_parameters;using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::last_frame;
    using BASE::stream_type;using BASE::frame_rate;
    RIGID_BODY<TV>* square1,*square2;
    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;

    ARB_EXAMPLE_3D(const STREAM_TYPE stream_type,int parameter);
    ~ARB_EXAMPLE_3D();
    void Initialize_Rigid_Bodies();
    void Two_Twins_Tumbling();
    void Three_Things_Throwing();
    void Four_Fools_Flailing();
    void Five_Phallic_Falling();
    void One_Wonder_Wobbling();
    void Diez_Draping_Diddlehoppers();
    void Twelve_Twirps_Twiddling();
    void Coney_Island();
};
}
#endif
