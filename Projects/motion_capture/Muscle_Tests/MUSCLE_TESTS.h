//#####################################################################
// Copyright 2008, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_TESTS
//#####################################################################
//   1. 
//#####################################################################
#ifndef __MUSCLE_TESTS__
#define __MUSCLE_TESTS__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class MUSCLE_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;

    MUSCLE_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection)
    {
    }

    ~MUSCLE_TESTS()
    {}
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
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
};
}
#endif
