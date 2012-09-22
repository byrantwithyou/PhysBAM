//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "ARB_Example/ARB_EXAMPLE.h"
//#include "Little_Man/LITTLE_MAN_EXAMPLE.h"
//#include "Simple_Muscle/SIMPLE_MUSCLE_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

int main(int argc,char** argv)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;

    bool opt_arb=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-arb",&opt_arb,"Use arb test");
    parse_args.Parse(true);

    if(opt_arb) example=new ARB_EXAMPLE<float,float>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(parse_args);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
