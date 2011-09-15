//#####################################################################
// Copyright 2007, Craig Schroeder, Jon Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_DRIVER.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    SOLIDS_EXAMPLE<TV>* example;

    example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(argc,argv);

    SOLIDS_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################
