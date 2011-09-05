//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;

    example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(argc,argv);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################
