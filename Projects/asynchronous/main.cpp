//#####################################################################
// Copyright 2008, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
    typedef float RW;
#else
    typedef float T;
    typedef float RW;
#endif
    typedef VECTOR<T,3> TV;

    STREAM_TYPE stream_type((RW()));

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;

    example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(argc,argv);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    dynamic_cast<STANDARD_TESTS<T>*>(example)->Set_Driver(&driver);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################
