//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    SOLIDS_EXAMPLE<TV>* example;

    PARSE_ARGS parse_args(argc,argv);
    example=new STANDARD_TESTS<T>(stream_type,parse_args);
    example->After_Construction();

    SOLIDS_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################

