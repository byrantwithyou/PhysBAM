//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    STANDARD_TESTS<T> example(stream_type,parse_args);
    //TEST_EXAMPLE<float> example(stream_type,example_number);
    //KINEMATIC_EXAMPLE<float> example(stream_type,example_number);
    //ARB_EXAMPLE<float> example(stream_type,example_number);
    example.mpi_world=new MPI_WORLD(parse_args);
    example.After_Construction();

    SOLIDS_DRIVER<TV> driver(example);
    driver.Execute_Main_Program();

    return 0;
}
