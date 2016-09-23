//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,1> TV;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    STANDARD_TESTS<T> example(stream_type,parse_args);
    //TEST_EXAMPLE<float> example(stream_type,example_number);
    //KINEMATIC_EXAMPLE<float> example(stream_type,example_number);
    //ARB_EXAMPLE<float> example(stream_type,example_number);
    example.mpi_world=new MPI_WORLD(parse_args);
    example.After_Construction();

    if(example.mpi_world->initialized) example.solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example.Adjust_Output_Directory_For_MPI(example.solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_DRIVER<TV> driver(example);
    driver.Execute_Main_Program();

    return 0;
}
