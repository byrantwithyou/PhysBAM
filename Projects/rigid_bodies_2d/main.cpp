//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    STANDARD_TESTS<T> example(stream_type);
    //TEST_EXAMPLE<float> example(stream_type,example_number);
    //KINEMATIC_EXAMPLE<float> example(stream_type,example_number);
    //ARB_EXAMPLE<float> example(stream_type,example_number);
    example.want_mpi_world=true;
    example.Parse(argc,argv);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(example);
    driver.Execute_Main_Program();

    return 0;
}
