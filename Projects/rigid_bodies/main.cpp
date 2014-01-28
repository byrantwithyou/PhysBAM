//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include "MPI_Example/MPI_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    EXAMPLE<TV>* example;

    bool opt_testmpi=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-testmpi",&opt_testmpi,"Use testmpi test");
    parse_args.Parse(true);

    if(opt_testmpi) example=new MPI_EXAMPLE<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(parse_args);

    SOLIDS_EXAMPLE<TV>* solid_fluid_example=dynamic_cast<SOLIDS_EXAMPLE<TV>*>(example);
    SOLIDS_DRIVER<TV> driver(*solid_fluid_example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
//#####################################################################
